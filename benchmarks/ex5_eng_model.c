#include "../benchmarks/ex5_eng_model.h"

#define M_PI 3.14159265358979323846
#define DIM 3

static bool eng_print = true;
static bool small_model = true;
static bool model_mod = true;
static bool diag = false;
static bool ghost_bc = true;

// RK4 temporary storage arrays
static double *k1,*k2,*k3,*k4;

static unsigned int Nghost_top = 0;
static unsigned int Nghost_bar = 0;
static unsigned int Nghost_core = 0;

static unsigned int n_bar;
static unsigned int n_connector_bar_center;
static unsigned int n_connector_heat_sink;
static unsigned int n_connector_inlet;
static unsigned int n_connector_top;
static unsigned int n_cross_1, n_cross_2;
static unsigned int n_heat_sink;
static unsigned int n_inlet;
static unsigned int n_top;

static unsigned int n_air_surf;
static unsigned int n_bar_surf;
static unsigned int n_con_surf;

static unsigned int n_top_inl_surf;
static unsigned int n_inl_top_surf;
static unsigned int n_bot_inl_surf;
static unsigned int n_inl_bot_surf;

static double VOL_INLET;
static double VOL_CON_SURF;
static double VOL_BAR_SURF;
static double VOL_AIR_SURF;
static double AREA_AIR_SURF;
static double AREA_SURF_TOP;
static double AREA_SURF_BOT;
static double VOL_SURF_TOP;
static double VOL_SURF_BOT;
static double VOL_CON_TOP;
static double VOL_CON_INL;

static bool is_in_bounding_box(particle* particles, unsigned int i, double xl, double xh, double yl, double yh, double zl, double zh) {
	double xp = particles[i].px;
	double yp = particles[i].py;
	double zp = particles[i].pz;
	return ((xp>xl && xp<xh) && (yp>yl && yp<yh) && (zp>zl && zp<zh));
}

static double compute_temp_inlet_top(particle* particles, unsigned int Ntot) {
	double sum_T_bnd = 0.0;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].inl_top_surf) continue;
		sum_T_bnd += particles[i].f;
	}
	return (sum_T_bnd/n_inl_top_surf);
}

static double compute_temp_inlet_bot(particle* particles, unsigned int Ntot) {
	double sum_T_bnd = 0.0;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].inl_bot_surf) continue;
		sum_T_bnd += particles[i].f;
	}
	return (sum_T_bnd/n_inl_bot_surf);
}

static double compute_temp_con_top(particle* particles, unsigned int Ntot) {
	double sum_T_bnd = 0.0;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].top_inl_surf) continue;
		sum_T_bnd += particles[i].f;
	}
	return (sum_T_bnd/n_top_inl_surf);
}

static double compute_temp_con_inlet(particle* particles, unsigned int Ntot) {
	double sum_T_bnd = 0.0;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bot_inl_surf) continue;
		sum_T_bnd += particles[i].f;
	}
	return (sum_T_bnd/n_bot_inl_surf);
}

static double compute_temp_bearing_surf(particle* particles, unsigned int Ntot) {
	// compute the average temperature of ball_bearing surfaces
	double sum_T_bnd = 0.0;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].con_surf) continue;
		sum_T_bnd += particles[i].f;
	}
	return (sum_T_bnd/n_con_surf);
}

static double compute_temp_bar_surf(particle* particles, unsigned int Ntot) {
	double sum_T_bnd = 0.0;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bar_surf) continue;
		sum_T_bnd += particles[i].f;
	}
	return (sum_T_bnd/n_bar_surf);
}

static void apply_boundary_inlet_top(particle* particles, double dt, double Qbnd, unsigned int Ntot) {
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].inl_top_surf) continue;
		particles[i].f -= dt*Qbnd;
	}
}

static void apply_boundary_inlet_bot(particle* particles, double dt, double Qbnd, unsigned int Ntot) {
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].inl_bot_surf) continue;
		particles[i].f -= dt*Qbnd;
	}
}

static void apply_boundary_connector_top(particle* particles, double dt, double Qbnd, unsigned int Ntot) {
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].top_inl_surf) continue;
		particles[i].f -= dt*Qbnd;
	}
}

static void apply_boundary_connector_inlet(particle* particles, double dt, double Qbnd, unsigned int Ntot) {
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bot_inl_surf) continue;
		particles[i].f -= dt*Qbnd;
	}
}

static void apply_boundary_bearing_bar(particle* particles, double dt, double Qbnd, unsigned int Ntot) {
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].con_surf) continue;
		particles[i].f -= dt*Qbnd;
	}
}

static void apply_boundary_bar_bearing(particle* particles, double dt, double Qbnd, unsigned int Ntot) {
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bar_surf) continue;
		particles[i].f -= dt*Qbnd;
	}
}

double perform_eng_heat_euler_3D(particle* particles, METHOD method, double dt, unsigned int step) {
	singleton_geometry *geometry = get_singleton_geometry();
	unsigned int Ntot = geometry->n;

	double Tb = 22.0;       // ambient temperature
	double ha = 5.0;        // convection coefficient (W/m2.k)
	double power = 1.3;     // inlet power in (W)
	double delta = 1e-4;    // temperature difference over boundaries
	double Rth = 1.0/0.46;  // thermal resistance (k/W)
	double Rth1 = 0.46;     // thermal conductance (W/k)
	double hc_top = 0.5*(KA_ALU+KA_PUR); // ALU-PUR conductance coefficient (W/m.k)
	double hc_bot = 0.5*(KA_ALU+KA_ALU); // ALU-ALU conductance coefficient (W/m.k)

	double totalT = 0.0;
	double Tmin = +DBL_MAX;
	double Tmax = -DBL_MAX;

	//***********************************
	// assume a constant uniform volume
	//***********************************
	double dy = pow(particles[0].m/particles[0].rho,1.0/3.0);
	double dV = particles[0].m/particles[0].rho;
	double dA = pow(dV,2.0/3.0);
	double dV1 = 1./dV;
	double dV_air_surf1 = 1./VOL_AIR_SURF;
	double dV_inl_top1 = 1./VOL_SURF_TOP;
	double dV_inl_bot1 = 1./VOL_SURF_BOT;
	double dV_con_top1 = 1./VOL_CON_TOP;
	double dV_con_inl1 = 1./VOL_CON_INL;
	double dV_inl1 = 1./VOL_INLET;
	double dV_con_bar_surf1 = 1./VOL_CON_SURF;
	double dV_bar_surf1 = 1./VOL_BAR_SURF;

	double H_air     = (ha*AREA_AIR_SURF)/(RHO_STL*CP_STL);
	double H_inl_top = (hc_top*AREA_SURF_TOP)/(dy);
	double H_inl_bot = (hc_bot*AREA_SURF_BOT)/(dy);
	double H_bearing = (Rth1)/(RHO_STL*CP_STL);
	//*********************************************************

	// 1. add INLET HEAT CHARGE
	double Q_charge = ((power)/VOL_INLET)/(RHO_ALU*CP_ALU); // dim: (k/sec)
	for (unsigned int i = 0; i < Ntot; i++) {
		if(particles[i].label!=INLET) continue;

		particles[i].f += dt*Q_charge;
	}

	// 2. compute AVERAGE
	double T_avg_inlet_top = compute_temp_inlet_top(particles,Ntot);
	double T_avg_inlet_bot = compute_temp_inlet_bot(particles,Ntot);
	double T_avg_con_top   = compute_temp_con_top(particles,Ntot);
	double T_avg_con_inlet = compute_temp_con_inlet(particles,Ntot);
	double T_avg_bearing   = compute_temp_bearing_surf(particles,Ntot);
	double T_avg_bar_surf  = compute_temp_bar_surf(particles,Ntot);

	// 3. compute HEAT FLUX
	double dTdt_inl_top = (1./(RHO_ALU*CP_ALU))*H_inl_top*dV_inl_top1*(T_avg_inlet_top - T_avg_con_top);
	double dTdt_top_inl = (1./(RHO_PUR*CP_PUR))*H_inl_top*dV_con_top1*(T_avg_con_top - T_avg_inlet_top);
	double dTdt_inl_bot = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_inl_bot1*(T_avg_inlet_bot - T_avg_con_inlet);
	double dTdt_bot_inl = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_con_inl1*(T_avg_con_inlet - T_avg_inlet_bot);
	double dTdt_bearing_bar =                   H_bearing*dV_con_bar_surf1*(T_avg_bearing - T_avg_bar_surf);
	double dTdt_bar_bearing =                   H_bearing*dV_bar_surf1*(T_avg_bar_surf - T_avg_bearing);

	// 4. apply BOUNDARY CONDITIONS
	apply_boundary_inlet_top(particles,dt,dTdt_inl_top,Ntot);
	apply_boundary_inlet_bot(particles,dt,dTdt_inl_bot,Ntot);
	apply_boundary_connector_top(particles,dt,dTdt_top_inl,Ntot);
	apply_boundary_connector_inlet(particles,dt,dTdt_bot_inl,Ntot);
	apply_boundary_bearing_bar(particles,dt,dTdt_bearing_bar,Ntot);
	apply_boundary_bar_bearing(particles,dt,dTdt_bar_bearing,Ntot);

	// 5. compute LAPLACIAN
	laplacian_meshfree_method(particles,method);

	// 6. integrate HEAT EQUATION
	for (unsigned int i = 0; i < Ntot; i++) {
		if(particles[i].bnd || particles[i].label==HEAT_SINK || particles[i].label==CONNECTOR_HEAT_SINK ||
				particles[i].label==CROSS_1 || particles[i].label==CROSS_2 || particles[i].label==TOP) continue;

		double Ta = particles[i].f;
		double Q_air = (particles[i].air_surf) ? H_air*dV_air_surf1*(Ta-Tb) : 0.0;

		// main heat equation-----------------------------------------------------------------------------------
		particles[i].f += dt*(particles[i].alpha*particles[i].LaplF - Q_air);

		assert(particles[i].f >= 21.0 && particles[i].f <= 26.0);

		Tmin = fmin(Tmin,particles[i].f);
		Tmax = fmax(Tmax,particles[i].f);

		// total concentration
		totalT += particles[i].f;
	}

	if (eng_print && step%10==0) {
		printf("step=%d t=%0.3f sum=%0.3f Tmin=%f Tmax=%f\n",step, step*dt, totalT, Tmin, Tmax);
		if(diag) {
			printf("T_avg_inlet_top=%0.3f T_avg_con_top=%0.3f\n",T_avg_inlet_top,T_avg_con_top);
			printf("T_avg_inlet_bot=%0.3f T_avg_con_inlet=%0.3f\n",T_avg_inlet_bot,T_avg_con_inlet);
			printf("T_avg_bearing=%0.3f T_avg_bar_surf=%0.3f\n",T_avg_bearing,T_avg_bar_surf);
			printf("dTdt_inl_top=%0.4f  dTdt_top_inl=%0.4f\n",dTdt_inl_top,dTdt_top_inl);
			printf("dTdt_inl_bot=%0.4f  dTdt_bot_inl=%0.4f\n",dTdt_inl_bot,dTdt_bot_inl);
			printf("dTdt_bearing_bar=%0.4f  dTdt_bar_bearing=%0.4f\n",dTdt_bearing_bar,dTdt_bar_bearing);
		}
	}

	return totalT;
}

double perform_eng_heat_rk4_3D(particle* particles, METHOD method, double dt, unsigned int step) {
	singleton_geometry *geometry = get_singleton_geometry();
	unsigned int Ntot = geometry->n;
	double dt2 = 0.5*dt;

	k1 = (double*) calloc(Ntot, sizeof(double));
	k2 = (double*) calloc(Ntot, sizeof(double));
	k3 = (double*) calloc(Ntot, sizeof(double));
	k4 = (double*) calloc(Ntot, sizeof(double));

	double Tb = 22.0;       // ambient temperature
	double ha = 5.0;        // convection coefficient (W/m2.k)
	double power = 1.3;     // inlet power in (W)
	double delta = 1e-4;    // temperature difference over boundaries
	double Rth = 1.0/0.46;  // thermal resistance (k/W)
	double Rth1 = 0.46;     // thermal conductance (W/k)
	double hc_top = 0.5*(KA_ALU+KA_PUR); // ALU-PUR conductance coefficient (W/m.k)
	double hc_bot = 0.5*(KA_ALU+KA_ALU); // ALU-ALU conductance coefficient (W/m.k)

	double totalT = 0.0;
	double Tmin = +DBL_MAX;
	double Tmax = -DBL_MAX;

	//***********************************
	// assume a constant uniform volume
	//***********************************
	double dy = pow(particles[0].m/particles[0].rho,1.0/3.0);
	double dV = particles[0].m/particles[0].rho;
	double dA = pow(dV,2.0/3.0);
	double dV1 = 1./dV;
	double dV_air_surf1 = 1./VOL_AIR_SURF;
	double dV_inl_top1 = 1./VOL_SURF_TOP;
	double dV_inl_bot1 = 1./VOL_SURF_BOT;
	double dV_con_top1 = 1./VOL_CON_TOP;
	double dV_con_inl1 = 1./VOL_CON_INL;
	double dV_inl1 = 1./VOL_INLET;
	double dV_con_bar_surf1 = 1./VOL_CON_SURF;
	double dV_bar_surf1 = 1./VOL_BAR_SURF;

	double H_air     = (ha*AREA_AIR_SURF)/(RHO_STL*CP_STL);
	double H_inl_top = (hc_top*AREA_SURF_TOP)/(dy);
	double H_inl_bot = (hc_bot*AREA_SURF_BOT)/(dy);
	double H_bearing = (Rth1)/(RHO_STL*CP_STL);
	//*********************************************************

	// 0. add INLET HEAT CHARGE
	double Q_charge = ((power)/VOL_INLET)/(RHO_ALU*CP_ALU); // dim: (k/sec)
	for (unsigned int i = 0; i < Ntot; i++) {
		if(particles[i].label!=INLET) continue;
		particles[i].f += dt*Q_charge;
	}

	// 1. pre-assumption for temporary field
	for (unsigned int i = 0; i < Ntot; i++) {
		particles[i].T = particles[i].f; // boundary treatmnent
		particles[i].F = particles[i].f; // RK4 temporary field
	}

	double T_avg_inlet_top,T_avg_inlet_bot,T_avg_con_top,T_avg_con_inlet,T_avg_bearing,T_avg_bar_surf;
	double dTdt_inl_top,dTdt_top_inl,dTdt_inl_bot,dTdt_bot_inl,dTdt_bearing_bar,dTdt_bar_bearing;

	/*
	 ************
	1st RK4 call
	 ************
	 */
	// 2. compute AVERAGE
	T_avg_inlet_top = compute_temp_inlet_top(particles,Ntot);
	T_avg_inlet_bot = compute_temp_inlet_bot(particles,Ntot);
	T_avg_con_top   = compute_temp_con_top(particles,Ntot);
	T_avg_con_inlet = compute_temp_con_inlet(particles,Ntot);
	T_avg_bearing   = compute_temp_bearing_surf(particles,Ntot);
	T_avg_bar_surf  = compute_temp_bar_surf(particles,Ntot);
	// 3. compute HEAT FLUX
	dTdt_inl_top = (1./(RHO_ALU*CP_ALU))*H_inl_top*dV_inl_top1*(T_avg_inlet_top-T_avg_con_top);
	dTdt_top_inl = (1./(RHO_PUR*CP_PUR))*H_inl_top*dV_con_top1*(T_avg_con_top-T_avg_inlet_top);
	dTdt_inl_bot = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_inl_bot1*(T_avg_inlet_bot-T_avg_con_inlet);
	dTdt_bot_inl = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_con_inl1*(T_avg_con_inlet-T_avg_inlet_bot);
	dTdt_bearing_bar = H_bearing*dV_con_bar_surf1*(T_avg_bearing-T_avg_bar_surf);
	dTdt_bar_bearing = H_bearing*dV_bar_surf1*(T_avg_bar_surf-T_avg_bearing);
	// 4. apply BOUNDARY CONDITIONS
	apply_boundary_inlet_top(particles,dt,dTdt_inl_top,Ntot);
	apply_boundary_inlet_bot(particles,dt,dTdt_inl_bot,Ntot);
	apply_boundary_connector_top(particles,dt,dTdt_top_inl,Ntot);
	apply_boundary_connector_inlet(particles,dt,dTdt_bot_inl,Ntot);
	apply_boundary_bearing_bar(particles,dt,dTdt_bearing_bar,Ntot);
	apply_boundary_bar_bearing(particles,dt,dTdt_bar_bearing,Ntot);
	// 5. compute LAPLACIAN
	laplacian_meshfree_method(particles,method);
	// 6. integrate HEAT EQUATION
	for (unsigned int i = 0; i < Ntot; i++) {
		if(particles[i].bnd || particles[i].label==HEAT_SINK || particles[i].label==TOP) continue;
		double Ta = particles[i].f;
		double Q_air = (particles[i].air_surf) ? H_air*dV_air_surf1*(Ta-Tb) : 0.0;
		k1[i] = (particles[i].alpha*particles[i].LaplF - Q_air);
		particles[i].f = particles[i].F + dt2*k1[i];
	}

	/*
	 ************
	2nd RK4 call
	 ************
	 */
	// 2. compute AVERAGE
	T_avg_inlet_top = compute_temp_inlet_top(particles,Ntot);
	T_avg_inlet_bot = compute_temp_inlet_bot(particles,Ntot);
	T_avg_con_top   = compute_temp_con_top(particles,Ntot);
	T_avg_con_inlet = compute_temp_con_inlet(particles,Ntot);
	T_avg_bearing   = compute_temp_bearing_surf(particles,Ntot);
	T_avg_bar_surf  = compute_temp_bar_surf(particles,Ntot);
	// 3. compute HEAT FLUX
	dTdt_inl_top = (1./(RHO_ALU*CP_ALU))*H_inl_top*dV_inl_top1*(T_avg_inlet_top-T_avg_con_top);
	dTdt_top_inl = (1./(RHO_PUR*CP_PUR))*H_inl_top*dV_con_top1*(T_avg_con_top-T_avg_inlet_top);
	dTdt_inl_bot = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_inl_bot1*(T_avg_inlet_bot-T_avg_con_inlet);
	dTdt_bot_inl = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_con_inl1*(T_avg_con_inlet-T_avg_inlet_bot);
	dTdt_bearing_bar = H_bearing*dV_con_bar_surf1*(T_avg_bearing-T_avg_bar_surf);
	dTdt_bar_bearing = H_bearing*dV_bar_surf1*(T_avg_bar_surf-T_avg_bearing);
	// 4. apply BOUNDARY CONDITIONS
	apply_boundary_inlet_top(particles,dt,dTdt_inl_top,Ntot);
	apply_boundary_inlet_bot(particles,dt,dTdt_inl_bot,Ntot);
	apply_boundary_connector_top(particles,dt,dTdt_top_inl,Ntot);
	apply_boundary_connector_inlet(particles,dt,dTdt_bot_inl,Ntot);
	apply_boundary_bearing_bar(particles,dt,dTdt_bearing_bar,Ntot);
	apply_boundary_bar_bearing(particles,dt,dTdt_bar_bearing,Ntot);
	// 5. compute LAPLACIAN
	laplacian_meshfree_method(particles,method);
	// 6. integrate HEAT EQUATION
	for (unsigned int i = 0; i < Ntot; i++) {
		if(particles[i].bnd || particles[i].label==HEAT_SINK || particles[i].label==TOP) continue;
		double Ta = particles[i].f;
		double Q_air = (particles[i].air_surf) ? H_air*dV_air_surf1*(Ta-Tb) : 0.0;
		k2[i] = (particles[i].alpha*particles[i].LaplF - Q_air);
		particles[i].f = particles[i].F + dt2*k2[i];
	}

	/*
	 ************
	3rd RK4 call
	 ************
	 */
	// 2. compute AVERAGE
	T_avg_inlet_top = compute_temp_inlet_top(particles,Ntot);
	T_avg_inlet_bot = compute_temp_inlet_bot(particles,Ntot);
	T_avg_con_top   = compute_temp_con_top(particles,Ntot);
	T_avg_con_inlet = compute_temp_con_inlet(particles,Ntot);
	T_avg_bearing   = compute_temp_bearing_surf(particles,Ntot);
	T_avg_bar_surf  = compute_temp_bar_surf(particles,Ntot);
	// 3. compute HEAT FLUX
	dTdt_inl_top = (1./(RHO_ALU*CP_ALU))*H_inl_top*dV_inl_top1*(T_avg_inlet_top-T_avg_con_top);
	dTdt_top_inl = (1./(RHO_PUR*CP_PUR))*H_inl_top*dV_con_top1*(T_avg_con_top-T_avg_inlet_top);
	dTdt_inl_bot = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_inl_bot1*(T_avg_inlet_bot-T_avg_con_inlet);
	dTdt_bot_inl = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_con_inl1*(T_avg_con_inlet-T_avg_inlet_bot);
	dTdt_bearing_bar = H_bearing*dV_con_bar_surf1*(T_avg_bearing-T_avg_bar_surf);
	dTdt_bar_bearing = H_bearing*dV_bar_surf1*(T_avg_bar_surf-T_avg_bearing);
	// 4. apply BOUNDARY CONDITIONS
	apply_boundary_inlet_top(particles,dt,dTdt_inl_top,Ntot);
	apply_boundary_inlet_bot(particles,dt,dTdt_inl_bot,Ntot);
	apply_boundary_connector_top(particles,dt,dTdt_top_inl,Ntot);
	apply_boundary_connector_inlet(particles,dt,dTdt_bot_inl,Ntot);
	apply_boundary_bearing_bar(particles,dt,dTdt_bearing_bar,Ntot);
	apply_boundary_bar_bearing(particles,dt,dTdt_bar_bearing,Ntot);
	//	// 5. compute LAPLACIAN
	laplacian_meshfree_method(particles,method);
	// 6. integrate HEAT EQUATION
	for (unsigned int i = 0; i < Ntot; i++) {
		if(particles[i].bnd || particles[i].label==HEAT_SINK || particles[i].label==TOP) continue;
		double Ta = particles[i].f;
		double Q_air = (particles[i].air_surf) ? H_air*dV_air_surf1*(Ta-Tb) : 0.0;
		k3[i] = (particles[i].alpha*particles[i].LaplF - Q_air);
		particles[i].f = particles[i].F + dt*k3[i];
	}

	/*
	 ************
	4th RK4 call
	 ************
	 */
	// 2. compute AVERAGE
	T_avg_inlet_top = compute_temp_inlet_top(particles,Ntot);
	T_avg_inlet_bot = compute_temp_inlet_bot(particles,Ntot);
	T_avg_con_top   = compute_temp_con_top(particles,Ntot);
	T_avg_con_inlet = compute_temp_con_inlet(particles,Ntot);
	T_avg_bearing   = compute_temp_bearing_surf(particles,Ntot);
	T_avg_bar_surf  = compute_temp_bar_surf(particles,Ntot);
	// 3. compute HEAT FLUX
	dTdt_inl_top = (1./(RHO_ALU*CP_ALU))*H_inl_top*dV_inl_top1*(T_avg_inlet_top-T_avg_con_top);
	dTdt_top_inl = (1./(RHO_PUR*CP_PUR))*H_inl_top*dV_con_top1*(T_avg_con_top-T_avg_inlet_top);
	dTdt_inl_bot = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_inl_bot1*(T_avg_inlet_bot-T_avg_con_inlet);
	dTdt_bot_inl = (1./(RHO_ALU*CP_ALU))*H_inl_bot*dV_con_inl1*(T_avg_con_inlet-T_avg_inlet_bot);
	dTdt_bearing_bar = H_bearing*dV_con_bar_surf1*(T_avg_bearing-T_avg_bar_surf);
	dTdt_bar_bearing = H_bearing*dV_bar_surf1*(T_avg_bar_surf-T_avg_bearing);
	// 4. apply BOUNDARY CONDITIONS
	apply_boundary_inlet_top(particles,dt,dTdt_inl_top,Ntot);
	apply_boundary_inlet_bot(particles,dt,dTdt_inl_bot,Ntot);
	apply_boundary_connector_top(particles,dt,dTdt_top_inl,Ntot);
	apply_boundary_connector_inlet(particles,dt,dTdt_bot_inl,Ntot);
	apply_boundary_bearing_bar(particles,dt,dTdt_bearing_bar,Ntot);
	apply_boundary_bar_bearing(particles,dt,dTdt_bar_bearing,Ntot);
	// 5. compute LAPLACIAN
	laplacian_meshfree_method(particles,method);
	// 6. integrate HEAT EQUATION
	for (unsigned int i = 0; i < Ntot; i++) {
		if(particles[i].bnd || particles[i].label==HEAT_SINK || particles[i].label==TOP) continue;
		double Ta = particles[i].f;
		double Q_air = (particles[i].air_surf) ? H_air*dV_air_surf1*(Ta-Tb) : 0.0;
		k4[i] = (particles[i].alpha*particles[i].LaplF - Q_air);
	}

	// update TEMPERATURE-------------------------------------------------------------------------
	for (unsigned int i = 0; i < Ntot; i++) {
		if(particles[i].bnd || particles[i].label==HEAT_SINK || particles[i].label==TOP) continue;

		// RK4 sum
		particles[i].f = particles[i].F + (1./6.)*dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

		assert(particles[i].f >= 21.0 && particles[i].f <= 26.0);

		Tmin = fmin(Tmin,particles[i].f);
		Tmax = fmax(Tmax,particles[i].f);

		// total concentration
		totalT += particles[i].f;
	}

	if (eng_print && step%10==0) {
		printf("step=%d t=%0.3f sum=%0.3f Tmin=%f Tmax=%f\n",step, step*dt, totalT, Tmin, Tmax);
		if(diag) {
			printf("T_avg_inlet_top=%0.3f T_avg_con_top=%0.3f\n",T_avg_inlet_top,T_avg_con_top);
			printf("T_avg_inlet_bot=%0.3f T_avg_con_inlet=%0.3f\n",T_avg_inlet_bot,T_avg_con_inlet);
			printf("T_avg_bearing=%0.3f T_avg_bar_surf=%0.3f\n",T_avg_bearing,T_avg_bar_surf);
			printf("dTdt_inl_top=%0.4f  dTdt_top_inl=%0.4f\n",dTdt_inl_top,dTdt_top_inl);
			printf("dTdt_inl_bot=%0.4f  dTdt_bot_inl=%0.4f\n",dTdt_inl_bot,dTdt_bot_inl);
			printf("dTdt_bearing_bar=%0.4f  dTdt_bar_bearing=%0.4f\n",dTdt_bearing_bar,dTdt_bar_bearing);
		}
	}

	free(k1);
	free(k2);
	free(k3);
	free(k4);

	return totalT;
}

static void init_top(particle *particles, double hdy, double *dx1) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/top_mod.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/top_mod.txt", "r");
	else
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/top.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/top.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = 0;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	double max_x = -DBL_MAX;
	double max_y = -DBL_MAX;
	double max_z = -DBL_MAX;
	double min_x = +DBL_MAX;
	double min_y = +DBL_MAX;
	double min_z = +DBL_MAX;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		// to find the surface particles of BAR in contact with ambient
		max_x  = fmax(max_x,particles[ip_start + lines].px);
		max_y  = fmax(max_y,particles[ip_start + lines].py);
		max_z  = fmax(max_z,particles[ip_start + lines].pz);
		min_x  = fmin(min_x,particles[ip_start + lines].px);
		min_y  = fmin(min_y,particles[ip_start + lines].py);
		min_z  = fmin(min_z,particles[ip_start + lines].pz);

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_top = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_top-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx1 = dy;

	unsigned int iter = 0;
	for (unsigned int i = 0; i < n_top; i++) {
		particles[ip_start + i].label = TOP;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_ALU;
		particles[ip_start + i].m = RHO_ALU*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_ALU/(RHO_ALU*CP_ALU);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;

		// to find the surface particles of TOP in contact with ambient
		if(particles[ip_start + i].label==TOP &&
				(particles[ip_start + i].px==min_x || particles[ip_start + i].px==max_x ||
						particles[ip_start + i].py==min_y || particles[ip_start + i].py==max_y ||
						particles[ip_start + i].pz==min_z || particles[ip_start + i].pz==max_z)) {

			particles[ip_start + i].air_surf = true;
			iter++;
		}

	}

	n_air_surf += iter;
}

static void init_connector_top(particle *particles, double hdy, double *dx2) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/connector_top_mod.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/connector_top_mod.txt", "r");
	else
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/connector_top.txt", "r");
		else		     file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/connector_top.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = n_top;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_connector_top = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_connector_top-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx2 = dy;

	double xl = +0.01771 - 0.5*dy;
	double xh = +0.05334 + 0.5*dy;
	double yl = -DBL_MAX;
	double yh = +DBL_MAX;
	double zl = -DBL_MAX;
	double zh = +0.08757 + 0.5*dy;

	double xlm = +0.01734 - 0.5*dy;
	double xhm = +0.05142 + 0.5*dy;
	double ylm = -DBL_MAX;
	double yhm = +DBL_MAX;
	double zlm = +0.06562 - 0.5*dy;
	double zhm = +0.08550 + 0.5*dy;

	unsigned int iter1 = 0;
	unsigned int iter2 = 0;

	for (unsigned int i = 0; i < n_connector_top; i++) {
		particles[ip_start + i].label = CONNECTOR_TOP;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_PUR;
		particles[ip_start + i].m = RHO_PUR*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_PUR/(RHO_PUR*CP_PUR);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;

		// original CAD model
		if(!model_mod && is_in_bounding_box(particles,ip_start+i,xl,xh,yl,yh,zl,zh)) {
			particles[ip_start + i].top_inl_surf = true;
			iter1++;
		}

		// modified CAD model
		if(model_mod && is_in_bounding_box(particles,ip_start+i,xlm,xhm,ylm,yhm,zlm,zhm)) {
			particles[ip_start + i].top_inl_surf = true;
			iter2++;
		}
	}

	n_top_inl_surf = (model_mod) ? iter2 : iter1;
	VOL_CON_TOP = n_top_inl_surf*dy*dy*dy;
}

static void init_inlet(particle *particles, double hdy, double *dx3) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/inlet_mod.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/inlet_mod.txt", "r");
	else
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/inlet.txt", "r");
		else 			 file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/inlet.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = n_top+n_connector_top;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_inlet = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_inlet-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx3 = dy;
	VOL_INLET = n_inlet*dy*dy*dy;

	double xl1 = +0.02217 + 0.5*dy;
	double xh1 = +0.04889 - 0.5*dy;
	double yl1 = -DBL_MAX;
	double yh1 = +DBL_MAX;
	double zl1 = +0.06754 - 0.5*dy;
	double zh1 = +0.08312 - 0.5*dy;

	double xl2 = -DBL_MAX;
	double xh2 = +DBL_MAX;
	double yl2 = -DBL_MAX;
	double yh2 = +DBL_MAX;
	double zl2 = -DBL_MAX;
	double zh2 = +0.06976 + 0.5*dy;

	double xlm = +0.02302 - 0.5*dy;
	double xhm = +0.04574 + 0.5*dy;
	double ylm = -DBL_MAX;
	double yhm = +DBL_MAX;
	double zlm = +0.06846 - 0.5*dy;
	double zhm = +0.07982 + 0.5*dy;

	unsigned int iter1 = 0;
	unsigned int iter2 = 0;
	unsigned int iter3 = 0;
	unsigned int iter4 = 0;

	double xmin = +DBL_MAX; double xmax = -DBL_MAX;
	double ymin = +DBL_MAX; double ymax = -DBL_MAX;
	double zmin = +DBL_MAX; double zmax = -DBL_MAX;

	for (unsigned int i = 0; i < n_inlet; i++) {
		particles[ip_start + i].label = INLET;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_ALU;
		particles[ip_start + i].m = RHO_ALU*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_ALU/(RHO_ALU*CP_ALU);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;

		xmin = fmin(xmin,particles[ip_start + i].px); xmax = fmax(xmax,particles[ip_start + i].px);
		ymin = fmin(ymin,particles[ip_start + i].py); ymax = fmax(ymax,particles[ip_start + i].py);
		zmin = fmin(zmin,particles[ip_start + i].pz); zmax = fmax(zmax,particles[ip_start + i].pz);

		// original CAD model
		if (!model_mod && !is_in_bounding_box(particles,ip_start+i,xl1,xh1,yl1,yh1,zl1,zh1)) {
			particles[ip_start + i].inl_top_surf = true;
			iter1++;
		}
		if (!model_mod && is_in_bounding_box(particles,ip_start+i,xl2,xh2,yl2,yh2,zl2,zh2)) {
			particles[ip_start + i].inl_bot_surf = true;
			iter2++;
		}

		// modified CAD model
		if (model_mod && !is_in_bounding_box(particles,ip_start+i,xlm,xhm,ylm,yhm,zlm,zhm)) {
			particles[ip_start + i].inl_top_surf = true;
			iter3++;
		}
		if (model_mod && particles[ip_start + i].pz<(0.06864 + 0.5*dy)) {
			particles[ip_start + i].inl_bot_surf = true;
			iter4++;
		}

	}

	if (model_mod) {
		n_inl_top_surf = iter3;
		n_inl_bot_surf = iter4;
	} else {
		n_inl_top_surf = iter1;
		n_inl_bot_surf = iter2;
	}

	assert(xmax>xmin && ymax>ymin && zmax>zmin);
	AREA_SURF_TOP = (xmax-xmin)*(ymax-ymin) + 2.0*(zmax-zmin)*(ymax-ymin);
	AREA_SURF_BOT = (xmax-xmin)*(ymax-ymin);
	VOL_SURF_TOP = n_inl_top_surf*dy*dy*dy;
	VOL_SURF_BOT = n_inl_bot_surf*dy*dy*dy;
}

static void init_connector_inlet(particle *particles, double hdy, double *dx4) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/connector_inlet_mod.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/connector_inlet_mod.txt", "r");
	else
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/connector_inlet.txt", "r");
		else 			 file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/connector_inlet.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = n_top+n_connector_top+n_inlet;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_connector_inlet = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_connector_inlet-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx4 = dy;

	double xl = +0.01994 - 0.5*dy;
	double xh = +0.05111 + 0.5*dy;
	double yl = -DBL_MAX;
	double yh = +DBL_MAX;
	double zl = +0.06531 - 0.5*dy;
	double zh = +DBL_MAX;

	double xlm = +0.02018 - 0.5*dy;
	double xhm = +0.04858 + 0.5*dy;
	double ylm = -DBL_MAX;
	double yhm = +DBL_MAX;
	double zlm = +0.06562 - 0.5*dy;
	double zhm = +0.06562 + 0.5*dy;

	unsigned int iter1 = 0;
	unsigned int iter2 = 0;

	for (unsigned int i = 0; i < n_connector_inlet; i++) {
		particles[ip_start + i].label = CONNECTOR_INLET;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_ALU;
		particles[ip_start + i].m = RHO_ALU*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_ALU/(RHO_ALU*CP_ALU);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;

		// original CAD model
		if(!model_mod && is_in_bounding_box(particles,ip_start+i,xl,xh,yl,yh,zl,zh)) {
			particles[ip_start + i].bot_inl_surf = true;
			iter1++;
		}

		// modified CAD model
		if(model_mod && is_in_bounding_box(particles,ip_start+i,xlm,xhm,ylm,yhm,zlm,zhm)) {
			particles[ip_start + i].bot_inl_surf = true;
			iter2++;
		}
	}

	n_bot_inl_surf = (model_mod) ? iter2 : iter1;
	VOL_CON_INL = n_bot_inl_surf*dy*dy*dy;
}

static void init_connector_bar_center(particle *particles, double hdy, double *dx5) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/connector_bar_center_mod.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/connector_bar_center_mod.txt", "r");
	else
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/connector_bar_center.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/connector_bar_center.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = n_top+n_connector_top+n_inlet+n_connector_inlet;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_connector_bar_center = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_connector_bar_center-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx5 = dy;

	double xl = -DBL_MAX;
	double xh = +DBL_MAX;
	double yl = +0.00940 - 0.5*dy;
	double yh = +0.02499 + 0.5*dy;
	double zl = +0.02300 - 0.5*dy;
	double zh = +0.03859 + 0.5*dy;

	double xlm = -DBL_MAX;
	double xhm = +DBL_MAX;
	double ylm = +0.00882 - 0.5*dy;
	double yhm = +0.02586 + 0.5*dy;
	double zlm = +0.02586 - 0.5*dy;
	double zhm = +0.03722 + 0.5*dy;

	unsigned int iter1 = 0;
	unsigned int iter2 = 0;

	for (unsigned int i = 0; i < n_connector_bar_center; i++) {
		particles[ip_start + i].label = CONNECTOR_BAR_CENTER;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_STL;
		particles[ip_start + i].m = RHO_STL*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_STL/(RHO_STL*CP_STL);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;

		// original CAD model
		if(!model_mod && is_in_bounding_box(particles,ip_start+i,xl,xh,yl,yh,zl,zh)) {
			particles[ip_start + i].con_surf = true;
			particles[ip_start + i].T = 22.0;
			iter1++;
		}

		// modified CAD model
		if(model_mod && is_in_bounding_box(particles,ip_start+i,xlm,xhm,ylm,yhm,zlm,zhm)) {
			particles[ip_start + i].con_surf = true;
			particles[ip_start + i].T = 22.0;
			iter2++;
		}

	}

	n_con_surf = (model_mod) ? iter2 : iter1;
	VOL_CON_SURF = n_con_surf*dy*dy*dy;
}

static void init_bar(particle *particles, double hdy, double *dx6) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/bar_mod.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/bar_mod.txt", "r");
	else
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/bar.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/bar.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = n_top+n_connector_top+n_inlet+n_connector_inlet+n_connector_bar_center;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_bar = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_bar-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx6 = dy;

	double xl1 = +0.02217 + 0.5*dy;
	double xh1 = +0.04889 - 0.5*dy;
	double yl1 = -DBL_MAX;
	double yh1 = +DBL_MAX;
	double zl1 = +0.06754 - 0.5*dy;
	double zh1 = +0.08312 - 0.5*dy;

	double xl2 = -DBL_MAX;
	double xh2 = +DBL_MAX;
	double yl2 = +0.01163 + 0.5*dy;
	double yh2 = +0.02276 - 0.5*dy;
	double zl2 = +0.02078 + 0.5*dy;
	double zh2 = +0.03636 - 0.5*dy;

	double xlm = -DBL_MAX;
	double xhm = +DBL_MAX;
	double ylm = +0.01450 + 0.5*dy;
	double yhm = +0.02018 - 0.5*dy;
	double zlm = +0.02018 + 0.5*dy;
	double zhm = +0.03438 - 0.5*dy;

	unsigned int iter1 = 0;
	unsigned int iter2 = 0;
	unsigned int iter3 = 0;
	unsigned int iter4 = 0;

	for (unsigned int i = 0; i < n_bar; i++) {
		particles[ip_start + i].label = BAR;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_STL;
		particles[ip_start + i].m = RHO_STL*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_STL/(RHO_STL*CP_STL);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;

		double xp = particles[ip_start + i].px;
		double zp = particles[ip_start + i].pz;

		// original CAD model
		if(!model_mod && !is_in_bounding_box(particles,ip_start+i,xl2,xh2,yl2,yh2,zl2,zh2)) {
			particles[ip_start + i].air_surf = true;
			iter1++;
		}
		if(!model_mod && particles[ip_start + i].air_surf && (xp>0.0132 && xp<0.0556) && zp>0.024) {
			particles[ip_start + i].bar_surf = true;
			iter2++;
		}

		// modified CAD model
		if(model_mod && !is_in_bounding_box(particles,ip_start+i,xlm,xhm,ylm,yhm,zlm,zhm)) {
			particles[ip_start + i].air_surf = true;
			iter3++;
		}
		if(model_mod && particles[ip_start + i].air_surf && (xp>0.01450 && xp<0.05426) && zp>0.024) {
			particles[ip_start + i].bar_surf = true;
			iter4++;
		}

	}

	if (!model_mod) {
		n_air_surf = iter1;
		n_bar_surf = iter2;
	} else {
		n_air_surf = iter3;
		n_bar_surf = iter4;
	}

	AREA_AIR_SURF = n_air_surf*dy*dy;
	VOL_AIR_SURF = n_air_surf*dy*dy*dy;
	VOL_BAR_SURF = n_bar_surf*dy*dy*dy;
}

static void init_cross_1(particle *particles, double hdy, double *dx7) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/cross_1_mod.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/cross_1_mod.txt", "r");
	else
		if (small_model)file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/cross_1.txt", "r");
		else 	    	file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/cross_1.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = n_top+n_connector_top+n_inlet+n_connector_inlet+n_connector_bar_center+n_bar;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_cross_1 = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_cross_1-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx7 = dy;

	for (unsigned int i = 0; i < n_cross_1; i++) {
		particles[ip_start + i].label = CROSS_1;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_ALU;
		particles[ip_start + i].m = RHO_ALU*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_ALU/(RHO_ALU*CP_ALU);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;
	}
}

static void init_cross_2(particle *particles, double hdy, double *dx8) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/cross_2_mod.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/cross_2_mod.txt", "r");
	else
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/cross_2.txt", "r");
		else 			 file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/cross_2.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = n_top+n_connector_top+n_inlet+n_connector_inlet+n_connector_bar_center+n_bar+n_cross_1;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_cross_2 = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_cross_2-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx8 = dy;

	for (unsigned int i = 0; i < n_cross_2; i++) {
		particles[ip_start + i].label = CROSS_2;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_ALU;
		particles[ip_start + i].m = RHO_ALU*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_ALU/(RHO_ALU*CP_ALU);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;
	}
}

static void init_connector_heat_sink(particle *particles, double hdy, double *dx9) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/connector_heat_sink_mod.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/connector_heat_sink_mod.txt", "r");
	else
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/connector_heat_sink.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/connector_heat_sink.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = n_top+n_connector_top+n_inlet+n_connector_inlet+n_connector_bar_center+
			n_bar+n_cross_1+n_cross_2;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_connector_heat_sink = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_connector_heat_sink-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx9 = dy;

	for (unsigned int i = 0; i < n_connector_heat_sink; i++) {
		particles[ip_start + i].label = CONNECTOR_HEAT_SINK;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_ALU;
		particles[ip_start + i].m = RHO_ALU*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_ALU/(RHO_ALU*CP_ALU);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;
	}
}

static void init_heat_sink(particle *particles, double hdy, double *dx10) {
	FILE *file;

	if (model_mod)
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_smallest_simpl_heatsink/heat_sink_mod_simpl.txt", "r");
		else             file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_checked_smaller/heat_sink_mod.txt", "r");
	else
		if (small_model) file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_small/heat_sink.txt", "r");
		else 			 file = fopen("/home/mamzi/Desktop/MY_PUBLICATIONS/Model_For_Engineering_Case/ascii_large/heat_sink.txt", "r");

	if (file == NULL) {
		printf("Error: input file cannot be found!");
		exit(1);
	}

	int d;
	unsigned int ip_start = n_top+n_connector_top+n_inlet+n_connector_inlet+n_connector_bar_center+
			n_bar+n_cross_1+n_cross_2+n_connector_heat_sink;
	unsigned int lines = 0;
	unsigned int ch = 0;
	float x,y,z;

	while(!feof(file)) {

		d = fscanf(file,"%f %f %f",&x,&y,&z);

		particles[ip_start + lines].px = x * 1e-3;
		particles[ip_start + lines].py = y * 1e-3;
		particles[ip_start + lines].pz = z * 1e-3;

		ch = fgetc(file);
		if(ch == '\n')
			lines++;
	}
	n_heat_sink = lines;
	fclose(file);

	// calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	for(unsigned int i = 0; i < (n_heat_sink-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dx10 = dy;

	for (unsigned int i = 0; i < n_heat_sink; i++) {
		particles[ip_start + i].label = HEAT_SINK;
		particles[ip_start + i].bnd = false;
		particles[ip_start + i].air_surf = false;
		particles[ip_start + i].con_surf = false;
		particles[ip_start + i].bar_surf = false;
		particles[ip_start + i].top_inlet = false;
		particles[ip_start + i].bot_inlet = false;
		particles[ip_start + i].bot_inl_surf = false;
		particles[ip_start + i].top_inl_surf = false;
		particles[ip_start + i].inl_bot_surf = false;
		particles[ip_start + i].inl_top_surf = false;
		particles[ip_start + i].rho = RHO_ALU;
		particles[ip_start + i].m = RHO_ALU*dy*dy*dy;
		particles[ip_start + i].h = hdy*dy;
		particles[ip_start + i].alpha = KA_ALU/(RHO_ALU*CP_ALU);
		particles[ip_start + i].f = 22.0;
		particles[ip_start + i].T = 0.0;
	}
}

static void make_ghost_particles_top(particle* particles, unsigned int Ntot, unsigned int Nbnd, double hdy, double *dxp) {
	// define dummy parameters
	double rho = 1.0;

	double *xghost,*yghost,*zghost;

	// 1. calculate the particle spacing (e.g. check y direction)
	unsigned int ip_start = 0;
	double dy = 0;
	for(unsigned int i = 0; i < (n_top-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dxp = dy;

	// 2. find the boundary coordinates
	double max_x = -DBL_MAX;
	double max_y = -DBL_MAX;
	double max_z = -DBL_MAX;
	double min_x = +DBL_MAX;
	double min_y = +DBL_MAX;
	double min_z = +DBL_MAX;
	for (unsigned int i = 0; i < n_top; i++) {
		max_x  = fmax(max_x,particles[ip_start + i].px);
		max_y  = fmax(max_y,particles[ip_start + i].py);
		max_z  = fmax(max_z,particles[ip_start + i].pz);
		min_x  = fmin(min_x,particles[ip_start + i].px);
		min_y  = fmin(min_y,particles[ip_start + i].py);
		min_z  = fmin(min_z,particles[ip_start + i].pz);
	}

	double lx = fabs(max_x - min_x);
	double ly = fabs(max_y - min_y);
	double lz = fabs(max_z - min_z);
	unsigned int nx = lx/dy + 1;
	unsigned int ny = ly/dy + 1;
	unsigned int nz = lz/dy + 1;
	unsigned int n_bulk = (nx+2*Nbnd)*(ny+2*Nbnd)*(nz+2*Nbnd);

	xghost = (double*) calloc(n_bulk, sizeof(double));
	yghost = (double*) calloc(n_bulk, sizeof(double));
	zghost = (double*) calloc(n_bulk, sizeof(double));

	// 3. generate bulk ghost particles
	double Lbnd = (Nbnd)*dy;
	for (unsigned int i = 0; i < (nx+2*Nbnd); i++) {
		for (unsigned int j = 0; j < (ny+2*Nbnd); j++) {
			for (unsigned int k = 0; k < (nz+2*Nbnd); k++) {
				unsigned int ID = i*(ny+2*Nbnd)*(nz+2*Nbnd) + j*(nz+2*Nbnd) + k;

				xghost[ID] = (min_x - Lbnd) + i*dy;
				yghost[ID] = (min_y - Lbnd) + j*dy;
				zghost[ID] = (min_z - Lbnd) + k*dy;
			}
		}
	}

	// 4. deduct from the main top-particles
	unsigned int ip = 0;
	for (unsigned int i = 0; i < n_bulk; i++) {
		particles[Ntot + ip].bnd = false;
		double r_min = +DBL_MAX;

		for (unsigned int j = 0; j < n_top; j++) {
			double xtop = particles[ip_start + j].px;
			double ytop = particles[ip_start + j].py;
			double ztop = particles[ip_start + j].pz;

			double rij = sqrt((xghost[i]-xtop)*(xghost[i]-xtop) + (yghost[i]-ytop)*(yghost[i]-ytop) + (zghost[i]-ztop)*(zghost[i]-ztop));

			r_min = fmin(r_min,rij);
		}

		if (r_min>=0.5*dy && !(zghost[i]<(0.1082-0.5*dy))) {
			particles[Ntot + ip].label = GHOST_TOP;
			particles[Ntot + ip].bnd = true;
			particles[Ntot + ip].air_surf = false;

			particles[Ntot + ip].px = xghost[i];
			particles[Ntot + ip].py = yghost[i];
			particles[Ntot + ip].pz = zghost[i];

			particles[Ntot + ip].rho = rho;
			particles[Ntot + ip].m = rho*dy*dy*dy;
			particles[Ntot + ip].h = hdy*dy;
			particles[Ntot + ip].alpha = 1.0;
			particles[Ntot + ip].f = 22.0;
			particles[Ntot + ip].T = 0.0;

			ip++;
		}
	}

	// 5. to find the surface particles of main body in contact with ambient
	for(unsigned int j = 0; j < n_top; j++) {
		double r_min = +DBL_MAX;

		double xtop = particles[ip_start + j].px;
		double ytop = particles[ip_start + j].py;
		double ztop = particles[ip_start + j].pz;

		for(unsigned int i = 0; i < ip; i++) {
			double xgh = particles[Ntot + i].px;
			double ygh = particles[Ntot + i].py;
			double zgh = particles[Ntot + i].pz;

			double rij = sqrt((xgh-xtop)*(xgh-xtop) + (ygh-ytop)*(ygh-ytop) + (zgh-ztop)*(zgh-ztop));

			r_min = fmin(r_min,rij);
		}

		if(r_min <= dy)
			particles[ip_start + j].air_surf = true;
	}

	free(xghost);
	free(yghost);
	free(zghost);

	Nghost_top = ip;
}

static void make_ghost_particles_core(particle* particles, unsigned int Ntot, unsigned int Nbnd, double hdy, double *dxp) {
	// define dummy parameters
	double rho = 1.0;

	double *xghost,*yghost,*zghost;

	// 1. calculate the particle spacing (e.g. check y direction)
	double dy = 0;
	unsigned int ip_start = n_top;
	unsigned int np = n_connector_top+n_inlet+n_connector_inlet+n_connector_bar_center;
	for(unsigned int i = 0; i < (np-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dxp = dy;

	// 2. find the boundary coordinates
	double max_x = -DBL_MAX;
	double max_y = -DBL_MAX;
	double max_z = -DBL_MAX;
	double min_x = +DBL_MAX;
	double min_y = +DBL_MAX;
	double min_z = +DBL_MAX;
	for (unsigned int i = 0; i < np; i++) {
		max_x  = fmax(max_x,particles[ip_start + i].px);
		max_y  = fmax(max_y,particles[ip_start + i].py);
		max_z  = fmax(max_z,particles[ip_start + i].pz);
		min_x  = fmin(min_x,particles[ip_start + i].px);
		min_y  = fmin(min_y,particles[ip_start + i].py);
		min_z  = fmin(min_z,particles[ip_start + i].pz);
	}

	double lx = fabs(max_x - min_x);
	double ly = fabs(max_y - min_y);
	double lz = fabs(max_z - min_z);
	unsigned int nx = lx/dy + 1;
	unsigned int ny = ly/dy + 1;
	unsigned int nz = lz/dy + 1;
	unsigned int n_bulk = (nx+2*Nbnd)*(ny+2*Nbnd)*(nz+2*Nbnd);

	xghost = (double*) calloc(n_bulk, sizeof(double));
	yghost = (double*) calloc(n_bulk, sizeof(double));
	zghost = (double*) calloc(n_bulk, sizeof(double));

	// 3. generate bulk ghost particles
	double Lbnd = (Nbnd)*dy;
	for (unsigned int i = 0; i < (nx+2*Nbnd); i++) {
		for (unsigned int j = 0; j < (ny+2*Nbnd); j++) {
			for (unsigned int k = 0; k < (nz+2*Nbnd); k++) {
				unsigned int ID = i*(ny+2*Nbnd)*(nz+2*Nbnd) + j*(nz+2*Nbnd) + k;

				xghost[ID] = (min_x - Lbnd) + i*dy;
				yghost[ID] = (min_y - Lbnd) + j*dy;
				zghost[ID] = (min_z - Lbnd) + k*dy;
			}
		}
	}

	// 4. deduct from the main core-particles
	unsigned int ip = 0;
	for (unsigned int i = 0; i < n_bulk; i++) {
		particles[Ntot + Nghost_top + ip].bnd = false;
		double r_min = +DBL_MAX;

		for (unsigned int j = 0; j < np; j++) {
			double xtop = particles[ip_start + j].px;
			double ytop = particles[ip_start + j].py;
			double ztop = particles[ip_start + j].pz;

			double rij = sqrt((xghost[i]-xtop)*(xghost[i]-xtop) + (yghost[i]-ytop)*(yghost[i]-ytop) + (zghost[i]-ztop)*(zghost[i]-ztop));

			r_min = fmin(r_min,rij);
		}

		if (r_min>=0.5*dy && !(zghost[i]>(0.1054+0.5*dy))) {
			particles[Ntot + Nghost_top + ip].label = GHOST_CORE;
			particles[Ntot + Nghost_top + ip].bnd = true;
			particles[Ntot + Nghost_top + ip].air_surf = false;

			particles[Ntot + Nghost_top + ip].px = xghost[i];
			particles[Ntot + Nghost_top + ip].py = yghost[i];
			particles[Ntot + Nghost_top + ip].pz = zghost[i];

			particles[Ntot + Nghost_top + ip].rho = rho;
			particles[Ntot + Nghost_top + ip].m = rho*dy*dy*dy;
			particles[Ntot + Nghost_top + ip].h = hdy*dy;
			particles[Ntot + Nghost_top + ip].alpha = 1.0;
			particles[Ntot + Nghost_top + ip].f = 22.0;
			particles[Ntot + Nghost_top + ip].T = 0.0;

			ip++;
		}
	}

	free(xghost);
	free(yghost);
	free(zghost);

	Nghost_core = ip;
}

static void make_ghost_particles_bar(particle* particles, unsigned int Ntot, unsigned int Nbnd, double hdy, double *dxp) {
	// define dummy parameters
	double rho = 1.0;

	double *xghost,*yghost,*zghost;

	// 1. calculate the particle spacing (e.g. check y direction)
	unsigned int ip_start = n_top+n_connector_top+n_inlet+n_connector_inlet+n_connector_bar_center;
	double dy = 0;
	for(unsigned int i = 0; i < (n_bar-1); i++) {
		dy = fabs(particles[ip_start + i].py - particles[ip_start + (i+1)].py);

		if (dy > 1e-3)
			break;
	}
	*dxp = dy;

	// 2. find the boundary coordinates
	double max_x = -DBL_MAX;
	double max_y = -DBL_MAX;
	double max_z = -DBL_MAX;
	double min_x = +DBL_MAX;
	double min_y = +DBL_MAX;
	double min_z = +DBL_MAX;
	for (unsigned int i = 0; i < n_bar; i++) {
		max_x  = fmax(max_x,particles[ip_start + i].px);
		max_y  = fmax(max_y,particles[ip_start + i].py);
		max_z  = fmax(max_z,particles[ip_start + i].pz);
		min_x  = fmin(min_x,particles[ip_start + i].px);
		min_y  = fmin(min_y,particles[ip_start + i].py);
		min_z  = fmin(min_z,particles[ip_start + i].pz);
	}

	double lx = fabs(max_x - min_x);
	double ly = fabs(max_y - min_y);
	double lz = fabs(max_z - min_z);
	unsigned int nx = lx/dy + 1;
	unsigned int ny = ly/dy + 1;
	unsigned int nz = lz/dy + 1;
	unsigned int n_bulk = (nx+2*Nbnd)*(ny+2*Nbnd)*(nz+2*Nbnd);

	xghost = (double*) calloc(n_bulk, sizeof(double));
	yghost = (double*) calloc(n_bulk, sizeof(double));
	zghost = (double*) calloc(n_bulk, sizeof(double));

	// 3. generate bulk ghost particles
	double Lbnd = (Nbnd)*dy;
	for (unsigned int i = 0; i < (nx+2*Nbnd); i++) {
		for (unsigned int j = 0; j < (ny+2*Nbnd); j++) {
			for (unsigned int k = 0; k < (nz+2*Nbnd); k++) {
				unsigned int ID = i*(ny+2*Nbnd)*(nz+2*Nbnd) + j*(nz+2*Nbnd) + k;

				xghost[ID] = (min_x - Lbnd) + i*dy;
				yghost[ID] = (min_y - Lbnd) + j*dy;
				zghost[ID] = (min_z - Lbnd) + k*dy;
			}
		}
	}

	// 4. deduct from the main bar-particles
	unsigned int ip = 0;
	for (unsigned int i = 0; i < n_bulk; i++) {
		particles[Ntot + Nghost_top + Nghost_core + ip].bnd = false;
		double r_min = +DBL_MAX;

		for (unsigned int j = 0; j < n_bar; j++) {
			double xbar = particles[ip_start + j].px;
			double ybar = particles[ip_start + j].py;
			double zbar = particles[ip_start + j].pz;

			double rij = sqrt((xghost[i]-xbar)*(xghost[i]-xbar) + (yghost[i]-ybar)*(yghost[i]-ybar) + (zghost[i]-zbar)*(zghost[i]-zbar));

			r_min = fmin(r_min,rij);
		}

		if (r_min>=0.5*dy) {
			particles[Ntot + Nghost_top + Nghost_core + ip].label = GHOST_BAR;
			particles[Ntot + Nghost_top + Nghost_core + ip].bnd = true;
			particles[Ntot + Nghost_top + Nghost_core + ip].air_surf = false;

			particles[Ntot + Nghost_top + Nghost_core + ip].px = xghost[i];
			particles[Ntot + Nghost_top + Nghost_core + ip].py = yghost[i];
			particles[Ntot + Nghost_top + Nghost_core + ip].pz = zghost[i];

			particles[Ntot + Nghost_top + Nghost_core + ip].rho = rho;
			particles[Ntot + Nghost_top + Nghost_core + ip].m = rho*dy*dy*dy;
			particles[Ntot + Nghost_top + Nghost_core + ip].h = hdy*dy;
			particles[Ntot + Nghost_top + Nghost_core + ip].alpha = 1.0;
			particles[Ntot + Nghost_top + Nghost_core + ip].f = 22.0;
			particles[Ntot + Nghost_top + Nghost_core + ip].T = 0.0;

			ip++;
		}
	}

	// 5. to find the surface particles of main body in contact with ambient
	for(unsigned int j = 0; j < n_bar; j++) {
		double r_min = +DBL_MAX;

		double xbar = particles[ip_start + j].px;
		double ybar = particles[ip_start + j].py;
		double zbar = particles[ip_start + j].pz;

		for(unsigned int i = 0; i < ip; i++) {
			double xgh = particles[Ntot + Nghost_top + Nghost_core + i].px;
			double ygh = particles[Ntot + Nghost_top + Nghost_core + i].py;
			double zgh = particles[Ntot + Nghost_top + Nghost_core + i].pz;

			double rij = sqrt((xgh-xbar)*(xgh-xbar) + (ygh-ybar)*(ygh-ybar) + (zgh-zbar)*(zgh-zbar));

			r_min = fmin(r_min,rij);
		}

		if(r_min <= dy)
			particles[ip_start + j].air_surf = true;
	}

	free(xghost);
	free(yghost);
	free(zghost);

	Nghost_bar = ip;
}

particle* eng_init3D(double hdy, unsigned int Nbnd, double *spacing) {
	singleton_geometry *geometry = get_singleton_geometry();
	particle* particles = (small_model) ? (particle*) calloc(100000, sizeof(particle)) : (particle*) calloc(200000,sizeof(particle));

	double dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,dx10;

	// build the model - downward
	init_top(particles,hdy,&dx1);
	init_connector_top(particles,hdy,&dx2);
	init_inlet(particles,hdy,&dx3);
	init_connector_inlet(particles,hdy,&dx4);
	init_connector_bar_center(particles,hdy,&dx5);
	init_bar(particles,hdy,&dx6);
	init_cross_1(particles,hdy,&dx7);
	init_cross_2(particles,hdy,&dx8);
	init_connector_heat_sink(particles,hdy,&dx9);
	init_heat_sink(particles,hdy,&dx10);

	*spacing = fmax(dx1,dx10);

	unsigned int n_tot = n_top+n_connector_top+n_inlet+n_connector_inlet+n_connector_bar_center+
			n_bar+n_cross_1+n_cross_2+n_connector_heat_sink+n_heat_sink;

	if(ghost_bc) {
		double dxp;
		make_ghost_particles_top(particles,n_tot,Nbnd,hdy,&dxp);
		make_ghost_particles_core(particles,n_tot,Nbnd,hdy,&dxp);
		make_ghost_particles_bar(particles,n_tot,Nbnd,hdy,&dxp);
	}
	unsigned int Nghost = Nghost_bar + Nghost_core + Nghost_top;

	printf("3D engineering simulation is running on %d physical particles + %d ghost particles\n", n_tot, Nghost);
	printf("=================================================================================================\n");

	geometry->n = n_tot + Nghost;

	return particles;
}
