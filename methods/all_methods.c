//Copyright ETH Zurich, IWF

//This file is part of iwf_thermal.

//iwf_thermal is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//iwf_thermal is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with iwf_thermal.  If not, see <http://www.gnu.org/licenses/>.

#include "all_methods.h"
#include "sph.h"
#include "pse.h"
#include "rkpm.h"
#include "cspm.h"
#include "nmfs.h"

void initialize_meshfree_method(particle* particles, METHOD method, kernel kernel_func, unsigned int Ntot) {
	switch (method) {
	case SPH:
		sph_init(kernel_func);
		break;
	case PSE:
		break;
	case RKPM:
		rkpm_init(Ntot,kernel_func);
		break;
	case CSPM:
		cspm_init(Ntot,kernel_func);
		break;
	case NMFS:
		nmfs_init(Ntot,kernel_func);
		break;
	}
}

void correction_terms_meshfree_method(particle* particles, METHOD method) {
	switch (method) {
	case SPH:
		break;
	case PSE:
		break;
	case RKPM:
		rkpm_compute_correction_terms(particles);
		break;
	case CSPM:
		cspm_compute_correction_terms(particles);
		break;
	case NMFS:
		nmfs_compute_correction_terms(particles);
		break;
	}
}

void laplacian_meshfree_method(particle* particles, METHOD method) {
	switch (method) {
	case SPH:
		perform_sph(particles);
		break;
	case PSE:
		perform_pse(particles);
		break;
	case RKPM:
		perform_rkpm(particles);
		break;
	case CSPM:
		perform_cspm(particles);
		break;
	case NMFS:
		perform_nmfs(particles);
		break;
	}
}

void methods_wipe_out() {
	rkpm_wipe_out();
	cspm_wipe_out();
	nmfs_wipe_out();
}
