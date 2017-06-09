
particle make_particle(unsigned int part_id, double px, double py, doube pz) {
	particle p;
	memset(&p,0,sizeof(particle));

	p.part_id = part_id;
	p.px = px;
	p.py = py;
	p.pz = pz;

	return p;
}
