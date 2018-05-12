/*
 * P3M.h
 *
 *  Created on: 9 May 2018
 *      Author: lorenzo
 */

#ifndef P3M_H_
#define P3M_H_

#include "defs.h"

#include <rfftw.h>
#include <complex>

class System;

class P3M {
public:
	P3M(System &syst, int NM, int assignment_OP);
	virtual ~P3M();

	void print_energy();

private:
	void _assign_dipole_density();
	int _cell_index(Eigen::Vector3i v);
	std::pair<number, number> _G_factors(vec3 &k);

	System &_syst;
	int _N_mesh_side;
	int _assignment_OP;
	bool _assignment_OP_is_odd;
	number _mesh_spacing;
	number _mesh_cell_volume;
	int _N_mesh;

	std::vector<fftw_real> _dipole_density[3];
	std::vector<fftw_complex> _transformed_dipole_density[3];
	std::vector<number> _green_function;
	rfftwnd_plan _fftw_plan;
};

#endif /* P3M_H_ */
