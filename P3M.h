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
	P3M(System &syst, uint NM, uint assignment_OP);
	virtual ~P3M();

	void print_energy();

private:
	void _assign_dipole_density();

	System &_syst;
	uint _NM;
	uint _assignment_OP;
	number _lattice_spacing;

	fftw_real *_dipole_density[3];
	rfftwnd_plan _fftw_plan;
};

#endif /* P3M_H_ */
