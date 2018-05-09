/*
 * Ewald.h
 *
 *  Created on: 20 Apr 2018
 *      Author: lorenzo
 */

#ifndef EWALD_H_
#define EWALD_H_

#include "defs.h"

#include <complex>

class System;

class Ewald {
public:
	Ewald(System &syst);
	virtual ~Ewald();

	void print_energy();

private:
	System &_syst;
	int _kcut;
	std::vector<vec3> _k_vectors;
	std::vector<number> _k_factors;
};

#endif /* EWALD_H_ */
