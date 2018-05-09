/*
 * P3M.cpp
 *
 *  Created on: 9 May 2018
 *      Author: lorenzo
 */

#include "P3M.h"

#include "System.h"

P3M::P3M(System &syst, uint NM, uint assignment_OP) :
				_syst(syst),
				_NM(NM),
				_assignment_OP(assignment_OP) {

	_dipole_density[0] = _dipole_density[1] = _dipole_density[2] = nullptr;
	uint partial = _NM;
	while(partial > 1) {
		if(partial % 2) {
			std::cerr << "NM should be a power of two" << std::endl;
			exit(1);
		}
		partial /= 2;
	}

	_lattice_spacing = _syst.box / _NM;

	_fftw_plan = rfftw3d_create_plan(_NM, _NM, _NM, FFTW_FORWARD, FFTW_ESTIMATE);

	for(int d = 0; d < 3; d++) {
		_dipole_density[d] = new fftw_real[CUB(_NM)];
	}
}

P3M::~P3M() {
	for(int d = 0; d < 3; d++) {
		if(_dipole_density[d] != nullptr) {
			delete[] _dipole_density[d];
		}
	}

	rfftwnd_destroy_plan(_fftw_plan);
}

void P3M::_assign_dipole_density() {

}

void P3M::print_energy() {
	number E_self = -2. * _syst.N() * CUB(_syst.alpha) / (3. * sqrt(M_PI));

	number E_r = 0.;
	for(uint p = 0; p < _syst.N(); p++) {
		for(uint q = 0; q < _syst.N(); q++) {
			if(p != q) {
				vec3 r = _syst.positions[q] - _syst.positions[p];
				r -= ((r.array() / _syst.box).round() * _syst.box).matrix();

				number r_sqr = r.dot(r);

				if(r_sqr <= _syst.rcut_sqr) {
					number r_mod = sqrt(r_sqr);

					number erfc_part = erfc(_syst.alpha * r_mod) / r_mod;
					number exp_part = M_2_SQRTPI * _syst.alpha * exp(-SQR(_syst.alpha) * r_sqr);

					vec3 &p_dip = _syst.dipoles[p];
					vec3 &q_dip = _syst.dipoles[q];

					E_r += (p_dip.dot(q_dip) * (erfc_part + exp_part) - p_dip.dot(r) * q_dip.dot(r) * (3. * erfc_part / r_sqr + (2. * SQR(_syst.alpha) + 3. / r_sqr) * exp_part)) / r_sqr;
				}
			}
		}
	}
	E_r /= 2.;

	number E_k = 0.;

	number E_tot = E_r + E_k + E_self;

	std::cout << "Real part: " << E_r << std::endl;
	std::cout << "Reciprocal part: " << E_k << std::endl;
	std::cout << "Self part: " << E_self << std::endl;
	std::cout << "Total: " << E_tot << " " << E_tot / _syst.N() << std::endl;
}
