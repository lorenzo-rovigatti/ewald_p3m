/*
 * Ewald.cpp
 *
 *  Created on: 20 Apr 2018
 *      Author: lorenzo
 */

#include "Ewald.h"

#include "System.h"

#include <cmath>

Ewald::Ewald(System &syst) :
				_syst(syst),
				_kcut(20) {

	number reciprocal = 2 * M_PI / _syst.box;
	for(int kx = -_kcut; kx <= _kcut; kx++) {
		for(int ky = -_kcut; ky <= _kcut; ky++) {
			for(int kz = -_kcut; kz <= _kcut; kz++) {
				if(kx != 0 || ky != 0 || kz != 0) {
					vec3 k_new = reciprocal * vec3(kx, ky, kz);
					_k_vectors.emplace_back(k_new);

					number k_sqr = k_new.dot(k_new);
					number k_factor = 4 * M_PI * exp(-k_sqr / (4 * SQR(_syst.alpha))) / (2. * k_sqr * CUB(_syst.box));
					_k_factors.emplace_back(k_factor);
				}
			}
		}
	}
}

Ewald::~Ewald() {

}

void Ewald::print_energy() {
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
	for(uint k = 0; k < _k_vectors.size(); k++) {
		std::complex<number> k_density(0., 0.);
		for(uint p = 0; p < _syst.N(); p++) {
			std::complex<number> cexp_arg(0, -_syst.positions[p].dot(_k_vectors[k]));
			k_density += std::exp(cexp_arg) * _syst.dipoles[p].dot(_k_vectors[k]);
		}

		E_k += (k_density * std::conj(k_density)).real() * _k_factors[k];
	}

	number E_tot = E_r + E_k + E_self;

	std::cout << "Real part: " << E_r << std::endl;
	std::cout << "Reciprocal part: " << E_k << std::endl;
	std::cout << "Self part: " << E_self << std::endl;
	std::cout << "Total: " << E_tot << " " << E_tot / _syst.N() << std::endl;
}
