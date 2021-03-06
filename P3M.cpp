/*
 * P3M.cpp
 *
 *  Created on: 9 May 2018
 *      Author: lorenzo
 */

#include "P3M.h"

#include "System.h"

P3M::P3M(System &syst, int NM, int assignment_OP) :
				_syst(syst),
				_N_mesh_side(NM),
				_assignment_OP(assignment_OP),
				_assignment_OP_is_odd(true) {

	int partial = _N_mesh_side;
	while(partial > 1) {
		if(partial % 2) {
			std::cerr << "NM should be a power of two" << std::endl;
			exit(1);
		}
		partial /= 2;
	}

	if(_assignment_OP == 0) {
		std::cerr << "The assignment order parameter must be > 0 (P == " << _assignment_OP << ")" << std::endl;
		exit(1);
	}

	if(_assignment_OP % 2 == 0) {
		_assignment_OP_is_odd = false;
		std::cerr << "The assignment order parameter should be an odd number (P == " << _assignment_OP << ")" << std::endl;
		exit(1);
	}

	_mesh_spacing = _syst.box / _N_mesh_side;
	_mesh_cell_volume = CUB(_mesh_spacing);
	_N_mesh = CUB(_N_mesh_side);

	_fftw_plan = rfftw3d_create_plan(_N_mesh_side, _N_mesh_side, _N_mesh_side, FFTW_FORWARD, FFTW_ESTIMATE);

	for(int d = 0; d < 3; d++) {
		_dipole_density[d].resize(_N_mesh);
		// see FFTW docs
		_transformed_dipole_density[d].resize(_N_mesh_side * _N_mesh_side * (_N_mesh_side / 2 + 1));
	}

	_green_function.resize(_N_mesh);

	_E_self_per_particle = -2. * CUB(_syst.alpha) / (3. * sqrt(M_PI));
	_E_corr_per_particle = 0.;
	vec3 k;
	vec3 base_ik(2. * M_PI / _syst.box, 2. * M_PI / _syst.box, 2. * M_PI / _syst.box);
	Eigen::Vector3i idx(0, 0, 0);
	for(idx[0] = 0; idx[0] < _N_mesh_side; idx[0]++) {
		k[0] = base_ik[0];
		k[0] *= (idx[0] > _N_mesh_side / 2) ? idx[0] - _N_mesh_side : idx[0];
		for(idx[1] = 0; idx[1] < _N_mesh_side; idx[1]++) {
			k[1] = base_ik[1];
			k[1] *= (idx[1] > _N_mesh_side / 2) ? idx[1] - _N_mesh_side : idx[1];
			for(idx[2] = 0; idx[2] < _N_mesh_side; idx[2]++) {
				k[2] = base_ik[2];
				k[2] *= (idx[2] > _N_mesh_side / 2) ? idx[2] - _N_mesh_side : idx[2];

				if(!(idx[0] == 0 && idx[1] == 0 && idx[2] == 0)) {
					int r_idx;
					if(idx[2] > _N_mesh_side / 2) {
						Eigen::Vector3i new_idx(idx);
						new_idx[0] = (idx[0] > _N_mesh_side / 2) ? _N_mesh_side - idx[0] : idx[0];
						new_idx[1] = (idx[1] > _N_mesh_side / 2) ? _N_mesh_side - idx[1] : idx[1];
						new_idx[2] = _N_mesh_side - idx[2];
						r_idx = _cell_index(new_idx);
					}
					else {
						r_idx = _cell_index(idx);
					}

					cvec3 Dk = k * complex(0., 1.);
					auto factors = _G_factors(k);
					_green_function[r_idx] = factors.first / SQR(factors.second * Dk.dot(Dk).real());

					_E_corr_per_particle += -k.dot(k) * _green_function[r_idx] * factors.second;

					//printf("U_sum %d %d %d = %lf\n", idx[0], idx[1], idx[2], factors.second);
				}
			}
		}
	}

	_E_corr_per_particle /= 6. * CUB(_syst.box);

	std::cout << "Average dipolar energy: " << _E_corr_per_particle * _syst.N() * CUB(_syst.box) << std::endl;;

	_E_corr_per_particle += 2. * M_PI / (3. * CUB(_syst.box)) + _E_self_per_particle;
	// E_corr is defined as - Delta U in the Cerda paper, meaning that we need to change the
	// sign of what we have computed till now, which is essentially Delta U
	_E_corr_per_particle *= -1;
}

P3M::~P3M() {
	rfftwnd_destroy_plan(_fftw_plan);
}

number sinc(number arg) {
#define epsi 0.1

#define c2 -0.1666666666667e-0
#define c4  0.8333333333333e-2
#define c6 -0.1984126984127e-3
#define c8  0.2755731922399e-5

	if(fabs(arg) > epsi) {
		return sin(arg) / arg;
	}
	else {
		number arg_sqr = SQR(arg);
		return 1.0 + arg_sqr * (c2 + arg_sqr * (c4 + arg_sqr * (c6 + arg_sqr * c8)));
	}
}

std::pair<number, number> P3M::_G_factors(vec3 &k) {
	number denominator_factor = 0;
	number numerator_factor = 0.;
	// taken from espresso
	int limit = 5;

	cvec3 Dk = k * complex(0., 1.);
	Eigen::Vector3i idx;
	for(idx[0] = -limit; idx[0] <= limit; idx[0]++) {
		number k_m_x = k[0] + 2. * M_PI / _mesh_spacing * idx[0];
		number Uk_1 = sinc(0.5 * k_m_x * _mesh_spacing);

		for(idx[1] = -limit; idx[1] <= limit; idx[1]++) {
			number k_m_y = k[1] + 2. * M_PI / _mesh_spacing * idx[1];
			number Uk_2 = Uk_1 * sinc(0.5 * k_m_y * _mesh_spacing);

			for(idx[2] = -limit; idx[2] <= limit; idx[2]++) {
				number k_m_z = k[2] + 2. * M_PI / _mesh_spacing * idx[2];
				number Uk = Uk_2 * sinc(0.5 * k_m_z * _mesh_spacing);

//				printf("%e %e %e -- %e\n", Uk_1, Uk_2, Uk, 0.5 * k_m_x * _mesh_spacing);
//				exit(1);

				Uk = pow(Uk, _assignment_OP);
				vec3 k_m(k_m_x, k_m_y, k_m_z);
				cvec3 Dk_m = k_m * complex(0., 1.);
				number k_m_sqr = k_m.dot(k_m);
				number phik_m = 4. * M_PI / k_m_sqr * exp(-k_m_sqr / (4. * SQR(_syst.alpha)));
				numerator_factor += SQR(Dk.dot(Dk_m).real()) * SQR(Uk) * phik_m;
				denominator_factor += SQR(Uk);
			}
		}
	}

	return std::pair<number, number>(numerator_factor, denominator_factor);
}

void P3M::_assign_dipole_density() {
	// reset the densities
	for(int d = 0; d < 3; d++) {
		std::fill(_dipole_density[d].begin(), _dipole_density[d].end(), 0.);
	}

	for(int i = 0; i < _syst.N(); i++) {
		vec3 &p_pos = _syst.positions[i];
		vec3 &p_dipole = _syst.dipoles[i];

		vec3 normalised_coords = ((p_pos / _syst.box).array() - (p_pos / _syst.box).array().floor()) * _N_mesh_side;
		vec3 reference_coords;
		Eigen::Vector3i closest_mesh_point;
		if(_assignment_OP_is_odd) {
			vec3 cmp = normalised_coords.array().round();
			// we once again use PBC to map the last edges on the first ones
			closest_mesh_point = cmp.cast<int>().unaryExpr([&](const int x) { return x % _N_mesh_side; });
			reference_coords = closest_mesh_point.cast<number>() * _mesh_spacing;
		}

		switch(_assignment_OP) {
		case 1: {
			int idx = _cell_index(closest_mesh_point);
			_dipole_density[0][idx] += p_dipole[0];
			_dipole_density[1][idx] += p_dipole[1];
			_dipole_density[2][idx] += p_dipole[2];
			break;
		}
		default:
			std::cerr << "Unsupported P value " << _assignment_OP << std::endl;
			exit(1);
		}
	}

	for(int i = 0; i < _N_mesh; i++) {
		_dipole_density[0][i] /= _mesh_cell_volume;
		_dipole_density[1][i] /= _mesh_cell_volume;
		_dipole_density[2][i] /= _mesh_cell_volume;
	}
}

int P3M::_cell_index(Eigen::Vector3i v) {
	return v[2] + _N_mesh_side * (v[1] + _N_mesh_side * v[0]);
}

void P3M::print_energy() {
	number E_self = _syst.N() * _E_self_per_particle;
	number E_corr = _syst.N() * _E_corr_per_particle;

	number E_r = 0.;
	for(int p = 0; p < _syst.N(); p++) {
		for(int q = 0; q < _syst.N(); q++) {
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
	_assign_dipole_density();
	rfftwnd_one_real_to_complex(_fftw_plan, _dipole_density[0].data(), _transformed_dipole_density[0].data());
	rfftwnd_one_real_to_complex(_fftw_plan, _dipole_density[1].data(), _transformed_dipole_density[1].data());
	rfftwnd_one_real_to_complex(_fftw_plan, _dipole_density[2].data(), _transformed_dipole_density[2].data());

	cvec3 Dk;
	cvec3 base_ik = complex(0., 2. * M_PI / _syst.box) * vec3(1., 1., 1.);
	Eigen::Vector3i idx(0, 0, 0);
	for(idx[0] = 0; idx[0] < _N_mesh_side; idx[0]++) {
		Dk[0] = base_ik[0];
		Dk[0] *= (idx[0] > _N_mesh_side / 2) ? idx[0] - _N_mesh_side : idx[0];
		for(idx[1] = 0; idx[1] < _N_mesh_side; idx[1]++) {
			Dk[1] = base_ik[1];
			Dk[1] *= (idx[1] > _N_mesh_side / 2) ? idx[1] - _N_mesh_side : idx[1];
			for(idx[2] = 0; idx[2] < _N_mesh_side; idx[2]++) {
				Dk[2] = base_ik[2];
				Dk[2] *= (idx[2] > _N_mesh_side / 2) ? idx[2] - _N_mesh_side : idx[2];

				if(!(idx[0] == 0 && idx[1] == 0 && idx[2] == 0)) {
					int r_idx;
					if(idx[2] > _N_mesh_side / 2) {
						Eigen::Vector3i new_idx(idx);
						new_idx[0] = (idx[0] > _N_mesh_side / 2) ? _N_mesh_side - idx[0] : idx[0];
						new_idx[1] = (idx[1] > _N_mesh_side / 2) ? _N_mesh_side - idx[1] : idx[1];
						new_idx[2] = _N_mesh_side - idx[2];
						r_idx = _cell_index(new_idx);
					}
					else {
						r_idx = _cell_index(idx);
					}

					// TODO fftw_complex -> std::complex can be made much more straightforward
					cvec3 dipole_density_k;
					dipole_density_k[0].real(_transformed_dipole_density[0][r_idx].re);
					dipole_density_k[0].imag(_transformed_dipole_density[0][r_idx].im);
					dipole_density_k[1].real(_transformed_dipole_density[1][r_idx].re);
					dipole_density_k[1].imag(_transformed_dipole_density[1][r_idx].im);
					dipole_density_k[2].real(_transformed_dipole_density[2][r_idx].re);
					dipole_density_k[2].imag(_transformed_dipole_density[2][r_idx].im);
					dipole_density_k *= _mesh_cell_volume;

					E_k += std::norm(dipole_density_k.dot(Dk)) * _green_function[r_idx];
				}
			}
		}
	}
	E_k /= 2 * CUB(_syst.box);

	number E_tot = E_r + E_k + E_self + E_corr;

	std::cout << "Real part: " << E_r << std::endl;
	std::cout << "Reciprocal part: " << E_k << std::endl;
	std::cout << "Self part: " << E_self << std::endl;
	std::cout << "P3M correction: " << E_corr << std::endl;
	std::cout << "Total: " << E_tot << " " << E_tot / _syst.N() << std::endl;
}
