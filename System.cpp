/*
 * System.cpp
 *
 *  Created on: 20 Apr 2018
 *      Author: lorenzo
 */

#include "System.h"

#include <fstream>

System::System(number mybox, int N) :
				box(mybox),
				alpha(2.0) {
	positions.resize(N);
	dipoles.resize(N);
	for(int i = 0; i < N; i++) {
		positions[i] = 0.5 * box * (vec3::Random() + vec3(1., 1., 1.));
		dipoles[i] = _random_vector_on_unit_sphere();
	}
}

System::~System() {

}

uint System::N() {
	return positions.size();
}

void System::print_conf(std::string filename) {
	std::ofstream conf(filename.c_str());

	conf << 0 << " " << 0 << " " << N() << " " << 0 << " " << 21 << std::endl;
	conf << box << " " << box << " " << box << " " << 0. << " " << 0. << " " << 0.2 << std::endl;

	for(int p = 0; p < N(); p++) {
		auto &dip = dipoles[p];
		auto &pos = positions[p];

		conf << dip[0] << " " << dip[1] << " " << dip[2] << std::endl;
		conf << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	}

	conf.close();
}

vec3 System::_random_vector_on_unit_sphere() {
	number ransq = 1.;
	number ran1, ran2;

	while(ransq >= 1) {
		ran1 = 1. - 2. * drand48();
		ran2 = 1. - 2. * drand48();
		ransq = ran1 * ran1 + ran2 * ran2;
	}

	number ranh = 2. * sqrt(1. - ransq);
	return vec3(ran1 * ranh, ran2 * ranh, 1. - 2. * ransq);
}
