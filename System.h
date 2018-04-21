/*
 * System.h
 *
 *  Created on: 20 Apr 2018
 *      Author: lorenzo
 */

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "defs.h"

#include <vector>

struct System {
public:
	System(number my_box, int N);
	virtual ~System();

	uint N();
	void print_conf(std::string filename);

	std::vector<vec3> positions;
	std::vector<vec3> dipoles;

	number box;
	number alpha;

private:
	vec3 _random_vector_on_unit_sphere();
};

#endif /* SYSTEM_H_ */
