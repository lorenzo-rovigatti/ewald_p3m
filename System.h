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

class System {
public:
	System();
	virtual ~System();

private:
	number _box;

	std::vector<vec3> _positions;
	std::vector<vec3> _dipoles;
};

#endif /* SYSTEM_H_ */
