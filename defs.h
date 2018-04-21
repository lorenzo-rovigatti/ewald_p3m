/*
 * defs.h
 *
 *  Created on: 20 Apr 2018
 *      Author: lorenzo
 */

#ifndef DEFS_H_
#define DEFS_H_

#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <iostream>

#ifdef EP_FLOAT
typedef float number;
using vec3 = Eigen::Vector3f;
using cvec3 = Eigen::Vector3cf;
using mat3 = Eigen::Matrix3f;
#else
typedef double number;
using vec3 = Eigen::Vector3d;
using cvec3 = Eigen::Vector3cd;
using mat3 = Eigen::Matrix3d;
#endif

typedef unsigned int uint;

#endif /* DEFS_H_ */
