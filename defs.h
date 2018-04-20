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

#ifdef EP_FLOAT
typedef float number;
using vec3 = Eigen::Vector3f;
using mat3 = Eigen::Matrix3f;
#else
typedef double number;
using vec3 = Eigen::Vector3d;
using mat3 = Eigen::Matrix3d;
#endif

#endif /* DEFS_H_ */
