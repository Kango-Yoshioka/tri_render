//
// Created by kango on 2022/07/09.
//

#ifndef TRI_RENDER_UTIL_H
#define TRI_RENDER_UTIL_H

#include <Eigen/Core>

void randInit();

double drand48();

double dClamp(const double &val, const double &min, const double &max);

Eigen::Vector3d rgbNormalize(const Eigen::Vector3d rgb);

Eigen::Vector3d coordinateTransformation(const Eigen::Vector3d &x,
                                         const Eigen::Vector3d &p,
                                         const Eigen::Vector3d &q,
                                         const Eigen::Vector3d &r);

Eigen::Vector3d plainToSolidAngle(const double &theta, const double &phi);

#endif //TRI_RENDER_UTIL_H
