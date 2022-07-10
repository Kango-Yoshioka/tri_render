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

#endif //TRI_RENDER_UTIL_H
