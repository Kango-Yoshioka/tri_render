//
// Created by kango on 2022/07/09.
//

#include <time.h>
#include <stdlib.h>
#include "Eigen/Dense"

void randInit() {
    srand((unsigned)time(NULL));
}

double drand48() {
    return (double) (rand()) / (RAND_MAX); /* RAND_MAX = 32767 */
}

/**
 * @param rgb 0~255
 * @param out_rgb 0~1
 */
Eigen::Vector3d rgbNormalize(const Eigen::Vector3d rgb) {
    Eigen::Vector3d out_rgb{
            rgb.x() / 255,
            rgb.y() / 255,
            rgb.z() / 255
    };

    return out_rgb;
}