//
// Created by kango on 2022/07/09.
//

#include <time.h>
#include <stdlib.h>
#include "Eigen/Dense"

/**
 * 乱数のシードの設定
 */
void randInit() {
    srand((unsigned)time(NULL));
}

/**
 * 0~1のdouble型乱数を返す
 * @return 0~1のdouble型乱数
 */
double drand48() {
    return (double) (rand()) / (RAND_MAX); /* RAND_MAX = 32767 */
}

/**
 * 与えられた値を指定範囲内の値で返す 範囲外の値を範囲内の値にして返す
 * @param val
 * @param min
 * @param max
 * @return
 */
double dClamp(const double &val, const double &min, const double &max) {
    double out_val = val;
    if (val < min) {
        out_val = min;
    } else if (val > max) {
        out_val = max;
    }

    return out_val;
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