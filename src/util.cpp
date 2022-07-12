//
// Created by kango on 2022/07/09.
//

#include <time.h>
#include <stdlib.h>
#include "Eigen/Dense"
#include "util.h"


/**
 * 乱数のシードの設定
 */
void randInit() {
    srand((unsigned) time(NULL));
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

/**
 * 任意の座標軸におけるxベクトルをデカルト座標軸上で表現したベクトルを返す
 * @param x 任意の座標軸におけるベクトル
 * @param p 任意のx軸
 * @param q 任意のy軸
 * @param r 任意のz軸
 * @return デカルト座標系で表現したベクトル
 */
Eigen::Vector3d coordinateTransformation(const Eigen::Vector3d &x,
                                         const Eigen::Vector3d &p,
                                         const Eigen::Vector3d &q,
                                         const Eigen::Vector3d &r) {
    const Eigen::Vector3d pn = p.normalized();
    const Eigen::Vector3d qn = q.normalized();
    const Eigen::Vector3d rn = r.normalized();

    return Eigen::Vector3d{
            x.x() * pn.x() + x.y() * qn.x() + x.z() * rn.x(),
            x.x() * pn.y() + x.y() * qn.y() + x.z() * rn.y(),
            x.x() * pn.z() + x.y() * qn.z() + x.z() * rn.z()
    };
}

/**
 * y軸からの回転をtheta,x軸からz軸への回転をphiとするときの
 * 立体角omegaの方向ベクトルを返す
 * @param theta
 * @param phi
 * @return
 */
Eigen::Vector3d plainToSolidAngle(const double &theta, const double &phi) {
    return Eigen::Vector3d{
            sin(theta) * cos(phi),
            cos(theta),
            sin(theta) * sin(phi)
    }.normalized();
}
