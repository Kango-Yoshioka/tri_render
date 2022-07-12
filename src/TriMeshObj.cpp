//
// Created by kango on 2022/07/11.
//

#include <iostream>
#include "TriMeshObj.h"

#include "util.h"

/**
 * TriMeshObjのコンストラクタ
 * TriMeshをまとめたもの
 * @param color255
 * @param triMeshType (TriMeshObjのタイプ。LIGHT,DIFFUSE)
 * @param kd
 * @param intensity (デフォルト引数=DEFAULT_INTENSITY)
 */
TriMeshObj::TriMeshObj(const Eigen::Vector3d &color255,
                       const TriMeshType &triMeshType, const double &kd,
                       const double &intensity)
        : color255(color255), triMeshType(triMeshType), kd(kd),
          intensity(intensity) {}

/**
 * 四角形メッシュを作成
 * @param center
 * @param width
 * @param height
 */
void TriMeshObj::setRectangleMesh(const Eigen::Vector3d &center,
                                  const Eigen::Vector3d &width,
                                  const Eigen::Vector3d &height) {
    clearTriMeshes();
    const Eigen::Vector3d v1 = center - width - height;
    const Eigen::Vector3d v2 = center + width - height;
    const Eigen::Vector3d v3 = center + width + height;
    const Eigen::Vector3d v4 = center - width + height;

    const TriMesh t1 = TriMesh(v1, v2, v4, color255, triMeshType, kd,
                               intensity);
    const TriMesh t2 = TriMesh(v3, v4, v2, color255, triMeshType, kd,
                               intensity);
    triMeshes.push_back(t1);
    triMeshes.push_back(t2);
}

/**
 * 八面体メッシュを作成
 * @param center
 * @param width
 * @param height
 * @param depth
 */
void TriMeshObj::setOctahedron(const Eigen::Vector3d &center,
                               const Eigen::Vector3d &width,
                               const Eigen::Vector3d &height,
                               const Eigen::Vector3d &depth) {
    clearTriMeshes();
    const Eigen::Vector3d v1 = center + height;
    const Eigen::Vector3d v2 = center + width + depth;
    const Eigen::Vector3d v3 = center + width - depth;
    const Eigen::Vector3d v4 = center - width - depth;
    const Eigen::Vector3d v5 = center - width + depth;
    const Eigen::Vector3d v6 = center - height;

    const TriMesh t1 = TriMesh(v1, v5, v2, color255, triMeshType, kd,
                               intensity);
    const TriMesh t2 = TriMesh(v1, v2, v3, color255, triMeshType, kd,
                               intensity);
    const TriMesh t3 = TriMesh(v1, v3, v4, color255, triMeshType, kd,
                               intensity);
    const TriMesh t4 = TriMesh(v1, v4, v5, color255, triMeshType, kd,
                               intensity);
    const TriMesh t5 = TriMesh(v6, v2, v5, color255, triMeshType, kd,
                               intensity);
    const TriMesh t6 = TriMesh(v6, v3, v2, color255, triMeshType, kd,
                               intensity);
    const TriMesh t7 = TriMesh(v6, v4, v3, color255, triMeshType, kd,
                               intensity);
    const TriMesh t8 = TriMesh(v6, v5, v4, color255, triMeshType, kd,
                               intensity);

    triMeshes.push_back(t1);
    triMeshes.push_back(t2);
    triMeshes.push_back(t3);
    triMeshes.push_back(t4);
    triMeshes.push_back(t5);
    triMeshes.push_back(t6);
    triMeshes.push_back(t7);
    triMeshes.push_back(t8);
}

/**
 * triMeshesに保存されているtriMeshを削除
 */
void TriMeshObj::clearTriMeshes() {
    triMeshes.clear();
}

/**
 * @param center 円錐底面の中心座標
 * @param radius 円錐底面の半径
 * @param height 円錐の向き、高さを示すベクトル
 */
void TriMeshObj::setCone(const Eigen::Vector3d &center, const double &radius,
                         const Eigen::Vector3d &height,
                         const int &division_num) {
    clearTriMeshes();

    const Eigen::Vector3d v0 = center + height;
    Eigen::Vector3d arm;
    if (height.x() == 0) {
        arm = Eigen::Vector3d{1.0, 0.0, 0.0};
    } else {
        arm = Eigen::Vector3d{(-height.y() - height.z()) / height.x(),
                              1.0, 1.0}.normalized();
    }
    const Eigen::Vector3d p = arm.normalized();
    const Eigen::Vector3d q = height.normalized();
    const Eigen::Vector3d r = (p.cross(q)).normalized();

    for (int i = 0; i < division_num; i++) {
        const double theta = 2 * M_PI / division_num * i;
        const double next_theta = 2 * M_PI / division_num * (i + 1);

        const Eigen::Vector3d color{
                color255.x() * sin(0.5 * theta),
                color255.y() * (1 - sin(0.5 * theta)),
                color255.z() * 0.5,
        };

        Eigen::Vector3d arm_theta{
                cos(theta),
                0.0,
                sin(theta)
        };

        Eigen::Vector3d arm_next_theta{
                cos(next_theta),
                0.0,
                sin(next_theta)
        };

        arm_theta = coordinateTransformation(arm_theta, p, q, r);
        arm_next_theta = coordinateTransformation(arm_next_theta, p, q, r);
        /**
         * 側面のメッシュ
         */
        const TriMesh t1 = TriMesh(v0, center + radius * arm_next_theta,
                                   center + radius * arm_theta,
                                   color, triMeshType, kd, intensity);
        /**
         * 底面のメッシュ
         */
        const TriMesh t2 = TriMesh(center, center + radius * arm_theta,
                                   center + radius * arm_next_theta,
                                   color, triMeshType, kd, intensity);
        triMeshes.push_back(t1);
        triMeshes.push_back(t2);
    }
}

/**
 * 球メッシュを作成
 * @param center
 * @param radius
 * @param latitude_division_num 緯度方向の分割数。経度方向の分割数はその倍になる。
 */
void TriMeshObj::setSphere(const Eigen::Vector3d &center, const double &radius,
                           const int &latitude_division_num) {
    clearTriMeshes();
    /**
     * 0 < theta < pi
     * 0 < phi < 2 * pi
     */
    const int longitude_division_num = 2 * latitude_division_num;
    for (int i = 0; i < latitude_division_num; i++) {
        const double theta = M_PI / latitude_division_num * i;
        const double next_theta = M_PI / latitude_division_num * (i + 1);

        for (int j = 0; j < longitude_division_num; j++) {
            const double phi = 2 * M_PI / longitude_division_num * j;
            const double next_phi =
                    2 * M_PI / longitude_division_num * (j + 1);

            const Eigen::Vector3d v1 = center +
                                       radius *
                                       plainToSolidAngle(theta, next_phi);
            const Eigen::Vector3d v2 = center +
                                       radius *
                                       plainToSolidAngle(next_theta, next_phi);
            const Eigen::Vector3d v3 = center +
                                       radius *
                                       plainToSolidAngle(next_theta, phi);
            const Eigen::Vector3d v4 =
                    center + radius * plainToSolidAngle(theta, phi);

            const TriMesh t1 = TriMesh(v1, v2, v4, color255, triMeshType, kd,
                                       intensity);
            const TriMesh t2 = TriMesh(v2, v3, v4, color255, triMeshType, kd,
                                       intensity);

            /**
             * もしi = 0ならt2のみ追加、i = 2 * latitude_division_num - 1ならt1のみ追加
             */
            if (i == 0) {
                triMeshes.push_back(t2);
            } else if (i == latitude_division_num - 1) {
                triMeshes.push_back(t1);
            } else {
                triMeshes.push_back(t1);
                triMeshes.push_back(t2);
            }
        }
    }
}



