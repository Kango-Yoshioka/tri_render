//
// Created by kango on 2022/07/09.
//

#include "TriMesh.h"
#include "util.h"
#include <stdio.h>
#include <iostream>

TriMesh::TriMesh() {}

TriMesh::TriMesh(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
                 const Eigen::Vector3d &v3, const Eigen::Vector3d &color255,
                 bool isLight, double kd) : v1(v1), v2(v2), v3(v3),
                                            is_light(isLight), kd(kd) {
    this->n = ((v2 - v1).cross(v3 - v1)).normalized();
    this->color = rgbNormalize(color255);

    /**
     * 面積の計算
     */
    this->A = 0.5 * ((v2 - v1).cross(v3 - v1)).norm();
}

/**
 * すでに存在するTriMeshオブジェクトのパラメータを再設定する
 * @param v1
 * @param v2
 * @param v3
 * @param color255
 * @param isLight
 * @param kd
 */
void TriMesh::setTriMesh(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
                         const Eigen::Vector3d &v3,
                         const Eigen::Vector3d &color255, bool isLight,
                         double kd) {
    this->v1 = v1;
    this->v2 = v2;
    this->v3 = v3;
    this->is_light = isLight;
    this->kd = kd;

    this->n = ((v2 - v1).cross(v3 - v1)).normalized();
    this->color = rgbNormalize(color255);
}

/**
 * TriMeshに現在格納されている情報を表示する
 */
void TriMesh::printInfo() {
    std::cout << "Name:\t" << this << std::endl;
    std::cout << "v1:\t{" << this->v1.x() << "," << this->v1.y() << "," <<
              this->v1.z() << "}" << std::endl;
    std::cout << "v2:\t{" << this->v2.x() << "," << this->v2.y() << "," <<
              this->v2.z() << "}" << std::endl;
    std::cout << "v3:\t{" << this->v3.x() << "," << this->v3.y() << "," <<
              this->v3.z() << "}" << std::endl;
    std::cout << "n:\t{" << this->n.x() << "," << this->n.y() << "," <<
              this->n.z() << "}" << std::endl;
    std::cout << "isLight:\t" << this->is_light << std::endl;
    std::cout << "color255:\t{" << this->color.x() << "," << this->color.y() <<
              "," << this->color.z() << "}" << std::endl;
    std::cout << "kd:\t" << this->kd << std::endl;
    std::cout << "-------------------------------" << std::endl;
}