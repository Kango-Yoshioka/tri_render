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
                 const TriMeshType &triMeshType, const double &kd, const
                 double &intensity) :
        v1(v1), v2(v2), v3(v3), triMeshType(triMeshType), kd(kd),
        intensity(intensity) {
    this->n = ((v2 - v1).cross(v3 - v1)).normalized();
    this->color = rgbNormalize(color255);

    /**
     * 面積の計算
     */
    this->A = 0.5 * ((v2 - v1).cross(v3 - v1)).norm();
}

/**
 * TriMeshに現在格納されている情報を表示する
 */
void TriMesh::printInfo() {
    std::string triMeshTypeStr = "UNDEFINED";
    switch (triMeshType) {
        case LIGHT:
            triMeshTypeStr = "LIGHT";
            break;
        case DIFFUSE:
            triMeshTypeStr = "DIFFUSE";
            break;
    }
    std::cout << "Name:\t" << this << std::endl;
    std::cout << "v1:\t{" << v1.x() << "," << v1.y() << "," << v1.z() << "}"
              << std::endl;
    std::cout << "v2:\t{" << v2.x() << "," << v2.y() << "," << v2.z() << "}"
              << std::endl;
    std::cout << "v3:\t{" << v3.x() << "," << v3.y() << "," <<
              v3.z() << "}" << std::endl;
    std::cout << "n:\t{" << n.x() << "," << n.y() << "," << n.z() << "}" <<
              std::endl;
    std::cout << "TriMeshType:\t" << triMeshTypeStr << std::endl;
    std::cout << "color:\t{" << color.x() << "," << color.y() << "," <<
              color.z() << "}" << std::endl;
    std::cout << "kd:\t" << kd << std::endl;
    std::cout << "intensity:\t" << intensity << std::endl;
    std::cout << "-------------------------------" << std::endl;
}
