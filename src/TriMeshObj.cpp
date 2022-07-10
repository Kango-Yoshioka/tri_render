//
// Created by kango on 2022/07/11.
//

#include "TriMeshObj.h"


TriMeshObj::TriMeshObj() {

}

TriMeshObj::TriMeshObj(const Eigen::Vector3d &color255,
                       const bool &isLight, const double &kd) : color255(
        color255), is_light(isLight), kd(kd) {}

void TriMeshObj::setRectangleMesh(const Eigen::Vector3d &center,
                                  const Eigen::Vector3d &width,
                                  const Eigen::Vector3d &height) {
    const Eigen::Vector3d v1 = center - width - height;
    const Eigen::Vector3d v2 = center + width - height;
    const Eigen::Vector3d v3 = center + width + height;
    const Eigen::Vector3d v4 = center - width + height;

    const TriMesh t1 = TriMesh(v1, v2, v4, color255, is_light, kd);
    const TriMesh t2 = TriMesh(v3, v4, v2, color255, is_light, kd);
    triMeshes.push_back(t1);
    triMeshes.push_back(t2);
}
