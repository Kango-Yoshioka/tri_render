//
// Created by kango on 2022/07/11.
//

#include "TriMeshObj.h"


TriMeshObj::TriMeshObj() {

}

TriMeshObj::TriMeshObj(const Eigen::Vector3d &color255,
                       const bool &isLight, const double &kd) : color255(
        color255), is_light(isLight), kd(kd) {}

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

    const TriMesh t1 = TriMesh(v1, v2, v4, color255, is_light, kd);
    const TriMesh t2 = TriMesh(v3, v4, v2, color255, is_light, kd);
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
void TriMeshObj::setOctahedron(const Eigen::Vector3d &center, const Eigen::Vector3d &width, const Eigen::Vector3d &height, const Eigen::Vector3d &depth) {
    clearTriMeshes();
    const Eigen::Vector3d v1 = center + height;
    const Eigen::Vector3d v2 = center + width + depth;
    const Eigen::Vector3d v3 = center + width - depth;
    const Eigen::Vector3d v4 = center - width - depth;
    const Eigen::Vector3d v5 = center - width + depth;
    const Eigen::Vector3d v6 = center - height;

    const TriMesh t1 = TriMesh(v1, v5, v2, color255, is_light, kd);
    const TriMesh t2 = TriMesh(v1, v2, v3, color255, is_light, kd);
    const TriMesh t3 = TriMesh(v1, v3, v4, color255, is_light, kd);
    const TriMesh t4 = TriMesh(v1, v4, v5, color255, is_light, kd);
    const TriMesh t5 = TriMesh(v6, v2, v5, color255, is_light, kd);
    const TriMesh t6 = TriMesh(v6, v3, v2, color255, is_light, kd);
    const TriMesh t7 = TriMesh(v6, v4, v3, color255, is_light, kd);
    const TriMesh t8 = TriMesh(v6, v5, v4, color255, is_light, kd);

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


