//
// Created by kango on 2022/07/09.
//

#ifndef TRI_RENDER_TRIMESH_H
#define TRI_RENDER_TRIMESH_H

#include "Eigen/Dense"

constexpr double DEFAULT_INTENSITY = 100;

enum TriMeshType {
    LIGHT,
    DIFFUSE
};

class TriMesh {

public:
    TriMesh();

    TriMesh(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
            const Eigen::Vector3d &v3, const Eigen::Vector3d &color255,
            const TriMeshType &triMeshType, const double &kd,
            const double &intensity = DEFAULT_INTENSITY);

    void printInfo();

    Eigen::Vector3d color;
    Eigen::Vector3d v1;
    Eigen::Vector3d v2;
    Eigen::Vector3d v3;
    Eigen::Vector3d n;
    TriMeshType triMeshType;
    double kd;
    double A;
    double intensity;
};


#endif //TRI_RENDER_TRIMESH_H
