//
// Created by kango on 2022/07/09.
//

#ifndef TRI_RENDER_TRIMESH_H
#define TRI_RENDER_TRIMESH_H

#include "Eigen/Dense"

class TriMesh {

public:
    TriMesh();

    TriMesh(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
            const Eigen::Vector3d &v3, const Eigen::Vector3d &color255,
            bool isLight, double kd);

    void setTriMesh(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
                    const Eigen::Vector3d &v3, const Eigen::Vector3d &color255,
                    bool isLight, double kd);

    void printInfo();

    Eigen::Vector3d color;
    Eigen::Vector3d v1;
    Eigen::Vector3d v2;
    Eigen::Vector3d v3;
    Eigen::Vector3d n;
    bool is_light;
    double kd;
};


#endif //TRI_RENDER_TRIMESH_H
