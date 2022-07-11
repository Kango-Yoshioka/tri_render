//
// Created by kango on 2022/07/11.
//

#ifndef TRI_RENDER_TRIMESHOBJ_H
#define TRI_RENDER_TRIMESHOBJ_H

#include <vector>
#include "TriMesh.h"
#include "Eigen/Dense"

class TriMeshObj {
public:
    std::vector<TriMesh> triMeshes;
    Eigen::Vector3d color255;
    double kd;
    bool is_light;

    TriMeshObj();

    TriMeshObj(const Eigen::Vector3d &color255, const bool
    &isLight, const double &kd);

    void clearTriMeshes();

    void setRectangleMesh(const Eigen::Vector3d &center, const
    Eigen::Vector3d &width, const Eigen::Vector3d &height);

    void setOctahedron(const Eigen::Vector3d &center, const
    Eigen::Vector3d &width, const Eigen::Vector3d &height, const Eigen::Vector3d &depth);
};


#endif //TRI_RENDER_TRIMESHOBJ_H
