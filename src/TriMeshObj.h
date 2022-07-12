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
    TriMeshType triMeshType;
    double intensity;

    TriMeshObj(const Eigen::Vector3d &color255, const TriMeshType &triMeshType,
               const double &kd, const double &intensity = DEFAULT_INTENSITY);

    void clearTriMeshes();

    void setRectangleMesh(const Eigen::Vector3d &center, const
    Eigen::Vector3d &width, const Eigen::Vector3d &height);

    void setOctahedron(const Eigen::Vector3d &center, const
    Eigen::Vector3d &width, const Eigen::Vector3d &height,
                       const Eigen::Vector3d &depth);

    void setCone(const Eigen::Vector3d &center, const double &radius, const
    Eigen::Vector3d &height, const int &division_num);

    void setSphere(const Eigen::Vector3d &center, const double &radius,
                   const int &latitude_division_num);
};


#endif //TRI_RENDER_TRIMESHOBJ_H
