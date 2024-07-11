//
//  tet_quality.h
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 12/28/23.
#pragma once

#include <mtet/mtet.h>
#include <mtet/io.h>
#include <valarray>
#include <nlohmann/json.hpp>
#include <fstream>

/// Labels for grid stats.
const std::array<std::string, 6> tet_metric_labels = {"total tet number: ",
    "active tet number: ",
    "minimum radius ratio among all tets: ",
    "minimum radius ratio amond active tets: ",
    "two functions' zero-crossing checks: ",
    "three functions' zero-crossing checks: "
};

/// calculate the radius ratio of the tetrahedra based on the vertices' locations. This value is scaled by 3 such that the perfect tet has the radius ratio of 1.
///@param[in] pts           4 arrays of three-tuples represent the coordinate of tet vertices. Each subarray represents a coordinate.
///
///@return          radius ratio of this tet.
///
double tet_radius_ratio(const std::array<std::valarray<double>,4> &pts);

/// saves the tet metrics to a file
/// @param[in] tet_metric_labels            The labels of the metrics.
/// @param[in] tet_metric           The values of the metrics, corresponding to the lables.
/// @param[in] filename            The name of the output file.
///
/// @return         Whether this saving procedure is successful.
bool save_metrics(const std::string& filename,
                  const std::array<std::string, 6>& tet_metric_labels,
                  const std::valarray<double>& tet_metric);

/// returns the dot product of input arrays `a` and `b`.
double dot(const std::valarray<double> &a, const std::valarray<double> &b);

/// returns the normalized vector of the input array `a` where it represents a vector in 3D.
std::valarray<double> normalize(const std::valarray<double> &a);

/// returns the norm of the input array `a` where it represents a vector in 3D.
double norm(const std::valarray<double> &a);

/// returns a perpendicular vector for the 2D vector `a`.
std::valarray<double> perp(const std::valarray<double> &a);

/// returns a cross product for the 3D vector `a` and `b`.
std::valarray<double> cross(const std::valarray<double> &a, const std::valarray<double> &b);
