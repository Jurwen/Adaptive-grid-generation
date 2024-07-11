//
//  tet_quality.cpp
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 6/20/24.
//

#include "tet_quality.h"

double dot(const std::valarray<double> &a, const std::valarray<double> &b) {
    return (a * b).sum();
}

std::valarray<double> normalize(const std::valarray<double> &a) {
    return a / sqrt(dot(a, a));
}

double norm(const std::valarray<double> &a) {
    return sqrt(dot(a, a));
}

std::valarray<double> perp(const std::valarray<double> &a){
    return {-a[1], a[0]};
}

std::valarray<double> cross(const std::valarray<double> &a, const std::valarray<double> &b) {
    std::valarray<double> c(3);
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

double tet_radius_ratio(const std::array<std::valarray<double>,4> &pts)
{
    
    // Determine side vectors
    std::array<std::valarray<double>, 6> side;
    for (int i = 0; i < 6; ++i)
    {
        side[i].resize(3);
    }
    side[0] = {pts[1][0] - pts[0][0], pts[1][1] - pts[0][1],
        pts[1][2] - pts[0][2]};
    
    side[1]={pts[2][0] - pts[1][0], pts[2][1] - pts[1][1],
        pts[2][2] - pts[1][2]};
    
    side[2] = {pts[0][0] - pts[2][0], pts[0][1] - pts[2][1],
        pts[0][2] - pts[2][2]};
    
    side[3] = {pts[3][0] - pts[0][0], pts[3][1] - pts[0][1],
        pts[3][2] - pts[0][2]};
    
    side[4] = {pts[3][0] - pts[1][0], pts[3][1] - pts[1][1],
        pts[3][2] - pts[1][2]};
    
    side[5] = {pts[3][0] - pts[2][0], pts[3][1] - pts[2][1],
        pts[3][2] - pts[2][2]};
    
    std::valarray<double> numerator = dot(side[3], side[3]) * cross(side[2],side[0]) +
    dot(side[2], side[2]) * cross(side[3] ,side[0]) + dot(side[0], side[0]) * cross(side[3] , side[2]);
    
    double area_sum;
    area_sum = (norm(cross(side[2], side[0]))
                + norm(cross(side[3], side[0]))
                + norm(cross(side[4], side[1]))
                + norm(cross(side[3], side[2]))) *
    0.5;
    
    double volume = dot(pts[0] - pts[3], cross(pts[1] - pts[3], pts[2] - pts[3]))/6;
    const double radius_ratio = (108 * volume * volume)/(norm(numerator) * area_sum)  ;
    return radius_ratio;
}

bool save_metrics(const std::string& filename,
                  const std::array<std::string, 6>& tet_metric_labels,
                  const std::valarray<double>& tet_metric)
{
    // assert stats_labels.size() == stats.size()
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    //fout.open(filename.c_str(),std::ios::app);
    json jOut;
    for (size_t i = 0; i < tet_metric.size(); ++i) {
        jOut[tet_metric_labels[i]] = tet_metric[i];
    }
    fout << jOut << std::endl;
    fout.close();
    return true;
}

