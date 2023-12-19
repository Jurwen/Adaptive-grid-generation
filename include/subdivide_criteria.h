#pragma once
#include <valarray>
#include <array>
#include <limits>

// =================== compute error of a tetrahedron ===========================

int get_sign(double x) {
    return (x > 0) ? 1 : -1;
}

double dot(const std::valarray<double> &a, const std::valarray<double> &b) {
    return (a * b).sum();
}

std::valarray<double> normalize(const std::valarray<double> &a) {
    return a / sqrt(dot(a, a));
}

double norm(const std::valarray<double> &a) {
    return sqrt(dot(a, a));
}

std::valarray<double> cross(const std::valarray<double> &a, const std::valarray<double> &b) {
    std::valarray<double> c(3);
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

double compute_tet_error(const std::array<std::valarray<double>,4> &pts,
                         const std::array<double,4> &vals,
                         const std::array<std::valarray<double>,4> &grads) {
    std::valarray<double> p0 = pts[0], p1 = pts[1], p2 = pts[2], p3 = pts[3];
    double v0 = vals[0], v1 = vals[1], v2 = vals[2], v3 = vals[3];
    const std::valarray<double> &g0 = grads[0], &g1 = grads[1], &g2 = grads[2], &g3 = grads[3];
    // Bezier control points
    std::valarray<double> v0s = {v0 + dot(g0, p1 - p0) / 3, v0 + dot(g0, p2 - p0) / 3, v0 + dot(g0, p3 - p0) / 3};
    std::valarray<double> v1s = {v1 + dot(g1, p2 - p1) / 3, v1 + dot(g1, p3 - p1) / 3, v1 + dot(g1, p0 - p1) / 3};
    std::valarray<double> v2s = {v2 + dot(g2, p3 - p2) / 3, v2 + dot(g2, p0 - p2) / 3, v2 + dot(g2, p1 - p2) / 3};
    std::valarray<double> v3s = {v3 + dot(g3, p0 - p3) / 3, v3 + dot(g3, p1 - p3) / 3, v3 + dot(g3, p2 - p3) / 3};
    double e0 = (v1s[0] + v1s[1] + v2s[0] + v2s[2] + v3s[1] + v3s[2]) / 6;
    double vMid0 = e0 + (e0 - (v1 + v2 + v3) / 3) / 2;
    double e1 = (v0s[1] + v0s[2] + v2s[0] + v2s[1] + v3s[0] + v3s[2]) / 6;
    double vMid1 = e1 + (e1 - (v0 + v2 + v3) / 3) / 2;
    double e2 = (v0s[0] + v0s[2] + v1s[1] + v1s[2] + v3s[0] + v3s[1]) / 6;
    double vMid2 = e2 + (e2 - (v0 + v1 + v3) / 3) / 2;
    double e3 = (v0s[0] + v0s[1] + v1s[0] + v1s[2] + v2s[1] + v2s[2]) / 6;
    double vMid3 = e3 + (e3 - (v0 + v1 + v2) / 3) / 2;
    std::valarray<double> all_Bezier_vals = {v0, v1, v2, v3, v0s[0], v0s[1], v0s[2], v1s[0], v1s[1], v1s[2], v2s[0], v2s[1], v2s[2],
                                        v3s[0], v3s[1], v3s[2], vMid0, vMid1, vMid2, vMid3};
    // return 0 if all values have the same sign
    if (get_sign(all_Bezier_vals.max()) == get_sign(all_Bezier_vals.min())) {
        return 0;
    }
    // compute level-set distance error
    // linear interpolation of v0, v1, v2, v3 in tet: w = a*x + b*y + c*z + d
    double a = (-(p2[2]*p3[1]*v0) + p2[1]*p3[2]*v0 + p0[2]*p2[1]*v1 - p0[1]*p2[2]*v1 -
                p0[2]*p3[1]*v1 + p2[2]*p3[1]*v1 + p0[1]*p3[2]*v1 - p2[1]*p3[2]*v1 +
                p0[2]*p3[1]*v2 - p0[1]*p3[2]*v2 - p0[2]*p2[1]*v3 + p0[1]*p2[2]*v3 +
                p1[1]*(p2[2]*v0 - p3[2]*v0 - p0[2]*v2 + p3[2]*v2 + p0[2]*v3 - p2[2]*v3) +
                p1[2]*(p3[1]*v0 + p0[1]*v2 - p3[1]*v2 - p0[1]*v3 + p2[1]*(-v0 + v3)))/
               (-(p0[0]*p1[2]*p2[1]) + p0[0]*p1[1]*p2[2] + p1[2]*p2[1]*p3[0] - p1[1]*p2[2]*p3[0] +
                p0[0]*p1[2]*p3[1] - p1[2]*p2[0]*p3[1] - p0[0]*p2[2]*p3[1] + p1[0]*p2[2]*p3[1] +
                p0[2]*(p1[0]*p2[1] - p2[1]*p3[0] + p1[1]*(-p2[0] + p3[0]) - p1[0]*p3[1] +
                       p2[0]*p3[1]) - p0[0]*p1[1]*p3[2] + p1[1]*p2[0]*p3[2] + p0[0]*p2[1]*p3[2] -
                p1[0]*p2[1]*p3[2] + p0[1]*(p1[2]*p2[0] - p1[0]*p2[2] - p1[2]*p3[0] + p2[2]*p3[0] +
                                           p1[0]*p3[2] - p2[0]*p3[2]));
    double b = (p2[2]*p3[0]*v0 - p2[0]*p3[2]*v0 - p0[2]*p2[0]*v1 + p0[0]*p2[2]*v1 + p0[2]*p3[0]*v1 -
                p2[2]*p3[0]*v1 - p0[0]*p3[2]*v1 + p2[0]*p3[2]*v1 - p0[2]*p3[0]*v2 +
                p0[0]*p3[2]*v2 + p0[2]*p2[0]*v3 - p0[0]*p2[2]*v3 +
                p1[2]*(-(p3[0]*v0) - p0[0]*v2 + p3[0]*v2 + p2[0]*(v0 - v3) + p0[0]*v3) +
                p1[0]*(-(p2[2]*v0) + p3[2]*v0 + p0[2]*v2 - p3[2]*v2 - p0[2]*v3 + p2[2]*v3))
               /(-(p0[0]*p1[2]*p2[1]) + p0[0]*p1[1]*p2[2] + p1[2]*p2[1]*p3[0] - p1[1]*p2[2]*p3[0] +
                 p0[0]*p1[2]*p3[1] - p1[2]*p2[0]*p3[1] - p0[0]*p2[2]*p3[1] + p1[0]*p2[2]*p3[1] +
                 p0[2]*(p1[0]*p2[1] - p2[1]*p3[0] + p1[1]*(-p2[0] + p3[0]) - p1[0]*p3[1] +
                        p2[0]*p3[1]) - p0[0]*p1[1]*p3[2] + p1[1]*p2[0]*p3[2] + p0[0]*p2[1]*p3[2] -
                 p1[0]*p2[1]*p3[2] + p0[1]*(p1[2]*p2[0] - p1[0]*p2[2] - p1[2]*p3[0] + p2[2]*p3[0] +
                                            p1[0]*p3[2] - p2[0]*p3[2]));
    double c = (-(p2[1]*p3[0]*v0) + p2[0]*p3[1]*v0 + p0[1]*p2[0]*v1 - p0[0]*p2[1]*v1 -
                p0[1]*p3[0]*v1 + p2[1]*p3[0]*v1 + p0[0]*p3[1]*v1 - p2[0]*p3[1]*v1 +
                p0[1]*p3[0]*v2 - p0[0]*p3[1]*v2 - p0[1]*p2[0]*v3 + p0[0]*p2[1]*v3 +
                p1[0]*(p2[1]*v0 - p3[1]*v0 - p0[1]*v2 + p3[1]*v2 + p0[1]*v3 - p2[1]*v3) +
                p1[1]*(p3[0]*v0 + p0[0]*v2 - p3[0]*v2 - p0[0]*v3 + p2[0]*(-v0 + v3)))/
               (-(p0[0]*p1[2]*p2[1]) + p0[0]*p1[1]*p2[2] + p1[2]*p2[1]*p3[0] - p1[1]*p2[2]*p3[0] +
                p0[0]*p1[2]*p3[1] - p1[2]*p2[0]*p3[1] - p0[0]*p2[2]*p3[1] + p1[0]*p2[2]*p3[1] +
                p0[2]*(p1[0]*p2[1] - p2[1]*p3[0] + p1[1]*(-p2[0] + p3[0]) - p1[0]*p3[1] +
                       p2[0]*p3[1]) - p0[0]*p1[1]*p3[2] + p1[1]*p2[0]*p3[2] + p0[0]*p2[1]*p3[2] -
                p1[0]*p2[1]*p3[2] + p0[1]*(p1[2]*p2[0] - p1[0]*p2[2] - p1[2]*p3[0] + p2[2]*p3[0] +
                                           p1[0]*p3[2] - p2[0]*p3[2]));
    // slope in tetrahedron
    double tan = sqrt(a*a + b*b + c*c);
    if (tan == 0) {
        return std::numeric_limits<double>::infinity();
    }
    std::valarray<double> normal = {a, b, c};
    normal /= tan;
    // slope on triangle 0
    double cos = dot(normalize(cross(p2-p1, p3-p1)), normal);
    double sin = abs(cos) >= 1 ? 0 : sqrt(1 - cos*cos);
    double tan0 = tan * sin;
    // slope on triangle 1
    cos = dot(normalize(cross(p2-p0, p3-p0)), normal);
    sin = abs(cos) >= 1 ? 0 : sqrt(1 - cos*cos);
    double tan1 = tan * sin;
    // slope on triangle 2
    cos = dot(normalize(cross(p1-p0, p3-p0)), normal);
    sin = abs(cos) >= 1 ? 0 : sqrt(1 - cos*cos);
    double tan2 = tan * sin;
    // slope on triangle 3
    cos = dot(normalize(cross(p1-p0, p2-p0)), normal);
    sin = abs(cos) >= 1 ? 0 : sqrt(1 - cos*cos);
    double tan3 = tan * sin;
    // slope along edges
    double tan01 = abs(v0 - v1) / norm(p0 - p1);
    double tan02 = abs(v0 - v2) / norm(p0 - p2);
    double tan03 = abs(v0 - v3) / norm(p0 - p3);
    double tan12 = abs(v1 - v2) / norm(p1 - p2);
    double tan13 = abs(v1 - v3) / norm(p1 - p3);
    double tan23 = abs(v2 - v3) / norm(p2 - p3);
    // if any slope is 0, return infinity
    if (tan0 == 0 || tan1 == 0 || tan2 == 0 || tan3 == 0 ||
        tan01 == 0 || tan02 == 0 || tan03 == 0 ||
        tan12 == 0 || tan13 == 0 || tan23 == 0) {
        return std::numeric_limits<double>::infinity();
    }
    // compute difference between Bezier values and linear values
    std::valarray<double> Bezier_vals = {v0s[0], v0s[1], v0s[2], v1s[0], v1s[1], v1s[2],
                                    v2s[0], v2s[1], v2s[2], v3s[0], v3s[1], v3s[2],
                                    vMid0, vMid1, vMid2, vMid3};
    std::vector<std::array<double, 4>> coeff = {
            {2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1}, {0, 2, 1, 0},
            {0, 2, 0, 1}, {1, 2, 0, 0}, {0, 0, 2, 1}, {1, 0, 2, 0},
            {0, 1, 2, 0}, {1, 0, 0, 2}, {0, 1, 0, 2}, {0, 0, 1, 2},
            {0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}
    };
    std::valarray<double> linear_vals(coeff.size());
    for (size_t i = 0; i < coeff.size(); ++i) {
        linear_vals[i] = v0 * coeff[i][0] + v1 * coeff[i][1] + v2 * coeff[i][2] + v3 * coeff[i][3];
    }
    linear_vals /= 3.0;
    double diff = abs(Bezier_vals - linear_vals).max();
    // distance error in tetrahedron
    double distance_error = diff / tan;
    // distance error on triangle 0
    std::valarray<double> Bezier_vals0 = {v1, v2, v3, Bezier_vals[3], Bezier_vals[4], Bezier_vals[6],
                                     Bezier_vals[8], Bezier_vals[10], Bezier_vals[11], Bezier_vals[12]};
    std::valarray<double> linear_vals0 = {v1, v2, v3, linear_vals[3], linear_vals[4], linear_vals[6],
                                     linear_vals[8], linear_vals[10], linear_vals[11], linear_vals[12]};
    if (get_sign(Bezier_vals0.max()) != get_sign(Bezier_vals0.min())) {
        double diff0 = abs(Bezier_vals0 - linear_vals0).max();
        distance_error = std::max(distance_error, diff0 / tan0);
    }
    // distance error on triangle 1
    std::valarray<double> Bezier_vals1 = {v0, v2, v3, Bezier_vals[1], Bezier_vals[2], Bezier_vals[6],
                                     Bezier_vals[7], Bezier_vals[9], Bezier_vals[11], Bezier_vals[13]};
    std::valarray<double> linear_vals1 = {v0, v2, v3, linear_vals[1], linear_vals[2], linear_vals[6],
                                     linear_vals[7], linear_vals[9], linear_vals[11], linear_vals[13]};
    if (get_sign(Bezier_vals1.max()) != get_sign(Bezier_vals1.min())) {
        double diff1 = abs(Bezier_vals1 - linear_vals1).max();
        distance_error = std::max(distance_error, diff1 / tan1);
    }
    // distance error on triangle 2
    std::valarray<double> Bezier_vals2 = {v0, v1, v3, Bezier_vals[0], Bezier_vals[2], Bezier_vals[4], Bezier_vals[5],
                                     Bezier_vals[9], Bezier_vals[10], Bezier_vals[14]};
    std::valarray<double> linear_vals2 = {v0, v1, v3, linear_vals[0], linear_vals[2], linear_vals[4], linear_vals[5],
                                     linear_vals[9], linear_vals[10], linear_vals[14]};
    if (get_sign(Bezier_vals2.max()) != get_sign(Bezier_vals2.min())) {
        double diff2 = abs(Bezier_vals2 - linear_vals2).max();
        distance_error = std::max(distance_error, diff2 / tan2);
    }
    // distance error on triangle 3
    std::valarray<double> Bezier_vals3 = {v0, v1, v2, Bezier_vals[0], Bezier_vals[1], Bezier_vals[3], Bezier_vals[5], \
Bezier_vals[7], Bezier_vals[8], Bezier_vals[15]};
    std::valarray<double> linear_vals3 = {v0, v1, v2, linear_vals[0], linear_vals[1], linear_vals[3], linear_vals[5], \
linear_vals[7], linear_vals[8], linear_vals[15]};
    if (get_sign(Bezier_vals3.max()) != get_sign(Bezier_vals3.min())) {
        double diff3 = abs(Bezier_vals3 - linear_vals3).max();
        distance_error = std::max(distance_error, diff3 / tan3);
    }
    // distance error on edge 01
    std::valarray<double> Bezier_vals01 = {v0, v1, Bezier_vals[0], Bezier_vals[5]};
    std::valarray<double> linear_vals01 = {v0, v1, linear_vals[0], linear_vals[5]};
    if (get_sign(Bezier_vals01.max()) != get_sign(Bezier_vals01.min())) {
        double diff01 = abs(Bezier_vals01 - linear_vals01).max();
        distance_error = std::max(distance_error, diff01 / tan01);
    }
    // distance error on edge 02
    std::valarray<double> Bezier_vals02 = {v0, v2, Bezier_vals[1], Bezier_vals[7]};
    std::valarray<double> linear_vals02 = {v0, v2, linear_vals[1], linear_vals[7]};
    if (get_sign(Bezier_vals02.max()) != get_sign(Bezier_vals02.min())) {
        double diff02 = abs(Bezier_vals02 - linear_vals02).max();
        distance_error = std::max(distance_error, diff02 / tan02);
    }
    // distance error on edge 03
    std::valarray<double> Bezier_vals03 = {v0, v3, Bezier_vals[2], Bezier_vals[9]};
    std::valarray<double> linear_vals03 = {v0, v3, linear_vals[2], linear_vals[9]};
    if (get_sign(Bezier_vals03.max()) != get_sign(Bezier_vals03.min())) {
        double diff03 = abs(Bezier_vals03 - linear_vals03).max();
        distance_error = std::max(distance_error, diff03 / tan03);
    }
    // distance error on edge 12
    std::valarray<double> Bezier_vals12 = {v1, v2, Bezier_vals[3], Bezier_vals[8]};
    std::valarray<double> linear_vals12 = {v1, v2, linear_vals[3], linear_vals[8]};
    if (get_sign(Bezier_vals12.max()) != get_sign(Bezier_vals12.min())) {
        double diff12 = abs(Bezier_vals12 - linear_vals12).max();
        distance_error = std::max(distance_error, diff12 / tan12);
    }
    // distance error on edge 13
    std::valarray<double> Bezier_vals13 = {v1, v3, Bezier_vals[4], Bezier_vals[10]};
    std::valarray<double> linear_vals13 = {v1, v3, linear_vals[4], linear_vals[10]};
    if (get_sign(Bezier_vals13.max()) != get_sign(Bezier_vals13.min())) {
        double diff13 = abs(Bezier_vals13 - linear_vals13).max();
        distance_error = std::max(distance_error, diff13 / tan13);
    }
    // distance error on edge 23
    std::valarray<double> Bezier_vals23 = {v2, v3, Bezier_vals[6], Bezier_vals[11]};
    std::valarray<double> linear_vals23 = {v2, v3, linear_vals[6], linear_vals[11]};
    if (get_sign(Bezier_vals23.max()) != get_sign(Bezier_vals23.min())) {
        double diff23 = abs(Bezier_vals23 - linear_vals23).max();
        distance_error = std::max(distance_error, diff23 / tan23);
    }
    //
    return distance_error;
}