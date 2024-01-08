//
//  main.cpp
//  tet_subdivision
//
//  Created by Yiwen Ju on 12/2/23.
//

#ifndef subdivide_multi_h
#define subdivide_multi_h


#include <iostream>
#include <string>
#include <valarray>
#include <array>
#include <SmallVector.h>
#include <timer.h>
#include <convex_hull_membership/contains.h>

using namespace std;


int sub_call_two = 0;
int sub_call_three = 0;

const std::array<std::array<double, 4>, 16> coeff = {{{2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1}, {0, 2, 1, 0},{0, 2, 0, 1}, {1, 2, 0, 0}, {0, 0, 2, 1}, {1, 0, 2, 0},{0, 1, 2, 0}, {1, 0, 0, 2}, {0, 1, 0, 2}, {0, 0, 1, 2},{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}}};
// constant matrix to build linear values with

int get_sign(double x) {
    return (x > 0) ? 1 : -1;
}

double dot(const std::array<double, 3> &a, const std::array<double, 3> &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double dot(const std::array<double, 2> &a, const std::array<double, 2> &b) {
    return a[0] * b[0] + a[1] * b[1];
}

double dot(const llvm_vecsmall::SmallVector<double, 20> &a, llvm_vecsmall::SmallVector<double, 20> &b) {
    double ret = 0;
    for (int i = 0; i < a.size(); i++){
        ret += a[i] * b[i];
    }
    return ret;
}

double norm(const std::array<double, 3> &a) {
    return sqrt(dot(a, a));
}

std::array<double, 3> getVec(const std::array<double, 3> &p1, const std::array<double, 3> &p2){
    std::array<double, 3> vec;
    for (int i = 0; i < 3; i++){
        vec[i] = p1[i] - p2[i];
    }
    return vec;
}

std::array<double, 3> vecPlus(const std::array<double, 3> &v1, const std::array<double, 3> &v2){
    std::array<double, 3> vec;
    for (int i = 0; i < 3; i++){
        vec[i] = v1[i] + v2[i];
    }
    return vec;
}

std::array<double, 2> perp(const std::array<double, 2> &a){
    return {-a[1], a[0]};
}

llvm_vecsmall::SmallVector<double, 20> cross(const std::array<double, 3> &a, const std::array<double, 3> &b) {
    llvm_vecsmall::SmallVector<double, 20> c(3);
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}
//
double det(const std::array<double, 3>& vec1,
           const std::array<double, 3>& vec2,
           const std::array<double, 3>& vec3) {
    return vec1[0] * vec2[1] * vec3[2] +
    vec1[1] * vec2[2] * vec3[0] +
    vec1[2] * vec2[0] * vec3[1] -
    vec1[2] * vec2[1] * vec3[0] -
    vec1[1] * vec2[0] * vec3[2] -
    vec1[0] * vec2[2] * vec3[1];
}
//
double det(const std::array<double, 2>& vec1,
           const std::array<double, 2>& vec2) {
    return vec1[0] * vec2[1] - vec1[1] * vec2[0];
}
//
std::array<double, 40> transpose2d(const std::array<std::array<double, 20>, 2>& matrix) {
    std::array<double, 40> transposed{};
    for (size_t i = 0; i < 20; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            transposed[i * 2 + j] = matrix[j][i];
        }
    }
    return transposed;
}

std::array<double, 60> transpose3d(const std::array<std::array<double, 20>, 3>& matrix) {
    std::array<double, 60> transposed{};
    for (size_t i = 0; i < 20; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            transposed[i * 3 + j] = matrix[j][i];
        }
    }
    return transposed;
}

llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> matrixMultiply(const llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> &matA, const llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> &matB) {
    size_t rowsA = matA.size();
    size_t colsA = matA[0].size();
    size_t colsB = matB[0].size();
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> result(rowsA, llvm_vecsmall::SmallVector<double, 20>(colsB));
    for (size_t i = 0; i < rowsA; ++i) {
        for (size_t j = 0; j < colsB; ++j) {
            for (size_t k = 0; k < colsA; ++k) {
                result[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }
    return result;
}


double subTet(std::array<std::array<double, 3>,4> &pts,
              const llvm_vecsmall::SmallVector<std::array<double,4>, 20> &vals,
              const llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> &grads, const double threshold, bool& active) {
    std::array<double, 3> p0 = pts[0], p1 = pts[1], p2 = pts[2], p3 = pts[3];
    std::array<double, 3> vec1, vec2, vec3;
    for (int i = 0; i < 3; i++){
        vec1[i] = p1[i] - p0[i];
        vec2[i] = p2[i] - p0[i];
        vec3[i] = p3[i] - p0[i];
    }
    double D = det(vec1, vec2, vec3);
    double sqD = D*D;
    double tetEdgeLen[] = {norm(getVec(p1, p0)),norm(getVec(p2,p0)), norm(getVec(p3,p0)), norm(getVec(p2,p1)), norm(getVec(p3,p1)), norm(getVec(p3,p2))};
    double score = *std::max_element(tetEdgeLen, tetEdgeLen + 6); // find the largest edge length using 6 edges.
    const size_t funcNum = vals.size();
    llvm_vecsmall::SmallVector<std::array<double, 3>, 20> gradList(funcNum);
    llvm_vecsmall::SmallVector<std::array<double, 20>, 20> valList(funcNum);
    llvm_vecsmall::SmallVector<std::array<double, 16>, 20> diffList(funcNum);
    llvm_vecsmall::SmallVector<double, 20> errorList(funcNum);
    llvm_vecsmall::SmallVector<bool, 20> activeTF(funcNum);
    int activeNum = 0;
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<bool, 20>, 20> zeroXResult(funcNum, llvm_vecsmall::SmallVector<bool, 20>(funcNum));
    //valarray<double> funcInt[funcNum];
    
    //single function linearity check:
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        double v0 = vals[funcIter][0], v1 = vals[funcIter][1], v2 = vals[funcIter][2], v3 = vals[funcIter][3];
        std::array<double, 3> g0 = grads[funcIter][0], g1 = grads[funcIter][1], g2 = grads[funcIter][2], g3 = grads[funcIter][3];
        //        std:: cout << v1 << std::endl;
        //        std:: cout << g0[0] << ", " << g0[1] << ", " << g0[2] << std::endl;
        double d1 = v1-v0, d2 = v2-v0, d3 = v3-v0;
        llvm_vecsmall::SmallVector<double, 20> unNormF = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({{d1, d2, d3}}), llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({cross(vec2, vec3), cross(vec3, vec1), cross(vec1, vec2)}))[0];
        //std:: cout << unNormF[0] << ", " << unNormF[1] << ", " << unNormF[2] << std::endl;
        gradList[funcIter] = {unNormF[0],unNormF[1],unNormF[2]};
        // Bezier control points
        std::array<double, 3> v0s = {v0 + dot(g0, getVec(p1 , p0)) / 3, v0 + dot(g0, getVec(p2 , p0)) / 3, v0 + dot(g0, getVec(p3 , p0)) / 3};
        std::array<double, 3> v1s = {v1 + dot(g1, getVec(p2 , p1)) / 3, v1 + dot(g1, getVec(p3 , p1)) / 3, v1 + dot(g1, getVec(p0 , p1)) / 3};
        std::array<double, 3> v2s = {v2 + dot(g2, getVec(p3 , p2)) / 3, v2 + dot(g2, getVec(p0 , p2)) / 3, v2 + dot(g2, getVec(p1 , p2)) / 3};
        std::array<double, 3> v3s = {v3 + dot(g3, getVec(p0 , p3)) / 3, v3 + dot(g3, getVec(p1 , p3)) / 3, v3 + dot(g3, getVec(p2 , p3)) / 3};
        double e0 = (v1s[0] + v1s[1] + v2s[0] + v2s[2] + v3s[1] + v3s[2]) / 6;
        double vMid0 = e0 + (e0 - (v1 + v2 + v3) / 3) / 2;
        double e1 = (v0s[1] + v0s[2] + v2s[0] + v2s[1] + v3s[0] + v3s[2]) / 6;
        double vMid1 = e1 + (e1 - (v0 + v2 + v3) / 3) / 2;
        double e2 = (v0s[0] + v0s[2] + v1s[1] + v1s[2] + v3s[0] + v3s[1]) / 6;
        double vMid2 = e2 + (e2 - (v0 + v1 + v3) / 3) / 2;
        double e3 = (v0s[0] + v0s[1] + v1s[0] + v1s[2] + v2s[1] + v2s[2]) / 6;
        double vMid3 = e3 + (e3 - (v0 + v1 + v2) / 3) / 2;

        //storing bezier and linear info for later linearity comparison
        valList[funcIter] = {v0, v1, v2, v3, v0s[0], v0s[1], v0s[2], v1s[0], v1s[1], v1s[2], v2s[0], v2s[1], v2s[2],
            v3s[0], v3s[1], v3s[2], vMid0, vMid1, vMid2, vMid3};
        std::array<double, 16> bezier_vals = {v0s[0], v0s[1], v0s[2], v1s[0], v1s[1], v1s[2],
            v2s[0], v2s[1], v2s[2], v3s[0], v3s[1], v3s[2],
            vMid0, vMid1, vMid2, vMid3};
        std::array<double, 16> linear_vals;
        for (size_t i = 0; i < coeff.size(); ++i) {
            linear_vals[i] = (v0 * coeff[i][0] + v1 * coeff[i][1] + v2 * coeff[i][2] + v3 * coeff[i][3]) / 3.0;
        }
        double diff = 0;
        std::array<double, 16> difflist_elements;
        for (int i = 0; i < 16; ++i) {
            difflist_elements[i] = bezier_vals[i] - linear_vals[i];
            diff = std::max(diff, std::abs(difflist_elements[i]));
        }
        diffList[funcIter] = difflist_elements;
        errorList[funcIter] = diff;

        Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        activeTF[funcIter] = get_sign(*std::max_element(valList[funcIter].begin(), valList[funcIter].end())) == get_sign(*std::min_element(valList[funcIter].begin(), valList[funcIter].end())) ? false : true;
        if (activeTF[funcIter]){
            if (!active){
                active = true;
            }
            activeNum++;
            double lhs = errorList[funcIter] * errorList[funcIter] * sqD;
            double rhs = threshold * threshold * dot(gradList[funcIter], gradList[funcIter]);
            if (lhs > rhs) {
                return score;
            }
        }
        single_timer.Stop();
    }
    
    Timer get_func_timer(getActiveMuti, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    if(activeNum < 2)
        return -1;
    llvm_vecsmall::SmallVector<int, 20> activeFunc(activeNum);
    int activeFuncIter = 0;
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        if (activeTF[funcIter]){
            activeFunc[activeFuncIter] = funcIter;
            activeFuncIter++;
        }
    }
    const int pairNum = activeNum * (activeNum-1)/2, triNum = activeNum * (activeNum-1) * (activeNum - 2)/ 6;
    llvm_vecsmall::SmallVector<array<int, 2>,40> pair(pairNum);
    llvm_vecsmall::SmallVector<array<int, 3>, 100> triple(triNum);
    int pairIter = 0, triIter = 0;
    for (int i = 0; i < activeNum - 1; i++){
        for (int j = i + 1; j < activeNum; j++){
            pair[pairIter] = {activeFunc[i], activeFunc[j]};
            pairIter ++;
            if (j < activeNum - 1){
                for (int k = j + 1; k < activeNum; k++){
                    triple[triIter] = {activeFunc[i], activeFunc[j], activeFunc[k]};
                    triIter ++;
                }
            }
        }
    }
    get_func_timer.Stop();
    // 2-function checks
    {
        Timer timer(twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int pairIter = 0; pairIter < pairNum; pairIter ++){
            
            std::array<double, 40> nPoints = transpose2d({valList[pair[pairIter][0]], valList[pair[pairIter][1]]});// X0, Y0, X1, Y1, ...
            std::array<double, 2> query = {0.0, 0.0}; // X, Y
            sub_call_two ++;
            Timer sub_timer(sub_twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<2, double>(nPoints, query);
            sub_timer.Stop();
            if (zeroX){
                zeroXResult[pair[pairIter][0]][pair[pairIter][1]] = true;
                zeroXResult[pair[pairIter][1]][pair[pairIter][0]] = true;

                // two function linearity test:
                std::array<double, 2> w1 = {dot(gradList[pair[pairIter][0]], gradList[pair[pairIter][0]]), dot(gradList[pair[pairIter][0]], gradList[pair[pairIter][1]])};
                std::array<double, 2> w2 = {dot(gradList[pair[pairIter][0]], gradList[pair[pairIter][1]]), dot(gradList[pair[pairIter][1]], gradList[pair[pairIter][1]])};
                double E = det(w1, w2);
                std::array<double, 2>invPerp_w2 = perp(w2);
                for(int i = 0; i < 2; i++){
                    invPerp_w2[i] *= -1;
                }
                std::array<double, 2> perp_w1 = perp(w1);
                llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> H = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({llvm_vecsmall::SmallVector<double, 20>(std::begin(invPerp_w2), std::end(invPerp_w2)), llvm_vecsmall::SmallVector<double, 20>(std::begin(perp_w1), std::end(perp_w1))}), llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({llvm_vecsmall::SmallVector<double, 20>(std::begin(gradList[pair[pairIter][0]]), std::end(gradList[pair[pairIter][0]])), llvm_vecsmall::SmallVector<double, 20>(std::begin(gradList[pair[pairIter][1]]), std::end(gradList[pair[pairIter][1]]))}));

                //find the largest max error (max squared gamma: the LHS of the equation) among all 16 bezier control points
                double maxGammaSq = 0;
                for (int i = 0; i < 16; i++){
                    //std::cout << diffList[pair[pairIter][0]][i] << " " << diffList[pair[pairIter][1]][i] << std::endl;
                    llvm_vecsmall::SmallVector<double, 20> unNormDis = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({{diffList[pair[pairIter][0]][i], diffList[pair[pairIter][1]][i]}}), llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>(H))[0];
                    double currError = sqD * dot(unNormDis, unNormDis);
                    if (maxGammaSq < currError)
                        maxGammaSq = currError;
                }
                if (maxGammaSq > threshold*threshold * E*E){
                    return score;
                }
            }
        }
        timer.Stop();
    }
    if(activeNum < 3)
        return -1;
    {
        Timer timer(threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int triIter = 0; triIter < triNum; triIter ++){
            if(!(zeroXResult[triple[triIter][0]][triple[triIter][1]]&&zeroXResult[triple[triIter][0]][triple[triIter][2]]&&zeroXResult[triple[triIter][1]][triple[triIter][2]]))
                continue;
            std::array<double, 60> nPoints = transpose3d({valList[pair[pairIter][0]], valList[pair[pairIter][1]], valList[triple[triIter][2]]});
            std::array<double, 3> query = {0.0, 0.0, 0.0}; // X, Y, Z
            sub_call_three ++;
            Timer sub_timer(sub_threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<3, double>(nPoints, query);
            sub_timer.Stop();
            if (zeroX){
                std::array<double, 3> fi = gradList[triple[triIter][0]];
                std::array<double, 3> fj = gradList[triple[triIter][1]];
                std::array<double, 3> fk = gradList[triple[triIter][2]];
                double E = det(fi, fj, fk);
                llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> H = {cross(fj, fk), cross(fk, fi), cross(fi, fj)};
                double maxGammaSq = 0;
                for (int i = 0; i < 16; i++){
                    llvm_vecsmall::SmallVector<double, 20> unNormDis = matrixMultiply({{diffList[triple[triIter][0]][i], diffList[triple[triIter][1]][i], diffList[triple[triIter][2]][i]}}, H)[0];
                    double currError = sqD * dot(unNormDis, unNormDis);
                    if (maxGammaSq < currError)
                        maxGammaSq = currError;
                }
                if (maxGammaSq > threshold*threshold * E*E)
                    return score;
            }
        }
        timer.Stop();
    }
    return -1;
}

#endif /* subdivide_multi_h */
