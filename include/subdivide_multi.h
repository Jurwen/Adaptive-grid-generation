//
//  main.cpp
//  tet_subdivision
//
//  Created by Yiwen Ju on 12/2/23.
//




#include <iostream>
#include <string>
#include <valarray>
#include <array>
#include "gurobi_c++.h"
#include "timer.h"

using namespace std;

//{total time, time spent on double functions, time spent on triple functions, time spent on zero crossing test}

int gurobi_call_two = 0;
int gurobi_call_three = 0;

const valarray<std::array<double, 4>> coeff = {
    {2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1}, {0, 2, 1, 0},
    {0, 2, 0, 1}, {1, 2, 0, 0}, {0, 0, 2, 1}, {1, 0, 2, 0},
    {0, 1, 2, 0}, {1, 0, 0, 2}, {0, 1, 0, 2}, {0, 0, 1, 2},
    {0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}
}; // constant matrix to build linear values with

int get_sign(double x) {
    return (x > 0) ? 1 : -1;
}

double dot(const valarray<double> &a, const valarray<double> &b) {
    return (a * b).sum();
}

valarray<double> normalize(const valarray<double> &a) {
    return a / sqrt(dot(a, a));
}

double norm(const valarray<double> &a) {
    return sqrt(dot(a, a));
}

valarray<double> perp(const valarray<double> &a){
    return {-a[1], a[0]};
}

valarray<double> cross(const valarray<double> &a, const valarray<double> &b) {
    valarray<double> c(3);
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

double det(const std::valarray<double>& vec1,
           const std::valarray<double>& vec2,
           const std::valarray<double>& vec3) {
    return vec1[0] * vec2[1] * vec3[2] +
    vec1[1] * vec2[2] * vec3[0] +
    vec1[2] * vec2[0] * vec3[1] -
    vec1[2] * vec2[1] * vec3[0] -
    vec1[1] * vec2[0] * vec3[2] -
    vec1[0] * vec2[2] * vec3[1];
}

double det(const std::valarray<double>& vec1,
           const std::valarray<double>& vec2) {
    return vec1[0] * vec2[1] -
    vec1[1] * vec2[0];
}

valarray<valarray<double>> transpose(const valarray<valarray<double>>& matrix) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    valarray<valarray<double>> transposed(valarray<double>(rows), cols);
    
    for (size_t i = 0; i < cols; ++i) {
        for (size_t j = 0; j < rows; ++j) {
            transposed[i][j] = matrix[j][i];
        }
    }
    return transposed;
}

valarray<valarray<double>> matrixMultiply(const valarray<valarray<double>> &matA, const valarray<valarray<double>> &matB) {
    size_t rowsA = matA.size();
    size_t colsA = matA[0].size();
    size_t colsB = matB[0].size();
    valarray<valarray<double>> result(valarray<double>(colsB),rowsA);
    
    for (size_t i = 0; i < rowsA; ++i) {
        for (size_t j = 0; j < colsB; ++j) {
            for (size_t k = 0; k < colsA; ++k) {
                result[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }
    return result;
}

double subTet(std::array<valarray<double>,4> &pts,
              const valarray<std::array<double,4>> &vals,
              const valarray<std::array<valarray<double>,4>> &grads, const double threshold, bool& active, GRBEnv& env) {
    valarray<double> p0 = pts[0], p1 = pts[1], p2 = pts[2], p3 = pts[3];
    valarray<double> vec1 = p1 - p0, vec2 = p2 - p0, vec3 = p3 - p0;
    double D = det(vec1, vec2, vec3);
    double sqD = D*D;
    double tetEdgeLen[] = {norm(p1-p0),norm(p2-p0), norm(p3-p0), norm(p2-p1), norm(p3-p1), norm(p3-p2)};
    double score = *std::max_element(tetEdgeLen, tetEdgeLen + 6); // find the largest edge length using 6 edges.
    const size_t funcNum = vals.size();
    valarray<double> gradList[funcNum];
    valarray<double> valList[funcNum];
    valarray<double> diffList[funcNum];
    double errorList[funcNum];
    bool activeTF[funcNum];
    int activeNum = 0;
    valarray<valarray<bool>> zeroXResult(valarray<bool>(funcNum), funcNum);
    //valarray<double> funcInt[funcNum];
    
    //single function linearity check:
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        double v0 = vals[funcIter][0], v1 = vals[funcIter][1], v2 = vals[funcIter][2], v3 = vals[funcIter][3];
        valarray<double> g0 = grads[funcIter][0], g1 = grads[funcIter][1], g2 = grads[funcIter][2], g3 = grads[funcIter][3];
        //        std:: cout << v1 << std::endl;
        //        std:: cout << g0[0] << ", " << g0[1] << ", " << g0[2] << std::endl;
        double d1 = v1-v0, d2 = v2-v0, d3 = v3-v0;
        valarray<double> unNormF = matrixMultiply({{d1, d2, d3}}, {cross(vec2, vec3), cross(vec3, vec1), cross(vec1, vec2)})[0];
        //std:: cout << unNormF[0] << ", " << unNormF[1] << ", " << unNormF[2] << std::endl;
        gradList[funcIter] = unNormF;
        
        // Bezier control points
        valarray<double> v0s = {v0 + dot(g0, p1 - p0) / 3, v0 + dot(g0, p2 - p0) / 3, v0 + dot(g0, p3 - p0) / 3};
        valarray<double> v1s = {v1 + dot(g1, p2 - p1) / 3, v1 + dot(g1, p3 - p1) / 3, v1 + dot(g1, p0 - p1) / 3};
        valarray<double> v2s = {v2 + dot(g2, p3 - p2) / 3, v2 + dot(g2, p0 - p2) / 3, v2 + dot(g2, p1 - p2) / 3};
        valarray<double> v3s = {v3 + dot(g3, p0 - p3) / 3, v3 + dot(g3, p1 - p3) / 3, v3 + dot(g3, p2 - p3) / 3};
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
        valarray<double> bezier_vals = {v0s[0], v0s[1], v0s[2], v1s[0], v1s[1], v1s[2],
            v2s[0], v2s[1], v2s[2], v3s[0], v3s[1], v3s[2],
            vMid0, vMid1, vMid2, vMid3};
        valarray<double> linear_vals(coeff.size());
        for (size_t i = 0; i < coeff.size(); ++i) {
            linear_vals[i] = v0 * coeff[i][0] + v1 * coeff[i][1] + v2 * coeff[i][2] + v3 * coeff[i][3];
        }
        linear_vals /= 3.0;
        double diff = abs(bezier_vals - linear_vals).max();
        diffList[funcIter] = bezier_vals - linear_vals;
        errorList[funcIter] = diff;
        
        activeTF[funcIter] = get_sign(valList[funcIter].max()) == get_sign(valList[funcIter].min()) ? false : true;
        
        
        if (activeTF[funcIter]){
            if (!active){
                valarray<double> valList = {v0, v1, v2,v3};
                if (get_sign(valList.max()) != get_sign(valList.min())){
                    active = true;
                }
            }
            activeNum++;
            double lhs = errorList[funcIter] * errorList[funcIter] * sqD;
            double rhs = threshold * threshold * dot(gradList[funcIter], gradList[funcIter]);
            if (lhs > rhs) {
                return score;
            }
        }
    }
    
    if(activeNum < 2)
        return -1;
    int activeFunc[activeNum];
    int activeFuncIter = 0;
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        if (activeTF[funcIter]){
            activeFunc[activeFuncIter] = funcIter;
            activeFuncIter++;
        }
    }
    
    const int pairNum = activeNum * (activeNum-1)/2, triNum = activeNum * (activeNum-1) * (activeNum - 2)/ 6;
    array<int, 2> pair[pairNum];
    array<int, 3> triple[triNum];
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
    
    // 2-function checks
    {
    bool zeroX;
        Timer timer(twoFunc, [&](auto profileResult){profileTimer += profileResult;});
        for (int pairIter = 0; pairIter < pairNum; pairIter ++){
            valarray<valarray<double>> nPoints = transpose({valList[pair[pairIter][0]], valList[pair[pairIter][1]]});
            valarray<double> cent(0.0, nPoints[0].size());
            for (size_t i = 0; i < nPoints.size(); i++){
                //std::cout << nPoints[i][0] << " " << nPoints[i][1] << std::endl;
                for (size_t j = 0; j < nPoints[0].size(); j++){
                    cent[j] += nPoints[i][j];
                }
            }
            cent /= 20;
            if(cent[0] == 0.0 & cent[1] == 0.0 ){
                zeroX = true;
            }else{
                Timer timer(gurobi_twoFunc, [&](auto profileResult){profileTimer += profileResult;});
                gurobi_call_two += 1;
                try {
                    // Create an environment
//                    GRBEnv env = GRBEnv(true);
//                    env.set(GRB_IntParam_OutputFlag, 0);
//                    env.set("LogFile", "");
//                    env.start();
//
                    // Create an empty model
                    GRBModel model = GRBModel(env);

                    // Create variables
                    GRBVar x = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x");
                    GRBVar y = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y");
                    GRBVar ep = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "epsilon");

                    model.setObjective(ep * 1.0, GRB_MAXIMIZE);
                    for (size_t i = 0; i < nPoints.size();i++){
                        model.addConstr(x * nPoints[i][0] + y * nPoints[i][1] - ep >= 0);
                    }
                    model.addConstr(cent[0] * x + cent[1] * y - 1 == 0);
                    // Optimize model
                    model.optimize();
                    //cout << "largest epsilon value is " << ep.get(GRB_DoubleAttr_X) << endl;
                    zeroX = ep.get(GRB_DoubleAttr_X) <= 0.0 ? true : false;
                    //cout << zeroX << endl;
                } catch(GRBException e) {
                    cout << "Error code = " << e.getErrorCode() << endl;
                    cout << e.getMessage() << endl;
                } catch(...) {
                    cout << "Exception during optimization" << endl;
                }
            }
            if (zeroX){
                zeroXResult[pair[pairIter][0]][pair[pairIter][1]] = true;
                zeroXResult[pair[pairIter][1]][pair[pairIter][0]] = true;
                
                // two function linearity test:
                valarray<double> w1 = {dot(gradList[pair[pairIter][0]], gradList[pair[pairIter][0]]), dot(gradList[pair[pairIter][0]], gradList[pair[pairIter][1]])};
                valarray<double> w2 = {dot(gradList[pair[pairIter][0]], gradList[pair[pairIter][1]]), dot(gradList[pair[pairIter][1]], gradList[pair[pairIter][1]])};
                double E = det(w1, w2);
                valarray<valarray<double>> H = matrixMultiply({-perp(w2), perp(w1)}, {gradList[pair[pairIter][0]], gradList[pair[pairIter][1]]});
                
                //find the largest max error (max squared gamma: the LHS of the equation) among all 16 bezier control points
                double maxGammaSq = 0;
                for (int i = 0; i < 16; i++){
                    //std::cout << diffList[pair[pairIter][0]][i] << " " << diffList[pair[pairIter][1]][i] << std::endl;
                    valarray<double> unNormDis = matrixMultiply({{diffList[pair[pairIter][0]][i], diffList[pair[pairIter][1]][i]}}, H)[0];
                    double currError = sqD * dot(unNormDis, unNormDis);
                    if (maxGammaSq < currError)
                        maxGammaSq = currError;
                }
                if (maxGammaSq > threshold*threshold * E*E){
                    return score;
                }
            }
        }
    }
    if(activeNum < 3)
        return -1;
    {
        Timer timer(threeFunc, [&](auto profileResult){profileTimer += profileResult;});
        bool zeroX;
        for (int triIter = 0; triIter < triNum; triIter ++){
            if(!(zeroXResult[triple[triIter][0]][triple[triIter][1]]&&zeroXResult[triple[triIter][0]][triple[triIter][2]]&&zeroXResult[triple[triIter][1]][triple[triIter][2]]))
                continue;
            valarray<valarray<double>> nPoints = transpose({valList[triple[triIter][0]], valList[triple[triIter][1]], valList[triple[triIter][2]]});
            valarray<double> cent(0.0, nPoints[0].size());
            for (size_t i = 0; i < nPoints.size(); i++){
                for (size_t j = 0; j < nPoints[0].size(); j++){
                    cent[j] += nPoints[i][j];
                }
            }
            cent /= 20;
            if(cent[0] == 0.0 && cent[1] == 0.0 && cent[2] == 0.0 ){
                zeroX = true;
            }else{
                Timer timer(gurobi_threeFunc, [&](auto profileResult){profileTimer += profileResult;});
                gurobi_call_three += 1;
                try {
                    // Create an empty model
                    GRBModel model = GRBModel(env);
                    // Create variables
                    GRBVar x = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x");
                    GRBVar y = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y");
                    GRBVar z = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z");
                    GRBVar ep = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "epsilon");

                    // Create an environment
                    model.setObjective(ep * 1.0, GRB_MAXIMIZE);
                    for (size_t i = 0; i < nPoints.size();i++){
                        model.addConstr(x * nPoints[i][0] + y * nPoints[i][1] + z * nPoints[i][2] - ep >= 0);
                    }
                    model.addConstr(cent[0] * x + cent[1] * y + cent[2] * z - 1 == 0);
                    // Optimize model
                    model.optimize();
                    //cout << "largest epsilon value is " << ep.get(GRB_DoubleAttr_X) << endl;
                    zeroX = ep.get(GRB_DoubleAttr_X) <= 0.0 ? true : false;
                    //cout << zeroX << endl;
                } catch(GRBException e) {
                    cout << "Error code = " << e.getErrorCode() << endl;
                    cout << e.getMessage() << endl;
                } catch(...) {
                    cout << "Exception during optimization" << endl;
                }
            }
            
            if (zeroX){
                valarray<double> fi = gradList[triple[triIter][0]];
                valarray<double> fj = gradList[triple[triIter][1]];
                valarray<double> fk = gradList[triple[triIter][2]];
                double E = det(fi, fj, fk);
                valarray<valarray<double>> H = {cross(fj, fk), cross(fk, fi), cross(fi, fj)};
                double maxGammaSq = 0;
                for (int i = 0; i < 16; i++){
                    valarray<double> unNormDis = matrixMultiply({{diffList[triple[triIter][0]][i], diffList[triple[triIter][1]][i], diffList[triple[triIter][2]][i]}}, H)[0];
                    double currError = sqD * dot(unNormDis, unNormDis);
                    if (maxGammaSq < currError)
                        maxGammaSq = currError;
                }
                if (maxGammaSq > threshold*threshold * E*E)
                    return score;
            }
        }
    }
    return -1;
}

