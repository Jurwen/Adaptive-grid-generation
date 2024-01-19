//
//  main.cpp
//  tet_subdivision
//
//  Created by Yiwen Ju on 12/2/23.
//

#ifndef subdivide_multi_h
#define subdivide_multi_h
//#define No_Multi_Check
//#define Only_ZeroX
//#define Only_Geometry
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

llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<array<int, 4>, 100>, 3>, 20> multiple_indices;

enum geo_obj {
    IA,
    CSG,
    MI
};

int GLOBAL_METHOD = IA;
bool curve_network = false;

enum csg_operations{
    Intersection,
    Union,
    Negation
};

struct csg_unit{
    int operation;
    std::array<int, 2> elements;
};

bool load_csgTree(const std::string filename, llvm_vecsmall::SmallVector<csg_unit, 20>& tree){
    using json = nlohmann::json;
    std::ifstream fin(filename.c_str());
    //llvm_vecsmall::SmallVector<csg_unit, 20> tree = {};
    if (!fin)
    {
        std::cout << "function file not exist!" << std::endl;
        return false;
    }
    json tree_data;
    fin >> tree_data;
    fin.close();
    //
    size_t n_units = tree_data.size();
    tree.resize(n_units);
    for (size_t j = 0 ; j < n_units; j++){
        std::string type = tree_data[j]["type"].get<std::string>();
        std::array<int, 2> elements;
        for (int i = 0; i < 2; i ++){
            elements[i] = tree_data[j]["elements"][i].get<int>();
        }
        if (type == "Intersection"){
            tree[j] = {0, elements};
        }else if (type == "Union"){
            tree[j] = {1, elements};
        }else if (type == "Negation"){
            tree[j] = {2, elements};
        }
    }
    return true;
}

llvm_vecsmall::SmallVector<csg_unit, 20> GLOBAL_CSGTREE ={{0, {2, 8}}, {1, {3, -9}}, {0,{-1,4}}, {0, {-2, 5}}, {0, {-3,6}}, {0, {-4, 7}}, {0, {-5, -6}},{2,{9, 0}}, {1, {-7, -8}}};

std::pair<array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> iterTree(const llvm_vecsmall::SmallVector<csg_unit, 20>csgTree,const int curNode,const llvm_vecsmall::SmallVector<array<double , 2>, 20> funcInt){
    csg_unit curUnit = csgTree[curNode - 1];
    array<double, 2> interval, childInt1, childInt2;
    llvm_vecsmall::SmallVector<int, 20> af(funcInt.size(), 1), childAF1(funcInt.size(), 1), childAF2(funcInt.size(), 1);
    if (curUnit.elements[0] > 0){
        std::pair<array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> child1 = iterTree(csgTree, curUnit.elements[0], funcInt);
        childInt1 = child1.first;
        childAF1 = child1.second;
    }else{
        childInt1 = funcInt[-curUnit.elements[0] - 1];
        childAF1[-curUnit.elements[0] - 1] = 0;
    }
    if (childInt1[0] * childInt1[1]>0){
        for (size_t i = 0; i < childAF1.size(); i++){
            childAF1[i] = 1;
        }
    }
    if (curUnit.operation != Negation){
        if (curUnit.elements[1] > 0){
            std::pair<array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> child2 = iterTree(csgTree, curUnit.elements[1], funcInt);
            childInt2 = child2.first;
            childAF2 = child2.second;
        }else{
            childInt2 = funcInt[-curUnit.elements[1] - 1];
            childAF2[-curUnit.elements[1] - 1] = 0;
        }
    }
    if (childInt2[0] * childInt2[1]>0){
        for (size_t i = 0; i < childAF2.size(); i++){
            childAF2[i] = 1;
        }
    }
    switch (curUnit.operation){
        case Intersection:
            interval = {std::max(childInt1[0], childInt2[0]), std::max(childInt1[1], childInt2[1])};
            if(interval[0]*interval[1] <= 0){
                for (int i = 0; i < funcInt.size(); i++){
                    af[i] = childAF1[i] * childAF2[i];
                }
            }
            break;
        case Union:
            interval = {std::min(childInt1[0], childInt2[0]), std::min(childInt1[1], childInt2[1])};
            if(interval[0]*interval[1] <= 0){
                for (int i = 0; i < funcInt.size(); i++){
                    af[i] = childAF1[i] * childAF2[i];
                }
            }
            break;
        case Negation:
            interval = {std::min(-childInt1[0], -childInt1[1]), std::max(-childInt1[0], -childInt1[1])};
            if(interval[0]*interval[1] <= 0)
                af = childAF1;
            break;
        default:
            std::cout << "not a valid CSG operation" << std::endl;
    }
    return std::pair(interval, af);
}

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


bool subTet(std::array<std::array<double, 3>,4> &pts,
            const llvm_vecsmall::SmallVector<std::array<double,4>, 20> &vals,
            const llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> &grads, const double threshold, bool& active) {
    std::array<double, 3> p0 = pts[0], p1 = pts[1], p2 = pts[2], p3 = pts[3];
    std::array<double, 3> vec1 = getVec(p1, p0), vec2 = getVec(p2, p0), vec3 = getVec(p3, p0),
    vec4 = getVec(p2, p1), vec5 = getVec(p3, p1), vec6 = getVec(p3, p2);
    double D = det(vec1, vec2, vec3);
    double sqD = D*D;
    bool score = true;
    //    double tetEdgeLen[] = {norm(getVec(p1, p0)),norm(getVec(p2,p0)), norm(getVec(p3,p0)), norm(getVec(p2,p1)), norm(getVec(p3,p1)), norm(getVec(p3,p2))};
    //    double score = *std::max_element(tetEdgeLen, tetEdgeLen + 6); // find the largest edge length using 6 edges.
    const size_t funcNum = vals.size();
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> crossMatrix = {cross(vec2, vec3), cross(vec3, vec1), cross(vec1, vec2)};
    llvm_vecsmall::SmallVector<std::array<double, 3>, 20> gradList(funcNum);
    llvm_vecsmall::SmallVector<std::array<double, 20>, 20> valList(funcNum);
    llvm_vecsmall::SmallVector<std::array<double, 16>, 20> diffList(funcNum);
    llvm_vecsmall::SmallVector<double, 20> errorList(funcNum);
    llvm_vecsmall::SmallVector<bool, 20> activeTF(funcNum);
    llvm_vecsmall::SmallVector<array<double , 2>, 20> funcInt(funcNum);
    int activeNum = 0;
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<bool, 20>, 20> zeroXResult(funcNum, llvm_vecsmall::SmallVector<bool, 20>(funcNum));
    
    //single function linearity check:
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        double v0 = vals[funcIter][0], v1 = vals[funcIter][1], v2 = vals[funcIter][2], v3 = vals[funcIter][3];
        std::array<double, 3> g0 = grads[funcIter][0], g1 = grads[funcIter][1], g2 = grads[funcIter][2], g3 = grads[funcIter][3];
        // Bezier control points
        std::array<double, 3> v0s = {v0 + dot(g0, vec1) / 3, v0 + dot(g0, vec2) / 3, v0 + dot(g0, vec3) / 3};
        std::array<double, 3> v1s = {v1 + dot(g1, vec4) / 3, v1 + dot(g1, vec5) / 3, v1 - dot(g1, vec1) / 3};
        std::array<double, 3> v2s = {v2 + dot(g2, vec6) / 3, v2 - dot(g2, vec2) / 3, v2 - dot(g2, vec4) / 3};
        std::array<double, 3> v3s = {v3 - dot(g3, vec3) / 3, v3 - dot(g3, vec5) / 3, v3 - dot(g3, vec6) / 3};
        //double e0 = (v1s[0] + v1s[1] + v2s[0] + v2s[2] + v3s[1] + v3s[2]) / 6;
        double vMid0 = (9 * (v1s[0] + v1s[1] + v2s[0] + v2s[2] + v3s[1] + v3s[2]) / 6 - v1 - v2 - v3)/ 6;
        //double e1 = (v0s[1] + v0s[2] + v2s[0] + v2s[1] + v3s[0] + v3s[2]) / 6;
        double vMid1 =(9 * (v0s[1] + v0s[2] + v2s[0] + v2s[1] + v3s[0] + v3s[2]) / 6 - v0 - v2 - v3)/ 6;
        //double e2 = (v0s[0] + v0s[2] + v1s[1] + v1s[2] + v3s[0] + v3s[1]) / 6;
        double vMid2 =(9 * (v0s[0] + v0s[2] + v1s[1] + v1s[2] + v3s[0] + v3s[1]) / 6 - v0 - v1 - v3)/ 6;
        //double e3 = (v0s[0] + v0s[1] + v1s[0] + v1s[2] + v2s[1] + v2s[2]) / 6;
        double vMid3 =(9 * (v0s[0] + v0s[1] + v1s[0] + v1s[2] + v2s[1] + v2s[2]) / 6 - v0 - v1 - v2)/ 6;
        
        //storing bezier and linear info for later linearity comparison
        valList[funcIter] = {v0, v1, v2, v3, v0s[0], v0s[1], v0s[2], v1s[0], v1s[1], v1s[2], v2s[0], v2s[1], v2s[2],
            v3s[0], v3s[1], v3s[2], vMid0, vMid1, vMid2, vMid3};
        if(GLOBAL_METHOD == IA){
            Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            activeTF[funcIter] = get_sign(*std::max_element(valList[funcIter].begin(), valList[funcIter].end())) == get_sign(*std::min_element(valList[funcIter].begin(), valList[funcIter].end())) ? false : true;
            single_timer.Stop();
            if (activeTF[funcIter]){
                if (GLOBAL_METHOD == IA){
                    if (!active){
                        active = true;
                    }
                }
                activeNum++;
                double d1 = v1-v0, d2 = v2-v0, d3 = v3-v0;
                llvm_vecsmall::SmallVector<double, 20> unNormF = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({{d1, d2, d3}}), crossMatrix)[0];
                gradList[funcIter] = {unNormF[0],unNormF[1],unNormF[2]};
                for (int i = 0; i < 16; ++i) {
                    diffList[funcIter][i] = valList[funcIter][i + 4] -
                    (v0 * coeff[i][0] + v1 * coeff[i][1] + v2 * coeff[i][2] + v3 * coeff[i][3]) / 3.0;
                }
                errorList[funcIter] = std::max(*max_element(diffList[funcIter].begin(), diffList[funcIter].end()), std::abs(*min_element(diffList[funcIter].begin(), diffList[funcIter].end())));
                Timer single2_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
                double lhs = errorList[funcIter] * errorList[funcIter] * sqD;
                double rhs;
                if (!curve_network){
                    rhs = threshold * threshold * dot(gradList[funcIter], gradList[funcIter]);
                }else{
                    rhs = numeric_limits<double>::infinity() * dot(gradList[funcIter], gradList[funcIter]);
                }
                if (lhs > rhs) {
                    single2_timer.Stop();
                    return score;
                }
                single2_timer.Stop();
            }
        }else{
            funcInt[funcIter] = {*std::min_element(valList[funcIter].begin(), valList[funcIter].end()), *std::max_element(valList[funcIter].begin(), valList[funcIter].end())};
        }
    }
    
    if (GLOBAL_METHOD == CSG){
        //Timer csg_timer(csgtree, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        std::pair<array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> csgResult = iterTree(GLOBAL_CSGTREE, 1, funcInt);
        //csg_timer.Stop();
        if(csgResult.first[0] * csgResult.first[1] > 0){
            return false;
        }else{
            for (size_t funcIter = 0; funcIter < funcNum; funcIter++){
                Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
                activeTF[funcIter] = !csgResult.second[funcIter];
                single_timer.Stop();
                if (activeTF[funcIter]){
                    if (!active){
                        active = true;
                    }
                    activeNum++;
                    double v0 = vals[funcIter][0], v1 = vals[funcIter][1], v2 = vals[funcIter][2], v3 = vals[funcIter][3];
                    double d1 = v1-v0, d2 = v2-v0, d3 = v3-v0;
                    llvm_vecsmall::SmallVector<double, 20> unNormF = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({{d1, d2, d3}}), crossMatrix)[0];
                    //std:: cout << unNormF[0] << ", " << unNormF[1] << ", " << unNormF[2] << std::endl;
                    gradList[funcIter] = {unNormF[0],unNormF[1],unNormF[2]};
                    for (int i = 0; i < 16; ++i) {
                        diffList[funcIter][i] = valList[funcIter][i + 4] -
                        (v0 * coeff[i][0] + v1 * coeff[i][1] + v2 * coeff[i][2] + v3 * coeff[i][3]) / 3.0;
                    }
                    errorList[funcIter] = std::max(*max_element(diffList[funcIter].begin(), diffList[funcIter].end()), std::abs(*min_element(diffList[funcIter].begin(), diffList[funcIter].end())));
                    Timer single2_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
                    double lhs = errorList[funcIter] * errorList[funcIter] * sqD;
                    double rhs;
                    if (!curve_network){
                        rhs = threshold * threshold * dot(gradList[funcIter], gradList[funcIter]);
                    }else{
                        rhs = numeric_limits<double>::infinity() * dot(gradList[funcIter], gradList[funcIter]);
                    }
                    if (lhs > rhs) {
                        single2_timer.Stop();
                        return score;
                    }
                    single2_timer.Stop();
                }
            }
        }
    }
    Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    if(activeNum < 2){
        single_timer.Stop();
        return false;
    }
    llvm_vecsmall::SmallVector<int, 20> activeFunc(activeNum);
    int activeFuncIter = 0;
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        if (activeTF[funcIter]){
            activeFunc[activeFuncIter] = funcIter;
            activeFuncIter++;
        }
    }
    single_timer.Stop();
    const int pairNum = activeNum * (activeNum-1)/2, triNum = activeNum * (activeNum-1) * (activeNum - 2)/ 6;
    llvm_vecsmall::SmallVector<array<int, 2>,40> pair(pairNum);
    llvm_vecsmall::SmallVector<array<int, 3>, 100> triple(triNum);
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<array<int, 4>, 100>, 3> multiples = multiple_indices[activeNum - 1];
//    int pairIt = 0, triIt = 0;
//    for (int i = 0; i < activeNum - 1; i++){
//        for (int j = i + 1; j < activeNum; j++){
//            pair[pairIt] = {activeFunc[i], activeFunc[j]};
//            pairIt ++;
//            if (j < activeNum - 1){
//                for (int k = j + 1; k < activeNum; k++){
//                    triple[triIt] = {activeFunc[i], activeFunc[j], activeFunc[k]};
//                    triIt ++;
//                }
//            }
//        }
//    }
    
    // 2-function checks
    int activeDouble_count = 0;
    {
        Timer timer(twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int pairIter = 0; pairIter < pairNum; pairIter ++){
            array<int, 2> pairIndices = {multiples[0][pairIter][0],multiples[0][pairIter][1]};
            pair[pairIter] = {activeFunc[pairIndices[0]], activeFunc[pairIndices[1]]};
            std::array<double, 40> nPoints = transpose2d({valList[pair[pairIter][0]], valList[pair[pairIter][1]]});// X0, Y0, X1, Y1, ...
            std::array<double, 2> query = {0.0, 0.0}; // X, Y
            sub_call_two ++;
            Timer sub_timer(sub_twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<2, double>(nPoints, query);
            sub_timer.Stop();
#ifdef No_Multi_Check
            zeroX = false;
#endif
#ifdef Only_Geometry
            zeroX = true;
#endif
#ifdef Only_ZeroX
            if (zeroX){
                return score;
            }else{
                return false;
            }
#endif
            if (zeroX){
                activeDouble_count++;
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
                    timer.Stop();
                    return score;
                }
            }
        }
        timer.Stop();
    }
    if(activeDouble_count < 3)
        return false;
    {
        Timer timer(threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int triIter = 0; triIter < triNum; triIter ++){
            array<int, 3> tripleIndices = {multiples[1][triIter][0], multiples[1][triIter][1], multiples[1][triIter][2]};
            triple[triIter] = {activeFunc[tripleIndices[0]], activeFunc[tripleIndices[1]], activeFunc[tripleIndices[2]]};
            if(!(zeroXResult[triple[triIter][0]][triple[triIter][1]]&&zeroXResult[triple[triIter][0]][triple[triIter][2]]&&zeroXResult[triple[triIter][1]][triple[triIter][2]]))
                continue;
            std::array<double, 60> nPoints = transpose3d({valList[triple[triIter][0]], valList[triple[triIter][1]], valList[triple[triIter][2]]});
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
                if (maxGammaSq > threshold*threshold * E*E){
                    timer.Stop();
                    return score;
                }
            }
        }
        timer.Stop();
    }
    return false;
}

bool subMI(std::array<std::array<double, 3>,4> &pts,
           const llvm_vecsmall::SmallVector<std::array<double,4>, 20> &vals,
           const llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> &grads, const double threshold, bool& active)
{
    std::array<double, 3> p0 = pts[0], p1 = pts[1], p2 = pts[2], p3 = pts[3];
    std::array<double, 3> vec1 = getVec(p1, p0), vec2 = getVec(p2, p0), vec3 = getVec(p3, p0),
    vec4 = getVec(p2, p1), vec5 = getVec(p3, p1), vec6 = getVec(p3, p2);
    double D = det(vec1, vec2, vec3);
    double sqD = D*D;
    bool score = true;
    //    double tetEdgeLen[] = {norm(getVec(p1, p0)),norm(getVec(p2,p0)), norm(getVec(p3,p0)), norm(getVec(p2,p1)), norm(getVec(p3,p1)), norm(getVec(p3,p2))};
    //    double score = *std::max_element(tetEdgeLen, tetEdgeLen + 6); // find the largest edge length using 6 edges.
    const size_t funcNum = vals.size();
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> crossMatrix = {cross(vec2, vec3), cross(vec3, vec1), cross(vec1, vec2)};
    llvm_vecsmall::SmallVector<std::array<double, 3>, 20> gradList(funcNum);
    llvm_vecsmall::SmallVector<std::array<double, 20>, 20> valList(funcNum);
    llvm_vecsmall::SmallVector<std::array<double, 16>, 20> diffList(funcNum);
    llvm_vecsmall::SmallVector<double, 20> errorList(funcNum);
    llvm_vecsmall::SmallVector<bool, 20> activeList(funcNum);
    llvm_vecsmall::SmallVector<array<double , 2>, 20> funcInt(funcNum);
    double maxLow = -1 * numeric_limits<double>::infinity();
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<bool, 20>, 20> activePair(funcNum, llvm_vecsmall::SmallVector<bool, 20>(false, funcNum));
    
    //single function linearity check:
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        double v0 = vals[funcIter][0], v1 = vals[funcIter][1], v2 = vals[funcIter][2], v3 = vals[funcIter][3];
        std::array<double, 3> g0 = grads[funcIter][0], g1 = grads[funcIter][1], g2 = grads[funcIter][2], g3 = grads[funcIter][3];
        // Bezier control points
        std::array<double, 3> v0s = {v0 + dot(g0, vec1) / 3, v0 + dot(g0, vec2) / 3, v0 + dot(g0, vec3) / 3};
        std::array<double, 3> v1s = {v1 + dot(g1, vec4) / 3, v1 + dot(g1, vec5) / 3, v1 - dot(g1, vec1) / 3};
        std::array<double, 3> v2s = {v2 + dot(g2, vec6) / 3, v2 - dot(g2, vec2) / 3, v2 - dot(g2, vec4) / 3};
        std::array<double, 3> v3s = {v3 - dot(g3, vec3) / 3, v3 - dot(g3, vec5) / 3, v3 - dot(g3, vec6) / 3};
        //double e0 = (v1s[0] + v1s[1] + v2s[0] + v2s[2] + v3s[1] + v3s[2]) / 6;
        double vMid0 = (9 * (v1s[0] + v1s[1] + v2s[0] + v2s[2] + v3s[1] + v3s[2]) / 6 - v1 - v2 - v3)/ 6;
        //double e1 = (v0s[1] + v0s[2] + v2s[0] + v2s[1] + v3s[0] + v3s[2]) / 6;
        double vMid1 =(9 * (v0s[1] + v0s[2] + v2s[0] + v2s[1] + v3s[0] + v3s[2]) / 6 - v0 - v2 - v3)/ 6;
        //double e2 = (v0s[0] + v0s[2] + v1s[1] + v1s[2] + v3s[0] + v3s[1]) / 6;
        double vMid2 =(9 * (v0s[0] + v0s[2] + v1s[1] + v1s[2] + v3s[0] + v3s[1]) / 6 - v0 - v1 - v3)/ 6;
        //double e3 = (v0s[0] + v0s[1] + v1s[0] + v1s[2] + v2s[1] + v2s[2]) / 6;
        double vMid3 =(9 * (v0s[0] + v0s[1] + v1s[0] + v1s[2] + v2s[1] + v2s[2]) / 6 - v0 - v1 - v2)/ 6;
        
        //storing bezier and linear info for later linearity comparison
        valList[funcIter] = {v0, v1, v2, v3, v0s[0], v0s[1], v0s[2], v1s[0], v1s[1], v1s[2], v2s[0], v2s[1], v2s[2],
            v3s[0], v3s[1], v3s[2], vMid0, vMid1, vMid2, vMid3};
        funcInt[funcIter] = {*std::min_element(valList[funcIter].begin(), valList[funcIter].end()), *std::max_element(valList[funcIter].begin(), valList[funcIter].end())};
        //cout << funcInt[funcIter][0] << " " << funcInt[funcIter][1] << endl;
        if (maxLow < funcInt[funcIter][0]){
            maxLow = funcInt[funcIter][0];
        }
    }
    //llvm_vecsmall::SmallVector<std::array<double, 20>, 20> active_valList;
    llvm_vecsmall::SmallVector<int, 20> activeFunc;
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        if(funcInt[funcIter][1] > maxLow){
            activeFunc.push_back(funcIter);
            //active_valList.push_back(valList[funcIter]);
        }
    }
    
    int activeNum = activeFunc.size();
    
    if(activeNum < 2)
        return false;
    
    //Timer get_func_timer(getActiveMuti, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    
    
    const int pairNum = activeNum * (activeNum-1)/2, triNum = activeNum * (activeNum-1) * (activeNum - 2)/ 6, quadNum = activeNum * (activeNum - 1) * (activeNum - 2) * (activeNum - 3)/ 24;
    llvm_vecsmall::SmallVector<array<int, 2>,40> pair(pairNum);
    llvm_vecsmall::SmallVector<array<int, 3>, 100> triple(triNum);
    llvm_vecsmall::SmallVector<array<int, 4>, 300> quad(quadNum);
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<array<int, 4>, 100>, 3> multiples = multiple_indices[activeNum - 1];
//    int pairIt = 0, triIt = 0, quadIt = 0;
//    for (int i = 0; i < activeNum - 1; i++){
//        for (int j = i + 1; j < activeNum; j++){
//            pair[pairIt] = {activeFunc[i], activeFunc[j]};
//            pairIt ++;
//            if (j < activeNum - 1){
//                for (int k = j + 1; k < activeNum; k++){
//                    triple[triIt] = {activeFunc[i], activeFunc[j], activeFunc[k]};
//                    triIt ++;
//                    if (k < activeNum - 1){
//                        for (int m = k + 1; m < activeNum; m++){
//                            quad[quadIt] = {activeFunc[i], activeFunc[j], activeFunc[k], activeFunc[m]};
//                            quadIt++;
//                        }
//                    }
//                }
//            }
//        }
//    }
//    get_func_timer.Stop();
    
    for (int pairIter = 0; pairIter < pairNum; pairIter ++){
        Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        array<int, 2> pairIndices = {multiples[0][pairIter][0],multiples[0][pairIter][1]};
        pair[pairIter] = {activeFunc[pairIndices[0]], activeFunc[pairIndices[1]]};
        array<double, 20> diff_at_point;
        int funcIndex1 = pair[pairIter][0];
        int funcIndex2 = pair[pairIter][1];
        for (int i = 0; i < 20; i++){
            diff_at_point[i] = valList[funcIndex1][i] - valList[funcIndex2][i];
        }
        bool activeTF = get_sign(*std::max_element(diff_at_point.begin(), diff_at_point.end())) == get_sign(*std::min_element(diff_at_point.begin(), diff_at_point.end())) ? false : true;
        single_timer.Stop();
        if (activeTF){
            if (!active){
                active = true;
            }
            activePair[pair[pairIter][0]][pair[pairIter][1]] = true;
            activePair[pair[pairIter][1]][pair[pairIter][0]] = true;
            if (!activeList[funcIndex1]){
                activeList[funcIndex1] = true;
                double d1 = valList[funcIndex1][1]-valList[funcIndex1][0], d2 = valList[funcIndex1][2]-valList[funcIndex1][0], d3 = valList[funcIndex1][3]-valList[funcIndex1][0];
                llvm_vecsmall::SmallVector<double, 20> unNormF = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({{d1, d2, d3}}), crossMatrix)[0];
                gradList[funcIndex1] = {unNormF[0],unNormF[1],unNormF[2]};
                
                for (int i = 0; i < 16; ++i) {
                    diffList[funcIndex1][i] = valList[funcIndex1][i + 4] -
                    (valList[funcIndex1][0] * coeff[i][0] + valList[funcIndex1][1] * coeff[i][1] + valList[funcIndex1][2] * coeff[i][2] + valList[funcIndex1][3] * coeff[i][3]) / 3.0;
                }
                //errorList[funcIndex1] = std::max(*max_element(diffList[funcIndex1].begin(), diffList[funcIndex1].end()), std::abs(*min_element(diffList[funcIndex1].begin(), diffList[funcIndex1].end())));
            }
            if (!activeList[funcIndex2]){
                activeList[funcIndex2] = true;
                double d1 = valList[funcIndex2][1]-valList[funcIndex2][0], d2 = valList[funcIndex2][2]-valList[funcIndex2][0], d3 = valList[funcIndex2][3]-valList[funcIndex2][0];
                llvm_vecsmall::SmallVector<double, 20> unNormF = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({{d1, d2, d3}}), crossMatrix)[0];
                gradList[funcIndex2] = {unNormF[0],unNormF[1],unNormF[2]};
                
                for (int i = 0; i < 16; ++i) {
                    diffList[funcIndex2][i] = valList[funcIndex2][i + 4] -
                    (valList[funcIndex2][0] * coeff[i][0] + valList[funcIndex2][1] * coeff[i][1] + valList[funcIndex2][2] * coeff[i][2] + valList[funcIndex2][3] * coeff[i][3]) / 3.0;
                }
                //errorList[funcIndex1] = std::max(*max_element(diffList[funcIndex2].begin(), diffList[funcIndex2].end()), std::abs(*min_element(diffList[funcIndex2].begin(), diffList[funcIndex2].end())));
            }
            array<double, 16> diff_twofunc;
            for (int i = 0; i < 16; ++i){
                diff_twofunc[i] = diffList[funcIndex1][i] - diffList[funcIndex2][i];
            }
            
            Timer single2_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            double error =std::max(*max_element(diff_twofunc.begin(), diff_twofunc.end()), std::abs(*min_element(diff_twofunc.begin(), diff_twofunc.end())));
            array<double, 3> grad;
            for (int i = 0; i < 3; ++i){
                grad[i] = gradList[funcIndex1][i] - gradList[funcIndex2][i];
            }
            double lhs = error * error * sqD;
            double rhs;
            if (!curve_network){
                rhs = threshold * threshold * dot(grad, grad);
            }else{
                rhs = numeric_limits<double>::infinity() * dot(grad, grad);
            }
            //cout << lhs << " " << rhs << endl;
            if (lhs > rhs) {
                single2_timer.Stop();
                return score;
            }
            single2_timer.Stop();
            
        }
    }
    
    // 2-function checks
    int activeTriple_count = 0;
    {
        Timer timer(twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int tripleIter = 0; tripleIter < triNum; tripleIter ++){
            array<int, 3> tripleIndices = {multiples[1][tripleIter][0], multiples[1][tripleIter][1], multiples[1][tripleIter][2]};
            triple[tripleIter] = {activeFunc[tripleIndices[0]], activeFunc[tripleIndices[1]], activeFunc[tripleIndices[2]]};
            int funcIndex1 = triple[tripleIter][0];
            int funcIndex2 = triple[tripleIter][1];
            int funcIndex3 = triple[tripleIter][2];
            if(!(activePair[funcIndex1][funcIndex2]&&activePair[funcIndex1][funcIndex3]&&activePair[funcIndex2][funcIndex3]))
                continue;
            array<double, 20> diffList1, diffList2;
            for (int i = 0; i < 20; ++i){
                diffList1[i] = valList[funcIndex1][i] - valList[funcIndex2][i];
                diffList2[i] = valList[funcIndex2][i] - valList[funcIndex3][i];
            }
            std::array<double, 40> nPoints = transpose2d({diffList1, diffList2});// X0, Y0, X1, Y1, ...
            std::array<double, 2> query = {0.0, 0.0}; // X, Y
            sub_call_two ++;
            Timer sub_timer(sub_twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<2, double>(nPoints, query);
            sub_timer.Stop();
            if (zeroX){
                activeTriple_count++;
                array<double, 16> diff_twofunc1, diff_twofunc2;
                for (int i = 0; i < 16; ++i){
                    diff_twofunc1[i] = diffList[funcIndex1][i] - diffList[funcIndex2][i];
                    diff_twofunc2[i] = diffList[funcIndex2][i] - diffList[funcIndex3][i];
                }
                array<double, 3> grad1, grad2;
                for (int i = 0; i < 3; ++i){
                    grad1[i] = gradList[funcIndex1][i] - gradList[funcIndex2][i];
                    grad2[i] = gradList[funcIndex2][i] - gradList[funcIndex3][i];
                }
                
                // two function linearity test:
                std::array<double, 2> w1 = {dot(grad1, grad1), dot(grad1, grad2)};
                std::array<double, 2> w2 = {dot(grad1, grad2), dot(grad2, grad2)};
                double E = det(w1, w2);
                std::array<double, 2>invPerp_w2 = perp(w2);
                for(int i = 0; i < 2; i++){
                    invPerp_w2[i] *= -1;
                }
                std::array<double, 2> perp_w1 = perp(w1);
                llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> H = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({llvm_vecsmall::SmallVector<double, 20>(std::begin(invPerp_w2), std::end(invPerp_w2)), llvm_vecsmall::SmallVector<double, 20>(std::begin(perp_w1), std::end(perp_w1))}), llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({llvm_vecsmall::SmallVector<double, 20>(std::begin(grad1), std::end(grad1)), llvm_vecsmall::SmallVector<double, 20>(std::begin(grad2), std::end(grad2))}));
                
                //find the largest max error (max squared gamma: the LHS of the equation) among all 16 bezier control points
                double maxGammaSq = 0;
                for (int i = 0; i < 16; i++){
                    //std::cout << diffList[pair[pairIter][0]][i] << " " << diffList[pair[pairIter][1]][i] << std::endl;
                    llvm_vecsmall::SmallVector<double, 20> unNormDis = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({{diff_twofunc1[i], diff_twofunc2[i]}}), llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>(H))[0];
                    double currError = sqD * dot(unNormDis, unNormDis);
                    if (maxGammaSq < currError)
                        maxGammaSq = currError;
                }
                if (maxGammaSq > threshold*threshold * E*E){
                    timer.Stop();
                    return score;
                }
            }
        }
        timer.Stop();
    }
    if(activeTriple_count < 4)
        return false;
    {
        Timer timer(threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int quadIter = 0; quadIter < quadNum; quadIter ++){
            array<int, 4> quadIndices = {multiples[2][quadIter][0], multiples[2][quadIter][1], multiples[2][quadIter][2],multiples[2][quadIter][3]};
            quad[quadIter] = {activeFunc[quadIndices[0]], activeFunc[quadIndices[1]], activeFunc[quadIndices[2]],activeFunc[quadIndices[3]]};
            int funcIndex1 = quad[quadIter][0];
            int funcIndex2 = quad[quadIter][1];
            int funcIndex3 = quad[quadIter][2];
            int funcIndex4 = quad[quadIter][3];
            if(!(activePair[funcIndex1][funcIndex2]&&activePair[funcIndex1][funcIndex3]&&activePair[funcIndex1][funcIndex4]&&activePair[funcIndex2][funcIndex3]&&activePair[funcIndex2][funcIndex4]&&activePair[funcIndex3][funcIndex4]))
                continue;
            array<double, 20> diffList1, diffList2, diffList3;
            for (int i = 0; i < 20; ++i){
                diffList1[i] = valList[funcIndex1][i] - valList[funcIndex2][i];
                diffList2[i] = valList[funcIndex2][i] - valList[funcIndex3][i];
                diffList3[i] = valList[funcIndex3][i] - valList[funcIndex4][i];
            }
            std::array<double, 60> nPoints = transpose3d({diffList1, diffList2, diffList3});
            std::array<double, 3> query = {0.0, 0.0, 0.0}; // X, Y, Z
            sub_call_three ++;
            Timer sub_timer(sub_threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<3, double>(nPoints, query);
            sub_timer.Stop();
            if (zeroX){
                array<double, 16> diff_twofunc1, diff_twofunc2, diff_twofunc3;
                for (int i = 0; i < 16; ++i){
                    diff_twofunc1[i] = diffList[funcIndex1][i] - diffList[funcIndex2][i];
                    diff_twofunc2[i] = diffList[funcIndex2][i] - diffList[funcIndex3][i];
                    diff_twofunc3[i] = diffList[funcIndex3][i] - diffList[funcIndex4][i];
                }
                array<double, 3> fi, fj, fk;
                for (int i = 0; i < 3; ++i){
                    fi[i] = gradList[funcIndex1][i] - gradList[funcIndex2][i];
                    fj[i] = gradList[funcIndex2][i] - gradList[funcIndex3][i];
                    fk[i] = gradList[funcIndex3][i] - gradList[funcIndex4][i];
                }
                double E = det(fi, fj, fk);
                llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> H = {cross(fj, fk), cross(fk, fi), cross(fi, fj)};
                double maxGammaSq = 0;
                for (int i = 0; i < 16; i++){
                    llvm_vecsmall::SmallVector<double, 20> unNormDis = matrixMultiply(llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20>({{diff_twofunc1[i], diff_twofunc2[i], diff_twofunc3[i]}}), H)[0];
                    double currError = sqD * dot(unNormDis, unNormDis);
                    if (maxGammaSq < currError)
                        maxGammaSq = currError;
                }
                if (maxGammaSq > threshold*threshold * E*E){
                    timer.Stop();
                    return score;
                }
            }
        }
        timer.Stop();
    }
    return false;
}



#endif /* subdivide_multi_h */
