//
//  subdivide_multi.cpp
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 6/20/24.
//

#include "refine_crit.h"

//using namespace std;

bool curve_network = false;
//llvm_vecsmall::SmallVector<csg_unit, 20> GLOBAL_CSGTREE = {};

/// Below are the variable/constant/function definitions that will only be used in `subdivide_multi.cpp`
///
///
/// Stores the index of permutations of n less than `funcNum` for the use in `subTet` and `subMI` functions.
///
llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<std::array<int, 4>, 100>, 3>, 20> multiple_indices;

std::array<double, 2> query_2d = {0.0, 0.0}; // X, Y
std::array<double, 3> query_3d = {0.0, 0.0, 0.0}; // X, Y, Z

/// Below are the local functions servicing `subTet` and `subMI`

/// returns a `bool` value that `true` represents positive and `false` represents negative of the input value `x`.
bool get_sign(double x) {
    return (x > 0) ? 1 : 0;
}

/// returns the dot product of input arrays `a` and `b`.
double dot(const std::array<double, 3> &a, const std::array<double, 3> &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/// returns a vector in space by subtracting a 3D coordinate `p2` from `p1`.
std::array<double, 3> getVec(const std::array<double, 3> &p1, const std::array<double, 3> &p2){
    return {p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]};
}

/// returns a vector in space by subtracting a 3D coordinate `p2` from `p1`. Output is in Eigen vector.
Eigen::Vector3d getEigenVec(const std::array<double, 3> &p1, const std::array<double, 3> &p2){
    return Eigen::Vector3d(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
}

/// transforms the input of errors at 20 bezier control points for two functions into the correct format that `convex_hull_membership` library can use.
std::array<double, 40> parse_convex_points2d(const Eigen::Matrix<double, 2, 20> valList) {
    std::array<double, 40> transposed;
    Eigen::MatrixXd::Map(transposed.data(), 1, 40) = valList.transpose();
    return transposed;
}

/// transforms the input of errors at 20 bezier control points for three functions into the correct format that `convex_hull_membership` library can use.
std::array<double, 60> parse_convex_points3d(const Eigen::Matrix<double, 3, 20> valList) {
    std::array<double, 60> transposed;
    Eigen::MatrixXd::Map(transposed.data(), 1, 60) = valList.transpose();
    return transposed;
}

Eigen::Vector<double, 20> bezierConstruct(const double v0,
                                       const double v1,
                                       const double v2,
                                       const double v3,
                                       const std::array<double, 3> g0,
                                       const std::array<double, 3> g1,
                                       const std::array<double, 3> g2,
                                       const std::array<double, 3> g3,
                                       const std::array<double, 3> vec1,
                                       const std::array<double, 3> vec2,
                                       const std::array<double, 3> vec3,
                                       const std::array<double, 3> vec4,
                                       const std::array<double, 3> vec5,
                                       const std::array<double, 3> vec6) {
    std::array<double, 3> v0s, v1s, v2s, v3s;
    v0s = {v0 + dot(g0, vec1) / 3, v0 + dot(g0, vec2) / 3, v0 + dot(g0, vec3) / 3};
    v1s = {v1 + dot(g1, vec4) / 3, v1 + dot(g1, vec5) / 3, v1 - dot(g1, vec1) / 3};
    v2s = {v2 + dot(g2, vec6) / 3, v2 - dot(g2, vec2) / 3, v2 - dot(g2, vec4) / 3};
    v3s = {v3 - dot(g3, vec3) / 3, v3 - dot(g3, vec5) / 3, v3 - dot(g3, vec6) / 3};
    double vMid0 = (9 * (v1s[0] + v1s[1] + v2s[0] + v2s[2] + v3s[1] + v3s[2]) / 6 - v1 - v2 - v3)/ 6;
    //double e1 = (v0s[1] + v0s[2] + v2s[0] + v2s[1] + v3s[0] + v3s[2]) / 6;
    double vMid1 =(9 * (v0s[1] + v0s[2] + v2s[0] + v2s[1] + v3s[0] + v3s[2]) / 6 - v0 - v2 - v3)/ 6;
    //double e2 = (v0s[0] + v0s[2] + v1s[1] + v1s[2] + v3s[0] + v3s[1]) / 6;
    double vMid2 =(9 * (v0s[0] + v0s[2] + v1s[1] + v1s[2] + v3s[0] + v3s[1]) / 6 - v0 - v1 - v3)/ 6;
    //double e3 = (v0s[0] + v0s[1] + v1s[0] + v1s[2] + v2s[1] + v2s[2]) / 6;
    double vMid3 =(9 * (v0s[0] + v0s[1] + v1s[0] + v1s[2] + v2s[1] + v2s[2]) / 6 - v0 - v1 - v2)/ 6;
    return Eigen::Vector<double, 20> {v0, v1, v2, v3, v0s[0], v0s[1], v0s[2], v1s[0], v1s[1], v1s[2], v2s[0], v2s[1], v2s[2],
        v3s[0], v3s[1], v3s[2], vMid0, vMid1, vMid2, vMid3};
}

Eigen::Vector<double, 16> bezierDiff(const Eigen::Vector<double,20> valList)
{
    /// Constant coefficient to obtain linear interpolated values at each bezier control points
    const Eigen::Matrix<double, 16, 4> linear_coeff {{2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1}, {0, 2, 1, 0},{0, 2, 0, 1}, {1, 2, 0, 0}, {0, 0, 2, 1}, {1, 0, 2, 0},{0, 1, 2, 0}, {1, 0, 0, 2}, {0, 1, 0, 2}, {0, 0, 1, 2},{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}};
    Eigen::Vector<double, 16> linear_val = (linear_coeff * valList.head(4)) / 3;
    return valList.tail(16) - linear_val;
}

void init_multi(const size_t funcNum,
                const int mode)
{
    multiple_indices.resize(funcNum);
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        multiple_indices[funcIter].resize(3);
        int activeNum = funcIter + 1;
        int pairNum = activeNum * (activeNum-1)/2, triNum = activeNum * (activeNum-1) * (activeNum - 2)/ 6;
        int quadNum = activeNum * (activeNum - 1) * (activeNum - 2) * (activeNum - 3)/ 24;
        llvm_vecsmall::SmallVector<std::array<int, 4>,100> pair(pairNum);
        llvm_vecsmall::SmallVector<std::array<int, 4>, 100> triple(triNum);
        llvm_vecsmall::SmallVector<std::array<int, 4>, 100> quad(quadNum);
        int pairIt = 0, triIt = 0, quadIt = 0;
        for (int i = 0; i < activeNum - 1; i++){
            for (int j = i + 1; j < activeNum; j++){
                pair[pairIt] = {i, j, 0, 0};
                pairIt ++;
                if (j < activeNum - 1){
                    for (int k = j + 1; k < activeNum; k++){
                        triple[triIt] = {i, j, k, 0};
                        triIt ++;
                        if (mode == MI){
                            if (k < activeNum - 1){
                                for (int m = k + 1; m < activeNum; m++){
                                    quad[quadIt] = {i, j, k, m};
                                    quadIt++;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (mode == MI){
            multiple_indices[funcIter] = {pair, triple, quad};
        }else{
            multiple_indices[funcIter] = {pair, triple};
        }
    }
}

bool two_func_check (Eigen::Matrix<double, 2, 3> grad,
                     const Eigen::Matrix<double, 16, 2> grad_matrix,
                     const double sqD,
                     const double threshold)
{
    Eigen::Matrix2d w;
    w << grad.row(0).squaredNorm(), grad.row(0).dot(grad.row(1)),
    grad.row(0).dot(grad.row(1)), grad.row(1).squaredNorm();
    double E = w.determinant();
    Eigen::Matrix<double, 2, 3> H = Eigen::Matrix2d({{w(1, 1), -w(1, 0)}, {-w(0, 1), w(0,0)}}) * grad;
    
    //find the largest max error (max squared gamma: the LHS of the equation) among all 16 bezier control points
    Eigen::Matrix<double, 16, 3> unNormDis = grad_matrix * H;
    Eigen::Vector<double, 16> dotProducts = sqD * unNormDis.cwiseProduct(unNormDis).rowwise().sum();
    return (dotProducts.maxCoeff() > threshold*threshold * E * E);
}

bool three_func_check (Eigen::Matrix<double, 3, 3> grad,
                     const Eigen::Matrix<double, 16, 3> grad_matrix,
                     const double sqD,
                     const double threshold)
{
    double E = grad.determinant();
    Eigen::Matrix<double, 3, 3> H;
    H << grad.row(1).cross(grad.row(2)),
    grad.row(2).cross(grad.row(0)),
    grad.row(0).cross(grad.row(1));
    Eigen::Matrix<double, 16, 3> unNormDis_eigen = grad_matrix * H;
    Eigen::Vector<double, 16> dotProducts = sqD * unNormDis_eigen.cwiseProduct(unNormDis_eigen).rowwise().sum();
    //double maxGammaSq = dotProducts.maxCoeff();
    return (dotProducts.maxCoeff() > threshold*threshold * E * E);
}

bool critIA(
            const std::array<std::array<double, 3>,4> &pts,
            const llvm_vecsmall::SmallVector<std::array<double,4>, 20> &vals,
            const llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> &grads,
            const double threshold,
            bool& active,
            int &sub_call_two,
            int &sub_call_three)
{
    std::array<double, 3> p0 = pts[0], p1 = pts[1], p2 = pts[2], p3 = pts[3];
    std::array<double, 3> vec1 = getVec(p1, p0), vec2 = getVec(p2, p0), vec3 = getVec(p3, p0),
    vec4 = getVec(p2, p1), vec5 = getVec(p3, p1), vec6 = getVec(p3, p2);
    const size_t funcNum = vals.size();
    
    Eigen::Matrix<double, Eigen::Dynamic, 20> valList (funcNum, 20);
    Eigen::Matrix<double, Eigen::Dynamic, 16> diffList(funcNum, 16);
    llvm_vecsmall::SmallVector<bool, 20> activeTF;
    Eigen::Matrix<double, 20, 3> gradList;
    Eigen::Vector3d eigenVec1 = getEigenVec(p1, p0), eigenVec2 = getEigenVec(p2, p0), eigenVec3 = getEigenVec(p3, p0);
    Eigen::Matrix3d vec;
    vec << eigenVec1, eigenVec2, eigenVec3;
    double D = vec.determinant();
    double sqD = D*D;
    Eigen::Matrix3d crossMatrix;
    crossMatrix << eigenVec2.cross(eigenVec3), eigenVec3.cross(eigenVec1), eigenVec1.cross(eigenVec2);
    
    int activeNum = 0;
    //single function linearity check:
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        //Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        double v0 = vals[funcIter][0], v1 = vals[funcIter][1], v2 = vals[funcIter][2], v3 = vals[funcIter][3];
        std::array<double, 3> g0 = grads[funcIter][0], g1 = grads[funcIter][1], g2 = grads[funcIter][2], g3 = grads[funcIter][3];
       //storing bezier and linear info for later linearity comparison
        valList.row(funcIter) = bezierConstruct(v0, v1, v2, v3, g0, g1, g2, g3, vec1, vec2, vec3, vec4, vec5, vec6);
        //Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        activeTF[funcIter] = get_sign(valList.row(funcIter).maxCoeff()) != get_sign(valList.row(funcIter).minCoeff());
        //single_timer.Stop();
        if (activeTF[funcIter]){
            if (!active){
                active = true;
            }
            activeNum++;
            Eigen::Vector3d unNormF = Eigen::RowVector3d(v1-v0, v2-v0, v3-v0) * crossMatrix.transpose();
            gradList.row(funcIter) = unNormF;
            diffList.row(funcIter) = bezierDiff(valList.row(funcIter));
            double error = std::max(diffList.row(funcIter).maxCoeff(), -diffList.row(funcIter).minCoeff());
            //Timer single2_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            double lhs = error * error * sqD;
            double rhs;
            if (!curve_network){
                rhs = threshold * threshold * gradList.row(funcIter).squaredNorm();
            }else{
                rhs = std::numeric_limits<double>::infinity() * gradList.row(funcIter).squaredNorm();
            }
            if (lhs > rhs) {
                //single2_timer.Stop();
                return true;
            }
            //single2_timer.Stop();
        }
    }
    //Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    if(activeNum < 2){
        //single_timer.Stop();
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
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<bool, 20>, 20> zeroXResult(funcNum, llvm_vecsmall::SmallVector<bool, 20>(funcNum));
    //single_timer.Stop();
    const int pairNum = activeNum * (activeNum-1)/2, triNum = activeNum * (activeNum-1) * (activeNum - 2)/ 6;
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<std::array<int, 4>, 100>, 3> multiples = multiple_indices[activeNum - 1];
    
    // 2-function checks
    int activeDouble_count = 0;
    {
        //Timer timer(twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int pairIter = 0; pairIter < pairNum; pairIter ++){
            std::array<int, 2> pairIndices = {activeFunc[multiples[0][pairIter][0]],activeFunc[multiples[0][pairIter][1]]};
            std::array<double, 40> nPoints = parse_convex_points2d(valList(pairIndices, Eigen::all));
            //Timer sub_timer(sub_twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<2, double>(nPoints, query_2d);
            //sub_timer.Stop();
            
            if (zeroX){
                activeDouble_count++;
                sub_call_two ++;
                zeroXResult[pairIndices[0]][pairIndices[1]] = true;
                zeroXResult[pairIndices[1]][pairIndices[0]] = true;
                Eigen::Matrix<double, 2, 3> grad = gradList(pairIndices, Eigen::all);
                Eigen::Matrix<double, 16, 2> grad_matrix = diffList(pairIndices, Eigen::all).transpose();
                // two function linearity test:
                if (two_func_check (grad, grad_matrix, sqD, threshold)){
                    //timer.Stop();
                    return true;
                }
            }
        }
        //timer.Stop();
    }
    if(activeDouble_count < 3)
        return false;
    {
        //Timer timer(threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int triIter = 0; triIter < triNum; triIter ++){
            std::array<int, 3> tripleIndices = {activeFunc[multiples[1][triIter][0]], activeFunc[multiples[1][triIter][1]], activeFunc[multiples[1][triIter][2]]};
            if(!(zeroXResult[tripleIndices[0]][tripleIndices[1]]&&zeroXResult[tripleIndices[0]][tripleIndices[2]]&&zeroXResult[tripleIndices[1]][tripleIndices[2]]))
                continue;
            std::array<double, 60> nPoints = parse_convex_points3d(valList(tripleIndices, Eigen::all));
            sub_call_three ++;
            //Timer sub_timer(sub_threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<3, double>(nPoints, query_3d);
            //sub_timer.Stop();

            if (zeroX){
                Eigen::Matrix<double, 3, 3> grad = gradList(tripleIndices, Eigen::all);
                Eigen::Matrix<double, 16, 3> grad_matrix = diffList(tripleIndices, Eigen::all).transpose();
                if (three_func_check (grad, grad_matrix, sqD, threshold)){
                    //timer.Stop();
                    return true;
                }
            }
        }
        //timer.Stop();
    }
    return false;
}

bool critCSG(
            const std::array<std::array<double, 3>,4> &pts,
            const llvm_vecsmall::SmallVector<std::array<double,4>, 20> &vals,
            const llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> &grads,const std::function<std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>>(llvm_vecsmall::SmallVector<std::array<double, 2>, 20>)> csg_func,
            const double threshold,
            bool& active,
            int &sub_call_two,
            int &sub_call_three)
{
    std::array<double, 3> p0 = pts[0], p1 = pts[1], p2 = pts[2], p3 = pts[3];
    std::array<double, 3> vec1 = getVec(p1, p0), vec2 = getVec(p2, p0), vec3 = getVec(p3, p0),
    vec4 = getVec(p2, p1), vec5 = getVec(p3, p1), vec6 = getVec(p3, p2);
    const size_t funcNum = vals.size();
    
    Eigen::Matrix<double, Eigen::Dynamic, 20> valList (funcNum, 20);
    Eigen::Matrix<double, Eigen::Dynamic, 16> diffList(funcNum, 16);
    llvm_vecsmall::SmallVector<bool, 20> activeTF;
    llvm_vecsmall::SmallVector<std::array<double , 2>, 20> funcInt(funcNum);
    Eigen::Matrix<double, 20, 3> gradList;
    Eigen::Vector3d eigenVec1 = getEigenVec(p1, p0), eigenVec2 = getEigenVec(p2, p0), eigenVec3 = getEigenVec(p3, p0);
    Eigen::Matrix3d vec;
    vec << eigenVec1, eigenVec2, eigenVec3;
    double D = vec.determinant();
    double sqD = D*D;
    Eigen::Matrix3d crossMatrix;
    crossMatrix << eigenVec2.cross(eigenVec3), eigenVec3.cross(eigenVec1), eigenVec1.cross(eigenVec2);
    
    int activeNum = 0;
    //single function linearity check:
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        //Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        double v0 = vals[funcIter][0], v1 = vals[funcIter][1], v2 = vals[funcIter][2], v3 = vals[funcIter][3];
        std::array<double, 3> g0 = grads[funcIter][0], g1 = grads[funcIter][1], g2 = grads[funcIter][2], g3 = grads[funcIter][3];
       //storing bezier and linear info for later linearity comparison
        valList.row(funcIter) = bezierConstruct(v0, v1, v2, v3, g0, g1, g2, g3, vec1, vec2, vec3, vec4, vec5, vec6);
        funcInt[funcIter] = {valList.row(funcIter).minCoeff(), valList.row(funcIter).maxCoeff()};
    }
    
        std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> csgResult = csg_func(funcInt);
        if(csgResult.first[0] * csgResult.first[1] > 0){
            return false;
        }else{
            for (size_t funcIter = 0; funcIter < funcNum; funcIter++){
                //Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
                activeTF[funcIter] = !csgResult.second[funcIter];
                //single_timer.Stop();
                if (activeTF[funcIter]){
                    if (!active){
                        active = true;
                    }
                    activeNum++;
                    double v0 = vals[funcIter][0], v1 = vals[funcIter][1], v2 = vals[funcIter][2], v3 = vals[funcIter][3];
                    Eigen::Vector3d unNormF = Eigen::RowVector3d(v1-v0, v2-v0, v3-v0) * crossMatrix.transpose();
                    gradList.row(funcIter) = unNormF;
                    diffList.row(funcIter) = bezierDiff(valList.row(funcIter));
                    double error = std::max(diffList.row(funcIter).maxCoeff(), -diffList.row(funcIter).minCoeff());
                    //Timer single2_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
                    double lhs = error * error * sqD;
                    double rhs;
                    if (!curve_network){
                        rhs = threshold * threshold * gradList.row(funcIter).squaredNorm();
                    }else{
                        rhs = std::numeric_limits<double>::infinity() * gradList.row(funcIter).squaredNorm();
                    }
                    if (lhs > rhs) {
                        //single2_timer.Stop();
                        return true;
                    }
                    //single2_timer.Stop();
                }
            }
        }
    //Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    if(activeNum < 2){
        //single_timer.Stop();
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
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<bool, 20>, 20> zeroXResult(funcNum, llvm_vecsmall::SmallVector<bool, 20>(funcNum));
    //single_timer.Stop();
    const int pairNum = activeNum * (activeNum-1)/2, triNum = activeNum * (activeNum-1) * (activeNum - 2)/ 6;
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<std::array<int, 4>, 100>, 3> multiples = multiple_indices[activeNum - 1];
    
    // 2-function checks
    int activeDouble_count = 0;
    {
        //Timer timer(twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int pairIter = 0; pairIter < pairNum; pairIter ++){
            std::array<int, 2> pairIndices = {activeFunc[multiples[0][pairIter][0]],activeFunc[multiples[0][pairIter][1]]};
            std::array<double, 40> nPoints = parse_convex_points2d(valList(pairIndices, Eigen::all));
            //Timer sub_timer(sub_twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<2, double>(nPoints, query_2d);
            //sub_timer.Stop();
            
            if (zeroX){
                activeDouble_count++;
                sub_call_two ++;
                zeroXResult[pairIndices[0]][pairIndices[1]] = true;
                zeroXResult[pairIndices[1]][pairIndices[0]] = true;
                Eigen::Matrix<double, 2, 3> grad = gradList(pairIndices, Eigen::all);
                Eigen::Matrix<double, 16, 2> grad_matrix = diffList(pairIndices, Eigen::all).transpose();
                // two function linearity test:
                if (two_func_check (grad, grad_matrix, sqD, threshold)){
                    //timer.Stop();
                    return true;
                }
            }
        }
        //timer.Stop();
    }
    if(activeDouble_count < 3)
        return false;
    {
        //Timer timer(threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int triIter = 0; triIter < triNum; triIter ++){
            std::array<int, 3> tripleIndices = {activeFunc[multiples[1][triIter][0]], activeFunc[multiples[1][triIter][1]], activeFunc[multiples[1][triIter][2]]};
            if(!(zeroXResult[tripleIndices[0]][tripleIndices[1]]&&zeroXResult[tripleIndices[0]][tripleIndices[2]]&&zeroXResult[tripleIndices[1]][tripleIndices[2]]))
                continue;
            std::array<double, 60> nPoints = parse_convex_points3d(valList(tripleIndices, Eigen::all));
            sub_call_three ++;
            //Timer sub_timer(sub_threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<3, double>(nPoints, query_3d);
            //sub_timer.Stop();

            if (zeroX){
                Eigen::Matrix<double, 3, 3> grad = gradList(tripleIndices, Eigen::all);
                Eigen::Matrix<double, 16, 3> grad_matrix = diffList(tripleIndices, Eigen::all).transpose();
                if (three_func_check (grad, grad_matrix, sqD, threshold)){
                    //timer.Stop();
                    return true;
                }
            }
        }
        //timer.Stop();
    }
    return false;
}
bool critMI(
           const std::array<std::array<double, 3>,4> &pts,
           const llvm_vecsmall::SmallVector<std::array<double,4>, 20> &vals,
           const llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> &grads, 
           const double threshold,
           bool& active,
           int &sub_call_two,
           int &sub_call_three)
{
    std::array<double, 3> p0 = pts[0], p1 = pts[1], p2 = pts[2], p3 = pts[3];
    std::array<double, 3> vec1 = getVec(p1, p0), vec2 = getVec(p2, p0), vec3 = getVec(p3, p0),
    vec4 = getVec(p2, p1), vec5 = getVec(p3, p1), vec6 = getVec(p3, p2);
    
    
    Eigen::Vector3d eigenVec1 = getEigenVec(p1, p0), eigenVec2 = getEigenVec(p2, p0), eigenVec3 = getEigenVec(p3, p0);
    Eigen::Matrix3d vec;
    vec << eigenVec1, eigenVec2, eigenVec3;
    double D = vec.determinant();
    double sqD = D*D;
    const size_t funcNum = vals.size();
    Eigen::Matrix3d crossMatrix_eigen;
    crossMatrix_eigen << eigenVec2.cross(eigenVec3), eigenVec3.cross(eigenVec1), eigenVec1.cross(eigenVec2);
    Eigen::Matrix<double, 20, 3> gradList_eigen;
    Eigen::Matrix<double, Eigen::Dynamic, 20> valList_eigen (funcNum, 20);
    Eigen::Matrix<double, Eigen::Dynamic, 16> diffList(funcNum, 16);
    
    llvm_vecsmall::SmallVector<bool, 20> activeList(funcNum);
    llvm_vecsmall::SmallVector<std::array<double , 2>, 20> funcInt(funcNum);
    double maxLow = -1 * std::numeric_limits<double>::infinity();
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<bool, 20>, 20> activePair(funcNum, llvm_vecsmall::SmallVector<bool, 20>(false, funcNum));
    //single function linearity check:
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        double v0 = vals[funcIter][0], v1 = vals[funcIter][1], v2 = vals[funcIter][2], v3 = vals[funcIter][3];
        std::array<double, 3> g0 = grads[funcIter][0], g1 = grads[funcIter][1], g2 = grads[funcIter][2], g3 = grads[funcIter][3];
        // Bezier control points
        valList_eigen.row(funcIter) = bezierConstruct(v0, v1, v2, v3, g0, g1, g2, g3, vec1, vec2, vec3, vec4, vec5, vec6);
        funcInt[funcIter] = {valList_eigen.row(funcIter).minCoeff(), valList_eigen.row(funcIter).maxCoeff()};
        //cout << funcInt[funcIter][0] << " " << funcInt[funcIter][1] << endl;
        if (maxLow < funcInt[funcIter][0]){
            maxLow = funcInt[funcIter][0];
        }
    }
    llvm_vecsmall::SmallVector<int, 20> activeFunc;
    for (int funcIter = 0; funcIter < funcNum; funcIter++){
        if(funcInt[funcIter][1] > maxLow){
            activeFunc.push_back(funcIter);
        }
    }
    size_t activeNum = activeFunc.size();
    if(activeNum < 2)
        return false;

    //Timer get_func_timer(getActiveMuti, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    const size_t pairNum = activeNum * (activeNum-1)/2, triNum = activeNum * (activeNum-1) * (activeNum - 2)/ 6, quadNum = activeNum * (activeNum - 1) * (activeNum - 2) * (activeNum - 3)/ 24;
    llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<std::array<int, 4>, 100>, 3> multiples = multiple_indices[activeNum - 1];
//    get_func_timer.Stop();
    
    for (int pairIter = 0; pairIter < pairNum; pairIter ++){
        //Timer single_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        std::array<int, 2> pairIndices = {activeFunc[multiples[0][pairIter][0]],activeFunc[multiples[0][pairIter][1]]};
        int funcIndex1 = pairIndices[0];
        int funcIndex2 = pairIndices[1];
        Eigen::Vector<double, 20> diff_at_point;
        diff_at_point = valList_eigen.row(funcIndex2) - valList_eigen.row(funcIndex1);
        bool activeTF = get_sign(diff_at_point.maxCoeff()) == get_sign(diff_at_point.minCoeff()) ? false : true;
        //single_timer.Stop();
        if (activeTF){
            if (!active){
                active = true;
            }
            activePair[pairIndices[0]][pairIndices[1]] = true;
            activePair[pairIndices[1]][pairIndices[0]] = true;
            if (!activeList[funcIndex1]){
                activeList[funcIndex1] = true;
                double v0 = valList_eigen(funcIndex1, 0), v1 = valList_eigen(funcIndex1, 1), v2 = valList_eigen(funcIndex1, 2), v3 = valList_eigen(funcIndex1, 3);
                Eigen::Vector3d unNormF_eigen = Eigen::RowVector3d(v1-v0, v2-v0, v3-v0) * crossMatrix_eigen.transpose();
                gradList_eigen.row(funcIndex1) = unNormF_eigen;
                
                diffList.row(funcIndex1) = bezierDiff(valList_eigen.row(funcIndex1));
            }
            if (!activeList[funcIndex2]){
                activeList[funcIndex2] = true;
                double v0 = valList_eigen(funcIndex2, 0), v1 = valList_eigen(funcIndex2, 1), v2 = valList_eigen(funcIndex2, 2), v3 = valList_eigen(funcIndex2, 3);
                Eigen::Vector3d unNormF_eigen = Eigen::RowVector3d(v1-v0, v2-v0, v3-v0) * crossMatrix_eigen.transpose();
                gradList_eigen.row(funcIndex2) = unNormF_eigen;
                diffList.row(funcIndex2) = bezierDiff(valList_eigen.row(funcIndex2));

            }
            Eigen::Vector<double, 16> diff_twofunc;
            diff_twofunc = diffList.row(funcIndex1) - diffList.row(funcIndex2);
            //Timer single2_timer(singleFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            double error = std::max(diff_twofunc.maxCoeff(), -diff_twofunc.minCoeff());
            Eigen::Vector3d grad_eigen;
            grad_eigen = gradList_eigen.row(funcIndex1) - gradList_eigen.row(funcIndex2);
            double lhs = error * error * sqD;
            double rhs;
            if (!curve_network){
                rhs = threshold * threshold * grad_eigen.squaredNorm();
            }else{
                rhs = std::numeric_limits<double>::infinity() * grad_eigen.squaredNorm();
            }
            if (lhs > rhs) {
                //single2_timer.Stop();
                return true;
            }
            //single2_timer.Stop();
            
        }
    }
    
    // 2-function checks
    int activeTriple_count = 0;
    {
        //Timer timer(twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int tripleIter = 0; tripleIter < triNum; tripleIter ++){
            std::array<int, 3> tripleIndices = {multiples[1][tripleIter][0], multiples[1][tripleIter][1], multiples[1][tripleIter][2]};
            int funcIndex1 = /*triple[tripleIter][0]*/activeFunc[multiples[1][tripleIter][0]];
            int funcIndex2 = activeFunc[multiples[1][tripleIter][1]];
            int funcIndex3 = activeFunc[multiples[1][tripleIter][2]];
            if(!(activePair[funcIndex1][funcIndex2]&&activePair[funcIndex1][funcIndex3]&&activePair[funcIndex2][funcIndex3]))
                continue;
            Eigen::Matrix<double,Eigen::Dynamic, 20> diff_mi(2, 20);
            diff_mi.row(0) = valList_eigen.row(funcIndex1) - valList_eigen.row(funcIndex2);
            diff_mi.row(1) =  valList_eigen.row(funcIndex2) - valList_eigen.row(funcIndex3);
            std::array<double, 40> nPoints = parse_convex_points2d(diff_mi);
            //Timer sub_timer(sub_twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<2, double>(nPoints, query_2d);
            //sub_timer.Stop();
            if (zeroX){
                sub_call_two ++;
                activeTriple_count++;
                Eigen::Matrix<double, 2, 3> grad(2, 3);
                grad.row(0) = gradList_eigen.row(funcIndex1) - gradList_eigen.row(funcIndex2);
                grad.row(1) = gradList_eigen.row(funcIndex2) - gradList_eigen.row(funcIndex3);
                Eigen::Matrix<double, 2, 16> grad_matrix(2, 16);
                grad_matrix.row(0) = diffList.row(funcIndex1) - diffList.row(funcIndex2);
                grad_matrix.row(1) = diffList.row(funcIndex2) - diffList.row(funcIndex3);
                if (two_func_check (grad, grad_matrix.transpose(), sqD, threshold)){
                    //timer.Stop();
                    return true;
                }
            }
        }
        //timer.Stop();
    }
    if(activeTriple_count < 4)
        return false;
    {
        //Timer timer(threeFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        bool zeroX;
        for (int quadIter = 0; quadIter < quadNum; quadIter ++){
            std::array<int, 4> quadIndices = {multiples[2][quadIter][0], multiples[2][quadIter][1], multiples[2][quadIter][2],multiples[2][quadIter][3]};
            int funcIndex1 = activeFunc[multiples[2][quadIter][0]];
            int funcIndex2 = activeFunc[multiples[2][quadIter][1]];
            int funcIndex3 = activeFunc[multiples[2][quadIter][2]];
            int funcIndex4 = activeFunc[multiples[2][quadIter][3]];
            if(!(activePair[funcIndex1][funcIndex2]&&activePair[funcIndex1][funcIndex3]&&activePair[funcIndex1][funcIndex4]&&activePair[funcIndex2][funcIndex3]&&activePair[funcIndex2][funcIndex4]&&activePair[funcIndex3][funcIndex4]))
                continue;
            
            Eigen::Matrix<double,Eigen::Dynamic, 20> diff_mi(3, 20);
            diff_mi.row(0) = valList_eigen.row(funcIndex1) - valList_eigen.row(funcIndex2);
            diff_mi.row(1) =  valList_eigen.row(funcIndex2) - valList_eigen.row(funcIndex3);
            diff_mi.row(2) =  valList_eigen.row(funcIndex3) - valList_eigen.row(funcIndex4);
            std::array<double, 60> nPoints = parse_convex_points3d(diff_mi);
            //Timer sub_timer(sub_twoFunc, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            zeroX = convex_hull_membership::contains<3, double>(nPoints, query_3d);
            //sub_timer.Stop();

            if (zeroX){
                sub_call_three ++;
                Eigen::Matrix<double, 3, 3> grad(3, 3);
                grad.row(0) = gradList_eigen.row(funcIndex1) - gradList_eigen.row(funcIndex2);
                grad.row(1) = gradList_eigen.row(funcIndex2) - gradList_eigen.row(funcIndex3);
                grad.row(2) = gradList_eigen.row(funcIndex3) - gradList_eigen.row(funcIndex4);
                Eigen::Matrix<double, 3, 16> grad_matrix(3, 16);
                grad_matrix.row(0) = diffList.row(funcIndex1) - diffList.row(funcIndex2);
                grad_matrix.row(1) = diffList.row(funcIndex2) - diffList.row(funcIndex3);
                grad_matrix.row(2) = diffList.row(funcIndex3) - diffList.row(funcIndex4);
                if (three_func_check (grad, grad_matrix.transpose(), sqD, threshold)){
                    //timer.Stop();
                    return true;
                }
            }
        }
        //timer.Stop();
    }
    return false;
}


