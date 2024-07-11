//
//  main.cpp
//  tet_subdivision
//
//  Created by Yiwen Ju on 12/2/23.
//
#pragma once

//#define Non_Robust_Test
//#define No_Multi_Check
//#define Only_ZeroX
//#define Only_Geometry
//#define No_Multi_Check_3
//#define Only_ZeroX_3
//#define Only_Geometry_3
//#define CSG_Base
//#define MI_Base
#include <iostream>
#include <string>
#include <valarray>
#include <array>
#include <SmallVector.h>
#include <timer.h>
#include <convex_hull_membership/contains.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>


/// Two checks for each grid element: zero-crossing check and distance check. This function performs these two checks under the setting of implicit arrangement(IA), constructive solid geometry(CSG) and their curve networks.
///
/// @param[in] pts          4 arrays of three-tuples represent the coordinate of tet vertices. Each subarray represents a coordinate.
/// @param[in] vals         functions' values at each vertex. Each subarray represents one function's values at four vertices.
/// @param[in] grads            functions' gradients at each vertices. Each subarray represents one functions' gradients at four vertices and in each axis-aligned partial derivative.
/// @param[in] threshold            The user-input epsilon value for the distance check.
/// @param[out] active          A `bool` represents whether the tet is containing part of the geometry, i.e., this tet passes the zero-crossing test.
///
/// @return         A `bool` represents whether the tet is "refinable", i.e., passing the zero-crossing test and contains error greater than `threshold`.
///
///
bool subTet(const std::array<std::array<double, 3>,4> &pts,
            const llvm_vecsmall::SmallVector<std::array<double,4>, 20> &vals,
            const llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> &grads,
            const double threshold,
            bool& active);
/// Two checks for each grid element: zero-crossing check and distance check. This function performs these two checks under the setting of material interface(MI) and its curve networks.
///
/// @param[in] pts          4 arrays of three-tuples represent the coordinate of tet vertices. Each subarray represents a coordinate.
/// @param[in] vals         functions' values at each vertex. Each subarray represents one function's values at four vertices.
/// @param[in] grads            functions' gradients at each vertices. Each subarray represents one functions' gradients at four vertices and in each axis-aligned partial derivative.
/// @param[in] threshold            The user-input epsilon value for the distance check.
/// @param[out] active          A `bool` represents whether the tet is containing part of the geometry, i.e., this tet passes the zero-crossing test.
///
/// @return         A `bool` represents whether the tet is "refinable", i.e., passing the zero-crossing test and contains error greater than `threshold`.
///
///
bool subMI(const std::array<std::array<double, 3>,4> &pts,
           const llvm_vecsmall::SmallVector<std::array<double,4>, 20> &vals,
           const llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> &grads,
           const double threshold,
           bool& active);

/// Enums for the current settings of implicit complexes
enum geo_obj {
    IA,
    CSG,
    MI
};

/// Enums for CSG boolean operations
enum csg_operations{
    Intersection,
    Union,
    Negation
};


///Defines the enumeration of CSG operations. In the data structure, Intersection is 0, Union is 1, and Negation is 2.
struct csg_unit{
    int operation;
    std::array<int, 2> elements;
};

///load the csg file
///
///@param[in] filename          The name of the input CSG tree
///@param[out] tree         The loaded tree structure
///
bool load_csgTree(const std::string filename, llvm_vecsmall::SmallVector<csg_unit, 20>& tree);

/// an iterative algortihm that traverses through the csg tree
///
///@param[in] csgTree           The CSG structure: a list of csg units
///@param[in] curNode           The current index in the csg structure
///@param[in] funcInt           The intervals of all the functions
///
///@param[out] std::pair
std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> iterTree(const llvm_vecsmall::SmallVector<csg_unit, 20>csgTree,const int curNode,const llvm_vecsmall::SmallVector<std::array<double , 2>, 20> funcInt);

/// Below are the variable/constant/function definitions that will only be used in `subdivide_multi.cpp`
///
///
/// Stores the index of permutations of n less than `funcNum` for the use in `subTet` and `subMI` functions.
///
extern llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<std::array<int, 4>, 100>, 3>, 20> multiple_indices;

/// Definations of varialbes before the computations to avoid repetitive defnitions.
extern int sub_call_two;
extern int sub_call_three;
extern int GLOBAL_METHOD;
extern bool curve_network;
extern llvm_vecsmall::SmallVector<csg_unit, 20> GLOBAL_CSGTREE;
extern llvm_vecsmall::SmallVector<std::array<double, 20>, 20> valList;
extern llvm_vecsmall::SmallVector<std::array<double, 16>, 20> diffList;
extern llvm_vecsmall::SmallVector<bool, 20> activeTF;
extern llvm_vecsmall::SmallVector<std::array<double , 2>, 20> funcInt;
extern llvm_vecsmall::SmallVector<std::array<double, 3>, 20> gradList;
extern size_t funcNum;

/// Constant coefficient to obtain linear interpolated values at each bezier control points
const std::array<std::array<double, 4>, 16> coeff = {{{2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1}, {0, 2, 1, 0},{0, 2, 0, 1}, {1, 2, 0, 0}, {0, 0, 2, 1}, {1, 0, 2, 0},{0, 1, 2, 0}, {1, 0, 0, 2}, {0, 1, 0, 2}, {0, 0, 1, 2},{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}}};

/// Below are the local functions servicing `subTet` and `subMI`

/// returns a `bool` value that `true` represents positive and `false` represents negative of the input value `x`.
bool get_sign(double x);

/// returns the dot product of input arrays `a` and `b`.
double dot(const std::array<double, 3> &a, const std::array<double, 3> &b);

/// returns the dot product of input arrays `a` and `b`.
double dot(const std::array<double, 2> &a, const std::array<double, 2> &b);

/// returns the dot product of input SmallVector `a` and `b`.
double dot(const llvm_vecsmall::SmallVector<double, 20> &a, llvm_vecsmall::SmallVector<double, 20> &b);

/// returns the norm of the input array `a` where it represents a vector in 3D.
double norm(const std::array<double, 3> &a);

/// returns a vector in space by subtracting a 3D coordinate `p2` from `p1`.
std::array<double, 3> getVec(const std::array<double, 3> &p1, const std::array<double, 3> &p2);


/// returns a vector in space by subtracting a 3D coordinate `p2` from `p1`. Output is in Eigen vector.
Eigen::Vector3d getEigenVec(const std::array<double, 3> &p1, const std::array<double, 3> &p2);

/// returns a perpendicular vector for the 2D vector `a`.
std::array<double, 2> perp(const std::array<double, 2> &a);

/// returns the cross product of 3D vectors `a` and `b`.
llvm_vecsmall::SmallVector<double, 20> cross(const std::array<double, 3> &a, const std::array<double, 3> &b);

/// returns the determinant of the matrix {vec1, vec2, vec3}. Under testing, such computation is faster than turning arrays into `Eigen::Vec3d`
double det(const std::array<double, 3>& vec1,
           const std::array<double, 3>& vec2,
           const std::array<double, 3>& vec3);
/// returns the determinant of the matrix {vec1, vec2}
double det(const std::array<double, 2>& vec1,
           const std::array<double, 2>& vec2);

/// transforms the input of errors at 20 bezier control points for two functions into the correct format that `convex_hull_membership` library can use.
std::array<double, 40> transpose2d(const std::array<std::array<double, 20>, 2>& matrix);

/// transforms the input of errors at 20 bezier control points for three functions into the correct format that `convex_hull_membership` library can use.
std::array<double, 60> transpose3d(const std::array<std::array<double, 20>, 3>& matrix);

/// previous matrix multiple function. Now it's replaced by `eigen` package.
llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> matrixMultiply(const llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> &matA, const llvm_vecsmall::SmallVector<llvm_vecsmall::SmallVector<double, 20>, 20> &matB);



//const Eigen::Matrix<double, 16, 4> linear_coeff {{2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1}, {0, 2, 1, 0},{0, 2, 0, 1}, {1, 2, 0, 0}, {0, 0, 2, 1}, {1, 0, 2, 0},{0, 1, 2, 0}, {1, 0, 0, 2}, {0, 1, 0, 2}, {0, 0, 1, 2},{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}};

//typedef Eigen::Triplet<double> T;
//std::vector<T> tripletList;
//tripletList.reserve(1);
//tripletList[0] = {0, 0, 0.0};
//
//Eigen::SparseMatrix<double> coeff_matrix(20, 20);
////coeff_matrix.insert(0,0) = 1;
//
////std::vector< Eigen::Triplet<double>> tripletList;
//
//void matrix_init(Eigen::SparseMatrix<double>& coeff_matrix){
//    coeff_matrix.insert(0,0) = 1;
//    coeff_matrix.insert(1,1) = 1;
//    coeff_matrix.insert(2,2) = 1;
//    coeff_matrix.insert(3,3) = 1;
//    
//    coeff_matrix.insert(4,0) = 1;
//    coeff_matrix.insert(4,4) = -1.0/3.0;
//    coeff_matrix.insert(4,5) = 1.0/3.0;
//    coeff_matrix.insert(5,0) = 1;
//    coeff_matrix.insert(5,4) = -1.0/3.0;
//    coeff_matrix.insert(5,6) = 1.0/3.0;
//    coeff_matrix.insert(6,0) = 1;
//    coeff_matrix.insert(6,4) = -1.0/3.0;
//    coeff_matrix.insert(6,7) = 1.0/3.0;
//    
//    coeff_matrix.insert(7,1) = 1;
//    coeff_matrix.insert(7,9) = -1.0/3.0;
//    coeff_matrix.insert(7,10) = 1.0/3.0;
//    coeff_matrix.insert(8,1) = 1;
//    coeff_matrix.insert(8,9) = -1.0/3.0;
//    coeff_matrix.insert(8,11) = 1.0/3.0;
//    coeff_matrix.insert(9,1) = 1;
//    coeff_matrix.insert(9,9) = -1.0/3.0;
//    coeff_matrix.insert(9,8) = 1.0/3.0;
//    
//    coeff_matrix.insert(10,2) = 1;
//    coeff_matrix.insert(10,14) = -1.0/3.0;
//    coeff_matrix.insert(10,15) = 1.0/3.0;
//    coeff_matrix.insert(11,2) = 1;
//    coeff_matrix.insert(11,14) = -1.0/3.0;
//    coeff_matrix.insert(11,12) = 1.0/3.0;
//    coeff_matrix.insert(12,2) = 1;
//    coeff_matrix.insert(12,14) = -1.0/3.0;
//    coeff_matrix.insert(12,13) = 1.0/3.0;
//    
//    coeff_matrix.insert(13,3) = 1;
//    coeff_matrix.insert(13,19) = -1.0/3.0;
//    coeff_matrix.insert(13,16) = 1.0/3.0;
//    coeff_matrix.insert(14,3) = 1;
//    coeff_matrix.insert(14,19) = -1.0/3.0;
//    coeff_matrix.insert(14,17) = 1.0/3.0;
//    coeff_matrix.insert(15,3) = 1;
//    coeff_matrix.insert(15,19) = -1.0/3.0;
//    coeff_matrix.insert(15,18) = 1.0/3.0;
//    
//    coeff_matrix.insert(16,1) = 1.0/3.0;
//    coeff_matrix.insert(16,2) = 1.0/3.0;
//    coeff_matrix.insert(16,3) = 1.0/3.0;
//    coeff_matrix.insert(16,9) = -1.0/6.0;
//    coeff_matrix.insert(16,10) = 1.0/12.0;
//    coeff_matrix.insert(16,11) = 1.0/12.0;
//    coeff_matrix.insert(16,13) = 1.0/12.0;
//    coeff_matrix.insert(16,14) = -1.0/6.0;
//    coeff_matrix.insert(16,15) = 1.0/12.0;
//    coeff_matrix.insert(16,17) = 1.0/12.0;
//    coeff_matrix.insert(16,18) = 1.0/12.0;
//    coeff_matrix.insert(16,19) = -1.0/6.0;
//    
//    coeff_matrix.insert(17,0) = 1.0/3.0;
//    coeff_matrix.insert(17,2) = 1.0/3.0;
//    coeff_matrix.insert(17,3) = 1.0/3.0;
//    coeff_matrix.insert(17,4) = -1.0/6.0;
//    coeff_matrix.insert(17,6) = 1.0/12.0;
//    coeff_matrix.insert(17,7) = 1.0/12.0;
//    coeff_matrix.insert(17,12) = 1.0/12.0;
//    coeff_matrix.insert(17,14) = -1.0/6.0;
//    coeff_matrix.insert(17,15) = 1.0/12.0;
//    coeff_matrix.insert(17,16) = 1.0/12.0;
//    coeff_matrix.insert(17,18) = 1.0/12.0;
//    coeff_matrix.insert(17,19) = -1.0/6.0;
//    
//    coeff_matrix.insert(18,0) = 1.0/3.0;
//    coeff_matrix.insert(18,1) = 1.0/3.0;
//    coeff_matrix.insert(18,3) = 1.0/3.0;
//    coeff_matrix.insert(18,4) = -1.0/6.0;
//    coeff_matrix.insert(18,5) = 1.0/12.0;
//    coeff_matrix.insert(18,7) = 1.0/12.0;
//    coeff_matrix.insert(18,8) = 1.0/12.0;
//    coeff_matrix.insert(18,9) = -1.0/6.0;
//    coeff_matrix.insert(18,11) = 1.0/12.0;
//    coeff_matrix.insert(18,16) = 1.0/12.0;
//    coeff_matrix.insert(18,17) = 1.0/12.0;
//    coeff_matrix.insert(18,19) = -1.0/6.0;
//    
//    coeff_matrix.insert(19,0) = 1.0/3.0;
//    coeff_matrix.insert(19,1) = 1.0/3.0;
//    coeff_matrix.insert(19,2) = 1.0/3.0;
//    coeff_matrix.insert(19,4) = -1.0/6.0;
//    coeff_matrix.insert(19,5) = 1.0/12.0;
//    coeff_matrix.insert(19,6) = 1.0/12.0;
//    coeff_matrix.insert(19,8) = 1.0/12.0;
//    coeff_matrix.insert(19,9) = -1.0/6.0;
//    coeff_matrix.insert(19,10) = 1.0/12.0;
//    coeff_matrix.insert(19,12) = 1.0/12.0;
//    coeff_matrix.insert(19,13) = 1.0/12.0;
//    coeff_matrix.insert(19,14) = -1.0/6.0;
//}
//const Eigen::Matrix<double, 20, 20> coeff_matrix{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//v0
//{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//v1
//{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//v2
//{0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//v3
//{1, 0, 0, 0, -1.0/3.0, 1.0/3.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//v0s_0
//{1, 0, 0, 0, -1.0/3.0, 0, 1.0/3.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//v0s_1
//{1, 0, 0, 0, -1.0/3.0, 0, 0, 1.0/3.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//v0s_2
//{0, 1, 0, 0, 0, 0, 0, 0, 0, -1.0/3.0, 1.0/3.0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//v1s_0
//{0, 1, 0, 0, 0, 0, 0, 0, 0, -1.0/3.0, 0, 1.0/3.0, 0, 0, 0, 0, 0, 0, 0, 0},//v1s_1
//{0, 1, 0, 0, 0, 0, 0, 0, 1.0/3.0, -1.0/3.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//v1s_2
//{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0/3.0, 1.0/3.0, 0, 0, 0, 0},//v2s_0
//{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0/3.0, 0, -1.0/3.0, 0, 0, 0, 0, 0},//v2s_1
//{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0/3.0, -1.0/3.0, 0, 0, 0, 0, 0},//v2s_2
//{0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0/3.0, 0, 0, -1.0/3.0},//v3s_0
//{0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0/3.0, 0, -1.0/3.0},//v3s_1
//{0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0/3.0, -1.0/3.0},//v3s_2
//{0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 0, 0, 0, 0, 0, -1.0/6.0, 1.0/12.0, 1.0/12.0, 0, 1.0/12.0, -1.0/6.0, 1.0/12.0, 0, 1.0/12.0, 1.0/12.0, -1.0/6.0},//vMid0
//{1.0/3.0, 0, 1.0/3.0, 1.0/3.0, -1.0/6.0, 0, 1.0/12.0, 1.0/12.0, 0, 0, 0, 0, 1.0/12.0, 0, -1.0/6.0, 1.0/12.0, 1.0/12.0, 0, 1.0/12.0, -1.0/6.0},//vMid1
//{1.0/3.0, 1.0/3.0, 0, 1.0/3.0, -1.0/6.0, 1.0/12.0, 0, 1.0/12.0, 1.0/12.0, -1.0/6.0, 0, 1.0/12.0, 0, 0, 0, 0, 1.0/12.0, 1.0/12.0, 0, -1.0/6.0},//vMid2
//{1.0/3.0, 1.0/3.0, 1.0/3.0, 0, -1.0/6.0, 1.0/12.0, 1.0/12.0, 0, 1.0/12.0, -1.0/6.0, 1.0/12.0, 0, 1.0/12.0, 1.0/12.0, -1.0/6.0, 0, 0, 0, 0, 0}};//vMid3

