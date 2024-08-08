//
//  csg.hpp
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 8/1/24.
//

#ifndef csg_h
#define csg_h
#include "adaptive_grid_gen.h"

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

#endif /* csg_h */
