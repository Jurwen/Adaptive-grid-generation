//
//  csg.cpp
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 8/1/24.
//

#include "csg.h"

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
            tree[j] = {Intersection, elements};
        }else if (type == "Union"){
            tree[j] = {Union, elements};
        }else if (type == "Negation"){
            tree[j] = {Negation, elements};
        }
    }
    return true;
}

std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> iterTree(const llvm_vecsmall::SmallVector<csg_unit, 20>csgTree,const int curNode,const llvm_vecsmall::SmallVector<std::array<double , 2>, 20> funcInt){
    csg_unit curUnit = csgTree[curNode - 1];
    std::array<double, 2> interval, childInt1, childInt2;
    llvm_vecsmall::SmallVector<int, 20> af(funcInt.size(), 1), childAF1(funcInt.size(), 1), childAF2(funcInt.size(), 1);
    if (curUnit.elements[0] > 0){
        std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> child1 = iterTree(csgTree, curUnit.elements[0], funcInt);
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
            std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> child2 = iterTree(csgTree, curUnit.elements[1], funcInt);
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
