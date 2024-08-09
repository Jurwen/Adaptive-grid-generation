//
//  mesh_refine.h
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 8/4/24.
//

#pragma once

#include <SmallVector.h>
#include <implicit_functions.h>

#include "adaptive_grid_gen.h"
#include "io.h"
#include "refine_crit.h"
#include "tet_quality.h"

using namespace mtet;

bool gridRefine(
                const std::string function_file,
                const int mode,
                const double threshold,
                const double alpha,
                const int max_elements,
                const std::function<std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>>(llvm_vecsmall::SmallVector<std::array<double, 2>, 20>)> csg_func,
                mtet::MTetMesh &grid,
                tet_metric &metric_list,
                std::array<double, timer_amount> profileTimer
                );
