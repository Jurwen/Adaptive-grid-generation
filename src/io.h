//
//  init.h
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 8/1/24.
//

#pragma once

#include <ankerl/unordered_dense.h>
#include "adaptive_grid_gen.h"
#include "timer.h"

struct tet_metric {
    size_t total_tet = 0;
    int active_tet = 0;
    double min_radius_ratio = 1;
    double active_radius_ratio = 1;
    int two_func_check = 0;
    int three_func_check = 0;
};

bool save_mesh_json(const std::string& filename,
                    const mtet::MTetMesh mesh);

bool save_function_json(const std::string& filename,
                        const mtet::MTetMesh mesh,
                        ankerl::unordered_dense::map<uint64_t, llvm_vecsmall::SmallVector<std::array<double, 4>, 20>> vertex_func_grad_map,
                        const size_t funcNum);

/// saves the timing profile to a file
/// @param[in] time_labels            The labels of the timings.
/// @param[in] timings           The values of the timings, corresponding to the lables.
/// @param[in] filename            The name of the output file.
///
/// @return         Whether this saving procedure is successful.
bool save_timings(const std::string& filename,
                  const std::array<std::string, timer_amount>& time_label,
                  const std::array<double, timer_amount>& timings);

/// saves the tet metrics to a file
/// @param[in] tet_metric_labels            The labels of the metrics.
/// @param[in] tet_metric           The values of the metrics, corresponding to the lables.
/// @param[in] filename            The name of the output file.
///
/// @return         Whether this saving procedure is successful.
bool save_metrics(const std::string& filename,
                  const std::array<std::string, 6>& tet_metric_labels,
                  const tet_metric metric_list);
