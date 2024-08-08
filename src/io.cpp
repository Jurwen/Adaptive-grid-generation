//
//  init.cpp
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 8/1/24.
//

#include "io.h"

using namespace mtet;

using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;

bool save_mesh_json(const std::string& filename,
                    const mtet::MTetMesh mesh)
{
    std::vector<std::array<double, 3>> vertices((int)mesh.get_num_vertices());
    std::vector<std::array<size_t, 4>> tets((int)mesh.get_num_tets());
    IndexMap vertex_tag_map;
    vertex_tag_map.reserve(mesh.get_num_vertices());
    int counter = 0;
    mesh.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data){
        size_t vertex_tag = vertex_tag_map.size() + 1;
        vertex_tag_map[value_of(vid)] = vertex_tag;
        vertices[counter] = {data[0], data[1], data[2]};
        counter ++;
    });
    counter = 0;
    mesh.seq_foreach_tet([&](TetId, std::span<const VertexId, 4> data) {
        tets[counter] = {vertex_tag_map[value_of(data[0])] - 1, vertex_tag_map[value_of(data[1])] - 1, vertex_tag_map[value_of(data[2])] - 1, vertex_tag_map[value_of(data[3])] - 1};
        counter ++;
    });
    if (std::filesystem::exists(filename.c_str())){
        std::filesystem::remove(filename.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    json jOut;
    jOut.push_back(json(vertices));
    jOut.push_back(json(tets));
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
    return true;
}

bool save_function_json(const std::string& filename,
                        const mtet::MTetMesh mesh,
                        ankerl::unordered_dense::map<uint64_t, llvm_vecsmall::SmallVector<std::array<double, 4>, 20>> vertex_func_grad_map,
                        const size_t funcNum)
{
    std::vector<std::vector<double>> values(funcNum);
    for (size_t funcIter = 0; funcIter <  funcNum; funcIter++){
        values[funcIter].reserve(((int)mesh.get_num_vertices()));
    }
    mesh.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data){
        llvm_vecsmall::SmallVector<std::array<double, 4>, 20> func_gradList(funcNum);
        func_gradList = vertex_func_grad_map[value_of(vid)];
        for (size_t funcIter = 0; funcIter < funcNum; funcIter++){
            values[funcIter].push_back(func_gradList[funcIter][0]);
        }
    });
    if (std::filesystem::exists(filename.c_str())){
        std::filesystem::remove(filename.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    json jOut;
    for (size_t funcIter = 0; funcIter <  funcNum; funcIter++){
        json jFunc;
        jFunc["type"] = "customized";
        jFunc["value"] = values[funcIter];
        jOut.push_back(jFunc);
    }
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
    return true;
}

bool save_timings(const std::string& filename,
                  const std::array<std::string, timer_amount>& time_label,
                  const std::array<double, timer_amount>& timings)
{
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    //fout.open(filename.c_str(),std::ios::app);
    json jOut;
    for (size_t i = 0; i < timings.size(); ++i) {
        jOut[time_label[i]] = timings[i];
    }
    //
    fout << jOut << std::endl;
    fout.close();
    return true;
}

bool save_metrics(const std::string& filename,
                  const std::array<std::string, 6>& tet_metric_labels,
                  const tet_metric metric_list)
{
    // assert stats_labels.size() == stats.size()
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    //fout.open(filename.c_str(),std::ios::app);
    json jOut;
    jOut[tet_metric_labels[0]] = metric_list.total_tet;
    jOut[tet_metric_labels[1]] = metric_list.active_tet;
    jOut[tet_metric_labels[2]] = metric_list.min_radius_ratio;
    jOut[tet_metric_labels[3]] = metric_list.active_radius_ratio;
    jOut[tet_metric_labels[4]] = metric_list.two_func_check;
    jOut[tet_metric_labels[5]] = metric_list.three_func_check;
//    for (size_t i = 0; i < tet_metric.size(); ++i) {
//        jOut[tet_metric_labels[i]] = tet_metric[i];
//    }
    fout << jOut << std::endl;
    fout.close();
    return true;
}
