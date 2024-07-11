//#define Check_Flip_Tets
#include <mtet/mtet.h>
#include <mtet/io.h>
#include <ankerl/unordered_dense.h>
#include <span>
#include <queue>
#include <optional>
#include <SmallVector.h>
#include <implicit_functions.h>
#include <CLI/CLI.hpp>

#include "subdivide_multi.h"
#include "tet_quality.h"
#include "timer.h"
#include "grid_mesh.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>


using namespace mtet;

bool save_mesh_json(const std::string& filename,
                    const mtet::MTetMesh mesh)
{
    std::vector<std::array<double, 3>> vertices((int)mesh.get_num_vertices());
    std::vector<std::array<size_t, 4>> tets((int)mesh.get_num_tets());
    using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;
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
//hash for mounting a boolean that represents the activeness to a tet
//since the tetid isn't const during the process, mount the boolean using vertexids of 4 corners.
uint64_t vertexHash(std::span<VertexId, 4>& x)
{
    ankerl::unordered_dense::hash<uint64_t> hash_fn;
    return hash_fn(value_of(x[0])) + hash_fn(value_of(x[1])) + hash_fn(value_of(x[2])) + hash_fn(value_of(x[3]));
}

int main(int argc, const char *argv[])
{
    struct
    {
        std::string mesh_file;
        std::string function_file;
        double threshold = 0.0001;
        double alpha = std::numeric_limits<double>::infinity();
        int max_elements = -1;
        double smallest_edge_length = 0;
        std::string method = "IA";
        std::string csg_file;
        bool bfs = false;
        bool dfs = false;
        bool curve_network = false;
        //bool analysis_mode = false;
    } args;
    CLI::App app{"Longest Edge Bisection Refinement"};
    app.add_option("mesh", args.mesh_file, "Initial mesh file")->required();
    app.add_option("function", args.function_file, "Implicit function file")->required();
    app.add_option("-t,--threshold", args.threshold, "Threshold value");
    app.add_option("-a,--alpha", args.alpha, "Alpha value");
    app.add_option("-o,--option", args.method, "Options of implicit manifold");
    app.add_option("--tree", args.csg_file, "CSG Tree file");
    app.add_option("-m,--max-elements", args.max_elements, "Maximum number of elements");
    app.add_option("-s,--shortest-edge", args.smallest_edge_length, "Shortest edge length");
    app.add_option("-b, --bfs", args.bfs, "Toggle BFS Mode");
    app.add_option("-d, --dfs", args.dfs, "Toggle DFS Mode");
    app.add_option("-c, --curve_network", args.curve_network, "Generate Curve Network only");
    CLI11_PARSE(app, argc, argv);
    // Read mesh
    mtet::MTetMesh mesh;
    if (args.mesh_file.find(".json") != std::string::npos){
        mesh = grid_mesh::load_tet_mesh(args.mesh_file);
        mtet::save_mesh("init.msh", mesh);
        mesh = mtet::load_mesh("init.msh");
    } else {
        mesh = mtet::load_mesh(args.mesh_file);
    }
    
    // Read implicit function
    std::vector<std::unique_ptr<ImplicitFunction<double>>> functions;
    load_functions(args.function_file, functions);
    funcNum = functions.size();
    valList.reserve(funcNum);
    diffList.reserve(funcNum);
    activeTF.reserve(funcNum);
    funcInt.reserve(funcNum);
    gradList.reserve(funcNum);
    // Read options
    if (args.max_elements < 0)
    {
        args.max_elements = std::numeric_limits<int>::max();
    }
    double threshold = args.threshold;
    double alpha = args.alpha;
    double smallest_edge_length = args.smallest_edge_length;
    if (args.method == "IA"){
        GLOBAL_METHOD = IA;
    }
    if (args.method == "CSG"){
        GLOBAL_METHOD = CSG;
        load_csgTree(args.csg_file, GLOBAL_CSGTREE);
    }
    if (args.method == "MI"){
        GLOBAL_METHOD = MI;
    }
    if (args.curve_network){
        curve_network = true;
    }
    
    //precomputing active multiples' indices:
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
                        if (GLOBAL_METHOD == MI){
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
        if (GLOBAL_METHOD == MI){
            multiple_indices[funcIter] = {pair, triple, quad};
        }else{
            multiple_indices[funcIter] = {pair, triple};
        }
    }
    int search_counter;
    if (args.bfs || args.dfs){
        search_counter = 0;
    }
    //matrix_init(coeff_matrix);
    // initialize vertex map: vertex index -> {{f_i, gx, gy, gz} | for all f_i in the function}
    using IndexMap = ankerl::unordered_dense::map<uint64_t, llvm_vecsmall::SmallVector<std::array<double, 4>, 20>>;
    IndexMap vertex_func_grad_map;
    vertex_func_grad_map.reserve(mesh.get_num_vertices());
    
    //initialize activeness map: four vertexids (v0, v1, v2, v3) -> hash(v0, v1, v2, v3) -> active boolean
    using activeMap = ankerl::unordered_dense::map<uint64_t, bool>;
    activeMap vertex_active_map;
    vertex_active_map.reserve(mesh.get_num_tets());
    
    mesh.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data)
                            {
        llvm_vecsmall::SmallVector<std::array<double, 4>, 20> func_gradList(funcNum);
        for(size_t funcIter = 0; funcIter < funcNum; funcIter++){
            auto &func = functions[funcIter];
            std::array<double, 4> func_grad;
            func_grad[0] = func->evaluate_gradient(data[0], data[1], data[2], func_grad[1], func_grad[2], func_grad[3]);
            func_gradList[funcIter] = func_grad;
        }
        vertex_func_grad_map[value_of(vid)] = func_gradList;});
    auto comp = [](std::pair<mtet::Scalar, mtet::EdgeId> e0,
                   std::pair<mtet::Scalar, mtet::EdgeId> e1)
    { return e0.first < e1.first; };
    std::vector<std::pair<mtet::Scalar, mtet::EdgeId>> Q;
    
    std::array<std::array<double, 3>, 4> pts;
    llvm_vecsmall::SmallVector<std::array<double, 4>, 20> vals(funcNum);
    llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> grads(funcNum);
    double activeTet = 0;
    auto push_longest_edge = [&](mtet::TetId tid)
    {
        std::span<VertexId, 4> vs = mesh.get_tet(tid);
        {
            //Timer eval_timer(evaluation, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            for (int i = 0; i < 4; ++i)
            {
                auto vid = vs[i];
                auto coords = mesh.get_vertex(vid);
                pts[i][0] = coords[0];
                pts[i][1] = coords[1];
                pts[i][2] = coords[2];
                llvm_vecsmall::SmallVector<std::array<double, 4>, 20> func_gradList(funcNum);
                std::array<double, 4> func_grad;
                if (!vertex_func_grad_map.contains(value_of(vid))) {
                    for(size_t funcIter = 0; funcIter < funcNum; funcIter++){
                        auto &func = functions[funcIter];
                        std::array<double, 4> func_grad;
                        func_grad[0] = func->evaluate_gradient(coords[0], coords[1], coords[2], func_grad[1], func_grad[2],
                                                               func_grad[3]);
                        func_gradList[funcIter] = func_grad;
                    }
                    vertex_func_grad_map[value_of(vid)] = func_gradList;
                }
                else {
                    func_gradList = vertex_func_grad_map[value_of(vid)];
                }
                for(size_t funcIter = 0; funcIter < funcNum; funcIter++){
                    vals[funcIter][i] = func_gradList[funcIter][0];
                    grads[funcIter][i][0] = func_gradList[funcIter][1];
                    grads[funcIter][i][1] = func_gradList[funcIter][2];
                    grads[funcIter][i][2] = func_gradList[funcIter][3];
                }
            }
            //eval_timer.Stop();
        }
        bool isActive = 0;
        bool subResult;
        {
            //Timer sub_timer(subdivision, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            if (GLOBAL_METHOD != MI){
                subResult = subTet(pts, vals, grads, threshold, isActive);
            }else{
                subResult = subMI(pts, vals, grads, threshold, isActive);
            }
            //sub_timer.Stop();
        }
        vertex_active_map[vertexHash(vs)] = isActive;
        //Timer eval_timer(evaluation, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        if (subResult)
        {
            mtet::EdgeId longest_edge;
            mtet::Scalar longest_edge_length = 0;
            mesh.foreach_edge_in_tet(tid, [&](mtet::EdgeId eid, mtet::VertexId v0, mtet::VertexId v1)
                                     {
                auto p0 = mesh.get_vertex(v0);
                auto p1 = mesh.get_vertex(v1);
                mtet::Scalar l = (p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                (p0[2] - p1[2]) * (p0[2] - p1[2]);
                if (l > longest_edge_length) {
                    longest_edge_length = l;
                    longest_edge = eid;
                } });
            if (args.bfs){
                Q.emplace_back(search_counter, longest_edge);
                search_counter--;
            }else if(args.dfs){
                Q.emplace_back(search_counter, longest_edge);
                search_counter++;
            }
            else{
                Q.emplace_back(longest_edge_length, longest_edge);
            }
            //eval_timer.Stop();
            return true;
        }
        //eval_timer.Stop();
        return false;
    };
    
    
    {
        Timer timer(total_time, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        
        // Initialize priority queue.
        mesh.seq_foreach_tet([&](mtet::TetId tid, [[maybe_unused]] std::span<const mtet::VertexId, 4> vs)
                             { push_longest_edge(tid); });
        std::make_heap(Q.begin(), Q.end(), comp);
        
        // Keep splitting the longest edge
        while (!Q.empty())
        {
            std::pop_heap(Q.begin(), Q.end(), comp);
            auto [edge_length, eid] = Q.back();
            if (!mesh.has_edge(eid)){
                Q.pop_back();
                continue;
            }
            //implement alpha value:
            mtet::Scalar comp_edge_length = alpha * edge_length;
            bool addedActive = false;
            mesh.foreach_tet_around_edge(eid,[&](mtet::TetId tid){
                std::span<VertexId, 4> vs = mesh.get_tet(tid);
                if(vertex_active_map.contains(vertexHash(vs))){
                    if (vertex_active_map[vertexHash(vs)]){
                        mtet::EdgeId longest_edge;
                        mtet::Scalar longest_edge_length = 0;
                        mesh.foreach_edge_in_tet(tid, [&](mtet::EdgeId eid_active, mtet::VertexId v0, mtet::VertexId v1)
                                                 {
                            auto p0 = mesh.get_vertex(v0);
                            auto p1 = mesh.get_vertex(v1);
                            mtet::Scalar l = (p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                            (p0[2] - p1[2]) * (p0[2] - p1[2]);
                            if (l > longest_edge_length) {
                                longest_edge_length = l;
                                longest_edge = eid_active;
                            }
                        });
                        if (longest_edge_length > comp_edge_length) {
                            Q.emplace_back(longest_edge_length, longest_edge);
                            addedActive = true;
                        }
                    }
                }
            });
            if(addedActive){
                std::push_heap(Q.begin(), Q.end(), comp);
                continue;
            }
            Q.pop_back();
            std::array<VertexId, 2> vs_old = mesh.get_edge_vertices(eid);
            //Timer split_timer(splitting, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            auto [vid, eid0, eid1] = mesh.split_edge(eid);
            //split_timer.Stop();
            //std::cout << "Number of elements: " << mesh.get_num_tets() << std::endl;
            if (mesh.get_num_tets() > args.max_elements) {
                break;
            }
            mesh.foreach_tet_around_edge(eid0, [&](mtet::TetId tid)
                                         {
                if (push_longest_edge(tid)) {
                    std::push_heap(Q.begin(), Q.end(), comp);
                } });
            mesh.foreach_tet_around_edge(eid1, [&](mtet::TetId tid)
                                         {
                if (push_longest_edge(tid)) {
                    std::push_heap(Q.begin(), Q.end(), comp);
                } });
#ifdef Check_Flip_Tets
            std::array<VertexId, 2> vs_new = mesh.get_edge_vertices(eid0);
            mesh.foreach_tet_around_edge(eid0, [&](mtet::TetId tid)
                                         {
                std::span<VertexId, 4> vs = mesh.get_tet(tid);
                std::vector<VertexId> parent(4);
                parent[0] = vs_old[0]; parent[1] = vs_old[1];
                int parentIndex = 2;
                for (auto vIter : vs){
                    if (vIter !=  vs_new[0] && vIter != vs_new[1]){
                        parent[parentIndex] = vIter;
                        parentIndex ++;
                    }
                }
                std::span<VertexId, 4> spanVec(parent);
                if (!vertex_active_map[vertexHash(spanVec)] && vertex_active_map[vertexHash(vs)]){
                    {
                        using json = nlohmann::json;
                        std::string filePath = "flip_tets.json";
                        //                        if (std::filesystem::exists(filePath)) {
                        //                            std::filesystem::remove(filePath);
                        //                        }
                        std::ofstream fout(filePath,std::ios::app);
                        json jOut;
                        std::array<std::array<double, 3>, 4> pts;
                        llvm_vecsmall::SmallVector<std::array<double, 4>, 20> vals(funcNum);
                        llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> grads(funcNum);
                        for (size_t i = 0; i < 4; ++i) {
                            auto coords = mesh.get_vertex(vs[i]);
                            pts[i][0] = coords[0];
                            pts[i][1] = coords[1];
                            pts[i][2] = coords[2];
                            auto func_gradList = vertex_func_grad_map[value_of(vs[i])];
                            for(size_t funcIter = 0; funcIter < funcNum; funcIter++){
                                vals[funcIter][i] = func_gradList[funcIter][0];
                                grads[funcIter][i][0] = func_gradList[funcIter][1];
                                grads[funcIter][i][1] = func_gradList[funcIter][2];
                                grads[funcIter][i][2] = func_gradList[funcIter][3];
                            }
                        }
                        jOut["vertices: "] = pts;
                        jOut["value: "] = vals;
                        jOut["gradient: "] = grads;
                        //
                        fout << jOut << std::endl;
                        fout.close();
                    }
                    {
                        using json = nlohmann::json;
                        std::string filePath = "flip_tets_parent.json";
                        //                        if (std::filesystem::exists(filePath)) {
                        //                            std::filesystem::remove(filePath);
                        //                        }
                        std::ofstream fout(filePath,std::ios::app);
                        json jOut;
                        std::array<std::array<double, 3>, 4> pts;
                        llvm_vecsmall::SmallVector<std::array<double, 4>, 20> vals(funcNum);
                        llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> grads(funcNum);
                        for (size_t i = 0; i < 4; ++i) {
                            auto coords = mesh.get_vertex(spanVec[i]);
                            pts[i][0] = coords[0];
                            pts[i][1] = coords[1];
                            pts[i][2] = coords[2];
                            auto func_gradList = vertex_func_grad_map[value_of(spanVec[i])];
                            for(size_t funcIter = 0; funcIter < funcNum; funcIter++){
                                vals[funcIter][i] = func_gradList[funcIter][0];
                                grads[funcIter][i][0] = func_gradList[funcIter][1];
                                grads[funcIter][i][1] = func_gradList[funcIter][2];
                                grads[funcIter][i][2] = func_gradList[funcIter][3];
                            }
                        }
                        jOut["vertices: "] = pts;
                        jOut["value: "] = vals;
                        jOut["gradient: "] = grads;
                        //
                        fout << jOut << std::endl;
                        fout.close();
                    }
                }
            });
            vs_new = mesh.get_edge_vertices(eid1);
            mesh.foreach_tet_around_edge(eid1, [&](mtet::TetId tid)
                                         {
                std::span<VertexId, 4> vs = mesh.get_tet(tid);
                std::vector<VertexId> parent(4);
                parent[0] = vs_old[0]; parent[1] = vs_old[1];
                int parentIndex = 2;
                for (auto vIter : vs){
                    if (vIter !=  vs_new[0] && vIter != vs_new[1]){
                        parent[parentIndex] = vIter;
                        parentIndex ++;
                    }
                }
                std::span<VertexId, 4> spanVec(parent);
                if (!vertex_active_map[vertexHash(spanVec)] && vertex_active_map[vertexHash(vs)]){
                    {
                        using json = nlohmann::json;
                        std::string filePath = "flip_tets.json";
                        //                        if (std::filesystem::exists(filePath)) {
                        //                            std::filesystem::remove(filePath);
                        //                        }
                        std::ofstream fout(filePath,std::ios::app);
                        json jOut;
                        std::array<std::array<double, 3>, 4> pts;
                        llvm_vecsmall::SmallVector<std::array<double, 4>, 20> vals(funcNum);
                        llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> grads(funcNum);
                        for (size_t i = 0; i < 4; ++i) {
                            auto coords = mesh.get_vertex(vs[i]);
                            pts[i][0] = coords[0];
                            pts[i][1] = coords[1];
                            pts[i][2] = coords[2];
                            auto func_gradList = vertex_func_grad_map[value_of(vs[i])];
                            for(size_t funcIter = 0; funcIter < funcNum; funcIter++){
                                vals[funcIter][i] = func_gradList[funcIter][0];
                                grads[funcIter][i][0] = func_gradList[funcIter][1];
                                grads[funcIter][i][1] = func_gradList[funcIter][2];
                                grads[funcIter][i][2] = func_gradList[funcIter][3];
                            }
                        }
                        jOut["vertices: "] = pts;
                        jOut["value: "] = vals;
                        jOut["gradient: "] = grads;
                        //
                        fout << jOut << std::endl;
                        fout.close();
                    }
                    {
                        using json = nlohmann::json;
                        std::string filePath = "flip_tets_parent.json";
                        //                        if (std::filesystem::exists(filePath)) {
                        //                            std::filesystem::remove(filePath);
                        //                        }
                        std::ofstream fout(filePath,std::ios::app);
                        json jOut;
                        std::array<std::array<double, 3>, 4> pts;
                        llvm_vecsmall::SmallVector<std::array<double, 4>, 20> vals(funcNum);
                        llvm_vecsmall::SmallVector<std::array<std::array<double, 3>,4>, 20> grads(funcNum);
                        for (size_t i = 0; i < 4; ++i) {
                            auto coords = mesh.get_vertex(spanVec[i]);
                            pts[i][0] = coords[0];
                            pts[i][1] = coords[1];
                            pts[i][2] = coords[2];
                            auto func_gradList = vertex_func_grad_map[value_of(spanVec[i])];
                            for(size_t funcIter = 0; funcIter < funcNum; funcIter++){
                                vals[funcIter][i] = func_gradList[funcIter][0];
                                grads[funcIter][i][0] = func_gradList[funcIter][1];
                                grads[funcIter][i][1] = func_gradList[funcIter][2];
                                grads[funcIter][i][2] = func_gradList[funcIter][3];
                            }
                        }
                        jOut["vertices: "] = pts;
                        jOut["value: "] = vals;
                        jOut["gradient: "] = grads;
                        //
                        fout << jOut << std::endl;
                        fout.close();
                    }
                }
            });
#endif
        }
        timer.Stop();
    }
    //profiled time(see details in time.h) and profiled number of calls to zero
    for (int i = 0; i < profileTimer.size(); i++){
        timeProfileName time_type = static_cast<timeProfileName>(i);
        std::cout << time_label[i] << ": " << profileTimer[i] << std::endl;
    }
    //    std::cout << profileTimer[0] << " "<< profileTimer[1] << " "<< profileTimer[2] << " "<< profileTimer[3] << " "<< profileTimer[4] << " "<< profileTimer[5] << " "<< profileTimer[6] << " "<< profileTimer[7] << " "<< profileTimer[8] << " "<< profileTimer[9] << " "<< sub_call_two << " "<< sub_call_three << std::endl;
    //std::cout << "sub two func calls: " << sub_call_two << std::endl;
    //std::cout << "sub three func calls: " << sub_call_three << std::endl;
    double min_rratio_all = 1;
    double min_rratio_active = 1;
    std::vector<mtet::TetId> activeTetId;
    mesh.seq_foreach_tet([&](mtet::TetId tid, std::span<const VertexId, 4> data) {
        std::span<VertexId, 4> vs = mesh.get_tet(tid);
        std::array<std::valarray<double>,4> vallPoints;
        for (int i = 0; i < 4; i++){
            vallPoints[i] = {0.0,0.0,0.0};
        }
        for (int i = 0; i < 4; i++){
            VertexId vid = vs[i];
            std::span<Scalar, 3> coords = mesh.get_vertex(vid);
            vallPoints[i][0] = coords[0];
            vallPoints[i][1] = coords[1];
            vallPoints[i][2] = coords[2];
        }
        double ratio = tet_radius_ratio(vallPoints);
        if (ratio < min_rratio_all){
            min_rratio_all = ratio;
            
        }
        if(vertex_active_map.contains(vertexHash(vs))){
            if (vertex_active_map[vertexHash(vs)]){
                activeTet++;
                activeTetId.push_back(tid);
                if (ratio < min_rratio_active){
                    min_rratio_active = ratio;
                }
            }
        }
    });
    // save timing records
    save_timings("timings.json",time_label, profileTimer);
    // save statistics
    save_metrics("stats.json", tet_metric_labels, {(double)mesh.get_num_tets(), activeTet, min_rratio_all, min_rratio_active, (double)sub_call_two, (double) sub_call_three});
    // save the mesh output for isosurfacing tool
    //save_mesh_json("mesh.json", mesh);
    // save the mesh output for isosurfacing tool
    //save_function_json("function_value.json", mesh, vertex_func_grad_map, funcNum);
    //write mesh and active tets
    mtet::save_mesh("tet_mesh.msh", mesh);
    mtet::save_mesh("active_tets.msh", mesh, std::span<mtet::TetId>(activeTetId));
    return 0;
}
