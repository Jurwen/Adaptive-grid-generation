#include <mtet/mtet.h>
#include <mtet/io.h>
#include <ankerl/unordered_dense.h>
#include <span>
#include <queue>
#include <optional>
#include <SmallVector.h>

#include <implicit_functions.h>
#include <subdivide_multi.h>
#include <CLI/CLI.hpp>
#include <tet_quality.h>
#include <timer.h>


using namespace mtet;

const bool GLOBAL_ANALYSIS_MODE = false;

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
        string mesh_file;
        string function_file;
        double threshold = 0.0001;
        int max_elements = -1;
        bool bfs = false;
        //bool analysis_mode = false;
    } args;
    CLI::App app{"Longest Edge Bisection Refinement"};
    app.add_option("mesh", args.mesh_file, "Initial mesh file")->required();
    app.add_option("function", args.function_file, "Implicit function file")->required();
    app.add_option("-t,--threshold", args.threshold, "Threshold value");
    app.add_option("-m,--max-elements", args.max_elements, "Maximum number of elements");
    app.add_option("-b, --bfs", args.bfs, "toggle BFS Mode");
    CLI11_PARSE(app, argc, argv);
    //extern const bool analysis = args.analysis_mode;
    if (args.max_elements < 0)
    {
        args.max_elements = numeric_limits<int>::max();
    }
    double threshold = args.threshold;
    
    // Read mesh
    mtet::MTetMesh mesh = mtet::load_mesh(args.mesh_file);
    
    // Read implicit function
    vector<unique_ptr<ImplicitFunction<double>>> functions;
    load_functions(args.function_file, functions);
    size_t funcNum = functions.size();
    long long int largeNumber;
    if (args.bfs){
        largeNumber = 99999;
        //std::numeric_limits<long long int>::max();
    }
    //setup gurobi env
//    GRBEnv env = GRBEnv(true);
//    env.set(GRB_IntParam_OutputFlag, 0);
//    env.set("LogFile", "");
//    env.start();
    
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
            array<double, 4> func_grad;
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
//    size_t sdim = 3;
//    for (int i = 0; i < 4; ++i)
//    {
//        pts[i].resize(sdim);
//        for (size_t j = 0; j < funcNum; j++){
//            grads[j][i].resize(sdim);
//        }
//    }
    
    auto push_longest_edge = [&](mtet::TetId tid)
    {
        std::span<VertexId, 4> vs = mesh.get_tet(tid);
        bool isActive = 0;
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
                    array<double, 4> func_grad;
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
        double subResult;
        {
            Timer sub_timer(subdivision, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            subResult = subTet(pts, vals, grads, threshold, isActive);
            sub_timer.Stop();
        }
        vertex_active_map[vertexHash(vs)] = isActive;

        if (subResult != -1)
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
                Q.emplace_back(largeNumber, longest_edge);
                largeNumber--;
            }else{
                Q.emplace_back(longest_edge_length, longest_edge);
            }
            return true;
        }
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
            Q.pop_back();
            if (!mesh.has_edge(eid))
                continue;
            Timer split_timer(splitting, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            auto [vid, eid0, eid1] = mesh.split_edge(eid);
            split_timer.Stop();
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
        }
        timer.Stop();
    }
    
    //profiled time(see details in time.h) and profiled number of calls to zero
    std::cout << profileTimer[0] << " " << profileTimer[1] << " " << profileTimer[2] << " " << profileTimer[3] << " " << profileTimer[4] << " " << profileTimer[5] << " " << profileTimer[6] << " " << profileTimer[7] << " " << sub_call_two << " " << sub_call_three << std::endl;
    double min_rratio_all = 1, min_rratio_active = 1;
    std::array<valarray<double>,4> valPoints;
    for (int i = 0; i < 4; i ++){
        valPoints[i].resize(3);
    }
    mesh.seq_foreach_tet([&](mtet::TetId tid, std::span<const VertexId, 4> data) {
        std::span<VertexId, 4> vs = mesh.get_tet(tid);
        for (int i = 0; i < 4; i++){
            auto vid = vs[i];
            auto coords = mesh.get_vertex(vid);
            valPoints[i][0] = coords[0];
            valPoints[i][1] = coords[1];
            valPoints[i][2] = coords[2];
        }
        double ratio = tet_radius_ratio(valPoints);
        if (ratio < min_rratio_all)
            min_rratio_all = ratio;
        if(vertex_active_map.contains(vertexHash(vs))){
            if (vertex_active_map[vertexHash(vs)]){
                activeTet++;
                if (ratio < min_rratio_active)
                    min_rratio_active = ratio;
            }
        }
    });
    // save timing records
    save_timings("timings.json",time_label, profileTimer);
        // save statistics
    save_metrics("stats.json", tet_metric_labels, {(double)mesh.get_num_tets(), activeTet, min_rratio_all, min_rratio_active, (double)sub_call_two, (double) sub_call_three});
    //write mesh
    mtet::save_mesh("output.msh", mesh);
    
    return 0;
}
