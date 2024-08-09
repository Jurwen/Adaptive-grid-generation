//#define Check_Flip_Tets
//#include <mtet/mtet.h>
//#include <mtet/io.h>
//#include <ankerl/unordered_dense.h>
#include <span>
#include <queue>
#include <optional>
//#include <implicit_functions.h>
#include <CLI/CLI.hpp>

#include "io.h"
#include "timer.h"
#include "csg.h"
#include "grid_mesh.h"
#include "grid_refine.h"


using namespace mtet;

int main(int argc, const char *argv[])
{
    struct
    {
        std::string grid_file;
        std::string function_file;
        double threshold;
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
    app.add_option("grid", args.grid_file, "Initial grid file")->required();
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
    // Read initial grid
    mtet::MTetMesh grid;
    if (args.grid_file.find(".json") != std::string::npos){
        grid = grid_mesh::load_tet_mesh(args.grid_file);
        mtet::save_mesh("init.msh", grid);
        grid = mtet::load_mesh("init.msh");
    } else {
        grid = mtet::load_mesh(args.grid_file);
    }
    
    double max_elements = args.max_elements;
    if (max_elements < 0)
    {
        max_elements = std::numeric_limits<int>::max();
    }
    std::string function_file = args.function_file;
    double threshold = args.threshold;
    double alpha = args.alpha;
    int mode;
    double smallest_edge_length = args.smallest_edge_length;
    llvm_vecsmall::SmallVector<csg_unit, 20> csg_tree = {};
    
    if (args.method == "IA"){
        mode = IA;
    }
    if (args.method == "CSG"){
        mode = CSG;
        load_csgTree(args.csg_file, csg_tree);
    }
    if (args.method == "MI"){
        mode = MI;
    }
    if (args.curve_network){
        curve_network = true;
    }
    
    
    auto implicit_func = [&](){
        
    }
    
    ///
    /// the lambda funciton for csg tree iteration/evaluation.
    /// @param[in] funcInt          Given an input of value range std::array<double, 2> for an arbitrary number of functions
    /// @return   A value range of this CSG operation in a form of `std::array<double, 2>` and a list of active function in a form of    `llvm_vecsmall::SmallVector<int, 20>>`
    ///
    auto csg_func = [&](llvm_vecsmall::SmallVector<std::array<double, 2>, 20> funcInt){
        if (args.csg_file == ""){
            throw std::runtime_error("ERROR: no csg file provided");
            std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> null_csg = {{},{}};
            return null_csg;
        }else{
            return iterTree(csg_tree, 1, funcInt);
        }
    };
    
    //perform main grid refinement algorithm:
    tet_metric metric_list;
    //an array of 8 {total time getting the multiple indices, total time,time spent on single function, time spent on double functions, time spent on triple functions time spent on double functions' zero crossing test, time spent on three functions' zero crossing test, total subdivision time, total evaluation time,total splitting time}
    std::array<double, timer_amount> profileTimer = {0,0,0,0,0,0,0,0,0,0};
    if (!gridRefine(function_file, mode, threshold, alpha, max_elements, csg_func, grid, metric_list, profileTimer))
    {
        throw std::runtime_error("ERROR: unsuccessful grid refinement");
    }
    // save timing records
    save_timings("timings.json",time_label, profileTimer);
    //profiled time(see details in time.h) and profiled number of calls to zero
    for (int i = 0; i < profileTimer.size(); i++){
        timeProfileName time_type = static_cast<timeProfileName>(i);
        std::cout << time_label[i] << ": " << profileTimer[i] << std::endl;
    }
    // save tet metrics
    save_metrics("stats.json", tet_metric_labels, metric_list);
    return 0;
}
