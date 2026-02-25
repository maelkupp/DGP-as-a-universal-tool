#include "geometry.h"
#include "dgp.h"
#include <chrono>


int main(int argc, char* argv[]){
    if(argc != 3){
        std::cerr << "Usage of script: " << argv[0] << " <dgp_file_name> <are_there_real_weights> (1 for real edge weights 0 for integer edge weights)\n";
        return 1;
    }

    bool real_edge_weights = (std::string(argv[2]) == "1");

    auto [edges, vertex_ids, vertex_id_to_name] = parse_dgp_instance_dat_file(argv[1]);
    auto start = std::chrono::high_resolution_clock::now();
    std::unordered_map<int, std::vector<Adjacency>> adj_list = create_adjacency_list_from_edges(edges, vertex_ids.size());
    std::vector<std::vector<Edge>> spanning_tree_edges = get_spanning_forest(adj_list);
    std::cerr << "number of connected components " << std::size(spanning_tree_edges) << "\n";
    std::cerr << "number of edges in first spanning tree " << std::size(spanning_tree_edges[0]) << "\n";
    auto greedy_cycles = get_greedy_edge_disjoint_cycles(spanning_tree_edges, adj_list);
    double min_errDGP_cycle_LB = compute_minErrDGP_cycle(greedy_cycles, real_edge_weights);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cerr << "got error " << min_errDGP_cycle_LB << "\n";

    std::cout << min_errDGP_cycle_LB << " " << elapsed.count() << std::endl;
}