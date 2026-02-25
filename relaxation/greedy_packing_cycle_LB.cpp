#include "geometry.h"
#include "dgp.h"
#include <chrono>


int main(int argc, char* argv[]){
    if(argc != 2){
        std::cerr << "Usage of script: " << argv[0] << " <dgp_file_name>\n";
        return 1;
    }

    auto [edges, vertex_ids, vertex_id_to_name] = parse_dgp_instance_dat_file(argv[1]);

    auto start = std::chrono::high_resolution_clock::now();
    std::unordered_map<int, std::vector<Adjacency>> adj_list = create_adjacency_list_from_edges(edges, vertex_ids.size());
    std::vector<std::vector<Edge>> spanning_tree_edges = get_spanning_forest(adj_list);
    std::cerr << "number of connected components " << std::size(spanning_tree_edges) << "\n";
    std::cerr << "number of edges in first spanning tree " << std::size(spanning_tree_edges[0]) << "\n";
    auto cycle_basis = get_cycle_basis(spanning_tree_edges, adj_list);


    std::vector<std::vector<double>> cycle_errors; //represents the error of each cycle
    cycle_errors.resize(cycle_basis.size());
    for(int i=0; i<static_cast<int>(cycle_basis.size()); ++i){
        cycle_errors[i].resize(cycle_basis[i].size());
        for(int j=0; j<static_cast<int>(cycle_basis[i].size()); ++j){
            std::vector<double> cycle_values {};
            for(const Edge& e: cycle_basis[i][j]){
                cycle_values.push_back(e.weight);
            }
            cycle_errors[i][j] = DP_cycle_error(cycle_values);
        }
    }

    double greedy_packing_LB = greedy_packing_cycle_err(cycle_basis, cycle_errors);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cerr << "got error " << greedy_packing_LB << "\n";
    std::cout << greedy_packing_LB << " " << elapsed.count() << std::endl;
}