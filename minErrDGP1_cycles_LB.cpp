#include "geometry.h"
#include "dgp.h"



int main(int argc, char* argv[]){
    if(argc != 2){
        std::cerr << "Usage of script: " << argv[0] << " <dgp_file_name> \n";
        return 1;
    }


    auto [edges, vertex_ids, vertex_id_to_name] = parse_dgp_instance_dat_file(argv[1]);
    std::unordered_map<int, std::vector<Adjacency>> adj_list = create_adjacency_list_from_edges(edges, vertex_ids.size());
    std::vector<std::vector<Edge>> spanning_tree_edges = get_spanning_forest(adj_list);
    std::vector<std::vector<std::vector<Edge>>> cycle_basis = get_cycle_basis(spanning_tree_edges, adj_list);

    double min_errDGP_cycle_LB = compute_minErrDGP_cycle_basis(cycle_basis);
    std::cerr << "got error " << min_errDGP_cycle_LB;
    std::cout << min_errDGP_cycle_LB;
}