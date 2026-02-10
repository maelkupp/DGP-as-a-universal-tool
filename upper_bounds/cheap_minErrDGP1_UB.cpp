#include "dgp.h"
#include "geometry.h"
#include <iostream>
#include <vector>




int main(int argc, char* argv[]){

    if(argc != 2){
        std::cout << "Usage: " << argv[0] << " dgp_filename.dat \n";
        std::cout << "This script is used to get as tight an upper bound as possible using discrete techniques on MinErrDGP1 of a DGP instance \n";
        return 1;
    }

    auto [edges, vertex_ids, vertex_id_to_name] = parse_dgp_instance_dat_file(argv[1]);
    double UB_optimized_error = optimized_projection_minErrDGP1_UB(edges, vertex_ids);
    std::cerr << "Optimized error " << UB_optimized_error << "\n";
    std::cout << UB_optimized_error;

    return 0;
}