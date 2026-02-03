#include "dgp.h"
#include "geometry.h"
#include <iostream>
#include <vector>




int main(int argc, char* argv[]){
    // n: number of points, so number of vertices in the graph
    // p: density of the graph, between 0 and 1 (does not make that much sense for high values of p)
    // A: [0,A]x[0,A] the points will be generated in this square

    if (argc != 4){
        std::cerr << "Usage: " << argv[0] << " <number of vertices> <density of edges> <height/width of square>" << std::endl;
        //return error code 1
        return 1;
    }
    int n = std::stoi(argv[1]);
    double p = std::stod(argv[2]);
    double A = std::stod(argv[3]);

    auto [edges, points] = generate_2d_dgp_instance(n, p, A);
    write_edges_to_dat_file(edges);

    double UB_error = cheap_rotation_minErrDGP1_UB(edges, points);
    std::cerr << "Got cheap UB: " << UB_error << "\n";
    std::cout << UB_error;

    return 0;
}