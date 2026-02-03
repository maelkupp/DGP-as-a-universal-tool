#include "geometry.h"
#include "dgp.h"
#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <tuple>
#include <cmath>
#include <format>
#include <string>
#include <limits>
#include <iomanip>



std::vector<std::string> generate_vertex_names(int n){
    std::vector<std::string> names;
    for(int i=0; i<n; ++i){
        names.push_back(std::format("x_{}", i)); //ignore
    }
    return names;
}



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


    double tot_err = simple_projection_relax(edges, points);
    double rot_err = rotation_projection_relax(edges, points);
    std::cerr << "tot_err: " << tot_err << std::endl;
    std::cerr << "rot_err: " << rot_err << std::endl;
    std::cout << std::fixed << std::setprecision(10)
    <<rot_err << "," << tot_err;
    return 0;
}
