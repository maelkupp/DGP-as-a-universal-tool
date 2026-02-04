#ifndef DGP_H
#define DGP_H
#include "geometry.h"
#include <vector>
#include <tuple>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <limits>

struct Edge{
    int u; //id of 1st vertex so if there is a directed edge this is always the source node
    int v; //if of 2nd vertex
    double weight; //the weight of the edge
    bool directed; //whether it is a directed edge


    Edge() = default;
    //simple constructor for now
    Edge(int u, int v, double weight, bool directed): u(u), v(v), weight(weight), directed(directed){}

    void display() const {
        std::string directed_string = "";
        if(directed){
            directed_string = "Directed";
        }
        std::cerr << directed_string << "Edge between vertices u: " << u << " and v: " << v << " of weight " << weight << "\n";
    }

};

struct Adjacency {
    int neighbourId;
    double dist;
    bool directed;   // is it a direction constraint? if directed if false we can ignore outgoing and do not initialize it for efficiency
    bool outgoing;   // is the edge outgoing (true) or incoming (false) with respect to the neighbour of this adjacent node

    //if there are three parameters to the contstructor we know it is directed
    Adjacency(int neighbourId, double dist, bool outgoing): neighbourId(neighbourId), dist(dist), directed(true), outgoing(outgoing){} 

    //if two parameters we know it is undirected edge and so do not initialise the outgoing attribute
    Adjacency(int neighbourId, double dist):  neighbourId(neighbourId), dist(dist), directed(false) {}

    void display(){
        std::string directed_string = "";
        if(directed){
            directed_string += "Directed";
        }
        std::cerr << directed_string << " Neighbour with id " << neighbourId << ", " <<  dist << " away \n";
    }

};

struct BacktrackState{
    double best_error;
    std::unordered_map<int, double> best_embedding;

    BacktrackState(): best_error(std::numeric_limits<double>::infinity()) {}
};

// struct that allows me to hash a pair of integers
struct PairHash {
    std::size_t operator()(const std::pair<int,int>& p) const {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};
// -------- functions to create/handle DGP instances ----------

std::tuple<std::vector<Edge>, std::set<int>, std::unordered_map<int, std::string>> parse_dgp_instance_dat_file(std::string file_path);

std::unordered_map<int, std::vector<Adjacency>> create_adjacency_list_from_edges(const std::vector<Edge>& edges,const int num_vertices);

std::vector<std::vector<Edge>> get_spanning_forest(std::unordered_map<int, std::vector<Adjacency>>& adj_list);

std::vector<std::vector<std::vector<Edge>>> get_cycle_basis(std::vector<std::vector<Edge>>& spanning_forest_edges, std::unordered_map<int, std::vector<Adjacency>>& adj_list);

double compute_minErrDGP_cycle_basis(std::vector<std::vector<std::vector<Edge>>>& cycle_basis, bool real_edge_weightsDP_cycle_error);

void display_1Dembedding(std::unordered_map<int, double> embedding, std::unordered_map<int, std::string> vertex_id_2_name);

std::set<int> dfs_through_adj_list(std::unordered_map<int, std::vector<Adjacency>>& adj_list, int id);

std::vector<std::set<int>> get_connected_components(std::unordered_map<int, std::vector<Adjacency>>& adj_list, const std::set<int>& vertex_ids);

double compute_vertex_embedding_error(int vertex_id_to_place, std::unordered_map<int, double> current_embedding, 
    std::unordered_map<int, std::vector<Adjacency>> adj_list, std::unordered_map<int, int>& pos);

double verbose_compute_1Dembedding_error(std::vector<Edge> edges, std::unordered_map<int, double> embedding, std::unordered_map<int, std::string> vertex_id_to_name);
void write_edges_to_dat_file(const std::vector<Edge>& edges);


// ------------- obtained Upper Bounds from geometry ----------

double cheap_rotation_minErrDGP1_UB(const std::vector<Edge>& edges, const std::vector<Point>& points);
double optimized_rotation_minErrDGP1_UB(const std::vector<Edge>& edges, const std::vector<Point>& points);


// ------------ things to do with 2D instances and generating them etc ------------
std::pair<std::vector<Edge>, std::vector<Point>> generate_2d_dgp_instance(int n, double p, double A);
std::vector<Point> generate_points(int n, double A);
std::vector<Edge> generate_edges_across_points(std::vector<Point>& points, double p);



// ------------ relaxations

double simple_projection_relax(const std::vector<Edge>&, const std::vector<Point>& points);
double rotation_projection_relax(const std::vector<Edge>&, std::vector<Point>& points);
#endif