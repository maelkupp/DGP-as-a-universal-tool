#include "dgp_bnb.h"
#include "random.h"
//#include "gurobi_c++.h"
#include <tuple>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>
#include <cmath>
#include <random>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <numeric>


// 1 ---------------------- Graph + Preprocessing

std::tuple<std::vector<Edge>, std::set<int>, std::unordered_map<int, std::string>> parse_dgp_instance_dat_file(std::string dgp_dat_file_path){
    std::vector<Edge> edges;
    std::set<int> vertex_ids;
    std::unordered_map<int, std::string> vertex_id_to_name;

    std::ifstream infile(dgp_dat_file_path);
    if (!infile) {
        std::cerr << "Error opening file: " << dgp_dat_file_path << "\n";
        return {edges, vertex_ids, vertex_id_to_name};
    }

    std::string line;
    bool in_edge_block = false;
    

    while (std::getline(infile, line)) {

        // trim leading whitespace
        auto first = line.find_first_not_of(" \t");
        if (first == std::string::npos) continue;
        line = line.substr(first);

        // skip full-line comments
        if (line.starts_with("#")) continue;

        // start of edge block
        if (line.starts_with("param : E")) {
            in_edge_block = true;
            continue;
        }

        // end of edge block
        if (in_edge_block && line == ";") {
            in_edge_block = false;
            continue;
        }

        if (!in_edge_block) {
            // ignore param lines for now
            continue;
        }

        // ----- EDGE LINE -----

        // split comment
        std::string data_part = line;
        std::string comment_part;

        size_t hash_pos = line.find('#');
        if (hash_pos != std::string::npos) {
            data_part = line.substr(0, hash_pos);
            comment_part = line.substr(hash_pos + 1);
        }

        // parse numeric part
        std::stringstream ss(data_part);

        int u, v, flag;
        double dist;

        if (!(ss >> u >> v >> dist >> flag)) {
            std::cerr << "Malformed edge line: " << line << "\n";
            continue;
        }

        edges.emplace_back(std::min(u,v), std::max(u,v), dist, flag);
        vertex_ids.insert(u);
        vertex_ids.insert(v);

        // parse names from comment if present
        if (!comment_part.empty()) {

            size_t lbracket = comment_part.find('[');
            size_t rbracket = comment_part.find(']');

            if (lbracket != std::string::npos && rbracket != std::string::npos) {

                std::string names =
                    comment_part.substr(lbracket + 1,
                                         rbracket - lbracket - 1);

                size_t comma = names.find(',');
                if (comma != std::string::npos) {

                    std::string u_name = names.substr(0, comma);
                    std::string v_name = names.substr(comma + 1);

                    if (!vertex_id_to_name.count(u))
                        vertex_id_to_name[u] = u_name;
                    if (!vertex_id_to_name.count(v))
                        vertex_id_to_name[v] = v_name;
                }
            }
        }
    }
    return {edges, vertex_ids, vertex_id_to_name};
}

std::unordered_map<int, std::vector<Adjacency>> create_adjacency_list_from_edges(const std::vector<Edge>& edges, const int num_vertices){
    //function that creates an adjacency list of a DGP instance from a vector of Edges
    std::cerr << "Creating an adjacency list for a DGP instance with " << edges.size() << "edges \n";
    std::unordered_map<int, std::vector<Adjacency>> adj_list;

    for(int i=1; i<=num_vertices; ++i){
        //initialise each adjacency list of the vertices to the empty list
        adj_list[i] = {};
    }

    for(const Edge& edge: edges){
        if(edge.directed){
            adj_list[edge.u].push_back(Adjacency(edge.v, edge.weight, true));
            adj_list[edge.v].push_back(Adjacency(edge.u, edge.weight, false));
        }else{
            adj_list[edge.u].push_back(Adjacency(edge.v, edge.weight));
            adj_list[edge.v].push_back(Adjacency(edge.u, edge.weight));    
        }
    }

    return adj_list;
}


std::vector<IndexedEdge> build_indexed_edges(
    const std::vector<Edge>& edges
){
    std::vector<IndexedEdge> indexed_edges;
    int n {static_cast<int>(edges.size())};
    indexed_edges.reserve(n);
    for(int i=0; i<n; ++i){
        //make the edges 0 indexed (may come back to this and change it)
        indexed_edges.push_back({i, edges[i].u, edges[i].v, edges[i].weight});
    }

    //creates a vector so that indexed_edges[i] has id i
    return indexed_edges;
}


std::pair<std::vector<Edge>, std::unordered_set<int>> dfs_for_spanning_tree(
    std::unordered_map<int, std::vector<Adjacency>>& adj_list, 
    int start_vertex
){
    //returns a pair, a vector of Edges in the spanning tree of the connected component and a set of all the vertices it visitied during the DFS
    std::unordered_set<int> visited;
    std::vector<Edge> tree_edges;

    std::stack<std::tuple<int, int, double>> stk; // (parent, current)
    stk.push({-1, start_vertex, 0.0});
    visited.insert(start_vertex);

    while (!stk.empty()) {
        auto [parent, u, dist] = stk.top();
        stk.pop();

        if (parent != -1) {
            // find distance from parent -> u, add edge to the vector of edges of the spanning tree
            tree_edges.push_back({parent, u, dist, 0});

        }
        
        //look through each neighbour in the adjacency of u, if we have not yet visited the node, we mark it as visited and and the edge to the stack
        for (const Adjacency& adj : adj_list[u]) {
            int v = adj.neighbourId;
            if (visited.find(v) == visited.end()) {
                visited.insert(v);
                stk.push({u, v, adj.dist});
            }
        }
    }

    return {tree_edges, visited};
}

std::vector<std::vector<Edge>> get_spanning_forest(std::unordered_map<int, std::vector<Adjacency>>& adj_list){
    //this function returns a vector of edges (all assumed to be undirected) corresponding to a spanning tree of a graph represented with its
    //adjacency list, this is done by performing a DFS through the tree
    //if there are multiple connected components in our graph this will return a spanning forest
    
    std::vector<std::vector<Edge>> spanning_forest {};
    std::unordered_set<int> visited {};

    for(const auto& [v_id, list]: adj_list){
        //looks if we have visited this vertex, if we have not, we are looking at a new connected component of the graph
        if(visited.find(v_id) == visited.end()){
            //we have not visited the ID yet
            auto tree_pair = dfs_for_spanning_tree(adj_list, v_id);
            std::vector<Edge> spanning_tree = tree_pair.first;
            std::unordered_set<int> visited_tree = tree_pair.second;


            //update the visited set, sets all the vertices in the connected component as visited
            spanning_forest.push_back(spanning_tree);
            for(const auto& visited_id: visited_tree){
                visited.insert(visited_id);
            }
        }
    }
    return spanning_forest;
}

//function that creates a standardised edge set for an edge (u,v) -> (min(u,v), max(u,v))
std::unordered_set<std::pair<int, int>, PairHash> create_tree_edge_set(const std::vector<Edge>& spanning_tree_edges){
    std::unordered_set<std::pair<int, int>, PairHash> tree_edges_set{};
    for(const auto& edge: spanning_tree_edges){
        tree_edges_set.insert({std::min(edge.u, edge.v), std::max(edge.u, edge.v)});
    }
    return tree_edges_set;
}

std::pair<std::unordered_map<int, Edge>,std::unordered_map<int, int>>get_parent_and_depth_maps_DFS(
    std::unordered_map<int, std::vector<Adjacency>>& adj_list,
    int root_id
){
    // parent_edge[v] = edge connecting v to its parent in the DFS tree
    std::unordered_map<int, Edge> parent_edge;
    //depth[v] = depth of the vertex v in the DFS tree
    std::unordered_map<int, int> depth;

    depth[root_id] = 0;
    // root has no parent edge

    std::stack<int> stk;
    std::unordered_set<int> visited;

    visited.insert(root_id);
    stk.push(root_id);

    while (!stk.empty()) {
        int u = stk.top();
        stk.pop();

        for (const auto& adj : adj_list[u]) {
            int v = adj.neighbourId;

            if (visited.find(v) == visited.end()) {
                visited.insert(v);
                stk.push(v);

                // store the actual parent edge
                parent_edge[v] = Edge(u, v, adj.dist, 0);
                depth[v] = depth[u] + 1;
            }
        }
    }

    return {parent_edge, depth};
}

std::vector<std::vector<std::vector<Edge>>> get_cycle_basis(std::vector<std::vector<Edge>>& spanning_forest_edges, std::unordered_map<int, std::vector<Adjacency>>& adj_list){
    //this function will take a spanning tree and the adjacency list of the entire graph and will return a cycle basis of the entire graph
   std::vector<std::vector<std::vector<Edge>>> cycle_basis;
   //adjacency list representing only the spanning tree


   for(const auto& spanning_tree_edges: spanning_forest_edges){
    //preprocessing for each tree in the forest
    std::unordered_set<std::pair<int, int>, PairHash> tree_edge_set = create_tree_edge_set(spanning_tree_edges);
    if(spanning_tree_edges.size() == 0) continue;
    
    //don't think I need this tree_adj_list, everything should work fine with the original adj_list, only need to be sure to put a vertex of the correct c.c
    //but that is done with spanning_tree_edges[0].u, need to come back to this and check it
    std::unordered_map<int, std::vector<Adjacency>> tree_adj_list = create_adjacency_list_from_edges(spanning_tree_edges, spanning_tree_edges.size());
    auto [parent_edge_tree, depth_tree] = get_parent_and_depth_maps_DFS(tree_adj_list, spanning_tree_edges[0].u); // get an arbitrary vertex id that is in the tree from the first element in the vector
    std::vector<std::vector<Edge>> component_cycle_basis;


    for(auto& it: adj_list){
        //iterating over all edges of the graph
        for(auto& adj: adj_list[it.first]){
            //iterating over each neighbour
            
            //saving steps and duplicates
            if(adj.neighbourId >= it.first) continue;
                
            // the edge is not in this connected component
            if(parent_edge_tree.find(it.first) == parent_edge_tree.end() || parent_edge_tree.find(adj.neighbourId) == parent_edge_tree.end()) continue;             
            
            //the edge is a tree edge so we skip it
            if(tree_edge_set.find({it.first, adj.neighbourId}) != tree_edge_set.end()) continue;
                
            //now dealing with a non-tree edge for the correct spanning tree of the connected component

            std::vector<Edge> cycle {}; // initialise an empty cycle
                int x = it.first;
                int y = adj.neighbourId;

                // Lift deeper vertex until depths match
                while (depth_tree[x] > depth_tree[y]) {
                    Edge e = parent_edge_tree[x];
                    cycle.push_back(e);
                    x = (e.u == x ? e.v : e.u);
                }

                while (depth_tree[y] > depth_tree[x]) {
                    Edge e = parent_edge_tree[y];
                    cycle.push_back(e);
                    y = (e.u == y ? e.v : e.u);
                }

                // Lift both simultaneously until lowest common ancestor is reached
                while (x != y) {
                    Edge ex = parent_edge_tree[x];
                    Edge ey = parent_edge_tree[y];
                    cycle.push_back(ex);
                    cycle.push_back(ey);
                    x = (ex.u == x ? ex.v : ex.u);
                    y = (ey.u == y ? ey.v : ey.u);
                }

                // Close the cycle with the non-tree edge ----
                cycle.push_back(Edge(it.first, adj.neighbourId, adj.dist, 0));

                component_cycle_basis.push_back(cycle);
        }
    }
    if(component_cycle_basis.size() > 0){
        cycle_basis.push_back(component_cycle_basis);

    }

   }
   return cycle_basis;
}

//don't really see the point of this function for now - will implement it later if I actually need it
std::vector<CycleID> convert_cycles_to_edge_ids(
    const std::vector<std::vector<std::vector<Edge>>>& cycle_basis,
    const std::vector<IndexedEdge>& indexed_edges
){
    std::vector<CycleID> indexed_cycles;
    std::unordered_map<std::pair<int, int>, int, PairHash> edge_vertices_to_index_map;

    //perform a unique mapping of the vertices of the edge to the edge index by doing (2^{edge.u} x 3^{edge.v} = edge.id), unique by fundemental theory of arithmetic
    for(const auto& ie: indexed_edges){
        std::pair<int, int> p_key = std::make_pair(std::min(ie.u, ie.v), std::max(ie.u, ie.v));
        edge_vertices_to_index_map[p_key] = ie.id;
    }

    std::vector<int> edge_ids {};
    //iterate over each connected component
    for(const auto& cc_basis: cycle_basis){
        //iterate of each cycle in the cycle basis for this connected component
        for(const auto& cycle: cc_basis){
            edge_ids.resize(static_cast<int>(cycle.size()));
            for(const auto& e: cycle){
                std::pair<int, int> p_key = std::make_pair(std::min(e.u, e.v), std::max(e.u, e.v));
                edge_ids.push_back(edge_vertices_to_index_map[p_key]);
            }
            indexed_cycles.push_back({edge_ids});
        }
    }

    return indexed_cycles;
}



// 2 ---------------- Relaxation + Lower Bound system ---------------------



double DP_cycle_error(const std::vector<double>& cycle) {
    int n = cycle.size();
    if (n == 0) return 0.0;

    int S = std::accumulate(cycle.begin(), cycle.end(), 0);

    std::vector<char> dp(S + 1, 0);
    dp[0] = 1;

    for (int w : cycle) {
        for (int j = S - w; j >= 0; --j) {
            if (dp[j]) dp[j + w] = 1;
        }
    }

    int target = S / 2;
    int best = 0;

    for (int p = 0; p <= S; ++p) {
        if (!dp[p]) continue; //if we cannot reach that value with the edge weight we continue
        if (std::abs(p - target) < std::abs(best - target)) {
            best = p;
        }
    }

    double cycle_error {std::fabs(2 * best - S)};
    std::cerr << "Got additional error on a cycle " << cycle_error << "\n";
    return cycle_error;
}

double DP_cycle_error(
    const CycleID& cycle,
    const std::vector<IndexedEdge>& i_edges
){
    int n = cycle.edge_ids.size();
    if (n == 0) return 0.0;
    //capturing i_edges only and by reference
    int S = std::accumulate(cycle.edge_ids.begin(), cycle.edge_ids.end(), 0,
            [&i_edges](int acc, int edge_id){
                return acc + i_edges[edge_id].weight;
            }
        );


    std::vector<char> dp(S + 1, 0);
    dp[0] = 1;

    for (int edge_id : cycle.edge_ids) {
        for (int j = S - i_edges[edge_id].weight; j >= 0; --j) {
            if (dp[j]) dp[j + i_edges[edge_id].weight] = 1;
        }
    }

    int target = S / 2;
    int best = 0;

    for (int p = 0; p <= S; ++p) {
        if (!dp[p]) continue; //if we cannot reach that value with the edge weight we continue
        if (std::abs(p - target) < std::abs(best - target)) {
            best = p;
        }
    }

    double cycle_error {std::fabs(2 * best - S)};
    std::cerr << "Got additional error on a cycle " << cycle_error << "\n";
    return cycle_error;
    

}



double DP_cycle_error_with_fixed_signs(
    const CycleID& cycle,
    const std::vector<int>& fixed_signs,
    const std::vector<IndexedEdge>& edges
){
    double fixed_sum {0.0};
    int n = static_cast<int>(fixed_signs.size());
    std::vector<double> free_weights {};


    for(int i=0; i<n; ++i){
        if(fixed_signs[i] != 0){
            fixed_sum += fixed_signs[i]*edges[i].weight;
        }else{
            free_weights.push_back(edges[i].weight);
        }
    }

    double S = std::accumulate(free_weights.begin(), free_weights.end(), 0);
    std::vector<char> dp(S + 1, 0);
    dp[0] = 1;

    for(auto w: free_weights){
        for(int j=S-w; j>=0; --j){
            if (dp[j]) dp[j + w] = 1;
        }
    }

    int target = S / 2;
    int best = 0;

    for (int p = 0; p <= S; ++p) {
        if (!dp[p]) continue; //if we cannot reach that value with the edge weight we continue
        if (std::abs(fixed_sum + p - target) < std::abs(fixed_sum + best - target)) {
            best = p;
        }
    }

    double cycle_error {std::fabs(2 * best - S)};
    std::cerr << "Got additional error on a cycle " << cycle_error << "\n";
    return cycle_error;
}

std::vector<double> get_cycle_basis_error(
    const std::vector<CycleID>& indexed_cycles,
    const std::vector<IndexedEdge>& i_edges
){
    std::vector<double> cycle_basis_error(static_cast<int>(indexed_cycles.size()), 0); //we already know the size
    for(int i=0; i<static_cast<int>(indexed_cycles.size()); ++i){
        cycle_basis_error[i] = DP_cycle_error(indexed_cycles[i], i_edges);
    }

    return cycle_basis_error;
}

//this function takes a cycle basis and the associated error (weight) of that cycle, the point is to obtain the best possible error
//by greedily selecting the cycle that are edge disjoint and have the greatest error
double greedy_packing_cycle_err(const std::vector<std::vector<std::vector<Edge>>>& cycle_basis,
    const std::vector<std::vector<double>>& cycle_errors){

    //first flatten the cycle into 1dimension
    std::vector<WeightedCycle> all_cycles;
    for(int i=0; i<static_cast<int>(cycle_basis.size()); ++i){
        for(int j=0; j<static_cast<int>(cycle_basis[i].size()); ++j){
            all_cycles.push_back({cycle_basis[i][j], cycle_errors[i][j]});
        }
    }

    //sort in decreasing order by weight (to have the cycles contributing to the error the most at the beginning)
    std::sort(all_cycles.begin(), all_cycles.end(),
            [](const WeightedCycle& a, const WeightedCycle& b) {
                return a.error > b.error;
            });

    std::set<std::pair<int,int>> used_edges;
    double total_lb = 0.0;

    for (const auto& cycle : all_cycles) {

        bool disjoint = true;

        for (const Edge& e : cycle.edges) {
            auto key = std::make_pair(std::min(e.u, e.v),
                                    std::max(e.u, e.v));
            if (used_edges.count(key)) {
                disjoint = false;
                break;
            }
        }

        if (disjoint) {
            total_lb += cycle.error;

            for (const Edge& e : cycle.edges) {
                std::pair<int, int> key = std::make_pair(std::min(e.u, e.v), std::max(e.u, e.v));
                used_edges.insert(key);
            }
        }
    }

    return total_lb;
    
}


double compute_cycle_packing_LB(
    const std::vector<CycleID>& cycles,
    const std::vector<int>& fixed_signs,
    const std::vector<IndexedEdge>& edges
){

    int n = static_cast<int>(cycles.size());
    std::vector<WeightedCycleID> vec_weighted_cycle;
    vec_weighted_cycle.resize(n);
    for(int i=0; i<n; ++i){
        vec_weighted_cycle[i] = {i, DP_cycle_error_with_fixed_signs(cycles[i], fixed_signs, edges), cycles[i].edge_ids};
    }

    //sort in decreasing order by weight (to have the cycles contributing to the error the most at the beginning)
    std::sort(vec_weighted_cycle.begin(), vec_weighted_cycle.end(),
            [](const WeightedCycleID& a, const WeightedCycleID& b) {
                return a.error > b.error;
            });

    std::set<std::pair<int, int>> used_edges;
    double total_lb = 0.0;

    for (const auto& cycle : vec_weighted_cycle) {

        bool disjoint = true;

        for (int e_id : cycle.edge_ids) {
            std::pair<int, int> p_key = std::make_pair(std::min(edges[e_id].u, edges[e_id].v), std::max(edges[e_id].u, edges[e_id].v));
            if (used_edges.count(p_key)) {
                disjoint = false;
                break;
            }
        }

        if (disjoint) {
            total_lb += cycle.error;

            for (int e_id : cycle.edge_ids) {
                std::pair<int, int> p_key = std::make_pair(std::min(edges[e_id].u, edges[e_id].v), std::max(edges[e_id].u, edges[e_id].v));
                used_edges.insert(p_key);
            }
        }
    }

    return total_lb;
}


// 3 ------------------ Upper bound system -----------------------------

double local_obj_abs(
    double x,
    int i_vertex_id_1based,
    const std::unordered_map<int, std::vector<Adjacency>>& adj_list,
    const std::vector<double>& t
){
    double sum = 0.0;
    const auto& nbrs = adj_list.at(i_vertex_id_1based);
    for (const auto& nbr : nbrs) {
        int j = nbr.neighbourId - 1;   // 0-based index into t
        sum += std::abs(std::abs(x - t[j]) - nbr.dist);
    }
    return sum;
}


double coefficient_descent_on_line(const std::vector<Edge>& edges,
    const std::unordered_map<int, std::vector<Adjacency>>& adj_list,
    std::vector<double>& t
){  

    //callable lambda to obtain the error of an embedding
    auto global_error = [&](const std::vector<double>& tt){
        double E = 0.0;
        for (const auto& e : edges) {
            E += std::abs(std::abs(tt[e.u - 1] - tt[e.v - 1]) - e.weight);
        }
        return E;
    };

    double prev_E = global_error(t);
    std::cerr << "Initial 1D error: " << prev_E << "\n";

    //best constants I have found with some testing, the testing what not very rigorous, can be improved upon
    const double lambda = 0.5;
    const double rel_tol = 5e-5;
    const int max_sweeps = 100; //still have a max sweep to prevent non - convergence


    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        // Jacobi-style buffer for stability
        std::vector<double> t_new = t;
        
        //have 1 -indexed for the vertex ids
        for (int i = 1; i < static_cast<int>(t.size()); ++i) {
            int vid = i-1;
            auto it = adj_list.find(i);
            if (it == adj_list.end() || it->second.empty()) continue; //if for some reason the vertex is not in the adjacency list keys or it does not have any neighbours continue to the next vertex

            // Build candidate set: {t_j ± d_ij} plus current t_i
            std::vector<double> cand;
            cand.reserve(2 * it->second.size() + 1);
            cand.push_back(t[i]);

            for (const auto& nbr : it->second) {
                int j = nbr.neighbourId - 1;
                double dij = nbr.dist;
                cand.push_back(t[j] + dij);
                cand.push_back(t[j] - dij);
            }

            // Pick best candidate for the TRUE local objective
            double best_x = cand[0];
            double best_val = local_obj_abs(best_x, vid, adj_list, t);
            for (int k = 1; k < (int)cand.size(); ++k) {
                double val = local_obj_abs(cand[k], vid, adj_list, t);
                if (val < best_val) {
                    best_val = val;
                    best_x = cand[k];
                }
            }

            // Damped update
            t_new[i] = (1.0 - lambda) * t[i] + lambda * best_x; //move part of the way towards the locally optimal position
            //if we jump all the way we might have strong reactions to conflicting constraints etc, often ping - pong back and forth
        }

        t.swap(t_new); // same as t = t_new but faster

        double curr_E = global_error(t);
        double rel_impr = std::abs(prev_E - curr_E) / std::max(1.0, prev_E);

        if (rel_impr < rel_tol) {
            std::cerr << "Converged after " << sweep + 1 << " sweeps \n";
            break;
        }

        prev_E = curr_E;


        std::cerr << "After sweep: " << curr_E << "\n";
    }

    double Efinal = global_error(t);
    std::cerr << "Optimized error: " << Efinal << "\n";
    return Efinal;
}


static inline double wrap_angle_pi(double theta) {
    // For line directions in 1D projection, angles are equivalent mod pi
    // because u and -u define the same line. Keep theta in [0, pi).
    double pi = PI;
    theta = std::fmod(theta, pi);
    if (theta < 0) theta += pi;
    return theta;
}

double optimized_projection_minErrDGP1_UB(
    const std::vector<Edge>& edges,
    const std::set<int>& vertex_ids,
    std::unordered_map<int, std::vector<Adjacency>>& adj_list
){
    //this function creates a UB for the MinerrDGP1 using random scatterings of the points onto the line and then performing the coefficient descent algorithm
    // Build adjacency list once

    const int n = static_cast<int>(vertex_ids.size());
    if (n == 0) return 0.0;

    // ---- estimate a reasonable scale from the edges ----
    double avg_d = 0.0;
    for (const auto& e : edges) avg_d += e.weight;
    avg_d /= std::max(1, (int)edges.size());

    const int num_trials = 10;          // number of random restarts
    const double spread = 2.0 * avg_d;  // scale of random scattering
    const double noise = 0.1 * avg_d;   // small perturbation

    double best_error = std::numeric_limits<double>::infinity();

    // Convert vertex_ids to vector for indexing
    std::vector<int> vids(vertex_ids.begin(), vertex_ids.end());

    for (int trial = 0; trial < num_trials; ++trial) {

        // ---- initialize t ----
        std::vector<double> t(n, 0.0);

        // Pick a random root
        int root_idx = Random::get(0, n - 1);
        int root_vid = vids[root_idx];
        t[root_idx] = 0.0;

        // Simple BFS-style initialization
        std::unordered_map<int, bool> visited;
        visited[root_vid] = true;

        std::vector<int> stack;
        stack.push_back(root_vid);

        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();

            int u_idx = std::distance(vids.begin(),
                                      std::find(vids.begin(), vids.end(), u));

            for (const auto& nbr : adj_list[u]) {
                int v = nbr.neighbourId;
                if (visited[v]) continue;

                int v_idx = std::distance(vids.begin(),
                                          std::find(vids.begin(), vids.end(), v));

                // random sign
                double sign = (Random::get(0, 1) == 0) ? -1.0 : 1.0;
                t[v_idx] = t[u_idx] + sign * nbr.dist;

                // small noise to avoid symmetry lock-in
                t[v_idx] += Random::get<double>(-noise, noise);

                visited[v] = true;
                stack.push_back(v);
            }
        }

        // For disconnected vertices, scatter randomly, should change this to scatter points of each connected component, here we are only doing it for one connected component
        for (int i = 0; i < n; ++i) {
            if (!visited[vids[i]]) {
                t[i] = Random::get<double>(-spread, spread);
            }
        }

        // ---- run the 1D solver ----
        double err = coefficient_descent_on_line(edges, adj_list, t);

        best_error = std::min(best_error, err);
    }

    return best_error;
}



// 4 ------------------ Exact embedding solver -------------

std::vector<std::vector<double>> build_incidence_matrix(
    const std::vector<IndexedEdge>& edges,
    int n
){

    int m = edges.size();

    std::vector<std::vector<double>> B(
        m,
        std::vector<double>(n, 0.0)
    );

    for (int e = 0; e < m; ++e) {
        int u = edges[e].u - 1;  // assuming 1-based input
        int v = edges[e].v - 1;

        B[e][u] = 1.0;
        B[e][v] = -1.0;
    }

    return B;
}



//wrapper that does all the steps in computing the exact embedding that minimizes the error for a fixed sign, is what guarantees the correctness of the BnB
double solve_exact_embedding(
    const std::vector<int>& fixed_signs,
    const std::vector<IndexedEdge>& edges,
    int n_vertices
){
    return 0.0;
}
/*
double solve_exact_embedding(
    const std::vector<int>& fixed_signs,
    const std::vector<IndexedEdge>& edges,
    int n_vertices
) {
    try {

        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, 0);
        env.start();

        GRBModel model(env);

        int m = edges.size();

        // ------------------------
        // Variables x
        // ------------------------

        std::vector<GRBVar> x(n_vertices);

        for (int i = 0; i < n_vertices; ++i) {
            x[i] = model.addVar(
                -GRB_INFINITY,
                GRB_INFINITY,
                0.0,
                GRB_CONTINUOUS
            );
        }

        // Fix translation
        model.addConstr(x[0] == 0.0);

        // ------------------------
        // Slack variables t_e
        // ------------------------

        std::vector<GRBVar> t(m);

        for (int e = 0; e < m; ++e) {
            t[e] = model.addVar(
                0.0,
                GRB_INFINITY,
                1.0, // objective coefficient
                GRB_CONTINUOUS
            );
        }

        model.update();

        // ------------------------
        // Constraints
        // ------------------------

        for (int e = 0; e < m; ++e) {

            int u = edges[e].u - 1;
            int v = edges[e].v - 1;

            double b = fixed_signs[e] * edges[e].weight;

            // t_e >= (x_u - x_v) - b
            model.addConstr(
                t[e] >= x[u] - x[v] - b
            );

            // t_e >= -(x_u - x_v) + b
            model.addConstr(
                t[e] >= -x[u] + x[v] + b
            );
        }

        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        model.optimize();

        if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
            return std::numeric_limits<double>::infinity();
        }

        return model.get(GRB_DoubleAttr_ObjVal);

    } catch (GRBException& e) {
        std::cerr << "Gurobi error: "
                  << e.getMessage() << "\n";
        return std::numeric_limits<double>::infinity();
    }
}
*/
int select_branching_edge(
    const std::vector<int>& fixed_signs
){
    //for now implement this simple strategy where we simply iterate over each edge, simply give the next edge in the vector
    for(int i=0; i<static_cast<int>(fixed_signs.size()); ++i){
        if(fixed_signs[i] == 0){
            return i;
        }
    }
    return -1;
}

void branch_and_bound(
    const std::vector<CycleID>& i_cycles,
    const std::vector<IndexedEdge>& i_edges,
    std::vector<int>& fixed_signs,
    double& bestUB,
    int num_vertices,
    int depth
){
    double LB = compute_cycle_packing_LB(i_cycles, fixed_signs, i_edges);
    if(LB >= bestUB){
        return;
    }

    int eid = select_branching_edge(fixed_signs);

    if(eid == -1){
        //have no more nodes to select, we have decided on all of them
        double err = solve_exact_embedding(fixed_signs, i_edges, num_vertices);
        if(err < bestUB){
            bestUB = err;
        }
        return;
    }

    //branch 1
    fixed_signs[eid] = 1;
    branch_and_bound(i_cycles, i_edges, fixed_signs, bestUB, num_vertices, depth-1);
    //branch2
    fixed_signs[eid] = -1;
    branch_and_bound(i_cycles, i_edges, fixed_signs, bestUB, num_vertices, depth+1);

    fixed_signs[eid] = 0;
}


double solve_minerr_dgp1(const std::string& filename){
    
    auto [edges, vertex_ids, vertex_id_to_name] = parse_dgp_instance_dat_file(filename);
    std::unordered_map<int, std::vector<Adjacency>> adj_list = create_adjacency_list_from_edges(edges, vertex_ids.size());
    std::vector<std::vector<Edge>> spanning_tree_edges = get_spanning_forest(adj_list);
    auto cycle_basis = get_cycle_basis(spanning_tree_edges, adj_list);


    std::vector<IndexedEdge> i_edges = build_indexed_edges(edges); //i_edges[id] has .id = id
    std::vector<CycleID> i_cycles =  convert_cycles_to_edge_ids(cycle_basis, i_edges);
    std::vector<double> error_cycle_basis = get_cycle_basis_error(i_cycles, i_edges); //precompute the error of each cycle in the basis, error_cycle_basis[i] has error of indexed_cycles[i], maybe dont need it right now

    
    std::vector<int> fixed_signs = std::vector<int>(static_cast<int>(edges.size()), 0); // all the edges are free
    int depth = 0;
    double bestUB = optimized_projection_minErrDGP1_UB(edges, vertex_ids, adj_list); //obtain a tight UB with our heuristic


    branch_and_bound(i_cycles, i_edges, fixed_signs, bestUB, vertex_ids.size(), 0);
    return bestUB;
}