#include "dgp.h"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <string>
#include <stack>
#include <queue>
#include <limits>
#include <cmath>
#include <set>



int count_placed_neighbours(
    int v,
    const std::unordered_set<int>& embedding,
    const std::unordered_map<int,std::vector<Adjacency>>& adj_list
){
    int count = 0;
    for(const auto& adj : adj_list.at(v)){
        if(embedding.count(adj.neighbourId)){
            count++;
        }
    }
    return count;
}


std::vector<int> create_single_ordering_for_directed_dgps(
    std::unordered_map<int, std::vector<Adjacency>>& adj_list,
    int root,
    std::unordered_map<int, std::vector<Adjacency>> directed_adj_list,
    int component_size){  
    /*
    this function returns an ordering of the vertices of the connected component so that each placed vertex has one of its neighbours already in the ordering (apart for the first vertex)
    if a vertex has directed neighbours, one of its directed neighbours is already in the ordering,
    for this to work we take the root node to have directed neighbours, so we took the first node of the component which is the A node if taking the BLP -> DGP reduction
    additionally among all valid choices, we will prioritse vertices with the highest number of already placed vertices, that way in the backtracking we will be able to explore as many options early on as possible
    */


    std::vector<int> ordering;
    ordering.reserve(component_size);

    std::unordered_set<int> placed;

    // Candidates that are currently allowed to be placed
    std::vector<int> good_frontier;
    std::vector<int> bad_frontier;

    // Place root
    ordering.push_back(root);
    placed.insert(root);

    // Initialize frontier from root
    for(const auto& adj : adj_list.at(root)){
        int u = adj.neighbourId;
        bool has_directed = !directed_adj_list.at(u).empty();

        if(!has_directed){
            good_frontier.push_back(u);
        }else{
            bad_frontier.push_back(u);
        }
    }

    // Grow ordering
    while(static_cast<int>(ordering.size()) < component_size){

        int v = -1;
        int best_score = -1;

        // --- 1) Try GOOD frontier, pick most constrained ---
        for(int cand : good_frontier){
            if(placed.count(cand)) continue;

            int score = count_placed_neighbours(cand, placed, adj_list);
            if(score > best_score){
                best_score = score;
                v = cand;
            }
        }

        // --- 2) If none, try BAD frontier (reluctantly) ---
        if(v == -1){
            for(int cand : bad_frontier){
                if(placed.count(cand)) continue;

                int score = count_placed_neighbours(cand, placed, adj_list);
                if(score > best_score){
                    best_score = score;
                    v = cand;
                }
            }
        }

        if(v == -1){
            throw std::runtime_error(
                "Failed to construct ordering: no valid frontier vertex"
            );
        }

        // Place v
        ordering.push_back(v);
        placed.insert(v);

        // Update frontiers using neighbours of v
        for(const auto& adj : adj_list.at(v)){
            int u = adj.neighbourId;
            if(placed.count(u)) continue;

            bool has_directed = !directed_adj_list.at(u).empty();

            if(!has_directed){
                good_frontier.push_back(u);
            }else{
                // Check if u is now activated by a placed directed neighbour
                bool activated = false;
                for(const auto& dir_adj : directed_adj_list.at(u)){
                    if(placed.count(dir_adj.neighbourId)){
                        activated = true;
                        break;
                    }
                }

                if(activated){
                    good_frontier.push_back(u);
                }else{
                    bad_frontier.push_back(u);
                }
            }
        }
    }

    return ordering;
}


std::unordered_map<int, std::vector<Adjacency>> get_directed_undirected_adj_list(std::unordered_map<int, std::vector<Adjacency>> adj_list, int num_vertices){
    std::unordered_map<int, std::vector<Adjacency>> directed_adjacency_list;
    for(int i=1; i<= num_vertices; ++i){
        directed_adjacency_list[i] = {};
    }

    for(const auto& [key, adj_list]: adj_list){
        for(const auto& adj: adj_list){
            if(adj.directed){
                directed_adjacency_list[key].push_back(adj);
            }
        }
    }

    return {directed_adjacency_list};

}

std::vector<std::vector<int>> create_vertex_orderings_for_directed_dgps(
    std::unordered_map<int, std::vector<Adjacency>>& adj_list,
    const std::vector<std::set<int>>& components,
    int num_vertices){
        /*
        this function takes the adjacency list and the vector sets of components and returns a vector of orderings of id per connected component
        the ordering is such that each vertex added in the ordering already has one of its neighbour placed earlier in the ordering
        */
        auto directed_adj_list = get_directed_undirected_adj_list(adj_list, num_vertices);
        std::vector<std::vector<int>> orderings;
        std::cout << "Create the vertex ordering for each of the " << components.size() << " connected components \n";
        for(auto& component: components){
            orderings.push_back(create_single_ordering_for_directed_dgps(adj_list, *std::next(component.begin(), 0),directed_adj_list, component.size()));
        }

        return orderings;
    }


void dir_backtrack_place(
    const std::vector<int>& ordering,
    int index,
    std::unordered_map<int, double>& current_embedding,
    double partial_error,
    BacktrackState& state,
    std::unordered_map<int, std::vector<Adjacency>>& adj_list,
    std::unordered_map<int,int>& pos,
    std::unordered_map<int, std::string> vertex_id_to_name){
    /*
    this function is a recursive backtracking algorithm that will find the embedding that minimises the error
    we carry with us a partial error store in the BacktrackState class so that if the partial error surpasses the
    best error we can backtrack as this cannot be the best embedding
    */

    if(index == 0){
        // place the first vertex at the origin
        current_embedding[ordering[index]] = 0.0;
        dir_backtrack_place(ordering, 1, current_embedding, partial_error, state, adj_list, pos, vertex_id_to_name);
        return;
    }

    if(index == 1){
        //when we are placing the second vertex, we do not need to branch as that will simply give us the mirror image of the of the solution which we do not care about
        //this saves a lot of time as we reduce the branching done early on
        for(const auto& adj: adj_list[ordering[1]]){
            if(adj.neighbourId == ordering[0]){
                if(adj.directed){
                    if(adj.outgoing){
                        current_embedding[ordering[1]] = current_embedding[ordering[0]] - adj.dist;
                    }else{
                        current_embedding[ordering[1]] = current_embedding[ordering[0]] + adj.dist;
                    }
                }else{
                    //otherwise per default we embed it to the right
                    current_embedding[ordering[1]] = current_embedding[ordering[0]] + adj.dist;
                }
                break;
            }
        }
        dir_backtrack_place(ordering, 2, current_embedding, partial_error, state, adj_list, pos, vertex_id_to_name);
        return;
    }

    if(partial_error > state.best_error){
        //this path cannot be optimal so we return
        return;
    }

    if(index == static_cast<int>(ordering.size())){
        //we have placed all the vertices so we compare the partial error with the best error so far
        if(partial_error < state.best_error){
            std::cout << "improving the minimum error from " << state.best_error << " to " << partial_error << "\n";
            std::cout << "This new best embedding is: \n";
            display_1Dembedding(current_embedding, vertex_id_to_name);
            state.best_error = partial_error;
            state.best_embedding = current_embedding;
        }
        return;
    }

    //need to change this so that we look for the first directed edge with which it is neighbour, if it does not have any then go on an look for its undirected neighbour
    int vertex_id_to_place = ordering[index];
    double delta_error {0.0};
    bool found_directed = false;
    double implied_pos;
    for(const auto& adj: adj_list[vertex_id_to_place]){
        //iterate through all of its neighbours to check where we can place the vertex, need to find the neighbour that has already been embedded
        if(current_embedding.find(adj.neighbourId) != current_embedding.end()){
            //this neighbour has already been placed
            if(adj.directed){
                //this is a directed edge so we only need to embed it it one direction, use the convention u->v means v is to the right of u
                //as it is a directed edge we can trust it and place it as is
                if(adj.outgoing){
                    current_embedding[vertex_id_to_place] = current_embedding[adj.neighbourId] - adj.dist;
                }else{
                    current_embedding[vertex_id_to_place] = current_embedding[adj.neighbourId] + adj.dist;
                }

                if(!found_directed){
                    implied_pos = current_embedding[vertex_id_to_place];
                    found_directed = true;
                }else if(std::abs(implied_pos - current_embedding[vertex_id_to_place]) > 1e-9){
                        return; //prune this branch as it cannote be correct
                }

            }
        }
    }


    if(found_directed){
        //we have found a directed edge and checked that it was valid for all other directed egdes
        current_embedding[vertex_id_to_place] = implied_pos;
        delta_error = compute_vertex_embedding_error(vertex_id_to_place, current_embedding, adj_list, pos);
        dir_backtrack_place(ordering, index+1, current_embedding, partial_error + delta_error, state, adj_list, pos, vertex_id_to_name);
        current_embedding.erase(vertex_id_to_place);
        return; //we return after this as we know that the vertex needs to be positionned
    }


    //if we end up here it means that the vertex has no directed neighbour so we try out all the different branching options for this edge by placing it in all possible positions
    for(const auto& adj: adj_list[vertex_id_to_place]){
        if(current_embedding.find(adj.neighbourId) != current_embedding.end()){
            current_embedding[vertex_id_to_place] = current_embedding[adj.neighbourId] + adj.dist;
            delta_error = compute_vertex_embedding_error(vertex_id_to_place, current_embedding, adj_list, pos);
            if(partial_error + delta_error > state.best_error){
                continue; //to stop unnecessary recursion
            }
            dir_backtrack_place(ordering, index+1, current_embedding, partial_error + delta_error, state, adj_list, pos, vertex_id_to_name);
            current_embedding.erase(vertex_id_to_place);

            current_embedding[vertex_id_to_place] = current_embedding[adj.neighbourId] - adj.dist;        
            delta_error = compute_vertex_embedding_error(vertex_id_to_place, current_embedding, adj_list, pos);
            if(partial_error + delta_error > state.best_error){
                continue; //to stop unnecessary recursion
            }
            dir_backtrack_place(ordering, index+1, current_embedding, partial_error + delta_error, state, adj_list, pos, vertex_id_to_name);
            current_embedding.erase(vertex_id_to_place);
        }
    }
}


std::unordered_map<int, double> compute_min_Err_dir_DGP1_UB(
    std::vector<std::vector<int>>& orderings,
    std::unordered_map<int, std::vector<Adjacency>>& adj_list, 
    double& min_error,
    std::unordered_map<int, std::string>& vertex_id_to_name){
    
    std::unordered_map<int, double> best_total_embedding;
    
    for(auto& ordering: orderings){
        //compute the best embedding for each distinct connected component, add the cumulative error and record the embeddings in one big hash map
        std::unordered_map<int, double> current_embedding;
        BacktrackState state;
        std::unordered_map<int,int> pos;
        for (int i = 0; i < (int)ordering.size(); ++i) {
            pos[ordering[i]] = i;
        }
        dir_backtrack_place(ordering, 0, current_embedding, 0, state, adj_list, pos, vertex_id_to_name);

        //accumulate the total error
        min_error += state.best_error;

        for(const auto& [vertex_id, coord]: state.best_embedding){
            best_total_embedding[vertex_id] = coord;
        }
    }
    return best_total_embedding;
}


int main(int argc, char* argv[]){
    if(argc != 2){
        std::cout << "Usage: " << argv[0] << "dgp_instance_file.dat \n";
        std::cout << "The point of this script is to calculate an efficient upper bound on directed DGP instances, this is useful when thinking about the reduction BLP -> DGP which can use directed edges \n";
        return 1;
    }
    auto [edges, vertex_ids, vertex_id_to_name] = parse_dgp_instance_dat_file(argv[1]);
    auto adj_list = create_adjacency_list_from_edges(edges, vertex_ids.size());


    auto components = get_connected_components(adj_list, vertex_ids);
    std::vector<std::vector<int>> orderings = create_vertex_orderings_for_directed_dgps(adj_list, components, vertex_ids.size());

    std::vector<int> ordering = orderings[0];
    for(auto& id: ordering){
        std::cout << id << ": " << vertex_id_to_name[id] << " ";
    }
    double min_error {0.0};
    std::unordered_map<int, double> best_embedding = compute_min_Err_dir_DGP1_UB(orderings, adj_list, min_error, vertex_id_to_name);

    std::cerr << "Have found an optimal embedding with error: " << min_error << "\n";
    std::cout << min_error;
    display_1Dembedding(best_embedding, vertex_id_to_name);

    verbose_compute_1Dembedding_error(edges, best_embedding, vertex_id_to_name);
}