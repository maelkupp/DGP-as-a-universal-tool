#include "dgp.h"
#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_set>

std::unordered_map<int, int> create_num_neighbours_map(std::unordered_map<int, std::vector<Adjacency>>&adj_list){
    std::unordered_map<int, int> num_neighbours_map;
    for(const auto& [vertex_id, adj_list]: adj_list){
        num_neighbours_map[vertex_id] = adj_list.size();
    }

    return num_neighbours_map;
}

std::vector<int> create_single_vertex_ordering(const std::unordered_map<int, std::vector<Adjacency>>& adj_list, const std::unordered_map<int, int>& degree, const std::set<int>& component)
{
    std::vector<int> ordering;
    std::unordered_set<int> placed;
    std::unordered_set<int> in_frontier;

    // find highest-degree vertex in component
    int root = *std::max_element(
        component.begin(), component.end(),
        [&](int a, int b){ return degree.at(a) < degree.at(b); }
    );

    ordering.push_back(root);
    placed.insert(root);

    // priority queue: max degree first
    std::priority_queue<std::pair<int,int>> frontier;

    for (const auto& adj : adj_list.at(root)) {
        if (component.count(adj.neighbourId)) {
            frontier.push({degree.at(adj.neighbourId), adj.neighbourId});
            in_frontier.insert(adj.neighbourId);
        }
    }

    while (!frontier.empty()) {
        auto [deg, v] = frontier.top();
        frontier.pop();

        if (placed.count(v)) continue;

        ordering.push_back(v);
        placed.insert(v);

        for (const auto& adj : adj_list.at(v)) {
            int u = adj.neighbourId;
            if (component.count(u) && !placed.count(u) && !in_frontier.count(u)) {
                frontier.push({degree.at(u), u});
                in_frontier.insert(u);
            }
        }
    }

    return ordering;
}

std::vector<std::vector<int>> create_vertex_orderings(std::unordered_map<int, std::vector<Adjacency>>& adj_list, std::vector<std::set<int>>& components, int num_vertices){
    std::unordered_map<int, int> num_neighbours_map = create_num_neighbours_map(adj_list);
    std::vector<std::vector<int>> orderings;
    for(const auto& component: components){
        orderings.push_back(create_single_vertex_ordering(adj_list, num_neighbours_map, component));
    }
    return orderings;
}




void backtrack_place(
    const std::vector<int>& ordering,
    int index,
    std::unordered_map<int, double>& current_embedding,
    double partial_error,
    BacktrackState& state,
    std::unordered_map<int, std::vector<Adjacency>>& adj_list,
    std::unordered_map<int,int>& pos,
    std::unordered_map<int, std::string> vertex_id_to_name){
    
    int v = ordering[index];
    if(index == 0){
        current_embedding[v] = 0.0;
        backtrack_place(ordering, 1, current_embedding, 0.0, state, adj_list, pos, vertex_id_to_name);
        current_embedding.erase(v);
        return;
    }

    if(index == 1){
        for(const auto& adj: adj_list[v]){
            //going through all the neighbours of the 2nd vertex in the ordering, we want to find the first vertex we placed
            if(adj.neighbourId == ordering[0]){
                //have found the 1st vertex of the ordering
                current_embedding[v] = adj.dist; // 0.0 + adj.dist = adj.dist we place the 2nd vertex to the right of the 1st vertex makes the algo go faster no branching on first edge
            }
            backtrack_place(ordering, 2, current_embedding, 0.0, state, adj_list, pos, vertex_id_to_name);
            current_embedding.erase(v);
            return;
        }
    }

    if(partial_error > state.best_error){
        return; // we prune this branch as it cannot be the best one
    }

    if(index == static_cast<int>(ordering.size())){
        //we have placed all the vertices so we compare the partial error with the best error so far
        std::cout << "have reached final level \n";
        if(partial_error < state.best_error){
            std::cout << "improving the minimum error from " << state.best_error << " to " << partial_error << "\n";
            std::cout << "This new best embedding is: \n";
            display_1Dembedding(current_embedding, vertex_id_to_name);
            state.best_error = partial_error;
            state.best_embedding = current_embedding;
        }
        return;
    }

    //----- the general case now that it is not the first two vertices that we are placing
    double delta_error {0.0};
    for(const auto& adj: adj_list[v]){
        //iterating through all the neighbours of v
        if(current_embedding.find(adj.neighbourId) != current_embedding.end()){
            //we have already placed the neighbour so we can branch on this edge
            current_embedding[v] = current_embedding[adj.neighbourId] + adj.dist;
            delta_error = compute_vertex_embedding_error(v, current_embedding, adj_list, pos);
            if(partial_error + delta_error > state.best_error){
                continue; //to stop unnecessary recursion
            }
            backtrack_place(ordering, index+1, current_embedding, partial_error + delta_error, state, adj_list, pos, vertex_id_to_name);
            current_embedding.erase(v);

            current_embedding[v] = current_embedding[adj.neighbourId] - adj.dist;        
            delta_error = compute_vertex_embedding_error(v, current_embedding, adj_list, pos);
            if(partial_error + delta_error > state.best_error){
                continue; //to stop unnecessary recursion
            }
            backtrack_place(ordering, index+1, current_embedding, partial_error + delta_error, state, adj_list, pos, vertex_id_to_name);
            current_embedding.erase(v);
        }
    }
}

std::unordered_map<int, double> compute_min_Err_DGP1_UB(std::vector<std::vector<int>>& orderings,
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
        backtrack_place(ordering, 0, current_embedding, 0, state, adj_list, pos, vertex_id_to_name);

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
        std::cout << "Usage: " << argv[0] << " dgp_filename.dat \n";
        std::cout << "This script is used to get as tight an upper bound as possible using discrete techniques on MinErrDGP1 of a DGP instance \n";
        return 1;
    }

    auto [edges, vertex_ids, vertex_id_to_name] = parse_dgp_instance_dat_file(argv[1]);
    auto adj_list = create_adjacency_list_from_edges(edges, vertex_ids.size());

    auto components = get_connected_components(adj_list, vertex_ids);
    std::vector<std::vector<int>> orderings = create_vertex_orderings(adj_list, components, vertex_ids.size());

    std::vector<int> ordering = orderings[0];
    for(auto& id: ordering){
        std::cout << id << ": " << vertex_id_to_name[id] << " ";
    }
    std::cout << "\n";

    double min_error {0.0};
    std::unordered_map<int, double> best_embedding = compute_min_Err_DGP1_UB(orderings, adj_list, min_error, vertex_id_to_name);

    std::cout << "Have found an optimal embedding with error: " << min_error << "\n";
    display_1Dembedding(best_embedding, vertex_id_to_name);

    verbose_compute_1Dembedding_error(edges, best_embedding, vertex_id_to_name);

}