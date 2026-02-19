import sys
import subprocess
import random
from upper_bounds.minErrDGP1 import min_err_dgp1_gurobi, min_err_dgp1_slab_gurobi, parse_dgp_dat_file

def generate_feasible_dat_graph(filename, n=8, m=10, A=10):
    #randomly position n points on the line at integer coordinates between 0 and A
    points = random.sample(range(A), n)
    points.sort()
    edges = set()
    while len(edges) < m:
        i = random.randint(1, n)
        j = random.randint(1, n)
        if i != j:
            if i < j:
                edges.add(tuple((i, j, abs(points[i-1] - points[j-1]))))
            else:
                edges.add(tuple((j, i, abs(points[i-1] - points[j-1]))))

    with open(filename, "w") as f:
        # header
        f.write("# DGP instance in .dat format for dgp.mod\n")
        f.write("# this graph is a DGP encoding of a random instance\n\n")

        # edge list
        f.write("param : E : c I :=\n")

        for (i, j, w) in edges:
            #the weight is not random here as it comes from the distance between the points in the 1D embedding
            f.write(f"  {i} {j}  {w:.3f} 0\n")

        f.write(";\n")

def generate_dat_graph(filename, n=8, m=12, weight_min=1, weight_max=5):
    edges = set()
    while len(edges) < m:
        i = random.randint(1, n)
        j = random.randint(1, n)
        if i != j:
            edges.add(tuple(sorted((i, j))))

    edges = list(edges)

    with open(filename, "w") as f:
        # header
        f.write("# DGP instance in .dat format for dgp.mod\n")
        f.write("# this graph is a DGP encoding of a random instance\n\n")

        # edge list
        f.write("param : E : c I :=\n")

        for (i, j) in edges:
            #here attribute a random weight to the edge between weight_min and weight_max
            w = random.randint(weight_min, weight_max)
            f.write(f"  {i} {j}  {w:.3f} 0\n")

        f.write(";\n")


def create_adjacency_list_from_edges(edges, n):
    adj_list = {i: [] for i in range(0, n)}
    for (i, j, w) in edges:
        adj_list[i].append((j, w)) 
        adj_list[j].append((i, w)) 

    return adj_list

def dfs_ordering(adj_list, root):
    visited = set()
    order = []

    stk = [root]
    while stk:
        node = stk.pop()
        if node not in visited:
            visited.add(node)
            order.append(node)
            for (nbr, w) in adj_list[node]:
                if nbr not in visited:
                    stk.append(nbr)
    return order

def feasibility_from_ordering(edges, minim_order):
    adj_list = create_adjacency_list_from_edges(edges, len(minim_order))
    #this is the dictionary that will store the position of each vertex on the line
    pos_map = {minim_order[0]: 0.0} #we can place the first vertex in the minimizing order at position 0 on the line
    index_map = {v: i for i, v in enumerate(minim_order)} #map from vertex id to its index in the minimizing order
    dfs_order = dfs_ordering(adj_list, minim_order[0]) #get a DFS ordering of the vertices starting from the first vertex in the minimizing order
    for i in range(1, len(dfs_order)):
        v = dfs_order[i]
        ngbrs = [nbr for nbr in adj_list[v] if nbr[0] in pos_map] #get the neighbors of v that have already been placed on the line
        if not ngbrs:
            #if v has no neighbors that have been placed on the line yet, we can place it anywhere, we choose to place it at position 0
            pos_map[v] = 0.0
        else:
            if index_map[ngbrs[0][0]] < index_map[v]:
                cand_pos = pos_map[ngbrs[0][0]] + ngbrs[0][1]
            else:
                cand_pos = pos_map[ngbrs[0][0]] - ngbrs[0][1] 
            for j in range(1, len(ngbrs)):
                (nbr, w) = ngbrs[j]
                new_cand_pos = None
                if index_map[nbr] < index_map[v]:
                    new_cand_pos = pos_map[nbr] + w
                else:
                    new_cand_pos = pos_map[nbr] - w
                if abs(new_cand_pos - cand_pos) > 1e-6:
                    #if the new candidate position is not consistent with the previous candidate position, then the ordering is not feasible
                    return False
            pos_map[v] = cand_pos
    return True



    #now that everything has been placed according to the minimizing order we just need to verify that all the distance constraints are satisfied
    return True

if __name__ == "__main__":
    #create a random DGP instance and obtain its edges and vertices
    if len(sys.argv) != 2:
        print("Usage: python -m orderings.obtain_minimising_ordering infeas/feas")
        sys.exit(1)
    
    if sys.argv[1] != "feas" and sys.argv[1] != "infeas":
        print("Usage: python -m orderings.obtain_minimising_ordering infeas/feas")
        sys.exit(1)
    
    if sys.argv[1] == "feas":
        generate_feasible_dat_graph("temp_graph.dat", n=5, m=8, A=10)
    else:
        #not necessarily an infeasible instance but highly likely to be one
        generate_dat_graph("temp_graph.dat", n=5, m=8, weight_min=1, weight_max=10)  

    edges, n, vertex_names = parse_dgp_dat_file("temp_graph.dat")

    #solve the DGP instance with Gurobi, to obtain an upper bound and the embedding of the vertices on the line
    ub, emb = min_err_dgp1_gurobi(edges, n, time_limit=60, mip_gap=0.05)
    print(f"Upper bound objective = {ub:.6f}", file=sys.stderr)
    print(f"emd: {emb} with error {ub}", file=sys.stderr) #sort the embedding by the x coordinate of the vertices

    minimizing_order = sorted(range(n), key=lambda i: emb[i])
    print(f"minimizing order: {minimizing_order}", file=sys.stderr)


    epsilons = [5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001]

    for epsilon in epsilons:
        print(f"Solving minErrDGP1 with slab of height {epsilon}", file=sys.stderr)
        ub, emb = min_err_dgp1_slab_gurobi(edges, n, epsilon=0.1, time_limit=60, mip_gap=0.05)
        print(f"Upper bound objective with slab of height {epsilon} = {ub:.6f}", file=sys.stderr)
        print(f"emd with slab of height {epsilon}: {emb} with error {ub}", file=sys.stderr)
        minim_order = sorted(range(n), key=lambda i: emb[i][0]) #sort wrt to the x coordinate, we do a projection
        print(f"minimizing order with slab of height {epsilon}: {minim_order}", file=sys.stderr)
        print(f"Same order as 1D embedding minimizing order? {minim_order == minimizing_order}", file=sys.stderr)
        if minim_order == minimizing_order:
            #exit the loop
            break


    #is_feas = feasibility_from_ordering(edges, minimizing_order)
    #print(f"Is the minimizing order feasible? {is_feas}", file=sys.stderr)





