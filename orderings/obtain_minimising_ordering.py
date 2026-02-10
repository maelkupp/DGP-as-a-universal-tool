import sys
import subprocess
import random
from upper_bounds.minErrDGP1 import min_err_dgp1_gurobi, parse_dgp_dat_file


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
            w = random.randint(weight_min, weight_max)
            f.write(f"  {i} {j}  {w:.3f} 0\n")

        f.write(";\n")


#need to create functions that generate feasible/infeasible instances of the DGP directly into the .dat format
# feasible instance will be generate by positionning n points on the line at integer coordinates and then creating edges 
#between them with weights equal to the distance between the points



if __name__ == "__main__":
    #create a random DGP instance and obtain its edges and vertices
    generate_dat_graph("temp_graph.dat", n=8, m=12, weight_min=1, weight_max=5)
    edges, n, vertex_names = parse_dgp_dat_file("temp_graph.dat")

    #solve the DGP instance with Gurobi, to obtain an upper bound and the embedding of the vertices on the line
    ub, emb = min_err_dgp1_gurobi(edges, n, time_limit=60, mip_gap=0.05)
    


