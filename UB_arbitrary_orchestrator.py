import sys
import subprocess
import time
import re
import random


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


if __name__ == "__main__":
    graph_filename = "temp_graph.dat"

    N = 25
    UB_ratios = []
    time_ratios = []
    print(f"Going to do {N} trials \n")
    for i in range(N):
        #generate random graph and save to temp_graph.dat
        generate_dat_graph(graph_filename)

        # ---- Cheap UB ----
        t0 = time.perf_counter()
        cheap_UB = subprocess.run(
            ["cheap_minErrDGP1_UB.exe", graph_filename],
            check=True,
            text=True,
            capture_output=True
        )
        t1 = time.perf_counter()
        cheap_time = t1 - t0

        # ---- Gurobi UB ----
        t0 = time.perf_counter()
        Gurobi_UB = subprocess.run(
            [sys.executable, "minErrDGP1.py", "temp_graph.dat"],
            check=True,
            text=True,
            capture_output=True
        )
        t1 = time.perf_counter()
        gurobi_time = t1 - t0

        # ---- Parse outputs ----
        UB_lines = Gurobi_UB.stdout.strip().splitlines()
        cheap_UB_value = float(cheap_UB.stdout.strip())
        Gurobi_UB_value = float(UB_lines[-1].split()[-1])

        UB_ratios.append(cheap_UB_value / Gurobi_UB_value)
        time_ratios.append(cheap_time / gurobi_time)

        print(
            f"Trial {i}: "
            f"Cheap UB = {cheap_UB_value:.6f} "
            f"(time = {cheap_time:.3f}s), "
            f"Gurobi UB = {Gurobi_UB_value:.6f} "
            f"(time = {gurobi_time:.3f}s)\n"
        )
        if Gurobi_UB_value == 0.0:
            print(f"Gurobi UB is zero {Gurobi_UB_value} for trial {i}, cheap UB value {cheap_UB_value}.\n")
        else:
            print(f"(UB ratio = {cheap_UB_value / Gurobi_UB_value:.4f})\n")
        print(f"(time ratio = {cheap_time / gurobi_time:.4f})\n")

    avg_UB_ratio = sum(UB_ratios) / N
    avg_time_ratio = sum(time_ratios) / N
    print(f"Average UB ratio (Cheap / Gurobi) over {N} trials: {avg_UB_ratio:.4f}")
