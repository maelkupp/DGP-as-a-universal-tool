import random
import datetime
import sys
import subprocess
import time
import numpy as np
import pandas as pd


def get_number_of_vertices():
    n = -1
    while n < 1:
        try:
            n = int(input("Enter the number of vertices (n > 0): "))
        except ValueError:
            print("Invalid input. Please enter a positive integer.")
    return n

def get_number_edges():
    m = -1
    while m < 0:
        try:
            m = int(input("Enter the number of edges (m >= 0): "))
        except ValueError:
            print("Invalid input. Please enter a non-negative integer.")
    return m

def get_minimum_weights():
    w_min = -1
    while w_min < 0:
        try:
            w_min = float(input("Enter the minimum edge weight (w_min >= 0): "))
        except ValueError:
            print("Invalid input. Please enter a non-negative number.")
    return w_min

def get_maximum_weights(min_weight):
    w_max = -1
    while w_max < min_weight:
        try:
            w_max = float(input(f"Enter the maximum edge weight (w_max >= {min_weight}): "))
        except ValueError:
            print("Invalid input. Please enter a number greater than or equal to the minimum weight.")
    return w_max


def generate_dat_graph(filename, n=10, m=15, weight_min=1, weight_max=5):
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

    if len(sys.argv) > 1:
        n = get_number_of_vertices() 
        m = get_number_edges()
        weight_min = get_minimum_weights()
        weight_max = get_maximum_weights()
    else:
        n = 50
        m = 75
        weight_min = 10
        weight_max = 15
    #generate random graph and save to temp_graph.dat

    n_values = [20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]
    density_levels = [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6]
    REPS = 20

    results = []
    for n in n_values:
        for density in density_levels:
            m = int(density * n*(n-1)/2)

            for weight_range in [(1,5), (5,10), (10,50)]:

                for rep in range(REPS):

                    generate_dat_graph("temp_graph.dat",
                                    n, m,
                                    weight_range[0],
                                    weight_range[1])

                    # --- Cycle relaxation ---
                    proc1 = subprocess.run(
                        ["relaxation\\minErrDGP1_cycles_LB.exe",
                        "temp_graph.dat", "0"],
                        check=True, text=True, capture_output=True)
                    cycle_val, cycle_time = map(float, proc1.stdout.split())

                    # --- Greedy relaxation ---
                    proc2 = subprocess.run(
                        ["relaxation\\greedy_packing_cycle_LB.exe",
                        "temp_graph.dat"],
                        check=True, text=True, capture_output=True)
                    greedy_val, greedy_time = map(float, proc2.stdout.split())

                    results.append({
                        "n": n,
                        "m": m,
                        "density": density,
                        "weight_min": weight_range[0],
                        "weight_max": weight_range[1],
                        "cycle_val": cycle_val,
                        "greedy_val": greedy_val,
                        "cycle_time": cycle_time,
                        "greedy_time": greedy_time,
                        "greedy_better": greedy_val > cycle_val
                    })

    df = pd.DataFrame(results)
    df.to_csv("experiment_results.csv", index=False)




