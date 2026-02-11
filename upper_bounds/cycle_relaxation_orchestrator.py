import random
import datetime
import sys
import subprocess

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
        n = 5
        m = 8
        weight_min = 1
        weight_max = 5
    #generate random graph and save to temp_graph.dat

    N = 50
    count = 0
    for i in range(N):
        
        if i%10 == 0:
            print(f"test number {i} have {count} / {i} UB >= Relaxation \n")
        #generate a random graph
        generate_dat_graph("temp_graph.dat", n, m, weight_min, weight_max)
        #implement the cycle based relaxation
        cycle_relax_proc = subprocess.run(["minErrDGP1_cycles_LB.exe", "temp_graph.dat", "0"],
                                           check=True, text=True, capture_output=True)

        #getting the gurobi UB value

        cycle_relax_value = float(cycle_relax_proc.stdout.strip())

        print(f"relax value {cycle_relax_value} \n")

        #using sys.executable to ensure the same Python interpreter is used
        UB_proc = subprocess.run([sys.executable, "minErrDGP1.py", "temp_graph.dat"],
                    check=True, text=True, capture_output=True)
        
        UB_lines = UB_proc.stdout.strip().splitlines()
        UB = float(UB_lines[-1])
        
        print(f" UB {UB} \n")
        if UB >= cycle_relax_value:
            count += 1
            print("UB >= Relaxation \n")
            print(f"UB value: {UB} \n Cycle Relaxation value: {cycle_relax_value} \n")
        else:
            print("UB < Relaxation \n")
            break
    print(f"Out of {N} tests, UB >= Relaxation in {count} cases.")




