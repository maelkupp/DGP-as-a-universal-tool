import subprocess
import sys
import time


def get_number_of_vertices():
    n = -1
    while n < 1:
        try:
            n = int(input("Enter the number of vertices (n > 0): "))
        except ValueError:
            print("Invalid input. Please enter a positive integer.")
    return n

def get_edge_density():
    p = -1.0
    while p < 0.0 or p > 1.0:
        try:
            p = float(input("Enter the edge density (0.0 <= p <= 1.0): "))
        except ValueError:
            print("Invalid input. Please enter a number between 0.0 and 1.0.")
    return p

def get_square_length():
    A = -1.0
    while A <= 0.0:
        try:
            A = float(input("Enter the square length (A > 0.0): "))
        except ValueError:
            print("Invalid input. Please enter a positive number.")
    return A

if __name__ == "__main__":
    # Run the relaxation orchestrator script using subprocess

    #quick if statement for now to speed up processing during testing
    if len(sys.argv) > 1:
        n = get_number_of_vertices()
        p = get_edge_density()
        A = get_square_length()
    else:
        n = 15
        p = 0.3
        A = 10.0

    N = 20
    print(f"Going to do {N} trials \n")
    for i in range(N):
        # ---- Cheap UB ----
        t0 = time.perf_counter()
        cheap_UB = subprocess.run(
            ["cheap_minErrDGP1_UB.exe", str(n), str(p), str(A)],
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

        print(
            f"Trial {i}: "
            f"Cheap UB = {cheap_UB_value:.6f} "
            f"(time = {cheap_time:.3f}s), "
            f"Gurobi UB = {Gurobi_UB_value:.6f} "
            f"(time = {gurobi_time:.3f}s)\n"
            f"(ratio = {cheap_UB_value / Gurobi_UB_value:.4f})\n"
        )

