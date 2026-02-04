import subprocess
import sys

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
        n = get_number_of_vertices()  # Assume this function is defined elsewhere
        p = get_edge_density()
        A = get_square_length()
    else:
        n = 10
        p = 0.3
        A = 10.0

    N = 50
    print(f"Going to do {N} tests \n")
    count = 0
    for i in range(N):
        if i%10 == 0:
            print(f"test number {i} have {count} / {i} UB >= Relaxation \n")

        
        relax_value = subprocess.run(["relax_gen2Dembedding.exe", str(n), str(p), str(A)],
                    check=True, text=True, capture_output=True)
        values = list(map(float, relax_value.stdout.strip().split(",")))
        print(values)
        rotate_relax_value = values[0]
        total_relax_value = values[1]
        print(f"rotational relaxation: {values[0]} \ntotal relaxation: {values[1]}\n")
        
        #using sys.executable to ensure the same Python interpreter is used
        UB_value = subprocess.run([sys.executable, "minErrDGP1.py", "temp_graph.dat"],
                    check=True, text=True, capture_output=True)
        UB_lines = UB_value.stdout.strip().splitlines()
        UB = float(UB_lines[-1])
        print(f"UB value: {UB}")
        if UB >= rotate_relax_value:
            count += 1
            print("UB >= Relaxation \n")
        else:
            print("UB < Relaxation \n")
            print(f"UB: {UB}, rot Relaxation: {rotate_relax_value} total: {total_relax_value} \n")
    print(f"Out of {N} tests, UB >= Relaxation in {count} cases.")




"""

formalize reduction proof (Mael)
think about how to turn BP into BPB (both of us)
find adversarial instances where Gurobi is slow and perhaps we are not (Mael)
fill in the proofs in the current LaTeX "notes" document (Mael)
start writing a journal paper (Leo)
"""