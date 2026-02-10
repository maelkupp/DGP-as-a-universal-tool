import sys
import subprocess
import os
import time
import gurobipy
from gurobipy import GRB
import tempfile
import re
import numpy as np

### GLOBAL CONSTANTS
myEps = 1e-2
minVal = -10
maxVal = 10
wrongPerturbRange = 3

##############################
# Instance Generation Functions
##############################

def generateSparseIntegerMatrix(m, n, s, min_val=minVal, max_val=maxVal):
    A = np.random.randint(min_val, max_val+1, size=(m, n))
    mask = np.random.random(size=(m, n))
    A[mask > s] = 0
    for i in range(m):
        if abs(np.sum(A[i, :])) < myEps:
            newindex = np.random.randint(0, n)
            A[i, newindex] = np.random.randint(min_val, max_val+1)
    return A

def generate_instance(n, m, s, feasType):
    A = generateSparseIntegerMatrix(m, n, s)
    b = np.zeros(m, dtype=int)
    perturb = None
    perturbColIndex = -1

    if feasType == "feas":
        xrnd = np.random.randint(0, 2, size=n)
        b = A.dot(xrnd)
    elif feasType == "rnd":
        b = np.random.randint(minVal, maxVal+1, size=m)
    elif feasType == "infeas":
        for i in range(m):
            b[i] = int(np.min(A[i, :]) - 1)
    elif feasType == "perturb":
        xrnd = np.random.randint(0, 2, size=n)
        b = A.dot(xrnd)
        perturb = np.random.randint(-wrongPerturbRange, wrongPerturbRange, size=m)
        perturbColIndex = np.random.randint(0, n)
        A[:, perturbColIndex] += perturb
    opb_lines = []
    for i in range(m):
        terms = []
        for j in range(n):
            if A[i, j] > myEps:
                terms.append("+{}*x{}".format(int(A[i,j]), j+1))
            elif A[i, j] < -myEps:
                terms.append("{}*x{}".format(int(A[i,j]), j+1))
        opb_lines.append(" ".join(terms) + " <= {};".format(int(b[i])))
    opb_str = "\n".join(opb_lines)
    return opb_str



def solve_Min_Err_BLP_with_Gurobi(blp_instance_file, time_limit=300, mip_gap=0.01):
    model = gurobipy.read(blp_instance_file)

    model.Params.TimeLimit = time_limit
    model.Params.MIPGap = mip_gap
    model.Params.OutputFlag = 0

    # Remove any existing objective
    model.setObjective(0.0, GRB.MINIMIZE)

    error_vars = []

    # Add one error variable per constraint
    for constr in model.getConstrs():
        e = model.addVar(lb=0.0, name=f"err_{constr.ConstrName}")
        error_vars.append(e)

        # a^T x - e <= b   <=>   a^T x <= b + e
        model.chgCoeff(constr, e, -1.0)

    # Minimize total error
    model.setObjective(gurobipy.quicksum(error_vars), GRB.MINIMIZE)

    model.update()

    start_time = time.perf_counter()
    model.optimize()
    end_time = time.perf_counter()

    if model.status in (GRB.OPTIMAL, GRB.TIME_LIMIT):
        obj_value = model.ObjVal
    else:
        obj_value = float('inf')

    solve_time = end_time - start_time
    return obj_value, solve_time


def algebraic_to_opb(src_path):
    with open(src_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    out_lines = []
    varset = set()
    constr_count = 0

    for line in lines:
        if not line.strip():
            continue

        # replace "a*xk" -> "a xk"
        line = re.sub(r'([+-]?\d+)\s*\*\s*(x\d+)', r'\1 \2', line)

        # collect variables
        for v in re.findall(r'(x\d+)', line):
            varset.add(v)

        out_lines.append(line.strip())
        constr_count += 1

    header = f"* #variable= {len(varset)} #constraint= {constr_count}\n"

    fd, tmp_path = tempfile.mkstemp(suffix=".opb")
    with os.fdopen(fd, "w", encoding="utf-8") as tmp:
        tmp.write(header)
        tmp.write("\n".join(out_lines))

    return tmp_path


"""
if __name__ == "__main__":
    #this python script takes a BLP instance and compares the time it takes for Gurobi to solve the minErr BLP versus
    #the time it takes for us to reduce the instance to a directed DGP instance and then solve the minErrDGP instance with our code
    if len(sys.argv) != 2:
        print("Usage: python blp_solver_comparator.py <blp_instance_file>")
        sys.exit(1)

    blp_file = sys.argv[1]
    new_blp_file = algebraic_to_opb(blp_file)
    # ---- Gurobi BLP Solver ----
    gurobi_obj, gurobi_time = solve_Min_Err_BLP_with_Gurobi(new_blp_file, time_limit=300, mip_gap=0.01)
    print(f"Gurobi MinErr BLP: Objective = {gurobi_obj:.6f}, Time = {gurobi_time:.3f}s\n")

    #first reduce our BLP instance to a directed DGP instance -- we do not take this into account for the timing comparison
    subprocess.run(
        [sys.executable, "blp2dgp_opt.py", blp_file],
        check=True,
        text=True,
        capture_output=True
    )
    dgpf = ''.join(blp_file.split('.')[:-1]) + "-dgp.dat"  # output DGP to AMPL .dat

    print(f"Reduced BLP instance to DGP instance: {dgpf}\n")
    # ---- Our MinErrDGP Solver ----
    t0 = time.perf_counter()
    dgp_ub = subprocess.run(
        ["dir_minErrDGP1_UB.exe", dgpf],
        check=True, capture_output=True, text=True
    )

    print(f"DGP dir solver output:\n{dgp_ub.stdout} in time: {time.perf_counter() - t0:.3f}s\n")
"""



if __name__ == "__main__":
    #this python script takes a BLP instance and compares the time it takes for Gurobi to solve the minErr BLP versus
    #the time it takes for us to reduce the instance to a directed DGP instance and then solve the minErrDGP instance with our code

    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <num_vars> <num_constraints> <density> infeas")
        sys.exit(1)
    
    n_vars = int(sys.argv[1])
    n_cons = int(sys.argv[2])
    density = float(sys.argv[3])

    N = 20
    for i in range(N):
        
        opb_str = generate_instance(n_vars, n_cons, density, "infeas")
        print(f"generate BLP instance: {opb_str}")
        blp_file = f"tmp_instance_{i}.opb"
        with open(blp_file, "w") as f:
            f.write(opb_str)
        new_blp_file = algebraic_to_opb(blp_file)
        # ---- Gurobi BLP Solver ----
        gurobi_obj, gurobi_time = solve_Min_Err_BLP_with_Gurobi(new_blp_file, time_limit=300, mip_gap=0.01)
        print(f"Gurobi MinErr BLP: Objective = {gurobi_obj:.6f}, Time = {gurobi_time:.3f}s\n")

        #first reduce our BLP instance to a directed DGP instance -- we do not take this into account for the timing comparison
        subprocess.run(
            [sys.executable, "blp2dgp_opt.py", blp_file],
            check=True,
            text=True,
            capture_output=True
        )
        dgpf = ''.join(blp_file.split('.')[:-1]) + "-dgp.dat"  # output DGP to AMPL .dat

        print(f"Reduced BLP instance to DGP instance: {dgpf}\n")
        # ---- Our MinErrDGP Solver ----
        t0 = time.perf_counter()
        dgp_ub = subprocess.run(
            ["cheap_minErrDGP1_UB.exe", dgpf],
            check=True, capture_output=True, text=True
        )

        print(f"DGP dir solver output:\n{float(dgp_ub.stdout)/5} in time: {time.perf_counter() - t0:.3f}s\n")




    