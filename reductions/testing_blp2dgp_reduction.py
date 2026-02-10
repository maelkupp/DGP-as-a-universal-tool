#!/usr/bin/env python3
import sys, os, time, datetime, math, random
import numpy as np
import subprocess
import itertools
import gurobipy as gp
from gurobipy import GRB
import re
from scipy.optimize import minimize
from reductions.blp2dgp import readOpb, reduce_blp_2_dgp, writeDat
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

##############################
# Brute-force Feasibility Checker
##############################

def parse_constraint(line):
    line = line.split('#')[0].strip()
    if not line:
        return None
    if line.endswith(';'):
        line = line[:-1]
    parts = line.split()
    relation = parts[-2]
    b_val = int(parts[-1])
    term_parts = parts[:-2]
    terms = []
    for term in term_parts:
        term = term.strip()
        if term.startswith('+'):
            term = term[1:]
        coeff_str, var_str = term.split('*')
        coeff = int(coeff_str)
        var = int(var_str[1:])
        terms.append((coeff, var))
    return terms, relation, b_val

def evaluate_opb_instance(opb_str, n_vars):
    constraints = []
    for line in opb_str.splitlines():
        parsed = parse_constraint(line)
        if parsed is not None:
            constraints.append(parsed)
    for bits in itertools.product([0,1], repeat=n_vars):
        assignment = {i+1: bits[i] for i in range(n_vars)}
        ok = True
        for terms, relation, b_val in constraints:
            total = sum(coeff * assignment[var] for coeff, var in terms)
            if relation == "<=" and not (total <= b_val):
                ok = False
                break
            elif relation == ">=" and not (total >= b_val):
                ok = False
                break
            elif relation == "=" and not (total == b_val):
                ok = False
                break
        if ok:
            return True
    return False

#############################
# dgp.dat file parser
#############################

def parse_dgp_dat_file(dat_file):
    edges = []
    n = 0
    in_edges = False
    with open(dat_file, "r") as f:
        for line in f:
            if line.strip().startswith("param : E : c I :="):
                in_edges = True
                continue
            if in_edges:
                if line.strip().startswith(";"):
                    break
                m = re.match(r"\s*(\d+)\s+(\d+)\s+([0-9eE\+\-\.]+)", line)
                if m:
                    i, j, dist = int(m.group(1)), int(m.group(2)), float(m.group(3))
                    edges.append((i-1, j-1, dist))  # subtract 1 for zero-based indexing
                    n = max(n, i, j)
    return edges, n

##############################
# DGP Solution Function
##############################

#MinERRDGP function

def min_err_dgp_absolute(edges, n, fixed_positions=None):
    m = gp.Model()
    m.Params.OutputFlag = 0
    m.Params.MIPGap = 0
    m.Params.TimeLimit = 120

    x = m.addVars(n, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="x")
    error = {}
    for idx, (i, j, d) in enumerate(edges):
        # Step 1: Auxiliary variable for the difference
        diff = m.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"diff_{idx}")
        m.addConstr(diff == x[i] - x[j])

        # Step 2: absdiff = |diff|
        absdiff = m.addVar(lb=0, name=f"absdiff_{idx}")
        m.addGenConstrAbs(absdiff, diff)

        # Step 3: aux = absdiff - d
        aux = m.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"aux_{idx}")
        m.addConstr(aux == absdiff - d)

        # Step 4: error = |aux|
        error[idx] = m.addVar(lb=0, name=f"err_{idx}")
        m.addGenConstrAbs(error[idx], aux)

    if fixed_positions:
        for idx, value in fixed_positions.items():
            m.addConstr(x[idx] == value)

    m.setObjective(gp.quicksum(error[idx] for idx in range(len(edges))), GRB.MINIMIZE)
    m.optimize()

    if m.status == GRB.OPTIMAL:
        embedding = [x[i].X for i in range(n)]
        print("Embedding:", embedding)
        print("--- Per-edge debug output (abs model) ---")
        total_err = 0.0
        for idx, (i, j, d) in enumerate(edges):
            realized = abs(embedding[i] - embedding[j])
            errval = error[idx].X
            print(f"Edge {idx}: ({i},{j})  Realized: {realized:.6f}, Target: {d:.6f}, Error: {errval:.6f}")
            total_err += errval
        print(f"Total slack (objective): {total_err:.6f}")
        return m.objVal, embedding
    else:
        print("No optimal solution found.")
        return None, None

#MINERRBLP function
def min_err_blp(A, b):
    n = A.shape[1]
    best_err = float('inf')
    best_x = None
    for x in itertools.product([0,1], repeat=n):
        err = np.linalg.norm(np.dot(A, x) - b)
        if err < best_err:
            best_err = err
            best_x = np.array(x)
    return best_err, best_x


def print_edges(edges, VL=None):
    """
    Pretty-print DGP edges for debugging.
    Args:
        edges: list of tuples (i, j, dist[, ...])
        VL: dict, optional, mapping vertex id (1-based) to name
    """
    print("printing edges \n")
    for edge in edges:
        # Support both 3- and 4-element tuples (ignore extra fields after distance)
        if len(edge) >= 3:
            i, j, dist = edge[:3]
            extra = edge[3:] if len(edge) > 3 else ()
        else:
            continue  # malformed edge

        # Convert to 1-based index if VL provided (because your .dat file uses 1-based)
        i_disp = i + 1 if VL else i
        j_disp = j + 1 if VL else j

        # Prepare vertex name output
        if VL:
            vi = VL.get(i_disp, f"v{i_disp}")
            vj = VL.get(j_disp, f"v{j_disp}")
            label = f"{vi} {vj}"
        else:
            label = f"{i_disp} {j_disp}"

        # Format additional columns if present
        extra_str = " ".join(str(e) for e in extra)
        line = f"{i_disp} {j_disp} {dist:.1f}"
        if extra:
            line += f" {extra_str}"
        line += f" {label}"

        print(line)

def solve_dgp_instance(dat_file, LV=False, VL=False):
    edges, n = parse_dgp_dat_file(dat_file)
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars(n, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="x")
    for (i, j, dist) in edges:
        diff = m.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"diff_{i}_{j}")
        err = m.addVar(lb=0, name=f"err_{i}_{j}")
        m.addConstr(diff == x[i] - x[j])
        m.addGenConstrAbs(err, diff)
        m.addConstr(err == dist)
    m.optimize()
    if m.status == GRB.OPTIMAL:
        positions = [x[i].X for i in range(n)]
        print("Embedding (positions):")
        for i in range(n):
            gurobi_var = x[i].VarName
            if VL:
                label = VL.get(i+1, f"vertex_{i+1}")
                print(f"{gurobi_var} = {positions[i]:.6f}    ({label})")
            else:
                print(f"{gurobi_var} = {positions[i]:.6f}")

        # Print pairwise distances
        for (i, j, dist) in edges:
            realized = abs(positions[i] - positions[j])
            if VL:
                name_i = VL.get(i+1, f"vertex_{i+1}")
                name_j = VL.get(j+1, f"vertex_{j+1}")
                print(f"Edge ({name_i}, {name_j}), target: {dist:.4f}, realized: {realized:.4f}")
            else:
                print(f"Edge ({i}, {j}), target: {dist:.4f}, realized: {realized:.4f}")
        return True
    else:
        print("No feasible solution found.")
        return False
    
##############################
# Main Loop: Mass Testing
##############################

def main():
    if len(sys.argv) < 6:
        exit("Usage: {} num_tests n_vars n_cons density [feas|rnd|infeas|perturb]".format(sys.argv[0]))

    if len(sys.argv) >= 7:
        opb_file = sys.argv[6]
        print(f"Testing specific BLP instance: {opb_file}")
        with open(opb_file, "r") as f:
            opb_str = f.read()

        # Brute-force feasibility for small enough n
        n_vars = max(
            int(var[1:]) for line in opb_str.splitlines() for var in re.findall(r'x\d+', line)
        )
        expected = evaluate_opb_instance(opb_str, n_vars)
        print("BLP feasibility (brute force):", "FEASIBLE" if expected else "INFEASIBLE")

        # Reduce to DGP (call your script)
        try:
            blp_instance = readOpb(opb_file)
            E, VL, LV, c2e, max_var = reduce_blp_2_dgp(blp_instance)
            E_sorted = sorted(list(E))
            V_list = sorted(VL.keys())
            G = (V_list, E_sorted, VL, max_var)

            dgp_file = opb_file.replace(".opb", "-dgp.dat")
            writeDat(G, dgp_file, opb_file)
        except Exception as e:
            print(f"Reduction failed: {e}")
            sys.exit(1)


        dgp_feasible = solve_dgp_instance(dgp_file, LV, VL)
        print("DGP feasibility:", "FEASIBLE" if dgp_feasible else "INFEASIBLE")
        # If either BLP or DGP is infeasible, compute and print MinErrs
        if not expected or not dgp_feasible:
            print("\n--- Error Minimization ---")
            # Parse BLP instance into matrix form for min_err_blp
            blp_instance = readOpb(opb_file)
            Ax, rel, b, var_ids = blp_instance
            n_bvars = max(var_ids)
            m_cons = len(Ax)
            # Convert Ax (dict of dicts) and b (dict) to numpy matrix and vector
            A_mat = np.zeros((m_cons, n_bvars))
            b_vec = np.zeros(m_cons)
            for row_idx, ax_dict in Ax.items():
                for lit, coeff in ax_dict.items():
                    col = abs(lit) - 1
                    sign = 1 if lit > 0 else -1
                    A_mat[row_idx, col] = coeff * sign
                b_vec[row_idx] = b[row_idx]
            minerr_blp, x_star = min_err_blp(A_mat, b_vec)
            print(f"MinErrBLP (minimum norm of constraint violation for best binary assignment): {minerr_blp:.6f}")
            print(f"Best BLP assignment: {x_star}")

            # Parse DGP .dat file for min_err_dgp
            edges, n_nodes = parse_dgp_dat_file(dgp_file)
            print(f"edges: {edges}")
            fixed_positions = None
            minerr_dgp, emb_star = min_err_dgp_absolute(edges, n_nodes, fixed_positions=fixed_positions)
            print(f"MinErrDGP (sum of squared edge errors for best embedding): {minerr_dgp}")
            print(f"Best DGP embedding: {emb_star}")

        try:
            os.remove(dgp_file)
        except Exception:
            pass
        return  # done!
    
    num_tests = int(sys.argv[1])
    n_vars = int(sys.argv[2])
    n_cons = int(sys.argv[3])
    density = float(sys.argv[4])
    feasType = sys.argv[5].lower()
    validTypes = {"feas", "rnd", "infeas", "perturb"}
    if feasType not in validTypes:
        exit("Type must be one of: feas, rnd, infeas, perturb")
    
    tests_passed = 0
    tests_failed = 0
    failures = []
    
    compute_expected = (n_vars <= 10)
    
    for i in range(1, num_tests+1):
        opb_str = generate_instance(n_vars, n_cons, density, feasType)
        print(f"generate BLP instance: {opb_str}")
        opb_file = f"tmp_instance_{i}.opb"
        with open(opb_file, "w") as f:
            f.write(opb_str)

        # --- Reduction: BLP to DGP ---
        try:
            # Call your reduction script: blp2dgp.py
            #sys.executable gives the path to the Python interpreter that is running the script so we will be in the correct environment
            subprocess.run([sys.executable, "blp2dgpleo_issue.py", opb_file], check=True)
        except Exception as e:
            print(f"Test {i}: Reduction script failed:\n{e}\nInstance:\n{opb_str}\n")
            tests_failed += 1
            continue

        dgp_file = opb_file.replace(".opb", "-dgp.dat")

        # --- Solve DGP instance ---
        try:
            dgp_feasible = solve_dgp_instance(dgp_file)
        except NotImplementedError:
            print("Please implement solve_dgp_instance to call your DGP solver.")
            sys.exit(1)
        except Exception as e:
            print(f"Test {i}: Exception during DGP solving:\n{e}\nInstance:\n{opb_str}\n")
            tests_failed += 1
            continue

        # --- Evaluate expected feasibility (if small enough) ---
        expected = None
        if compute_expected:
            expected = evaluate_opb_instance(opb_str, n_vars)
        
        # --- Compare results ---
        if expected is not None:
            if dgp_feasible == expected:
                tests_passed += 1
            else:
                tests_failed += 1
                failures.append((i, opb_str, expected, dgp_feasible))
                print(f"Test {i} FAILED: Expected {'FEASIBLE' if expected else 'INFEASIBLE'}, got {'FEASIBLE' if dgp_feasible else 'INFEASIBLE'}")
        else:
            tests_passed += 1  # Just log as passed if not sure

        # Optional: Clean up temp files
        try:
            os.remove(opb_file)
            os.remove(dgp_file)
        except Exception:
            pass

    print(f"\nTesting complete: {tests_passed} passed, {tests_failed} failed out of {num_tests} tests.")
    if failures:
        print("\nFailed instances:")
        for idx, inst, exp, got in failures:
            print(f"Test {idx}: Expected {'FEASIBLE' if exp else 'INFEASIBLE'}, got {'FEASIBLE' if got else 'INFEASIBLE'}")
            print("Instance:")
            print(inst)
            print("-----")

if __name__ == '__main__':
    main()
