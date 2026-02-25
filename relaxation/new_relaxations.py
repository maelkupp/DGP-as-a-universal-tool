import re
import gurobipy as gp
from gurobipy import GRB
import sys
import math
import random
import subprocess
import numpy as np
from collections import defaultdict
import pandas as pd


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


def parse_dgp_dat_file(dat_file):
    edges = []
    n = 0
    in_edges = False

    vertex_names = {}  # id (0-based) -> name

    with open(dat_file, "r") as f:
        for line in f:
            if line.strip().startswith("param : E : c I :="):
                in_edges = True
                continue

            if in_edges:
                if line.strip().startswith(";"):
                    break

                # match edge data
                m = re.match(r"\s*(\d+)\s+(\d+)\s+([0-9eE\+\-\.]+)", line)
                if not m:
                    continue

                i, j, dist = int(m.group(1)), int(m.group(2)), float(m.group(3))
                i0, j0 = i - 1, j - 1
                edges.append((i0, j0, dist))
                n = max(n, i, j)

                # extract names after # [name_i,name_j]
                name_match = re.search(r"#\s*\[([^\]]+)\]", line)
                if name_match:
                    names = name_match.group(1).split(",")
                    if len(names) == 2:
                        name_i, name_j = names[0].strip(), names[1].strip()
                        vertex_names.setdefault(i0, name_i)
                        vertex_names.setdefault(j0, name_j)

    return edges, n, vertex_names

def get_oriented_edge_list(edges):
    oriented_edges = []
    for (i, j, w) in edges:
        #orient our edges so that tail = min(i,j) and head = max(i,j) to break symmetry and make it easier to interpret the output embedding
        if i < j:
            oriented_edges.append((i, j, w))
        else:
            oriented_edges.append((j, i, w))
    return oriented_edges


def build_incidence_matrix(edge_list, n):
    m = len(edge_list)
    B = np.zeros((m, n))

    for e, (tail, head, _) in enumerate(edge_list):
        B[e, tail] = -1
        B[e, head] = +1

    return B

def create_adjacency_list_from_edges(edges, n):
    adj_list = {i: [] for i in range(0, n)}
    for (i, j, w) in edges:
        adj_list[i].append((j, w)) 
        adj_list[j].append((i, w)) 

    return adj_list

import numpy as np
from collections import defaultdict

def build_adj(edge_list):
    adj = defaultdict(list)
    for e, (u, v, _) in enumerate(edge_list):
        adj[u].append((v, e))
        adj[v].append((u, e))
    return adj


def fundamental_cycle_basis(edge_list, n, B=None, verify=True):
    """
    edge_list: list of (tail, head, d)
    n: number of vertices
    B: incidence matrix (optional, for verification)
    verify: whether to check B^T c = 0
    """

    m = len(edge_list)
    adj = build_adj(edge_list)

    # --- Build spanning tree ---
    parent = [-1]*n
    parent_edge = [-1]*n
    visited = [False]*n
    tree_edges = set()

    for start in range(n):
        if visited[start]:
            continue

        stack = [start]
        visited[start] = True

        while stack:
            u = stack.pop()
            for v, e in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    parent[v] = u
                    parent_edge[v] = e
                    tree_edges.add(e)
                    stack.append(v)

    # --- Helper: recover tree path from a to b ---
    def tree_path(a, b):
        """
        Returns list of (edge_index, traversal_direction)
        traversal_direction = +1 if traversing tail→head,
                              -1 if head→tail
        """
        # compute ancestors of a
        ancestors = {}
        x = a
        while x != -1:
            ancestors[x] = True
            x = parent[x]

        # find LCA
        x = b
        while x not in ancestors:
            x = parent[x]
        lca = x

        path = []

        # path from a up to LCA
        x = a
        while x != lca:
            e = parent_edge[x]
            tail, head, _ = edge_list[e]

            # traversal is x → parent[x]
            if tail == x and head == parent[x]:
                direction = +1
            else:
                direction = -1

            path.append((e, direction))
            x = parent[x]

        # path from b up to LCA (reverse later)
        temp = []
        x = b
        while x != lca:
            e = parent_edge[x]
            tail, head, _ = edge_list[e]

            # traversal is parent[x] → x (reverse direction)
            if tail == parent[x] and head == x:
                direction = +1
            else:
                direction = -1

            temp.append((e, direction))
            x = parent[x]

        path.extend(temp)

        return path

    # --- Build fundamental cycles ---
    C_rows = []

    for e, (u, v, _) in enumerate(edge_list):
        if e in tree_edges:
            continue  # skip tree edges

        c = np.zeros(m)

        # 1) traverse chord u → v
        tail, head, _ = edge_list[e]
        if tail == u and head == v:
            c[e] = +1
        elif tail == v and head == u:
            c[e] = -1
        else:
            raise RuntimeError("Edge orientation mismatch.")

        # 2) traverse tree path v → u
        path = tree_path(v, u)

        for f, direction in path:
            c[f] += direction

        # Verification
        if verify and B is not None:
            residual = B.T @ c
            if not np.allclose(residual, 0):
                raise ValueError(
                    f"Cycle row does not satisfy B^T c = 0.\nResidual={residual}"
                )

        C_rows.append(c)

    C = np.array(C_rows)
    if len(C_rows) == 0:
        C = np.zeros((0, m))
    else:
        C = np.vstack(C_rows)
    # Final global verification
    if verify and B is not None:
        if not np.allclose(C @ B, 0):
            raise ValueError("C B ≠ 0 — cycle basis invalid.")

    return C

def compute_phi(B, v):
    m, n = B.shape
    model = gp.Model("phi_lp")

    x = model.addVars(n, lb=-GRB.INFINITY, name="x")
    t = model.addVars(m, lb=0.0, name="t")

    for e in range(m):
        expr = gp.LinExpr()
        for vtx in range(n):
            if B[e, vtx] != 0:
                expr += B[e, vtx] * x[vtx]

        model.addConstr(t[e] >=  expr - v[e])
        model.addConstr(t[e] >= -expr + v[e])

    model.setObjective(gp.quicksum(t[e] for e in range(m)), GRB.MINIMIZE)
    model.optimize()

    return model.ObjVal


def cycle_max_relaxation(C, d, verbose=False):
    """
    C: r x m cycle matrix
    d: length m numpy array
    Returns:
        lower bound value
        optimal s vector
    """
    if C.shape[0] == 0:
        return 0.0, np.ones(len(d))

    r, m = C.shape

    model = gp.Model("cycle_max_relax")
    if not verbose:
        model.setParam("OutputFlag", 0)

    # Binary encoding: s_e = 2 y_e - 1
    y = model.addVars(m, vtype=GRB.BINARY, name="y")
    t = model.addVar(lb=0.0, name="t")

    # Build constraints
    for k in range(r):
        expr = gp.LinExpr()
        for e in range(m):
            if C[k, e] != 0:
                coeff = C[k, e] * d[e]
                expr += coeff * (2*y[e] - 1)

        model.addConstr(t >=  expr)
        model.addConstr(t >= -expr)

    model.setObjective(t, GRB.MINIMIZE)
    model.optimize()

    if model.status != GRB.OPTIMAL:
        raise RuntimeError("Cycle-max relaxation did not solve optimally.")

    s_opt = np.array([2*y[e].X - 1 for e in range(m)])

    return model.ObjVal, s_opt

def defect_norm_relaxation(C, d, verbose=False):
    """
    C: r x m cycle matrix
    d: length m numpy array
    Returns:
        lower bound value
        optimal s vector
    """
    if C.shape[0] == 0:
        return 0.0, np.ones(len(d))

    r, m = C.shape

    model = gp.Model("defect_norm_relax")
    if not verbose:
        model.setParam("OutputFlag", 0)

    y = model.addVars(m, vtype=GRB.BINARY, name="y")
    z = model.addVars(r, lb=-GRB.INFINITY, name="z")
    u = model.addVars(r, lb=0.0, name="u")

    # Define z_k
    for k in range(r):
        expr = gp.LinExpr()
        for e in range(m):
            if C[k, e] != 0:
                expr += C[k, e] * d[e] * (2*y[e] - 1)

        model.addConstr(z[k] == expr)

    # Absolute value linearization
    for k in range(r):
        model.addConstr(u[k] >=  z[k])
        model.addConstr(u[k] >= -z[k])

    model.setObjective(gp.quicksum(u[k] for k in range(r)), GRB.MINIMIZE)
    model.optimize()

    if model.status != GRB.OPTIMAL:
        #still need to return a value
        LB = model.ObjVal if model.ObjVal is not None else 0
        return LB, np.array([2*y[e].X - 1 if y[e].X is not None else 1 for e in range(m)])

    # operator norm ||C||_{1→1}
    norm_C = np.max(np.sum(np.abs(C), axis=0))

    s_opt = np.array([2*y[e].X - 1 for e in range(m)])

    LB = model.ObjVal / norm_C

    return LB, s_opt


def min_err_dgp1_gurobi(edges, n, time_limit=240, mip_gap=0.01):
    m = gp.Model("MinErrDGP1")
    m.Params.OutputFlag = 0

    m.Params.TimeLimit = time_limit
    m.Params.MIPGap    = mip_gap
    m.Params.Heuristics = 0.7
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions  = 10

    x = m.addVars(n, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="x")
    
    #add the following constraints to speed up processing by removing symmetries
    m.addConstr(x[0] == 0) #fix the first vertex at 0 to remove translation symmetry
    m.addConstr(x[1] >= 0) #fix the second vertex to be non-negative to remove reflection symmetrys


    err = {}
    for idx, (i, j, d) in enumerate(edges):
        #print(f"Looking at edge idx {idx} between vertex {i} and vertex {j} with distance {d}", file=sys.stderr)
        diff = m.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY)
        m.addConstr(diff == x[i] - x[j])

        absdiff = m.addVar(lb=0)
        m.addGenConstrAbs(absdiff, diff)

        aux = m.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY)
        m.addConstr(aux == absdiff - d)

        err[idx] = m.addVar(lb=0)
        m.addGenConstrAbs(err[idx], aux)

    m.setObjective(gp.quicksum(err[idx] for idx in err), GRB.MINIMIZE)
    m.optimize()

    if m.Status in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL):
        embedding = [x[i].X for i in range(n)]
        total_err = sum(abs(abs(embedding[i] - embedding[j]) - d)
                        for i, j, d in edges)
        return total_err, embedding
    else:
        return None, None

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python minErrDGP1.py <dgp_instance.dat>", file=sys.stderr)
        sys.exit(1)

    N = 15
    weights = [(1,5), (5,10), (10,50)]
    n_vertices = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
    densities = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    
    results = []
    for num in n_vertices:
        for density in densities:
            m = int(density * num*(num-1)/2)
            for weight_min, weight_max in weights:
                for _ in range(N):
                    generate_dat_graph(sys.argv[1], n=num, m=m, weight_min=weight_min, weight_max=weight_max)

                    edges, n_parsed, vertex_names = parse_dgp_dat_file(sys.argv[1])

                    oriented_edges = get_oriented_edge_list(edges)
                    B = build_incidence_matrix(oriented_edges, num)
                    C = fundamental_cycle_basis(oriented_edges, num, B=B, verify=True)
                    d = np.array([edge[2] for edge in oriented_edges])

                    lb1, s1 = cycle_max_relaxation(C, d)
                    print("Cycle-max LB:", lb1)

                    lb2, s2 = defect_norm_relaxation(C, d)
                    print("Defect-norm LB:", lb2)

                    lb3_proc = subprocess.run(
                        ["relaxation\\minErrDGP1_cycles_LB.exe", sys.argv[1], "0"],
                        check=True, text=True, capture_output=True)
                    cycle_val, cycle_time = map(float, lb3_proc.stdout.split())
                    print("Greedy edge disjoint cycle:", cycle_val)

                    lb4_proc = subprocess.run(
                        ["relaxation\\greedy_packing_cycle_LB.exe", sys.argv[1]],
                        check=True, text=True, capture_output=True)
                    greedy_val, greedy_time = map(float, lb4_proc.stdout.split())
                    print("Greedy cycle packing:", greedy_val)
                    print(f"\n")

                    results.append({
                        "n": num,
                        "m": m,
                        "density": density,
                        "weight_min": weight_min,
                        "weight_max": weight_max,
                        "cycle_max": lb1,
                        "defect_norm": lb2,
                        "cycle_val": cycle_val,
                        "greedy_val": greedy_val,
                        })


    df = pd.DataFrame(results)
    df.to_csv("lower_bounds.csv", index=False)
