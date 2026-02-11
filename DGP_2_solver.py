import re
import gurobipy as gp
from gurobipy import GRB
import sys

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


def solve_dgp2_gurobi(dat_file, time_limit=60, verbose=False):
    """
    Input:
      dat_file: path to .dat DGP instance

    Output:
      (feasible, embedding)
      embedding: dict i -> (x_i, y_i) if feasible, else None
    """
    edges, n, vertex_names = parse_dgp_dat_file(dat_file)

    model = gp.Model("DGP2")
    model.Params.NonConvex = 2
    model.Params.OutputFlag = 1 if verbose else 0
    model.Params.TimeLimit = time_limit

    # variables
    x = model.addVars(n, lb=-GRB.INFINITY, name="x")
    y = model.addVars(n, lb=-GRB.INFINITY, name="y")

    # gauge fixing using first edge
    i0, j0, d0 = edges[0]
    model.addConstr(x[i0] == 0.0)
    model.addConstr(y[i0] == 0.0)
    model.addConstr(y[j0] == 0.0)
    model.addConstr(x[j0] == d0)

    # distance constraints
    for i, j, d in edges:
        model.addQConstr(
            (x[i] - x[j]) * (x[i] - x[j]) +
            (y[i] - y[j]) * (y[i] - y[j]) == d * d
        )

    model.setObjective(0.0, GRB.MINIMIZE)
    model.optimize()

    if model.Status not in (GRB.OPTIMAL, GRB.SUBOPTIMAL):
        return False, None

    embedding = {i: (x[i].X, y[i].X) for i in range(n)}
    return True, embedding


if __name__ == "__main__":
    dat_file = sys.argv[1]
    feasible, embedding = solve_dgp2_gurobi(dat_file, time_limit=300, verbose=True)
    if feasible:
        print("Feasible embedding found:")
        for i in range(len(embedding)):
            print(f"Vertex {i}: (x={embedding[i][0]}, y={embedding[i][1]})")
    else:
        print("No feasible embedding found within the time limit.")

