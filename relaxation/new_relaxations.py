import re
import gurobipy as gp
from gurobipy import GRB
import sys
import math
import random
import numpy as np

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
