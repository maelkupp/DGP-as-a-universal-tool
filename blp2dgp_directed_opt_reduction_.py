################ IMPORTS ##################

import sys
import os
import time
import subprocess
import datetime
import mimetypes
import re
import pickle
import math
import numpy as np
from collections import defaultdict
import math

####################### GLOBALS ###############################

myEps = 1e-2
myZero = 1e-9

satHeader = "p cnf"
verbose = 1 
verbose = 0 # quiet

ts = time.time()
timestamp = datetime.datetime.fromtimestamp(ts).isoformat().replace('-','').replace('T','').replace(':','').split('.')[0][2:]
exn = sys.argv[0]

##################### FUNCTIONS ###############################3
## read an OPB file Ax rel b (assume integer coeffs)
def readOpb(opbf):
    b = {}    # RHS
    rel = {}  # Relations
    Ax = {}   # LHS
    var_ids = set()
    with open(opbf, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line[0] == '#':
                continue
            l = line.split()
            if l[-1] == ';':
                b[i] = int(l[-2])
                rel[i] = l[-3]
                l = l[:-3]
            else:
                b[i] = int(l[-1][:-1])
                rel[i] = l[-2]
                l = l[:-2]
            Ax[i] = {}
            for term in l:
                coeff, var = term.split('*')
                j = int(var[1:])  # 'x3' → 3
                coeff = int(coeff)
                Ax[i][j] = coeff
                var_ids.add(abs(j))
    return (Ax, rel, b, var_ids)


## write an AMPL .dat file encoding the DGP by list of edges
def writeDat(G, dgpf, opbf, e=exn, t=timestamp):
    (V,E, VL, nlits) = G
    n = max(V)
    with open(dgpf, "w") as f:
        print("# DGP instance in .dat format for dgp.mod", file=f)
        print("# written by {} on {}".format(e, t), file=f)
        print("# this graph is a DGP encoding of the BLP instance", file=f)
        print("#   found in {}".format(opbf), file=f)
        print("param n := {};".format(n), file=f)
        print("param nlits := {};".format(nlits), file=f)
        print("param Kdim := 1;", file=f)
        print("param : E : c I :=", file=f)

        for el in E:
            print(el[0], el[1], el[2], el[3], VL[el[0]], VL[el[1]])
            print("  {} {}  {:.3f} {}  # [{},{}]".format(el[0],el[1],el[2],el[3],VL[el[0]],VL[el[1]]), file=f)
        print(";", file=f)
        #print("# vertex map", file=f)
        #for i in V:
        #    print("# vertex {} is {}".format(i, VL[i]), file=f)
    return

def normalizing_single_constraint(Ax_i, b, rel):
    """
    Normalizes a single BLP constraint so that:
      - The inequality is in the form "≤".
      - All coefficients are positive.
    If the original relation is "≥", the function multiplies both sides by -1.
    Additionally, if any coefficient is negative, it flips the corresponding literal
    and adds the absolute value of the coefficient to the bound.
    """

    norm_Ax = {}
    norm_b = b
    relation = rel

    if relation == ">=":
        # Multiply both sides by -1
        relation = "<="
        for var, coeff in Ax_i.items():
            norm_Ax[var] = -coeff
        norm_b = -norm_b
    else:
        norm_Ax = dict(Ax_i)  # copy

    #make sure all the coefficients are positive
    norm_Ax_final = {}
    for xj, c in norm_Ax.items():
        if c < 0:
            # Flip variable sign if coeff is negative
            c = abs(c)
            xj = -xj #make the literal a 'not' literal
            
            norm_b += c #add the coefficient over to the other side
        norm_Ax_final[xj] = c
    
    return (norm_Ax_final, norm_b, relation)

def normalize_blp(blp_instance):
    #obtain a blp_instance of form (Ax, relation, b, var)
    #return an equivalent blp_instance where all the coefficients are positive and all the relations are <=
    # the equality constraint is written as <= and a >= constraint
    (Ax, relation, b, var) = blp_instance
    Ax_norm = {}
    rel_norm = {}
    b_norm = {}
    
    new_index = 0  # Use a new index to store normalized constraints

    for i in Ax:
        norm_Axi, norm_b, norm_rel = normalizing_single_constraint(Ax[i], b[i], relation[i])

        if norm_rel == "=":
            # 1. ∑ aᵢ·xᵢ ≤ b
            Ax_norm[new_index] = dict(norm_Axi)
            rel_norm[new_index] = "<="
            b_norm[new_index] = norm_b
            new_index += 1

            # 2. ∑ (-aᵢ)·xᵢ ≤ -b, with proper literal sign flipping
            flipped_Ax = {}
            flipped_b = -norm_b
            for var, coeff in norm_Axi.items():
                flipped_coeff = -coeff
                if flipped_coeff < 0:
                    #print('landed here')
                    flipped_coeff = abs(flipped_coeff)
                    var = -var
                    flipped_b += flipped_coeff
                flipped_Ax[var] = flipped_coeff
                #print(f"here, flipped_Ax: {flipped_Ax}")
            Ax_norm[new_index] = flipped_Ax
            rel_norm[new_index] = "<="
            b_norm[new_index] = flipped_b
            new_index += 1

        else:
            # Otherwise, just store the constraint as-is
            Ax_norm[new_index] = norm_Axi
            rel_norm[new_index] = norm_rel
            b_norm[new_index] = norm_b
            new_index += 1
    
    return (Ax_norm, rel_norm, b_norm, var)

#################### Main ###########################
opbf = sys.argv[1] #the first argument given to the program
opbbn = os.path.basename(opbf)  # basename
dgpf = ''.join(opbbn.split('.')[:-1]) + "-dgp.dat"  # output DGP to AMPL .dat
blp2dgpf = ''.join(opbbn.split('.')[:-1]) + "-blp2dgp.pkl"  # output DGP to AMPL .dat

blp_instance = readOpb(opbf)
(Ax, rel, b, var) = normalize_blp(blp_instance)

max_var = max(var) # here we obtain the number of literals in our BLP instance
nconstraints = len(Ax) #the number of DGP contraints we have to map to

c2e = {}
E = set([])  # weighted edges of the DGP graph
VL = {} #vertex_id to literal name
LV = {} #literal name to vertex_id
vertex_id = 1
## literal vertices
for i in range(1,max_var+1):
    VL[vertex_id] = "X+{}".format(i)
    LV[VL[vertex_id]] = vertex_id
    vertex_id += 1
## A vertices
VL[vertex_id] = "A"
LV[VL[vertex_id]] = vertex_id
vertex_id += 1

E.update([(LV["A"], LV["X+{}".format(i)], 1.0, 0) for i in range(1,max_var+1)]) # add the weight from the vertex A to each vertex encoding the literals

for index in Ax:
    #iterating over each constraint in the BLP
    print(Ax[index])

    #################### ENCODING THE PARTIAL SUM #############################
    VL[vertex_id] = "P_0^{}".format(index)
    LV[VL[vertex_id]] = vertex_id
    vertex_id += 1
    edges_to_add = [(LV["A"], LV["P_0^{}".format(index)], 1, 0)] # initiating our edges to add list with the edge from A to the start of the partial sum for this BLP 
    VL[vertex_id] = "M_1^{}".format(index)
    LV[VL[vertex_id]] = vertex_id
    vertex_id += 1
    edges_to_add.append((LV["P_0^{}".format(index)], LV["M_1^{}".format(index)], (Ax[index].get(1, 0))/2, 1))

    for i in range(1, max_var):
        VL[vertex_id] = "P_{}^{}".format(i, index)
        LV[VL[vertex_id]] = vertex_id
        vertex_id += 1

        VL[vertex_id] = "M_{}^{}".format(i+1, index)
        LV[VL[vertex_id]] = vertex_id
        vertex_id += 1
        edges_to_add.append((LV["M_{}^{}".format(i, index)], LV["P_{}^{}".format(i,index)], Ax[index].get(i, 0)/2, 0))
        edges_to_add.append((LV["P_{}^{}".format(i, index)], LV["M_{}^{}".format(i+1,index)], Ax[index].get(i+1, 0)/2, 1))

    VL[vertex_id] = "P_{}^{}".format(max_var, index)
    LV[VL[vertex_id]] = vertex_id
    vertex_id += 1
    edges_to_add.append((LV["M_{}^{}".format(max_var, index)], LV["P_{}^{}".format(max_var,index)], Ax[index].get(max_var, 0)/2, 0))
    #here we have created the system of partial sums we now need to encode the buffer

    ############ ENCODING THE BUFFER ##########################################
    l = math.floor(math.log2(b[index])) #obtain the floor of log base 2 of our bound to encode the buffer
    VL[vertex_id] = "R_{}^{}".format(0, index)
    LV[VL[vertex_id]] = vertex_id
    vertex_id += 1

    edges_to_add.append((LV["P_{}^{}".format(max_var, index)], LV["R_{}^{}".format(0,index)], 1/2, 1))

    for i in range(0, l+1):
        VL[vertex_id] = "S_{}^{}".format(i, index)
        LV[VL[vertex_id]] = vertex_id
        vertex_id += 1

        VL[vertex_id] = "R_{}^{}".format(i+1, index)
        LV[VL[vertex_id]] = vertex_id
        vertex_id += 1

        edges_to_add.append((LV["R_{}^{}".format(i, index)], LV["S_{}^{}".format(i,index)], 2**(i-1), 0))
        edges_to_add.append((LV["S_{}^{}".format(i, index)], LV["R_{}^{}".format(i+1,index)], 2**i, 1))
    
    VL[vertex_id] = "S_{}^{}".format(l, index)
    LV[VL[vertex_id]] = vertex_id
    vertex_id += 1
    edges_to_add.append((LV["R_{}^{}".format(l, index)], LV["S_{}^{}".format(l,index)], 2**(l-1), 0))
    #finished encoding the buffer

    ########### adding the auxiliary variables #############
    ############ connect the x_1 variable to the partial sum #####################
    edges_to_add.append((LV["A"], LV["P_{}^{}".format(1, index)], Ax[index].get(1, 0), 0))
    edges_to_add.append((LV["X+{}".format(1)], LV["P_{}^{}".format(1,index)], 2, 0))

    #connect the rest of the literal
    for j in range(2, max_var):
        if Ax[index].get(j, 0) != 0:
            #only do this for the edges that have a weight larger than 1
            VL[vertex_id] = "B_{}^{}".format(j, index)
            LV[VL[vertex_id]] = vertex_id
            vertex_id += 1

            VL[vertex_id] = "x_{}^{}".format(j, index)
            LV[VL[vertex_id]] = vertex_id
            vertex_id += 1
            
            ######### mechanism is 1D that forces an error ###################
            edges_to_add.append((LV["B_{}^{}".format(j, index)], LV["P_{}^{}".format(j-1,index)], 1, 1))
            edges_to_add.append((LV["B_{}^{}".format(j, index)], LV["P_{}^{}".format(j,index)], Ax[index].get(j, 0), 0))

            edges_to_add.append((LV["x_{}^{}".format(j, index)], LV["B_{}^{}".format(j,index)], 1, 0))
            edges_to_add.append((LV["x_{}^{}".format(j, index)], LV["P_{}^{}".format(j,index)], 2, 0))

        
        ########### ENCODING THE LITERAL EXTENSION MECHANISM ##################################
            VL[vertex_id] = "N_({}, {})^{}".format(j, 1, index)
            LV[VL[vertex_id]] = vertex_id
            vertex_id += 1

            edges_to_add.append((LV["X+{}".format(j)], LV["N_({}, {})^{}".format(j, 1,index)], Ax[index].get(1, 0)/2, 1))
            for k in range(1, j-1):
                VL[vertex_id] = "Q_({}, {})^{}".format(j, k, index)
                LV[VL[vertex_id]] = vertex_id
                vertex_id += 1

                VL[vertex_id] = "N_({}, {})^{}".format(j, k+1, index)
                LV[VL[vertex_id]] = vertex_id
                vertex_id += 1

                edges_to_add.append((LV["N_({}, {})^{}".format(j, k, index)], LV["Q_({}, {})^{}".format(j, k,index)], Ax[index].get(k, 0)/2, 0))
                edges_to_add.append((LV["Q_({}, {})^{}".format(j, k, index)], LV["N_({}, {})^{}".format(j, k+1,index)], Ax[index].get(k+1, 0)/2, 1))
            
            edges_to_add.append((LV["N_({}, {})^{}".format(j, j-1, index)], LV["x_{}^{}".format(j,index)], Ax[index].get(j-1, 0)/2, 0))

        print(edges_to_add)
        E.update(edges_to_add)
        c2e[index] = edges_to_add



## make a graph G=(V,E)
E = sorted(list(E))
V = [vid for vid in VL]
V.sort()
G = (V,E,VL,max_var)

## write DGP1 instance
writeDat(G, dgpf, opbf)

## write c2e dictionary
with open(blp2dgpf, 'wb') as f:
    pickle.dump(c2e, f)

