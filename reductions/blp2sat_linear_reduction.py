
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
#from pysat.solvers import Solver



################ GLOBALS ##################



################ Classes ###################
class CNFEncoder:
    def __init__(self):
        self.var_counter = 0
        self.var_map = {}  # maps (name, id) → variable number
        self.clauses = []
        #list of clauses where a clause is a list of integers
        #[1, -2, 3] means x1 \/ -x2 \/ x3

    def new_var(self, name="v"):
        self.var_counter += 1
        return self.var_counter

    def add_clause(self, *literals):
        #adding a clause to the CNF formula
        #to represent x1 \/ -x2, do encoder.add_clause(1, -2)
        self.clauses.append(list(literals))

    def add_and(self, out, in1, in2):
        # out <=> in1 ∧ in2
        self.add_clause(-in1, -in2, out)
        self.add_clause(in1, -out)
        self.add_clause(in2, -out)

    def add_or(self, out, in1, in2):
        # out <=> in1 ∨ in2
        self.add_clause(in1, in2, -out)
        self.add_clause(-in1, out)
        self.add_clause(-in2, out)

    def add_eq(self, out, var):
        if out == var:
            # Don't add any clauses if identical, since out ↔ var is trivially true
            return
        if out == -var:
            # out ↔ ¬var: must explicitly add clauses
            self.add_clause(-out, -var)
            self.add_clause(out, var)
            return
        self.add_clause(-out, var)
        self.add_clause(out, -var)


################ FUNCTIONS ##################

#function helpful for debugging when I was still testing the code

## write a DIMACS-formatted .sat file
def writeSat(sat, fn):
    with open(fn, "w") as f:
        nlits = max(max(abs(l) for l in c) for c in sat)
        nclauses = len(sat)
        print("p cnf {} {}".format(nlits, nclauses), file=f)
        for clause in sat:
            for lit in clause:
                print("{} ".format(lit), end='', file=f)
            print('0', file=f)
    return


def format_clause(clause):
    """
    Formats a clause like [1, -2, 3] into a string: x1 ∨ ¬x2 ∨ x3
    """
    literals = []
    for lit in clause:
        var_name = f"x{abs(lit)}"
        if lit < 0:
            literals.append(f"¬{var_name}")
        else:
            literals.append(var_name)
    return " ∨ ".join(literals)

#function that prints which clauses have been encoded for which initial BLP constraint
def print_clauses_by_constraint(encoder, c2c):
    """
    Prints each group of CNF clauses by BLP constraint ID using the encoder's clause list.
    """
    print("\n===== CNF Clauses Grouped by BLP Constraint =====\n")
    for cid in sorted(c2c):
        print(f"Constraint {cid}:")
        for clause_idx in c2c[cid]:
            clause = encoder.clauses[clause_idx]
            print("   ", format_clause(clause))
        print()

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


##################### Functions to encode the BLP into a SAT

#obtained normalized but equivalent blp instance
#break coefficients into binary bits, create bit label variable for each a_j * x_j term
#sum them recursively, add constraints that enforce sum <= b
def xor2(encoder, a, b, out):
    #implementing a Xor b into the CNF
    encoder.add_clause(-a, -b, -out)
    encoder.add_clause(-a, b, out)
    encoder.add_clause(a, -b, out)
    encoder.add_clause(a, b, -out)

def xor3(encoder, a, b, c, out):
    #implement a Xor b Xor c into out in the CNF
    tmp = encoder.new_var("tmp_xor")
    xor2(encoder, a, b, tmp)
    xor2(encoder, tmp, c, out)

def majority3(encoder, a, b, c, carry):
    # carry is true if at least two of a, b, c are true
    encoder.add_clause(-a, -b, carry)
    encoder.add_clause(-a, -c, carry)
    encoder.add_clause(-b, -c, carry)
    encoder.add_clause(a, b, -carry)
    encoder.add_clause(a, c, -carry)
    encoder.add_clause(b, c, -carry)



def warners_leq(weights, vars_, bound, encoder):
    """
    Encode the inequality: sum_i (weights[i] * vars_[i]) <= bound
    using Warner’s method.
    
    We can assume that:
      - The BLP has been normalized so that all coefficients (weights) are positive.
      - Negative coefficients have been handled by flipping the corresponding literal and shifting the bound.
      - The 'vars_' list contains signed literals in the standard SAT convention meaning that -n corresponds to (not)x_n
    
    Returns:
      The list of CNF clauses (appended into encoder) representing the inequality.
    """

    # If the bound is negative, the constraint is unsatisfiable.
    
    if bound < 0:
        # Compute constant K as the sum of weights.
        K = sum(weights)
        # Shift the bound to a non-negative number.
        new_bound = bound + K
        # Create a new variable that will act as a constant 1.
        const_var = encoder.new_var("const")
        # Force this variable to be true.
        encoder.add_clause(const_var)
        # Append the constant term to the weights and the corresponding literal.
        weights.append(K)
        vars_.append(const_var)
        # Update bound to the shifted value.
        bound = new_bound


    # If the sum of coefficients is <= bound, the constraint is trivially satisfied.
    if sum(weights) <= bound:
        return encoder.clauses

    max_sum = sum(weights)
    M = max_sum.bit_length()  # Number of bits needed to represent the sum

    # ----------------------------------------------------------------
    # Step 1. Build the binary representation of each term.
    # For each term weight * literal, create new auxiliary bit variables.
    # The auxiliary variable for bit position k will be true if the kth bit
    # of the term is 1.
    # ----------------------------------------------------------------
    sums_by_level = defaultdict(list)
    for weight, lit in zip(weights, vars_):
        for k in range(M):
            if (weight >> k) & 1:
                # Create a new variable to represent the kth bit of this term.
                # We use a name that encodes the original literal and the bit position.
                bit_var = encoder.new_var(f"bit_{lit}_{k}")
                # Enforce that bit_var is equivalent to the original literal.
                # (Assuming that the literal already encodes any necessary negation.)
                encoder.add_eq(bit_var, lit)
                sums_by_level[k].append(bit_var)

    # ----------------------------------------------------------------
    # Step 2. Compute the binary sum of all terms.
    # We add the bits level by level, using pair–wise (binary) addition.
    # To avoid duplicate handling (which might otherwise lead to clauses like
    # (x ∨ x) or (¬x ∨ ¬x)), we “cancel out” duplicate occurrences.
    # ----------------------------------------------------------------
    carry_bits = defaultdict(list)
    sum_bits = {}  # maps bit position to auxiliary variable representing that bit of the total sum

    # For each bit level k, combine the bits from the binary representations of the terms
    # along with any carry bits from lower levels. We pair identical bits (since x XOR x = 0)
    # and add their AND as a carry to the next level.
    for k in range(M + int(math.ceil(math.log2(len(weights)))) + 2):
        current_bits = sums_by_level[k] + carry_bits[k]

        # Remove duplicates: if the same literal appears multiple times, process in pairs.
        bit_counts = defaultdict(int)
        unique_bits = []
        for bit in current_bits:
            bit_counts[bit] += 1
        for bit, count in bit_counts.items():
            pairs, remainder = divmod(count, 2)
            # For each pair: XOR(bit, bit) = 0, and AND(bit, bit) = bit gets carried.
            for _ in range(pairs):
                carry_bits[k+1].append(bit)
            if remainder:
                unique_bits.append(bit)
        current_bits = unique_bits

        # Now, while more than one bit remains, add them pairwise.
        while len(current_bits) > 1:
            a = current_bits.pop()
            b = current_bits.pop()
            sum_var = encoder.new_var(f"sum_{k}_{len(current_bits)}")
            carry_var = encoder.new_var(f"carry_{k}_{len(current_bits)}")
            # sum_var = a XOR b, carried bit = a AND b
            xor2(encoder, a, b, sum_var)
            encoder.add_and(carry_var, a, b)
            current_bits.append(sum_var)
            carry_bits[k+1].append(carry_var)
        if current_bits:
            sum_bits[k] = current_bits[0]
        else:
            #if all the bits cancel out explicitly create a zero bit at that position
            zero_bit = encoder.new_var(f"zero_bit_{k}")
            encoder.add_clause(-zero_bit)  # Force zero_bit to be false.
            sum_bits[k] = zero_bit      

    # ----------------------------------------------------------------
    # Step 3. Enforce the inequality on the binary sum.
    # Let bound_bits be the binary representation of the bound.
    # For each bit position k, if the kth bit of the bound is 0,
    # then if all higher bits exactly match, the kth sum bit must be 0.
    # Additionally, for any sum bit positions above the bound's length,
    # force those bits to be false.
    # ----------------------------------------------------------------
    max_level = max(sum_bits.keys())
    bound_bits = [(bound >> k) & 1 for k in range(max_level + 1)]
    higher_bits_match = []
    for k in reversed(range(len(bound_bits))):
        sum_bit = sum_bits.get(k)
        b_bit = bound_bits[k]
        if sum_bit is None:
            continue
        if b_bit == 0:
            if higher_bits_match:
                match_var = encoder.new_var(f"match_above_{k}")
                for hb in higher_bits_match:
                    encoder.add_clause(-match_var, hb)
                encoder.add_clause(match_var, *[-hb for hb in higher_bits_match])
                encoder.add_clause(-match_var, -sum_bit)
            else:
                # If there are no higher bits to check, then force this sum bit to 0.
                encoder.add_clause(-sum_bit)
        # For each bit position, create a “match” variable that is true if the sum bit equals the bound bit.
        match_current = encoder.new_var(f"match_bit_{k}")
        if b_bit == 1:
            encoder.add_eq(match_current, sum_bit)
        else:
            encoder.add_eq(match_current, -sum_bit)
        higher_bits_match.append(match_current)
    # Force any sum bits beyond the bound's length to be false.
    #for k in range(len(bound_bits), max(sum_bits.keys()) + 1):
    #    if k in sum_bits:
    #        encoder.add_clause(-sum_bits[k])

            
    return encoder.clauses


def encode_constraints(Ax, rel, b, encoder):
    """
    Encodes each normalized constraint (already in ≤ form with positive coefficients and signed literals)
    using Warner's method. Appends CNF clauses into the encoder.
    """
    contraint2clause = {}
    for i in Ax:
        coeffs = list(Ax[i].values())     # [a₁, a₂, ...]
        vars_ = list(Ax[i].keys())        # [±x₁, ±x₂, ...]
        print(coeffs)
        print(vars_)
        
        if rel[i] != "<=":
            raise ValueError(f"Expected <= constraint at index {i}, got: {rel[i]}")
        
        # Encode the constraint: ∑ aᵢ * xᵢ ≤ b
        #print(f"i: {i}, Ax[i]: {Ax[i]}")
        start_idx = len(encoder.clauses)
        warners_leq(coeffs, vars_, b[i], encoder)
        end_idx = len(encoder.clauses)

        contraint2clause[i] = list(range(start_idx, end_idx))
    return contraint2clause

#################### Main ###########################
opbf = sys.argv[1] #the first argument given to the program
satbn = os.path.basename(opbf)  # basename
sat3f = '.'.join(satbn.split('.')[:-1]) + ".3sat"  # 3sat name
blp_instance = readOpb(opbf) #obtain initial blp instance from the file
(Ax, rel, b, var) = normalize_blp(blp_instance) 

print(f"Normalized Ax: {Ax} \n")
print(f"Normalized bounds: {b} \n")
print(f"max_var: {max(var)}")
           



encoder = CNFEncoder()
# Reserve original variable indices before adding auxiliary variables.
original_vars = set(abs(v) for Ax_i in Ax.values() for v in Ax_i)

if original_vars:
    encoder.var_counter = max(original_vars)
constraint2clause = encode_constraints(Ax, rel, b, encoder)
writeSat(encoder.clauses, sat3f)

#print_clauses_by_constraint(encoder, constraint2clause)



'''
########## main code for solving the SAT instance
with Solver(name='glucose3') as solver:
    for clause in encoder.clauses:
        #print(f"Adding clause: {clause}")
        solver.add_clause(clause)

    satisfiable = solver.solve()
    model = solver.get_model() if satisfiable else None

print(model)'''

#To run this script call:
#python3 blp2satmael.py opbfile.opb
