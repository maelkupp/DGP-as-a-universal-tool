
# 🧠 DGP as a Universal Computational Framework

Welcome to this repository, which contains all the scripts developed during my research project on using the **Distance Geometry Problem (DGP)** as a universal tool for computation. The goal of this research was to investigate whether the DGP could be leveraged to model and solve **Binary Linear Programming (BLP)** problems via geometric transformations.

---

## 📌 Project Goals

- Explore whether the **BLP problem** can be reduced to the **DGP problem**.
- Leverage the **geometric structure** of DGP to apply **relaxations** and **continuous optimization techniques**.
- Investigate **error mapping** from DGP back to BLP in approximate/infeasible cases.

---

## 📘 Definitions

### 🔷 Binary Linear Programming (BLP)
> Given a cost vector `c`, matrix `A`, and bound vector `b`, the feasibility version of the BLP is defined as:
```
Find x ∈ {0,1}^n  such that A·x ≤ b
```
This problem is **NP-hard** and appears in various optimization settings like scheduling, logistics, and allocation.

### 🔷 Distance Geometry Problem (DGP)
> Given a graph G = (V, E, d) and dimension `K`, determine whether it's possible to assign coordinates to each vertex such that:
```
∀ {u, v} ∈ E,  ||y_u − y_v||² = d_uv²
```
DGP naturally allows **geometric relaxations** by embedding in higher-dimensional spaces.

---

## ⚙️ What’s in This Repository

This repository includes:
- 🔁 A **reduction from BLP to SAT**, based on the **linear-time method by J.P. Warners** ([6] in report).
- 🔗 A **reduction from 3SAT to DGP**, based on **Saxe’s 1979 construction** ([5] in report), along with my formal proof of correctness.
- 🧩 A **direct reduction from BLP to DGP**, developed during this research, with several variants implemented.
- 🧠 A **dynamic programming solution to the Subset Sum problem**, used in early modeling experiments.

---

## 🧪 Implemented Reductions

### ✅ BLP → SAT
- Implements Warner’s [6] linear-time reduction.
- Converts BLP constraints to CNF using normalization, binary decomposition, and logical encoding.

### ✅ SAT → 3SAT → DGP
- Uses **Saxe's graph gadgets** to encode 3SAT clauses as DGP constraints.
- I provide a **formal embeddability proof** to complete the theoretical foundation.

### ✅ Direct BLP → DGP
- Designed a **custom graph encoding** for BLP constraints directly into DGP instances.
- Supports **error-based relaxation** and provides a path toward computing `MinErrBLP` via `MinErrDGP`.

### 🔄 BLP → DGP2 (relaxed)
- A second variant lifts some DGP embeddings into higher dimensions to avoid infeasibility artifacts.
- Models constraints **without error**, using additional degrees of freedom (e.g., DGP2 embeddings with partial fixing).

---

## 📂 Repository Structure

```
/scripts/
├── blp_to_sat.py              # Warner’s reduction
├── sat_to_dgp.py              # DGP graph construction from 3SAT
├── blp_to_dgp_direct.py       # Direct encoding of BLP into DGP
├── subset_sum_dp.py           # Dynamic programming solver
├── utils/                     # Graph utilities, CNF encoders, etc.
```

---

## 📄 References

1. **Karp, R. M.** "Reducibility Among Combinatorial Problems", 1972.  
2. **Saxe, J.** "Embeddability of Weighted Graphs in k-Space is Strongly NP-Hard", 1979.  
3. **Warners, J.P.** "A Linear-Time Transformation of Linear Inequalities into CNF", *Information Processing Letters*, 1998.  
4. See full list in the [report](./PRL_report.pdf) 📄.

---

## 📈 Future Work

- Adapt and test the **Branch and Prune algorithm** for solving generated DGP instances.
- Develop robust **error mapping techniques** for translating MinErrDGP back to BLP.
- Explore more compact graph gadgets to enable **efficient reductions without instance blow-up**.

---

## 🤝 Acknowledgments

This project was conducted under the supervision of **Prof. Leo Liberti** at **École Polytechnique**, as part of my undergraduate research thesis.
