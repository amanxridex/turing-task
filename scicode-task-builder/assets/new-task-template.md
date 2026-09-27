# New SciCode Task — Planning Worksheet (Path B)

Optional. Fill this in yourself if you want to think through the architecture before handing the paper to Claude — Claude can also design this directly from the paper without it. Once filled in, paste it into the conversation along with the paper.

---

## Paper Information

**Paper Title:**

**Authors:**

**Publication Year:**

**Link/Citation:**

**Open Access?** ☐ Yes ☐ No

---

## Task Concept

**Core physics problem** (2–3 sentences):

**Why suitable for SciCode?**
- ☐ Requires numerical computation (not analytically solvable)
- ☐ Has clear input data (e.g., eigenvalues, measurements)
- ☐ Has measurable output (e.g., parameter extraction, comparison)
- ☐ Has a natural 2–3 subproblem decomposition
- ☐ Reflects modern research (paper is recent, cutting-edge)

**Physics domain:** (e.g. Statistical Physics, Wave Physics, Quantum Mechanics, Electromagnetism, Thermal Physics, Materials Science, Fluid Dynamics, Astrophysics — or other)

---

## Task Architecture (3 Functions)

### Subproblem 1 (SP1)
- **Function name:**
- **Role:** (data processing / feature extraction / simulation-output processing / other)
- **What it does:**
- **Input parameters** (name, type, units, range):
- **Output** (dict keys, types, physical meaning):
- **Why SP1 comes first / paper reference:**

### Subproblem 2 (SP2)
- **Function name:**
- **Role:** (fitting/optimization / parameter extraction / analytical model / comparison / other)
- **What it does:**
- **Input parameters:**
- **Output:**
- **Why SP2 builds on SP1 / paper reference:**

### Main Problem (MP)
- **Function name:**
- **Role:** (synthesize SP1+SP2 / compare atomistic vs analytical / validate model / other)
- **What it does — and why it's NOT just "call SP1 then SP2":**
- **Input parameters:**
- **Output:**
- **Paper reference for the synthesis/comparison:**

---

## Difficulty Lever Design

See `references/difficulty-levers.md` in this skill for ready-to-adapt implementations (return-shape traps, n-vs-(n-2) denominators, boundary-inclusion traps, complex-phase traps, aggregate-vs-count traps, the zero-cost-correlation trap, None-handling, self-contained MP).

- **Which lever(s), and in which section(s)?**
- **Natural mistake a naive implementer would make:**
- **Why this lever fits this task's physics:**

---

## Target Pass Rate

- ☐ 2–3/8 (sweet spot) ☐ 1–5/8 (acceptable band) ☐ other, with reasoning:

---

## Pre-Design Checklist

- ☐ Problem requires numerical work (not analytically solvable)
- ☐ Problem is NOT just a recipe ("do this, then do that")
- ☐ Parameters introduce real physical context, not bare symbols
- ☐ Main problem adds real synthesis, not just chaining SP1+SP2
- ☐ SP1 and SP2 are each smaller in scope than MP
- ☐ Paper is recent and cutting-edge

## Key Equations

**Equation 1:** (what does it compute?)

**Equation 2:**

**Key concept the task should teach:**

## Example Input/Output Values (for writing tests later)

- Simple/edge case:
- Realistic case:
- Boundary case:
