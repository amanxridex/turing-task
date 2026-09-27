---
name: "scicode-task-builder"
description: >
  Build or fix SciCode-style physics benchmark tasks (Subproblem 1, Subproblem 2, Main Problem) for the CTP platform, delivered as Colab-ready code boxes. Trigger this whenever the user hands over an arXiv paper, PDF, or paper title and asks to build, write, or draft a task from it — even without the word "SciCode," e.g. "turn this paper into subproblems," "build a CTP task from this," "write the notebook sections for this arXiv paper," "here's the paper, make me a task." Also trigger when the user hands over an existing task plus reviewer or pod-lead feedback and wants it fixed or reworked, for designing difficulty levers or discriminative tests, for debugging a CTP validator error, or for interpreting or tuning a Gemini pass@k calibration result (e.g. "getting 8/8," "0/8," "reviewer says SP1 is too easy," "model keeps passing every test").
---

# SciCode Task Builder

Builds a complete SciCode-style benchmark task for the CTP platform — Subproblem 1, Subproblem 2, Main Problem — from a research paper, or fixes an existing task against reviewer feedback or a pass@k result. Either way the output is Colab-ready code boxes.

## Which path?

- **Have a paper, no task yet** → **Path B: New Task Creation**
- **Have an existing task (all 9 sections) plus reviewer comments or a pass@k number** → **Path A: Task Reworking**
- **Both at once** (e.g. reworking with the original paper on hand for grounding) → treat it as Path A; use the paper the same way Path B would, for grounding the fix.

Both paths share the pod-lead rules, notebook structure, test-writing standards, difficulty levers, validator fixes, and pass@k table below. Path A doesn't relax any of them — it just applies them to the specific sections that are flagged, instead of writing all nine from a blank page.

## Pod-lead rules (apply throughout, both paths, not just at the end)

These are non-negotiable constraints on the task itself, independent of the phase-by-phase mechanics below:

- No instructional or implementation-style prompts — the prompt states the physics problem, not a recipe for solving it
- Physical parameters must be introduced with real context, never as bare unexplained symbols
- No over-simplified toy models
- The problem must not be analytically solvable — it should require numerical work
- Test cases must be complete, cover the relevant edge cases, and never leave an undefined `ValueError` check dangling
- Subproblems must not be conceptually bigger than the main problem
- The main problem must not be *just* "call SP1 then call SP2" — it needs its own physical synthesis
- Every task is built from scratch and physically motivated — no shortcuts, no generic filler
- Target recent, cutting-edge papers

## Notebook structure

```
Metadata
Title
# Subproblem 1
  ## Prompt / ## Background / ## Testing Template / ## Solution
# Subproblem 2
  ## Prompt / ## Background / ## Testing Template / ## Solution
# Main Problem
  ## Prompt / ## Background / ## Testing Template / ## Solution
```

---

## Path B — New Task Creation

Work through all four phases below in order and hand back each notebook section in its own code box.

### Phase 1 — Task Architecture

**Step 1 — Identify the 3 functions**, one per section:

| Section | Role | Example |
|---|---|---|
| SP1 | Data processing / binning / transformation | Bin phonon modes into spectral density |
| SP2 | Fitting / optimization / parameter extraction | Fit analytical model via `least_squares` |
| MP | Combine SP1 + SP2 → physical observable | Compare atomistic vs analytical bandwidth |

**Step 2 — Design at least one difficulty lever per section — SP1, SP2, AND MP, all three, not just the two subproblems.** A single well-chosen feature where the natural, spec-compliant implementation is subtly wrong is often enough, but it's normal — and sometimes necessary to land in band — to layer two or three independent traps into one section (e.g. a data-processing subproblem stacking a boundary-inclusion trap, a complex-phase trap, and an aggregate-vs-count trap all at once). Reach for the proven discriminators in `references/difficulty-levers.md` first (return-shape traps, `n`-vs-`n-2`/chi-square-factor denominators, boundary-inclusion traps, complex-phase traps, aggregate-vs-count traps, scale-invariant edge cases, `None`-handling exceptions, self-containment for MP); invent a new paper-specific one only if none of those fit the physics naturally. MP's lever is usually self-containment (Lever E) — forcing the model to embed correct copies of SP1 and SP2 from memory — but add a numerical trap on top of that whenever the physics naturally supports one; self-containment alone doesn't always land MP in band on its own.

Before moving to Phase 2, write down a one-line estimated pass rate for each of the three sections individually. If any one of them looks like it'll land outside the acceptable band (below 1/8 or above 5/8), fix that section's lever now — don't let SP1 and SP2 look good while MP is left too easy or too hard.

**Step 3 — Target pass rates.** These bands are nested, not competing — aim for the tightest one and treat the wider ones as "still acceptable, don't touch it further":

| Band | Range | Meaning |
|---|---|---|
| Sweet spot | 2–3 / 8 | What to aim for per section |
| Acceptable | SP1 & SP2: 2–4/8, MP: 3–5/8 | Fine as-is, no further tightening needed |
| Outer still-in-band | 1–5 / 8 | Don't touch; only Phase 4's 0/8 or 6–8/8 rows call for a fix |

If the user hands over a paper without much architecture detail already worked out, design this phase yourself directly from the paper — don't stop to make them fill out a worksheet first. `assets/new-task-template.md` is available if they'd rather plan on paper before handing you the material, but it's optional scaffolding, not a prerequisite.

### Phase 2 — Writing the Notebook

#### Prompt

- LaTeX: `\( inline \)` and `\[ display \]`
- State units for EVERY parameter
- Full docstring: Inputs, Outputs, Raises
- The docstrings in the function templates of the prompts must EXACTLY match those in the corresponding solutions
- The prompt section must contain NO MORE THAN ONE Python code string/block
- Return stub must EXACTLY match the solution's return — copy-paste the return block between them
- Explicitly state `None`/null rejection, default values, edge cases

#### Background

- 1–2 paragraphs of prose — no bullets, no LaTeX, no headers
- Explain WHY each design choice matters physically
- Must read like a domain expert wrote it, not like an AI — don't just repeat the prompt, add complementary context

#### Tests

Tests are not a fixed-count checklist. Accepted tasks typically run **6–10 test cases per section — not some round number** — and the right count follows from what the physics actually needs to guard against, not from a template. Padding a suite out to a fixed target produces repetitive, formulaic tests that read as AI-written and rarely survive review.

**Illustrative shape, from an accepted task** (an illustration of proportions, not a layout to copy):

| Section | # tests | Rough breakdown |
|---|---|---|
| SP1 | 9 | 3 core-correctness · 1 type/count check · 2 boundary edge cases (each paired with its own discriminative negative assertion) · 1 zero-coupling sanity check · 1 bundled invalid-input test covering 4 scenarios |
| SP2 | 8 | 2 core-correctness (clean fit + noisy fit) · 1 type check · 2 discriminative checks · 1 covariance cross-check against finite differences · 1 bundled invalid-input test covering 3 scenarios |
| MP | 6 | 1 structural+range check · 2 exact-value regression tests on different data · 1 bundled parameter-validation test · 2 physically-motivated `ValueError` tests |

Cover these categories, folding several into one test function wherever they belong together:

- **Core correctness** — inputs simple enough to verify the expected value by hand or in a couple of lines of numpy.
- **Structural / convention checks** — return type, dict keys, tuple length, sign convention. Usually folded into a correctness test rather than given its own slot.
- **Boundary/edge behavior** — whatever the prompt calls out as a special case (an inclusive endpoint, a value exactly on a boundary, a zero-cost/degenerate limit).
- **Discriminative check(s)** — the lever(s) from Phase 1 Step 2. Usually an *extra assertion added to a correctness test* — compute the right answer, assert it matches, then assert it does NOT match the natural-mistake value — rather than a separately numbered test. One section can carry two or three independent discriminative assertions spread across different tests.
- **Invalid-input handling** — bundle every `ValueError`/`None`/out-of-range case the prompt specifies into a **single** test that loops over a list of bad-parameter dicts, or chains several `try/except (ValueError, TypeError): ... assert raised` blocks. Don't spend one `test_case_N` per invalid scenario.
- **(MP only) physically-motivated failure modes** — ways the combination can legitimately break even when SP1 and SP2 are each correct (a degenerate all-zero curve, a peak with no decaying branch). One test per failure mode.

Mechanical rules that always hold, whatever the count:

1. Every test starts with `def test_case_N():` as the first character (no leading blank lines)
2. Every test does `import numpy as np` inside the function body
3. The test's own `assert` statements sit at the top level of `test_case_N`, never inside a nested `def`
4. Every test references the main function by name at the top level
5. A small nested `def` that only computes an expected/reference value — then is called and asserted on at the top level — is fine and common; that's different from nesting the asserts themselves, and several accepted tasks reach for a helper `def` here rather than a `lambda`
6. Each test is fully self-contained: own data, own imports, no shared fixtures
7. `None`/invalid-input tests catch `except (ValueError, TypeError)` — never just `except ValueError`
8. MP tests never call the SP1/SP2 functions directly to derive an expected value. Instead, either hand-derive it inline, or — more common for MP, since its output is usually too composite to hand-derive — hardcode a high-precision value (e.g. `rtol=1e-7`) obtained by actually running the correct reference solution once
9. Discriminative tests check exact numerical values that differ between the correct implementation and the natural mistake
10. Every discriminative test carries a negative assertion too: `assert not np.isclose(out['x'], wrong_value)`

#### Solution

- `if X is None or not np.isfinite(X):` — the `None` check goes BEFORE `isfinite`
- Cast every output explicitly: `float(...)`, `int(...)`
- Solution Ordering: Unchanged SP1 function, followed by unchanged SP2 function, and then the Main Problem function that calls and uses both
- The MP solution embeds complete copies of the SP1 and SP2 functions
- Embedded copies must match the SP1/SP2 solutions exactly (copy-paste, don't retype)

### Phase 3 — Quality Control (QC1 & QC2 Validation)

When syncing the notebook to the CTP workbench, Tier 1 and Tier 2 validation are triggered automatically. Always monitor the validation sequence in the Activity tab on the right side.

#### QC1: Deterministic and Structural Checks (Non-bypassable)
Resolve every valid QC1 flag directly in the Colab notebook.

| QC1 Error | Root Cause & Fix |
|---|---|
| `Missing Assert Statements` in test_case_N | Asserts are trapped inside a nested `def` — move them back to `test_case_N`'s top level. |
| `Tests Do Not Reference Declared Function` | Move the main function call to the top level, out of any nested def. |
| `prompt return does not match solution` | Copy-paste the return block so both are identical, whitespace included. |
| `docstring does not match solution` | Copy-paste the docstring from the solution into the prompt template. |
| `Multiple code blocks in prompt` | Keep strictly ONE Python code string/block in the prompt section. |
| `None/null handling not covered` | Add a test with `except (ValueError, TypeError)` for None inputs. |
| `IndentationError line 1` | Remove ALL blank lines and comments before `def test_case_N():`. |
| `test block exec failed` | Make the test self-contained; remove any cross-block dependency. |

#### QC2: Evaluation & Defensible Flags
QC2 evaluates scientific correctness, test-case discriminativeness, and determinism.
- **Genuine Issue:** Fix in Colab, save, sync, and re-validate.
- **False Positive:** Provide a detailed, factual, non-vague defense in the platform's response field explaining why the task is scientifically correct and why the automated check does not apply.

### Phase 4 — Model Pass-Rate Evaluation Ladder (Gemini 3.1 Pro → GPT-5.6 → Claude Opus 4.8)

Automated pass-rate evaluation is run 8 times per problem component (SP1, SP2, MP).

#### Model Sequence & Escalation Rules
1. **Base Model:** Start on **Gemini 3.1 Pro** for all problems.
2. **Acceptance Threshold:** All components (SP1, SP2, MP) must score between **1/8 and 5/8 (inclusive) on the SAME model**.
3. **Immediate Rejection Threshold:** Any score of **6/8, 7/8, or 8/8 on ANY component** immediately rejects the task (too easy). Rework the difficulty levers, save, sync, pass QC, and rerun.
4. **Progression (0/8 rule):** A score of **0/8 on any component** (with no component $\ge 6/8$) triggers escalation to the next stronger model:
   $$\text{Gemini 3.1 Pro} \longrightarrow \text{GPT-5.6} \longrightarrow \text{Claude Opus 4.8}$$
5. **Task Not Submittable:** If Claude Opus 4.8 still yields 0/8 or fails to get all parts into 1–5/8 without exceeding 5/8, the task is discarded.

#### CTP UI Controls: Unit Re-run vs Header Continue
- **Per-Problem Unit Re-run (Lightning icon on SP1 / SP2 / MP):**
  - Use after changing that problem's code/tests and syncing Colab.
  - Restarts Gemini for that problem only (live k). This is your calibration loop.
- **Header Continue ("Continue GPT on all problems"):**
  - Use ONLY when Gemini has finished on every problem, is in progression (has 0/8 with none $\ge 6/8$), and the UI indicates escalation.
  - Runs the next model on ALL problems simultaneously.
  - *Do NOT use Continue* if you just synced a notebook fix, or if any problem scored $\ge 6/8$.
- **Note:** Syncing Colab alone does NOT clear old pass@k scores; clicking Re-run does.

#### Calibration & Debugging
- **Debugging 0/8:** Read model response transcript. Common blockers: `TypeError: ufunc 'isfinite' not supported` (widen to `except (ValueError, TypeError)`), `NameError` in MP (inline reference values instead of calling SP), or overly strict tuple/dict shape requirements.
- **Debugging 6–8/8 (Too Easy):** Add a subtle difficulty trap from `references/difficulty-levers.md` (e.g. $n$ vs $n-2$ degrees-of-freedom denominator, factor of 2 in $\chi^2_{\text{red}}$, boundary inclusive/exclusive check). Verify that reference solution passes while the natural mistake fails.

---

## Path A — Task Reworking

Use this when the user already has a task (all 9 sections: Prompt/Background/Tests/Solution × SP1, SP2, MP) and wants it fixed — because a reviewer flagged issues, a pass@k run came back out of band, or a validator error is blocking submission. It applies the exact same standards as Path B, targeted at whichever sections are actually broken, rather than rewriting all nine from scratch.

### Phase 0R — Parse & Diagnose

1. Get the current task's 9 sections and whatever feedback exists — reviewer comments, a pass@k number, or a validator error message. If any of this is missing or vague, ask for it rather than guessing at what's wrong.
2. Sort each piece of feedback by what it's actually about: physics grounding, a specific test-writing standard, the difficulty lever, or implementation/validator mechanics — a single comment can point at more than one of these.
3. Map every issue to the specific section (SP1/SP2/MP) and layer (Prompt/Background/Tests/Solution) it lives in. This map is what keeps Phase 2R targeted instead of a full rewrite.
4. Re-check the current task against the pod-lead rules above, independent of what the reviewer said — feedback can miss a violation that's still there.

### Phase 1R — Architecture Review

1. Confirm SP1/SP2/MP still have the right conceptual sizes relative to each other per the pod-lead rules — feedback sometimes reveals that a "subproblem" quietly grew bigger than the main problem, or that MP is just chaining SP1→SP2.
2. Re-derive by hand whether the difficulty lever(s) in each flagged section actually discriminate: does the natural mistake really fail the current discriminative test? A common rework issue is a lever that looks right on paper but was never actually verified against the tests.
3. If there's a pass@k number, use the Phase 4 table above to tell whether the issue is "too hard" (fix a blocker, don't touch the lever) or "too easy" (needs a new or stronger lever from `references/difficulty-levers.md`).

### Phase 2R — Notebook Rework

Rework only the sections Phase 0R flagged, applying the same Prompt/Background/Tests/Solution standards as Path B's Phase 2 — same mechanical test rules, same "reads like a domain expert, not an AI" bar for Background, same copy-paste-identical return stubs. Don't touch sections that weren't flagged unless a change elsewhere forces it (e.g. reworking SP1's return shape means MP's embedded copy of SP1 has to match it).

### Phase 3R — Validator

Same table as Path B's Phase 3 — a rework request is often triggered by exactly one of those rows.

### Phase 4R (optional) — Pass@k Recalibration

Same table and debugging notes as Path B's Phase 4.

### Handing it back

Hand back only the sections you actually changed, each in its own code box, plus a short note at the top: what was flagged, which pod-lead rule / test standard / lever it violated, and what you changed. List which of the 9 sections are untouched so the user knows what to leave alone. Self-check the same way as Path B before calling it done: reference solution passes all tests, the natural mistake fails the discriminators, and the estimated pass rate for every section you touched lands in band.

If the user pastes the current task and feedback directly into chat, work from that — don't make them fill out a worksheet first. `assets/rework-input-template.md` is there if they'd rather organize their notes before handing you the material.

---

## Reference material

`references/difficulty-levers.md` has the proven, ready-to-adapt discriminator implementations (return-shape traps, `n`-vs-`n-2`/chi-square-factor traps, boundary-inclusion traps, complex-phase traps, aggregate-vs-count traps, the zero-cost correlation trap, `None`-handling test pattern, the self-contained-MP-test pattern, and the golden-regression-value pattern for MP) with working code. Used by both paths — read it while designing Phase 1 Step 2, while debugging an 8/8 result in Phase 4, or while picking a stronger lever in Phase 1R.

`assets/new-task-template.md` and `assets/rework-input-template.md` are optional fill-in worksheets for the user's own planning, one per path — point to the relevant one only if the user wants structure before handing over the material.

## Output format

Always hand back each notebook section in its own separate code box, ready to paste directly into the Colab notebook — never as one long undifferentiated block. Before handing anything over, self-check separately for every section you're delivering: reference solution passes all tests, the natural/naive implementation fails the discriminative tests, and the estimated pass rate lands in band per the Phase 1 / Phase 1R table above. Don't declare a task — or a rework — done until every delivered section clears this check.
