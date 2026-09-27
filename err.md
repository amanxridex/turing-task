# Turing CTP SciCode Skill Compliance & Gap Analysis (`err.md`)

This document compares the skill package (`scicode-task-builder`) against the **Turing Central Tasking Platform (CTP) SOP** and the **12-Point Workflow Checklist** for the `Nvidia_STEM_SciCode_3k` project.

---

## Executive Verdict

| Area | Match Level | Summary |
|---|---|---|
| **Scientific & Problem Architecture** | **High (95%)** | Aligns with subproblem decomposition (SP1, SP2, MP), numerical requirements, domain-expert context, and avoiding toy models. |
| **Code, Prompts & LaTeX Standards** | **High (90%)** | Enforces LaTeX formatting, units, docstrings, exact return stubs, and copy-paste return matching. |
| **Test Suite & Difficulty Levers** | **Very High (95%)** | Follows the 6–10 test structure, negative assertions (`assert not np.isclose`), bundled invalid tests, and rich traps in `difficulty-levers.md`. |
| **Model Pass@k Calibration Ladder** | **Partial (40%)** | **Critical Gap**: The skill only documents **Gemini 3.1 Pro**, omitting the multi-model cascade to **GPT-5.6** and **Claude Opus 4.8**, as well as the platform UI rules for "Unit Re-run" vs "Header Continue". |
| **CTP Platform & QC Mechanics** | **Low (30%)** | **Gap**: Does not cover CTP queue rules, metadata fields, "no more than 1 code block in prompt", Colab Sync mechanics, QC1 vs QC2 defense writing, or L1/L2 delivery workflows. |

---

## Part 1: What is PRESENT in the Skill

The skill files (`SKILL.md`, `difficulty-levers.md`, `new-task-template.md`, `rework-input-template.md`) cover the core problem authoring and testing requirements:

### 1. 12-Point Checklist Alignments
- **Checklist Item 1 (Paper selection & decomposition):** Enforces selecting cutting-edge papers, avoiding analytically solvable problems, and decomposing into SP1 (data processing), SP2 (fitting/modeling), and MP (physical synthesis).
- **Checklist Item 2 (Scientific accuracy & consistency):** Grounding in physics research, realistic parameters, and avoiding over-simplified toy models.
- **Checklist Item 3 & 4 (Well-posed, self-contained prompts):** Prompts must state physical context, specify units for every parameter, define mathematical symbols, and avoid instructional recipe-style hints.
- **Checklist Item 5 (LaTeX & phrasing):** Strict rules for LaTeX (`\( inline \)` and `\[ display \]`) and domain-expert phrasing for the Background section (no bullet points, no LLM clichés).
- **Checklist Item 6 (Comprehensive test coverage):** 6–10 tests per section covering core correctness, boundaries, bundled `ValueError`/`TypeError` handling, and discriminative negative assertions.
- **Checklist Item 7 & 11 (Assert-based tests):** Tests are purely assert-based, starting with `def test_case_N():` with top-level `assert` statements.
- **Checklist Item 9 (Main Problem structure):** Explicitly requires MP solution to embed identical copy-pasted versions of SP1 and SP2, and call both to synthesize the final result.
- **Checklist Item 10 (Return matching):** Enforces exact copy-paste matching of return statements between function templates and solutions.
- **Checklist Item 12 (Direct function validation):** Mandatory mechanical rule that every test calls the corresponding declared function at top level.

### 2. Difficulty Levers & Calibrated Traps (`references/difficulty-levers.md`)
- Multi-part return with mixed types (`float` vs `int`).
- $n$ vs $n-2$ denominator and missing factor-of-2 degrees-of-freedom traps in $\chi^2_{\text{red}}$.
- Boundary inclusion/exclusion traps ($<$ vs $\le$).
- Complex-phase conventions and quadrant traps.
- Self-contained regression tests with golden reference values.

### 3. Basic Validator Handling
- Table of common Tier 1 & Tier 2 validator syntax errors: `Missing Assert Statements`, `Tests Do Not Reference Declared Function`, `prompt return does not match solution`, `None/null handling not covered`, `IndentationError line 1`.

---

## Part 2: What is NOT PRESENT / Missing in the Skill

The following requirements from the CTP SOP and workflow are **missing or incomplete** in the skill:

### 1. Multi-Model Evaluation Ladder (CRITICAL)
- **Missing Models:** The skill assumes evaluation is solely on **Gemini 3.1 Pro**. The CTP SOP mandates a sequential 3-model cascade:
  $$\text{Gemini 3.1 Pro} \longrightarrow \text{GPT-5.6} \longrightarrow \text{Claude Opus 4.8}$$
- **Progression Logic:**
  - A score of `0/8` on any component (with no scores $\ge 6/8$) does **not** always mean the task is broken; in CTP, it triggers progression to the next stronger model in the sequence.
  - Immediate Rejection: Any score $\ge 6/8$ (`6/8`, `7/8`, `8/8`) on **any** component on the current model triggers instant rejection.
  - Acceptance Rule: All components (SP1, SP2, MP) must simultaneously land in **1/8 to 5/8 on the same model**.
  - Exhaustion Rule: If the task still has a `0/8` on Claude Opus 4.8, it is un-submittable and discarded.

### 2. CTP UI Controls: "Unit Re-run" vs "Header Continue"
- **Unit Re-run (Lightning icon on SP1 / SP2 / MP):** Must be used when code/tests for that specific problem are edited and synced from Colab. Restarts Gemini only for that problem.
- **Header Continue ("Continue [Next Model] on all problems"):** Only used when Gemini is complete across all problems and the UI prompts escalation. Runs the next model on **all** problems simultaneously.
- **Warning Rule:** Saving/syncing Colab does **not** clear old pass@k scores; only running/re-running does.

### 3. QC1 vs QC2 Rubrics & Formal Defense
- The skill groups all validator errors into a simple bug-fix table. It lacks:
  - **QC1 (Deterministic/Structural Checks):** Non-bypassable hard checks (formatting, required fields, code placement, reference links).
  - **QC2 (Rubric Evaluation & Defensible Flags):** Automated checks for discriminativeness and correctness.
  - **Defense Writing:** Missing instructions on how to write a factual, non-vague defense in the platform's response field when QC2 produces a false-positive flag.
  - **QC Wait Times:** SOP notes QC takes 5–30 minutes, during which contributors should claim another task.

### 4. Notebook Constraints & Metadata Requirements
- **Subproblem Count:** The skill hardcodes a 3-part layout (SP1, SP2, MP). The CTP SOP specifies that while 2 subproblems is standard, tasks **may use 3 or 4 subproblems** (max 4) when justified by decomposition.
- **Code Block Limit in Prompt:** CTP strictly mandates that the prompt section must contain **no more than one Python code string/block**.
- **Metadata Fields:** Missing CTP-specific fields: `Nvidia_STEM_SciCode_3k`, Batch domain (e.g., `chem_started_01_sep`), Paper Reference Link, and Subdomain tagging.
- **Full Docstring Matching:** Checklist Item 8 requires the docstring in the prompt function template to **exactly match** the docstring in the solution; the skill currently only enforces matching the return statement.

### 5. Task Queues, Review & Delivery (L1/L2)
- Queue navigation in CTP (`Batches` $\rightarrow$ claim domain $\rightarrow$ `My Tasks`).
- L1 reviewer evaluation and resolution of reviewer rubric feedback.
- Escalation rule: contacting the domain lead if a task remains pending review longer than 1 week.

---

## Actionable Recommendations to Upgrade the Skill

1. **Update `SKILL.md` Phase 4:** Add the complete model cascade table (Gemini 3.1 Pro $\rightarrow$ GPT-5.6 $\rightarrow$ Claude Opus 4.8) and the UI execution rules (Unit Re-run vs Header Continue).
2. **Add QC1/QC2 Defense Guidelines:** Include a section on writing valid, evidence-based rebuttal defenses for false-positive QC2 flags.
3. **Add CTP Colab Formatting Constraints:** Enforce the "max 1 python code string per prompt" rule, metadata block requirements, and exact docstring matching between prompt template and solution.
