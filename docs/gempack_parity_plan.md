# GEMPACK TABLO Feature-Parity Plan

Status: DRAFT (2026-07-24). Companion to ROADMAP Phase 6 ("GEMPACK parity
program — performance-vetoed"). Refines/expands ROADMAP item 6.4 ("Statement
completion … full IF semantics, POSTSIM, others"). Not yet committed.

Spec source: `teems-docs/solver_optimization/GEMPACK User Manual.html`
(Release 12.2+). Ground truth: solver parser (`teems-solver/src/*.c`) + R
front-end (`teems-R/R/tablo_process.R`, `chk_statements.R`,
`data-raw/teems_data-raw.R`). Test corpus (per ROADMAP): `original_tabs/`,
`teems_archive/tab_files/` (GDYNv3.6.tab, gtapep.tab) — NOT the reduced TABs
in `teems-R/inst/models/`.

Standing rules (ROADMAP): full parity is the goal; **performance is the only
veto**; every feature must be **zero-cost when unused**; the manual is the
spec; model semantics live solver-side (single-source-of-truth), the R side
only pre-flight-validates and may do declarative rewrites (as it already does
for IF and netcut).

## Legend
- Status: ✅ full · ◐ partial/inert · ✗ absent
- Call: **PARITY** (implement) · **FIX** (bug/inert) · **QUESTION** (scrutinize scope/obsolescence) · **SKIP** (obsolete — parse-and-ignore at most)
- 🔍 = needs an explicit design/scope decision before work starts

---

## Tier 0 — POSTSIM  (the headline addition)  — Call: PARITY

**Current:** ✗ solver (no stage concept at all); R ignores it with a warning
(`ignored_state` in `teems_data-raw.R:41`; `tablo_process.R:165` literally
comments "potentially handle postsim here"). ROADMAP 6.4 parks it, unspecced.

**What it is (manual ch.12, §10.18):** a section of the model's own TAB file
whose statements run **once, after the simulation solves**, to compute and
write result-dependent numbers (decompositions, summary tables, welfare, GDP
aggregates). It replaces the old multi-program chore (SLTOHT + a separate
report TAB reading initial data + updated data + SOL file, e.g.
DECOMP.TAB/GTAPVOL.TAB).

**Syntax:**
- Section: `POSTSIM (BEGIN); … POSTSIM (END);` (multiple sections allowed;
  all conceptually concatenated at end of file).
- Single-statement qualifier: `(POSTSIM)` on a `Write` or `Display`.
- CMF side: `XPostSim (Begin)…(End)` with `XSet/XSubset/XFile/XWrite/XDisplay`;
  run-time off switch `PostSim = no;`.

**Semantics (the subtle part — this is what makes it more than Formula+Write):**
- Inside a PostSim section you may access **simulation results as if the
  Variables were Coefficients** (their %-change / change values).
- An **ordinary coefficient on a RHS carries its POST-simulation (updated)
  value**, not its pre-sim value. To use a pre-sim value you must have
  captured it earlier with `Formula (Initial) SAVEORIG = SAVE;` into a
  separate (parameter) coefficient. → depends on real `(parameter)`/`Initial`
  semantics (Tier 1/3).
- LHS of a PostSim Formula may only be a **PostSim Coefficient** — never a
  Variable, never an ordinary Coefficient.
- **Scope isolation:** PostSim sets/subsets/mappings/coefficients cannot be
  referenced from the ordinary part.
- INITIAL/ALWAYS and PARAMETER/NON_PARAMETER qualifiers are **ignored** in
  PostSim (every PostSim formula runs exactly once); `Default` statements are
  disallowed there. ZeroDivide state **resets** at the start of the PostSim
  pass.
- Allowed inside: **Set, Subset, Coefficient, File, Read, Formula, ZeroDivide,
  Mapping, Write, Display, Assertion.** Forbidden inside: Variable, Equation,
  Update, Transfer, Omit/Substitute/Backsolve, Complementarity.
- Execution (manual §12.5): (1) re-run the ordinary Reads+Formulas so normal
  coefficients hold **updated** values; (2) organise PostSim sets/subsets/
  mappings; (3) run PostSim Reads/Formulas/Assertions/Writes/Displays in file
  order. PostSim sets go to the Solution file but not the Equations file.

**Prerequisites (must land before or with PostSim):**
1. **General WRITE** (Tier 2) — PostSim's whole point is writing result tables
   (HAR + text + terminal). Today's solver Write is a narrow CMF-tied form.
2. **Post-sim value exposure** — the solver already computes UPDATE results, so
   updated coefficient values exist post-solve; PostSim needs them surfaced,
   plus **variables-as-coefficients** (the %-change/change results readable by
   formula). This is the main new solver plumbing.
3. **`(parameter)` semantics + `Formula (Initial)` capture** (Tier 1/3) — the
   pre-sim SAVEORIG idiom depends on it.
4. **Intrinsic functions** (Tier 3) — decomposition/report formulas routinely
   use EXP/LOG/ID01/ID0V.
5. **Minimal DISPLAY** (Tier 2) if `Display (postsim)` support is wanted.

**Split:** engine is solver-side. R must stop treating Postsim as
`ignored_state`, pass the section through, and surface its written headers via
`ems_compose`/`ems_deploy`. 🔍 decision: HAR vs text output routing through the
container.

**🔍 Scrutiny:** the 3-pass "re-run ordinary reads+formulas to rebuild updated
coefficients" model is the load-bearing design choice — confirm TEEMS can
replay the update state cheaply (it already retains updated levels for
multi-step). Scope the first cut to the common case (Write summary tables from
post-sim coefficients + variable results); defer PostSim Mapping and PostSim
Read-into-coefficient to a second pass.

---

## Tier 1 — Correctness gaps in already-"supported" features (do first; cheap, high-value)

| # | Feature | Status | Manual role | Call | Notes |
|---|---|---|---|---|---|
|1.1|**ASSERTION evaluation**|◐ solver recognises the keyword in preprocessing but **never evaluates it** (no `tab_next_statement("assertion")` consumer); R passes it through|CORE, "use many assertions"|**FIX/PARITY**|Silent-confidence bug: models with assertions run but nothing is checked. Implement condition evaluation (reuse the formula/condition engine) with abort-on-fail + failing-element report; honour CMF `Assertions = no/warn`. High value, self-contained.|
|1.2|**`(parameter)` / `(non_parameter)` semantics**|◐ `(parameter)`/`(integer)` stripped, no distinct type; `(non_parameter)` unmatched|Parameters are constant through a sim|**FIX + QUESTION**|Two problems: (a) the parser strips the substring `"parameter"`, which **corrupts `non_parameter` → `non_`** (real bug — fix regardless). (b) Whether the missing parameter/non-parameter distinction actually changes solve results is 🔍 — it does matter for **PostSim** (post-sim value rules) and for the SAVEORIG idiom.|
|1.3|**Conditional quantifier `(all,i,SET:cond)` in EQUATIONs**|✗ solver|Manual §11: allowed **only** in FORMULA(ALWAYS) and UPDATE — **forbidden in equations**|**NO ACTION (spec-correct)**|Resolves the solver-agent's flagged "gap": TEEMS is already parity-correct. Document it so it isn't re-opened.|
|1.4|**Variable `ORIG_LEVEL=`**|◐ tolerated, skipped, never tracked|Reports levels values for a linear variable|**QUESTION**|Output-fidelity nicety. 🔍 assess whether target models rely on orig-level reporting; if so, track and surface via compose.|

---

## Tier 2 — Missing statements (parity candidates)

| # | Statement | Status | Manual role | Call | Notes / scrutiny |
|---|---|---|---|---|---|
|2.1|**WRITE (general)**|◐ solver handles only `Write X to file <log> header "…"` bound to a CMF outdata logname; R **drops** Write entirely (`tablo_process.R:162`)|CORE — HAR/text/terminal, `(SET)`/`(ALLSETS)`, LONGNAME, BY_ELEMENTS, `(POSTSIM)`|**PARITY**|Prerequisite for PostSim. Implement HAR + text + terminal targets and the set-writing forms. R must stop dropping it.|
|2.2|**FORMULA & EQUATION (combined)**|✗|COMMON shorthand = FORMULA(INITIAL)+EQUATION(LEVELS)|**PARITY (cheap)**|Pure sugar over two already-supported statements. 🔍 could be an **R-side declarative rewrite** (like the existing IF/netcut rewrites), keeping the solver unchanged — preferred if it fits single-source-of-truth.|
|2.3|**DISPLAY**|✗ (R hard-aborts)|Manual calls it "old-fashioned"; SLC/CVL preferred, but used with PostSim|**QUESTION**|Implement the **`Display (postsim)`** path minimally for PostSim; treat standalone Display as low-priority legacy. 🔍|
|2.4|**MAPPING … FROM … TO**|✗ (R hard-aborts)|COMMON/ADVANCED; used in index expressions; `(ONTO)`/`(PROJECT)`|**QUESTION**|🔍 Assess prevalence in the target corpus (GTAP standard uses little/no mapping). Real cost (index-expression evaluation). Defer unless a corpus TAB needs it.|
|2.5|**COMPLEMENTARITY / MCP**|✗ (R hard-aborts)|ADVANCED (ch.51); bounds, NO_SPLIT internals, condensation/subtotal interplay|**QUESTION (major)**|🔍 Highest-cost item. Standard GTAP models don't use MCP. Decide scope: likely **defer** until a target model requires it; if adopted, it's a project of its own (bounds engine + NO_SPLIT variables + Newton path).|
|2.6|**TRANSFER**|✗|ADVANCED; copy HAR headers between files|**QUESTION → likely SKIP**|🔍 Pure data-plumbing. TEEMS's R orchestration + container data flow likely already covers header movement; a `Transfer` is arguably obsoleted by the R I/O layer. Confirm no corpus TAB depends on it, then parse-and-reject with a clear message.|
|2.7|**LOOP / BREAK / CYCLE**|✗ (R hard-aborts)|ADVANCED (Release 12+); mandatory `(all,…)` driver|**QUESTION → defer**|🔍 Niche control flow; rare in standard models; per-iteration overhead. Defer decision pending a real use case.|

---

## Tier 3 — Qualifier / function / expression completeness

| # | Feature | Status | Manual role | Call | Notes |
|---|---|---|---|---|---|
|3.1|**Intrinsic functions** (EXP, LOGE, SQRT, ABS, **ID01/ID0V**, ROUND/TRUNC, MAX/MIN, …)|✗ formula compiler has only LOAD/±/×/÷/POW/IF_*|CORE expression functions|**PARITY (important)**|🔍 Real models use EXP/LOG in CES/CET nests; **ID01/ID0V** are the manual's recommended replacement for Zerodivide (and work in equations). Pick the subset first: EXP, LOGE, SQRT, ABS, ID01, ID0V, ROUND. Prerequisite-ish for PostSim reports.|
|3.2|**ZERODIVIDE completeness** (`ZERO_BY_ZERO`/`NONZERO_BY_ZERO`, `OFF`)|◐ default form only|COMMON but manual calls it error-prone, recommends ID01/ID0V|**QUESTION**|Prefer investing in ID01/ID0V (3.1) as the modern path; implement Zerodivide's stateful qualifiers only for legacy-TAB parity. 🔍|
|3.3|**Set `MAXIMUM SIZE`**|✗|**OBSOLETE** — modern GEMPACK silently ignores it (§11.7.2)|**SKIP (parse-and-ignore)**|Accept the keyword and discard it so legacy TABs parse; implement no semantics. Trivial, removes a hard failure.|
|3.4|**UPDATE `(product)` qualifier**|◐ product is TEEMS's internal default; keyword not parsed|PRODUCT is the default update qualifier|**PARITY (trivial)**|Accept `(product)` as an explicit no-op default so TABs that write it parse.|
|3.5|**`Keyword (DEFAULT=…)` statements**|◐ `default=levels/change` handled; other forms partial|COMMON|**PARITY**|🔍 Ensure all forms parse: COEFFICIENT/VARIABLE/FORMULA/EQUATION default resets, and default LOWER/UPPER_BOUND. Low cost, prevents surprise aborts.|
|3.6|**Variable `LINEAR_NAME=` / `LINEAR_VAR=`**|✗|Override the auto `p_`/`c_`; associate with an existing linear var|**PARITY (LINEAR_NAME) / QUESTION (LINEAR_VAR)**|LINEAR_NAME is cheap naming. LINEAR_VAR (share one linear var across levels vars) is 🔍 — assess need.|
|3.7|**Variable `NO_SPLIT` / `VPQType`**|✗|NO_SPLIT internal (MCP/Newton); VPQType for homogeneity testing (ch.57)|**QUESTION → defer**|🔍 NO_SPLIT couples to MCP (2.5). VPQType is a testing-only annotation. Both defer.|
|3.8|**FILE output qualifiers** (TEXT/GAMS, ROW/COL/SPREADSHEET/SSE, SEPARATOR, FOR_UPDATES)|◐ binding is CMF-driven, TAB File qualifiers mostly ignored|Mixed core/niche|**QUESTION**|🔍 Implement **text/SSE output** only if PostSim Write-to-text needs it; GAMS export is niche → SKIP. FOR_UPDATES relates to updated-data files.|
|3.9|**READ `IfHeaderExists` / `BY_ELEMENTS`**|◐|Conditional read / set-mapping element read|**PARITY (minor)**|IfHeaderExists is a small robustness feature; BY_ELEMENTS pairs with Mapping (2.4).|
|3.10|**Subset `(BY_NUMBERS)`**|✗|Rare/obsolete (positional)|**SKIP**|Element-name subsets (default) are the modern form; don't implement.|

---

## Tier 4 — Simulation-stage (closure / shocks) — already at parity

Closure (`exogenous/endogenous/rest`, swaps) and shocks (scalar, component,
`uniform`, `select from` HAR, file) are ✅ via the closure/shock files
(`closure_read`, `shocks_read`). **SKIP** the legacy forms the manual itself
deprecates: component-number closure/shock addressing, text-file shocks,
interactive/STI closure.

---

## Obsolete / cumbersome — recommend QUESTION→SKIP (parse-and-ignore or reject cleanly)

Per the "question obsoleted/cumbersome features" directive, these are called
out as **not for reimplementation** (the manual itself deprecates them):

- Set `SIZE <int>` / `MAXIMUM SIZE + SIZE(integer)` **unnamed** sets (§10.1) —
  discouraged for new models.
- `MAXIMUM SIZE` on read sets (§11.7.2) — silently ignored by GEMPACK → 3.3.
- `SUBSET (BY_NUMBERS)` (§10.2) — rare/obsolete → 3.10.
- **STI-file condensation** and interactive TABLO `o/s/b/e` dialogue — TEEMS's
  in-TAB Omit/Substitute/Backsolve (ROADMAP 6.2, done) already supersedes it.
- Closure/shocks by **component number**, **interactive/STI** closure,
  **text-file shocks** — superseded by CMF + element names + HAR (Tier 4).
- TABLO **stop/restart stage** files (F1/F2/F3, `.TBR`/`.TBT`, ch.64) —
  explicitly deprecated in the manual.
- **UDC files** (§12.4.3) — made redundant by PostSim/SLC.
- **GAMS text export** `FILE(NEW,GAMS)` and hard-wired actual filenames — niche.
- `UNION` (vs disjoint `+`) and `-` (vs relative `\`) — support for legacy
  parsing (already done), but prefer the modern operators in any TEEMS-authored
  TABs.

---

## Scrutiny register (decisions needed before work)

1. **PostSim value-replay design** — confirm cheap reconstruction of updated
   coefficients for the post-sim pass (Tier 0).
2. **PostSim output routing** — HAR vs text through the container; how
   `ems_compose` surfaces PostSim headers.
3. **Parameter/non-parameter** — does the missing distinction affect any solve
   result today, or only PostSim? (1.2)
4. **Intrinsic-function subset** — first tranche (EXP/LOGE/SQRT/ABS/ID01/ID0V/
   ROUND?) (3.1).
5. **MCP scope** — implement or defer? Standard GTAP doesn't need it (2.5).
6. **Mapping prevalence** in the target corpus (2.4).
7. **Transfer** — obsoleted by the R/container I/O layer? (2.6).
8. **FORMULA&EQUATION** — R-side rewrite vs solver support (2.2).
9. **Text/SSE file output** — needed for PostSim? (3.8).

---

## Recommended sequencing

1. **Quick correctness wins (independent, cheap):** ASSERTION evaluation
   (1.1), the `non_parameter` strip bug (1.2a), `MAXIMUM SIZE` ignore (3.3),
   `(product)` no-op (3.4), Default-statement forms (3.5). Each zero-cost when
   unused; each removes a silent-wrong or a spurious abort.
2. **PostSim foundation:** general WRITE (2.1) + intrinsic functions (3.1) +
   `(parameter)` semantics (1.2b) + post-sim value / variable-as-coefficient
   exposure.
3. **POSTSIM engine** (Tier 0), incl. minimal `Display (postsim)` (2.3) and R
   passthrough.
4. **Legacy-parity qualifiers** as the corpus demands: FORMULA&EQUATION (2.2),
   LINEAR_NAME (3.6), Read IfHeaderExists (3.9), Zerodivide qualifiers (3.2).
5. **Scrutiny-gated majors** — decide per target-model need: Mapping (2.4),
   Complementarity (2.5), Loop/Break/Cycle (2.7), Transfer (2.6).

Validate every step against `original_tabs/` + `teems_archive/tab_files/`
(GDYNv3.6, gtapep) under the standing bit-identity + ASan gates; each feature
zero-cost when unused.
