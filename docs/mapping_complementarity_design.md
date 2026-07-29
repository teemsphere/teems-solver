# Mapping + Complementarity implementation design

Status: DESIGN (2026-07-28). Decision record: the step-5 corpus survey
(parity plan progress log, 74f9f7d) showed Mapping load-bearing in
GTAP-E/-EP/-AEZ and Complementarity in GMig2; user decision 2026-07-28:
**GTAP-E/-EP/-AEZ and GMig2 are targets — implement both statements.**
Transfer and Loop/Break/Cycle remain reject-cleanly. Spec sources:
GEMPACK manual 11.9 (Mappings), 10.17/11.14 (Complementarity
statement/semantics), ch. 51 (solution algorithm) — extracted from
`teems-docs/solver_optimization/GEMPACK User Manual.html`.

## 1. Mapping (manual 11.9)

A mapping `MAPPING [(ONTO)] M from S1 to S2;` is a total function from
source-set elements to target-set elements (many-to-one). Storage is an
integer vector over S1 with values = element positions in S2 (exactly
an integer coefficient in shape; manual 11.9.10 writes it as such).

**Corpus-required surface** (everything the target models use):

| Piece | Form in targets | Manual |
|---|---|---|
| Declaration | `Mapping MCOM from DCOMM to COMM;` (no `(onto)` uses in corpus; support the qualifier + its check anyway — cheap and "strongly recommended" by the manual) | 11.9 |
| Population | `Read (by_elements) M from file F header "H";` — CHARACTER data: for each s1, an element NAME of S2 (≤12 chars); converted to positions at read; range/onto checked before first use | 11.9.1(a), 11.9.2, 11.9.3 |
| Index use | `COEF(M(i))`, `var(M(i))` in formulas AND equations (`NCTAXB(REGTOBLOC(r))`, `del_nctaxb(REGTOBLOC(r))`) | 11.9.4, 11.9.5 |
| Condition use | `sum{r,REG: M(r) = b, expr}` in formulas and equations; also `sum(...)` paren form (GMig2) and `M(cc) EQ c` spelling | 11.9, 11.4.11 |

**Deferred (clean fatal, named message + manual section)** — zero
corpus uses: integer Read (11.9.1c), formula-assigned mappings
(11.9.1 b/d/e, incl. `$POS`-sourced — subsumed by R-side aggregation),
partial reads (11.9.11), LHS mappings (11.9.8), composition
`M1(M2(i))` (11.9.6), projection mappings (11.9.1f), index offsets on
mapping args, Write of mappings (11.9.10).

**Validation semantics** (11.9.2/11.9.3): values checked before first
use / end of preliminary pass, even for unused mappings; unknown
target element name in by_elements data = fatal; `(onto)` failure =
fatal; a mapping is immutable once used. Argument-set typing (11.9.7):
`COEF(M(i))` legal iff M's codomain equals (or is a subset of) COEF's
argument set at that position — enforce at compile.

**teems-R side**: `Mapping` is an unknown keyword to the R pipeline
today (folds into the preceding statement → `stmt_missing_equals`
abort). R must parse the declaration (name, from-set, to-set, onto
flag), accept `Read (by_elements)` targeting a declared mapping
(read-target check currently requires a coefficient), carry both
through `.finalize_tab()` unchanged, and validate: sets declared
(S-row machinery), name uniqueness (mapping names join the 11.2.1
namespace), by_elements read has a header. Conditional sums pass
through the R condensation/IF-rewrite untouched (verify — the `:`
condition inside `sum{}` must not trigger the binary-switch or IF
paths).

## 2. Complementarity (manual 10.17/11.14, ch. 51)

Statement: `Complementarity (Variable = X, Lower_bound = L[,
Upper_bound = U]) name (quantifiers) Expression;` — X a LEVELS
variable, L/U each a levels variable, `Coefficient(parameter)`, or
real constant; name ≤ 10 chars; quantifier sets must match X's (and
any variable/parameter bound's) argument sets position-wise (equal or
same-ordered subset). Semantics: X ≥ L, Expression ≥ 0, and
complementary slackness (X > L ⟹ Expr = 0; Expr > 0 ⟹ X = L);
mirrored for upper bounds. Three states (51.2): 1 = X at lower bound,
2 = Expression zero (X interior), 3 = X at upper bound.

**GEMPACK's algorithm** (51.1, 51.7):
1. **Approximate run**: many-step Euler with a derived non-smooth
   equation active. Derived statements per complementarity (11.14,
   51.7.2): levels change variable `comp@E` = Expression (via a new
   formula+equation); linear dummy `$comp@D` (switches the derived
   equation on/off); no-split `$del_Comp`; the linear equation
   `E_$comp` with 3 IF branches keyed on the CURRENT state, each
   branch holding its state's invariant through the step, with a
   Newton-correction term pulling bound overshoots back (51.1.2).
   Steps may be redone on state flips (51.7.3).
2. **Accurate run** (51.7.1): from the approximate run's final
   states, automatically modify closure+shocks — state 2: `comp@E`
   exogenized, shocked to zero (shock = −presim value); state 1/3
   with known bound (constant/parameter/exogenous levels var): X
   exogenized, shocked to the bound; state 1/3 with endogenous
   levels-variable bound: the auto-introduced difference variable
   (X − bound) exogenized, shocked to zero. ALL `$comp@D` components
   endogenized (turns the non-smooth equation off). Solve with the
   normal extrapolated method; verify post-sim states match
   (51.5.4/51.7.5); with subintervals + automatic accuracy, redo a
   subinterval on mismatch (51.7.4).

**Prerequisite discovered in the target model**: GMig2's
complementarity sits in a levels module — `Variable (change,levels)`,
`Equation (levels)` (E_RMIGSP_C), and `Formula&Equation` statements.
TEEMS rejects `Equation (levels)` and `Formula&Equation` by design
today. So Complementarity phasing must include the MINIMAL levels
subset: linearize `Equation (levels) LHS = RHS;` by automatic
differentiation into the existing linear pipeline (GMig2's levels
equations are simple sum/product forms), which also unblocks
`Formula&Equation` (its expansion needs exactly this). Full levels
parity stays out of scope.

**Closure/count semantics** (11.14): for closure purposes the
complementarity counts as one equation of its quantifier size.
Condensation guards (11.14.1): the complementarity variable must not
be condensed/backsolved/omitted; bound variables must not be
substituted/omitted (backsolve allowed) — enforce in the condensation
validators (R side + solver backsolve checks).

## 3. Phasing (each phase gated: verify.sh bit-identical, warnings
baseline, kit legs; feature zero-cost when unused)

- **M1** — Mapping declaration parse + storage + `Read (by_elements)`
  + validation (range/onto/immutability) + deferred-form fatals.
- **M2** — mapping-valued index compile in formulas and equations
  (incl. variable columns / Jacobian).
- **M3** — mapping-equality conditions in sums (formula + equation),
  `EQ`/`=` spellings, brace and paren sum forms.
- **M4** — teems-R passthrough + validation + e2e (GTAP-E-shaped
  fixture through deploy).
- **C0** — minimal levels-equation linearization (auto-diff of
  `Equation (levels)`) + `Formula&Equation` expansion riding on it.
  Full design in section 5 below (2026-07-29).
- **C1** — Complementarity parse + validation (qualifier/set
  matching, bound typing, name length, condensation guards) + derived
  statements (comp@E formula+equation, dummy, difference vars).
- **C2** — approximate-run state machinery (per-step state
  evaluation, 3-branch equation, Newton correction, step redo).
- **C3** — accurate-run automation (closure/shock modification,
  state-match verification, subinterval redo) + reporting
  (51.5.3-style state-change log lines).

Integration points (solver internals): see the survey appendix below.

## 4. Solver integration survey (appendix, 2026-07-28)

Key findings from the internals survey (file:line refs at survey
time, tree @ 74f9f7d):

**Index machinery.** Every operand index is `arSet[j].indx`
optionally routed through the static `superset_pos` table via
`dim_addr{ADims,SupSet,SSIndx,leadlag}` (teems_solver.h:361). The
per-tuple index loop is duplicated ~29x across formula.c:731-1224
plus jacobian.c:300-308 and jacobian.c:2992-2999 —
**prerequisite refactor: factor into one inline helper before adding
a mapping mode.** A mapping is then a fifth dim_addr mode (`MapId`):
`superset_pos`-like per-element int lookup, but through a
mapping-owned `int*` table (dedicated array, NOT a float coefficient
— `store_real` is float, exactness cap 2^24).

**Parsing hazard.** All index-list tokenizers (`strtok` on
`{`,`,`,`}` — formula_bind_operand formula.c:45+, eq_linvar_read
jacobian.c:455-565) shred nested braces, so `coef{c,map{r}}` cannot
survive normalize. **Design: pre-lower `MAP(i)` inside index lists
to a flat token (`map@i`) in a preprocessing pass** before
formula_normalize; binders recognize the `@` form (mapping id +
domain index name). No tokenizer surgery.

**Conditions today.** Formula quantifier conditions EXIST
(formula.c:1714-1782, 1922-1980, 2028-2045) but RHS is a numeric
constant only. Conditional SUMS exist nowhere. Equations/updates
have NO condition support and FAIL SILENTLY (`strtok(NULL,")")`
yields "SET: cond", set lookup misses, calloc'd setid=0 → equation
expanded over sets[0]; sum_parse leaves sumsetid=0 and the body
binder reads elem_vals[0] — no error path). **M3 must add guards on
these silent paths regardless** (a `:` reaching an equation
quantifier or sum without support = named fatal). M3 implementation
shape: build-time boolean mask over the quantifier tuple per
statement (mappings are static after reads), consulted identically
in jacobian_preallocate's counting loop and stmt_prog_execute's fill
loop (jacobian.c:293-321) — condition on a sum prunes columns only;
equation-quantifier conditions (row-pruning, changes VecSize/
eq_addr/closure squareness) are NOT needed by the corpus and stay
fatal.

**Jacobian specifics for mapped variable indices** (`del_x(MAP(r))`):
column index at jacobian.c:300-320 (`li3` loop) needs the mapping
mode; the SAME computation must run in jacobian_preallocate or
dnnz/onnz miscount; dimsetnames recovery (textual `(all,<idx>,`
search, jacobian.c:466-492) needs a mapping-aware source (codomain
set from the declaration); dcountdim3 (jacobian.c:772-779) needs a
domain-position+map-id mode; DBBD/NDBBD border classification must
mark mapped references border (`ele_inter`) since they can cross
blocks invisibly; backsolve row↔element bijection check breaks under
mapped indices (reject backsolving variables referenced through
mappings).

**Reads.** `datafile_read_labels` (cmf_io.c:86-127) is a complete
string reader (the GTAPSETS element reader) — by_elements mapping
reads reuse it; resolve names by the same strcmp matching
subset_map_build uses (tab_parse.c:4287). Read dispatch: keyword at
cmf_io.c:625-628; `data_read_files` must skip `(by_elements)`
statements, which route to a new mapping-values pass. Call order:
after set-element resolution, before the sets/set_elems broadcast
(main.c:964-1043), so mappings exist before coefficients_read and
jacobian_preallocate; broadcast the mapping tables.

**Statement dispatch.** No statement loop — keyword cascade in
tab_preprocess (cmf_io.c:589-677) with STICKY commsyntax: an
unknown keyword is prefixed with the LAST recognized keyword and
mis-parsed silently. Add `mapping` (+ later `complementarity`)
recognition there; add to the names_validate reserved list
(tab_parse.c:1620); reader modeled on sets_read; PostSim gate
already lists `complementarity` (cmf_io.c:1433) with no
implementation behind it.

**Levels.** No machinery: `level_par` only suppresses the `p_`
rewrite (cmf_io.c:1061); `orig_level=` parsed and discarded.
C0 is real new work (auto-diff of levels equations into the linear
pipeline).

**Closure two-pass.** `closure_vals` is read-only after the
broadcast (zero writers in the drivers); exo_index numbering
(main.c:1405-1594) determines VecSize/eq_addr/preallocation/block
ordering — a state flip is a pipeline re-entry from main.c:1405
through the driver dispatch, NOT an in-place edit. The teardown
helpers needed between passes already exist and are collective-safe
(jacobian_cache_free, backsolve_cache_free, *_fastrefac_free —
main.c:1944-1950), and postsim_expose_results (main.c:2031) is the
existing "copy solution back into elem_vals and re-run the formula
engine" template for seeding pass 2. → C3 = wrap numbering→driver
in a loop with the closure-modification hook between passes;
backsolve_read's post-hoc closure_vals marking pass
(tab_parse.c:4650-4670) is the mutation template.

## 5. C0 design: minimal levels-equation linearization (2026-07-29)

**STATUS: SHIPPED solver-side 2026-07-30** (src/levels.c
tab_levels_transform + the three reader adjustments below;
.audit/levels-test-kit 13 legs; verify.sh 14/14 bit-identical,
warnings 102, all kits green; parity-plan progress log has the
record). teems-R side = the C0-R phase (lift
formula_equation/default_levels rejections, mirror bookkeeping).
Deferred with named fatals: functions in levels equations,
conditional sums, p_/c_-leading levels names (joint follow-on with
GMig2's user-declared c_* linear variables — C1 prerequisite).

Spec sources: manual 9.2 (gpd2.2.2, the TABLO conversion), 18.1/18.2.1
(differentiation rules), 10.9.1 (Formula&Equation), 11.4.8 (levels-
equation operand legality: ONLY levels variables, parameters, constants;
non-parameter coefficients and linear variables are illegal; functions
limited to SQRT/EXP/LOGE/LOG10). Corpus contract = GMig2 v3.0
(original_tabs/3617 gtap.tab): 4x `Equation (levels)` + 7x
`Formula&Equation`; expression surface `+ - * / ^ ( ) [ ]`, `sum()`,
constants, parameter coefficients, levels variables (both `(levels)`
percent and `(change,levels)` change forms; some Read from file, some
set by `Formula (initial)`; values also referenced by ordinary formulas
and by ordinary linear equations, e.g. `MIGNOSP00 = 100/MIGNOSP`,
`qops = 100/LFNOSP * c_LFNOSP`).

**Change differentiation only (TABLO's ACD mode).** Manual 9.2.6:
TABLO uses change differentiation when either side is zero, any change
variable appears, the top-level operator is +/-, or a SUM occurs —
every GMig2 levels equation and F&E meets at least one of these, so
change-diff-only reproduces TABLO exactly on the whole target corpus,
and ACD ("Always use Change Differentiation") is a sanctioned global
TABLO option, so the mode is spec-legal for any other TAB. Divergence
(TABs where TABLO would default to percent-change differentiation)
affects Johansen one-step results and convergence path, not the
extrapolated solution; recorded, not guarded.

**Rules (18.1), change mode:** d(A+B)=dA+dB; d(A-B)=dA-dB;
d(A*B)=B*dA+A*dB; d(A/B)=(1/B)*dA-(A/B^2)*dB; d(A^B)=B*A^(B-1)*dA
+LOGE(A)*A^B*dB; d(sum(j,S,A))=sum(j,S,dA); d(param)=d(const)=0.
Terminals: dX -> p_X (change var, solution already in change units) or
(X/100)*p_X (percent var). Implementation: recursive AST differentiation
producing a flat term list [(coef-expr, levels-var-ref, sum-context)];
each side emitted as `sum(j,S,...)`-wrapped `(coef)*p_x(args)` terms so
every linear-variable factor is preceded by `*`/`+`/`-`/`(` (the
equation-side LinVar scan's accepted positions; a variable in a
denominator or exponent is invisible to it). All-parameter side emits
literal `0` (constant sides compile fine).

**The naming crux and its resolution ("auto-pair").** GEMPACK converts
each levels Variable X into Coefficient(non_parameter) X + linear
Variable p_X + Update (9.2.2). TEEMS convention is variables declared
BARE and referenced as p_X/c_X, with the binder stripping `p_` and
searching coefficients FIRST — so the literal GEMPACK pair is today
forbidden (names_validate coef-vs-var + the main.c coefficient-vs-p_var
guard) and would mis-bind update RHSs. Resolution: declare BOTH as `x`
— coefficient x (value store) + variable x (level_par=true, change_real
per qualifier) — with three principled reader adjustments, each
bit-neutral on any legal current TAB because 11.2.1 name-uniqueness
makes today's lookups unambiguous:
  1. names_validate: coefficient X + variable X is legal iff the
     variable has level_par (this combination IS the GEMPACK associated
     pair; hand-written pairs get the same semantics).
  2. formula_bind_operand: tokens carrying a p_ prefix search variables
     BEFORE coefficients (a p_ prefix explicitly denotes a linear
     variable); bare tokens keep coefficient-first (value reference).
  3. data_read_files (both branches): a level_par variable match defers
     to the paired coefficient (a levels read means "read the levels
     values").
Everything else already lines up: level_par's single existing behavior
(suppress the p_-insert on bare occurrences in equation/update lines,
cmf_io.c tab_write_variables) is exactly the value-reference semantics
the pair needs; `c_X` references become `p_X` in preprocess and resolve
to bare x; closure/shock/solution names are the bare variable name,
unchanged from every other TEEMS variable.

**New pass `tab_levels_transform` (new src/levels.c),** called on the
preprocessed one-statement-per-line stream right after
tab_postsim_split, before sets_read (declarations must exist before the
readers; postsim sections already forbid variables/equations). Scan A
collects levels variables (name, change/percent, quantifier text,
argument list) honoring explicit qualifiers AND positional
`Variable (default=levels|linear|change|percent_change)` statements
(mirrors variables_read), plus parameter/non-parameter coefficient
names and linear variable names for legality checks. Zero-cost gate: no
levels variables, no `&equation`, no `equation (levels` -> return
without rewriting (goldens bit-identical by construction). Pass B
rewrites to a tmp file + rename (tab_postsim_split model):
  - levels Variable statement: passed through VERBATIM (variables_read
    re-derives level_par/change_real from the same qualifiers/defaults)
    + emit `coefficient (non_parameter) <quants> x(args) ;` + emit
    `update <quants> x(args) = p_x(args) ;` (percent) or
    `update (change) <quants> x(args) = p_x(args) ;` (change — emitted
    p_-form directly since the preprocess c_->p_ rewrite already ran;
    updates_apply's (change) branch adds the RHS, its default branch
    produces x*(1+p_x/100), both correct for the single-variable RHS).
  - `Formula&Equation [quals] NAME [#lbl#] <quants> LHS = RHS ;` ->
    `formula (initial) <quants> LHS = RHS ;` + `equation (levels) NAME
    <quants> LHS = RHS ;` (10.9.1; only (initial)/(levels) qualifiers
    accepted), the latter then differentiated like any levels equation.
    The A7 fatal in formulas_execute stays as a dead-man backstop.
  - `Equation (levels) NAME <quants> LHS = RHS ;` -> parse both sides
    (recursive descent: numbers, refs, + - * / ^, parens, sum(idx,SET,
    body)), enforce 11.4.8 legality with named fatals (non-parameter
    coefficient / linear variable / unknown name / function — SQRT,
    EXP, LOGE, LOG10 deferred with a named fatal, zero corpus uses),
    differentiate, emit `equation NAME <quants> dLHS = dRHS ;`.
  - Backstop: `(levels` surviving on any equation line after the pass
    is a named fatal (today it would mis-parse silently).
Generated-name collision: none (no new names are minted — the pair
shares the levels variable's own name; sum indices are carried
verbatim, gen_* space untouched).

**What stays fatal (named):** functions in levels equations (the four
legal ones deferred); PROD(); `Equation (default=levels)` (statement-
level (levels) is the corpus form); LINEAR_NAME=/LINEAR_VAR= (A4,
zero corpus levels TABs use them... GMig2 included); NO_SPLIT.

**Recorded pre-existing gap (C1 blocker, not C0):** user-DECLARED
c_-prefixed variables (GMig2: c_shiftlf, c_migin, c_smigin, ...) are
mangled by the preprocess c_->p_ rewrite (reference becomes p_shiftlf,
declared name stays c_shiftlf, LinVar scan strcmp misses and BREAKS the
whole scan silently -> dropped columns). Fix belongs with C1 (GMig2
e2e): declaration-time normalization (strip c_, set change_real) or
scan-side fallback; either way add the M3-style named fatal on the
scan-break path.

**Other recorded interactions:** levels values ride coefficient storage
(store_real float — same precision as all data); the levels variable's
var_store/elem_vals slots hold ordinary per-step solution data (updates
and PostSim exposure unchanged); PostSim references of a levels name
bind the coefficient = post-update levels, the GEMPACK semantics;
division in generated derivative coefficients (quotient/power rules)
rides the existing zerodivide machinery; equation sides with zero
linear variables after differentiation are legal text (`0`) — the
existing "no diagnostic on an empty row" behavior is unchanged and the
closure counts them like GEMPACK does (the levels equation contributes
one row of its quantifier size regardless).

**Gates:** verify.sh 14/14 bit-identical + warnings 102 (goldens have
no levels statements — the zero-cost gate must hold); new
.audit/levels-test-kit: manual-pinned linearizations (9.2.4 examples 1
change/percent product forms, 9.2.5 changes-form sum), GMig2-shaped
legs (param*X+Y change pair = E_RMIGSP_C shape; X*Y-X*Y products =
E_RYnmvsPPP; ratio^levels-exponent with the LOGE term = E_LFNOSPendog;
sum over quantified index = E_POPNO F&E), F&E expansion leg, mixed
percent+change, value-pinned multi-step solve with per-step levels
updates asserted via PostSim, and named-fatal legs (non-parameter
coefficient, linear variable, function, unknown name, surviving
(levels)); ASan on kit legs; all prior kits re-run.
