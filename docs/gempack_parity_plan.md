# GEMPACK TABLO Feature-Parity Plan

Status: ACTIVE (committed 2026-07-25; drafted 2026-07-24). Companion to
ROADMAP Phase 6 ("GEMPACK parity program — performance-vetoed").
Refines/expands ROADMAP item 6.4 ("Statement completion … full IF
semantics, POSTSIM, others").

## Progress log

- **2026-08-04 — C3 Complementarity accurate run landed (79770af;
  ch. 51 parity approximate + accurate end to end)**: after the
  approximate pass, the 51.7.1 closure/shock modification runs on
  pre-simulation-restored values (state 2: comp@e exogenized shocked
  to zero; state 1/3 constant/parameter bound: X shocked to the
  bound, change/percent forms; levels-variable bound: the @l/@u
  difference shocked to zero; all dummies endogenous; one-for-one
  swap so the system size is unchanged), then main re-enters the
  closure-dependent pipeline (numbering through preallocation) and
  the REQUESTED method solves the accurate run. Post-accurate
  51.5.4/51.7.5 verification (state match + bounds) is fatal by
  default; CMF statements do_approx_run/do_acc_run/state_bound_error
  join steps_approx_run/redo_steps/redo_step_min_fraction (51.6).
  Remaining ch. 51 deferrals: per-subinterval approximate+accurate
  pairs under automatic accuracy (51.7.4; the approximate run treats
  the simulation as one interval and warns) and subtotals of
  \$del_Comp (29.1).

- **2026-08-03 — C2 Complementarity approximate-run state machinery
  landed (ed2b069; design doc section 8; C1 parse/validation/derived
  statements landed same day @ 298c0f1 with the C1a p_-leading
  levels-name prerequisite)**: the 11.14 E_$comp equation is emitted
  as static TAB text over per-component weight coefficients written
  by C code from the current levels values each step (the runtime
  3-state IF lives in the values; preallocation counts structurally
  so the static text preallocates the union of all branches). New
  solve_comp_approx driver: forward Euler with per-step state
  evaluation (51.7.5 whole-plane division), the $del_Comp Newton
  correction shocked 1 in full per step (NO_SPLIT), and step redo at
  the crossing fraction (51.7.3). CMF statements `complementarity
  steps_approx_run/redo_steps/redo_step_min_fraction` (51.6);
  pre/post-sim exactness warnings (51.7.5) + state-change report
  lines (51.5.3). Closure balance per component: endogenous
  complementarity variable = active (dummy auto-exogenous, counted
  as one equation element per 11.14); exogenous = inert (endogenous
  dummy absorbs the row). teems-R mirrors the counting in
  .check_system_square, exposes composed comp@e/@l/@u, and solves an
  active quota model e2e. Remaining for full ch. 51 parity: the C3
  accurate run (51.7.1 closure/shock auto-modification + two-pass
  re-entry + subinterval redo) and the do_approx_run/do_acc_run/
  state_bound_error CMF switches.

- **2026-07-30 — C0 levels-equation linearization landed (mapping/
  complementarity design doc section 5; prerequisite for GMig2
  Complementarity)**: new `tab_levels_transform` pass (src/levels.c)
  on the preprocessed statement stream — each levels Variable becomes
  the GEMPACK associated triple (9.2.2: declaration passes through,
  value Coefficient(non_parameter) + Update appended); Formula&Equation
  expands per 10.9.1; Equation (levels) is linearized by CHANGE
  differentiation (18.1/18.2.1 — TABLO's sanctioned ACD mode; every
  corpus levels equation meets a 9.2.6 change-diff trigger anyway).
  Auto-pair naming: coefficient and variable share the levels name —
  bare tokens bind the coefficient (value), p_/c_ tokens the variable
  (column) — via three bit-neutral reader adjustments (names_validate
  level_par exemption; binder variables-first for p_ tokens;
  data_read_files level_par deferral to the pair coefficient), each
  unambiguous on legal TABs by 11.2.1 uniqueness. Ride-along: the
  no-space `Formula&Equation` spelling no longer falls through the
  sticky-keyword prepend in tab_preprocess. Named fatals: non-parameter
  coefficient / linear variable / unknown name in a levels equation
  (11.4.8), functions (SQRT/EXP/LOGE/LOG10 deferred, zero corpus uses),
  conditional sums, p_/c_-leading levels names (the p_-substring
  scanners cannot carry them — recorded with GMig2's user-declared c_*
  linear variables as ONE naming-normalization follow-on, a C1
  prerequisite). Kit `.audit/levels-test-kit` 13 legs: value-pinned
  Johansen-exact forms (percent/change products, quotient, sum,
  parameter multiplier, mixed, F&E) + Gragg-extrapolated X^E (4.84
  pinned, loge dB term compiled in-equation) + failing-assertion
  canary + 6 named fatals; legacyq feq leg repointed to the 11.4.8
  fatal. Gates: verify.sh 14/14 bit-identical x2, warnings 102, all
  kits green, ASan sweep. Observed while testing (not C0 scope): a
  bare-name shock of a dimensioned variable (`Shock lv = uniform 1;`)
  silently shocks a single element — fold into the fail-fast sweep.

- **2026-07-28 — step-5 corpus survey (decision evidence for the
  scrutiny register)**: statement-leading keyword census over 71
  unique TABs (54 archive `tab_files/` + 19 extracted from the 11
  `original_tabs/` zips, deduped by content hash; comments and
  `# labels #` stripped before matching). **Transfer: 0 uses.
  Loop/Break/Cycle: 0 uses** — recommend reject-cleanly dispositions
  for 2.6 and 2.7. **Complementarity: exactly 1 model** — GMig2 v3.0
  (migration bound `Complementarity (variable=NMIGSP, lower_bound=0)`);
  2.5 needed only if migration modeling becomes a target.
  **Mapping: 19 statements across 7 files, all GTAP-family, and
  load-bearing where present**: GTAP-E/-EP use `REGTOBLOC` (REG→BLOC)
  20x each in the emissions-trading machinery — as an index function
  (`del_nctaxb(REGTOBLOC(r))`) and in mapping-conditional sums
  (`sum{r,REG: REGTOBLOC(r) = b, CO2T(r)}`); GTAP-AEZ uses
  MACTS/MAPFORSLCOV/MAPGRAZLCOV for land-cover; classic gtap.tab has
  the same REGTOBLOC pattern; GMig2 has MAPPROD. The upstream
  gtapv7.tab core's UCOM2COMM/UACT2ACTS are the user→model
  aggregation identity (`UCOM2COMM(c) = $POS(c)` feeding the CNT
  welfare accounting) — subsumed by TEEMS's R-side aggregation
  (UCOM≡COMM here) and coupled to the rejected `$POS` intrinsic, so
  the core model needs nothing. **The named validation targets
  GDYNv3.6 and gtapep use NONE of the four majors.** Net: 2.4 Mapping
  is the only major with real demand, and the demand is specifically
  GTAP-E/-EP/-AEZ (energy/emissions/land-use satellites); its hard
  part is the mapping-valued index + `MAP(r) = b` set-condition forms,
  not the declaration.
  UNIQUENESS rule, not a resolution order (manual 11.2.1)**: "the
  duplication of a name for two different purposes is not allowed …
  you cannot use 'X1' to denote a coefficient and 'x1' to be a
  variable" (names case-insensitive). So the POP-vs-pop PostSim
  finding (PS formula RHS admits both coefficients and variables; the
  binder silently took the variable = solution zeros) is an ILLEGAL
  TAB, and the fix is enforcement, not a search order. Corpus survey
  first: zero collisions across coefficient/variable/set and zero
  same-class duplicates in the full corpus — but 553 logical-FILE
  name shares (the write-all pattern names output Files after the
  coefficients they dump) and one real coefficient/equation share
  (LULC tmng), so FILE and EQUATION names are exempt from enforcement
  (neither namespace cross-binds in this engine). `names_validate`
  (tab_parse.c, once after variables_read): fatal named errors for
  cross-class collisions (coef/var, coef/set, var/set), same-class
  duplicates, and RESERVED WORDS as identifiers (sum/if/max/min/
  id0v/… — the A4 tokenizer and 3.1 normalize would misparse them;
  corpus-clean, verified). PostSim-vs-ordinary collisions are caught
  even earlier by the Tier-0 scope-isolation pass at the split (its
  12.2.1 message cites the colliding line) — layering recorded.
  Validated: verify.sh 14/14 bit-identical, warnings 102;
  `.audit/names-test-kit/run_names_tests.sh` 12/12 (the original
  popvar case fatal+named via scope isolation; coef/var, coef/set,
  var/set, duplicate, reserved each fatal+named via names_validate);
  all eight kits green; ASan johansen-lu + re-lu clean. The "use
  distinct names in PostSim expressions" workaround note is now
  superseded by enforcement.

- **2026-07-26 — step-4 legacy qualifiers closed: IfHeaderExists
  implemented; 2.2/3.6 dispositioned; 3.2 satisfied by A1**:
  spec-first from manual 10.9.1/11.11.8. **3.9 Read (IfHeaderExists)**
  (the one implementable item): qualifier stripped per statement in
  data_read_files; an absent header skips the read silently at all
  four not-found sites (bare/quantified × coefficient/variable);
  KEY SEMANTIC FINDING — the engine runs ALL reads before ALL
  formulas, so plain skip semantics would invert the manual's idiom
  (default Formula then conditional Read: the formula would overwrite
  the read every pass). A satisfied conditional read therefore MARKS
  its coefficient (`teems_coef_ifhdr` parallel array) and
  formulas_execute SUPERSEDES formulas assigning it (nloops=0 at LHS
  resolution — statement plumbing intact, nothing evaluated),
  reproducing GEMPACK's file-order outcome for the documented idiom.
  Divergence recorded: a formula placed AFTER the conditional read
  (would win in GEMPACK) also gets superseded — pathological pattern,
  zero corpus uses of the qualifier at all. **2.2 FORMULA & EQUATION**:
  expansion needs EQUATION(LEVELS) (10.9.1) — blocked on the levels
  major, fatal-by-design with a remedy message (was mislabeled "pure
  sugar" in the tier row). **3.6 LINEAR_NAME/LINEAR_VAR**: levels-only,
  zero corpus levels-variable presence — stays A4-fatal, deferred with
  the levels machinery. **3.2 Zerodivide qualifiers**: satisfied by
  A1's flag-gated dual-class state; adoption decision open. WART
  RECORDED: fully-quantified reads (`Read (all,r,REG) X(r) from ...`)
  silently no-op even unqualified — pre-existing, not a corpus form
  (bare + partial-with-elements are), fold into the fail-fast sweep.
  Ride-along: data_read_files failure now MPI_Aborts (was `return 0`,
  fourth exit-0 wart). Validated: verify.sh 14/14 bit-identical,
  warnings 102; `.audit/legacyq-test-kit/run_legacyq_tests.sh` 8/8
  (present leg: two conditional reads land AND supersede their
  default formulas, assertion-checked sums; absent leg: no error,
  defaults stand; headerless-form still errors; F&E fatal named);
  all seven kits green; ASan+UBSan johansen-lu + re-lu clean plus a
  sanitized present-leg run.

- **2026-07-26 — 3.1 remaining intrinsics landed: ID0V, MAX, MIN,
  TRUNC0, TRUNCB; $POS deferred**: spec-first from manual
  11.5/11.5.1. Zero corpus uses (pure forward parity). **Multi-arg
  machinery** (the new piece): formula_normalize classifies `max(`,
  `min(`, `id0v(` as function parens (both switches — operator-
  preceded and expression-start); formula_compile splits the innermost
  group on brace-depth-0 commas (indices ride inside {}), compiles
  each argument through the pow/muldiv/addsub pipeline under a FRESH
  group number (gen_* temp names embed the group id — same-id
  arguments would cross-bind), then folds pairwise with new op codes
  OP_MAXF/OP_MINF/OP_ID0VF (61-63; operands always compiled temps).
  MAX/MIN take 2+ args (fold), ID0V exactly 2 (x if x≠0 else v — the
  manual's recommended zerodivide guard, and it works in equations
  where ZERODIVIDE cannot). Arity violations are named fatal errors.
  ops-array sizing gains a top-level-comma slack term at the five
  allocation sites (per-arg loads + folds exceeded the 2-per-group
  budget). **Unary**: TRUNC0 (toward zero) / TRUNCB (floor) as
  transform types 48/49, same reversed-name+boundary recognition as
  ROUND. **$POS deferred**: needs per-tuple index-position machinery
  and mostly matters inside IF conditions the engine can't express;
  clean fatal at formula_compile entry (compile failures MPI_Abort in
  formula/update executors; assertions warn-skip per policy).
  Validated: verify.sh 14/14 bit-identical, warnings 102;
  `.audit/intrinsics-test-kit/run_intrinsics_tests.sh` 7/7 with 15
  in-solver assertions (ID0V const/expression/denominator-guard forms;
  MAX/MIN 2-arg, 4-arg fold, expression args, nested, indexed
  coefficient args; TRUNC0/TRUNCB both signs — all values exactly
  representable); all six kits 73/73; ASan+UBSan johansen-lu + re-lu
  clean plus a sanitized run of the values leg. Still open from 11.5:
  RANDOM (needs the randomize CMF switch + reproducibility design),
  the statistical functions, RAS_MATRIX — none in corpus, none
  planned.

- **2026-07-26 — PostSim residuals landed (Tier 0 now feature-complete
  solver-side)**: spec-first from manual 12.2.1-12.2.3. (1) **PostSim
  Read**: the split routes `read` lines to the _ps companion (`read
  elements` stays with the declarations); `postsim_reads_execute`
  (tab_parse.c) runs them once after the solve at the PS hook —
  targets validated (PostSim Coefficient only: ordinary coefficient /
  variable / unknown are named fatal errors per 12.2.3), then
  `data_read_files` reuses the whole format machinery against scratch
  stores and only the read coefficients' element ranges are copied
  into the live elem_vals, preserving every post-simulation value.
  12.2.3's same-file rule enforced in the split (normal + PostSim
  reads of one logical file = fatal). (2) **PS-LHS validation**
  (12.2.2): `teems_ps_pass` gates a check at formulas_execute's LHS
  resolution — assigning a variable or an ordinary coefficient in a
  PostSim formula is a named fatal error; PostSim coefficients are
  known via split-recorded names (`teems_ps_coefnames`) marked onto
  a `teems_coef_is_ps` parallel array after coefficients_read; the F2
  (parameter) warning is suppressed on the PS pass (qualifiers
  ignored there per 12.2.4). (3) **Scope isolation** (12.2.1): split
  pass A collects every section-declared name (set/subset/
  coefficient/file), pass B fatals on any ordinary statement
  referencing one (word-boundary match; also catches ordinary/PostSim
  redeclaration collisions — the POP-vs-pop shadowing class).
  Ride-along: the split's -1 now MPI_Aborts (was `return 0`, third
  exit-0 wart caught by kits). Validated: verify.sh 14/14
  bit-identical, warnings 102; `.audit/postsim-test-kit/
  run_postsim_tests.sh` 14/14 (positive leg: Read+formula+assertion
  over a dedicated PS data file sums to 60 in-solver AND ordinary
  outputs stay bit-identical to the golden shk2d manifest AND
  flipping the data file flips the assertion — values provably come
  from the file; six named fatal legs); all five kits 66/66;
  ASan+UBSan johansen-lu + re-lu + shk2d clean plus a sanitized
  PS-positive run through the read/copy/assert path. Tier-0 remains:
  execution order is per-statement-kind (reads, then formulas, then
  assertions; file order within each kind) — interleavings that
  depend on a read AFTER a formula would misorder (recorded); PostSim
  ZeroDivide statements ride the _ps scan with fresh state (A1).

- **2026-07-26 — A9 second bound slot landed**: a declaration may carry
  one lower (GE/GT) and one upper (LE/LT) bound (manual 10.19.1) —
  the GDYNv3.6.tab:1771 `(ge 0, le 10)` case. Slot 1 stays in
  `array_def.gltype/glval` (the struct is fwritten verbatim to sol.var
  at main.c ~1982 = binary-locked, same reason as F2's is_param);
  slot 2 rides the new `teems_coef_gltype2/teems_coef_glval2` parallel
  arrays (calloc'd in coefficients_read, NULL elsewhere — ranks
  without the read skip slot-2 checks, matching their prior slot-1-
  only behavior). Tokenizer: first bound → slot 1, opposite-direction
  second → slot 2, same-direction duplicate → fatal named error
  (replaces A4's interim keep-last warning). Enforcement refactored
  into one `coef_range_check()` helper (formula.c) called per slot in
  all three regions (formulas + both update executors) under the
  25.4.4 CMF modes; message text unchanged. Variables: both bounds
  parse cleanly into discard slots (variable bounds were never
  enforced solver-side — pre-existing gap, recorded). Ride-along:
  coefficients_read failure now MPI_Aborts (was `return 0` — another
  exit-0 wart, kit-caught). Validated: verify.sh 14/14 bit-identical,
  warnings 102; `.audit/bounds-test-kit/run_bounds_tests.sh` 15/15
  (slot-2 upper + slot-1 lower + reversed-order + in-range + CMF-yes
  fatal + duplicate fatal); all four kits green; ASan+UBSan
  johansen-lu + re-lu + gmshk clean. A4 kit's dblbound leg updated to
  the two-slot semantics. Audit A-items now ALL closed (A1-A10).

- **2026-07-26 — A6 Default-statement rework landed**: spec-first from
  manual 10.19/10.19.1. `tab_default_value` (tab_parse.c) detects a
  Default statement and extracts its space-squeezed value (spaced
  `(default = x)` forms now work — old bug 1); the readers apply
  semantics POSITIONALLY and reversibly: variables_read handles all
  four values (linear/levels/change/percent_change; the old sticky
  only knew levels/change one-way), coefficients_read gains the real-
  coefficient PARAMETER/NON_PARAMETER default (explicit token wins,
  INTEGER→PARAMETER per 10.3), formulas_execute handles initial AND
  always. `variables_read_defaults` DELETED — it pre-applied
  `default=levels/change` to ALL variables including ones declared
  before the statement (old bug 2, the audit's position-wrong).
  `tab_defaults_validate` (cmf_io.c, once after preprocess, MPI_Abort
  on failure) closes old bug 3 (defaults silently skipped): accepts
  the supported values plus Equation linear/not_add_homotopy as
  no-ops, and FATALS with named errors on Equation levels/
  add_homotopy (linearized-only solver), Coefficient lower_bound/
  upper_bound defaults (single bound slot, A9), and any unknown
  keyword or value. PostSim sections already reject Default
  statements (12.2.4) via the split's forbidden list. Nuance recorded:
  integer-LHS formulas nominally stay INITIAL under a
  `Formula (default=always)` (10.19) — not distinguished, corpus has
  zero Default uses. Validated: verify.sh 14/14 bit-identical,
  warnings 102; `.audit/defaults-test-kit/run_defaults_tests.sh`
  12/12 — equiv leg (mid-file explicit GEMPACK defaults incl. spaced
  forms + trailing defaults after the last declaration, which the old
  pre-pass would have mis-applied) and rewrite leg (existing
  `Variable (change)` declaration expressed as default=change / bare
  decl / default=percent_change) both BIT-IDENTICAL to the golden
  manifest; four fatal legs named. A1/A4 kits re-run green. ASan+UBSan
  johansen-lu + re-lu clean. HARNESS NOTE: deliberately-aborted kit
  runs leak MPICH shm segments; the container's 64M /dev/shm filled
  and later MPI inits died with a pre-main SIGBUS in
  MPIDI_POSIX_comm_bootstrap — all kits now clear
  /dev/shm/mpich_shm_* before each run.

- **2026-07-26 — A4 tokenized qualifier parsing landed + range-test CMF
  modes (25.4.4); batch-1 residual (b) closed**: `tab_qualifiers_parse`
  (tab_parse.c) tokenizes every LEADING `(qualifier,...)` group of
  VARIABLE/COEFFICIENT declarations — the first `(all,...)` quantifier
  ends the region — replacing the enumerated strstr/strip matching in
  `variables_read`/`coefficients_read` (8 bound blocks + 6 variable
  combos + the coefficient bare-word strips deleted; the
  `(orig_level` paren-count hack removed). Semantics: change/
  percent_change, linear/levels, real/integer, parameter/non_parameter
  (token-driven now — coefficient names containing "change"/
  "parameter"/"integer" no longer corrupted = residual (b); INTEGER
  defaults to PARAMETER per 10.3), ge/gt/le/lt bounds in both `(ge 0)`
  and list-member form. Unknown tokens and unsupported-semantics
  qualifiers (NO_SPLIT, LINEAR_NAME=/LINEAR_VAR=) are FATAL with named
  errors instead of silently corrupting the quantifier walk;
  ORIG_LEVEL=/VPQTYPE= parse-and-ignore (reporting/homogeneity
  metadata). Double bounds warn loudly and keep the last (A9's second
  slot still open, no longer silent). RIDE-ALONG (kit-exposed): the
  declared-range checks were print-and-continue always; now CMF-mode-
  driven per manual 25.4.4 — `range test initial values = yes|no|warn`
  (initial = IsIni formulas pass) and `range test updated values =
  updated|extrapolated|both|no|warn` (update executors + later
  formulas passes; the three selectors all map to fatal at update-time
  checks), default WARN both legs = GEMPACK's non-auto-accuracy
  default, restoring golden behavior (the gmshk matrix-shock run
  drives VTMFSD (ge 0) slightly negative — previously an invisible
  "Error:" print, now a spec-correct warning). Validated: verify.sh
  14/14 bit-identical, warnings 102; `.audit/quals-test-kit/
  run_quals_tests.sh` 16/16 (equiv leg: explicit-default qualifiers on
  existing golden declarations leave outputs BIT-IDENTICAL to the
  golden manifest; names leg: EXCHANGE/PARAMETERS/INTEGERS coefficients
  survive with assertions; fatal legs named; range-test warn/yes/no
  matrix); ASan+UBSan clean johansen-lu + re-lu + gmshk.

- **2026-07-25 — A1 ZERODIVIDE dual-class machinery landed (flag-gated)
  + batch-1 residual (a) closed**: spec-first from manual 10.11/10.11.1.
  Mechanism: `zdiv_state` (zbz/nbz value + on/off each) tracked
  positionally by `tab_next_statement_resolved_raw` on a line copy
  (legacy shared-default branches untouched → flag-off is byte-identical);
  `formulas_execute`/`assertions_execute` reset the state per rescan
  (which also gives the PostSim pass its fresh state, manual 12.2.4),
  capture it per statement, and disable it on exit; UPDATE executors
  disable at entry (10.11.1: zerodivide never applies in equations or
  updates — the jacobian leg therefore structurally never sees it
  enabled). Gate: `-gpzerodivide 1` (default 0 = legacy conflated
  default; adoption is a re-anchor-class change). The single
  `formula_eval` division site branches on numerator==0: class default
  if on, named fatal error if off; the pow(0,neg) site stays legacy
  (not a division). Initial state per the manual's convention:
  zero-by-zero default 0, nonzero-by-zero off. Residual (a): assertion
  `# message #` is captured in `tab_preprocess` BEFORE `#`-strip +
  lowercasing and reinserted (assertions only) at write-out, so failure
  reports print the verbatim mixed-case label; leading/trailing spaces
  trimmed. Validated: verify.sh 14/14 bit-identical ×2, warnings 102;
  `.audit/zdiv-test-kit/run_zdiv_tests.sh` 9/9 (dual defaults hold
  flag-on; legacy conflation reproduced and caught by an in-solver
  assertion whose mixed-case message prints; both `off` classes abort
  flag-on with named errors and complete under legacy); ASan+UBSan
  clean on johansen-lu + re-lu (flag-off) and on all three kit
  variants under `-gpzerodivide 1` (dual-state eval + both abort legs).

- **2026-07-25 — Tier 0 POSTSIM engine (first cut) landed**:
  `POSTSIM (BEGIN)/(END)` sections work end-to-end. Mechanism: the
  phase-1 preprocessor recognizes `postsim` as a statement keyword (so
  the sticky-keyword mechanism cannot mangle the markers into
  `write postsim (begin);`), and a post-preprocess splitter
  (`tab_postsim_split`) routes section contents — declarations
  (Set/Subset/Coefficient/File) STAY in the ordinary file (single
  namespace: TEEMS separates PostSim by EXECUTION, not by symbol
  table), executables (Formula/Assertion/ZeroDivide) move to a `_ps`
  companion consumed once after the solve with the F3 exposure in
  effect (IsIni=true so initial/always are both honored-as-ignored per
  12.2.4; per-pass zerodivide state gives the separate PostSim
  defaults for free), Write/Display in sections are dropped (outputs
  ride the write-all dump), PostSim Read and the twelve forbidden
  statement kinds are fatal at split time. CMF `PostSim = no;`
  honored. Zero cost without sections (one marker scan). Validated on
  golden-shk2d: a PS section computing a PostSim coefficient from the
  exposed solution (`sum of aoall shocks` = 3.6102) passes its
  identity assertion post-solve — the same section evaluated
  pre-solve (and failed) on the pre-Tier-0 build; wrong identity
  fails rc=1 with the %% report; `PostSim = no` skips; a forbidden
  Update in a section dies cleanly at split. **teems-R side landed same day** (teems-R
  506c269): Postsim is a supported statement; the pipeline tags section
  contents (declarations by name from a raw-text scan, executables by
  marker region), validates 12.2.1 (forbidden kinds abort), re-wraps
  executables in a trailing section at finalize, routes PostSim
  coefficients to `out/postsim/` outdata entries, and `ems_compose`
  surfaces them with type `postsim`. E2e-validated: GTAPv7 + section →
  deploy → solve (teems:dev) → compose row matching R's own sum to
  storage precision. FINDING recorded: coefficient/variable name
  collisions (GTAP `POP` coefficient vs `pop` variable, case-collapsed
  by preprocessing) resolve to the VARIABLE in PostSim formulas after
  solution exposure — use distinct names in PostSim expressions; a
  spec pass on GEMPACK's resolution order (12.2.2) is queued. REMAINING for full
  Tier 0: PostSim Read (fatal-with-remedy today), PS Formula
  LHS-must-be-PS-coefficient validation (single namespace accepts
  ordinary LHS — GEMPACK-illegal, currently permitted), scope
  isolation (ordinary→PS references not rejected), teems-R side
  (pass sections through instead of `ignored_state`, generate
  Write+outdata pairs for PostSim coefficients so they land in
  `ems_compose`), and `(postsim)` Write/Display qualifiers outside
  sections (dropped by R today anyway).

- **2026-07-25 — PostSim foundation (F2+F3) landed**: general WRITE
  reconsidered OFF the critical path (see the 2.1 entry — TEEMS writes
  all coefficients via R-generated Write+outdata pairs; PostSim outputs
  ride the same dump). **F2 `(parameter)` tracking**: parallel
  `teems_coef_is_param` array captured in `coefficients_read`
  (`array_def` itself is binary-locked to the R-side sol.var parser);
  Formula (always) assignment to a parameter WARNS once per run (not
  fatal — TEEMS never enforced the distinction; escalation with 1.2b
  full semantics). Update-of-parameter validation deferred to the same.
  **F3 post-sim value exposure**: after the solve, coefficient slots
  already hold post-simulation values (the updated-data writes read
  them); `postsim_expose_results` copies the composed solution `xcf`
  into the variables' elem_vals slots so post-solve statements read
  simulation results as if variables were coefficients (manual ch.12).
  First consumer landed with it: **`(postsim)`-qualified ASSERTIONS**
  (legal inside PostSim per 12.2.4) — skipped during ordinary passes,
  evaluated once post-solve with initial/always ignored; zero cost when
  absent (one statement scan). Validated on golden-shk2d: a
  sum-of-shocks assertion that would fail pre-solve passes only because
  it is skipped there, holds post-solve at `gt 3.6` and fails at
  `gt 3.7` (rc=1 with the %% report) — proving skip, exposure and
  evaluation against genuine solution values. Remaining for Tier 0:
  POSTSIM section scoping (BEGIN/END, isolated namespace), the
  ordinary-Reads+Formulas replay decision (likely unnecessary — values
  are already post-sim at the hook), PostSim Formula/Coefficient
  statements, and the R-side Write+outdata generation for PostSim
  coefficients.

- **2026-07-25 — audit findings batch (safe subset) landed**: A2 (set
  equality no longer copies the source header — the `=` readele drives
  `set_equality_build` for every source kind; validated end-to-end with
  an equality copy of a read-from-file set, size-asserted), A3 (unknown
  equality RHS — incl. set products — and undeclared sets in the
  expression route are now FATAL instead of silently producing empty
  sets), A5 (`Coefficient (real)` stripped; `(non_intertemporal)` no
  longer substring-routes into the intertemporal branch), A7
  (`Formula & Equation` fails loudly until plan 2.2), A8 (READ
  `(IfHeaderExists)` / `from terminal` / headerless forms get clean
  errors instead of strcpy(NULL) crashes), A10 (element range
  abbreviation in explicit lists is a fatal error — legal element names
  cannot contain `-`; the intertemporal `p[a]-p[b]` branch is
  unaffected). DEFERRED with design notes: A1 ZERODIVIDE (correct
  semantics change results the goldens rely on — flag-gated machinery
  as its own follow-up, adoption at a re-anchor slot), A4 (tokenized
  variable-qualifier parsing), A6 (Default-statement rework), A9
  (second bound slot). NOTE: the headerless-READ error path exits 0 on
  some callers (the known return-0 fail-fast wart — covered by the
  standing fail-fast sweep item).
- **2026-07-25 — 3.1 intrinsics tranche landed** (see commit c7ea3a5):
  EXP/SQRT/LOG10/ROUND + the paren-free-argument and
  loge-at-expression-start whitelist fixes.

- **2026-07-25 — batch 1 (quick correctness wins) landed**: 1.1 ASSERTION
  evaluation (`assertions_execute` in formula.c: residual `lhs-(rhs)`
  through the ordinary formula engine per quantifier tuple; `(initial)`/
  `(always)` honoured; word and symbol relops; per-element failure report
  in the manual's `%%` format capped at 10 details; CMF
  `Assertions = yes|no|warn` via `cmf_assertions_mode`, abort =
  MPI_Abort; rides every `formulas_execute` site incl. per-step; zero
  cost when no assertions — one `tab_count_statements` scan). 1.2a
  `non_parameter` strip-order fix. 3.3 `maximum size <n>` excision in
  the set readers. 3.4 `Update (product)` accepted at both scans.
  3.5 trimmed from the batch: the corpus shows ZERO `(default=…)` uses,
  so the remaining Default forms move with their features (e.g.
  `default=parameter` with 1.2b). RESIDUALS RECORDED: (a) assertion
  `# message #` text is stripped by `tab_preprocess` before the
  evaluator, so failure reports label with the condition text instead —
  message preservation is a follow-up; (b) the bare-word qualifier
  strips (`parameter`/`change`/`integer`) can still corrupt lowercase
  coefficient names containing those words — proper qualifier-token
  parsing is the fix, deferred; (c) unsupported assertion condition
  forms (conditional quantifiers, uncompilable functions) warn-and-skip
  rather than abort, so pre-batch models keep running; the warning
  closes once 3.1 intrinsics land. Validated: golden-re's 10 live
  `ge`/`lt` assertions pass on golden data across all verify configs;
  inverted-assertion negative test aborts rc=1 with named elements;
  warn/no modes behave; `maximum size`+`(product)`+`(non_parameter)`
  insertions parse and reproduce bit-identical outputs.

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

## Declaration-surface audit (2026-07-25, user-requested)

Full manual-vs-parser-vs-corpus audit of the ALREADY-implemented
statements (SET/SUBSET/COEFFICIENT/VARIABLE/FILE/READ/FORMULA/EQUATION/
UPDATE/ZERODIVIDE). Conclusion: everything the 54-file corpus actually
exercises is implemented consistently with the manual; the genuinely
wrong-results-today exposure is concentrated in ZERODIVIDE. Ranked
findings (A# = audit item; corpus-use counts in brackets):

- **A1 ZERODIVIDE qualifier conflation** — one shared default serves
  both `(ZERO_BY_ZERO)` and `(NONZERO_BY_ZERO)` [116 stmts]; unqualified
  `Zerodivide off;` [18] is a silent no-op (stale default persists);
  qualified `off` resets the default to 0 instead of restoring
  divide-is-an-error; equations/updates silently use the formula
  default where GEMPACK always errors. REAL numeric divergence in
  shipped TABs. Fix: two tracked defaults + off/error state per
  qualifier. ⚠ fixing may shift results of models that relied on the
  conflation — gate carefully, possible re-anchor-class change.
- **A2 `SET B = A;` broken for read-from-file sources** — `header`
  copied without `fileid`, so elements re-read from `iodata[0]` and
  `set_equality_build` (incl. the two implied SUBSETs) never runs [0
  uses today]. Fix: don't copy the header; let the equality path run.
- **A3 silent set forms** — set product `= A x B` and data-dependent
  `= (all,i,S: cond)` produce empty/garbage sets while the run
  CONTINUES [0 uses, but standard GTAP uses data-dependent ENDWS_COMM].
  Fix: fatal unsupported-statement errors.
- **A4 Variable qualifier combos enumerated, not parsed** — only
  `(linear,change)`/`(change,linear)` handled; `(levels,change)`,
  `(percent_change,…)`, `LINEAR_NAME=`, `NO_SPLIT`, `VPQTYPE=` corrupt
  the quantifier walk [0 uses]. Fix: tokenized qualifier parsing
  (also closes the batch-1 bare-word-strip residual).
- **A5 legal no-op qualifiers hard-error** — `Coefficient (REAL)`;
  `Set (NON_INTERTEMPORAL)` substring-matches `intertemporal` and dies
  in the intertemporal branch [0 uses]. Fix: strip/boundary-match.
- **A6 Default statements unreliable ×3** — spaced `(default = x)`
  not matched; `variables_read_defaults` applies position-wrong
  (affects variables declared BEFORE the statement); coefficient/
  equation defaults skipped entirely [0 uses].
- **A7 `FORMULA & EQUATION` silently drops the equation half** and
  runs the formula as ALWAYS [0 uses]. Fix: fatal until 2.2 lands.
- **A8 READ minor forms crash** — `FROM TERMINAL`, headerless text
  `FROM FILE`, `(IfHeaderExists)` reach strcpy(NULL)/misroute [0
  uses]. Fix: NULL-guards + clean unsupported errors.
- **A9 double bound `(GE 0, LE 10)` silently drops the LE** [1 use:
  GDYNv3.6.tab:1771] — single gltype slot. Fix: second bound slot.
- **A10 element range `(c1-c5)` not expanded** — becomes ONE element,
  silently [0 uses]. Fix: fatal or implement expansion.
- Error-path robustness (printf-and-continue in subsets_read/sets_read
  expression route; `read elements` with unknown logname indexes
  iodata[niodata]) — fold into the fail-fast convention.

Also corrected by the audit: this plan's Tier-3.1 premise — the
formula engine ALREADY had ABS/LOGE/ID01 (operand-transform types
41-43); 3.1 therefore extends that mechanism (EXP/SQRT/LOG10/ROUND
landed 2026-07-25 — including fixing normalize's whitelist so
function calls with paren-free arguments work, and the missing
loge-at-expression-start case) rather than building from scratch.
Remaining 3.1 candidates: ID0V, MAX/MIN, TRUNC0/TRUNCB, $POS [all 0
corpus uses — implement with PostSim if its report formulas need them].

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
|2.1|**WRITE (general)**|◐ solver handles only `Write X to file <log> header "…"` bound to a CMF outdata logname; R **drops** Write entirely (`tablo_process.R:162`)|CORE — HAR/text/terminal, `(SET)`/`(ALLSETS)`, LONGNAME, BY_ELEMENTS, `(POSTSIM)`|**RECONSIDERED 2026-07-25 (user)**|TEEMS handles writes by **writing all coefficients**: teems-R generates a `Write X to file <log> header` line per coefficient plus a matching CMF `outdata` entry (288 in the golden-re CMF), and `outputs_write_csv` dumps every one as a CSV that `ems_compose` surfaces. User-authored Write statements are architecturally superseded — R decides what is written (currently: everything). Therefore general WRITE is NOT a PostSim prerequisite: PostSim coefficients ride the same dump (R generates the Write+outdata pairs for the PostSim section exactly as for the ordinary section). GEMPACK file-targeting forms (HAR/text/terminal, `(SET)`, LONGNAME routing) stay unimplemented; user Writes in imported TABs are replaced by the write-everything scheme. Resolves scrutiny-register items 2 and 9. Residual niceties (implement only on demand): `Write ... to terminal` as a log print; a compose-side flag marking which outputs the TAB's own Writes selected.|
|2.2|**FORMULA & EQUATION (combined)**|✅ (2026-07-30, C0)|COMMON shorthand = FORMULA(INITIAL)+EQUATION(LEVELS)|**DONE**|tab_levels_transform expands per 10.9.1 and linearizes the EQUATION(LEVELS) half by change differentiation (design doc section 5); the A7 execute-time fatal stays as a backstop for PostSim sections. Operand legality (11.4.8) enforced with named fatals.|
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
|3.2|**ZERODIVIDE completeness** (`ZERO_BY_ZERO`/`NONZERO_BY_ZERO`, `OFF`)|✓ flag-gated (A1, 32297ec)|COMMON but manual calls it error-prone, recommends ID01/ID0V|**DONE (adoption open)**|Full dual-class semantics behind `-gpzerodivide 1`; default 0 keeps the bit-anchored legacy conflation. Flipping the default is a re-anchor-class decision. ID01/ID0V (3.1) are the modern path and are implemented.|
|3.3|**Set `MAXIMUM SIZE`**|✗|**OBSOLETE** — modern GEMPACK silently ignores it (§11.7.2)|**SKIP (parse-and-ignore)**|Accept the keyword and discard it so legacy TABs parse; implement no semantics. Trivial, removes a hard failure.|
|3.4|**UPDATE `(product)` qualifier**|◐ product is TEEMS's internal default; keyword not parsed|PRODUCT is the default update qualifier|**PARITY (trivial)**|Accept `(product)` as an explicit no-op default so TABs that write it parse.|
|3.5|**`Keyword (DEFAULT=…)` statements**|◐ `default=levels/change` handled; other forms partial|COMMON|**PARITY**|🔍 Ensure all forms parse: COEFFICIENT/VARIABLE/FORMULA/EQUATION default resets, and default LOWER/UPPER_BOUND. Low cost, prevents surprise aborts.|
|3.6|**Variable `LINEAR_NAME=` / `LINEAR_VAR=`**|✗ fatal (A4)|Override the auto `p_`/`c_`; associate with an existing linear var|**DEFERRED (levels variables)**|Levels-only qualifiers with ZERO corpus presence of levels VARIABLES at all; the naming override only matters once real levels-variable machinery exists. The A4 tokenizer fatals with a named remedy.|
|3.7|**Variable `NO_SPLIT` / `VPQType`**|✗|NO_SPLIT internal (MCP/Newton); VPQType for homogeneity testing (ch.57)|**QUESTION → defer**|🔍 NO_SPLIT couples to MCP (2.5). VPQType is a testing-only annotation. Both defer.|
|3.8|**FILE output qualifiers** (TEXT/GAMS, ROW/COL/SPREADSHEET/SSE, SEPARATOR, FOR_UPDATES)|◐ binding is CMF-driven, TAB File qualifiers mostly ignored|Mixed core/niche|**QUESTION**|🔍 Implement **text/SSE output** only if PostSim Write-to-text needs it; GAMS export is niche → SKIP. FOR_UPDATES relates to updated-data files.|
|3.9|**READ `IfHeaderExists` / `BY_ELEMENTS`**|✓ IfHeaderExists (2026-07-26); BY_ELEMENTS ✗|Conditional read / set-mapping element read|**IfHeaderExists DONE**|Bare header-form reads: absent header skips silently; a satisfied read SUPERSEDES formulas assigning that coefficient (reproduces GEMPACK file-order for the default+conditional-read idiom under pass-ordered execution). BY_ELEMENTS pairs with Mapping (2.4).|
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
