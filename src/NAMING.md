# Naming map (phase 4)

Terminology sources:
- **[GM]** GEMPACK Manual, Horridge et al. 2024 (teems-docs/solver_optimization)
- **[HK16]** Ha & Kompas 2016, *Economic Modelling* 52 — SBBD ordering (teems-docs/sbbd.pdf)
- **[KH19]** Kompas & Ha 2019, *Economic Modelling* 80 — NDBBD (teems-docs/ndbbd.pdf)
- **[D20]** Dixon Handbook ch. 20, *Solution Software* (teems-docs/dixon_handbook)

Conventions: snake_case; module prefix instead of `ha_`/`hcge_`/`hnew_`/`hlin_`;
the papers' mathematical symbols (`bivi` = B_i·V_i, `netcut`, interface matrix)
are kept and documented rather than paraphrased.

## Files

| current | new | rationale |
|---|---|---|
| ha_cgeglobal.h | teems_solver.h | umbrella header (split per-module in 4.2) |
| ha_cgefparse.c | str_util.c | case-insensitive search/replace helpers |
| ha_cgeiof.c | cmf_io.c | command (CMF) file reading, output staging [GM "command file"] |
| ha_cgetab.c | tab_parse.c | TAB-language parsing [GM "TABLO input file"] |
| ha_newfparse.c | formula.c | FORMULA compile/eval, UPDATE application [GM] |
| ha_newmfparse.c | (split, 4.2) jacobian.c, block_order.c, block_solve.c, solve_drivers.c | [HK16]/[KH19] architecture |
| ha_mp48.f90 | hsl_kernels.f90 | wrappers around HSL_MP48/MA48/MA51/MC66 |
| globals.c, main.c | (keep) | |

Binary `hsl` → `teems-solver` is deferred: path is hardcoded in teems-R
(`cmd_construct.R`) and the Dockerfiles; needs a coordinated change.

## Typedefs and structs (teems_solver.h)

| current | new | rationale |
|---|---|---|
| ha_cgetype | solve_real (double) | solution precision [GM real(8)] |
| ha_floattype | store_real (float) | coefficient storage precision |
| uvdim | dim_t | set size / dimension counter |
| uvadd | offset_t | element offset (long) |
| uvaddshort | exo_idx_t | exogenous-variable index width |
| forint | fortran_int | Fortran interop integer |
| hcge_iodata | cmf_file_entry | logical name → path from CMF [GM "file statement"] |
| ha_cgeset | set_def | SET statement [GM] |
| ha_cgesetele | set_element (.setsh → .superset_pos) | element + position within each superset [GM "subset mapping"] |
| hcge_cof | array_def (.antidims → .strides, .begadd → .offset, .matsize → .nelem) | COEFFICIENT/VARIABLE declaration [GM] |
| hcge_sumcof | sum_def | SUM quantifier instance [GM] |
| ha_cgevar | elem_value (.varval → .value, .var0 → .initial, .csolpupd → .substep_base) | per-element value record |
| ha_cgesumele | sum_value | partial-sum element |
| ha_cgecofele | elem_store | float-valued element (read staging) |
| ha_cgeexovar | closure_entry (.ShockId → .is_exogenous, .ShockVal → .shock_value, .ExoIndx → .exo_index) | closure + shock per element [GM "closure", "shock"] |
| ha_cgesetindx | quantifier (.arIndx → .index_name) | (all,r,REG) quantifier [GM] |
| hcge_linvars | eq_var_ref | linear-variable reference in an equation, with lead/lag |
| hcge_calvars | formula_op | one op of the compiled formula program |
| ha_cgemvar1 | datafile_labels | set-element labels from data file header |

## Globals

| current | new | rationale |
|---|---|---|
| temdir | scratch_dir | |
| mymaxnumthrd | max_threads | |
| smallthreads | section_threads | thread cap for selected OMP sections |
| step1,step2,step3 | steps1,steps2,steps3 | Gragg step counts 2-4-8 [GM "multi-step"] |
| kindx1,kindx2 | step_ratio2,step_ratio3 | steps2/steps1, steps3/steps1 |
| kval1,kval2,kval3 | extrap_w1,extrap_w2,extrap_w3 | Richardson extrapolation weights [GM "extrapolation", D20] |
| HA_COMM,HA1_COMM | node_comm,node_tail_comm | per-node split; ranks-last-on-node group |

## Enums replacing magic codes

```c
enum matrix_method  { MM_LU=0, MM_SBBD=1, MM_DBBD=2, MM_NDBBD=3 };   /* [HK16][KH19]; matches teems-R */
enum solution_method{ SM_MODIFIED_MIDPOINT=1, SM_JOHANSEN=10, SM_NOSOLVE=100 }; /* [GM "Gragg"; D20; Pearson 1991] */
enum bound_type     { BT_GE=1, BT_GT=2, BT_LE=3, BT_LT=4 };          /* hcge_cof.gltype */
/* formula_op.op and operand type codes get named constants likewise */
```

## Functions

### str_util.c (was ha_cgefparse.c)
| current | new |
|---|---|
| ha_chrfrall | str_replace_char_all |
| ha_cgefrchr | str_replace_char |
| ha_revstrpbrk | str_rfind_any |
| ha_cgerevfind | str_rfind_ci |
| ha_cgenchf | str_count_char |
| ha_cgenfind | str_count_ci |
| ha_cgeeqfind | str_rfind_toplevel |
| ha_cgefind | str_find_ci (from tab_parse) |
| ha_cgefrstr / ha_cgefrstr1 | str_replace_all / str_replace_first (from tab_parse) |
| ha_cgefrstrvbz / vbz1 | str_replace_all_bounded / str_replace_first_bounded |
| ha_cgedrcmt | str_strip_comment |

### cmf_io.c (was ha_cgeiof.c)
| current | new |
|---|---|
| hcge_niodata | cmf_count_files |
| hcge_rcmd | cmf_read |
| ha_cgerdvar1 | datafile_read_header_info |
| ha_cgermvar1 | datafile_read_labels |
| hcge_rsetname | tab_read_set_name |
| hcge_wtab | tab_preprocess |
| ha_csumindx | sum_dedup_indices |
| hcge_wvar | tab_write_variables |
| hcge_wdata | outputs_write_csv |

### tab_parse.c (was ha_cgetab.c)
| current | new |
|---|---|
| ha_cgenset / ha_cgerset | sets_count / sets_read |
| hcge_rinterset | sets_read_intertemporal |
| ha_setunion / ha_setplus / ha_setminus | set_union_named / set_union_op / set_difference |
| ha_cgersubset | subsets_read |
| ha_cgesubsetchck | subset_map_build |
| ha_cgeralltime | set_find_alltime |
| ha_cgencof | tab_count_statements |
| hcge_rcof / hcge_rvar | coefficients_read / variables_read |
| hcge_defvar | variables_read_defaults |
| hcge_rcofele | coef_resolve_sets |
| hcge_rexo | closure_read |
| hcge_rshock | shocks_read |
| hcge_readff | data_read_files |
| ha_cgercls | closure_next_statement |
| ha_cgertabl / ha_cgertabl1 | tab_next_statement / tab_next_statement_resolved |
| hcge_dsum / hcge_nsum | sum_parse / sum_count |
| ha_cgecutsum | sum_extract |
| ha_cgerecovar | formula_normalize |
| hcge_repllin / hcge_rlinzero | eq_replace_linvar / eq_zero_linvar |

### formula.c (was ha_newfparse.c)
| current | new |
|---|---|
| ha_newfparse | formula_compile |
| ha_newfppow / ha_newfpmuldiv / ha_newfpplumin / ha_newfpif | formula_compile_pow / _muldiv / _addsub / _if |
| ha_newfpcal | formula_eval |
| hnew_calcff | formulas_execute |
| hnew_calsum | sum_eval |
| hnew_simplrpl | formula_subst_scalar |
| hnew_varrepl | formula_bind_operand |
| hnew_intrpl | leadlag_encode |
| hnew_arset | parse_index_leadlag |
| hnew_update | updates_apply (midpoint flag) [GM "update"] |
| hnew_gupd | updates_apply_product |

### jacobian.c / block_order.c / block_solve.c / solve_drivers.c (was ha_newmfparse.c)
| current | new | rationale |
|---|---|---|
| HaNewMatVal | jacobian_fill | first-order derivative matrix [HK16 §1] |
| NewMatreadele | jacobian_preallocate | nnz counts for PETSc preallocation |
| NewMatvarRead | equation_order_read | variable/equation ordering [HK16 §5] |
| NestedMatvarRead | equation_order_read_nested | [KH19] |
| hlin_dsum / hlin_calsum / hlin_replsum / hlin_rlinone | eq_sum_parse / eq_sum_eval / eq_sum_replace / eq_linvar_read | linearized-equation terms |
| HaDBBDMatOder | dbbd_order | doubly bordered block diagonal [KH19] |
| HaNDBBDMatOderPre / HaNDBBDMatOder | ndbbd_order_presolve / ndbbd_order | [KH19] |
| HaDBBDParSol | dbbd_solve | |
| HaNDBBDParPre / HaNDBBDParSol | ndbbd_presolve / ndbbd_solve | pre-solve/back-solve [KH19 steps 1–5] |
| NDBBD_sol / NDBBD_sol_nread1 | ndbbd_block_solve / ndbbd_block_solve_mem | disk-staged vs in-memory block solve |
| HaReduce / HaReduceNoComp | reduce_to_rank / reduce_to_rank_nocompress | chunked MPI reduction of B_i·V_i products |
| Johansen | solve_johansen | [GM; D20; Johansen 1960] |
| ModMidPoint | solve_modified_midpoint | Gragg's modified midpoint [GM "Gragg"; Pearson 1991] |

Variables named after the papers' symbols are kept with documentation:
`vecbivi` (B_i·V_i product, [KH19] step 1c), `netcut` [HK16 §4],
`biviindx`, `insize` (solver metadata block), `counteq*`, `countvarintra*`.

All ha_*/hcge_*/hnew_*/hlin_* identifiers are gone, including parameters,
locals, and the generated-name prefixes (now gen_sum/gen_par/gen_pow/
gen_mul/gen_add). Pure loop temporaries (i1, j2, l2...) remain.
they are renamed only where a function is otherwise touched.
