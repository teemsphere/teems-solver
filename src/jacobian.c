#include <teems_solver.h>
#include <hsl_kernels.h>

/* ===== compiled equation programs (roadmap 6.2, small-model lever) =====
 *
 * jacobian_fill used to re-derive every equation's parse artifacts --
 * quantifier bindings, linear-variable references, SUM bodies, the
 * per-occurrence coefficient ops programs -- from the TAB text on every
 * per-step refill.  The structures below hold those artifacts explicitly,
 * splitting each statement into a build phase (parse + formula_compile,
 * step-invariant) and an execute phase (SUM evaluation + element loop +
 * MatSetValues, per-step values).  Ops programs bind coefficients,
 * variables and sums as offsets into elem_vals/sum_vals (see
 * formula_bind_operand), so a built program remains valid when values
 * change between steps.
 */

/* one SUM(...) body: compiled ops + its quantifier frame; evaluation
   writes the summed values into sum_vals[base .. base+nloops) */
typedef struct
{
  formula_op *ops;
  int nops;
  quantifier *arSet;                /* nouter carried dims + summed index */
  dim_t nouter;                     /* carried dims (= sum_cof size) */
  offset_t nloops;                  /* product of carried dim sizes */
  offset_t dcountdim[4*MAXVARDIM];  /* carried-dim strides */
  offset_t sumset_size;             /* elements of the summed set */
  offset_t base;                    /* first target slot in sum_vals */
} sum_prog;

/* one linear-variable occurrence: the compiled coefficient program */
typedef struct
{
  bool built;                       /* false = rows outside this rank */
  offset_t LinVarIndx;
  dim_t fdimlin;
  offset_t nloopsfac;
  quantifier *arSet;                /* fdimlin quantifiers */
  offset_t dcountdim2[4*MAXVARDIM];
  dim_t dcountdim3[MAXVARDIM];      /* variable dim -> arSet position */
  dim_t supset[MAXSUPSET];          /* superset_pos column per dim, 0 = direct */
  int dcountmap[MAXVARDIM];         /* >0: dim routes via teems_maps[id-1] (M2b) */
  dim_t dimleadlag[MAXVARDIM];
  sum_prog *sums;                   /* coefficient-expression sums */
  int nsums;
  formula_op *ops;
  int nops;
} linvar_prog;

/* one equation statement */
typedef struct
{
  bool inproc;                      /* any row of the block on this rank */
  offset_t nloops;                  /* rows in this equation block */
  solve_real zerodivide;            /* refreshed per fill (may reference a scalar coefficient) */
  sum_prog *eqsums;                 /* linear-variable-free sums */
  int neqsums;
  sum_value *sum_vals;              /* shared value store for all sums */
  linvar_prog *lv;
  int nlv;
} stmt_prog;

/* Build one innermost SUM of formulain into *out: compile the body,
   record the quantifier frame and target slots, and replace the SUM text
   with its generated-name reference.  This is the build half of the
   former eq_sum_eval/sum_eval pair; evaluation happens per step in
   sum_prog_eval.  skip_linvar_sums skips bodies referencing linear
   variables (the equation-level pass; those sums are handled per
   occurrence).  Returns 1 when a sum was built, 0 when none remain. */
static int sum_prog_build(char *formulain, char *commsyntax, bool skip_linvar_sums, sum_def *sum_cof, int j, set_def *sets, array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar, offset_t ncofele, int totalsum, formula_op *ops, int *sumindx, sum_prog *out) {
  char *readitem,*p;
  char interchar[NAMESIZE],line[TABREADLINE],line1[TABREADLINE];
  int nops,length,k=0,k1=0,i=0;
  dim_t dcount,l;
  offset_t nloops;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=str_find_ci(readitem,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0||formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      sum_extract(line);
      k1=str_find_ci(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else if (skip_linvar_sums&&(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1)) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        strcpy(line1,line);
        p=strtok(line,",");
        p=strtok(NULL,",");
        p=strtok(NULL,"\0");
        p[strlen(p)-1]='\0';
        out->nouter=sum_cof[j].size;
        out->arSet= (quantifier *) calloc (sum_cof[j].size+1,sizeof(quantifier));
        for (l=0; l<sum_cof[j].size; l++) {
          out->arSet[l].setid=sum_cof[j].setid[l];
          strcpy(out->arSet[l].index_name,sum_cof[j].dimnames[l]);
        }
        nloops=1;
        for (l=0; l<sum_cof[j].size; l++) {
          nloops=nloops*sets[out->arSet[l].setid].size;
          dcount=sum_cof[j].size-l;
          if(dcount==sum_cof[j].size) {
            out->dcountdim[dcount-1]=1;
          }
          else {
            out->dcountdim[dcount-1]=out->dcountdim[dcount]*sets[out->arSet[dcount].setid].size;
          }
        }
        out->nloops=nloops;
        out->arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
        strcpy(out->arSet[sum_cof[j].size].index_name,sum_cof[j].sumindx);
        out->sumset_size=sets[sum_cof[j].sumsetid].size;
        nops=0;
        if(!formula_compile(p,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,out->arSet,(dim_t)(sum_cof[j].size+1)))MPI_Abort(PETSC_COMM_WORLD,1);
        out->ops= (formula_op *) malloc (nops*sizeof(formula_op));
        memcpy(out->ops,ops,nops*sizeof(formula_op));
        out->nops=nops;
        out->base=*sumindx;
        *sumindx=*sumindx+nloops;
        strcpy(interchar,sum_cof[j].sumname);
        strcat(interchar,"{");
        for (l=0; l<sum_cof[j].size; l++) {
          strcat(interchar,sum_cof[j].dimnames[l]);
          strcat(interchar,",");
        }
        if (interchar[strlen(interchar)-1]==',') {
          interchar[strlen(interchar)-1]='}';
        }
        else {
          if (interchar[strlen(interchar)-1]=='{') {
            interchar[strlen(interchar)-1]='\0';
          }
          else {
            strcat(interchar,"}");
          }
        }
        while(str_replace_all(formulain,line1,interchar)!=NULL);
        return 1;
      }
    }
    else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

/* Evaluate one built SUM into sum_vals: same OpenMP layout, thread-copy
   pattern and per-element arithmetic order as the former inline loops. */
static void sum_prog_eval(sum_prog *sp, set_def *sets, set_element *set_elems, elem_value *elem_vals, sum_value *sum_vals, solve_real zerodivide) {
  dim_t dcount;
  offset_t l,l1,l2,superset_pos;
  solve_real vval;
  quantifier *arSet2=NULL;
  formula_op *ops1=NULL;
  #pragma omp parallel private(l,l1,l2,dcount,superset_pos,vval,arSet2,ops1) shared(elem_vals,sum_vals,sp)
  {
  if(omp_get_thread_num()!=0){
    arSet2=malloc((sp->nouter+1)*sizeof(quantifier));
    memcpy(arSet2,sp->arSet,(sp->nouter+1)*sizeof(quantifier));
    ops1=malloc(sp->nops*sizeof(formula_op));
    memcpy(ops1,sp->ops,sp->nops*sizeof(formula_op));
  }else{
    ops1=sp->ops;
    arSet2=sp->arSet;
  }
    /* workers must finish copying the shared scratch (master aliases it)
       before the master starts mutating per-iteration state */
    #pragma omp barrier
  #pragma omp for
  for (l=0; l<sp->nloops; l++) {
    l2=l;
    for (dcount=0; dcount<sp->nouter; dcount++) {
      superset_pos=l2/sp->dcountdim[dcount];
      arSet2[dcount].indx=superset_pos;
      l2=l2-superset_pos*sp->dcountdim[dcount];
    }
    vval=0;
    for (l1=0; l1<sp->sumset_size; l1++) {
      arSet2[sp->nouter].indx=l1;
      vval+=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,sp->nops,arSet2,(dim_t)(sp->nouter+1),zerodivide);
    }
    sum_vals[sp->base+l].value=vval;
  }
  if(omp_get_thread_num()!=0){
    free(arSet2);
    arSet2=NULL;
    free(ops1);
    ops1=NULL;
  }else{
    ops1=NULL;
    arSet2=NULL;
  }
  }
}

/* Free one statement's built programs. */
static void stmt_prog_free(stmt_prog *st) {
  int i,s;
  for (s=0; s<st->neqsums; s++) {
    free(st->eqsums[s].ops);
    free(st->eqsums[s].arSet);
  }
  free(st->eqsums);
  for (i=0; i<st->nlv; i++) {
    for (s=0; s<st->lv[i].nsums; s++) {
      free(st->lv[i].sums[s].ops);
      free(st->lv[i].sums[s].arSet);
    }
    free(st->lv[i].sums);
    free(st->lv[i].ops);
    free(st->lv[i].arSet);
  }
  free(st->lv);
  free(st->sum_vals);
  memset(st,0,sizeof(stmt_prog));
}

/* Per-rank cache of the built statement programs.  The equation
   structure is fixed for the whole solve -- only coefficient values
   change between steps -- so the first jacobian_fill call builds every
   statement's programs and later calls only re-evaluate them.  The
   cache is keyed to the matrix row ownership range; if a caller ever
   presents a different layout the cache rebuilds from scratch. */
static stmt_prog *stmt_cache=NULL;
static int stmt_cache_n=0;
static bool stmt_cache_built=false;
static PetscInt stmt_cache_Istart=-1,stmt_cache_Iend=-1;

void jacobian_cache_free(void) {
  int i;
  for (i=0; i<stmt_cache_n; i++) stmt_prog_free(&stmt_cache[i]);
  free(stmt_cache);
  stmt_cache=NULL;
  stmt_cache_n=0;
  stmt_cache_built=false;
}

/* Execute one built statement: evaluate its SUM programs in build order
   (equation-level first, then per occurrence -- the order they were
   evaluated inline before the split), then run each occurrence's element
   loop filling A (endogenous columns) and B (exogenous columns).  Loop
   layout, OpenMP scheduling and MatSetValues pattern are unchanged. */
static void stmt_prog_execute(stmt_prog *st, offset_t matrow, offset_t *eq_addr, offset_t nloops, set_def *sets, set_element *set_elems, elem_value *elem_vals, closure_entry *closure_vals, array_def *vars, PetscInt Istart1, PetscInt Iend1, Mat A, Mat B) {
  int s,i;
  dim_t dcount;
  offset_t lj,l1,l2,li3,i3,i5,sj;
  PetscInt Iindx,Jindx;
  solve_real vval;
  quantifier *arSet1=NULL;
  formula_op *ops1=NULL;
  if(!st->inproc) return;
  for (s=0; s<st->neqsums; s++) sum_prog_eval(&st->eqsums[s],sets,set_elems,elem_vals,st->sum_vals,st->zerodivide);
  for (i=0; i<st->nlv; i++) {
    linvar_prog *lv=&st->lv[i];
    if(!lv->built) continue;
    for (s=0; s<lv->nsums; s++) sum_prog_eval(&lv->sums[s],sets,set_elems,elem_vals,st->sum_vals,st->zerodivide);
    #pragma omp parallel private(lj,Jindx,i3,sj,i5,l2,dcount,l1,li3,Iindx,arSet1,ops1,vval) shared(elem_vals,st,lv,closure_vals,vars,eq_addr)
    {
    if(omp_get_thread_num()!=0){
      arSet1=malloc(lv->fdimlin*sizeof(quantifier));
      memcpy(arSet1,lv->arSet,lv->fdimlin*sizeof(quantifier));
      ops1=malloc(lv->nops*sizeof(formula_op));
      memcpy(ops1,lv->ops,lv->nops*sizeof(formula_op));
    }else{
      ops1=lv->ops;
      arSet1=lv->arSet;
    }
    /* workers must finish copying the shared scratch (master aliases it)
       before the master starts mutating per-iteration state */
    #pragma omp barrier
      solve_real *value= (solve_real *) calloc (lv->nloopsfac,sizeof(solve_real));
      PetscInt *jcn= (PetscInt *) calloc (lv->nloopsfac,sizeof(PetscInt));
      solve_real *valueb= (solve_real *) calloc (lv->nloopsfac,sizeof(solve_real));
      PetscInt *jcnb= (PetscInt *) calloc (lv->nloopsfac,sizeof(PetscInt));
    #pragma omp for
      for (i5=0; i5<nloops; i5++) {
        Jindx=eq_addr[matrow+i5];
        if(Jindx>=Istart1&&Jindx<Iend1) {
          i3=0;
          sj=0;
          for (lj=i5*lv->nloopsfac; lj<(i5+1)*lv->nloopsfac; lj++) {
          l2=lj;
          for (dcount=0; dcount<lv->fdimlin; dcount++) {
            l1=(offset_t) l2/lv->dcountdim2[dcount];
            arSet1[dcount].indx=l1;
            l2=l2-l1*lv->dcountdim2[dcount];
          }
          li3=0;
          for (dcount=0; dcount<vars[lv->LinVarIndx].size; dcount++) {
            if(lv->dcountmap[dcount]>0) {
              /* mapped dim (M2b): the loop runs the domain; the column
                 sits at the mapping's codomain position */
              li3=li3+((offset_t)teems_maps[lv->dcountmap[dcount]-1].values[arSet1[lv->dcountdim3[dcount]].indx])*vars[lv->LinVarIndx].strides[dcount];
            }
            else if(lv->supset[dcount]==0) {
              li3=li3+(arSet1[lv->dcountdim3[dcount]].indx+lv->dimleadlag[dcount])*vars[lv->LinVarIndx].strides[dcount];
            }
            else {
              li3=li3+(set_elems[sets[arSet1[lv->dcountdim3[dcount]].setid].offset+arSet1[lv->dcountdim3[dcount]].indx].superset_pos[lv->supset[dcount]]+lv->dimleadlag[dcount])*vars[lv->LinVarIndx].strides[dcount];
            }
          }
          vval=formula_eval(elem_vals,sets,set_elems,st->sum_vals,ops1,lv->nops,arSet1,lv->fdimlin,st->zerodivide);
          Iindx=closure_vals[vars[lv->LinVarIndx].offset+li3].exo_index;
          if (!closure_vals[vars[lv->LinVarIndx].offset+li3].is_exogenous&&vval!=0) {
            value[i3]=vval;
            jcn[i3]=Iindx;
            i3++;
          }
          if (closure_vals[vars[lv->LinVarIndx].offset+li3].is_exogenous&&vval!=0) {
            valueb[sj]=-vval;
            jcnb[sj]=Iindx;
            sj++;
          }
        }
        if(i3>0)MatSetValues(A,1,&Jindx,i3,jcn,value,ADD_VALUES);
        if(sj>0)MatSetValues(B,1,&Jindx,sj,jcnb,valueb,ADD_VALUES);
        }
      }
      free(value);
      free(jcn);
      free(valueb);
      free(jcnb);
    if(omp_get_thread_num()!=0){
      free(arSet1);
      arSet1=NULL;
      free(ops1);
      ops1=NULL;
    }else{
      ops1=NULL;
      arSet1=NULL;
    }
    }
  }
}

/* One index token of a linear-variable reference (shared by the
   statement builder, the preallocation scan and the two ordering
   scans, whose per-dim recovery blocks were verbatim copies): decode
   the lead/lag, then recover the set the dim loops over -- from the
   (all,idx,SET) quantifier when idx is quantified, else from the
   textually nearest enclosing sum(idx,SET,...) before the reference
   (lvar = the reference's offset in the statement text, the original
   occurrence disambiguation).  A lowered mapping token map@idx
   (design doc M2b) splits instead of searching: in split mode -- the
   builder/preallocation, whose loops iterate the DOMAIN and route the
   column offset through the mapping's value table -- dimnames gets
   the bare index, dimsetnames the domain set and dimmapid the
   mapping id; in keep mode -- the ordering scans -- the whole token
   is kept, so it matches no quantifier and the reference classifies
   as crossing blocks (the conservative direction), and dimsetnames
   gets the codomain so border marking resolves it to the declared
   set's full range. */
static void linvar_dim_read(char *p, char *linecopy, offset_t lvar,
                            eq_var_ref *ref, dim_t d, bool split_mapped,
                            set_def *sets) {
  char lintmp[TABREADLINE],*p1;
  offset_t l1,lvar3;
  int leadlag=0,lvar1,lvar2,lvar4,mp=0;
  if(p==NULL)p="";
  parse_index_leadlag(p,&leadlag);
  ref->dimleadlag[d]=leadlag;
  ref->dimmapid[d]=0;
  if(strchr(p,'@')!=NULL) {
    char *idx=mapping_token_split(p,&mp);
    if(split_mapped) {
      /* the loop set is whatever the domain index ranges over -- fall
         through to the textual recovery with the bare index so a sum
         over a SUBSET of the domain resolves to the subset and trips
         the exact-domain fatal in linvar_map_dim_check rather than
         silently looping the full domain (M2c) */
      ref->dimmapid[d]=mp;
      p=idx;
    }
    else {
      *(idx-1)='@';
      strcpy(ref->dimnames[d],p);
      strcpy(ref->dimsetnames[d],sets[teems_maps[mp-1].toset].setname);
      return;
    }
  }
  strcpy(ref->dimnames[d],p);
  strcpy(lintmp,"(all,");
  strcat(lintmp,p);
  strcat(lintmp,",");
  l1=str_find_ci(linecopy,lintmp);
  if (l1>-1) {
    p1=&linecopy[0]+l1;
    strncpy(ref->dimsetnames[d],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
  }
  else {
    strcpy(lintmp,"sum(");
    strcat(lintmp,p);
    strcat(lintmp,",");
    lvar1=str_count_ci(linecopy,lintmp);
    lvar3=str_find_ci(linecopy,lintmp);
    if (mp>0&&lvar3==-1) {
      /* mapped index with no quantifier or enclosing sum in sight:
         fall back to the declared domain (cannot crash on the miss) */
      strcpy(ref->dimsetnames[d],sets[teems_maps[mp-1].fromset].setname);
      return;
    }
    if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
        lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
        if (lvar4>-1&&lvar4<lvar) {
          lvar3=lvar3+lvar4+4;
        }
        else {
          break;
        }
      }
    p1=&linecopy[0]+lvar3;
    p1=p1+strlen(lintmp);
    strncpy(ref->dimsetnames[d],p1,strchr(p1,',')-p1);
  }
}

/* Validate one mapped dimension of a linear-variable reference and
   disable superset routing for it (the column offset goes through the
   mapping's value table instead): the loop index must range over the
   mapping's domain set exactly and the declared dim set must be the
   codomain exactly -- the same contract map_dim_bind enforces on the
   coefficient side -- and a lead/lag offset through a mapping has no
   meaning.  Named fatals, not mis-binds. */
static void linvar_map_dim_check(eq_var_ref *ref, dim_t d, offset_t frame_setid,
                                 array_def *vars) {
  map_def *md=&teems_maps[ref->dimmapid[d]-1];
  if((offset_t)md->fromset!=frame_setid) {
    printf("Error: the index of mapping %s does not range over its domain set (in %s); subset routing around a mapped argument is not supported\n",md->mapname,ref->LinVarName);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  if((offset_t)md->toset!=vars[ref->LinVarIndx].setid[d]) {
    printf("Error: mapping %s does not map into the argument set at that position of %s (manual 11.9.7)\n",md->mapname,ref->LinVarName);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  if(ref->dimleadlag[d]!=0) {
    printf("Error: a lead/lag offset on the mapped index of %s (mapping %s) is not supported\n",ref->LinVarName,md->mapname);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  md->used=true;
}

/* Build one equation statement's compiled programs into *stp (the build
   half of the former jacobian_fill loop body; the statement text in line
   is consumed).  force_all treats every row of the block as owned --
   backsolve recovery programs are evaluated outside the matrix ownership
   split -- while the fill path keeps the exact eq_addr/[Istart1,Iend1)
   pruning it always had. */
static void stmt_prog_build_one(char *line, stmt_prog *stp, char *commsyntax,
                                set_def *sets, offset_t nset, array_def *coefs, offset_t ncof,
                                array_def *vars, offset_t nvar, offset_t ncofele,
                                offset_t *eq_addr, offset_t matrow,
                                PetscInt Istart1, PetscInt Iend1, PetscMPIInt mpisize1,
                                bool force_all) {
  char tline[TABREADLINE],line1[TABREADLINE],leftline[TABREADLINE],linecopy[TABREADLINE];
  char vname[TABREADLINE],sumsyntax[NAMESIZE],lintmp[TABREADLINE];
  char *readitem=NULL,*p=NULL,*p1=NULL;
  PetscInt Jindx=0;
  bool isinproc;
  dim_t fdim,np,dcount,fdimlin=0,i4,sup,supset[MAXSUPSET];
  int totalsum,sumcount=1,sumcount1=0,lvar,lvar1,lvar2,lvar3,lvar4;
  offset_t lj,l1,i1=0,sumbegadd,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],dcountdim3[4*MAXVARDIM],nloops,nloopslin,nloopsfac,li3,nsumele,nsumele1,l2;
  int sumindx,npow,npar,nmul,nplu,ndiv,nmin,nops=0,nlinvars,leadlag,varindx1,varindx2;
  offset_t sj,l,i3,i,arsetdim=0,nops_alloc=0;

      str_replace_first(line, commsyntax, "");
      str_replace_first(line, "(linear)", "");
      while (str_replace_all(line,"  ", " "));
      while (str_replace_char(line, '[', '('));
      while (str_replace_char(line, ']', ')'));
      while (str_replace_char(line, '{', '('));
      while (str_replace_char(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=str_count_ci(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        i=str_rfind_ci(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (str_replace_all(readitem," ", ""));
      while (formula_normalize(readitem)==1);
      leadlag_encode(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      npow=str_count_char(readitem, '^');
      nmul=str_count_char(readitem, '*');
      ndiv=str_count_char(readitem, '/');
      nmul=nmul+ndiv;
      nplu=str_count_char(readitem, '+');
      nmin=str_count_char(readitem, '-');
      nplu=nplu+nmin;
      npar=str_count_char(readitem, '(')+str_count_char(readitem, ',');/* comma slack for multi-arg intrinsics (plan 3.1) */

      strcpy(line,line1);
      readitem=line;
      np=str_count_ci(readitem,"p_");
      nops_alloc=2*(npow+nmul+nplu+npar+1);
      formula_op *ops= (formula_op *) calloc (nops_alloc,sizeof(formula_op));
      eq_var_ref *LinVars= (eq_var_ref *) calloc (np,sizeof(eq_var_ref));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=str_find_ci(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;

        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          if (l>=nvar) break;
          p=strpbrk(readitem,"}+*-/^)");
          if (p==NULL) break;
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (vars[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            linvar_dim_read(p,linecopy,lvar,&LinVars[i3],0,true,sets);
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<vars[l].size-1; i4++) {
              p = strtok(NULL,",");
              linvar_dim_read(p,linecopy,lvar,&LinVars[i3],i4,true,sets);
            }
            p = strtok(NULL,"}");
            linvar_dim_read(p,linecopy,lvar,&LinVars[i3],i4,true,sets);
            break;
          }
          i3++;
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      i3=0;
      for (i=0; i<nlinvars; i++) i3+=vars[LinVars[i].LinVarIndx].size;
      arsetdim=fdim+i3;
      quantifier *arSet= (quantifier *) calloc (arsetdim,sizeof(quantifier));

      strcpy(line,linecopy);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].index_name,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
              arSet[i].setid=i4;
              break;
            }
          nloops=nloops*sets[arSet[i].setid].size;
        }
        dcountdim1[fdim-1]=1;
        for (dcount=fdim-2; dcount>-1; dcount--) {
          dcountdim1[dcount]=sets[arSet[dcount+1].setid].size*dcountdim1[dcount+1];
        }
      }
      isinproc=false;
      if(force_all) isinproc=true;
      else if(mpisize1>1) {
        for (lj=0; lj<nloops; lj++) {
          Jindx=eq_addr[matrow+lj];
          if(Jindx>=Istart1&&Jindx<Iend1) {
            isinproc=true;
            break;
          }
        }
      }
      else isinproc=true;
      stp->inproc=isinproc;
      if(isinproc) {
        strcpy(line,line1);
        readitem=line;
        strcpy(sumsyntax,"sum(");
        totalsum=sum_count(readitem,sumsyntax);
        sum_def *sum_cof= (sum_def *) calloc (totalsum*nlinvars+1,sizeof(sum_def));
        sumcount=0;
        while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim+1,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        li3=0;
        for (i4=0; i4<totalsum; i4++) {
          i1=1;
          for(sj=0; sj<sum_cof[i4].size; sj++) {
            i1=i1*sets[sum_cof[i4].setid[sj]].size;
          }
          sum_cof[i4].offset=li3;
          li3=li3+i1;
        }
        nsumele=li3;
        stp->sum_vals= (sum_value *) calloc (nsumele*nlinvars,sizeof(sum_value));
        strcpy(line,line1);
        readitem=line;
        sumcount=0;
        while (eq_sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim+1,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        li3=0;
        for (i4=0; i4<totalsum; i4++) {
          i1=1;
          for(sj=0; sj<sum_cof[i4].size; sj++) {
            i1=i1*sets[sum_cof[i4].setid[sj]].size;
          }
          sum_cof[i4].offset=li3;
          li3=li3+i1;
        }
        nsumele=li3;
        for (i3=0; i3<totalsum; i3++) {
          i1=1;
          sum_cof[i3].strides[sum_cof[i3].size-1]=1;
          for(sj=sum_cof[i3].size-2; sj>-1; sj--) {
            sum_cof[i3].strides[sj]=sum_cof[i3].strides[sj+1]*sets[sum_cof[i3].setid[sj+1]].size;
          }
        }
        sumcount=0;
        strcpy(line,line1);
        readitem=line;
        sumindx=0;
        stp->eqsums= (sum_prog *) calloc (totalsum+1,sizeof(sum_prog));
        while (sum_prog_build(readitem,sumsyntax,true,sum_cof,sumcount,sets,coefs,ncof,vars,nvar,ncofele,totalsum,ops,&sumindx,&stp->eqsums[stp->neqsums])==1) {
          sumcount++;
          stp->neqsums++;
        }
        strcpy(line1,readitem);
        sumbegadd=nsumele;
        sumcount1=sumcount;
        stp->lv= (linvar_prog *) calloc (nlinvars,sizeof(linvar_prog));
        stp->nlv=nlinvars;
        for (i=0; i<nlinvars; i++) {
          if(!force_all) {
            Jindx=eq_addr[matrow];
            if(Jindx>=Iend1)continue;
          }
          i3=0;
          nloopslin=nloops;
          if (fdim==0) {
            for (l2=0; l2<vars[LinVars[i].LinVarIndx].size; l2++) {
              strcpy(arSet[i3].index_name,LinVars[i].dimnames[l2]);
              for (i4=0; i4<nset; i4++) if(strcmp(sets[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                  break;
                }
              arSet[i3].setid=i4;
              nloopslin=nloopslin*sets[arSet[i3].setid].size;
              i3++;
            }
          }
          else for (l2=0; l2<vars[LinVars[i].LinVarIndx].size; l2++) {
              for (sj=0; sj<fdim; sj++) {
                if(strcmp(arSet[sj].index_name,LinVars[i].dimnames[l2])==0) {
                  break;
                }
                else if (sj==fdim-1) {
                  strcpy(arSet[fdim+i3].index_name,LinVars[i].dimnames[l2]);
                  for (i4=0; i4<nset; i4++) if(strcmp(sets[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                      break;
                    }
                  arSet[fdim+i3].setid=i4;
                  nloopslin=nloopslin*sets[arSet[fdim+i3].setid].size;
                  i3++;
                }
              }
            }
          nloopsfac=(offset_t)nloopslin/nloops;
          if(!force_all) {
            Jindx=eq_addr[matrow+(offset_t)(nloopslin-1)/nloopsfac];
            if(Jindx<Istart1)continue;
          }
          
          fdimlin=fdim+i3;
          if (i3>0) {
            dcountdim2[fdimlin-1]=1;
            for (dcount=fdimlin-2; dcount>-1; dcount--) {
              dcountdim2[dcount]=sets[arSet[dcount+1].setid].size*dcountdim2[dcount+1];
            }
          }
          else for (dcount=0; dcount<fdim; dcount++) {
              dcountdim2[dcount]=dcountdim1[dcount];
            }

          strcpy(leftline,line1);
          eq_zero_linvar(leftline,i);
          eq_replace_linvar(leftline,0);
          strcpy(sumsyntax,"sum(");
          eq_sum_replace(leftline,sumsyntax,i,LinVars,vars);
          eq_linvar_read(leftline,LinVars,i,vars);
          strcpy(line,leftline);
          readitem=line;
          sumcount1=sumcount;
          while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdimlin+1,sumcount1)==1) {
            sumcount1++;
          }
          totalsum=sumcount1;
          for (i4=sumcount; i4<totalsum; i4++) {
            i1=1;
            for(sj=0; sj<sum_cof[i4].size; sj++) {
              i1=i1*sets[sum_cof[i4].setid[sj]].size;
            }
            sum_cof[i4].offset=sumbegadd;
            sumbegadd=sumbegadd+i1;
          }
          nsumele1=sumbegadd+nsumele;
          for (i3=sumcount; i3<totalsum; i3++) {
            i1=1;
            sum_cof[i3].strides[sum_cof[i3].size-1]=1;
            for(sj=sum_cof[i3].size-2; sj>-1; sj--) {
              sum_cof[i3].strides[sj]=sum_cof[i3].strides[sj+1]*sets[sum_cof[i3].setid[sj+1]].size;
            }
          }
          strcpy(line,leftline);
          readitem=line;
          stp->lv[i].sums= (sum_prog *) calloc (totalsum-sumcount+1,sizeof(sum_prog));
          while (sum_prog_build(readitem,sumsyntax,false,sum_cof,sumcount,sets,coefs,ncof,vars,nvar,ncofele,totalsum,ops,&sumindx,&stp->lv[i].sums[stp->lv[i].nsums])==1) {
            sumcount++;
            stp->lv[i].nsums++;
          }
          nops=0;
          if(!formula_compile(readitem,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdimlin))MPI_Abort(PETSC_COMM_WORLD,1);
          for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
            for (i4=0; i4<fdimlin; i4++) {
              if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].index_name)==0) {
                dcountdim3[dcount]=i4;
                break;
              }
            }
          }
          for(dcount=0; dcount<MAXSUPSET; dcount++)supset[dcount]=0;
          for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
            if(LinVars[i].dimmapid[dcount]>0) {
              linvar_map_dim_check(&LinVars[i],dcount,arSet[dcountdim3[dcount]].setid,vars);
              supset[dcount]=0;
              continue;
            }
            if(sets[vars[LinVars[i].LinVarIndx].setid[dcount]].size!=sets[arSet[dcountdim3[dcount]].setid].size) {
              for(sup=1; sup<MAXSUPSET; sup++)if(vars[LinVars[i].LinVarIndx].setid[dcount]==sets[arSet[dcountdim3[dcount]].setid].subsetid[sup]) {
                  supset[dcount]=sup;
                  break;
                }
            }
            else supset[dcount]=0;
          }
          stp->lv[i].ops= (formula_op *) malloc (nops*sizeof(formula_op));
          memcpy(stp->lv[i].ops,ops,nops*sizeof(formula_op));
          stp->lv[i].nops=nops;
          stp->lv[i].LinVarIndx=LinVars[i].LinVarIndx;
          stp->lv[i].fdimlin=fdimlin;
          stp->lv[i].nloopsfac=nloopsfac;
          stp->lv[i].arSet= (quantifier *) malloc (fdimlin*sizeof(quantifier));
          memcpy(stp->lv[i].arSet,arSet,fdimlin*sizeof(quantifier));
          for (dcount=0; dcount<fdimlin; dcount++) stp->lv[i].dcountdim2[dcount]=dcountdim2[dcount];
          for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
            stp->lv[i].dcountdim3[dcount]=dcountdim3[dcount];
            stp->lv[i].supset[dcount]=supset[dcount];
            stp->lv[i].dcountmap[dcount]=LinVars[i].dimmapid[dcount];
            stp->lv[i].dimleadlag[dcount]=LinVars[i].dimleadlag[dcount];
          }
          stp->lv[i].built=true;
        }
        free(sum_cof);
      }
      stp->nloops=nloops;

      free(ops);
      free(LinVars);
      free(arSet);
}

int jacobian_fill(char *fname, char *commsyntax,set_def *sets,offset_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,offset_t *eq_addr,offset_t *counteq,offset_t nintraeq,Mat A,Mat B) {
  FILE * filehandle;
  char line[TABREADLINE];
  PetscInt Istart1,Iend1,matrow;
  PetscErrorCode ierr;
  solve_real zerodivide=0;
  PetscMPIInt  mpisize1;
  stmt_prog *stp=NULL;
  offset_t eqindx=0;
  ierr = MatGetOwnershipRange(A,&Istart1,&Iend1);
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize1);
  CHKERRQ(ierr);

  filehandle = fopen(fname,"r");
  matrow=0;

  if (stmt_cache_built&&(Istart1!=stmt_cache_Istart||Iend1!=stmt_cache_Iend)) jacobian_cache_free();
  if (stmt_cache_built) {
    /* fast path: statements already compiled; walk the file only to
       refresh zerodivide defaults (they may reference scalar
       coefficients whose values change between steps) */
    int sidx=0;
    while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
      if (strstr(line,"(default")==NULL) {
        stmt_cache[sidx].zerodivide=zerodivide;
        stmt_prog_execute(&stmt_cache[sidx],matrow,eq_addr,stmt_cache[sidx].nloops,sets,set_elems,elem_vals,closure_vals,vars,Istart1,Iend1,A,B);
        matrow+=stmt_cache[sidx].nloops;
        sidx++;
      }
    }
    fclose(filehandle);
    return 1;
  }

  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      stmt_cache=realloc(stmt_cache,(stmt_cache_n+1)*sizeof(stmt_prog));
      stp=&stmt_cache[stmt_cache_n];
      stmt_cache_n++;
      memset(stp,0,sizeof(stmt_prog));
      stp->zerodivide=zerodivide;
      stmt_prog_build_one(line,stp,commsyntax,sets,nset,coefs,ncof,vars,nvar,ncofele,eq_addr,matrow,Istart1,Iend1,mpisize1,false);
      stmt_prog_execute(stp,matrow,eq_addr,stp->nloops,sets,set_elems,elem_vals,closure_vals,vars,Istart1,Iend1,A,B);
      matrow+=stp->nloops;
      eqindx++;
    }
  }
  stmt_cache_built=true;
  stmt_cache_Istart=Istart1;
  stmt_cache_Iend=Iend1;
  fclose(filehandle);
  return 1;
}

/* ===== backsolve recovery (GEMPACK manual 14.1.3, roadmap 6.2) ========
 *
 * A backsolved variable and its nominated defining equation are excluded
 * from the condensed system (the "equation" scan filter in tab_parse.c
 * hides the equation; main.c assigns the variable no matrix column).
 * After each step's solve, before the data updates, the retained defining
 * equation is evaluated row by row with the step's known changes:
 *
 *     sum over occurrences of  a_v * z_v  =  0
 *     =>  x_bs = -(sum over surviving occurrences) / a_pivot
 *
 * where z_v is x[column] for endogenous survivors and the captured vece
 * value for exogenous ones.  The rows ride the same compiled statement
 * programs as jacobian_fill (stmt_prog_build_one with force_all: no
 * ownership pruning -- recovery is a per-rank computation wherever the
 * update loops run).  teems-R's condensation guarantees each defining
 * equation references surviving variables plus exactly one plain
 * occurrence of its own backsolved variable, index-bijective with the
 * equation's domain; the build checks those properties and the execute
 * aborts on a zero pivot. */
typedef struct
{
  stmt_prog st;
  int pair;        /* index into backsolves[] */
  bool checked;    /* row->element bijection verified on first execute */
} bs_prog;
static bs_prog *bs_cache=NULL;
static int bs_cache_n=0;
static bool bs_cache_built=false;

void backsolve_cache_free(void) {
  int i;
  for (i=0; i<bs_cache_n; i++) stmt_prog_free(&bs_cache[i].st);
  free(bs_cache);
  bs_cache=NULL;
  bs_cache_n=0;
  bs_cache_built=false;
}

/* Format element li3 of variable v as "name(ele1,ele2,...)" for error
   messages. */
static void bs_element_label(array_def *vars, offset_t v, offset_t li3, set_def *sets, set_element *set_elems, char *out) {
  dim_t d;
  offset_t idx;
  strcpy(out,vars[v].cofname);
  if (vars[v].size==0) return;
  strcat(out,"(");
  for (d=0; d<vars[v].size; d++) {
    idx=li3/vars[v].strides[d];
    li3-=idx*vars[v].strides[d];
    strcat(out,set_elems[sets[vars[v].setid[d]].offset+idx].setele);
    if (d<vars[v].size-1) strcat(out,",");
  }
  strcat(out,")");
}

/* Evaluate one defining equation's rows and write the recovered per-step
   changes into bsvals at the pair's elem_base.  Aborts the run on a zero
   pivot or a defective row->element mapping: continuing would write a
   wrong solution. */
static void bs_prog_execute(bs_prog *bp, set_def *sets, set_element *set_elems,
                            elem_value *elem_vals, closure_entry *closure_vals,
                            array_def *vars, solve_real *x, solve_real *exo_z,
                            solve_real *bsvals) {
  stmt_prog *st=&bp->st;
  backsolve_def *bd=&backsolves[bp->pair];
  char label[TABREADLINE];
  int s,i;
  dim_t dcount;
  offset_t lj,l1,l2,li3,i5,gidx;
  solve_real vval;
  quantifier *arSet1=NULL;
  formula_op *ops1=NULL;
  int pivbad=0;
  solve_real *acc= (solve_real *) calloc (st->nloops,sizeof(solve_real));
  solve_real *piv= (solve_real *) calloc (st->nloops,sizeof(solve_real));
  offset_t *pivelem= (offset_t *) malloc (st->nloops*sizeof(offset_t));
  for (i5=0; i5<st->nloops; i5++) pivelem[i5]=-1;
  for (s=0; s<st->neqsums; s++) sum_prog_eval(&st->eqsums[s],sets,set_elems,elem_vals,st->sum_vals,st->zerodivide);
  for (i=0; i<st->nlv; i++) {
    linvar_prog *lv=&st->lv[i];
    if(!lv->built) continue;
    for (s=0; s<lv->nsums; s++) sum_prog_eval(&lv->sums[s],sets,set_elems,elem_vals,st->sum_vals,st->zerodivide);
    #pragma omp parallel private(lj,i5,l2,dcount,l1,li3,gidx,arSet1,ops1,vval) shared(elem_vals,st,lv,closure_vals,vars,acc,piv,pivelem,x,exo_z,bp,bd,pivbad)
    {
    if(omp_get_thread_num()!=0){
      arSet1=malloc(lv->fdimlin*sizeof(quantifier));
      memcpy(arSet1,lv->arSet,lv->fdimlin*sizeof(quantifier));
      ops1=malloc(lv->nops*sizeof(formula_op));
      memcpy(ops1,lv->ops,lv->nops*sizeof(formula_op));
    }else{
      ops1=lv->ops;
      arSet1=lv->arSet;
    }
    /* workers must finish copying the shared scratch (master aliases it)
       before the master starts mutating per-iteration state */
    #pragma omp barrier
    #pragma omp for
      for (i5=0; i5<st->nloops; i5++) {
        for (lj=i5*lv->nloopsfac; lj<(i5+1)*lv->nloopsfac; lj++) {
          l2=lj;
          for (dcount=0; dcount<lv->fdimlin; dcount++) {
            l1=(offset_t) l2/lv->dcountdim2[dcount];
            arSet1[dcount].indx=l1;
            l2=l2-l1*lv->dcountdim2[dcount];
          }
          li3=0;
          for (dcount=0; dcount<vars[lv->LinVarIndx].size; dcount++) {
            if(lv->dcountmap[dcount]>0) {
              li3=li3+((offset_t)teems_maps[lv->dcountmap[dcount]-1].values[arSet1[lv->dcountdim3[dcount]].indx])*vars[lv->LinVarIndx].strides[dcount];
            }
            else if(lv->supset[dcount]==0) {
              li3=li3+(arSet1[lv->dcountdim3[dcount]].indx+lv->dimleadlag[dcount])*vars[lv->LinVarIndx].strides[dcount];
            }
            else {
              li3=li3+(set_elems[sets[arSet1[lv->dcountdim3[dcount]].setid].offset+arSet1[lv->dcountdim3[dcount]].indx].superset_pos[lv->supset[dcount]]+lv->dimleadlag[dcount])*vars[lv->LinVarIndx].strides[dcount];
            }
          }
          vval=formula_eval(elem_vals,sets,set_elems,st->sum_vals,ops1,lv->nops,arSet1,lv->fdimlin,st->zerodivide);
          gidx=vars[lv->LinVarIndx].offset+li3;
          if (lv->LinVarIndx==bd->varindx) {
            /* pivot: same-pattern occurrences (possibly on both sides)
               combine into the algebraic pivot coefficient */
            piv[i5]+=vval;
            if (pivelem[i5]==-1) pivelem[i5]=li3;
            else if (pivelem[i5]!=li3) pivbad=1;
          }
          else if (closure_vals[gidx].is_exogenous) {
            acc[i5]+=vval*exo_z[gidx];
          }
          else {
            acc[i5]+=vval*x[closure_vals[gidx].exo_index];
          }
        }
      }
    if(omp_get_thread_num()!=0){
      free(arSet1);
      arSet1=NULL;
      free(ops1);
      ops1=NULL;
    }else{
      ops1=NULL;
      arSet1=NULL;
    }
    }
  }
  if (pivbad) {
    printf("Error: the occurrences of %s in defining equation %s map to different elements in the same row; the condensation must combine same-pattern terms before deployment (GEMPACK manual 14.1.10)\n",vars[bd->varindx].cofname,bd->eqname);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  if (!bp->checked) {
    /* one-time bijection check: every element of the backsolved variable
       must be recovered by exactly one row of the defining equation */
    char *seen= (char *) calloc (vars[bd->varindx].nelem,sizeof(char));
    for (i5=0; i5<st->nloops; i5++) {
      if (pivelem[i5]<0||pivelem[i5]>=vars[bd->varindx].nelem||seen[pivelem[i5]]) {
        printf("Error: defining equation %s does not map one-to-one onto backsolved variable %s; redeploy the model so the condensation validates this backsolve\n",bd->eqname,vars[bd->varindx].cofname);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      seen[pivelem[i5]]=1;
    }
    free(seen);
    bp->checked=true;
  }
  for (i5=0; i5<st->nloops; i5++) {
    if (piv[i5]==0) {
      bs_element_label(vars,bd->varindx,pivelem[i5],sets,set_elems,label);
      printf("Error: zero pivot backsolving %s from equation %s at this step; the defining equation cannot determine the variable here (GEMPACK would report the same singularity) -- nominate a different equation or leave the variable in the system\n",label,bd->eqname);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    bsvals[bd->elem_base+pivelem[i5]]=-acc[i5]/piv[i5];
  }
  free(acc);
  free(piv);
  free(pivelem);
}

int backsolve_recover(char *fname, char *commsyntax,set_def *sets,offset_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofele,closure_entry *closure_vals,solve_real *x,solve_real *exo_z,solve_real *bsvals) {
  FILE *filehandle;
  char line[TABREADLINE],eqname[NAMESIZE];
  solve_real zerodivide=0;
  int i,pr,sidx=0;
  offset_t d;
  bs_prog *bp=NULL;
  if (nbacksolve==0) return 0;
  backsolve_scan_mode=BS_SCAN_ONLY;
  filehandle=fopen(fname,"r");
  if (bs_cache_built) {
    /* refresh zerodivide defaults, as in the jacobian_fill fast path */
    while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
      if (strstr(line,"(default")==NULL) {
        bs_cache[sidx].st.zerodivide=zerodivide;
        sidx++;
      }
    }
  }
  else {
    while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
      if (strstr(line,"(default")!=NULL) continue;
      tab_equation_name(line,eqname);
      pr=-1;
      for (i=0; i<nbacksolve; i++) {
        if (strcmp(backsolves[i].eqname,eqname)==0) {
          pr=i;
          break;
        }
      }
      if (pr==-1) continue; /* cannot happen: the ONLY filter admits nominated equations */
      bs_cache=realloc(bs_cache,(bs_cache_n+1)*sizeof(bs_prog));
      bp=&bs_cache[bs_cache_n];
      bs_cache_n++;
      memset(bp,0,sizeof(bs_prog));
      bp->pair=pr;
      bp->st.zerodivide=zerodivide;
      stmt_prog_build_one(line,&bp->st,commsyntax,sets,nset,coefs,ncof,vars,nvar,ncofele,NULL,0,0,0,1,true);
      if (bp->st.nloops!=vars[backsolves[pr].varindx].nelem) {
        printf("Error: defining equation %s has %ld rows but backsolved variable %s has %ld elements; they must match one-to-one\n",backsolves[pr].eqname,bp->st.nloops,vars[backsolves[pr].varindx].cofname,vars[backsolves[pr].varindx].nelem);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      {
        int npivot=0;
        for (i=0; i<bp->st.nlv; i++) {
          if (bp->st.lv[i].LinVarIndx!=backsolves[pr].varindx) continue;
          npivot++;
          if (bp->st.lv[i].nloopsfac!=1) {
            printf("Error: an occurrence of %s in defining equation %s carries indices outside the equation's quantifiers; it cannot be backsolved from this equation (GEMPACK manual 14.1.10)\n",vars[backsolves[pr].varindx].cofname,backsolves[pr].eqname);
            MPI_Abort(PETSC_COMM_WORLD,1);
          }
          for (d=0; d<vars[backsolves[pr].varindx].size; d++) {
            if (bp->st.lv[i].dimleadlag[d]!=0) {
              printf("Error: an occurrence of %s in defining equation %s carries a lead/lag offset; it cannot be backsolved from this equation (GEMPACK manual 14.1.10)\n",vars[backsolves[pr].varindx].cofname,backsolves[pr].eqname);
              MPI_Abort(PETSC_COMM_WORLD,1);
            }
          }
        }
        if (npivot==0) {
          printf("Error: defining equation %s does not reference backsolved variable %s\n",backsolves[pr].eqname,vars[backsolves[pr].varindx].cofname);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
      }
    }
    if (bs_cache_n!=nbacksolve) {
      printf("Error: %d backsolve statements but %d nominated defining equations found in the TAB file\n",nbacksolve,bs_cache_n);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    bs_cache_built=true;
  }
  fclose(filehandle);
  backsolve_scan_mode=BS_SCAN_SKIP;
  for (i=0; i<bs_cache_n; i++) {
    bs_prog_execute(&bs_cache[i],sets,set_elems,elem_vals,closure_vals,vars,x,exo_z,bsvals);
  }
  return 1;
}

int eq_sum_parse(char *formulain, char *commsyntax, sum_def *sum_cof,quantifier *arSet,set_def *sets,dim_t nset,dim_t fdim,int j) {
  char *readitem,*p,*p1,*p2,interchar2[NAMESIZE],argu[TABREADLINE],tempname[NAMESIZE];//,line5[TABREADLINE]
  char interchar[NAMESIZE],interchar1[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],line4[TABREADLINE];
  dim_t l,l1,l2,l3,l4,l5,l6,l7;
  int ncur=0,ncuri,i=0,k=0,k1=0,length;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=str_find_ci(readitem,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      sum_extract(line);
      k1=str_find_ci(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
          strcpy(line2,line);
          strcpy(line4,line);
          sprintf(interchar1, "%d", j);
          interchar[0]='\0';
          if (j<10) {
            strcat(interchar,"gen_sum000");
          }
          if (9<j&&j<100) {
            strcat(interchar,"gen_sum00");
          }
          if (99<j&&j<1000) {
            strcat(interchar,"gen_sum0");
          }
          if (999<j&&j<10000) {
            strcat(interchar,"gen_sum");
          }
          if (i>10000) {
            strcat(interchar,"gen_sum");
            printf("Error: too many sum() terms in one statement\n");
          }
          strcat(interchar,interchar1);
          strcpy(sum_cof[j].sumname,interchar);
          strcat(interchar,"{");
          strcpy(interchar2,interchar);
          strcpy(line,line1);
          p = strtok(line,"(");
          p = strtok(NULL,",");
          strcpy(sum_cof[j].sumindx,p);
          p = strtok(NULL,",");
          for (l7=0; l7<nset; l7++) if(strcmp(p,sets[l7].setname)==0) {
              sum_cof[j].sumsetid=l7;
              break;
            }

          ncur=str_count_ci(line2, "{");
          l2=0;
          l3=0;
          for (ncuri=0; ncuri<ncur; ncuri++) {
            strcpy(line2,line4);
            p = &line2[l2];
            p = strtok(p,"{");
            p = strtok(NULL,"}");
            l2=p-line2;
            strcpy(argu,p);
            strcat(argu,",");
            l=str_count_ci(argu, ",");
            if (l<2) {
              if(strcmp(sum_dim_identity(p),sum_cof[j].sumindx)!=0) {
                strcat(interchar,sum_cof[j].sumindx);
                for (l4=0; l4<l3; l4++) {
                  if(strcmp(sum_dim_identity(p),sum_cof[j].dimnames[l4])==0) {
                    break;
                  }
                }
                if (l4==l3||l3==0) {
                  strcpy(sum_cof[j].dimnames[l3],sum_dim_identity(p));
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(sum_dim_identity(p),arSet[l5].index_name)==0) {
                      sum_cof[j].setid[l3]=l5;
                      l6++;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,sum_dim_identity(p));
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=str_rfind_ci(line3,interchar1);
                    p1=&line3[l7+2];
                    p1 = strtok(p1,",");
                    for (l7=0; l7<nset; l7++) if(strcmp(p1,sets[l7].setname)==0) {
                        sum_cof[j].setid[l3]=l7;
                        break;
                      }
                  }
                  l3++;
                }
              }
            }
            else {
              for (l1=0; l1<l; l1++) {
                if(l1==0) {
                  p = strtok(argu,",");
                }
                else {
                  p = strtok(NULL,",");
                }
                if(strcmp(sum_dim_identity(p),sum_cof[j].sumindx)!=0) {
                  for (l4=0; l4<l3; l4++) {
                    if(strcmp(sum_dim_identity(p),sum_cof[j].dimnames[l4])==0) {
                      break;
                    }
                  }
                  if (l4==l3||l3==0) {
                    strcpy(sum_cof[j].dimnames[l3],sum_dim_identity(p));
                    strcat(interchar,sum_cof[j].dimnames[l3]);
                    l6=0;
                    for (l5=0; l5<fdim-1; l5++) if(strcmp(sum_dim_identity(p),arSet[l5].index_name)==0) {
                        sum_cof[j].setid[l3]=l5;
                        l6++;
                        break;
                      }
                    if (l6==0) {
                      interchar1[0]='\0';
                      strcat(interchar1,"sum(");
                      strcat(interchar1,sum_dim_identity(p));
                      strcpy(line3,formulain);
                      line3[readitem-formulain]='\0';
                      l7=str_rfind_ci(line3,interchar1);
                      p1=&line3[l7+2];
                      p2=strchr(p1,',');
                      strncpy(tempname,p1,p2-p1);
                      tempname[p2-p1]='\0';
                      for (l7=0; l7<nset; l7++) if(strcmp(tempname,sets[l7].setname)==0) {
                          sum_cof[j].setid[l3]=l7;
                          break;
                        }
                    }
                    l3++;
                  }
                }
              }
            }
          }
          if (interchar[strlen(interchar)-1]==',') {
            interchar[strlen(interchar)-1]='}';
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(str_replace_all(formulain,line4,interchar));
          sum_cof[j].size=l3;
          return 1;
        }
      }

    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      sum_extract(line);
      k1=str_find_ci(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
          strcpy(line2,line);
          strcpy(line4,line);
          sprintf(interchar1, "%d", j);
          interchar[0]='\0';
          if (j<10) {
            strcat(interchar,"gen_sum000");
          }
          if (9<j&&j<100) {
            strcat(interchar,"gen_sum00");
          }
          if (99<j&&j<1000) {
            strcat(interchar,"gen_sum0");
          }
          if (999<j&&j<10000) {
            strcat(interchar,"gen_sum");
          }
          if (i>10000) {
            strcat(interchar,"gen_sum");
            printf("Error: too many sum() terms in one statement\n");
          }
          strcat(interchar,interchar1);
          strcpy(sum_cof[j].sumname,interchar);
          strcat(interchar,"{");
          strcpy(interchar2,interchar);
          strcpy(line,line1);
          p = strtok(line,"(");
          p = strtok(NULL,",");
          strcpy(sum_cof[j].sumindx,p);
          p = strtok(NULL,",");
          for (l7=0; l7<nset; l7++) if(strcmp(p,sets[l7].setname)==0) {
              sum_cof[j].sumsetid=l7;
              break;
            }

          ncur=str_count_ci(line2, "{");
          l2=0;
          l3=0;
          for (ncuri=0; ncuri<ncur; ncuri++) {
            strcpy(line2,line4);
            p = &line2[l2];
            p = strtok(p,"{");
            p = strtok(NULL,"}");
            l2=p-line2;
            strcpy(argu,p);
            strcat(argu,",");
            l=str_count_ci(argu, ",");
            if (l<2) {
              if(strcmp(sum_dim_identity(p),sum_cof[j].sumindx)!=0) {
                for (l4=0; l4<l3; l4++) {
                  if(strcmp(sum_dim_identity(p),sum_cof[j].dimnames[l4])==0) {
                    break;
                  }
                }
                if (l4==l3||l3==0) {
                  strcpy(sum_cof[j].dimnames[l3],sum_dim_identity(p));
                  strcat(interchar,sum_cof[j].dimnames[l3]);
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(sum_dim_identity(p),arSet[l5].index_name)==0) {
                      sum_cof[j].setid[l3]=l5;
                      l6++;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,sum_dim_identity(p));
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=str_rfind_ci(line3,interchar1);
                    p1=&line3[l7+2];
                    p1 = strtok(p1,",");
                    for (l7=0; l7<nset; l7++) if(strcmp(p1,sets[l7].setname)==0) {
                        sum_cof[j].setid[l3]=l7;
                        break;
                      }
                  }
                  l3++;
                }
              }
            }
            else {
              for (l1=0; l1<l; l1++) {
                if(l1==0) {
                  p = strtok(argu,",");
                }
                else {
                  p = strtok(NULL,",");
                }
                if(strcmp(sum_dim_identity(p),sum_cof[j].sumindx)!=0) {
                  for (l4=0; l4<l3; l4++) {
                    if(strcmp(sum_dim_identity(p),sum_cof[j].dimnames[l4])==0) {
                      break;
                    }
                  }
                  if (l4==l3||l3==0) {
                    strcpy(sum_cof[j].dimnames[l3],sum_dim_identity(p));
                    strcat(interchar,sum_cof[j].dimnames[l3]);
                    strcat(interchar,",");
                    l6=0;
                    for (l5=0; l5<fdim-1; l5++) if(strcmp(sum_dim_identity(p),arSet[l5].index_name)==0) {
                        sum_cof[j].setid[l3]=arSet[l5].setid;//l5;
                        l6++;
                      }
                    if (l6==0) {
                      interchar1[0]='\0';
                      strcat(interchar1,"sum(");
                      strcat(interchar1,sum_dim_identity(p));
                      strcpy(line3,formulain);
                      line3[readitem-formulain]='\0';
                      l7=str_rfind_ci(line3,interchar1);
                      p1=&line3[l7+2];
                      p2=strchr(p1,',');
                      strncpy(tempname,p1,p2-p1);
                      tempname[p2-p1]='\0';
                      for (l7=0; l7<nset; l7++) if(strcmp(tempname,sets[l7].setname)==0) {
                          sum_cof[j].setid[l3]=l7;
                          break;
                        }
                    }
                    l3++;
                  }
                }
              }
            }
          }
          if (interchar[strlen(interchar)-1]==',') {
            interchar[strlen(interchar)-1]='}';
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(str_replace_all(formulain,line4,interchar));
          sum_cof[j].size=l3;
          return 1;
        }
      }

    }
    else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

int eq_sum_replace(char *formulain, char *commsyntax,int LinIndx, eq_var_ref *LinVars,array_def *vars) {
  char *readitem,line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],argu[TABREADLINE],*p0,*p1;
  int i=0,k=0,k1=0,k2,k3,length,tsum,j,l;
  bool IsRemSum;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=str_find_ci(formulain+i,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      sum_extract(line);
      strcpy(line3,line);
      if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
        tsum=sum_count(line,commsyntax);
        IsRemSum=false;
        for (j=tsum; j>1; j--) {
          k2=1;
          k1=0;
          while(k2<j) {
            k3=k1;
            k1=str_find_ci(line+k1+3,commsyntax);
            k1=k3+k1+3;//printf("k1 %d line %s\n",k1,line+k1-1);
            if (line[k1-1]=='+'||line[k1-1]=='-'||line[k1-1]=='*'||line[k1-1]=='/'||line[k1-1]=='^'||line[k1-1]=='('||line[k1-1]==',') k2++;//printf("k1 %d line %s\n",k1,line+k1-1);}
          }
          strcpy(line1,&line[k1]);
          sum_extract(line1);
          p0=strchr(line1,'(');
          p0++;
          p1=strchr(line1,',');
          strncpy(argu,p0,p1-p0);
          argu[p1-p0]='\0';
          for(l=0; l<vars[LinVars[LinIndx].LinVarIndx].size; l++) {
            if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
              p1++;
              p1=strchr(p1,',');
              p1++;
              strncpy(line2,p1,strlen(line1)-(p1-line1)-1);
              line2[strlen(line1)-(p1-line1)-1]='\0';
              str_replace_all(line,line1,line2);
              IsRemSum=true;
              break;
            }
          }
        }
        p0=strchr(line,'(');
        p0++;
        p1=strchr(line,',');
        strncpy(argu,p0,p1-p0);
        argu[p1-p0]='\0';
        for(l=0; l<vars[LinVars[LinIndx].LinVarIndx].size; l++) {
          if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
            p1++;
            p1=strchr(p1,',');
            p1++;
            line2[0]='(';
            strncpy(&line2[1],p1,strlen(line)-(p1-line)-1);
            line2[strlen(line)-(p1-line)]=')';
            line2[strlen(line)-(p1-line)+1]='\0';
            strcpy(line,line2);
            IsRemSum=true;
            break;
          }
        }
        if(IsRemSum) {
          str_replace_all(formulain,line3,line);
          i=i+k+strlen(line);
        }
        else i=i+k+strlen(line3);
      }
      else {
        str_replace_all(formulain,line,"0");
        i=i+k+1;
      }
    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      sum_extract(line);
      strcpy(line3,line);
      if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
        tsum=sum_count(line,commsyntax);
        IsRemSum=false;
        for (j=tsum; j>1; j--) {
          k2=1;
          k1=0;
          while(k2<j) {
            k3=k1;
            k1=str_find_ci(line+k1+3,commsyntax);
            k1=k3+k1+3;//printf("k1 %d line %s\n",k1,line+k1-1);
            if (line[k1-1]=='+'||line[k1-1]=='-'||line[k1-1]=='*'||line[k1-1]=='/'||line[k1-1]=='^'||line[k1-1]=='('||line[k1-1]==',') k2++;//printf("k1 %d line %s\n",k1,line+k1-1);}
          }
          strcpy(line1,&line[k1]);
          sum_extract(line1);
          p0=strchr(line1,'(');
          p0++;
          p1=strchr(line1,',');
          strncpy(argu,p0,p1-p0);
          argu[p1-p0]='\0';
          for(l=0; l<vars[LinVars[LinIndx].LinVarIndx].size; l++) {
            if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
              p1++;
              p1=strchr(p1,',');
              p1++;
              strncpy(line2,p1,strlen(line1)-(p1-line1)-1);
              line2[strlen(line1)-(p1-line1)-1]='\0';
              str_replace_all(line,line1,line2);
              IsRemSum=true;
              break;
            }
          }
        }
        p0=strchr(line,'(');
        p0++;
        p1=strchr(line,',');
        strncpy(argu,p0,p1-p0);
        argu[p1-p0]='\0';
        for(l=0; l<vars[LinVars[LinIndx].LinVarIndx].size; l++) {
          if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
            p1++;
            p1=strchr(p1,',');
            p1++;
            line2[0]='(';
            strncpy(&line2[1],p1,strlen(line)-(p1-line)-1);
            line2[strlen(line)-(p1-line)]=')';
            line2[strlen(line)-(p1-line)+1]='\0';
            strcpy(line,line2);
            IsRemSum=true;
            break;
          }
        }
        if(IsRemSum) {
          str_replace_all(formulain,line3,line);
          i=i+k+strlen(line);
        }
        else i=i+k+strlen(line3);
      }
      else {
        str_replace_all(formulain,line,"0");
        i=i+k+1;
      }
    }
    else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

int eq_linvar_read(char *formulain,eq_var_ref *LinVars,int linindx,array_def *vars) {
  char line[TABREADLINE];
  int i1,i2,i3,l;//,l,d;
  line[0]='p';
  line[1]='_';
  line[2]='\0';
  l=strlen(formulain);
  strcat(line,LinVars[linindx].LinVarName);
  i3=0;
  i1=0;
  while (i3==0) {
    i1=str_find_ci(formulain+i1,line);
    if (i1==-1) break;
    if (i1==0) {
      if(vars[LinVars[linindx].LinVarIndx].size==0) {
        i2=l-1;
      }
      else
        i2=str_find_ci(formulain,"}");
      formulain[0]='1';
      memmove(formulain+1,formulain+i2+1,l-i2);
      i3=1;
    }
    else {
      if (formulain[i1-1]=='+'||formulain[i1-1]=='-'||formulain[i1-1]=='*'||formulain[i1-1]=='/'||formulain[i1-1]=='^'||formulain[i1-1]=='('||formulain[i1-1]==',') {
        if(vars[LinVars[linindx].LinVarIndx].size==0) {
          i2=strlen(line)-1;
        }
        else i2=str_find_ci(formulain+i1,"}");
        formulain[i1]='1';
        memmove(formulain+i1+1,formulain+i1+i2+1,l-i2-i1);
        i3=1;
      }
      else i1++;
    }
  }
  return 1;
}

/* Element-level border marking for one variable reference (6.5 E3).
   Called for references that are not block-safe on the chain dimension
   (a lead/lag offset, or an index other than the equation's chain
   quantifier). Each dim narrows the mark to the positions its index
   set occupies inside the variable's declared set, shifted by that
   dim's lead/lag offset and clipped to the declared range
   (tab_preprocess turns quoted elements into singleton subsets, so
   fixed dims arrive here as 1-element sets); unresolvable dims stay
   full. A mark spanning the whole variable collapses to the
   whole-variable flag. */
static void border_mark_ref(eq_var_ref *ref, array_def *vars, set_def *sets, dim_t nset,
                            set_element *set_elems, offset_t alltimeset,
                            bool *var_inter, bool *ele_inter, bool nested,
                            dim_t *orderintra, dim_t *orderreg) {
  offset_t l=ref->LinVarIndx;
  offset_t poscnt[MAXVARDIM],mstr[MAXVARDIM],*poslist[MAXVARDIM];
  offset_t s,li,l2,i1,addr,total,dsize,pos;
  dim_t nd=vars[l].size,d,sup;
  int k;
  bool allfull=true;
  if(var_inter[l])return;
  if(alltimeset<0||orderintra[l]<0||nd==0) {
    var_inter[l]=true;
    if(nested) {
      orderintra[l]=-1;
      orderreg[l]=-1;
    }
    return;
  }
  for(d=0; d<nd; d++) {
    poslist[d]=NULL;
    dsize=sets[vars[l].setid[d]].size;
    poscnt[d]=dsize;
    k=ref->dimleadlag[d];
    for(s=0; s<nset; s++)if(strcmp(sets[s].setname,ref->dimsetnames[d])==0)break;
    if(s==nset||s==vars[l].setid[d]||sets[s].size==sets[vars[l].setid[d]].size) {
      /* index runs over the declared set itself (or cannot be
         resolved): full range unless shifted */
      if(k==0)continue;
      poslist[d]=(offset_t *) malloc (dsize*sizeof(offset_t));
      poscnt[d]=0;
      for(li=0; li<dsize; li++) {
        pos=li+k;
        if(pos>=0&&pos<dsize)poslist[d][poscnt[d]++]=pos;
      }
    }
    else {
      for(sup=1; sup<MAXSUPSET; sup++)if(sets[s].subsetid[sup]==vars[l].setid[d])break;
      if(sup==MAXSUPSET) { /* not a resolvable subset of the declared set */
        if(k==0)continue;
        poslist[d]=(offset_t *) malloc (dsize*sizeof(offset_t));
        poscnt[d]=0;
        for(li=0; li<dsize; li++) {
          pos=li+k;
          if(pos>=0&&pos<dsize)poslist[d][poscnt[d]++]=pos;
        }
      }
      else {
        poslist[d]=(offset_t *) malloc (sets[s].size*sizeof(offset_t));
        poscnt[d]=0;
        for(li=0; li<sets[s].size; li++) {
          pos=set_elems[sets[s].offset+li].superset_pos[sup]+k;
          if(pos>=0&&pos<dsize)poslist[d][poscnt[d]++]=pos;
        }
      }
    }
    if(poscnt[d]<dsize)allfull=false;
    if(poscnt[d]==0) { /* every position clipped out: nothing referenced */
      for(li=0; li<=d; li++)free(poslist[li]);
      return;
    }
  }
  if(allfull) {
    for(d=0; d<nd; d++)free(poslist[d]);
    var_inter[l]=true;
    if(nested) {
      orderintra[l]=-1;
      orderreg[l]=-1;
    }
    return;
  }
  mstr[nd-1]=1;
  for(d=nd-1; d>0; d--)mstr[d-1]=mstr[d]*poscnt[d];
  total=mstr[0]*poscnt[0];
  for(li=0; li<total; li++) {
    l2=li;
    addr=0;
    for(d=0; d<nd; d++) {
      i1=l2/mstr[d];
      l2-=i1*mstr[d];
      addr+=(poslist[d]==NULL?i1:poslist[d][i1])*vars[l].strides[d];
    }
    ele_inter[vars[l].offset+addr]=true;
  }
  for(d=0; d<nd; d++)free(poslist[d]);
}

/* The equation's chain quantifier: index name of the last (all,...)
   quantifier over an intertemporal set — the same selection
   equation_order_read[_nested] use for eq_orderintra — or NULL when the
   equation has none. */
static const char *eq_chain_index(quantifier *arSet, dim_t fdim, set_def *sets) {
  const char *q=NULL;
  dim_t i;
  for(i=0; i<fdim; i++)if(sets[arSet[i].setid].intertemp)q=arSet[i].index_name;
  return q;
}

int equation_order_read(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,bool *var_inter,bool *ele_inter,array_def *eq_defs,bool *eq_intertemp,dim_t *eq_orderintra,dim_t *eq_orderreg,offset_t allregset,offset_t alltimeset,dim_t *orderintra,dim_t *orderreg) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  solve_real zerodivide=0;
  dim_t fdim,np,i4;
  offset_t rowindx,j,l,l1,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],nloops,nloopslin,li3,l2,matrow,eqindx=0,nelem,l01,j01;//,sizelinvars,totlinvars,templinvars
  int sj,i,i3,nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,leadlag,varindx1,varindx2,tempint;
  solve_real dimmat[4*MAXVARDIM];

  filehandle = fopen(fname,"r");
  if (filehandle==NULL) return 0;
  matrow=0;

  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      str_replace_first(line, commsyntax, "");
      str_replace_first(line, "(linear)", "");
      while (str_replace_all(line,"  ", " "));
      while (str_replace_char(line, '[', '('));
      while (str_replace_char(line, ']', ')'));
      while (str_replace_char(line, '{', '('));
      while (str_replace_char(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=str_count_ci(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        if (readitem==NULL||strlen(readitem)>=sizeof(eq_defs[eqindx].cofname)) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcpy(eq_defs[eqindx].cofname,readitem);
        readitem = strtok(NULL,"=");
        if (readitem==NULL) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        if (readitem==NULL||(readitem-line)+strlen(readitem)+strlen(vname)+4>=sizeof(line)) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        readitem = strtok(line+1,"(");
        if (readitem==NULL||strlen(readitem)>=sizeof(eq_defs[eqindx].cofname)) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcpy(eq_defs[eqindx].cofname,readitem);
        strcpy(line,linecopy);
        i=str_rfind_ci(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        if (readitem==NULL) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        if (readitem==NULL||(readitem-line)+strlen(readitem)+strlen(vname)+4>=sizeof(line)) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (str_replace_all(readitem," ", ""));
      while (formula_normalize(readitem)==1);
      leadlag_encode(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=str_count_ci(readitem,"p_");
      eq_var_ref *LinVars= (eq_var_ref *) calloc (np+1,sizeof(eq_var_ref));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=str_find_ci(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          if (l>=nvar) break;
          p=strpbrk(readitem,"}+*-/^)");
          if (p==NULL) break;
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (vars[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            linvar_dim_read(p,linecopy,lvar,&LinVars[i3],0,false,sets);
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<vars[l].size-1; i4++) {
              p = strtok(NULL,",");
              linvar_dim_read(p,linecopy,lvar,&LinVars[i3],i4,false,sets);
            }
            p = strtok(NULL,"}");
            linvar_dim_read(p,linecopy,lvar,&LinVars[i3],i4,false,sets);
            break;
          }
          i3++;
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=vars[LinVars[i].LinVarIndx].size;
      quantifier *arSet= (quantifier *) calloc (fdim+i3+1,sizeof(quantifier));
      bool *linvarrcount= (bool *) calloc (nlinvars,sizeof(bool));

      strcpy(line,linecopy);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].index_name,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
              arSet[i].setid=i4;
              eq_defs[eqindx].setid[i]=i4;
              break;
            }
          nloops=nloops*sets[arSet[i].setid].size;
        }
      }
      eq_defs[eqindx].size=fdim;
      eq_defs[eqindx].nelem=nloops;
      if(alltimeset>=0) {
        if(eqindx==0) {
          for(i=0; i<nvar; i++) {
            if(var_inter[i]) continue;
            j=0;
            for(i4=0; i4<vars[i].size; i4++)if(sets[vars[i].setid[i4]].intertemp) {
                orderintra[i]=i4;
                j++;
                break;
              }
            if(j==0)var_inter[i]=true;
          }
        }
        j=0;
        for (i=0; i<fdim; i++)if(sets[eq_defs[eqindx].setid[i]].intertemp) {
            j++;
            eq_orderintra[eqindx]=i;
          }
        if(j==0) {
          if(allregset<0)for (i4=0; i4<nlinvars; i4++) var_inter[LinVars[i4].LinVarIndx]=true;
          eq_intertemp[eqindx]=true;
        }

      }
      /* 6.5 E3: references that are not block-safe on the chain
         dimension — a lead/lag offset, or (when the equation has a
         chain quantifier) a chain-dim index other than that quantifier
         (sum indices, fixed elements, other quantifiers) — border at
         element level; must run before the partition rules below,
         which read var_inter */
      {
        const char *eqchain=alltimeset>=0?eq_chain_index(arSet,fdim,sets):NULL;
        for (i4=0; i4<nlinvars; i4++) {
          l=LinVars[i4].LinVarIndx;
          for(i=0; i<vars[l].size; i++)if(LinVars[i4].dimleadlag[i]!=0)break;
          if(i==vars[l].size) {
            if(eqchain==NULL||var_inter[l]||orderintra[l]<0)continue;
            if(strcmp(LinVars[i4].dimnames[orderintra[l]],eqchain)==0)continue;
          }
          border_mark_ref(&LinVars[i4],vars,sets,nset,set_elems,alltimeset,var_inter,ele_inter,false,orderintra,orderreg);
        }
      }
      if(allregset>=0) {
        if(eqindx==0) {
          for(i=0; i<nvar; i++) {
            j=0;
            for(i4=0; i4<vars[i].size; i4++)if(sets[vars[i].setid[i4]].regional) {
                orderreg[i]=i4;
                j++;
              }
            if(j==0)var_inter[i]=true;
          }
        }
        j=0;
        for(i4=0; i4<nlinvars; i4++)for(i=0; i<vars[LinVars[i4].LinVarIndx].size; i++) if(sets[vars[LinVars[i4].LinVarIndx].setid[i]].regional)j++;
        if(j==0)eq_intertemp[eqindx]=true;
        j=0;
        for (i=0; i<fdim; i++)if(sets[eq_defs[eqindx].setid[i]].regional) {
            j++;
            eq_orderreg[eqindx]=i;
          }
        if(j==0) {
          eq_intertemp[eqindx]=true;
        }
        if(j>=2) {
          tempint=eq_orderreg[eqindx];
          eq_orderreg[eqindx]=-1;
          for (i=0; i<fdim; i++)dimmat[i]=-1;
          for (i=0; i<fdim; i++) {
            memset(linvarrcount,0,nlinvars*sizeof(bool));
            for(i4=0; i4<nlinvars; i4++) {
              for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;//{printf("here!!! linvarrcount[lvar4] %d\n",linvarrcount[lvar4]);
              if(lvar4<i4&&linvarrcount[lvar4]) {
                linvarrcount[i4]=true;
                continue;
              }
              for(i3=0; i3<vars[LinVars[i4].LinVarIndx].size; i3++) {
                if(strcmp(LinVars[i4].dimnames[i3],arSet[i].index_name)==0&&sets[eq_defs[eqindx].setid[i]].regional&&orderreg[LinVars[i4].LinVarIndx]==i3) { //ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i3]].regional) {
                  dimmat[i]+=vars[LinVars[i4].LinVarIndx].nelem;
                  l01=0;
                  for(j01=0; j01<vars[LinVars[i4].LinVarIndx].nelem; j01++)if(closure_vals[vars[LinVars[i4].LinVarIndx].offset+j01].is_exogenous)l01++;
                  dimmat[i]-=l01;
                  linvarrcount[i4]=true;
                  if(lvar4<i4)linvarrcount[lvar4]=true;
                }
              }
            }
          }
          nelem=0;
          for (i=0; i<fdim; i++)if(dimmat[i]>nelem)nelem=dimmat[i];
          for (i=0; i<fdim; i++)if(dimmat[i]==nelem)break;
          if(nelem==0) {
            eq_orderreg[eqindx]=tempint;
            eq_intertemp[eqindx]=true;
          }
          else {
            eq_orderreg[eqindx]=i;
          }

        }


        if(eq_orderreg[eqindx]>-1) {
          nelem=0;
          memset(linvarrcount,0,nlinvars*sizeof(bool));
          for(i4=0; i4<nlinvars; i4++) {
            for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
            if(lvar4<i4&&linvarrcount[lvar4]) {
              linvarrcount[i4]=true;
              continue;
            }
            if((orderreg[LinVars[i4].LinVarIndx]<0||strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].index_name)!=0)&&!var_inter[LinVars[i4].LinVarIndx]) {
              nelem+=vars[LinVars[i4].LinVarIndx].nelem;
              l01=0;
              for(j01=0; j01<vars[LinVars[i4].LinVarIndx].nelem; j01++)if(closure_vals[vars[LinVars[i4].LinVarIndx].offset+j01].is_exogenous)l01++;
              nelem+=vars[LinVars[i4].LinVarIndx].nelem-l01;
              linvarrcount[i4]=true;
              if(lvar4<i4)linvarrcount[lvar4]=true;
            }
          }
          if(nelem>=eq_defs[eqindx].nelem)eq_intertemp[eqindx]=true;
          else {
            for(i4=0; i4<nlinvars; i4++)
              if((orderreg[LinVars[i4].LinVarIndx]<0||strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].index_name)!=0)&&!var_inter[LinVars[i4].LinVarIndx]) {
                var_inter[LinVars[i4].LinVarIndx]=true;
              }
          }
        }
      }
      free(arSet);
      free(linvarrcount);
      free(LinVars);
      eqindx++;
    }
  }
  fclose(filehandle);
  return 1;
}

int equation_order_read_nested(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,bool *var_inter,bool *ele_inter,array_def *eq_defs,bool *eq_intertemp,dim_t *eq_orderintra,dim_t *eq_orderreg,offset_t allregset,offset_t alltimeset,dim_t *orderintra,dim_t *orderreg) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  solve_real zerodivide=0;
  dim_t fdim,np;
  offset_t rowindx,j,j01,l,l01,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],nloops,nloopslin,li3,l2,matrow,eqindx=0,nelem;//,sizelinvars,totlinvars,templinvars
  offset_t sj,i,i3,i4,nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,varindx1,varindx2,tempint;
  int leadlag;
  solve_real dimmat[4*MAXVARDIM];

  filehandle = fopen(fname,"r");
  if (filehandle==NULL) return 0;
  matrow=0;

  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      str_replace_first(line, commsyntax, "");
      str_replace_first(line, "(linear)", "");
      while (str_replace_all(line,"  ", " "));
      while (str_replace_char(line, '[', '('));
      while (str_replace_char(line, ']', ')'));
      while (str_replace_char(line, '{', '('));
      while (str_replace_char(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=str_count_ci(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        if (readitem==NULL||strlen(readitem)>=sizeof(eq_defs[eqindx].cofname)) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcpy(eq_defs[eqindx].cofname,readitem);
        readitem = strtok(NULL,"=");
        if (readitem==NULL) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        if (readitem==NULL||(readitem-line)+strlen(readitem)+strlen(vname)+4>=sizeof(line)) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        readitem = strtok(line+1,"(");
        if (readitem==NULL||strlen(readitem)>=sizeof(eq_defs[eqindx].cofname)) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcpy(eq_defs[eqindx].cofname,readitem);
        strcpy(line,linecopy);
        i=str_rfind_ci(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        if (readitem==NULL) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        if (readitem==NULL||(readitem-line)+strlen(readitem)+strlen(vname)+4>=sizeof(line)) {
          printf("Error: malformed %s statement in TAB file: %s\n",commsyntax,linecopy);
          fclose(filehandle);
          return 0;
        }
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (str_replace_all(readitem," ", ""));
      while (formula_normalize(readitem)==1);
      leadlag_encode(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=str_count_ci(readitem,"p_");
      eq_var_ref *LinVars= (eq_var_ref *) calloc (np+1,sizeof(eq_var_ref));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=str_find_ci(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          if (l>=nvar) break;
          p=strpbrk(readitem,"}+*-/^)");
          if (p==NULL) break;
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (vars[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            linvar_dim_read(p,linecopy,lvar,&LinVars[i3],0,false,sets);
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<vars[l].size-1; i4++) {
              p = strtok(NULL,",");
              linvar_dim_read(p,linecopy,lvar,&LinVars[i3],i4,false,sets);
            }
            p = strtok(NULL,"}");
            linvar_dim_read(p,linecopy,lvar,&LinVars[i3],i4,false,sets);
            break;
          }
          i3++;
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=vars[LinVars[i].LinVarIndx].size;
      quantifier *arSet= (quantifier *) calloc (fdim+i3+1,sizeof(quantifier));
      bool *linvarrcount= (bool *) calloc (nlinvars,sizeof(bool));

      strcpy(line,linecopy);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].index_name,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
              arSet[i].setid=i4;
              eq_defs[eqindx].setid[i]=i4;
              break;
            }
          nloops=nloops*sets[arSet[i].setid].size;
        }
      }
      eq_defs[eqindx].size=fdim;
      eq_defs[eqindx].nelem=nloops;


      if(eqindx==0) {
        for(i=0; i<nvar; i++) {
          if(var_inter[i]) continue;
          j=0;
          for(i4=0; i4<vars[i].size; i4++)if(sets[vars[i].setid[i4]].intertemp) {
              orderintra[i]=i4;
              j++;
              break;
            }
          if(j==0) {
            var_inter[i]=true;
            orderreg[i]=-1;
            orderintra[i]=-1;
          }
          else {
            j01=0;
            for(i4=0; i4<vars[i].size; i4++)if(sets[vars[i].setid[i4]].regional) {
                orderreg[i]=i4;
                j01++;
              }
          }
        }
      }

      /* 6.5 E3: references that are not block-safe on the chain
         dimension border at element level (see the non-nested reader);
         must run before the rules below, which read var_inter (a
         collapse to the whole-variable flag also resets the order dims,
         as the inline flagging used to) */
      {
        const char *eqchain=alltimeset>=0?eq_chain_index(arSet,fdim,sets):NULL;
        for (i4=0; i4<nlinvars; i4++) {
          l=LinVars[i4].LinVarIndx;
          for(i=0; i<vars[l].size; i++)if(LinVars[i4].dimleadlag[i]!=0)break;
          if(i==vars[l].size) {
            if(eqchain==NULL||var_inter[l]||orderintra[l]<0)continue;
            if(strcmp(LinVars[i4].dimnames[orderintra[l]],eqchain)==0)continue;
          }
          border_mark_ref(&LinVars[i4],vars,sets,nset,set_elems,alltimeset,var_inter,ele_inter,true,orderintra,orderreg);
        }
      }

      j=0;
      for(i4=0; i4<nlinvars; i4++)for(i=0; i<vars[LinVars[i4].LinVarIndx].size; i++) if(sets[vars[LinVars[i4].LinVarIndx].setid[i]].intertemp)j++;
      if(j==0) {
        eq_intertemp[eqindx]=true;//eq_orderreg[eqindx]=-1;eq_orderintra[eqindx]=-1;}
      }
      j=0;
      for (i=0; i<fdim; i++)if(sets[eq_defs[eqindx].setid[i]].regional) {
          j++;
          eq_orderreg[eqindx]=i;
        }
      if(j>=2) {
        for (i=0; i<fdim; i++)dimmat[i]=-1;
        for (i=0; i<fdim; i++) {
          memset(linvarrcount,0,nlinvars*sizeof(bool));
          for(i4=0; i4<nlinvars; i4++) {
            for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
            if(lvar4<i4&&linvarrcount[lvar4]) {
              linvarrcount[i4]=true;
              continue;
            }
            for(i3=0; i3<vars[LinVars[i4].LinVarIndx].size; i3++) {
              if(strcmp(LinVars[i4].dimnames[i3],arSet[i].index_name)==0&&sets[eq_defs[eqindx].setid[i]].regional&&orderreg[LinVars[i4].LinVarIndx]==i3) {
                dimmat[i]+=vars[LinVars[i4].LinVarIndx].nelem;
                linvarrcount[i4]=true;
                if(lvar4<i4)linvarrcount[lvar4]=true;
                l01=0;
                for(j01=0; j01<vars[LinVars[i4].LinVarIndx].nelem; j01++)if(closure_vals[vars[LinVars[i4].LinVarIndx].offset+j01].is_exogenous)l01++;
                dimmat[i]-=l01;
              }
            }
          }
        }
        nelem=0;
        for (i=0; i<fdim; i++)if(dimmat[i]>nelem)nelem=dimmat[i];
        for (i=0; i<fdim; i++)if(dimmat[i]==nelem)break;
        if(nelem==0) {
          eq_orderreg[eqindx]=-1;
        }
        else {
          eq_orderreg[eqindx]=i;
        }
      }
      if(eq_orderreg[eqindx]>-1) {
        nelem=0;
        memset(linvarrcount,0,nlinvars*sizeof(bool));
        for(i4=0; i4<nlinvars; i4++) {
          for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
          if(lvar4<i4&&linvarrcount[lvar4]) {
            linvarrcount[i4]=true;
            continue;
          }
          if(orderreg[LinVars[i4].LinVarIndx]!=-1&&strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].index_name)!=0) {
            l01=0;
            for(j01=0; j01<vars[LinVars[i4].LinVarIndx].nelem; j01++)if(closure_vals[vars[LinVars[i4].LinVarIndx].offset+j01].is_exogenous)l01++;
            nelem+=vars[LinVars[i4].LinVarIndx].nelem-l01;
            linvarrcount[i4]=true;
            if(lvar4<i4)linvarrcount[lvar4]=true;
          }
        }
        if(nelem>=eq_defs[eqindx].nelem)eq_orderreg[eqindx]=-1;//printf("OOOOOOOOOOOOOOOOO\n");}
        else {
          for(i4=0; i4<nlinvars; i4++)
            if(orderreg[LinVars[i4].LinVarIndx]!=-1&&strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].index_name)!=0) {
              orderreg[LinVars[i4].LinVarIndx]=-1;//orderintra[LinVars[i4].LinVarIndx]=-1;
            }
        }
      }


      j=0;
      for (i=0; i<fdim; i++)if(sets[eq_defs[eqindx].setid[i]].intertemp) {
          j++;
          eq_orderintra[eqindx]=i;
        }
      if(j==0) {
        eq_intertemp[eqindx]=true;
      }

      if(eq_intertemp[eqindx]==true) {
        eq_orderreg[eqindx]=-1;
        eq_orderintra[eqindx]=-1;
      }

      free(arSet);
      free(LinVars);
      free(linvarrcount);
      eqindx++;
    }
  }
  fclose(filehandle);
  return 1;
}

int jacobian_preallocate(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele, offset_t nexo,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,bool *eq_intertemp,offset_t *eq_addr,dim_t *eq_time,dim_t *eq_reg,offset_t *counteq,offset_t nintraeq,bool *sbbd_overrid,PetscInt Istart,PetscInt Iend,PetscInt *dnz,PetscInt *dnnz,PetscInt *onz,PetscInt *onnz,PetscInt *dnzB,PetscInt *dnnzB,PetscInt *onzB,PetscInt *onnzB,int nesteddbbd,eq_probe_meta *eqmeta,offset_t *neqmeta) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  PetscInt Iindx=0,Jindx;
  solve_real zerodivide=0;
  dim_t fdim=0,np,i4,sup,supset[MAXSUPSET];
  offset_t rowindx,rowindxorg,l,l1,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],dcountdim3[4*MAXVARDIM],dcountdim4[4*MAXVARDIM],dcountdim5[4*MAXVARDIM],nloops,nloopslin,nloopsfac,li3,l2,matrow,matroworg,ltime,lreg,leq=0,eqindx=0;//,sizelinvars,totlinvars,templinvars
  offset_t nreg=0,nint=0,sj,i,i3;
  if(allregset>-1)nreg=sets[allregset].size;
  if(alltimeset>-1)nint=sets[alltimeset].size;
  int nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,leadlag,varindx1,varindx2;
  offset_t *counteq1= (offset_t *) calloc (ndblock,sizeof(offset_t));
  for(i=0; i<ndblock; i++)counteq1[i]=counteq[i];
  filehandle = fopen(fname,"r");
  matrow=nintraeq;
  matroworg=0;
  Jindx=0;

  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      if(eqmeta!=NULL)tab_equation_name(line,eqmeta[eqindx].eqname);
      str_replace_first(line, commsyntax, "");
      str_replace_first(line, "(linear)", "");
      while (str_replace_all(line,"  ", " "));
      while (str_replace_char(line, '[', '('));
      while (str_replace_char(line, ']', ')'));
      while (str_replace_char(line, '{', '('));
      while (str_replace_char(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=str_count_ci(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        i=str_rfind_ci(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (str_replace_all(readitem," ", ""));
      while (formula_normalize(readitem)==1);
      leadlag_encode(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=str_count_ci(readitem,"p_");
      eq_var_ref *LinVars= (eq_var_ref *) calloc (np+1,sizeof(eq_var_ref));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=str_find_ci(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          if (l>=nvar) break;
          p=strpbrk(readitem,"}+*-/^)");
          if (p==NULL) break;
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (vars[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            linvar_dim_read(p,linecopy,lvar,&LinVars[i3],0,true,sets);
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<vars[l].size-1; i4++) {
              p = strtok(NULL,",");
              linvar_dim_read(p,linecopy,lvar,&LinVars[i3],i4,true,sets);
            }
            p = strtok(NULL,"}");
            linvar_dim_read(p,linecopy,lvar,&LinVars[i3],i4,true,sets);
            break;
          }
          i3++;
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=vars[LinVars[i].LinVarIndx].size;
      quantifier *arSet= (quantifier *) calloc (fdim+i3+1,sizeof(quantifier));

      strcpy(line,linecopy);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].index_name,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
              arSet[i].setid=i4;
              if(sets[i4].intertemp) *sbbd_overrid=true;
              break;
            }
          nloops=nloops*sets[arSet[i].setid].size;
        }
        //*****Adjust ha_eqreg*******************//
        //*****End Adjust ha_eqreg*******************//
        dcountdim1[fdim-1]=1;
        for (dcount=fdim-2; dcount>-1; dcount--) {
          dcountdim1[dcount]=sets[arSet[dcount+1].setid].size*dcountdim1[dcount+1];
        }
        if(alltimeset>=0&&allregset<0) {
          for (dcount=eq_time[eqindx]-1; dcount>-1; dcount--) {
            dcountdim3[dcount]=dcountdim1[dcount]/sets[arSet[eq_time[eqindx]].setid].size;
          }
          /* dcount>-1: eq_time==-1 ran one extra iteration writing
             dcountdim3[-1] (stack scribble) */
          for (dcount=fdim-1; dcount>eq_time[eqindx]-1&&dcount>-1; dcount--) {
            dcountdim3[dcount]=dcountdim1[dcount];
          }
        }
        if(alltimeset<0&&allregset>=0) {
          for (dcount=eq_reg[eqindx]-1; dcount>-1; dcount--) {
            dcountdim4[dcount]=dcountdim1[dcount]/sets[arSet[eq_reg[eqindx]].setid].size;
          }
          /* dcount>-1: eq_reg==-1 ran one extra iteration writing
             dcountdim4[-1] (stack scribble; ASan-confirmed on nsub-dbbd) */
          for (dcount=fdim-1; dcount>eq_reg[eqindx]-1&&dcount>-1; dcount--) {
            dcountdim4[dcount]=dcountdim1[dcount];
          }
        }
        if(alltimeset>=0&&allregset>=0) {
          i3=1;
          dcountdim4[fdim-1]=dcountdim1[fdim-1];
          for (dcount=fdim-2; dcount>-1; dcount--) {
            if(dcount+1==eq_reg[eqindx])i3*=sets[arSet[eq_reg[eqindx]].setid].size;
            if(dcount+1==eq_time[eqindx])i3*=sets[arSet[eq_time[eqindx]].setid].size;
            dcountdim4[dcount]=dcountdim1[dcount]/i3;
          }
        }
      }
      for (i=0; i<nlinvars; i++) {
        i3=0;
        nloopslin=nloops;
        if (fdim==0) {
          for (l2=0; l2<vars[LinVars[i].LinVarIndx].size; l2++) {
            strcpy(arSet[i3].index_name,LinVars[i].dimnames[l2]);
            for (i4=0; i4<nset; i4++) if(strcmp(sets[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                break;
              }
            arSet[i3].setid=i4;
            nloopslin=nloopslin*sets[arSet[i3].setid].size;
            if(sets[i4].intertemp)*sbbd_overrid=true;
            i3++;
          }
        }
        else for (l2=0; l2<vars[LinVars[i].LinVarIndx].size; l2++) {
            for (sj=0; sj<fdim; sj++) {
              if(strcmp(arSet[sj].index_name,LinVars[i].dimnames[l2])==0) {
                break;
              }
              else if (sj==fdim-1) {
                strcpy(arSet[fdim+i3].index_name,LinVars[i].dimnames[l2]);
                for (i4=0; i4<nset; i4++) if(strcmp(sets[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                    break;
                  }
                arSet[fdim+i3].setid=i4;
                nloopslin=nloopslin*sets[arSet[fdim+i3].setid].size;
                if(sets[i4].intertemp)*sbbd_overrid=true;
                i3++;
              }
            }
          }
        fdimlin=fdim+i3;
        if (i3>0) {
          dcountdim2[fdimlin-1]=1;
          for (dcount=fdimlin-2; dcount>-1; dcount--) {
            dcountdim2[dcount]=sets[arSet[dcount+1].setid].size*dcountdim2[dcount+1];
          }
        }
        else for (dcount=0; dcount<fdim; dcount++) {
            dcountdim2[dcount]=dcountdim1[dcount];
          }
        nloopsfac=(offset_t)nloopslin/nloops;
        if(eqmeta!=NULL) { /* probe structure capture: one update per reference, merged by variable */
          for(i4=0; i4<eqmeta[eqindx].nvars_ref; i4++)if(eqmeta[eqindx].var_ref[i4]==LinVars[i].LinVarIndx)break;
          if(i4<eqmeta[eqindx].nvars_ref)eqmeta[eqindx].var_w[i4]+=nloopslin;
          else if(eqmeta[eqindx].nvars_ref<PROBE_MAXEQVARS) {
            eqmeta[eqindx].var_ref[eqmeta[eqindx].nvars_ref]=LinVars[i].LinVarIndx;
            eqmeta[eqindx].var_w[eqmeta[eqindx].nvars_ref]=nloopslin;
            eqmeta[eqindx].nvars_ref++;
          }
        }
        for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
          for (i4=0; i4<fdimlin; i4++) {
            if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].index_name)==0) {
              dcountdim5[dcount]=i4;
              break;
            }
          }
        }
        for(dcount=0; dcount<MAXSUPSET; dcount++)supset[dcount]=0;
        for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
          if(LinVars[i].dimmapid[dcount]>0) {
            linvar_map_dim_check(&LinVars[i],dcount,arSet[dcountdim5[dcount]].setid,vars);
            supset[dcount]=0;
            continue;
          }
          if(sets[vars[LinVars[i].LinVarIndx].setid[dcount]].size!=sets[arSet[dcountdim5[dcount]].setid].size) {
            for(sup=1; sup<MAXSUPSET; sup++)if(vars[LinVars[i].LinVarIndx].setid[dcount]==sets[arSet[dcountdim5[dcount]].setid].subsetid[sup]) {
                supset[dcount]=sup;
                break;
              }
          }
          else supset[dcount]=0;
        }
        for (lj=0; lj<nloopslin; lj++) {
          l2=lj;
          rowindx=0;
          rowindxorg=0;
          for (dcount=0; dcount<fdimlin; dcount++) {
            l1=(offset_t) l2/dcountdim2[dcount];
            arSet[dcount].indx=l1;
            if(eq_intertemp[eqindx]) {
              if(alltimeset>=0&&allregset<0) {
                if (dcount<fdim&&dcount!=eq_time[eqindx]) {
                  rowindx=rowindx+l1*dcountdim3[dcount];
                }
                if(dcount==eq_time[eqindx])ltime=l1;
              }
              if(alltimeset<0&&allregset>=0) {
                if (dcount<fdim&&dcount!=eq_reg[eqindx]) {
                  rowindx=rowindx+l1*dcountdim4[dcount];
                }
                if(dcount==eq_reg[eqindx])ltime=l1;
              }
              if(alltimeset>=0&&allregset>=0) {
                if(nesteddbbd==1) {
                  if(eq_reg[eqindx]>-1) {
                    if (dcount<fdim&&dcount!=eq_time[eqindx]&&dcount!=eq_reg[eqindx])rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  else {
                    if (dcount<fdim&&dcount!=eq_time[eqindx])rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  if(dcount==eq_time[eqindx])ltime=l1;
                  if(dcount==eq_reg[eqindx])lreg=l1;
                }
                else {
                  if (dcount<fdim&&dcount!=eq_time[eqindx]&&dcount!=eq_reg[eqindx]) {
                    rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  if(dcount==eq_time[eqindx])ltime=l1;
                  if(dcount==eq_reg[eqindx])lreg=l1;
                }
              }
            }
            else {
              if (dcount<fdim) {
                rowindx=rowindx+l1*dcountdim1[dcount];
              }
            }
            l2=l2-l1*dcountdim2[dcount];
          }
          li3=0;
          for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
            if(LinVars[i].dimmapid[dcount]>0) {
              /* must mirror the fill loop exactly or dnnz/onnz miscount */
              li3=li3+((offset_t)teems_maps[LinVars[i].dimmapid[dcount]-1].values[arSet[dcountdim5[dcount]].indx])*vars[LinVars[i].LinVarIndx].strides[dcount];
            }
            else if(supset[dcount]==0) {
              li3=li3+(arSet[dcountdim5[dcount]].indx+LinVars[i].dimleadlag[dcount])*vars[LinVars[i].LinVarIndx].strides[dcount];
            }
            else {
              li3=li3+(set_elems[sets[arSet[dcountdim5[dcount]].setid].offset+arSet[dcountdim5[dcount]].indx].superset_pos[supset[dcount]]+LinVars[i].dimleadlag[dcount])*vars[LinVars[i].LinVarIndx].strides[dcount];
            }
          }

          if(eq_intertemp[eqindx]) {
            if(alltimeset>=0&&allregset>=0) {
              if(nesteddbbd==1) {
                if(eq_reg[eqindx]>-1)Jindx=counteq1[(set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[sets[arSet[eq_time[eqindx]].setid].intsup])*(nreg+1)+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+lreg].superset_pos[sets[arSet[eq_reg[eqindx]].setid].regsup]]+rowindx;
                else Jindx=counteq1[(set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[sets[arSet[eq_time[eqindx]].setid].intsup])*(nreg+1)+nreg]+rowindx;
              }
              else {
                if(arSet[eq_time[eqindx]].setid==alltimeset)
                  Jindx=counteq1[(set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[0])*sets[arSet[eq_reg[eqindx]].setid].size+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+lreg].superset_pos[0]]+rowindx;
                else {
                  for(i4=1; i4<MAXSUPSET; i4++)if(sets[arSet[eq_time[eqindx]].setid].subsetid[i4]==alltimeset)break;
                  Jindx=counteq1[(set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[i4])*sets[arSet[eq_reg[eqindx]].setid].size+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+lreg].superset_pos[0]]+rowindx;
                }
              }
            }
            if(alltimeset>=0&&allregset<0) {
              if(arSet[eq_time[eqindx]].setid==alltimeset)
                Jindx=counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[0]]+rowindx;
              else {
                for(i4=1; i4<MAXSUPSET; i4++)if(sets[arSet[eq_time[eqindx]].setid].subsetid[i4]==alltimeset)break;
                Jindx=counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[i4]]+rowindx;
              }
            }
            if(alltimeset<0&&allregset>=0) Jindx=counteq1[set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+ltime].superset_pos[0]]+rowindx;
          }
          else Jindx=matrow+rowindx;
          leq=matroworg+(offset_t)lj/nloopsfac;//matroworg+rowindxorg;
          eq_addr[leq]=Jindx;
          Iindx=closure_vals[vars[LinVars[i].LinVarIndx].offset+li3].exo_index;
          if(Istart<=Jindx&&Jindx<Iend) {
            if (!closure_vals[vars[LinVars[i].LinVarIndx].offset+li3].is_exogenous) {
              if (Istart<=Iindx&&Iindx<Iend) {
                dnnz[Jindx-Istart]=dnnz[Jindx-Istart]+1;
              }
              else {
                onnz[Jindx-Istart]=onnz[Jindx-Istart]+1;
              }
            }
            if (closure_vals[vars[LinVars[i].LinVarIndx].offset+li3].is_exogenous) {
              if (Istart<=Iindx&&Iindx<Iend) { //&&Iindx<nexo
                dnnzB[Jindx-Istart]=dnnzB[Jindx-Istart]+1;
              }
              else {
                onnzB[Jindx-Istart]=onnzB[Jindx-Istart]+1;
              }
            }
          }
        }
      }
      if(eq_intertemp[eqindx]) {
        if(alltimeset>=0&&allregset>=0) {
          if(nesteddbbd==1) {
            if(eq_reg[eqindx]>-1)for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)for(l2=0; l2<sets[arSet[eq_reg[eqindx]].setid].size; l2++) {
                  counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[sets[arSet[eq_time[eqindx]].setid].intsup]*(nreg+1)+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+l2].superset_pos[sets[arSet[eq_reg[eqindx]].setid].regsup]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size/sets[arSet[eq_reg[eqindx]].setid].size;
                  if(sets[arSet[eq_reg[eqindx]].setid].regsup!=0)printf("Error: subsets of the regional partition set are not supported in NDBBD ordering\n");
                }
            else for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)
                counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[sets[arSet[eq_time[eqindx]].setid].intsup]*(nreg+1)+nreg]+=nloops/sets[arSet[eq_time[eqindx]].setid].size;
          }
          else {
            if(arSet[eq_time[eqindx]].setid==alltimeset)
              for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)for(l2=0; l2<sets[arSet[eq_reg[eqindx]].setid].size; l2++)
                  counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[0]*sets[arSet[eq_reg[eqindx]].setid].size+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+l2].superset_pos[0]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size/sets[arSet[eq_reg[eqindx]].setid].size;
            else {
              for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)for(l2=0; l2<sets[arSet[eq_reg[eqindx]].setid].size; l2++) {
                  for(i4=1; i4<MAXSUPSET; i4++)if(sets[arSet[eq_time[eqindx]].setid].subsetid[i4]==alltimeset)break;
                  counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[i4]*sets[arSet[eq_reg[eqindx]].setid].size+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+l2].superset_pos[0]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size/sets[arSet[eq_reg[eqindx]].setid].size;
                }
            }
          }
        }
        if(alltimeset>=0&&allregset<0) {
          if(arSet[eq_time[eqindx]].setid==alltimeset)
            for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[0]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size;
          else {
            for(i4=1; i4<MAXSUPSET; i4++)if(sets[arSet[eq_time[eqindx]].setid].subsetid[i4]==alltimeset)break;
            for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[i4]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size;
          }
        }
        if(alltimeset<0&&allregset>=0)for(lj=0; lj<sets[arSet[eq_reg[eqindx]].setid].size; lj++)counteq1[set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+lj].superset_pos[0]]+=nloops/sets[arSet[eq_reg[eqindx]].setid].size;
      }
      else matrow+=nloops;
      if(eqmeta!=NULL) {
        eqmeta[eqindx].base=matroworg;
        eqmeta[eqindx].nrows=nloops;
        eqmeta[eqindx].fdim=(fdim<4*MAXVARDIM)?fdim:4*MAXVARDIM;
        for(i4=0; i4<eqmeta[eqindx].fdim; i4++)eqmeta[eqindx].setid[i4]=arSet[i4].setid;
      }
      matroworg+=nloops;
      eqindx++;
      free(LinVars);
      free(arSet);
    }
  }
  if(neqmeta!=NULL)*neqmeta=eqindx;
  free(counteq1);
  fclose(filehandle);
  return 1;
}

