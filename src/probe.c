#include <teems_solver.h>
#include <sys/time.h>
#include "hsl_mc79i.h"

/* -solmed probe structural diagnosis (HSL_MC79, ROADMAP 6.7).
   The probe run has executed the full pre-solve pipeline (data,
   formulas, closure, ordering, preallocation) but no solve method has
   assembled a Jacobian.  This leg assembles the condensed Jacobian
   sequentially, then runs the MC79 maximum matching on
     (a) the full stored pattern   -> structural closure validity, and
     (b) the numerically realized pattern (entries nonzero at base
         data) -> the zero-flow singularity class that today only
         surfaces as a mid-solve rank message with no names.
   Solvability of the condensed system = a perfect matching on the
   element-level incidence graph, so unmatched columns are exactly the
   under-determined variable elements and unmatched rows the
   over-constrained equations; both are translated back to names via
   the closure exo_index map (columns) and the eq_addr map plus the
   statement metadata captured by jacobian_preallocate (rows).  On
   defects the coarse Dulmage-Mendelsohn decomposition localizes the
   entangled region, with the defect and entangled sets additionally
   aggregated by variable / equation statement so the report stays
   readable at any scale; -probefine adds the fine decomposition's
   strongly connected components (the model's irreducible simultaneous
   cores) with the composition of the largest cores.  The report
   (version 2) also carries the statement-level equation-system
   structure: per statement, its quantifier sets, row count and the
   variables it references with element-level incidence weights.
   Results go to the log and to <solfiles>.probe.json. */

#define PROBE_JSON_LIST_MAX 200
#define PROBE_LOG_LIST_MAX 20
#define PROBE_AGG_LOG_MAX 10
#define PROBE_TOP_CORES 20
#define PROBE_CORE_NAMES_MAX 50

/* 1-based CSC of A's stored pattern (realized!=0 drops entries whose
   stored value is zero).  CSR row-major order makes each column's rows
   come out sorted, as MC79 requires.  Returns kept entry count. */
static PetscInt probe_csc_build(Mat A,PetscInt n,int realized,int **ptr_out,int **row_out) {
  Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
  PetscInt i,k,nkeep=0;
  int *ptr= (int *) calloc (n+1,sizeof(int));
  for(k=0; k<aa->nz; k++) {
    if(realized&&aa->a[k]==0.0)continue;
    ptr[aa->j[k]+1]++;
    nkeep++;
  }
  ptr[0]=1;
  for(i=1; i<=n; i++)ptr[i]+=ptr[i-1];
  int *row= (int *) calloc ((nkeep>0)?nkeep:1,sizeof(int));
  for(i=0; i<n; i++)for(k=aa->i[i]; k<aa->i[i+1]; k++) {
      if(realized&&aa->a[k]==0.0)continue;
      row[ptr[aa->j[k]]-1]=i+1;
      ptr[aa->j[k]]++;
    }
  for(i=n; i>0; i--)ptr[i]=ptr[i-1];
  ptr[0]=1;
  *ptr_out=ptr;
  *row_out=row;
  return nkeep;
}

/* variable index owning closure position j5 (vars[] offsets ascend) */
static offset_t probe_var_of_ele(offset_t j5,array_def *vars,offset_t nvar) {
  offset_t lo=0,hi=nvar-1,mid;
  while(lo<=hi) {
    mid=(lo+hi)/2;
    if(j5<vars[mid].offset)hi=mid-1;
    else if(j5>=vars[mid].offset+vars[mid].nelem)lo=mid+1;
    else return mid;
  }
  return -1;
}

/* statement index owning eq_addr position leq (bases ascend) */
static offset_t probe_stmt_of_leq(offset_t leq,eq_probe_meta *eqmeta,offset_t neqmeta) {
  offset_t lo=0,hi=neqmeta-1,mid;
  if(eqmeta==NULL)return -1;
  while(lo<=hi) {
    mid=(lo+hi)/2;
    if(leq<eqmeta[mid].base)hi=mid-1;
    else if(leq>=eqmeta[mid].base+eqmeta[mid].nrows)lo=mid+1;
    else return mid;
  }
  return -1;
}

/* "name(ele,ele)" for the variable element at closure position j5 */
static void probe_var_ele_name(offset_t j5,array_def *vars,offset_t nvar,set_def *sets,set_element *set_elems,char *out) {
  offset_t i,j,idx;
  dim_t d;
  i=probe_var_of_ele(j5,vars,nvar);
  if(i<0) {
    sprintf(out,"variable element %ld",j5);
    return;
  }
  strcpy(out,vars[i].cofname);
  if(vars[i].size>0) {
    strcat(out,"(");
    j=j5-vars[i].offset;
    for(d=0; d<vars[i].size; d++) {
      idx=j/vars[i].strides[d];
      j-=idx*vars[i].strides[d];
      if(d>0)strcat(out,",");
      strcat(out,set_elems[sets[vars[i].setid[d]].offset+idx].setele);
    }
    strcat(out,")");
  }
}

/* "eqname(ele,ele)" for the equation element at eq_addr position leq */
static void probe_eq_ele_name(offset_t leq,eq_probe_meta *eqmeta,offset_t neqmeta,set_def *sets,set_element *set_elems,char *out) {
  offset_t k,off,stride,idx;
  dim_t d,d1;
  k=probe_stmt_of_leq(leq,eqmeta,neqmeta);
  if(k<0) {
    sprintf(out,"equation row element %ld",leq);
    return;
  }
  strcpy(out,eqmeta[k].eqname);
  if(eqmeta[k].fdim>0) {
    strcat(out,"(");
    off=leq-eqmeta[k].base;
    for(d=0; d<eqmeta[k].fdim; d++) {
      stride=1;
      for(d1=d+1; d1<eqmeta[k].fdim; d1++)stride*=sets[eqmeta[k].setid[d1]].size;
      idx=off/stride;
      off-=idx*stride;
      if(d>0)strcat(out,",");
      strcat(out,set_elems[sets[eqmeta[k].setid[d]].offset+idx].setele);
    }
    strcat(out,")");
  }
}

/* ---- aggregation by variable / statement ------------------------------- */

typedef struct {
  offset_t idx;
  long c;
} probe_cnt ;

static int probe_cnt_cmp(const void *a,const void *b) {
  long d=((const probe_cnt*)b)->c-((const probe_cnt*)a)->c;
  return (d>0)-(d<0);
}

/* count matrix columns per variable (which=0) or rows per statement
   (which=1) over the 0-based indices in items[]; returns the nonzero
   (idx,count) pairs sorted descending (caller frees), nnz in *nnz_out */
static probe_cnt *probe_agg(PetscInt *items,PetscInt nitems,int which,
                            offset_t *col2ele,offset_t *row2leq,
                            array_def *vars,offset_t nvar,
                            eq_probe_meta *eqmeta,offset_t neqmeta,
                            offset_t *nnz_out) {
  offset_t n=(which==0)?nvar:neqmeta,k,nnz=0;
  PetscInt i;
  long *counts= (long *) calloc ((n>0)?n:1,sizeof(long));
  for(i=0; i<nitems; i++) {
    if(which==0)k=probe_var_of_ele(col2ele[items[i]],vars,nvar);
    else k=probe_stmt_of_leq(row2leq[items[i]],eqmeta,neqmeta);
    if(k>=0)counts[k]++;
  }
  for(k=0; k<n; k++)if(counts[k]>0)nnz++;
  probe_cnt *sorted= (probe_cnt *) malloc (((nnz>0)?nnz:1)*sizeof(probe_cnt));
  nnz=0;
  for(k=0; k<n; k++)if(counts[k]>0) {
      sorted[nnz].idx=k;
      sorted[nnz].c=counts[k];
      nnz++;
    }
  qsort(sorted,nnz,sizeof(probe_cnt),probe_cnt_cmp);
  free(counts);
  *nnz_out=nnz;
  return sorted;
}

static const char *probe_agg_name(probe_cnt *e,int which,array_def *vars,eq_probe_meta *eqmeta) {
  return (which==0)?vars[e->idx].cofname:eqmeta[e->idx].eqname;
}

/* JSON: [{"n": "name", "c": N}, ...] (no key, no trailing comma) */
static void probe_agg_json(FILE *fp,probe_cnt *sorted,offset_t nnz,int which,array_def *vars,eq_probe_meta *eqmeta,offset_t cap) {
  offset_t k;
  fprintf(fp,"[");
  for(k=0; k<nnz&&k<cap; k++)fprintf(fp,"%s{\"n\": \"%s\", \"c\": %ld}",(k>0)?", ":"",probe_agg_name(&sorted[k],which,vars,eqmeta),sorted[k].c);
  fprintf(fp,"]");
}

static void probe_agg_log(const char *label,probe_cnt *sorted,offset_t nnz,int which,array_def *vars,eq_probe_meta *eqmeta) {
  offset_t k;
  printf("probe:   %s:",label);
  for(k=0; k<nnz&&k<PROBE_AGG_LOG_MAX; k++)printf(" %s x%ld",probe_agg_name(&sorted[k],which,vars,eqmeta),sorted[k].c);
  if(nnz>PROBE_AGG_LOG_MAX)printf(" ... (%ld more)",(long)(nnz-PROBE_AGG_LOG_MAX));
  printf("\n");
}

/* aggregate + emit one keyed JSON line and one log line, then free */
static void probe_agg_report(FILE *fp,const char *jsonkey,const char *label,
                             PetscInt *items,PetscInt nitems,int which,
                             offset_t *col2ele,offset_t *row2leq,
                             array_def *vars,offset_t nvar,
                             eq_probe_meta *eqmeta,offset_t neqmeta) {
  offset_t nnz=0;
  probe_cnt *sorted=probe_agg(items,nitems,which,col2ele,row2leq,vars,nvar,eqmeta,neqmeta,&nnz);
  probe_agg_log(label,sorted,nnz,which,vars,eqmeta);
  if(fp!=NULL) {
    fprintf(fp,"    \"%s\": ",jsonkey);
    probe_agg_json(fp,sorted,nnz,which,vars,eqmeta,PROBE_JSON_LIST_MAX);
    fprintf(fp,",\n");
  }
  free(sorted);
}

/* log (capped) + JSON array (capped) of the defect columns or rows.
   which: 0 = columns (variables), 1 = rows (equations) */
static void probe_defect_report(FILE *fp,const char *jsonkey,const char *logline,int *match,PetscInt n,int which,
                                offset_t *col2ele,offset_t *row2leq,
                                array_def *vars,offset_t nvar,eq_probe_meta *eqmeta,offset_t neqmeta,
                                set_def *sets,set_element *set_elems) {
  PetscInt i,cnt=0;
  char name[TABREADLINE];
  printf("probe:   %s:",logline);
  for(i=0; i<n&&cnt<PROBE_LOG_LIST_MAX; i++)if(match[i]==0) {
      if(which==0)probe_var_ele_name(col2ele[i],vars,nvar,sets,set_elems,name);
      else probe_eq_ele_name(row2leq[i],eqmeta,neqmeta,sets,set_elems,name);
      printf(" %s",name);
      cnt++;
    }
  for(; i<n; i++)if(match[i]==0)cnt++;
  if(cnt>PROBE_LOG_LIST_MAX)printf(" ... (%d more)",(int)(cnt-PROBE_LOG_LIST_MAX));
  printf("\n");
  if(fp==NULL)return;
  fprintf(fp,"    \"%s\": [",jsonkey);
  cnt=0;
  for(i=0; i<n&&cnt<PROBE_JSON_LIST_MAX; i++)if(match[i]==0) {
      if(which==0)probe_var_ele_name(col2ele[i],vars,nvar,sets,set_elems,name);
      else probe_eq_ele_name(row2leq[i],eqmeta,neqmeta,sets,set_elems,name);
      fprintf(fp,"%s\"%s\"",(cnt>0)?", ":"",name);
      cnt++;
    }
  fprintf(fp,"],\n");
}

/* ---- on-failure diagnosis (adoption-shape part 3) ---------------------
   The factorization kernels hard-STOP on MA48B/BD failure; just before
   the STOP they call teems_onfail_diag_() (the Fortran-visible entry
   below).  The C solve paths register the system they are about to
   factorize -- the live PETSc matrix plus how its local rows/columns
   map back to condensed-system positions -- so the diagnosis can run
   the MC79 matching at the failure point, on the failing (sub)system
   and its current values, and name the defective elements.  When both
   patterns still have full structural rank the failure is numerical
   (near-singularity), which is reported as such -- the recorded MC79
   limitation.  Everything here is failure-path-only: registration is
   a few pointer stores per factorization. */

static struct {
  set_def *sets;
  set_element *set_elems;
  array_def *vars;
  offset_t nvar;
  eq_probe_meta *eqmeta;
  offset_t neqmeta;
  PetscInt VecSize;
  /* condensed position -> element, built eagerly at registration:
     closure_vals is spilled and freed mid-solve in disk mode, so it
     cannot be consulted at failure time */
  int *col2ele;
  int *row2leq;
  int set;
} onfail_ctx;

typedef struct {
  Mat A;                        /* live matrix (NULL = label-only note) */
  PetscInt m,n;
  const char *label;
  int block_id;                 /* >=0: printed with the label */
  int *row_order,*col_order;    /* NULL = identity map */
  offset_t row_base,col_base;   /* index offset into the order arrays */
  offset_t row_add,col_add;     /* added to the order-array value
                                   (DBBD: the base again; NDBBD: 0) */
  int set;
} probe_onfail_scope_t;
static __thread probe_onfail_scope_t onfail_scope;

void probe_onfail_context(set_def *sets,set_element *set_elems,array_def *vars,offset_t nvar,closure_entry *closure_vals,offset_t nvarele,offset_t *eq_addr,eq_probe_meta *eqmeta,offset_t neqmeta,PetscInt VecSize) {
  PetscInt i;
  offset_t j5;
  onfail_ctx.sets=sets;
  onfail_ctx.set_elems=set_elems;
  onfail_ctx.vars=vars;
  onfail_ctx.nvar=nvar;
  onfail_ctx.eqmeta=eqmeta;
  onfail_ctx.neqmeta=neqmeta;
  onfail_ctx.VecSize=VecSize;
  onfail_ctx.col2ele= (int *) malloc (VecSize*sizeof(int));
  onfail_ctx.row2leq= (int *) malloc (VecSize*sizeof(int));
  for(i=0; i<VecSize; i++) {
    onfail_ctx.col2ele[i]=-1;
    onfail_ctx.row2leq[i]=-1;
  }
  for(j5=0; j5<nvarele; j5++)if(!closure_vals[j5].is_exogenous&&!closure_vals[j5].is_backsolved) {
      if(closure_vals[j5].exo_index>=0&&closure_vals[j5].exo_index<VecSize)onfail_ctx.col2ele[closure_vals[j5].exo_index]=(int)j5;
    }
  for(j5=0; j5<VecSize; j5++)if(eq_addr[j5]>=0&&eq_addr[j5]<VecSize)onfail_ctx.row2leq[eq_addr[j5]]=(int)j5;
  onfail_ctx.set=1;
}

void probe_onfail_scope_set(Mat A,PetscInt m,PetscInt n,const char *label,int block_id,int *row_order,int *col_order,offset_t row_base,offset_t col_base,offset_t row_add,offset_t col_add) {
  onfail_scope.A=A;
  onfail_scope.m=m;
  onfail_scope.n=n;
  onfail_scope.label=label;
  onfail_scope.block_id=block_id;
  onfail_scope.row_order=row_order;
  onfail_scope.col_order=col_order;
  onfail_scope.row_base=row_base;
  onfail_scope.col_base=col_base;
  onfail_scope.row_add=row_add;
  onfail_scope.col_add=col_add;
  onfail_scope.set=1;
}

void probe_onfail_scope_clear(void) {
  onfail_scope.set=0;
  onfail_scope.A=NULL;
}

/* Fortran entry: collective abort after a diagnosed failure — a plain
   STOP on one rank leaves the others waiting in the next collective */
void teems_onfail_abort_(void) {
  MPI_Abort(PETSC_COMM_WORLD,1);
}

/* Fortran entry: CALL TEEMS_ONFAIL_DIAG(INFO(1)) before the failure
   aborts.  info1 is MA48's error code: only singularity-class failures
   get the (potentially expensive) matching — a workspace failure
   (INFO(1) = -3, LA too small) is not a structural question and at
   large scale the matching would add minutes to a failure whose remedy
   is just a bigger -laA. */
void teems_onfail_diag_(int *info1) {
  #pragma omp critical(teems_onfail)
  {
    probe_onfail_scope_t *s=&onfail_scope;
    int code=(info1!=NULL)?*info1:0;
    if(code==-3) {
      printf("probe: MA48 workspace failure (INFO(1) = -3): not a singularity — increase -laA (no structural diagnosis run)\n");
      goto onfail_done;
    }
    PetscInt i,nz,pass;
    offset_t g;
    int *ptr=NULL,*row=NULL;
    struct mc79_control_i control;
    struct mc79_info_i info;
    if(!s->set||!onfail_ctx.set) {
      printf("probe: factorization failed in an unregistered context — no structural diagnosis available\n");
    }
    else {
      if(s->block_id>=0)printf("probe: diagnosing the failed factorization: %s %d (%ld x %ld)\n",s->label,s->block_id,(long)s->m,(long)s->n);
      else printf("probe: diagnosing the failed factorization: %s (%ld x %ld)\n",s->label,(long)s->m,(long)s->n);
      if(s->A==NULL) {
        printf("probe: structural diagnosis is not wired for this system yet; run -solmed probe for the deploy-level diagnosis\n");
      }
      else if(s->m!=s->n) {
        printf("probe: non-square scope (%ld x %ld), diagnosis skipped\n",(long)s->m,(long)s->n);
      }
      else {
        /* local matrix position -> condensed position -> element */
        offset_t *col2ele_l= (offset_t *) malloc (s->n*sizeof(offset_t));
        offset_t *row2leq_l= (offset_t *) malloc (s->m*sizeof(offset_t));
        for(i=0; i<s->n; i++) {
          g=(s->col_order!=NULL)?(offset_t)s->col_order[s->col_base+i]+s->col_add:(offset_t)i;
          col2ele_l[i]=(g>=0&&g<onfail_ctx.VecSize)?(offset_t)onfail_ctx.col2ele[g]:-1;
        }
        for(i=0; i<s->m; i++) {
          g=(s->row_order!=NULL)?(offset_t)s->row_order[s->row_base+i]+s->row_add:(offset_t)i;
          row2leq_l[i]=(g>=0&&g<onfail_ctx.VecSize)?(offset_t)onfail_ctx.row2leq[g]:-1;
        }
        int *rowmatch= (int *) malloc (s->m*sizeof(int));
        int *colmatch= (int *) malloc (s->n*sizeof(int));
        PetscInt *items= (PetscInt *) malloc (s->n*sizeof(PetscInt));
        PetscInt nitems;
        mc79_default_control_i(&control);
        control.f_arrays=1;
        int defect=0;
        for(pass=0; pass<2; pass++) {
          const char *pname=(pass==0)?"stored pattern":"currently nonzero pattern";
          nz=probe_csc_build(s->A,s->n,pass,&ptr,&row);
          mc79_matching_i(s->m,s->n,ptr,row,rowmatch,colmatch,&control,&info);
          if(info.flag<0) {
            printf("probe: MC79 matching failed (%s, flag %d)\n",pname,info.flag);
          }
          else if(info.mbar==0&&info.nbar==0) {
            printf("probe: %s: full structural rank %ld of %ld\n",pname,(long)s->n,(long)s->n);
          }
          else {
            defect=1;
            printf("probe: %s: structurally singular — rank %ld of %ld (%d unmatched equations, %d unmatched variables)\n",pname,(long)(s->n-info.nbar),(long)s->n,info.mbar,info.nbar);
            probe_defect_report(NULL,"","under-determined variable elements",colmatch,s->n,0,col2ele_l,row2leq_l,onfail_ctx.vars,onfail_ctx.nvar,onfail_ctx.eqmeta,onfail_ctx.neqmeta,onfail_ctx.sets,onfail_ctx.set_elems);
            probe_defect_report(NULL,"","over-constrained equation elements",rowmatch,s->m,1,col2ele_l,row2leq_l,onfail_ctx.vars,onfail_ctx.nvar,onfail_ctx.eqmeta,onfail_ctx.neqmeta,onfail_ctx.sets,onfail_ctx.set_elems);
            nitems=0;
            for(i=0; i<s->n; i++)if(colmatch[i]==0)items[nitems++]=i;
            probe_agg_report(NULL,"","under-determined by variable",items,nitems,0,col2ele_l,row2leq_l,onfail_ctx.vars,onfail_ctx.nvar,onfail_ctx.eqmeta,onfail_ctx.neqmeta);
            nitems=0;
            for(i=0; i<s->m; i++)if(rowmatch[i]==0)items[nitems++]=i;
            probe_agg_report(NULL,"","over-constrained by equation",items,nitems,1,col2ele_l,row2leq_l,onfail_ctx.vars,onfail_ctx.nvar,onfail_ctx.eqmeta,onfail_ctx.neqmeta);
            mc79_coarse_i(s->m,s->n,ptr,row,rowmatch,colmatch,&control,&info);
            if(info.flag>=0)printf("probe:   DM localization: under %d x %d, well %d x %d, over %d x %d\n",info.m1,info.n1,info.m2,info.n2,info.m3,info.n3);
            free(ptr);
            free(row);
            ptr=NULL;
            row=NULL;
            break; /* structural verdict found; no need for the second pass */
          }
          free(ptr);
          free(row);
          ptr=NULL;
          row=NULL;
        }
        if(!defect)printf("probe: both patterns have full structural rank at the failure state — the failure is numerical (near-singular values / ill-conditioning), not structural; the MA48 rank output above is the primary evidence\n");
        free(rowmatch);
        free(colmatch);
        free(items);
        free(col2ele_l);
        free(row2leq_l);
      }
    }
    onfail_done: ;
  }
}

int probe_structural(PetscInt VecSize,offset_t nvarele,offset_t ncofele,PetscInt dnz,PetscInt *dnnz,PetscInt dnzB,PetscInt *dnnzB,char *tabfile,char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,offset_t *eq_addr,offset_t *counteq,offset_t nintraeq,eq_probe_meta *eqmeta,offset_t neqmeta,cmf_file_entry *iodata,int niodata,int noutdata,int nsoldata,int probefine,PetscInt mpisize,PetscInt rank) {
  Mat A,B;
  struct timeval t0,t1;
  offset_t j5,k,d;
  PetscInt i,pass,nz;
  int *ptr=NULL,*row=NULL,*rowmatch=NULL,*colmatch=NULL;
  struct mc79_control_i control;
  struct mc79_info_i info;
  char name[TABREADLINE];
  if(rank!=0)return 0;
  if(mpisize>1) {
    printf("probe: the MC79 structural diagnosis runs on a single rank; rerun with -n 1 to enable it\n");
    return 0;
  }
  if(VecSize<1) {
    printf("probe: empty system, structural diagnosis skipped\n");
    return 0;
  }
  gettimeofday(&t0,NULL);
  PetscInt BSize=(PetscInt)(nvarele-VecSize-nbselems);
  BSize=(BSize>VecSize)?BSize:VecSize;
  /* always a sequential copy: the probe's diagnosis is serial and the
     preallocation counts cover the whole system on this rank */
  MatCreate(PETSC_COMM_SELF,&A);
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
  MatSetType(A,MATSEQAIJ);
  MatSeqAIJSetPreallocation(A,dnz,dnnz);
  MatSetOption(A,MAT_SYMMETRIC,PETSC_FALSE);
  MatCreate(PETSC_COMM_SELF,&B);
  MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,BSize);
  MatSetType(B,MATSEQAIJ);
  MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
  jacobian_fill(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,ndblock,alltimeset,allregset,eq_addr,counteq,nintraeq,A,B);
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
  gettimeofday(&t1,NULL);
  logmsg(1,"probe: Jacobian assembled for structural diagnosis (%d x %d, %ld stored entries) in %.2f s\n",VecSize,VecSize,(long)((Mat_SeqAIJ*)A->data)->nz,(t1.tv_sec-t0.tv_sec)+((double)(t1.tv_usec-t0.tv_usec))/1000000);

  /* column -> variable element, row -> eq_addr position */
  offset_t *col2ele= (offset_t *) malloc (VecSize*sizeof(offset_t));
  offset_t *row2leq= (offset_t *) malloc (VecSize*sizeof(offset_t));
  for(i=0; i<VecSize; i++) {
    col2ele[i]=-1;
    row2leq[i]=-1;
  }
  for(j5=0; j5<nvarele; j5++)if(!closure_vals[j5].is_exogenous&&!closure_vals[j5].is_backsolved) {
      if(closure_vals[j5].exo_index>=0&&closure_vals[j5].exo_index<VecSize)col2ele[closure_vals[j5].exo_index]=j5;
    }
  for(j5=0; j5<VecSize; j5++)if(eq_addr[j5]>=0&&eq_addr[j5]<VecSize)row2leq[eq_addr[j5]]=j5;

  mc79_default_control_i(&control);
  control.f_arrays=1;
  rowmatch= (int *) malloc (VecSize*sizeof(int));
  colmatch= (int *) malloc (VecSize*sizeof(int));
  PetscInt *items= (PetscInt *) malloc (VecSize*sizeof(PetscInt));
  PetscInt nitems;

  /* <solfiles>.probe.json next to <solfiles>.stats.json */
  char probepath[TABREADLINE+16];
  int fidx;
  for(fidx=niodata+noutdata; fidx<niodata+noutdata+nsoldata; fidx++)if(strcmp("solfiles",iodata[fidx].logname)==0)break;
  if(fidx<niodata+noutdata+nsoldata)strcpy(probepath,iodata[fidx].filname);
  else strcpy(probepath,"solution");
  strcat(probepath,".probe.json");
  FILE *fp=fopen(probepath,"w");
  if(fp==NULL)printf("Warning: cannot write probe report %s\n",probepath);
  else {
    fprintf(fp,"{\n");
    fprintf(fp,"  \"version\": 2,\n");
    fprintf(fp,"  \"vecsize\": %ld,\n",(long)VecSize);
  }

  /* statement-level equation-system structure: per statement its
     quantifier sets, row count and referenced variables with
     element-level incidence weights (captured by preallocation) */
  if(fp!=NULL&&eqmeta!=NULL&&neqmeta>0) {
    fprintf(fp,"  \"statements\": [\n");
    for(k=0; k<neqmeta; k++) {
      fprintf(fp,"    {\"eq\": \"%s\", \"rows\": %ld, \"dims\": [",eqmeta[k].eqname,(long)eqmeta[k].nrows);
      for(d=0; d<eqmeta[k].fdim; d++)fprintf(fp,"%s\"%s\"",(d>0)?", ":"",sets[eqmeta[k].setid[d]].setname);
      fprintf(fp,"], \"vars\": [");
      for(d=0; d<eqmeta[k].nvars_ref; d++)fprintf(fp,"%s{\"v\": \"%s\", \"w\": %ld}",(d>0)?", ":"",vars[eqmeta[k].var_ref[d]].cofname,(long)eqmeta[k].var_w[d]);
      fprintf(fp,"]}%s\n",(k<neqmeta-1)?",":"");
    }
    fprintf(fp,"  ],\n");
  }

  int defects_total=0;
  for(pass=0; pass<2; pass++) {
    const char *pname=(pass==0)?"structural":"realized";
    gettimeofday(&t0,NULL);
    nz=probe_csc_build(A,VecSize,pass,&ptr,&row);
    mc79_matching_i(VecSize,VecSize,ptr,row,rowmatch,colmatch,&control,&info);
    gettimeofday(&t1,NULL);
    if(info.flag<0) {
      printf("probe: MC79 matching failed on the %s pattern (flag %d) — diagnosis skipped\n",pname,info.flag);
      if(fp!=NULL)fprintf(fp,"  \"%s\": {\"flag\": %d},\n",pname,info.flag);
      free(ptr);
      free(row);
      ptr=NULL;
      row=NULL;
      continue;
    }
    if(pass==0)logmsg(1,"probe: %s pattern: %ld entries, matching %.2f s\n",pname,(long)nz,(t1.tv_sec-t0.tv_sec)+((double)(t1.tv_usec-t0.tv_usec))/1000000);
    else logmsg(1,"probe: %s pattern (nonzero at base data): %ld entries, matching %.2f s\n",pname,(long)nz,(t1.tv_sec-t0.tv_sec)+((double)(t1.tv_usec-t0.tv_usec))/1000000);
    if(fp!=NULL) {
      fprintf(fp,"  \"%s\": {\n",pname);
      fprintf(fp,"    \"flag\": %d,\n",info.flag);
      fprintf(fp,"    \"entries\": %ld,\n",(long)nz);
      fprintf(fp,"    \"rank\": %ld,\n",(long)(VecSize-info.nbar));
      fprintf(fp,"    \"unmatched_rows\": %d,\n",info.mbar);
      fprintf(fp,"    \"unmatched_cols\": %d,\n",info.nbar);
    }
    if(info.mbar==0&&info.nbar==0) {
      if(pass==0)printf("probe: structural pattern: full structural rank %ld of %ld — the closure is structurally valid\n",(long)VecSize,(long)VecSize);
      else printf("probe: realized pattern: full structural rank %ld of %ld — no zero-flow undeterminable elements at base data\n",(long)VecSize,(long)VecSize);
      if(fp!=NULL)fprintf(fp,"    \"defective\": false\n  },\n");
    }
    else {
      defects_total+=info.mbar+info.nbar;
      printf("probe: %s pattern: system is structurally singular — rank %ld of %ld (%d unmatched equations, %d unmatched variables)\n",pname,(long)(VecSize-info.nbar),(long)VecSize,info.mbar,info.nbar);
      if(pass==1)printf("probe:   (defects on the realized pattern only are the zero-flow class: structurally present but zero-valued at base data)\n");
      probe_defect_report(fp,"under_determined_vars","under-determined variable elements",colmatch,VecSize,0,col2ele,row2leq,vars,nvar,eqmeta,neqmeta,sets,set_elems);
      probe_defect_report(fp,"over_constrained_eqs","over-constrained equation elements",rowmatch,VecSize,1,col2ele,row2leq,vars,nvar,eqmeta,neqmeta,sets,set_elems);
      /* aggregation by variable / statement: never truncates the way
         the element lists must, so it stays readable at any scale */
      nitems=0;
      for(i=0; i<VecSize; i++)if(colmatch[i]==0)items[nitems++]=i;
      probe_agg_report(fp,"under_determined_by_var","under-determined by variable",items,nitems,0,col2ele,row2leq,vars,nvar,eqmeta,neqmeta);
      nitems=0;
      for(i=0; i<VecSize; i++)if(rowmatch[i]==0)items[nitems++]=i;
      probe_agg_report(fp,"over_constrained_by_eq","over-constrained by equation",items,nitems,1,col2ele,row2leq,vars,nvar,eqmeta,neqmeta);
      /* coarse DM: localize the entangled region */
      int *rowperm= (int *) malloc (VecSize*sizeof(int));
      int *colperm= (int *) malloc (VecSize*sizeof(int));
      mc79_coarse_i(VecSize,VecSize,ptr,row,rowperm,colperm,&control,&info);
      if(info.flag>=0) {
        printf("probe:   DM localization: under-determined block %d x %d, well-determined block %d x %d, over-determined block %d x %d\n",info.m1,info.n1,info.m2,info.n2,info.m3,info.n3);
        if(info.n1>0&&info.n1<=PROBE_LOG_LIST_MAX) {
          printf("probe:   entangled under-determined columns:");
          for(i=0; i<info.n1; i++) {
            probe_var_ele_name(col2ele[colperm[i]-1],vars,nvar,sets,set_elems,name);
            printf(" %s",name);
          }
          printf("\n");
        }
        if(info.m3>0&&info.m3<=PROBE_LOG_LIST_MAX) {
          printf("probe:   entangled over-determined rows:");
          for(i=VecSize-info.m3; i<VecSize; i++) {
            probe_eq_ele_name(row2leq[rowperm[i]-1],eqmeta,neqmeta,sets,set_elems,name);
            printf(" %s",name);
          }
          printf("\n");
        }
        if(fp!=NULL)fprintf(fp,"    \"dm\": {\"m1\": %d, \"n1\": %d, \"m2\": %d, \"n2\": %d, \"m3\": %d, \"n3\": %d},\n",info.m1,info.n1,info.m2,info.n2,info.m3,info.n3);
        if(info.n1>0) {
          for(i=0; i<info.n1; i++)items[i]=colperm[i]-1;
          probe_agg_report(fp,"dm_under_by_var","DM under-determined block by variable",items,info.n1,0,col2ele,row2leq,vars,nvar,eqmeta,neqmeta);
        }
        if(info.m3>0) {
          nitems=0;
          for(i=VecSize-info.m3; i<VecSize; i++)items[nitems++]=rowperm[i]-1;
          probe_agg_report(fp,"dm_over_by_eq","DM over-determined block by equation",items,nitems,1,col2ele,row2leq,vars,nvar,eqmeta,neqmeta);
        }
      }
      free(rowperm);
      free(colperm);
      if(fp!=NULL)fprintf(fp,"    \"defective\": true\n  },\n");
    }
    /* fine DM on the realized pattern: SCC structure = irreducible
       simultaneous cores vs recursive structure (opt-in, -probefine) */
    if(pass==1&&probefine) {
      int *rowperm= (int *) malloc (VecSize*sizeof(int));
      int *colperm= (int *) malloc (VecSize*sizeof(int));
      int *rowptr= (int *) malloc ((VecSize+2)*sizeof(int));
      int *colptr= (int *) malloc ((VecSize+2)*sizeof(int));
      gettimeofday(&t0,NULL);
      mc79_fine_i(VecSize,VecSize,ptr,row,rowperm,colperm,rowptr,colptr,&control,&info);
      gettimeofday(&t1,NULL);
      if(info.flag<0)printf("probe: MC79 fine decomposition failed (flag %d)\n",info.flag);
      else {
        PetscInt maxblk=0,nontriv=0,blk;
        offset_t nnz;
        for(k=info.hz_comps; k<info.hz_comps+info.sq_comps; k++) {
          blk=colptr[k+1]-colptr[k];
          if(blk>maxblk)maxblk=blk;
          if(blk>1)nontriv++;
        }
        printf("probe: fine DM (%.2f s): %d strongly connected components in the well-determined block (%ld rows) — %ld simultaneous cores (>1 element), largest %ld, %ld recursive (single-element)\n",(t1.tv_sec-t0.tv_sec)+((double)(t1.tv_usec-t0.tv_usec))/1000000,info.sq_comps,(long)info.m2,(long)nontriv,(long)maxblk,(long)(info.sq_comps-nontriv));
        if(fp!=NULL) {
          fprintf(fp,"  \"fine\": {\"flag\": %d, \"sq_comps\": %d, \"cores_gt1\": %ld, \"largest_core\": %ld,\n",info.flag,info.sq_comps,(long)nontriv,(long)maxblk);
          /* core-size histogram */
          PetscInt *bsz= (PetscInt *) malloc (((info.sq_comps>0)?info.sq_comps:1)*sizeof(PetscInt));
          for(k=0; k<info.sq_comps; k++)bsz[k]=colptr[info.hz_comps+k+1]-colptr[info.hz_comps+k];
          /* sort ascending for run-length emission */
          for(i=1; i<info.sq_comps; i++) {
            PetscInt v=bsz[i];
            PetscInt j2=i-1;
            while(j2>=0&&bsz[j2]>v) {
              bsz[j2+1]=bsz[j2];
              j2--;
            }
            bsz[j2+1]=v;
          }
          fprintf(fp,"    \"core_sizes\": [");
          i=0;
          int first=1;
          while(i<info.sq_comps) {
            PetscInt v=bsz[i],c=0;
            while(i<info.sq_comps&&bsz[i]==v) {
              c++;
              i++;
            }
            fprintf(fp,"%s{\"size\": %ld, \"count\": %ld}",first?"":", ",(long)v,(long)c);
            first=0;
          }
          fprintf(fp,"],\n");
          free(bsz);
          /* composition of the largest nontrivial cores, aggregated by
             statement (rows) and variable (columns) */
          probe_cnt *blocks= (probe_cnt *) malloc (((info.sq_comps>0)?info.sq_comps:1)*sizeof(probe_cnt));
          for(k=0; k<info.sq_comps; k++) {
            blocks[k].idx=info.hz_comps+k;
            blocks[k].c=colptr[info.hz_comps+k+1]-colptr[info.hz_comps+k];
          }
          qsort(blocks,info.sq_comps,sizeof(probe_cnt),probe_cnt_cmp);
          fprintf(fp,"    \"top_cores\": [");
          int ncores=0;
          for(k=0; k<info.sq_comps&&ncores<PROBE_TOP_CORES; k++) {
            if(blocks[k].c<2)break;
            offset_t b=blocks[k].idx;
            fprintf(fp,"%s\n      {\"size\": %ld, \"eqs\": ",(ncores>0)?",":"",blocks[k].c);
            nitems=0;
            for(i=rowptr[b]; i<rowptr[b+1]; i++)items[nitems++]=rowperm[i-1]-1;
            probe_cnt *sorted=probe_agg(items,nitems,1,col2ele,row2leq,vars,nvar,eqmeta,neqmeta,&nnz);
            probe_agg_json(fp,sorted,nnz,1,vars,eqmeta,PROBE_CORE_NAMES_MAX);
            free(sorted);
            fprintf(fp,", \"vars\": ");
            nitems=0;
            for(i=colptr[b]; i<colptr[b+1]; i++)items[nitems++]=colperm[i-1]-1;
            sorted=probe_agg(items,nitems,0,col2ele,row2leq,vars,nvar,eqmeta,neqmeta,&nnz);
            probe_agg_json(fp,sorted,nnz,0,vars,eqmeta,PROBE_CORE_NAMES_MAX);
            free(sorted);
            fprintf(fp,"}");
            ncores++;
          }
          fprintf(fp,"]},\n");
          /* largest-core composition summary to the log */
          if(info.sq_comps>0&&blocks[0].c>1) {
            offset_t b=blocks[0].idx;
            nitems=0;
            for(i=rowptr[b]; i<rowptr[b+1]; i++)items[nitems++]=rowperm[i-1]-1;
            probe_cnt *sorted=probe_agg(items,nitems,1,col2ele,row2leq,vars,nvar,eqmeta,neqmeta,&nnz);
            probe_agg_log("largest core by equation",sorted,nnz,1,vars,eqmeta);
            free(sorted);
          }
          free(blocks);
        }
      }
      free(rowperm);
      free(colperm);
      free(rowptr);
      free(colptr);
    }
    free(ptr);
    free(row);
    ptr=NULL;
    row=NULL;
  }
  if(fp!=NULL) {
    fprintf(fp,"  \"defective\": %s\n",(defects_total>0)?"true":"false");
    fprintf(fp,"}\n");
    fclose(fp);
    logmsg(1,"probe: report written to %s\n",probepath);
  }
  free(items);
  free(rowmatch);
  free(colmatch);
  free(col2ele);
  free(row2leq);
  MatDestroy(&A);
  MatDestroy(&B);
  return (defects_total==0)?0:1;
}
