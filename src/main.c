#include <teems_solver.h>
#include <strings.h>
/* -solmed probe: the structure-detection lines carry the probe prefix so
   the probe log reads as one diagnosis; solve logs are unchanged */
static int teems_probe_mode=0;
#define PROBE_PFX (teems_probe_mode?"probe: ":"")

/* what a run wrote, as one line at verbosity 2: the per-file "Wrote
   <name>" lines were console noise under teems-R */
/* teems-R's condensation appends its substitution coefficients as
   CSUB<n> (Coefficient + Formula pairs); they are derived, not declared
   by the model, so the tally reports them apart */
static int cof_is_csub(const char *name) {
  const char *p=name+4;
  if(strncasecmp(name,"csub",4)!=0||*p=='\0')return 0;
  for(; *p; p++)if(*p<'0'||*p>'9')return 0;
  return 1;
}

static void outputs_summary_log(int cofdumped,array_def *coefs,offset_t ncof,offset_t ncofele,long nsets,long nother,long nskipped,int outputs_on,int solwritten,offset_t nvar,offset_t nvarele) {
  char buf[512];
  int n=0,first=1;
  if(!cofdumped&&!outputs_on&&!solwritten)return;
  n+=snprintf(buf+n,sizeof(buf)-n,"Wrote ");
  if(cofdumped) {
    offset_t i,ncsub=0,ncsubele=0;
    for(i=0; i<ncof; i++)if(cof_is_csub(coefs[i].cofname)) { ncsub++; ncsubele+=coefs[i].nelem; }
    if(ncsub>0)n+=snprintf(buf+n,sizeof(buf)-n,"%ld declared coefficients (%ld elements) and %ld condensation coefficients CSUB* (%ld elements)",(long)(ncof-ncsub),(long)(ncofele-ncsubele),(long)ncsub,(long)ncsubele);
    else n+=snprintf(buf+n,sizeof(buf)-n,"%ld coefficients (%ld elements)",(long)ncof,(long)ncofele);
    first=0;
  }
  if(outputs_on) {
    n+=snprintf(buf+n,sizeof(buf)-n,"%s%ld set%s",first?"":", ",nsets,nsets==1?"":"s"); first=0;
    if(nother>0)n+=snprintf(buf+n,sizeof(buf)-n,", %ld coefficient file%s",nother,nother==1?"":"s");
  }
  if(solwritten)n+=snprintf(buf+n,sizeof(buf)-n,"%s%ld variables (%ld elements)",first?"":", ",(long)nvar,(long)nvarele);
  if(nskipped>0)n+=snprintf(buf+n,sizeof(buf)-n,"; %ld output%s skipped (PostSim pass not run)",nskipped,nskipped==1?"":"s");
  logmsg(2,"%s\n",buf);
}


static char help[] = "teems-solver " TEEMS_SOLVER_VERSION ": solves a CGE model (TABLO/CMF) in parallel.\n\
  -version              print the solver version and exit 0\n\
  -cmdfile <path>       CMF file manifest (default ./reg.cmf)\n\
  -matsol {0,1,2,3}     matrix method LU/SBBD/DBBD/NDBBD\n\
  -solmed <name>        Gragg|Euler|RK2|Heun|RK4|BoSha32|DoPri54|Johansen|probe|nosim\n\
  -jacdump {1,2}        write the base-point Jacobian (<solfiles>.jac + .jac.json); 2 = and stop before the solve\n\
  -step1/-step2/-step3  step counts; -nsubints n; -verbosity {0,1,2}\n\
  Full option table: docs/solver-reference.md section 11.\n\
  Exit status: 0 = completed with no Error line; 1 = an Error line was\n\
  printed (named abort or reported failure); other = crash/kill/MPI.\n\n";

/* Coefficient dump (<stem>.cof + <stem>.cbin): the coefficient twin of
   sol.var + sol.bin, written after PostSim so PostSim coefficients carry
   their computed values.  .cof = int64 header {version, ncof, ncofele,
   reserved} + array_def x ncof (the sol.var struct, binary-locked) +
   uint8 kind x ncof (bit0 postsim, bit1 parameter).  .cbin = solve_real
   x ncofele, the value slice of elem_vals in offset order (file formats
   carry solve_real; store_real is in-memory only).  Replaces the
   per-coefficient CSV dump as teems-R's coefficient transport (ROADMAP
   6.13). */
static int coefficients_dump(const char *stem, array_def *coefs, offset_t ncof, offset_t ncofele, elem_value *elem_vals) {
  char path[TABREADLINE];
  FILE *fp;
  offset_t i,hdr[4];
  strcpy(path,stem);
  strcat(path,".cof");
  if((fp=fopen(path,"wb"))==NULL) {
    errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",path,strerror(errno),(int)getuid());
    return 1;
  }
  hdr[0]=1; hdr[1]=ncof; hdr[2]=ncofele; hdr[3]=0;
  fwrite(hdr,sizeof(offset_t),4,fp);
  for(i=0; i<ncof; i++) {
    /* a p_NAME coefficient runs as p@NAME (pcoef_rename): report the
       declared name */
    array_def rec=coefs[i];
    if(rec.cofname[0]=='p'&&rec.cofname[1]=='@')rec.cofname[1]='_';
    fwrite(&rec,sizeof(array_def),1,fp);
  }
  for(i=0; i<ncof; i++) {
    unsigned char kind=0;
    if(teems_coef_is_ps!=NULL&&teems_coef_is_ps[i])kind|=1;
    if(teems_coef_is_param!=NULL&&teems_coef_is_param[i])kind|=2;
    if(teems_sui&&teems_coef_is_fini!=NULL&&teems_coef_is_fini[i])kind|=4;
    fwrite(&kind,1,1,fp);
  }
  fclose(fp);
  outputs_note(path,"coefficient_declarations",1);
  strcpy(path,stem);
  strcat(path,".cbin");
  if((fp=fopen(path,"wb"))==NULL) {
    errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",path,strerror(errno),(int)getuid());
    return 1;
  }
  {
    /* stream through a fixed buffer: elem_value is a struct, so the value
       field is not contiguous */
    enum { CHUNK=1<<14 };
    solve_real buf[CHUNK];
    offset_t k=0,n;
    while(k<ncofele) {
      n=ncofele-k; if(n>CHUNK)n=CHUNK;
      for(i=0; i<n; i++)buf[i]=(solve_real)elem_vals[k+i].value;
      fwrite(buf,sizeof(solve_real),n,fp);
      k+=n;
    }
  }
  fclose(fp);
  outputs_note(path,"coefficients",0);
  return 0;
}

/* -ma48_cntl2 pivot threshold as given (<=0 = library defaults) */
static double teems_ma48_cntl2_opt=-1.0;
/* -ma48_cntl4: MA48 CNTL(4) in the DBBD/NDBBD rank probes */
static double teems_ma48_cntl4=1e-4;

/* Effective BLAS kernel family. The runtime image pins OPENBLAS_CORETYPE
   (docker/expedited_build/Dockerfile) so the same image gives the same
   bits on every host it supports, but OpenBLAS ignores a core name that
   is absent from its dispatch table WITHOUT a diagnostic -- so the run
   records what the library actually selected, never what was asked for.
   The symbol is weak: a build against a non-OpenBLAS BLAS, and the fuzz
   drivers' stub cluster, link and report "unknown". */
extern char *openblas_get_corename(void) __attribute__((weak));

static const char *blas_corename(void) {
  static const char *cached=NULL;
  if(cached==NULL) {
    const char *n=openblas_get_corename?openblas_get_corename():NULL;
    cached=(n!=NULL&&n[0]!='\0')?n:"unknown";
  }
  return cached;
}

/* Command-line options are strict (Tier B1). PETSc keeps every -name it
   is given and nothing checked for leftovers, so a flag the image does
   not know (a newer teems driving an older image, or a typo) ran as if
   it had not been given and the run exited 0. Every option must be one
   the solver reads, a legacy flag that an earlier teems sent and the
   solver has since dropped (accepted and ignored), or a PETSc runtime
   option. The check runs once, before any work. */
static const char *const cli_teems_flags[]={
  "adaptive","assertions","cmdfile","cofdump",
  "comp_do_acc","comp_do_approx","comp_redo","comp_redo_min_frac",
  "comp_sberr_warn","comp_steps","condest","convrule","epstol","fastrefac","fhtest",
  "inmemory","jacdump","laA","laD","laDi","ma48_cntl2","ma48_cntl4","matsol",
  "maxretries","maxthreads","ndcutcache","nowrites","nsbbdblocks",
  "nsubints","postsim","probefine","probepattern","random_seed",
  "range_test_initial","range_test_updated","refine","residcheck","retryadj","rkchart",
  "rkctrl","rkguard","rk_h0","rknorm","rkscope","single_run","sui","sup","two_run",
  "smllthreads","solmed","step1","step2","step3","tempdir","verbosity","zdivshift",
  "withmc66","version",NULL
};
/* sent by earlier teems releases or read by earlier solvers, no longer
   read: -regset, -enable_time, -presol, -nesteddbbd, -ndbbd_bl_rank
   (teems-R), -enable_iter, -isLinux, -medthreads, -nestfile, -no_hsl,
   -stoiter (solver) */
static const char *const cli_legacy_flags[]={
  "regset","enable_time","presol","nesteddbbd","ndbbd_bl_rank",
  "enable_iter","islinux","medthreads","nestfile","no_hsl","stoiter",NULL
};
/* PETSc runtime options a run may carry (the ASan/TSan gates pass
   -no_signal_handler; teems-R always passes -nox) */
static const char *const cli_petsc_flags[]={
  "nox","nox_warning","display","help","h","options_left","options_view",
  "options_monitor","options_file","skip_petscrc","log_view","log_trace",
  "info","malloc_debug","malloc_dump","malloc_view","malloc_test",
  "memory_view","on_error_abort","on_error_attach_debugger",
  "start_in_debugger","debugger_ranks","stop_for_debugger",
  "no_signal_handler","fp_trap","check_pointer_intensity",
  "mpi_linebuffer","objects_dump","history",NULL
};

static int cli_in(const char *const *list,const char *name) {
  int i;
  for(i=0; list[i]!=NULL; i++)if(strcasecmp(list[i],name)==0)return 1;
  return 0;
}

static int cli_options_check(PetscInt rank) {
  PetscInt n=0,i;
  char **names=NULL,**values=NULL;
  int nbad=0;
  if(PetscOptionsLeftGet(NULL,&n,&names,&values)!=0)return 0;
  for(i=0; i<n; i++) {
    const char *nm=names[i];
    if(nm==NULL)continue;
    while(*nm=='-')nm++;
    if(cli_in(cli_teems_flags,nm)||cli_in(cli_petsc_flags,nm))continue;
    if(cli_in(cli_legacy_flags,nm)) {
      if(rank==0)logmsg(2,"legacy option -%s accepted and ignored\n",nm);
      continue;
    }
    if(rank==0)errmsg("Error: unknown command-line option -%s: teems-solver %s does not read it (teems and its solver image are released together; use the image that matches this teems version, and see docs/solver-reference.md section 11 for the options)\n",nm,TEEMS_SOLVER_VERSION);
    nbad++;
  }
  PetscOptionsLeftRestore(NULL,&n,&names,&values);
  return nbad;
}

/* Side-cars that a later reader could take for this run's (Tier B2):
   removed before any work, so a run that fails leaves none of them and
   no completion marker behind. The locked five and the .cof/.cbin dump
   keep their existing rules. */
static const char *const sidecar_exts[]={".outputs.json",".outputs.json.tmp",".cols",".cols.json",".cbin0",".xac",".ud5",".ud6",".ud7",".jac",".jac.json",NULL};

static void sidecars_clear_stale(const char *cmfname) {
  FILE *f;
  char line[TABREADLINE],*a,*b,*c,*d,path[TABREADLINE+32];
  int i;
  strcpy(teems_sol_stem,"solution");
  f=fopen(cmfname,"r");
  if(f!=NULL) {
    while(tab_next_statement("soldata",f,line,TABREADLINE)) {
      a=strchr(line,'"');
      b=a?strchr(a+1,'"'):NULL;
      c=b?strchr(b+1,'"'):NULL;
      d=c?strchr(c+1,'"'):NULL;
      if(d==NULL||b-a-1!=8||strncasecmp(a+1,"solfiles",8)!=0)continue;
      if((size_t)(d-c-1)>=sizeof(teems_sol_stem))continue;
      memcpy(teems_sol_stem,c+1,d-c-1);
      teems_sol_stem[d-c-1]='\0';
      break;
    }
    fclose(f);
  }
  for(i=0; sidecar_exts[i]!=NULL; i++) {
    snprintf(path,sizeof(path),"%s%s",teems_sol_stem,sidecar_exts[i]);
    remove(path);
  }
}

/* Pre-simulation coefficient values (<stem>.cbin0, Tier B2): the values
   after the initial Reads and Formulas, before any step. Same layout as
   the post-simulation .cbin, behind the .cof header form: int64
   {version 1, ncof, ncofele, phase 0 = pre-simulation}, then ncofele
   doubles in begadd order; the declarations are the run's .cof. */
int coefficients_dump_phase(const char *stem, const char *ext, offset_t phase, offset_t ncof, offset_t ncofele, elem_value *elem_vals) {
  char path[TABREADLINE];
  FILE *fp;
  offset_t i,hdr[4];
  enum { CHUNK=1<<14 };
  solve_real buf[CHUNK];
  offset_t k=0,n;
  snprintf(path,sizeof(path),"%s%s",stem,ext);
  if((fp=fopen(path,"wb"))==NULL) {
    errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",path,strerror(errno),(int)getuid());
    return 1;
  }
  hdr[0]=1; hdr[1]=ncof; hdr[2]=ncofele; hdr[3]=phase;
  fwrite(hdr,sizeof(offset_t),4,fp);
  while(k<ncofele) {
    n=ncofele-k; if(n>CHUNK)n=CHUNK;
    for(i=0; i<n; i++)buf[i]=(solve_real)elem_vals[k+i].value;
    fwrite(buf,sizeof(solve_real),n,fp);
    k+=n;
  }
  if(fclose(fp)!=0) {
    errmsg("Error: cannot write %s: %s\n",path,strerror(errno));
    return 1;
  }
  outputs_note(path,phase==0?"coefficients_presim":"coefficients_updated_pass",1);
  return 0;
}

/* Extra solution columns (<stem>.cols + <stem>.cols.json, Tier B2): one
   container for every per-element column beyond the solution itself.
   int64 header {version 1, ncol, nrow, kind}, nrow int64 row offsets
   (element offsets in .bin order), then ncol x nrow doubles, column by
   column. kind: 0 mixed (per column in the JSON), 1 subtotal,
   2 sagem_individual, 3 approx_cumulative, 4 pass_solution. Written
   for the complementarity approximate run (kind 3) and the shock-group
   subtotals (kind 1; kind 2 for a one-step Johansen run, where one
   group per shock is GEMPACK's SAGEM individual column). The JSON
   "shocked" field is the group's item list (subtotals) or "all". */
static const char *const cols_kind_names[]={"mixed","subtotal","sagem_individual","approx_cumulative","pass_solution"};

static void json_string(FILE *fp,const char *s) {
  fputc('"',fp);
  for(; *s!='\0'; s++) {
    unsigned char c=(unsigned char)*s;
    if(c=='"'||c=='\\')fprintf(fp,"\\%c",c);
    else if(c<0x20)fprintf(fp,"\\u%04x",c);
    else fputc(c,fp);
  }
  fputc('"',fp);
}

static int cols_write(const char *stem, offset_t kind, offset_t ncol, offset_t nrow, const solve_real *const *cols, const char *const *labels, const char *const *shocked) {
  char path[TABREADLINE+16];
  FILE *fp;
  offset_t hdr[4],i,c;
  snprintf(path,sizeof(path),"%s.cols",stem);
  if((fp=fopen(path,"wb"))==NULL) {
    errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",path,strerror(errno),(int)getuid());
    return 1;
  }
  hdr[0]=1; hdr[1]=ncol; hdr[2]=nrow; hdr[3]=kind;
  fwrite(hdr,sizeof(offset_t),4,fp);
  for(i=0; i<nrow; i++)fwrite(&i,sizeof(offset_t),1,fp);
  for(c=0; c<ncol; c++)fwrite(cols[c],sizeof(solve_real),nrow,fp);
  if(fclose(fp)!=0) {
    errmsg("Error: cannot write %s: %s\n",path,strerror(errno));
    return 1;
  }
  outputs_note(path,"cols",1);
  snprintf(path,sizeof(path),"%s.cols.json",stem);
  if((fp=fopen(path,"w"))==NULL) {
    errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",path,strerror(errno),(int)getuid());
    return 1;
  }
  fprintf(fp,"{\n  \"version\": 1,\n  \"ncol\": %ld,\n  \"nrow\": %ld,\n  \"kind\": \"%s\",\n  \"rows\": \"all\",\n  \"columns\": [",(long)ncol,(long)nrow,cols_kind_names[kind]);
  for(c=0; c<ncol; c++) {
    fprintf(fp,"%s\n    {\"index\": %ld, \"kind\": \"%s\", \"label\": ",c?",":"",(long)c,cols_kind_names[kind]);
    json_string(fp,labels[c]);
    fprintf(fp,", \"shocked\": ");
    json_string(fp,shocked!=NULL?shocked[c]:"all");
    fprintf(fp,"}");
  }
  fprintf(fp,"\n  ]\n}\n");
  if(fclose(fp)!=0) {
    errmsg("Error: cannot write %s: %s\n",path,strerror(errno));
    return 1;
  }
  outputs_note(path,"cols_index",1);
  return 0;
}

/* The model-description files every run leaves next to its solution
   (the locked .var/.set/.sel/.mds of teems-R's parse_solution.cpp):
   variable declarations, set declarations, set elements and the four
   counts. Written by the solve path and by the no-simulation path. */
static int structure_files_write(const char *stem, array_def *vars, offset_t nvar, set_def *sets, dim_t nset, set_element *set_elems, offset_t nsetspace, offset_t nvarele) {
  char solchar[TABREADLINE+8];
  FILE *solution;
  const char *ext[3]={".var",".set",".sel"};
  const char *kind[3]={"variable_declarations","set_declarations","set_elements"};
  offset_t f;
  /* a p_NAME linear variable runs as p@NAME (pcoef_rename): the
     declarations carry the declared name */
  array_def *vout=(array_def *)malloc((nvar>0?nvar:1)*sizeof(array_def));
  if (vout==NULL) return 1;
  memcpy(vout,vars,(size_t)nvar*sizeof(array_def));
  for (f=0; f<nvar; f++) if (vout[f].cofname[0]=='p'&&vout[f].cofname[1]=='@') vout[f].cofname[1]='_';
  const void *data[3]={vout,sets,set_elems};
  size_t size[3]={sizeof(array_def),sizeof(set_def),sizeof(set_element)};
  size_t count[3]={(size_t)nvar,(size_t)nset,(size_t)nsetspace};
  offset_t modeldes[4];
  for(f=0; f<3; f++) {
    snprintf(solchar,sizeof(solchar),"%s%s",stem,ext[f]);
    logmsg(2,"solchar %s\n",solchar);
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",solchar,strerror(errno),(int)getuid());
      free(vout);
      return 1;
    }
    fwrite(data[f],size[f],count[f],solution);
    fclose(solution);
    outputs_note(solchar,kind[f],0);
  }
  free(vout);
  modeldes[0]=nsetspace;
  modeldes[1]=nvar;
  modeldes[2]=nvarele;
  modeldes[3]=(offset_t)nset;
  snprintf(solchar,sizeof(solchar),"%s.mds",stem);
  if ( (solution = fopen(solchar, "wb")) == NULL ) {
    errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",solchar,strerror(errno),(int)getuid());
    return 1;
  }
  fwrite(modeldes, sizeof(offset_t),4, solution);
  fclose(solution);
  outputs_note(solchar,"model_description",0);
  return 0;
}

/* Per-run ordering statistics (<solfiles>.stats.json): netcut, border
   sizes and per-block variable/equation counts. Written before the
   solve so failed runs still record their ordering; consumed by
   teems-R for matrix_method calibration. Fields are null when no
   bordered ordering was built (no chain/partition dimension).
   chain_source/partition_source say how each dimension was determined
   ("structural" detection or "none"); auto_json
   carries the partition-candidate table (NULL when nothing was
   probed).
   The "options" object is the run's EFFECTIVE-configuration record
   (posterity/reproducibility; the CMF stays a file manifest by
   design): resolved values after defaults, validation and forced
   changes -- e.g. the complementarity step-sum default and the
   fastrefac force-clear -- not merely what the caller passed.
   teems-R renders it into model_diagnostics.txt after each solve. */
typedef struct {
  int postsim_on;
  int cofdump;
  long subints;
  int adaptive;                 /* 0 fixed, 1 retry-on-check, 2 accuracy-only */
  double epstol;
  int rk_chart, rk_norm, rk_ctrl, rk_scope;
  double rk_h0;
  int maxretries;
  double retryadj;
  long laA, laDi, laD;
  int comp_steps, comp_redo, comp_do_approx, comp_do_acc, comp_sberr_warn;
  int comp_pass;                /* 0 approximate, 1 accurate (write of that pass) */
  double comp_minfrac;
} stats_run_options;

static void ordering_stats_write(cmf_file_entry *iodata, int niodata, int noutdata, int nsoldata,
                                 long VecSize, offset_t nvarele, offset_t nexo,
                                 dim_t matsol, char *solmed, dim_t nesteddbbd, long mpisize, dim_t mc66,
                                 offset_t alltimeset, offset_t allregset, set_def *sets,
                                 offset_t ntime, offset_t nreg, offset_t ndblock,
                                 offset_t netcut, offset_t nintraeq,
                                 offset_t *countvarintra1, offset_t *counteqnoadd,
                                 const char *chain_source, const char *partition_source, const char *auto_json,
                                 const stats_run_options *ropt) {
  static const char *matsol_names[] = {"LU","SBBD","DBBD","NDBBD"};
  char statspath[TABREADLINE+16];
  int i;
  offset_t j;
  for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
    if (strcmp("solfiles",iodata[i].logname)==0)break;
  }
  if(i<niodata+noutdata+nsoldata)strcpy(statspath,iodata[i].filname);
  else strcpy(statspath,"solution");
  strcat(statspath,".stats.json");
  FILE *fp=fopen(statspath,"w");
  if (fp==NULL) {
    printf("Warning: cannot write ordering stats %s\n",statspath);
    return;
  }
  outputs_note(statspath,"stats",2);
  bool bordered=alltimeset>=0||allregset>=0;
  fprintf(fp,"{\n");
  fprintf(fp,"  \"version\": 2,\n");
  fprintf(fp,"  \"solver_version\": \"%s\",\n",TEEMS_SOLVER_VERSION);
  fprintf(fp,"  \"vecsize\": %ld,\n",VecSize);
  fprintf(fp,"  \"nvarele\": %ld,\n",nvarele);
  fprintf(fp,"  \"nexo\": %ld,\n",nexo);
  /* condensation fields only when active, so uncondensed runs stay
     byte-identical to the version-2 layout */
  if(nbacksolve>0)fprintf(fp,"  \"nbacksolve\": %d,\n  \"nbselems\": %ld,\n",nbacksolve,nbselems);
  fprintf(fp,"  \"matrix_method\": \"%s\",\n",(matsol>=MM_LU&&matsol<=MM_NDBBD)?matsol_names[matsol]:"unknown");
  fprintf(fp,"  \"solution_method\": \"%s\",\n",solmed);
  fprintf(fp,"  \"nested_dbbd\": %d,\n",(int)nesteddbbd);
  fprintf(fp,"  \"mc66\": %d,\n",(int)mc66);
  fprintf(fp,"  \"mpi_size\": %ld,\n",mpisize);
  fprintf(fp,"  \"bordered\": %s,\n",bordered?"true":"false");
  fprintf(fp,"  \"chain_source\": \"%s\",\n",chain_source);
  fprintf(fp,"  \"partition_source\": \"%s\",\n",partition_source);
  if(alltimeset>=0)fprintf(fp,"  \"chain_set\": \"%s\",\n  \"time_set\": \"%s\",\n  \"ntime\": %ld,\n",sets[alltimeset].setname,sets[alltimeset].setname,ntime);
  else fprintf(fp,"  \"chain_set\": null,\n  \"time_set\": null,\n  \"ntime\": null,\n");
  if(allregset>=0)fprintf(fp,"  \"partition_set\": \"%s\",\n  \"reg_set\": \"%s\",\n  \"nreg\": %ld,\n",sets[allregset].setname,sets[allregset].setname,nreg);
  else fprintf(fp,"  \"partition_set\": null,\n  \"reg_set\": null,\n  \"nreg\": null,\n");
  if(auto_json!=NULL)fprintf(fp,"%s",auto_json);
  fprintf(fp,"  \"ndblock\": %ld,\n",ndblock);
  if(bordered) {
    fprintf(fp,"  \"netcut\": %ld,\n",netcut);
    fprintf(fp,"  \"nintraeq\": %ld,\n",nintraeq);
    fprintf(fp,"  \"border_neq\": %ld,\n",VecSize-nintraeq);
  }
  else {
    fprintf(fp,"  \"netcut\": null,\n  \"nintraeq\": null,\n  \"border_neq\": null,\n");
  }
  fprintf(fp,"  \"block_nvar\": [");
  for (j=0; j<ndblock; j++)fprintf(fp,"%s%ld",j?",":"",countvarintra1[j+1]-countvarintra1[j]);
  fprintf(fp,"],\n  \"block_neq\": [");
  for (j=0; j<ndblock; j++)fprintf(fp,"%s%ld",j?",":"",counteqnoadd[j]);
  fprintf(fp,"],\n");
  /* effective run-configuration record (see the header comment) */
  {
    static const char *mode_names[] = {"off","warn","fatal"};
    dim_t frchk=0;
    int isrk=(strcmp(solmed,"RK2")==0||strcmp(solmed,"Heun")==0||strcmp(solmed,"RK4")==0||strcmp(solmed,"BoSha32")==0||strcmp(solmed,"DoPri54")==0);
    PetscOptionsGetInt(NULL,NULL,"-fastrefac",&frchk,NULL); /* post force-clear = effective */
    fprintf(fp,"  \"options\": {\n");
    int multistep=(strcmp(solmed,"Gragg")==0||strcmp(solmed,"Midpoint")==0||strcmp(solmed,"Euler")==0);
    if(multistep&&teems_single_run)
      fprintf(fp,"    \"steps\": [%d],\n",steps1);
    else if(multistep&&teems_two_run)
      fprintf(fp,"    \"steps\": [%d,%d],\n",steps1,(int)llround(steps1*step_ratio2));
    else if(multistep)
      fprintf(fp,"    \"steps\": [%d,%d,%d],\n",steps1,(int)llround(steps1*step_ratio2),(int)llround(steps1*step_ratio3));
    else if(strcmp(solmed,"Johansen")==0||strcmp(solmed,"probe")==0)
      fprintf(fp,"    \"steps\": null,\n");
    else fprintf(fp,"    \"steps\": [%d],\n",steps1);
    fprintf(fp,"    \"subintervals\": %ld,\n",ropt->subints);
    fprintf(fp,"    \"single_run\": %s,\n",teems_single_run?"true":"false");
    {
      dim_t cr=0;
      PetscOptionsGetInt(NULL,NULL,"-convrule",&cr,NULL);
      fprintf(fp,"    \"convergence_rule\": %s,\n",cr?"true":"false");
    }
    fprintf(fp,"    \"save_updated_passes\": \"%s\",\n",teems_sup==2?"all":(teems_sup==1?"last":"none"));
    fprintf(fp,"    \"random_seed\": %ld,\n",teems_random_seed);
    if(isrk)fprintf(fp,"    \"adaptive\": %d,\n    \"eps_tolerance\": %g,\n    \"max_retries\": %d,\n    \"retry_adjust\": %g,\n    \"rk_chart\": \"%s\",\n    \"rk_norm\": \"%s\",\n    \"rk_controller\": \"%s\",\n    \"rk_scope\": \"%s\",\n    \"rk_h0\": %g,\n",ropt->adaptive,ropt->epstol,ropt->maxretries,ropt->retryadj,ropt->rk_chart==RK_CHART_LOG?"log":"percent",ropt->rk_norm==RK_NORM_RMS?"rms":"max",ropt->rk_ctrl==RK_CTRL_PI?"pi":"std",ropt->rk_scope==RK_SCOPE_ALL?"all":"pct",ropt->rk_h0);
    else fprintf(fp,"    \"adaptive\": null,\n    \"eps_tolerance\": null,\n    \"max_retries\": null,\n    \"retry_adjust\": null,\n    \"rk_chart\": null,\n    \"rk_norm\": null,\n    \"rk_controller\": null,\n    \"rk_scope\": null,\n    \"rk_h0\": null,\n");
    fprintf(fp,"    \"laA\": %ld,\n    \"laDi\": %ld,\n    \"laD\": %ld,\n",ropt->laA,ropt->laDi,ropt->laD);
    fprintf(fp,"    \"max_threads\": %d,\n",(int)max_threads);
    fprintf(fp,"    \"store_precision\": \"%s\",\n",TEEMS_STORE_PRECISION);
    fprintf(fp,"    \"blas_core\": \"%s\",\n",blas_corename());
    fprintf(fp,"    \"fastrefac\": %s,\n",frchk?"true":"false");
    fprintf(fp,"    \"refine\": %s,\n",(teems_refine&&matsol==MM_DBBD)?"true":"false");
    if(teems_ma48_cntl2_opt>0)fprintf(fp,"    \"ma48_cntl2\": %g,\n",teems_ma48_cntl2_opt);
    else fprintf(fp,"    \"ma48_cntl2\": null,\n");
    fprintf(fp,"    \"ma48_cntl4\": %g,\n",teems_ma48_cntl4);
    fprintf(fp,"    \"ndcutcache\": %d,\n",teems_ndcutcache);
    fprintf(fp,"    \"condest\": %s,\n",teems_condest?"true":"false");
    fprintf(fp,"    \"assertions\": \"%s\",\n",mode_names[teems_assertions_mode>=0&&teems_assertions_mode<=2?teems_assertions_mode:2]);
    fprintf(fp,"    \"range_test_initial\": \"%s\",\n",mode_names[teems_range_test_initial>=0&&teems_range_test_initial<=2?teems_range_test_initial:1]);
    fprintf(fp,"    \"range_test_updated\": \"%s\",\n",mode_names[teems_range_test_updated>=0&&teems_range_test_updated<=2?teems_range_test_updated:1]);
    fprintf(fp,"    \"postsim\": %s,\n",ropt->postsim_on?"true":"false");
    fprintf(fp,"    \"cofdump\": %s",ropt->cofdump?"true":"false");
    if(teems_comp_active>0) {
      fprintf(fp,",\n    \"complementarity\": {\n");
      fprintf(fp,"      \"active_components\": %ld,\n",(long)teems_comp_active);
      fprintf(fp,"      \"steps_approx_run\": %d,\n",ropt->comp_steps);
      fprintf(fp,"      \"redo_steps\": %s,\n",ropt->comp_redo?"true":"false");
      fprintf(fp,"      \"redo_step_min_fraction\": %g,\n",ropt->comp_minfrac);
      fprintf(fp,"      \"do_approx_run\": %s,\n",ropt->comp_do_approx?"true":"false");
      fprintf(fp,"      \"do_acc_run\": %s,\n",ropt->comp_do_acc?"true":"false");
      fprintf(fp,"      \"state_bound_error\": \"%s\",\n",ropt->comp_sberr_warn?"warn":"fatal");
      fprintf(fp,"      \"pass\": \"%s\"\n",ropt->comp_pass?"accurate":"approximate");
      fprintf(fp,"    }\n");
    }
    else fprintf(fp,"\n");
    fprintf(fp,"  }\n");
  }
  fprintf(fp,"}\n");
  fclose(fp);
}

/* Post-solve companion to ordering_stats_write: MA48 workspace growth
   during the solve can raise the effective -la* sizes above the
   recorded options, so patch a top-level "la_used" object (equivalent
   percents; the configured value when nothing grew) into the existing
   stats.json.  Feeds warm-start -la* defaults. */
static void stats_la_used_patch(cmf_file_entry *iodata, int niodata, int noutdata, int nsoldata,
                                long laA_used, long laDi_used, long laD_used) {
  char statspath[TABREADLINE+16];
  int i;
  for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
    if (strcmp("solfiles",iodata[i].logname)==0)break;
  }
  if(i<niodata+noutdata+nsoldata)strcpy(statspath,iodata[i].filname);
  else strcpy(statspath,"solution");
  strcat(statspath,".stats.json");
  FILE *fp=teems_fopen_opt(statspath,"r");
  if (fp==NULL)return;
  fseek(fp,0,SEEK_END);
  long len=ftell(fp);
  fseek(fp,0,SEEK_SET);
  char *buf=(char *) malloc (len+1);
  size_t rd=fread(buf,1,len,fp);
  fclose(fp);
  if((long)rd!=len) {
    free(buf);
    return;
  }
  buf[len]='\0';
  char *end=strrchr(buf,'}');
  if(end==NULL) {
    free(buf);
    return;
  }
  *end='\0';
  fp=fopen(statspath,"w");
  if (fp==NULL) {
    free(buf);
    return;
  }
  fprintf(fp,"%s,\n  \"la_used\": {\"laA\": %ld, \"laDi\": %ld, \"laD\": %ld}\n}\n",buf,laA_used,laDi_used,laD_used);
  fclose(fp);
  free(buf);
}

/* Runge-Kutta run record: stage-solve economy and the retry census,
   patched into stats.json like la_used */
static void stats_rk_patch(cmf_file_entry *iodata, int niodata, int noutdata, int nsoldata) {
  char statspath[TABREADLINE+16];
  int i;
  for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
    if (strcmp("solfiles",iodata[i].logname)==0)break;
  }
  if(i<niodata+noutdata+nsoldata)strcpy(statspath,iodata[i].filname);
  else strcpy(statspath,"solution");
  strcat(statspath,".stats.json");
  FILE *fp=teems_fopen_opt(statspath,"r");
  if (fp==NULL)return;
  fseek(fp,0,SEEK_END);
  long len=ftell(fp);
  fseek(fp,0,SEEK_SET);
  char *buf=(char *) malloc (len+1);
  size_t rd=fread(buf,1,len,fp);
  fclose(fp);
  if((long)rd!=len) {
    free(buf);
    return;
  }
  buf[len]='\0';
  char *end=strrchr(buf,'}');
  if(end==NULL) {
    free(buf);
    return;
  }
  *end='\0';
  fp=fopen(statspath,"w");
  if (fp==NULL) {
    free(buf);
    return;
  }
  fprintf(fp,"%s,\n  \"runge_kutta\": {\"steps\": %ld, \"stage_solves\": %ld, \"stage_solves_reused\": %ld, \"rejects_accuracy\": %ld, \"rejects_crossed\": %ld, \"rejects_range\": %ld, \"rejects_assertion\": %ld, \"rejects_guard\": %ld, \"rejects_singular\": %ld, \"h_min\": %g, \"h_max\": %g, \"worst_step_metric\": %g, \"worst_estimated_metric\": %g}\n}\n",
          buf,teems_rk_stats.steps,teems_rk_stats.stage_solves,teems_rk_stats.stage_solves_reused,
          teems_rk_stats.rejects_accuracy,teems_rk_stats.rejects_crossed,teems_rk_stats.rejects_range,
          teems_rk_stats.rejects_assert,teems_rk_stats.rejects_guard,teems_rk_stats.rejects_singular,
          teems_rk_stats.h_min,teems_rk_stats.h_max,teems_rk_stats.worst_step_metric,teems_rk_stats.worst_est_metric);
  fclose(fp);
  free(buf);
}

/* NDBBD thread-budget record: the team size each capped region ran
   with on the last step (-maxthreads / -smllthreads asked, the budget
   chose), patched into stats.json like la_used; absent unless NDBBD ran */
static void stats_ndbbd_threads_patch(cmf_file_entry *iodata, int niodata, int noutdata, int nsoldata) {
  char statspath[TABREADLINE+16];
  int i;
  if(ndbbd_threads_used[0]<=0)return;
  for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
    if (strcmp("solfiles",iodata[i].logname)==0)break;
  }
  if(i<niodata+noutdata+nsoldata)strcpy(statspath,iodata[i].filname);
  else strcpy(statspath,"solution");
  strcat(statspath,".stats.json");
  FILE *fp=teems_fopen_opt(statspath,"r");
  if (fp==NULL)return;
  fseek(fp,0,SEEK_END);
  long len=ftell(fp);
  fseek(fp,0,SEEK_SET);
  char *buf=(char *) malloc (len+1);
  size_t rd=fread(buf,1,len,fp);
  fclose(fp);
  if((long)rd!=len) {
    free(buf);
    return;
  }
  buf[len]='\0';
  char *end=strrchr(buf,'}');
  if(end==NULL) {
    free(buf);
    return;
  }
  *end='\0';
  fp=fopen(statspath,"w");
  if (fp==NULL) {
    free(buf);
    return;
  }
  fprintf(fp,"%s,\n  \"ndbbd_threads\": {\"asked\": %d, \"presolve\": %d, \"interface_rank\": %d, \"interface_factor\": %d, \"schur\": %d}\n}\n",
          buf,(int)max_threads,ndbbd_threads_used[0],ndbbd_threads_used[1],ndbbd_threads_used[2],ndbbd_threads_used[3]);
  fclose(fp);
  free(buf);
}

/* phase resident-memory record (6.16(a)): per phase the max-over-ranks
   and sum-over-ranks resident set in GB and the run's high-water mark,
   patched into stats.json like la_used (same append-at-tail flow) */
static void stats_rss_patch(cmf_file_entry *iodata, int niodata, int noutdata, int nsoldata) {
  char statspath[TABREADLINE+16];
  int i,j;
  double hwm_max=0.0,hwm_sum=0.0;
  for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
    if (strcmp("solfiles",iodata[i].logname)==0)break;
  }
  if(i<niodata+noutdata+nsoldata)strcpy(statspath,iodata[i].filname);
  else strcpy(statspath,"solution");
  strcat(statspath,".stats.json");
  FILE *fp=teems_fopen_opt(statspath,"r");
  if (fp==NULL)return;
  fseek(fp,0,SEEK_END);
  long len=ftell(fp);
  fseek(fp,0,SEEK_SET);
  char *buf=(char *) malloc (len+1);
  size_t rd=fread(buf,1,len,fp);
  fclose(fp);
  if((long)rd!=len) {
    free(buf);
    return;
  }
  buf[len]='\0';
  char *end=strrchr(buf,'}');
  if(end==NULL) {
    free(buf);
    return;
  }
  *end='\0';
  fp=fopen(statspath,"w");
  if (fp==NULL) {
    free(buf);
    return;
  }
  fprintf(fp,"%s,\n  \"rss_gb\": {\n",buf);
  for(i=0; i<teems_nrss; i++) {
    fprintf(fp,"    \"");
    for(j=0; teems_rss[i].phase[j]!='\0'; j++)fputc(teems_rss[i].phase[j]==' '?'_':teems_rss[i].phase[j],fp);
    fprintf(fp,"\": {\"max\": %.3f, \"sum\": %.3f, \"probes\": %ld},\n",teems_rss[i].rss_max,teems_rss[i].rss_sum,teems_rss[i].count);
    if(teems_rss[i].hwm_max>hwm_max)hwm_max=teems_rss[i].hwm_max;
    if(teems_rss[i].hwm_sum>hwm_sum)hwm_sum=teems_rss[i].hwm_sum;
  }
  fprintf(fp,"    \"peak\": {\"max\": %.3f, \"sum\": %.3f}\n  }\n}\n",hwm_max,hwm_sum);
  fclose(fp);
  free(buf);
}

/* -condest record: per-run maxima of the MA60 solve-quality measures,
   patched into stats.json like la_used (same append-at-tail flow) */
static void stats_condest_patch(cmf_file_entry *iodata, int niodata, int noutdata, int nsoldata,
                                double kw1max, double kw2max, double omegamax,
                                long solves, long skips) {
  char statspath[TABREADLINE+16];
  int i;
  for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
    if (strcmp("solfiles",iodata[i].logname)==0)break;
  }
  if(i<niodata+noutdata+nsoldata)strcpy(statspath,iodata[i].filname);
  else strcpy(statspath,"solution");
  strcat(statspath,".stats.json");
  FILE *fp=teems_fopen_opt(statspath,"r");
  if (fp==NULL)return;
  fseek(fp,0,SEEK_END);
  long len=ftell(fp);
  fseek(fp,0,SEEK_SET);
  char *buf=(char *) malloc (len+1);
  size_t rd=fread(buf,1,len,fp);
  fclose(fp);
  if((long)rd!=len) {
    free(buf);
    return;
  }
  buf[len]='\0';
  char *end=strrchr(buf,'}');
  if(end==NULL) {
    free(buf);
    return;
  }
  *end='\0';
  fp=fopen(statspath,"w");
  if (fp==NULL) {
    free(buf);
    return;
  }
  fprintf(fp,"%s,\n  \"condest\": {\"kappa_w1_max\": %.6e, \"kappa_w2_max\": %.6e, \"omega_max\": %.6e, \"solves\": %ld, \"zero_rhs_skips\": %ld}\n}\n",
          buf,kw1max,kw2max,omegamax,solves,skips);
  fclose(fp);
  free(buf);
}

/* solve-accuracy record (manual 30.1.5): maximum residual ratio over
   every checked solve, patched into stats.json like la_used; null when
   the matrix method releases A before the solve (no solve checked) */
static void stats_resid_patch(cmf_file_entry *iodata, int niodata, int noutdata, int nsoldata,
                              double rmax, long solves, long warns, long skipped) {
  char statspath[TABREADLINE+16];
  int i;
  for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
    if (strcmp("solfiles",iodata[i].logname)==0)break;
  }
  if(i<niodata+noutdata+nsoldata)strcpy(statspath,iodata[i].filname);
  else strcpy(statspath,"solution");
  strcat(statspath,".stats.json");
  FILE *fp=teems_fopen_opt(statspath,"r");
  if (fp==NULL)return;
  fseek(fp,0,SEEK_END);
  long len=ftell(fp);
  fseek(fp,0,SEEK_SET);
  char *buf=(char *) malloc (len+1);
  size_t rd=fread(buf,1,len,fp);
  fclose(fp);
  if((long)rd!=len) {
    free(buf);
    return;
  }
  buf[len]='\0';
  char *end=strrchr(buf,'}');
  if(end==NULL) {
    free(buf);
    return;
  }
  *end='\0';
  fp=fopen(statspath,"w");
  if (fp==NULL) {
    free(buf);
    return;
  }
  if(solves>0)fprintf(fp,"%s,\n  \"residual\": {\"max_residual_ratio\": %.6e, \"solves\": %ld, \"warnings\": %ld, \"unchecked_solves\": %ld}\n}\n",buf,rmax,solves,warns,skipped);
  else fprintf(fp,"%s,\n  \"residual\": {\"max_residual_ratio\": null, \"solves\": 0, \"warnings\": 0, \"unchecked_solves\": %ld}\n}\n",buf,skipped);
  fclose(fp);
  free(buf);
}

/* DBBD refinement record (-refine): solves refined, the worst residual
   ratio before and after the step, and the seconds it took, patched into
   stats.json like the residual record */
static void stats_refine_patch(cmf_file_entry *iodata, int niodata, int noutdata, int nsoldata,
                               long solves, double before, double after, double secs) {
  char statspath[TABREADLINE+16];
  int i;
  for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
    if (strcmp("solfiles",iodata[i].logname)==0)break;
  }
  if(i<niodata+noutdata+nsoldata)strcpy(statspath,iodata[i].filname);
  else strcpy(statspath,"solution");
  strcat(statspath,".stats.json");
  FILE *fp=teems_fopen_opt(statspath,"r");
  if (fp==NULL)return;
  fseek(fp,0,SEEK_END);
  long len=ftell(fp);
  fseek(fp,0,SEEK_SET);
  char *buf=(char *) malloc (len+1);
  size_t rd=fread(buf,1,len,fp);
  fclose(fp);
  if((long)rd!=len) {
    free(buf);
    return;
  }
  buf[len]='\0';
  char *end=strrchr(buf,'}');
  if(end==NULL) {
    free(buf);
    return;
  }
  *end='\0';
  fp=fopen(statspath,"w");
  if (fp==NULL) {
    free(buf);
    return;
  }
  fprintf(fp,"%s,\n  \"refine\": {\"solves\": %ld, \"residual_ratio_before_max\": %.6e, \"residual_ratio_after_max\": %.6e, \"seconds\": %.3f}\n}\n",buf,solves,before,after,secs);
  fclose(fp);
  free(buf);
}

/* ====================================================================
   Structural partition detection

   The bordered orderings need two inputs: the chain dimension (the set
   whose elements the equations couple through lead/lag index offsets;
   conventionally time) and the diagonal-block partition set
   (conventionally regions). Neither is a property of a set's name or
   declaration: both follow from how the system of equations references
   each set, so the routines below derive them from the equations.
   (The transitional -enable_time/-regset overrides were removed once
   detection reproduced their results bit-identically; stale flags in
   old commands are ignored by the options parser.)
   ==================================================================== */

/* Mark intsup on every set nested inside a qualified chain set, so the
   ordering can map subset elements onto chain-block positions. Extracted
   verbatim from the inline marking (both the explicit-flag path and the
   structural path must produce identical maps). */
static void chain_flags_apply(set_def *sets, dim_t nset, offset_t chain) {
  offset_t i,j;
  for(i=0; i<nset; i++) {
    for(j=1; j<MAXSUPSET; j++) {
      if(sets[i].subsetid[j]>-1) {
        if(sets[i].subsetid[j]==chain) {
          sets[i].intsup=j;
          break;
        }
      }
      else break;
    }
  }
}

/* Mark the chosen partition set and every set nested inside it
   (regsup gives the superset slot whose element positions index the
   diagonal blocks). Extracted verbatim from the inline marking. */
static void partition_flags_apply(set_def *sets, dim_t nset, offset_t partset) {
  offset_t i,j;
  sets[partset].regional=true;
  for(i=0; i<nset; i++) {
    for(j=1; j<MAXSUPSET; j++) {
      if(sets[i].subsetid[j]>-1) {
        if(sets[sets[i].subsetid[j]].regional) {
          sets[i].regional=true;
          sets[i].regsup=j;
          break;
        }
      }
      else break;
    }
  }
}

/* Reset the partition marks between probe candidates (regional/regsup
   are zero from calloc until a partition is applied). */
static void partition_flags_clear(set_def *sets, dim_t nset) {
  offset_t i;
  for(i=0; i<nset; i++) {
    sets[i].regional=false;
    sets[i].regsup=0;
  }
}

/* First ordering pass, shared by the live ordering below and the
   partition probe: count the endogenous variable elements that fall
   inside each diagonal block under the current chain/partition marks
   (border variables/elements and exogenous elements excluded).
   countvar must be zeroed, length ndblock. Extracted verbatim from the
   three inline counting passes. */
static void block_var_count(array_def *vars, offset_t nvar, set_def *sets, set_element *set_elems,
                            closure_entry *closure_vals, bool *var_inter, bool *ele_inter,
                            dim_t *orderintra, dim_t *orderreg,
                            offset_t alltimeset, offset_t allregset, offset_t nreg, dim_t nesteddbbd,
                            offset_t *countvar) {
  offset_t i,j,j0,j1,j2,j4;
  offset_t j3=0;
  if(nesteddbbd==1) {
    for (i=0; i<nvar; i++) {
      for (j=0; j<vars[i].nelem; j++) {
        if(!CL_EXO(j3+j)&&!CL_BS(j3+j)) {
          if(!var_inter[i]&&!ele_inter[j3+j]) {
            j0=j;
            j2=-1;
            for(j1=0; j1<orderintra[i]+1; j1++) {
              j2=j0/vars[i].strides[j1];
              j0-=j2*vars[i].strides[j1];
            }
            j0=j;
            j4=-1;
            for(j1=0; j1<orderreg[i]+1; j1++) {
              j4=j0/vars[i].strides[j1];
              j0-=j4*vars[i].strides[j1];
            }
            if(j4>-1)if(sets[vars[i].setid[orderreg[i]]].regsup>0)j4=set_elems[sets[vars[i].setid[orderreg[i]]].offset+j4].superset_pos[sets[vars[i].setid[orderreg[i]]].regsup];
            if(sets[vars[i].setid[orderintra[i]]].intsup>0)j2=set_elems[sets[vars[i].setid[orderintra[i]]].offset+j2].superset_pos[sets[vars[i].setid[orderintra[i]]].intsup];
            if(orderreg[i]>-1)countvar[j2*(nreg+1)+j4]++;
            else countvar[j2*(nreg+1)+nreg]++;
          }
        }
      }
      j3+=vars[i].nelem;
    }
  }
  else if(alltimeset>=0) {
    for (i=0; i<nvar; i++) {
      for (j=0; j<vars[i].nelem; j++) {
        if(!CL_EXO(j3+j)&&!CL_BS(j3+j)) {
          if(!var_inter[i]&&!ele_inter[j3+j]) {
            j0=j;
            j2=-1;
            for(j1=0; j1<orderintra[i]+1; j1++) {
              j2=j0/vars[i].strides[j1];
              j0-=j2*vars[i].strides[j1];
            }
            if(allregset>=0) {
              j0=j;
              j4=-1;
              for(j1=0; j1<orderreg[i]+1; j1++) {
                j4=j0/vars[i].strides[j1];
                j0-=j4*vars[i].strides[j1];
              }
              if(sets[vars[i].setid[orderreg[i]]].regsup>0)j4=set_elems[sets[vars[i].setid[orderreg[i]]].offset+j4].superset_pos[sets[vars[i].setid[orderreg[i]]].regsup];
              countvar[j2*nreg+j4]++;
            }
            else {
              countvar[j2]++;
            }
          }
        }
      }
      j3+=vars[i].nelem;
    }
  }
  else if(allregset>=0) {
    for (i=0; i<nvar; i++) {
      for (j=0; j<vars[i].nelem; j++) {
        if(!CL_EXO(j3+j)&&!CL_BS(j3+j)) {
          if(!var_inter[i]&&!ele_inter[j3+j]) {
            j0=j;
            j4=-1;
            for(j1=0; j1<orderreg[i]+1; j1++) {
              j4=j0/vars[i].strides[j1];
              j0-=j4*vars[i].strides[j1];
            }
            if(sets[vars[i].setid[orderreg[i]]].regsup>0)j4=set_elems[sets[vars[i].setid[orderreg[i]]].offset+j4].superset_pos[sets[vars[i].setid[orderreg[i]]].regsup];
            countvar[j4]++;
          }
        }
      }
      j3+=vars[i].nelem;
    }
  }
}

/* Structural chain scan: walk the equation statements and count, per
   declared dimension set, the variable references carrying a lead/lag
   index offset (x{t+1} after normalization/encoding). The offsets ARE
   the chain structure; the (intertemporal) qualifier only licenses the
   syntax and is cross-checked by the caller. Statement preprocessing
   and reference matching mirror equation_order_read. */
static int chain_refs_scan(char *fname, set_def *sets, dim_t nset, array_def *coefs, offset_t ncof,
                           array_def *vars, offset_t nvar, elem_value *elem_vals,
                           offset_t *refcount) {
  FILE *filehandle;
  char line[TABREADLINE],linecopy[TABREADLINE],idxlist[TABREADLINE],vname[TABREADLINE];
  char commsyntax[NAMESIZE];
  char *readitem=NULL,*p=NULL,*pend=NULL,*tok=NULL;
  solve_real zerodivide=0;
  dim_t fdim,i4;
  offset_t l;
  int i,np,varindx1,varindx2,leadlag;
  strcpy(commsyntax,"equation");
  filehandle=teems_fopen(fname,"r");
  if(filehandle==NULL)return 0;
  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")!=NULL)continue;
    str_replace_first(line, commsyntax, "");
    str_replace_first(line, "(linear)", "");
    while (str_replace_all(line,"  ", " "));
    while (str_replace_char(line, '[', '('));
    while (str_replace_char(line, ']', ')'));
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
      readitem = strtok(line+1,"(");
      strcpy(line,linecopy);
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
    str_delete_char(readitem,' ');
    str_sign_runs_collapse(readitem);
    while (formula_normalize(readitem)==1);
    leadlag_encode(readitem);
    np=str_count_ci(readitem,"p_");
    for (i=0; i<np; i++) {
      varindx2=0;
      while(-1<0) {
        varindx1=str_find_ci(readitem+varindx2,"p_");
        if(varindx1==-1) break;
        varindx2=varindx2+varindx1;
        if(varindx2==0) {
          /* genuine start only before the first token (section 6) */
          if(i==0) break;
          varindx2++;
        }
        else if(readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
        else varindx2++;
      }
      if(varindx1==-1) break;
      readitem=readitem+varindx2+2;
      p=strpbrk(readitem,"{+*-/^)");
      if(p==NULL||*p!='{')continue;/* reference without indices */
      strncpy(vname,readitem,p-readitem);
      vname[p-readitem]='\0';
      {
        offset_t lr=linvar_resolve(vname,vars,nvar);
        if(lr<0)continue;/* refcount scan stays lenient; the equation
                            builders fatal undeclared references */
        l=lr;
      }
      pend=strchr(p,'}');
      if(pend==NULL)continue;
      strncpy(idxlist,p+1,pend-p-1);
      idxlist[pend-p-1]='\0';
      i4=0;
      tok=strtok(idxlist,",");
      while(tok!=NULL&&i4<vars[l].size) {
        leadlag=0;
        parse_index_leadlag(tok,&leadlag);
        if(leadlag!=0)refcount[vars[l].setid[i4]]++;
        i4++;
        tok=strtok(NULL,",");
      }
    }
  }
  fclose(filehandle);
  return 1;
}

/* Aggregate the per-dimension lead/lag reference counts onto top-level
   sets and pick the chain dimension (most-referenced top-level set; -1
   when the equations use no offsets). The (intertemporal) qualifier is
   cross-checked and mismatches reported; structure decides. */
static offset_t chain_set_select(set_def *sets, dim_t nset, offset_t *refcount, PetscInt rank) {
  offset_t i,top,chain=-1,nchain=0;
  offset_t j;
  offset_t *topcount=(offset_t *) calloc (nset,sizeof(offset_t));
  for(i=0; i<nset; i++) {
    if(refcount[i]==0)continue;
    top=i;
    if(sets[i].subsetid[1]!=-1) {
      for(j=1; j<MAXSUPSET; j++) {
        if(sets[i].subsetid[j]==-1)break;
        if(sets[sets[i].subsetid[j]].subsetid[1]==-1) {
          top=sets[i].subsetid[j];
          break;
        }
      }
    }
    topcount[top]+=refcount[i];
  }
  for(i=0; i<nset; i++) {
    if(topcount[i]>0) {
      nchain++;
      if(chain<0||topcount[i]>topcount[chain])chain=i;
    }
  }
  if(rank==0) {
    if(nchain>1) {
      printf("Warning: lead/lag references span %ld top-level sets (",nchain);
      for(i=0; i<nset; i++)if(topcount[i]>0)printf(" %s:%ld",sets[i].setname,topcount[i]);
      printf(" ); using %s as the chain dimension, the others are ordered as ordinary sets\n",sets[chain].setname);
    }
    if(chain>=0&&!sets[chain].intertemp)
      printf("Warning: equations apply lead/lag offsets to set %s, which lacks the (intertemporal) qualifier\n",sets[chain].setname);
    if(chain<0) {
      i=set_find_alltime(sets,nset);
      if(i>=0)logmsg(1,"Set %s is declared (intertemporal) but no equation uses lead/lag offsets; no chain ordering applied\n",sets[i].setname);
    }
    if(chain>=0&&sets[chain].size<2)logmsg(1,"Set %s carries the lead/lag offsets but has %d element(s); no chain ordering applied (one period orders as a static model)\n",sets[chain].setname,sets[chain].size);
    else if(chain>=0)logmsg(1,"%sChain dimension detected structurally: set %s (size %d), %ld lead/lag references\n",PROBE_PFX,sets[chain].setname,sets[chain].size,topcount[chain]);
  }
  if(chain>=0&&sets[chain].size<2)chain=-1;
  free(topcount);
  return chain;
}

/* Probe one candidate partition set: apply its marks, run the
   pre-Jacobian ordering scan and the first counting pass, and report
   the border size and per-block extremes that partition would produce.
   Everything touched is restored or freed; the run itself is
   unchanged. In nested mode the per-chain-block interface column is
   excluded from the block extremes (it is not a parallel block). */
static int partition_probe(char *tabfile, set_def *sets, dim_t nset, set_element *set_elems,
                           array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar,
                           elem_value *elem_vals, offset_t ncofele, offset_t nvarele,
                           closure_entry *closure_vals, offset_t neq, offset_t VecSize,
                           offset_t alltimeset, offset_t ntime, dim_t nesteddbbd, offset_t cand,
                           offset_t *netcut_out, offset_t *nblocks_out, offset_t *blkmin_out, offset_t *blkmax_out) {
  char commsyntax[NAMESIZE];
  offset_t i,nreg,ndblock,total,bmin,bmax;
  strcpy(commsyntax,"equation");
  nreg=sets[cand].size;
  if(nesteddbbd==1)ndblock=ntime*(nreg+1);
  else if(alltimeset>=0)ndblock=ntime*nreg;
  else ndblock=nreg;
  bool *var_inter=(bool *) calloc (nvar,sizeof(bool));
  bool *ele_inter=(bool *) calloc (nvarele,sizeof(bool));
  dim_t *orderintra=(dim_t *) malloc (nvar*sizeof(dim_t));
  dim_t *orderreg=(dim_t *) malloc (nvar*sizeof(dim_t));
  for(i=0; i<nvar; i++) {
    orderintra[i]=-1;
    orderreg[i]=-1;
  }
  array_def *eq_defs=(array_def *) calloc (neq,sizeof(array_def));
  bool *eq_intertemp=(bool *) calloc (neq,sizeof(bool));
  dim_t *eq_time=(dim_t *) malloc (neq*sizeof(dim_t));
  dim_t *eq_reg=(dim_t *) malloc (neq*sizeof(dim_t));
  for(i=0; i<neq; i++) {
    eq_time[i]=-1;
    eq_reg[i]=-1;
  }
  offset_t *countvar=(offset_t *) calloc (ndblock,sizeof(offset_t));
  partition_flags_apply(sets,nset,cand);
  if(nesteddbbd==1) {
    if(!equation_order_read_nested(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,var_inter,ele_inter,eq_defs,eq_intertemp,eq_time,eq_reg,cand,alltimeset,orderintra,orderreg))MPI_Abort(PETSC_COMM_WORLD,1);
  }
  else if(!equation_order_read(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,var_inter,ele_inter,eq_defs,eq_intertemp,eq_time,eq_reg,cand,alltimeset,orderintra,orderreg))MPI_Abort(PETSC_COMM_WORLD,1);
  block_var_count(vars,nvar,sets,set_elems,closure_vals,var_inter,ele_inter,orderintra,orderreg,alltimeset,cand,nreg,nesteddbbd,countvar);
  partition_flags_clear(sets,nset);
  total=0;
  bmin=-1;
  bmax=0;
  for(i=0; i<ndblock; i++) {
    /* in nested mode the per-chain-block interface column is the local
       border the interface solve pays for — score it as border, not as
       an intra block */
    if(nesteddbbd==1&&i%(nreg+1)==nreg)continue;
    total+=countvar[i];
    if(bmin<0||countvar[i]<bmin)bmin=countvar[i];
    if(countvar[i]>bmax)bmax=countvar[i];
  }
  *netcut_out=VecSize-total;
  *nblocks_out=ndblock;
  *blkmin_out=bmin<0?0:bmin;
  *blkmax_out=bmax;
  free(var_inter);
  free(ele_inter);
  free(orderintra);
  free(orderreg);
  free(eq_defs);
  free(eq_intertemp);
  free(eq_time);
  free(eq_reg);
  free(countvar);
  return 1;
}

/* Structural partition selection: probe every candidate set (two or
   more elements, outside the chain family, indexing at least one
   variable dimension directly or through a subset) and choose the one
   with the smallest border, breaking near-ties (2%) by block balance.
   A candidate is viable when it yields at least n_tasks blocks, a
   border below half the system, and no empty block. Returns the set id
   or -1 when no candidate is viable. Deterministic integer arithmetic
   throughout, so every rank reaches the same answer independently.
   *auto_json receives the stats.json candidate table (malloc'd; caller
   frees). */
static offset_t partition_auto_select(char *tabfile, set_def *sets, dim_t nset, set_element *set_elems,
                                      array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar,
                                      elem_value *elem_vals, offset_t ncofele, offset_t nvarele,
                                      closure_entry *closure_vals, offset_t neq, offset_t VecSize,
                                      offset_t alltimeset, offset_t ntime, dim_t nesteddbbd,
                                      long mpisize, PetscInt rank, char **auto_json) {
  offset_t i,s,d,chosen=-1,bestcut=-1,ncand=0;
  offset_t j;
  bool *isdim=(bool *) calloc (nset,sizeof(bool));
  bool *iscand=(bool *) calloc (nset,sizeof(bool));
  bool *viable=(bool *) calloc (nset,sizeof(bool));
  offset_t *cnetcut=(offset_t *) calloc (nset,sizeof(offset_t));
  offset_t *cnblocks=(offset_t *) calloc (nset,sizeof(offset_t));
  offset_t *cmin=(offset_t *) calloc (nset,sizeof(offset_t));
  offset_t *cmax=(offset_t *) calloc (nset,sizeof(offset_t));
  for(i=0; i<nvar; i++)for(d=0; d<vars[i].size; d++)isdim[vars[i].setid[d]]=true;
  for(s=0; s<nset; s++) {
    if(sets[s].size<2)continue;
    if(sets[s].intertemp)continue;/* chain-family syntax carrier */
    if(alltimeset>=0) {
      if(s==alltimeset)continue;
      for(j=1; j<MAXSUPSET; j++) {
        if(sets[s].subsetid[j]==-1)break;
        if(sets[s].subsetid[j]==alltimeset)break;
      }
      if(j<MAXSUPSET&&sets[s].subsetid[j]==alltimeset)continue;/* nested in the chain set */
    }
    if(!isdim[s]) {/* usable when some variable dimension maps into s */
      for(d=0; d<nset; d++) {
        if(!isdim[d])continue;
        for(j=1; j<MAXSUPSET; j++) {
          if(sets[d].subsetid[j]==-1)break;
          if(sets[d].subsetid[j]==s)break;
        }
        if(j<MAXSUPSET&&sets[d].subsetid[j]==s)break;
      }
      if(d==nset)continue;
    }
    iscand[s]=true;
    ncand++;
  }
  if(rank==0)logmsg(1,"%sPartition detection: probing %ld candidate sets against the equation structure\n",PROBE_PFX,ncand);
  for(s=0; s<nset; s++) {
    if(!iscand[s])continue;
    partition_probe(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,nvarele,closure_vals,neq,VecSize,alltimeset,ntime,nesteddbbd,s,&cnetcut[s],&cnblocks[s],&cmin[s],&cmax[s]);
    viable[s]=cnblocks[s]>=mpisize&&2*cnetcut[s]<VecSize&&cmin[s]>0;
    if(rank==0)logmsg(1,"%s  %-14s blocks %6ld  border %8ld (%5.1f%%)  block min/max %ld/%ld  %s\n",
                        PROBE_PFX,sets[s].setname,cnblocks[s],cnetcut[s],VecSize>0?100.0*cnetcut[s]/VecSize:0.0,
                        cmin[s],cmax[s],viable[s]?"viable":"not viable");
  }
  for(s=0; s<nset; s++) {
    if(!viable[s])continue;
    if(bestcut<0||cnetcut[s]<bestcut)bestcut=cnetcut[s];
  }
  for(s=0; s<nset; s++) {
    if(!viable[s])continue;
    if(50*cnetcut[s]>51*bestcut)continue;/* outside the 2% near-tie band */
    if(chosen<0)chosen=s;
    else if(cmin[s]*cmax[chosen]>cmin[chosen]*cmax[s])chosen=s;/* better balance */
  }
  if(rank==0) {
    if(chosen>=0)logmsg(1,"%sPartition detection: selected set %s (%ld blocks, border %.1f%% of system)\n",
                          PROBE_PFX,sets[chosen].setname,cnblocks[chosen],VecSize>0?100.0*cnetcut[chosen]/VecSize:0.0);
    else if(ncand>0)logmsg(1,"%sPartition detection: no viable candidate (need >=%ld nonempty blocks and a border below half the system)\n",PROBE_PFX,mpisize);
    else logmsg(1,"%sPartition detection: no candidate sets to probe\n",PROBE_PFX);
  }
  {
    size_t cap=256+(size_t)ncand*256;
    char *js=(char *) malloc (cap);
    size_t off=0;
    off+=snprintf(js+off,cap-off,"  \"partition_auto\": {\n    \"candidates\": [");
    j=0;
    for(s=0; s<nset; s++) {
      if(!iscand[s])continue;
      off+=snprintf(js+off,cap-off,"%s\n      {\"set\": \"%s\", \"size\": %d, \"nblocks\": %ld, \"netcut\": %ld, \"block_min\": %ld, \"block_max\": %ld, \"viable\": %s}",
                    j?",":"",sets[s].setname,sets[s].size,cnblocks[s],cnetcut[s],cmin[s],cmax[s],viable[s]?"true":"false");
      j++;
      if(off>=cap-256)break;
    }
    if(chosen>=0)off+=snprintf(js+off,cap-off,"\n    ],\n    \"chosen\": \"%s\"\n  },\n",sets[chosen].setname);
    else off+=snprintf(js+off,cap-off,"\n    ],\n    \"chosen\": null\n  },\n");
    *auto_json=js;
  }
  free(isdim);
  free(iscand);
  free(viable);
  free(cnetcut);
  free(cnblocks);
  free(cmin);
  free(cmax);
  return chosen;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args) {
  Vec vecb,vece;  /* approx solution, RHS, exact solution */
  Mat      A,B;    /* linear system matrix */
  PetscInt  rank=0,mpisize,rank_hsl=0;
  PetscInt VecSize=0,Istart=0,Iend=0,dnz=0,onz=0,dnzB=0,onzB=0,*onnz,*dnnz,*onnzB,*dnnzB;
  PetscErrorCode ierr;
  PetscBool   flg;
  struct timeval begintime,endtime,probe_begin;
  offset_t i,j;
  offset_t j2=0,j1=0,j0=0,j3,j4,j5,j6;
  /* -version answers before MPI/PETSc start: one line, exit 0, no
     PETSc banner (PETSc owns -version too and would print its own).
     A pre-1.1 image gives no such line: no response IS the old-image
     detection the R side relies on. */
  for(int a=1; a<argc; a++) {
    if(strcmp(args[a],"-version")==0||strcmp(args[a],"--version")==0) {
      printf("teems-solver %s\n",TEEMS_SOLVER_VERSION);
      return 0;
    }
  }
  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize);

  char processor_name[MPI_MAX_PROCESSOR_NAME+1];
  int name_len,name_len_max,name_beg,class_size,color=0,group_size,node_rank;  /* color=0 only for flow analysis: the match loop below always assigns it, every rank's own name is in the class list */
  MPI_Get_processor_name(processor_name, &name_len);
  logmsg(2,"rank %d name len %d proc name %s\n",rank,name_len,processor_name);
  MPI_Allreduce(&name_len,&name_len_max,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
  name_len_max++;
  char *vec_pr_name=(char *) calloc (mpisize*name_len_max,sizeof(char));
  char *vec_pr_sname=(char *) calloc (mpisize*name_len_max,sizeof(char));
  name_beg=rank*name_len_max;
  for(i=name_beg; i<name_len+name_beg; i++)vec_pr_name[i]=processor_name[i-name_beg];
  vec_pr_name[i]='\0';
  MPI_Allreduce(vec_pr_name,vec_pr_sname,mpisize*name_len_max,MPI_CHAR,MPI_SUM,PETSC_COMM_WORLD);
  for(i=0; i<name_len_max; i++)vec_pr_name[i]=vec_pr_sname[i];
  j=1;
  for(i=0; i<mpisize; i++) {
    for(j1=0; j1<j; j1++) {
      if(strncmp(vec_pr_sname+i*name_len_max,vec_pr_name+j1*name_len_max,name_len_max)==0) {
        break;
      }
    }
    if(j1==j) {
      for(j6=0; j6<name_len_max; j6++)vec_pr_name[j*name_len_max+j6]=vec_pr_sname[j1*name_len_max+j6];
      j++;
    }
  }
  class_size=j;
  for(i=0; i<class_size; i++) {
    if(strncmp(processor_name,vec_pr_name+i*name_len_max,name_len_max)==0) {
      color=i;
      break;
    }
  }
  free(vec_pr_name);
  free(vec_pr_sname);
  MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&node_comm);
  MPI_Comm_rank( node_comm, &node_rank);
  MPI_Comm_size(node_comm,&group_size);
  if(node_rank==group_size-1)color=1;
  else color=0;
  MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&node_tail_comm);

  gettimeofday(&begintime, NULL);
  probe_begin=begintime;
  bool sbbd_overrid=false;
  PetscBool nohsl=false;
  if(rank==0) {
    /* ungated provenance line: every archived solver log records its
       producer (same constant as -version and stats.json) */
    printf("teems-solver %s\n",TEEMS_SOLVER_VERSION);
    fflush(stdout);
    logmsg(1,"Coefficient storage: %s precision\n",TEEMS_STORE_PRECISION);
    logmsg(1,"BLAS kernels: %s\n",blas_corename());
    logmsg(2,"Notes:\n  Shock statement values follow GEMPACK ordering (first subscript varies fastest).\n  Declare intertemporal variables with minimum dimension to minimise the net cut,\n  e.g. capital(REG,TIME)=qo(\"capital\",REG,TIME) rather than shocking qo(COM,REG,TIME).\n  laA/laDi control solver workspace sizes; use the smallest that solves.\n  Beware CRLF line endings in model text files.\n");
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  char run_id[64];
  solve_real *comp_approx_col=NULL; /* complementarity approximate-run solution, for <stem>.cols (kept across the accurate-run re-entry) */
  snprintf(run_id,sizeof(run_id),"%ld.%06ld-%ld",(long)begintime.tv_sec,(long)begintime.tv_usec,(long)getpid());
  MPI_Bcast(run_id,sizeof(run_id),MPI_CHAR,0,PETSC_COMM_WORLD);
  if(rank==0) {
    char cmfopt[TABREADLINE];
    PetscBool cmfflg=PETSC_FALSE;
    PetscOptionsGetString(NULL,NULL,"-cmdfile",cmfopt,TABREADLINE,&cmfflg);
    sidecars_clear_stale(cmfflg?cmfopt:"./reg.cmf");
  }
  if(cli_options_check(rank)>0) {
    PetscFinalize();
    return 1;
  }
  //**************************************************************************************
  //****************************** READ SET ELEMENT***************************************
  //**************************************************************************************
  char tabfile[TABREADLINE],newtabfile[TABREADLINE]="_temp_tab_file",newtabfile1[TABREADLINE]="_temp_tab_new_file",closure[TABREADLINE]="",shock[TABREADLINE]="",filename[TABREADLINE],longname[TABREADLINE],vname[NAMESIZE],copyline[TABREADLINE];
  char psfile[TABREADLINE]="\0";
  int npostsim=0,postsim_on=1;
  char tempchar[255],solmed[NAMESIZE],solchar[255];
  int niodata=0,nj,noutdata=0,nsoldata=0,nowrites=0,cofdump=1;
  offset_t nsetspace=0,ndblock=0,netcut=0,ndblock1,nreg=0,ntime=0;
  dim_t nset=0,vsize,dim1,nlength=0,matsol=0,laA=2,laDi=2,laD=2,nsbbdblocks=2,nesteddbbd=0,mc66=0,subints=1;
  offset_t alltimeset=-1,allregset=-1;
  map_def *maps=NULL;
  dim_t nmap=0;
  if(rank<10) {
    strcat(newtabfile,"000");
    strcat(newtabfile1,"000");
  }
  if(rank<100&&rank>9) {
    strcat(newtabfile,"00");
    strcat(newtabfile1,"00");
  }
  if(rank<1000&&rank>99) {
    strcat(newtabfile,"0");
    strcat(newtabfile1,"0");
  }
  sprintf(tempchar, "%d",rank);
  strcat(newtabfile,tempchar);
  strcat(newtabfile1,tempchar);
  strcat(newtabfile,".tab");
  strcat(newtabfile1,".tab");
  PetscOptionsGetInt(NULL,NULL,"-matsol",&matsol,NULL);/* enum matrix_method; >=MM_SBBD needs a regional or time set (first reg set orders variables) */
  if(matsol==MM_DBBD)nohsl=true;
  if(matsol==MM_NDBBD)nohsl=true;
  if(matsol==MM_NDBBD) {
    /* measured 2026-09-06 on the 4.36M-eq GTAP-RE rig, NDBBD-2 Gragg:
       the persistent refactorization is slower under NDBBD (19.0 vs
       17.4 s per step) and holds +0.4 GB per rank of resident factors;
       the flag is ignored here so a caller's LU/DBBD habit cannot
       cost the memory-bound method its reason to exist */
    dim_t frchk=0;
    PetscOptionsGetInt(NULL,NULL,"-fastrefac",&frchk,NULL);
    if(frchk) {
      if(rank==0)printf("Warning: -fastrefac is ignored under NDBBD (measured slower, +0.4 GB per rank); running without it\n");
      PetscOptionsSetValue(NULL,"-fastrefac","0");
    }
  }
  PetscOptionsGetInt(NULL,NULL,"-laA",&laA,NULL);
  if(laA==0)laA=2;
  PetscOptionsGetInt(NULL,NULL,"-laD",&laD,NULL);
  if(laD==0)laD=2;
  PetscOptionsGetInt(NULL,NULL,"-laDi",&laDi,NULL);
  if(laDi==0)laDi=2;
  PetscOptionsGetInt(NULL,NULL,"-withmc66",&mc66,NULL);
  /* Tier B3 (factor once, solve many): -residcheck 1 keeps the LHS
     matrix of one-shot SBBD/DBBD/NDBBD solves so the residual check
     (manual 30.1.5) covers every solve, at the cost of holding the
     matrix through the factorization; -fhtest 1 (kit hook) solves each
     step's right-hand side again with the kept factorization and
     compares. A method that cannot keep its factors stops here, before
     any work. */
  PetscOptionsGetInt(NULL,NULL,"-residcheck",&teems_resid_all,NULL);
  teems_resid_all=(teems_resid_all!=0);
  PetscOptionsGetInt(NULL,NULL,"-fhtest",&teems_fh_selftest,NULL);
  teems_fh_selftest=(teems_fh_selftest!=0);
  fh_request_check(matsol,mc66,rank);
  /* -refine {0,1} (default 1): one step of iterative refinement after
     every DBBD solve, with the kept factors and the kept LHS matrix
     (solve_drivers.c fh_refine_cols). The other methods are not refined:
     LU and SBBD are backward-stable already, and NDBBD cannot keep its
     factors yet. */
  {
    PetscBool rfset=PETSC_FALSE;
    PetscOptionsGetInt(NULL,NULL,"-refine",&teems_refine,&rfset);
    teems_refine=(teems_refine!=0);
    if(rfset&&teems_refine&&matsol!=MM_DBBD&&rank==0) {
      if(matsol==MM_NDBBD)logmsg(1,"Note: -refine 1 is ignored under NDBBD, which cannot keep its factorization for a refinement step yet; its solves are not refined\n");
      else logmsg(1,"Note: -refine applies to matrix_method DBBD only; this run's solves are not refined\n");
    }
  }
  /* -zdivshift x: shift every zero-divide default by x*(1+|default|),
     so a second run shows what depends on one (formula.c zdiv_shifted) */
  {
    PetscReal zs=0;
    PetscOptionsGetReal(NULL,NULL,"-zdivshift",&zs,NULL);
    teems_zdiv_shift=(double)zs;
  }
  PetscOptionsGetInt(NULL,NULL,"-step1",&steps1,NULL);
  if(steps1==0)steps1=2;
  PetscOptionsGetInt(NULL,NULL,"-step2",&steps2,NULL);
  if(steps2==0)steps2=4;
  PetscOptionsGetInt(NULL,NULL,"-step3",&steps3,NULL);
  if(steps3==0)steps3=8;
  PetscOptionsGetInt(NULL,NULL,"-single_run",&teems_single_run,NULL);
  teems_single_run=(teems_single_run!=0);
  PetscOptionsGetInt(NULL,NULL,"-two_run",&teems_two_run,NULL);
  teems_two_run=(teems_two_run!=0);
  /* -random_seed N: seed of RANDOM(a,b) (manual 11.5.2); the same seed
     reproduces every draw */
  {
    PetscInt rs=1;
    PetscBool rset=PETSC_FALSE;
    PetscOptionsGetInt(NULL,NULL,"-random_seed",&rs,&rset);
    if(rset&&rs<0) {
      errmsg("Error: -random_seed must be a non-negative integer, got %ld\n",(long)rs);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    teems_random_seed=(long)rs;
  }
  section_threads=0;
  max_threads=1;
  PetscOptionsGetInt(NULL,NULL,"-maxthreads",&max_threads,NULL);
  if(max_threads>1&&max_threads<=omp_get_max_threads( )){
    omp_set_num_threads(max_threads);
  }else{
    if(max_threads>1)printf("Warning: -maxthreads %d exceeds the OpenMP thread limit %d; using 1 thread\n",max_threads,omp_get_max_threads());
    max_threads=1;
    omp_set_num_threads(max_threads);
  }
  PetscOptionsGetInt(NULL,NULL,"-smllthreads",&section_threads,NULL);
  if(section_threads==0)section_threads=max_threads;
  PetscOptionsGetInt(NULL,NULL,"-nsubints",&subints,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nsbbdblocks",&nsbbdblocks,NULL);
  nesteddbbd=(matsol==MM_NDBBD)?1:0;/* nested ordering exists only for the NDBBD solve; any other pairing is broken */
  PetscOptionsGetInt(NULL,NULL,"-nowrites",&nowrites,NULL);
  PetscOptionsGetInt(NULL,NULL,"-cofdump",&cofdump,NULL); /* <stem>.cof/.cbin coefficient dump, default on */
  {
    PetscInt verb=verbosity;
    PetscOptionsGetInt(NULL,NULL,"-verbosity",&verb,NULL);/* 0 errors/results, 1 progress (default), 2 debug */
    verbosity=(int)verb;
    char verbstr[8];
    sprintf(verbstr,"%d",verbosity);
    setenv("TEEMS_VERBOSITY",verbstr,1);/* for the Fortran kernels (hsl_kernels.f90) */
  }
  {
    /* -ma48_cntl2 <x>: MA48/HSL_MP48 pivot threshold (CNTL(2), 0 < x <= 1).
       Absent = each library's own default (MA48 0.1, HSL_MP48 0.01),
       bit-identical to pre-option builds; given = applied at every
       MA48ID/MP48AD initialisation, the DBBD/NDBBD rank probes included,
       via TEEMS_MA48_CNTL2 (hsl_kernels.f90). Recorded in stats.json
       options as ma48_cntl2 (null when default). */
    PetscReal u=0;
    PetscBool uflg=PETSC_FALSE;
    PetscOptionsGetReal(NULL,NULL,"-ma48_cntl2",&u,&uflg);
    if(uflg) {
      if(!(u>0&&u<=1)) {
        if(rank==0)errmsg("Error: -ma48_cntl2 must be in (0,1] (MA48/HSL_MP48 pivot threshold CNTL(2); MA48 default 0.1, HSL_MP48 default 0.01), got %g\n",(double)u);
        PetscFinalize();
        return 1;
      }
      char ustr[32];
      snprintf(ustr,sizeof(ustr),"%.17g",(double)u);
      setenv("TEEMS_MA48_CNTL2",ustr,1);
      teems_ma48_cntl2_opt=(double)u;
    }
    else unsetenv("TEEMS_MA48_CNTL2");
  }
  {
    /* -ma48_cntl4 <x>: MA48 CNTL(4), pivots below x treated as zero, in
       the MA48/MA51 rank probes that size the DBBD/NDBBD blocks and the
       NDBBD interface cut (x >= 0; default 1e-4; 0 = MA48's own default,
       exact zeros only). Recorded in stats.json options as ma48_cntl4. */
    PetscReal c4=teems_ma48_cntl4;
    PetscOptionsGetReal(NULL,NULL,"-ma48_cntl4",&c4,NULL);
    if(!(c4>=0)) {
      if(rank==0)errmsg("Error: -ma48_cntl4 must be >= 0 (MA48 CNTL(4) rank-probe zero tolerance; default 1e-4), got %g\n",(double)c4);
      PetscFinalize();
      return 1;
    }
    teems_ma48_cntl4=(double)c4;
  }
  inmemory=-1;
  PetscOptionsGetInt(NULL,NULL,"-inmemory",&inmemory,NULL);/* keep value arrays resident instead of spilling to scratch */
  PetscOptionsGetInt(NULL,NULL,"-ndcutcache",&teems_ndcutcache,NULL);/* NDBBD: reuse the rank cuts across steps (1) or re-probe every step (0) */
  if(inmemory<0)inmemory=(matsol==MM_NDBBD)?0:1;/* default: resident except NDBBD, whose factor-file I/O wants the page cache that spilling frees (measured; LU/SBBD gain, DBBD neutral on both domains) */
  {
    /* Scratch directory for solver temp files: -tempdir option, else
       TMPDIR, else the compiled-in default (/tmp/). Capped at 200
       chars so appended file names always fit the 255/256-byte path
       buffers used throughout the solver. */
    char tmpopt[NAMESIZE];
    PetscBool tmpflg=PETSC_FALSE;
    PetscOptionsGetString(NULL,NULL,"-tempdir",tmpopt,NAMESIZE,&tmpflg);
    if(!tmpflg) {
      char *envtmp=getenv("TMPDIR");
      if(envtmp!=NULL&&envtmp[0]!='\0') {
        strncpy(tmpopt,envtmp,NAMESIZE-1);
        tmpopt[NAMESIZE-1]='\0';
        tmpflg=PETSC_TRUE;
      }
    }
    if(tmpflg) {
      size_t tlen=strlen(tmpopt);
      if(tlen==0||tlen>200) {
        if(rank==0)errmsg("Error: -tempdir must be 1-200 characters: %s\n",tmpopt);
        PetscFinalize();
        return 1;
      }
      if(tmpopt[tlen-1]!='/') {
        tmpopt[tlen]='/';
        tmpopt[tlen+1]='\0';
      }
      strcpy(scratch_dir,tmpopt);
    }
    /* Under -inmemory, and unless the user pinned a scratch location,
       place scratch on tmpfs so the block-factor handoff files written
       by the Fortran kernels never touch disk. */
    if(inmemory&&!tmpflg&&access("/dev/shm/",W_OK)==0) {
      strcpy(scratch_dir,"/dev/shm/");
      if(rank==0)logmsg(1,"inmemory: scratch on tmpfs (%s)\n",scratch_dir);
    }
    /* export for the Fortran kernels, which build their factor-file
       paths themselves (hsl_kernels.f90) */
    setenv("TEEMS_SCRATCH",scratch_dir,1);
    if(access(scratch_dir,W_OK)!=0) {
      if(rank==0)errmsg("Error: -tempdir is not a writable directory: %s\n",scratch_dir);
      PetscFinalize();
      return 1;
    }
  }
  /* Per-time-block interface ranks for NDBBD: the ordering presolve
     records min(nrow,ncol) per block, the rank-revealing factorization
     tightens it to the true numerical rank (allocated once ntime is
     known). */
  PetscInt *ndbbddrank1=NULL;
  logmsg(2,"matsol %d\n",matsol);
  PetscOptionsGetString(NULL,NULL,"-cmdfile",filename,TABREADLINE,&flg);
  if (!flg) {
    strcpy(filename,"./reg.cmf");//orani03.cmf");
  }
  PetscOptionsGetString(NULL,NULL,"-solmed",solmed,NAMESIZE,&flg);
  if (!flg) {
    strcpy(solmed,"Gragg");
  }
  if(strcmp(solmed,"Mmid")==0) {/* transitional alias: the multi-step method has always been Gragg's smoothed modified midpoint (Pearson 1991) */
    if(rank==0)printf("Warning: -solmed Mmid is deprecated and runs Gragg's method (smoothed modified midpoint) — use -solmed Gragg, or -solmed Midpoint for the midpoint method without the smoothing pass (GEMPACK manual 30.2)\n");
    strcpy(solmed,"Gragg");
  }
  if(strcmp(solmed,"NoSol")==0) {/* transitional alias: renamed — it is a structure probe, not a degenerate solve */
    if(rank==0)printf("Warning: -solmed NoSol is deprecated — use -solmed probe\n");
    strcpy(solmed,"probe");
  }
  int solmethod=0;
  if(strcmp(solmed,"Gragg")==0)solmethod=SM_GRAGG;
  if(strcmp(solmed,"Euler")==0)solmethod=SM_EULER;
  if(strcmp(solmed,"Midpoint")==0)solmethod=SM_MIDPOINT;
  if(strcmp(solmed,"RK2")==0)solmethod=SM_RK2;
  if(strcmp(solmed,"Heun")==0)solmethod=SM_HEUN;
  if(strcmp(solmed,"RK4")==0)solmethod=SM_RK4;
  if(strcmp(solmed,"BoSha32")==0)solmethod=SM_BOSHA32;
  if(strcmp(solmed,"DoPri54")==0)solmethod=SM_DOPRI54;
  if(strcmp(solmed,"Johansen")==0)solmethod=SM_JOHANSEN;
  if(strcmp(solmed,"probe")==0)solmethod=SM_PROBE;
  teems_probe_mode=(solmethod==SM_PROBE)?1:0;
  if(strcmp(solmed,"nosim")==0)solmethod=SM_NOSIM;
  if(solmethod==0) {
    if(rank==0)errmsg("Error: unknown -solmed %s (valid: Gragg, Midpoint, Euler, RK2, Heun, RK4, BoSha32, DoPri54, Johansen, probe, nosim)\n",solmed);
    PetscFinalize();
    return 1;
  }
  logmsg(2,"Sol med %d\n",solmethod);
  /* -solmed probe measures the system's structure irrespective of the
     requested -matsol: both structural detections run (chain dimension
     and diagonal-block partition, nested geometry when a chain exists),
     the method-vs-structure aborts below are skipped (nothing is
     solved), and the ordering rides the nohsl path like DBBD/NDBBD so
     the nested layout is valid. The stats.json this produces is what
     teems-R's matrix_method = "auto" consults (ROADMAP 6.10); the MC79
     diagnosis itself is ordering-invariant. */
  if(solmethod==SM_PROBE)nohsl=true;
  /* -single_run 1 (GEMPACK "method = euler|midpoint|gragg; steps = N;"): one pass
     of -step1 steps with no extrapolation. Johansen is a single step
     already; the Runge-Kutta family has its own step control. */
  if(teems_single_run) {
    if(solmethod==SM_RK2||solmethod==SM_HEUN||solmethod==SM_RK4||solmethod==SM_BOSHA32||solmethod==SM_DOPRI54) {
      errmsg("Error: -single_run 1 applies to -solmed Euler, Midpoint or Gragg; the Runge-Kutta methods (%s) are single-pass already and take their steps from -step1/-adaptive\n",solmed);
      PetscFinalize();
      return 1;
    }
    if(solmethod==SM_JOHANSEN||solmethod==SM_PROBE||solmethod==SM_NOSIM) teems_single_run=0;
    else if(rank==0) printf("Single-pass %s run: %d steps, no extrapolation (-single_run 1)\n",solmed,(int)steps1);
  }
  /* -two_run 1 (GEMPACK "steps = n1 n2;"): extrapolation from two
     multi-step solutions, -step1 and -step2, with no accuracy estimate
     (manual 26.1.2, 26.2.3) */
  if(teems_two_run) {
    if(teems_single_run) {
      errmsg("Error: -two_run 1 and -single_run 1 exclude each other\n");
      PetscFinalize();
      return 1;
    }
    if(solmethod!=SM_GRAGG&&solmethod!=SM_MIDPOINT&&solmethod!=SM_EULER) {
      if(solmethod==SM_JOHANSEN||solmethod==SM_PROBE||solmethod==SM_NOSIM) teems_two_run=0;
      else {
        errmsg("Error: -two_run 1 applies to -solmed Euler, Midpoint or Gragg (%s takes its steps from -step1/-adaptive)\n",solmed);
        PetscFinalize();
        return 1;
      }
    }
    if(teems_two_run) {
      steps3=steps2;
      if(rank==0) printf("Two-solution %s run: %d and %d steps extrapolated, no accuracy estimate (-two_run 1)\n",solmed,(int)steps1,(int)steps2);
    }
  }
  /* -sup 1|2 (manual 26.8.1 SUP = last|all): the updated data after
     the last or every separate multi-step solution, <stem>.ud5..ud7 in
     the .cbin0 layout (phase 5..7), indexed by the run's .cof. GEMPACK
     allows it with neither subintervals nor automatic accuracy. */
  PetscOptionsGetInt(NULL,NULL,"-sup",&teems_sup,NULL);
  PetscOptionsGetInt(NULL,NULL,"-sui",&teems_sui,NULL);
  teems_sui=(teems_sui!=0);
  if(teems_sui&&!cofdump) {
    errmsg("Error: -sui marks Formula (Initial) coefficients in the coefficient dump, so it needs -cofdump 1\n");
    PetscFinalize();
    return 1;
  }
  if(teems_sup<0||teems_sup>2) {
    errmsg("Error: -sup takes 0 (none), 1 (last) or 2 (all)\n");
    PetscFinalize();
    return 1;
  }
  if(teems_sup&&solmethod!=SM_GRAGG&&solmethod!=SM_MIDPOINT&&solmethod!=SM_EULER) {
    errmsg("Error: -sup saves the updated data of separate multi-step solutions, so it applies to -solmed Euler, Midpoint or Gragg (GEMPACK manual 26.8.1)\n");
    PetscFinalize();
    return 1;
  }
  if(teems_sup&&subints>1) {
    errmsg("Error: -sup cannot be used with more than one subinterval (GEMPACK manual 26.8.1)\n");
    PetscFinalize();
    return 1;
  }
  if(teems_sup&&!cofdump) {
    errmsg("Error: -sup needs the coefficient dump (-cofdump 1): the .cof declarations index the updated-data files\n");
    PetscFinalize();
    return 1;
  }
  /* -probefine: with -solmed probe, add the MC79 fine Dulmage-Mendelsohn
     report (strongly connected components of the well-determined block) */
  dim_t probefine=0;
  PetscOptionsGetInt(NULL,NULL,"-probefine",&probefine,NULL);
  /* -jacdump 1 (Tier B5): the base-point Jacobian of the condensed
     system over every variable element, written before the first step
     (jacobian.c jacobian_dump); the run then goes on as asked. -jacdump
     2 stops after the export, as a probe without its diagnosis: the
     homogeneity check needs the matrix, not a solution */
  dim_t jacdump=0;
  bool jac_only=false;
  PetscOptionsGetInt(NULL,NULL,"-jacdump",&jacdump,NULL);
  if(jacdump<0||jacdump>2) {
    if(rank==0)errmsg("Error: -jacdump must be 0 (off), 1 (write <solfiles>.jac) or 2 (write it and stop before the solve), got %ld\n",(long)jacdump);
    PetscFinalize();
    return 1;
  }
  /* -condest: opt-in per-solve quality diagnostics (MA60 iterative
     refinement + Arioli-Demmel-Duff backward/forward error and scaled
     condition numbers) on the sequential LU path.  Diagnostic-only:
     solutions are untouched.  Not implemented for the bordered
     methods (condition estimation needs transpose solves on the
     composed system, which only the LU path has). */
  {
    dim_t condest=0;
    PetscOptionsGetInt(NULL,NULL,"-condest",&condest,NULL);
    teems_condest=(condest!=0);
    if(teems_condest&&matsol!=MM_LU) {
      if(rank==0)printf("Warning: -condest is implemented for the sequential LU path (-matsol 0) only; ignored for this run\n");
      teems_condest=0;
    }
  }
  bool isrk=(solmethod==SM_RK2||solmethod==SM_HEUN||solmethod==SM_RK4||solmethod==SM_BOSHA32||solmethod==SM_DOPRI54);
  bool isrk_embedded=(solmethod==SM_BOSHA32||solmethod==SM_DOPRI54);
  /* Runge-Kutta controls (GEMPACK 26.5.1/26.5.2): -adaptive
     no|yes|accuracy-only (embedded flavors only; accuracy-only skips
     the check-failure retries), -epstol the per-step error-metric
     bound, -retryadj/-maxretries the check-failure retry policy */
  int adaptive=0,maxretries=3;
  PetscReal epstol=0.1,retryadj=0.5;
  rk_options rko={RK_CHART_LOG,RK_NORM_MAX,RK_CTRL_STD,RK_SCOPE_PCT,0.0,log(1e30)};
  {
    char adaptbuf[NAMESIZE];
    PetscOptionsGetString(NULL,NULL,"-adaptive",adaptbuf,NAMESIZE,&flg);
    if(flg) {
      if(strcmp(adaptbuf,"no")==0)adaptive=0;
      else if(strcmp(adaptbuf,"yes")==0)adaptive=1;
      else if(strcmp(adaptbuf,"accuracy-only")==0)adaptive=2;
      else {
        if(rank==0)errmsg("Error: unknown -adaptive %s (valid: no, yes, accuracy-only)\n",adaptbuf);
        PetscFinalize();
        return 1;
      }
      if(adaptive&&!isrk_embedded) {
        if(rank==0)errmsg("Error: -adaptive requires an embedded Runge-Kutta method (-solmed BoSha32 or DoPri54)\n");
        PetscFinalize();
        return 1;
      }
    }
    PetscOptionsGetReal(NULL,NULL,"-epstol",&epstol,NULL);
    PetscOptionsGetReal(NULL,NULL,"-retryadj",&retryadj,NULL);
    PetscOptionsGetInt(NULL,NULL,"-maxretries",&maxretries,NULL);
    /* -rkchart log|percent: the state chart the stages combine in (log
       levels: positivity unconditional, Munthe-Kaas 1999; percent: the
       GEMPACK-orientation arithmetic); -rknorm max|rms: the per-step
       error-metric norm the accept test uses; -rkctrl std|pi: the
       step-size controller; -rk_h0: initial step length (default 1/step1
       capped by the initial gradient); -rkguard: log-chart level ratio
       beyond which a stage is rejected */
    {
      char optbuf[NAMESIZE];
      PetscReal rtmp=0;
      PetscOptionsGetString(NULL,NULL,"-rkchart",optbuf,NAMESIZE,&flg);
      if(flg) {
        if(strcmp(optbuf,"log")==0)rko.chart=RK_CHART_LOG;
        else if(strcmp(optbuf,"percent")==0)rko.chart=RK_CHART_PERCENT;
        else {
          if(rank==0)errmsg("Error: unknown -rkchart %s (valid: log, percent)\n",optbuf);
          PetscFinalize();
          return 1;
        }
      }
      PetscOptionsGetString(NULL,NULL,"-rknorm",optbuf,NAMESIZE,&flg);
      if(flg) {
        if(strcmp(optbuf,"max")==0)rko.norm=RK_NORM_MAX;
        else if(strcmp(optbuf,"rms")==0)rko.norm=RK_NORM_RMS;
        else {
          if(rank==0)errmsg("Error: unknown -rknorm %s (valid: max, rms)\n",optbuf);
          PetscFinalize();
          return 1;
        }
      }
      PetscOptionsGetString(NULL,NULL,"-rkctrl",optbuf,NAMESIZE,&flg);
      if(flg) {
        if(strcmp(optbuf,"std")==0)rko.ctrl=RK_CTRL_STD;
        else if(strcmp(optbuf,"pi")==0)rko.ctrl=RK_CTRL_PI;
        else {
          if(rank==0)errmsg("Error: unknown -rkctrl %s (valid: std, pi)\n",optbuf);
          PetscFinalize();
          return 1;
        }
      }
      PetscOptionsGetString(NULL,NULL,"-rkscope",optbuf,NAMESIZE,&flg);
      if(flg) {
        if(strcmp(optbuf,"pct")==0)rko.scope=RK_SCOPE_PCT;
        else if(strcmp(optbuf,"all")==0)rko.scope=RK_SCOPE_ALL;
        else {
          if(rank==0)errmsg("Error: unknown -rkscope %s (valid: pct, all)\n",optbuf);
          PetscFinalize();
          return 1;
        }
      }
      PetscOptionsGetReal(NULL,NULL,"-rk_h0",&rtmp,&flg);
      if(flg) {
        if(rtmp<=0||rtmp>1) {
          if(rank==0)errmsg("Error: -rk_h0 must lie in (0, 1] (got %g)\n",(double)rtmp);
          PetscFinalize();
          return 1;
        }
        rko.h0=(double)rtmp;
      }
      PetscOptionsGetReal(NULL,NULL,"-rkguard",&rtmp,&flg);
      if(flg) {
        if(rtmp<=1) {
          if(rank==0)errmsg("Error: -rkguard must exceed 1 (a level ratio; got %g)\n",(double)rtmp);
          PetscFinalize();
          return 1;
        }
        rko.guard=log((double)rtmp);
      }
    }
    if(isrk) {
      if(steps1<1) {
        if(rank==0)errmsg("Error: -step1 must be at least 1 for Runge-Kutta methods (got %d)\n",steps1);
        PetscFinalize();
        return 1;
      }
      if(epstol<=0||retryadj<=0||retryadj>=1||maxretries<1) {
        if(rank==0)errmsg("Error: -epstol and -retryadj must be positive (-retryadj below 1) and -maxretries at least 1\n");
        PetscFinalize();
        return 1;
      }
      if(adaptive&&epstol<0.005) {
        if(rank==0)printf("Warning: -epstol %g is below 0.005; tolerances this tight are hard to achieve numerically and may reject many steps\n",(double)epstol);
      }
    }
  }
  if((solmethod==SM_GRAGG||solmethod==SM_MIDPOINT||solmethod==SM_EULER)&&steps1<1) {
    if(rank==0)errmsg("Error: -step1 must be at least 1 (got %d)\n",steps1);
    PetscFinalize();
    return 1;
  }
  if(subints<1) {
    if(rank==0)errmsg("Error: -nsubints must be at least 1 (got %d)\n",subints);
    PetscFinalize();
    return 1;
  }
  if((solmethod==SM_GRAGG||solmethod==SM_MIDPOINT)&&!teems_single_run) {
    /* the h^2 error expansion of Gragg and the midpoint method holds
       within one step parity (Pearson 1991 Thm 6.1; manual 26.1.2,
       30.2); mixed parity also breaks the shared-power extrapolation */
    step_ratio2=steps1/(double)2;
    i=(offset_t)steps1/2;
    if(step_ratio2!=i){
      //odd
      step_ratio2=steps2/(double)2;
      i=(offset_t)steps2/2;
      if(step_ratio2==i){
        errmsg("Error: -step1/-step2/-step3 must be all odd or all even (got %d %d %d)\n",steps1,steps2,steps3);
        PetscFinalize();
        return 1;
      }
      step_ratio2=steps3/(double)2;
      i=(offset_t)steps3/2;
      if(step_ratio2==i){
        errmsg("Error: -step1/-step2/-step3 must be all odd or all even (got %d %d %d)\n",steps1,steps2,steps3);
        PetscFinalize();
        return 1;
      }
    }else{
      //even
      step_ratio2=steps2/(double)2;
      i=(offset_t)steps2/2;
      if(step_ratio2!=i){
        errmsg("Error: -step1/-step2/-step3 must be all odd or all even (got %d %d %d)\n",steps1,steps2,steps3);
        PetscFinalize();
        return 1;
      }
      step_ratio2=steps3/(double)2;
      i=(offset_t)steps3/2;
      if(step_ratio2!=i){
        errmsg("Error: -step1/-step2/-step3 must be all odd or all even (got %d %d %d)\n",steps1,steps2,steps3);
        PetscFinalize();
        return 1;
      }
    }
  }
  if((solmethod==SM_GRAGG||solmethod==SM_MIDPOINT||solmethod==SM_EULER)&&teems_two_run&&!(steps1<steps2)) {
    errmsg("Error: -step1/-step2 must be strictly increasing (got %d %d)\n",steps1,steps2);
    PetscFinalize();
    return 1;
  }
  if((solmethod==SM_GRAGG||solmethod==SM_MIDPOINT||solmethod==SM_EULER)&&!teems_single_run&&!teems_two_run&&!(steps1<steps2&&steps2<steps3)) {
    /* extrapolation needs three distinct step sizes (GEMPACK: i<j<k) */
    errmsg("Error: -step1/-step2/-step3 must be strictly increasing (got %d %d %d)\n",steps1,steps2,steps3);
    PetscFinalize();
    return 1;
  }
  step_ratio2=steps2/(double)steps1;
  steps2=(PetscInt)steps2/steps1;
  step_ratio3=steps3/(double)steps1;
  steps3=(PetscInt)steps3/steps1;

  #pragma omp parallel private(i)
  {
  i=0;
  }

  if(nohsl) {
    rank_hsl=rank;
  }
  else {
    rank_hsl=0;
  }
  char *readitem=NULL;
  char subfile[TABREADLINE]="";
  int has_sub=0;
  if(rank==0) {
    if(manifest_check(filename)<0)MPI_Abort(PETSC_COMM_WORLD,1);
    has_sub=cmf_subtotals_file(filename,subfile);
    niodata=cmf_count_files(filename,"iodata");
    if(niodata==-1)MPI_Abort(PETSC_COMM_WORLD,1);/* rank 0 only: the other ranks wait in the Bcast below */
    noutdata=cmf_count_files(filename,"outdata");
    nsoldata=cmf_count_files(filename,"soldata");
  }
  /* Tier B3 subtotals: the combinations the step accumulation does not
     cover stop here, before any work */
  MPI_Bcast(&has_sub,1,MPI_INT,0,PETSC_COMM_WORLD);
  if(has_sub) {
    const char *why=NULL;
    if(matsol==MM_NDBBD)why="matrix_method NDBBD, which cannot keep its factorization for more solves yet";
    else if(matsol==MM_SBBD&&mc66!=0)why="SBBD under -withmc66 1, which cannot keep its factorization for more solves";
    else if(solmethod==SM_NOSIM)why="-solmed nosim, which runs no simulation";
    else if(solmethod!=SM_JOHANSEN&&solmethod!=SM_GRAGG&&solmethod!=SM_MIDPOINT&&solmethod!=SM_EULER&&solmethod!=SM_PROBE)why="a Runge-Kutta method, whose stage combination in the log chart has no settled subtotal convention yet";
    if(why!=NULL) {
      if(rank==0)errmsg("Error: subtotals (GEMPACK manual 29) are not available with %s; use matrix_method LU, SBBD or DBBD with the Johansen, Euler, midpoint or Gragg method\n",why);
      MPI_Barrier(PETSC_COMM_WORLD);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
  if(nohsl) {
    MPI_Bcast(&niodata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&noutdata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nsoldata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  cmf_file_entry *iodata= (cmf_file_entry *) calloc (niodata+noutdata+nsoldata,sizeof(cmf_file_entry));
  /* run-mode switches travel on the invocation like every other run
     control (the CMF is a file manifest by design); parsed on every
     rank from the shared argv, recorded in stats.json's options
     object. -assertions / -range_test_initial / -range_test_updated:
     0 = off, 1 = warn, 2 = fatal. -postsim: 0|1. */
  {
    dim_t mopt;
    mopt=teems_assertions_mode;
    PetscOptionsGetInt(NULL,NULL,"-assertions",&mopt,NULL);
    if(mopt<0||mopt>2) {
      if(rank==0)errmsg("Error: -assertions must be 0 (off), 1 (warn) or 2 (fatal)\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    teems_assertions_mode=(int)mopt;
    mopt=teems_range_test_initial;
    PetscOptionsGetInt(NULL,NULL,"-range_test_initial",&mopt,NULL);
    if(mopt<0||mopt>2) {
      if(rank==0)errmsg("Error: -range_test_initial must be 0 (off), 1 (warn) or 2 (fatal)\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    teems_range_test_initial=(int)mopt;
    mopt=teems_range_test_updated;
    PetscOptionsGetInt(NULL,NULL,"-range_test_updated",&mopt,NULL);
    if(mopt<0||mopt>2) {
      if(rank==0)errmsg("Error: -range_test_updated must be 0 (off), 1 (warn) or 2 (fatal)\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    teems_range_test_updated=(int)mopt;
    mopt=postsim_on;
    PetscOptionsGetInt(NULL,NULL,"-postsim",&mopt,NULL);
    postsim_on=mopt?1:0;
  }
  if(rank==rank_hsl) {
    cmf_read(filename,niodata,iodata,tabfile,closure,shock);
    for (nj=0; nj<niodata+noutdata+nsoldata; nj++) logmsg(2,"rank %d logname %s fname %s\n",rank,iodata[nj].logname,iodata[nj].filname);
    if(tab_preprocess(tabfile,newtabfile)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    /* audit A6: fail fast on unsupported/unknown Default statements
       before any reader applies them positionally */
    if(tab_defaults_validate(newtabfile)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    /* Tier 0: route POSTSIM sections to the _ps companion, consumed
       once after the solve */
    strcpy(psfile,newtabfile);
    strcat(psfile,"_ps");
    npostsim=tab_postsim_split(newtabfile,psfile);
    /* fail-fast: exit 0 here used to mask section errors */
    if(npostsim<0)MPI_Abort(PETSC_COMM_WORLD,1);
    /* conditional set builders (manual 10.1.3): data-dependent
       conditions evaluated from the input files, statements rewritten
       into explicit lists + subset relations BEFORE set resolution */
    if(tab_setbuilder_transform(newtabfile,iodata,niodata)<0)MPI_Abort(PETSC_COMM_WORLD,1);
    /* C1: parse Complementarity statements, validate (10.17/11.14)
       and replace them with their derived levels statements (51.7.2)
       BEFORE the levels transform expands those; untouched when the
       TAB has none (design doc section 7) */
    if(tab_complementarity_transform(newtabfile)<0)MPI_Abort(PETSC_COMM_WORLD,1);
    /* C0: expand Formula&Equation, pair levels variables with value
       coefficients + updates, linearize Equation (levels) by change
       differentiation; untouched when the TAB has no levels
       statements (design doc section 5) */
    if(tab_levels_transform(newtabfile)<0)MPI_Abort(PETSC_COMM_WORLD,1);
    /* loops (manual 11.18): marker statements, body statements over
       the loop, BREAK/CYCLE as probe assertions; a PostSim section's
       loops are in the companion file */
    if(tab_loop_transform(newtabfile)<0)MPI_Abort(PETSC_COMM_WORLD,1);
    if(npostsim>0&&tab_loop_transform(psfile)<0)MPI_Abort(PETSC_COMM_WORLD,1);
  }

  strcpy(tabfile,newtabfile);
  if(rank==0)nset=sets_count(tabfile)+teems_loop_nsyn; /* + one lp@<id> set per general loop */
  /* Tier B4: a TAB without equations is a data program (manual 5.1.2,
     6.3; GEMPACK runs it without a simulation), and -solmed nosim runs
     any TAB that way (CMF "simulation = no;", 25.1.8): reads, formulas
     and assertions, then the writes and dumps, with no closure, shocks
     or solve */
  {
    long neqstmt=0;
    if(rank==0)neqstmt=(long)tab_count_statements(tabfile,"equation");
    MPI_Bcast(&neqstmt,1,MPI_LONG,0,PETSC_COMM_WORLD);
    if(neqstmt==0&&solmethod!=SM_NOSIM) {
      if(rank==0)logmsg(1,"The TAB has no equations: running it as a data program (reads, formulas, assertions and writes; no closure, shocks or solve; GEMPACK manual 5.1.2)\n");
      solmethod=SM_NOSIM;
      strcpy(solmed,"nosim");
    }
    if(solmethod==SM_NOSIM) {
      const char *why=NULL;
      if(has_sub)why="subtotals";
      else if(jacdump)why="-jacdump";
      if(why!=NULL) {
        if(rank==0)errmsg("Error: %s needs a simulation, and this run has none (%s); drop %s or give the TAB its equations and a closure\n",why,neqstmt==0?"the TAB has no equations":"-solmed nosim",why);
        MPI_Barrier(PETSC_COMM_WORLD);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      if(rank==0&&neqstmt>0)logmsg(1,"No-simulation run (-solmed nosim): %ld equation statement(s) are not solved; reads, formulas, assertions and writes only (GEMPACK manual 25.1.8)\n",neqstmt);
    }
  }
  if(nohsl) {
    MPI_Bcast(iodata,niodata*sizeof(cmf_file_entry), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(closure,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(shock,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nset,sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  set_def *sets= (set_def *) calloc (nset,sizeof(set_def));
  teems_set_isprod= (bool *) calloc (nset>0?nset:1,sizeof(bool));
  teems_set_prod1= (dim_t *) calloc (nset>0?nset:1,sizeof(dim_t));
  teems_set_prod2= (dim_t *) calloc (nset>0?nset:1,sizeof(dim_t));
  teems_set_itstem= (char (*)[NAMESIZE]) calloc (nset>0?nset:1,NAMESIZE);
  teems_set_itfirst= (int *) calloc (nset>0?nset:1,sizeof(int));
  for(i=0; i<nset; i++) {
    sets[i].subsetid[0]=i;
    for(j=0; j<MAXSUPSET; j++)sets[i].subsetid[j]=-1;
  }
  if(rank==0) {
    if(sets_read(tabfile,niodata,iodata, sets,nset)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    if(sets_read_intertemporal(tabfile,niodata,iodata, sets,nset)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    for (i=0; i<teems_loop_nsyn; i++) {
      set_def *ls=&sets[nset-teems_loop_nsyn+i];
      dim_t par;
      snprintf(ls->setname,sizeof(ls->setname),"lp@%d",teems_loop_syns[i].id);
      snprintf(ls->readele,sizeof(ls->readele),"~%s",teems_loop_syns[i].parent);
      for (par=0; par<nset-teems_loop_nsyn; par++) if (strcmp(sets[par].setname,teems_loop_syns[i].parent)==0) break;
      if (par==nset-teems_loop_nsyn) {
        errmsg("Error: Loop over %s, which is not a declared set (GEMPACK manual 11.18)\n",teems_loop_syns[i].parent);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      ls->size=(sets[par].size>0)?1:0;
    }
    for (i=0; i<nset; i++) {
      if(sets[i].size<0){
        errmsg("Error: set %s has a negative size in TAB file\n",sets[i].setname);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      if(sets[i].size>MAXSETSIZE){
        errmsg("Error: set %s declares %ld elements, above the %d limit; check the element count in its data-file header\n",sets[i].setname,(long)sets[i].size,MAXSETSIZE);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      sets[i].offset=nsetspace;
      nsetspace=nsetspace+sets[i].size;
    }
  }
  if(nohsl) {
    MPI_Bcast(sets,nset*sizeof(set_def), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nsetspace,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(teems_set_itstem,(nset>0?nset:1)*NAMESIZE, MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(teems_set_itfirst,(nset>0?nset:1)*sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  set_element *set_elems= (set_element *) calloc (nsetspace,sizeof(set_element));
  /* set sizes come from the data-file headers (sets_read), so a header
     declaring an absurd count makes nsetspace absurd and this calloc
     fail -- the loop below then wrote through NULL (fuzz batch 13
     data-file probe: a count of 2147483647 segfaulted). Only size<0 was
     checked above. */
  if (nsetspace>0&&set_elems==NULL) {
    errmsg("Error: cannot allocate %ld set elements; check the element counts the data-file headers declare\n",(long)nsetspace);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  teems_sets=sets;
  teems_nset=nset;
  teems_set_elems=set_elems;
  logmsg(2,"nset %d nsetspace %ld\n",nset,nsetspace);
  for (i=0; i<nsetspace; i++)for (j=0; j<MAXSUPSET; j++)set_elems[i].superset_pos[j]=-1;
  if(rank==0) {
    for (i=0; i<nset; i++) {
      if (sets[i].readele[0]=='~') continue; /* loop set: tab_loop_sets_link */
      if (sets[i].readele[0]=='#') { /* unnamed set (manual 11.7.1): elements 1..n */
        for (j=0; j<sets[i].size; j++) {
          snprintf(set_elems[j+sets[i].offset].setele,NAMESIZE,"%ld",(long)(j+1));
          set_elems[j+sets[i].offset].superset_pos[0]=j;
        }
        continue;
      }
      strcpy(vname,sets[i].header);
      nlength=0;
      while (vname[nlength] != '\0') {
        nlength++;
      }
      if (nlength>0) {
        datafile_read_header_info(vname,iodata[sets[i].fileid].filname,&vsize,longname,&dim1);
        /* dim1 is the element count the DATA FILE declares; set_elems has
           only sets[i].size slots at this set's offset (allocated above
           from the TAB's declared sizes), so a larger count wrote past
           the array -- a header count of 9999 against a 3-element set
           segfaulted (fuzz batch 13 data-file probe). The TAB's size is
           authoritative, as the element-list path below assumes. */
        if (dim1<0||dim1>sets[i].size) {
          errmsg("Error: header \"%s\" in %s supplies %ld elements but set %s is declared with %ld in the TAB file\n",
                 vname,iodata[sets[i].fileid].filname,(long)dim1,sets[i].setname,(long)sets[i].size);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
        datafile_labels *matvar1= (datafile_labels *) calloc (dim1,sizeof(datafile_labels));
        if (dim1>0&&matvar1==NULL) {
          errmsg("Error: out of memory reading the elements of set %s\n",sets[i].setname);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
        datafile_read_labels(vname,iodata[sets[i].fileid].filname,dim1,matvar1);
        for (j=0; j<dim1; j++) {
          nj=0;
          while(matvar1[j].ch[nj]!='\0') {
            matvar1[j].ch[nj]=tolower((int) matvar1[j].ch[nj]);
            nj++;
          }
          strcpy(set_elems[j+sets[i].offset].setele,matvar1[j].ch);
          set_elems[j+sets[i].offset].superset_pos[0]=j;
        }
        free(matvar1);
      }
      else {
        if (sets[i].readele[0]=='@') {
          set_expr_build(set_elems, sets,nset,i);
        }
        else if (sets[i].readele[0]=='-'&&sets[i].readele[1]==',') {
          set_difference(set_elems, sets,nset,i);
        }
        else {
          if (sets[i].readele[0]=='+'&&sets[i].readele[1]==',') {
            set_union_op(set_elems, sets,nset,i);
          }
          else {
            if (sets[i].readele[0]=='^'&&sets[i].readele[1]==',') {
              set_union_named(set_elems, sets,nset,i);
            }
            else {
              if(sets[i].readele[0]=='=') {
                set_equality_build(set_elems, sets, i);
              }
              else {
                dim1=sets[i].size;
                if (dim1<=0) continue; /* nothing allocated to populate; last-set offset==nsetspace */
                strcpy(copyline,sets[i].readele);
                strcat(copyline,",");
                str_delete_char(copyline,' ');
                readitem = strtok(copyline,",");
                if (readitem==NULL||strlen(readitem)>=NAMESIZE) {
                  errmsg("Error: malformed element list for set %s\n",sets[i].setname);
                  MPI_Abort(PETSC_COMM_WORLD,1);
                }
                strcpy(set_elems[sets[i].offset].setele,readitem);
                set_elems[sets[i].offset].superset_pos[0]=0;
                for (j=1; j<dim1; j++) {
                  readitem = strtok(NULL,",");
                  if (readitem==NULL||strlen(readitem)>=NAMESIZE) {
                    errmsg("Error: malformed element list for set %s\n",sets[i].setname);
                    MPI_Abort(PETSC_COMM_WORLD,1);
                  }
                  strcpy(set_elems[j+sets[i].offset].setele,readitem);
                  set_elems[j+sets[i].offset].superset_pos[0]=j;
                }
              }
            }
          }
        }
      }

    }
    subsets_read(tabfile, set_elems, sets,nset);
    j2=1;
    while(j2==1)for(i=1; i<MAXSUPSET; i++)subset_map_build(set_elems,sets,nset,&j2); //printf("check %d\n",i);}
    if(tab_loop_sets_link(set_elems,sets,nset)<0)MPI_Abort(PETSC_COMM_WORLD,1);
    ndblock1=ndblock;
    /* MAPPING statements (manual 11.9): declarations + by_elements
       values, resolved against the set elements built above so they
       exist before any coefficient/equation machinery runs */
    nmap=tab_count_statements(tabfile,"mapping");
    if(nmap>0) {
      maps= (map_def *) calloc (nmap+MAP_SYNTH_MAX,sizeof(map_def));
      mappings_read(tabfile,maps,nmap,sets,nset);
      mapping_values_read(tabfile,niodata,iodata,maps,nmap,sets,nset,set_elems);
      /* formula-assigned mappings (manual 10.13.1) get their values
         when their Formula executes (main passes or PostSim), so the
         completeness check moves to first use for those */
      mapping_formula_scan(tabfile,maps,nmap);
      if(npostsim>0)mapping_formula_scan(psfile,maps,nmap);
      mappings_validate(maps,nmap,sets,set_elems);
    }
    /* C1: complementarity quantifier sets vs the variable's and
       bounds' argument sets (11.14 points 2-3) -- needs the set
       elements built above for the ordered-subset walk */
    if(complementarities_validate(sets,nset,set_elems)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  if(nohsl) {
    MPI_Bcast(sets,nset*sizeof(set_def), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(set_elems,nsetspace*sizeof(set_element), MPI_BYTE,0, PETSC_COMM_WORLD);
    /* every rank that runs formulas repoints its own loop sets */
    MPI_Bcast(&teems_loop_nsyn,1,MPI_INT,0,PETSC_COMM_WORLD);
    if(teems_loop_nsyn>0) {
      if(rank!=0) {
        free(teems_loop_syns); /* nohsl ranks ran the transform too */
        teems_loop_syns=(teems_loop_syn *) calloc (teems_loop_nsyn,sizeof(teems_loop_syn));
      }
      MPI_Bcast(teems_loop_syns,teems_loop_nsyn*sizeof(teems_loop_syn),MPI_BYTE,0,PETSC_COMM_WORLD);
    }
    MPI_Bcast(&nmap,sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    if(nmap>0) {
      if(rank!=0) maps= (map_def *) calloc (nmap+MAP_SYNTH_MAX,sizeof(map_def));
      MPI_Bcast(maps,nmap*sizeof(map_def), MPI_BYTE,0, PETSC_COMM_WORLD);
      for(i=0; i<nmap; i++) {
        if(rank!=0) maps[i].values= (dim_t *) calloc (sets[maps[i].fromset].size>0?sets[maps[i].fromset].size:1,sizeof(dim_t));
        if(rank!=0) maps[i].assigned= (unsigned char *) calloc (sets[maps[i].fromset].size>0?sets[maps[i].fromset].size:1,sizeof(unsigned char));
        MPI_Bcast(maps[i].values,sets[maps[i].fromset].size*sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
    }
    /* C2: comp records leave rank 0 (fixed-size PODs; the state
       machinery and the closure/driver dispatch read them) */
    MPI_Bcast(&teems_ncomp,sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    if(teems_ncomp>0) {
      if(rank!=0) {
        /* under nohsl every rank ran the transform and already holds
           its own records; rank 0's are canonical -- drop ours before
           the replacement or the transform's array leaks */
        free(teems_comps);
        teems_comps= (comp_def *) calloc (teems_ncomp,sizeof(comp_def));
      }
      MPI_Bcast(teems_comps,teems_ncomp*sizeof(comp_def), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
  }
  teems_maps=maps;
  teems_nmap=nmap;
  teems_nmap_user=nmap;
  /* The chain dimension and the diagonal-block partition are derived
     structurally just before the ordering, once the equations are
     readable. Only the bordered methods consume these dimensions. */
  bool structural_time=(matsol==MM_SBBD||matsol==MM_DBBD||matsol==MM_NDBBD||solmethod==SM_PROBE);
  bool structural_reg=(matsol==MM_DBBD||matsol==MM_NDBBD||solmethod==SM_PROBE);
  ndblock=ndblock1;

  //**************************************************************************************
  //****************************** END READ SET ELEMENT***********************************
  //**************************************************************************************

  //**************************************************************************************
  //****************************** READ COEFFICIENT NAME**********************************
  //**************************************************************************************
  char commsyntax[NAMESIZE];
  strcpy(commsyntax,"coefficient");
  offset_t ncof=0,ncofele=0,ncof1,ncofele1;
  if(rank==0) {
    logmsg(2,"tabfile %s\n",tabfile);
    ncof=tab_count_statements(tabfile,commsyntax);
    logmsg(2,"tabfile %s ncof %ld\n",tabfile,ncof);
    ncof1=ncof;
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  logmsg(2,"rank %d ncof %ld\n",rank,ncof);
  if(nohsl)MPI_Bcast(&ncof1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  ncof=ncof1;
  array_def *coefs= (array_def *) calloc (ncof,sizeof(array_def));//recycle ha_cgeset
  if(rank==0) {
    ncofele=coefficients_read(tabfile,commsyntax,coefs,ncof,sets,nset);
    /* fail-fast: exit 0 here used to mask declaration errors */
    if(ncofele==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    /* PostSim residuals: mark split-recorded PostSim coefficients for
       the 12.2.2 LHS and 12.2.3 Read-target rules */
    postsim_mark_coefs(coefs,ncof);
    ncofele1=ncofele;
  }
  if(nohsl) {
    MPI_Bcast(&ncofele1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(coefs,ncof*sizeof(array_def), MPI_BYTE,0, PETSC_COMM_WORLD);
    /* every rank runs the step formulas under nohsl: the integer flags
       decide their INITIAL default (manual 11.6.3) */
    MPI_Bcast(&teems_n_int_coefs,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    if(rank!=0) {
      free(teems_coef_is_int);
      teems_coef_is_int= (bool *) calloc (ncof+1,sizeof(bool));
    }
    MPI_Bcast(teems_coef_is_int,(ncof+1)*sizeof(bool), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  ncofele=ncofele1;
  logmsg(2,"rank %d ncofele %ld\n",rank,ncofele);
  //**************************************************************************************
  //****************************** END READ COEFFICIENT NAME******************************
  //**************************************************************************************

  //**************************************************************************************
  //****************************** READ VARIABLE NAME*************************************
  //**************************************************************************************
  strcpy(commsyntax,"variable");
  offset_t nvar=0,nvarele=0,nvar1,nvarele1;
  if(rank==0) {
    nvar=tab_count_statements(tabfile,commsyntax);
    nvar1=nvar;
  }
  if(nohsl)MPI_Bcast(&nvar1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  nvar=nvar1;
  array_def *vars= (array_def *) calloc (nvar,sizeof(array_def));//recycle ha_cgeset
  bool *var_inter= (bool *) calloc (nvar,sizeof(bool));//recycle ha_cgeset
  logmsg(2,"nvarele %ld\n",nvarele);
  if(rank==0) {
    nvarele=variables_read(tabfile,commsyntax,vars,nvar,sets,nset);
    if(nvarele==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    /* 11.2.1 name uniqueness (the 12.2.2 name-resolution spec pass):
       a coefficient/variable/set name collision would silently bind
       whichever list is searched first */
    if(names_validate(sets,nset,coefs,ncof,vars,nvar,maps,nmap)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    if(nmap>0)mapping_use_guards(tabfile,maps,nmap);
    nvarele1=nvarele;
  }
  logmsg(2,"nvarele %ld\n",nvarele);
  if(nohsl) {
    MPI_Bcast(&nvarele1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(vars,nvar*sizeof(array_def), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  nvarele=nvarele1;
  /* coefficient X + variable p_X/c_X is the hand-linearized pair
     idiom (GTAP-AEZ YIELD/p_YIELD) and is now resolved unambiguously:
     p_-prefixed tokens bind variables first (incl. the declared
     p_-/c_-leading name via linvar_resolve), bare tokens bind
     coefficients first. The old guard here fataled the idiom; the
     genuine ambiguity (variable X + variable p_X/c_X) is fataled in
     names_validate (design doc section 6). */
  logmsg(2,"nvarele %ld\n",nvarele);
  elem_value *elem_vals= (elem_value *) calloc ((ncofele+nvarele),sizeof(elem_value));
  elem_store *coef_store= (elem_store *) calloc (ncofele,sizeof(elem_store));
  elem_store *var_store= (elem_store *) calloc (nvarele,sizeof(elem_store));
  /* element counts are products of the set sizes, and set sizes come from
     the data-file headers, so an oversized header count makes these
     allocations fail; coef_resolve_sets and data_read_files then wrote
     through NULL (fuzz batch 14 data-file driver: a REG count of 99999
     segfaulted in data_read_files). */
  if ((ncofele+nvarele)>0&&elem_vals==NULL) {
    errmsg("Error: cannot allocate %ld coefficient and variable elements; check the element counts the data-file headers declare\n",(long)(ncofele+nvarele));
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  if ((ncofele>0&&coef_store==NULL)||(nvarele>0&&var_store==NULL)) {
    errmsg("Error: cannot allocate the coefficient (%ld) and variable (%ld) value stores; check the element counts the data-file headers declare\n",(long)ncofele,(long)nvarele);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  logmsg(2,"rankasd %d nvar %ld\n",rank,nvar);
  if(rank==0) {
    coef_resolve_sets(coefs,ncof,sets,nset,coef_store);
  }
  if(nohsl) {
    if(ncofele*sizeof(elem_store)>1500000000) {
      j1=1500000000/sizeof(elem_store);
      i=ncofele/j1;
      for(j=0; j<i; j++) {
        MPI_Bcast(coef_store+j*j1,j1*sizeof(elem_store), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
      i=ncofele-j*j1;
      MPI_Bcast(coef_store+j*j1,i*sizeof(elem_store), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(coef_store,ncofele*sizeof(elem_store), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
  }
  logmsg(2,"rank %d ncofele %ld\n",rank,ncofele);

  if(rank==0)coef_resolve_sets(vars,nvar,sets,nset,var_store);

  if(rank==rank_hsl) {
    tab_write_variables(tabfile,newtabfile1,vars,nvar);
  }
  strcpy(tabfile,newtabfile1);
  if(nohsl) {
    MPI_Bcast(sets,nset*sizeof(set_def), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  //**************************************************************************************
  //****************************** END READ VARIABLE NAME*********************************
  //**************************************************************************************

  //**************************************************************************************
  //********************* READ VARIABLE, COEFFICIENT VALUE FROM FILE**********************
  //**************************************************************************************
  if(rank==rank_hsl)ord_plan_build(tabfile,0,coefs,ncof,vars,nvar);
  ord_io_set(niodata,iodata,nvarele);
  if(rank==0) {
    strcpy(commsyntax,"read");
    /* fail-fast: exit 0 here used to mask read errors */
    if(data_read_files(tabfile,niodata,iodata,commsyntax,sets,nset,set_elems,coefs,ncof,coef_store,ncofele,vars,nvar,var_store,nvarele)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
  }
  //**************************************************************************************
  //********************* END READ VARIABLE, COEFFICIENT VALUE FROM FILE******************
  //**************************************************************************************

  //**************************************************************************************
  //***************** CALCULATE VARIABLE, COEFFICIENT VALUE FROM FORMULA******************
  //**************************************************************************************

  offset_t nexo=0,nexo1;
  if(rank==0) {
    for (i=0; i<ncofele; i++) {
      elem_vals[i].value=coef_store[i].value;
    }
    for (i=ncofele; i<nvarele+ncofele; i++) {
      elem_vals[i].value=var_store[i-ncofele].value;
    }
  }
  logmsg(2,"rank %d ncofvar %ld\n",rank,ncofele+nvarele);
  free(coef_store);
  free(var_store);
  closure_entry *closure_vals= (closure_entry *) calloc (nvarele,sizeof(closure_entry));
  teems_cl_flags= (unsigned char *) calloc (nvarele,sizeof(unsigned char));
  teems_cl_shock= (store_real *) calloc (nvarele,sizeof(store_real));
  /* element-level border marks (6.5 E3), populated by the ordering scan
     alongside var_inter on the same ranks */
  bool *ele_inter= (bool *) calloc (nvarele,sizeof(bool));
  if(rank==0&&solmethod!=SM_NOSIM) {
    strcpy(commsyntax,"exogenous");
    nexo=closure_read(closure,commsyntax,closure_vals,vars,nvar,sets,nset,set_elems);
    if(nexo==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    nexo1=nexo;
    strcpy(commsyntax,"shock");
    if(shocks_read(shock,commsyntax,closure_vals,nvarele,vars,nvar,sets,nset,set_elems,subints)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    if(shocks_check_floor(vars,nvar,subints,solmethod)<0)MPI_Abort(PETSC_COMM_WORLD,1);
    if(has_sub) {
      if(subtotals_read(subfile,vars,nvar,sets,nset,set_elems,nvarele)<0)MPI_Abort(PETSC_COMM_WORLD,1);
      logmsg(1,"Subtotals: %d shock group(s) from %s (GEMPACK manual 29)\n",teems_nsub,subfile);
    }
    /* backsolve statements: mark the eliminated elements (the flags ride
       the closure broadcast) and check the condensed system's references
       before any equation scan runs with the filter active */
    if(backsolve_read(tabfile,vars,nvar,closure_vals)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    if(backsolve_validate_refs(tabfile,vars)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    if(nbacksolve>0)logmsg(1,"Backsolving %d variables (%ld elements) from retained defining equations\n",nbacksolve,nbselems);
    /* C1/C2: auto-exogenize del_comp@ and the dummies of ACTIVE
       (X-endogenous) complementarity components; inert components
       keep their dummy endogenous so it absorbs the E_$comp row
       (51.7.2 (c)/(d); design doc sections 7-8). The marks and the
       nexo adjustment ride the closure broadcast below; the 11.14.1
       backsolve guards stay fatal. */
    if(comp_closure_check(closure_vals,vars,nvar,&nexo,sets,nset,set_elems)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    nexo1=nexo;
  }
  /* every rank dispatches on active-mode complementarities (C2) */
  MPI_Bcast(&teems_comp_active,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  if(has_sub) {
    int hc=(teems_ncomp>0);
    MPI_Bcast(&hc,1,MPI_INT,0,PETSC_COMM_WORLD);
    if(hc) {
      if(rank==0)errmsg("Error: subtotals (GEMPACK manual 29) are not available in a model with complementarities yet: the approximate and accurate runs change the closure and the states between steps, so the step right-hand sides are not the shocks alone (GEMPACK manual 52)\n");
      MPI_Barrier(PETSC_COMM_WORLD);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
  /* approximate/accurate-run controls (manual 51.6 semantics, TEEMS
     command-line flags -- teems-R passes them from
     ems_complementarity(); no CMF statements). Defaults: step count
     = the accurate method's step sum, redo on, both runs on, errors
     fatal. PetscOptions are identical on every rank, so no
     broadcasts. comp_acc_phase: 0 = approximate pass (or first
     pass), 1 = the accurate pass after the 51.7.1 closure
     modification (C3). */
  int comp_steps=0,comp_redo=1,comp_do_approx=1,comp_do_acc=1,comp_sberr_warn=0,comp_acc_phase=0;
  /* 51.7.4: with several subintervals each one runs its own approximate
     and accurate pair. The original closure and per-subinterval shares
     are kept to restore between pairs; results compound across them. */
  int comp_sub=0,comp_next_phase=0;
  unsigned char *comp_flags0=NULL;
  store_real *comp_shock0=NULL;
  solve_real *comp_tot=NULL,*comp_apx_tot=NULL,*comp_est_tot=NULL;
  double comp_minfrac=0.005;
  if(teems_comp_active>0) {
    PetscReal minfrac_opt=comp_minfrac;
    dim_t iopt;
    /* 51.6 default = the accurate run's step sum; steps2/steps3 were
       folded into ratios above, so rebuild the three counts the
       multistep driver will use */
    /* -single_run runs steps1 only, -two_run steps1 and steps2 */
    comp_steps=(solmethod==SM_GRAGG||solmethod==SM_MIDPOINT||solmethod==SM_EULER)
      ?steps1+(teems_single_run?0:(int)llround(steps1*step_ratio2))+((teems_single_run||teems_two_run)?0:(int)llround(steps1*step_ratio3))
      :steps1;
    if(comp_steps<1)comp_steps=10;
    iopt=comp_steps;
    PetscOptionsGetInt(NULL,NULL,"-comp_steps",&iopt,NULL);
    if(iopt<1) {
      if(rank==0)errmsg("Error: -comp_steps must be a positive Euler step count (GEMPACK manual 51.6)\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    comp_steps=(int)iopt;
    iopt=comp_redo;
    PetscOptionsGetInt(NULL,NULL,"-comp_redo",&iopt,NULL);
    comp_redo=iopt?1:0;
    PetscOptionsGetReal(NULL,NULL,"-comp_redo_min_frac",&minfrac_opt,NULL);
    if(minfrac_opt<=0||minfrac_opt>1) {
      if(rank==0)errmsg("Error: -comp_redo_min_frac must lie in (0,1] (GEMPACK manual 51.6)\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    comp_minfrac=(double)minfrac_opt;
    iopt=comp_do_approx;
    PetscOptionsGetInt(NULL,NULL,"-comp_do_approx",&iopt,NULL);
    comp_do_approx=iopt?1:0;
    iopt=comp_do_acc;
    PetscOptionsGetInt(NULL,NULL,"-comp_do_acc",&iopt,NULL);
    comp_do_acc=iopt?1:0;
    teems_comp_no_acc=!comp_do_acc;
    iopt=comp_sberr_warn;
    PetscOptionsGetInt(NULL,NULL,"-comp_sberr_warn",&iopt,NULL);
    comp_sberr_warn=iopt?1:0;
    if(!comp_do_approx&&!comp_do_acc) {
      if(rank==0)errmsg("Error: -comp_do_approx and -comp_do_acc cannot both be 0\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    /* state flips change the filtered nonzero pattern between steps,
       which breaks the persistent-pivot refactorization */
    {
      dim_t frchk=0;
      PetscOptionsGetInt(NULL,NULL,"-fastrefac",&frchk,NULL);
      if(frchk) {
        if(rank==0)printf("Warning: -fastrefac is disabled for the complementarity approximate run (state changes alter the nonzero pattern between steps)\n");
        PetscOptionsSetValue(NULL,"-fastrefac","0");
      }
    }
  }
  /* a Runge-Kutta run marches the whole interval under its own step
     control, so subintervals only mean something where they re-predict
     complementarity states: an approximate and an accurate pair per
     subinterval (51.7.4), the accurate run by the Runge-Kutta method
     with -step1 steps in each. Decided here, once the closure says
     which complementarities are active. */
  if(isrk&&subints>1&&!(teems_comp_active>0&&comp_do_acc)) {
    if(rank==0)errmsg("Error: subintervals are not available with Runge-Kutta methods (got -nsubints %d) except in complementarity runs with an accurate run; increase -step1 instead\n",subints);
    PetscFinalize();
    return 1;
  }
  if(nohsl) {
    if(nvarele*sizeof(closure_entry)>1500000000) {
      j1=1500000000/sizeof(closure_entry);
      i=nvarele/j1;
      for(j=0; j<i; j++) {
        MPI_Bcast(closure_vals+j*j1,j1*sizeof(closure_entry), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
      i=nvarele-j*j1;
      MPI_Bcast(closure_vals+j*j1,i*sizeof(closure_entry), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(closure_vals,nvarele*sizeof(closure_entry), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    /* the closure side arrays travel with closure_vals (same chunking) */
    j1=1500000000/sizeof(store_real);
    for(j=0; j*j1<nvarele; j++) {
      i=(nvarele-j*j1<j1)?nvarele-j*j1:j1;
      MPI_Bcast(teems_cl_flags+j*j1,i*sizeof(unsigned char), MPI_BYTE,0, PETSC_COMM_WORLD);
      MPI_Bcast(teems_cl_shock+j*j1,i*sizeof(store_real), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    MPI_Bcast(&nexo1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nbacksolve,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nbselems,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    if(nbacksolve>0) {
      if(rank!=0)backsolves=realloc(backsolves,nbacksolve*sizeof(backsolve_def));
      MPI_Bcast(backsolves,nbacksolve*sizeof(backsolve_def), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
  }
  nexo=nexo1;
  strcpy(commsyntax,"formula");
  bool IsIni=true;
  if(rank==0) {
    statements_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,0);
    if(cofdump&&elem_vals!=NULL)coefficients_dump_phase(teems_sol_stem,".cbin0",0,ncof,ncofele,elem_vals); /* reported, not fatal: the run goes on and exits 1 without a completion marker */
  }
  /* formula-assigned mappings (manual 10.13.1) got their values on rank
     0 during the pass above; the equation side on every rank routes
     mapped indices through the same tables */
  if(nohsl&&nmap>0) {
    for(i=0; i<nmap; i++) if(maps[i].formula_assigned) {
      MPI_Bcast(&maps[i].has_values,sizeof(bool), MPI_BYTE,0, PETSC_COMM_WORLD);
      MPI_Bcast(&maps[i].nassigned,sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
      MPI_Bcast(maps[i].values,sets[maps[i].fromset].size*sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
      MPI_Bcast(maps[i].assigned,sets[maps[i].fromset].size*sizeof(unsigned char), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
  }
  /* Tier B4 no-simulation finish: the values after the Reads, Formulas
     and Assertions are the run's values, written the way a solve writes
     its post-simulation ones (structure files, the .cof/.cbin dump, the
     TAB's Writes); no .bin, no stats.json. PostSim needs a simulation
     and is not run. */
  if(solmethod==SM_NOSIM) {
    if(rank==0) {
      if(npostsim>0||tab_has_postsim_assertions(tabfile))
        printf("Warning: PostSim statements and assertions are not run without a simulation (GEMPACK manual 12.1); %d PostSim statement(s) skipped\n",npostsim);
      if(structure_files_write(teems_sol_stem,vars,nvar,sets,nset,set_elems,nsetspace,nvarele)==0) {
        int cofdumped=0,wr;
        long nsets_w=0,nother_w=0,nskip_w=0;
        if(cofdump&&elem_vals!=NULL) {
          if(coefficients_dump(teems_sol_stem,coefs,ncof,ncofele,elem_vals)==0)cofdumped=1;
        }
        if(nowrites==0)for(i=0; i<noutdata; i++) {
          if(postsim_write_skipped(iodata[i+niodata].logname)) { nskip_w++; continue; }
          wr=outputs_write_csv(tabfile,iodata[i+niodata].logname,iodata[i+niodata].filname,sets,nset,set_elems,coefs,ncof,ncofele,vars,nvar,nvarele,elem_vals);
          if(wr!=-1) {
            if(wr==1)nsets_w++; else nother_w++;
            outputs_note(iodata[i+niodata].filname,"csv",0);
          }
        }
        outputs_summary_log(cofdumped,coefs,ncof,ncofele,nsets_w,nother_w,nskip_w,nowrites==0,0,nvar,nvarele);
      }
      gettimeofday(&endtime, NULL);
      logmsg(1,"No-simulation run complete in %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    }
    {
      int nerr=teems_error_count,nerrsum=0;
      MPI_Allreduce(&nerr,&nerrsum,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
      if(rank==0&&nerrsum==0)outputs_json_write(teems_sol_stem,run_id);
    }
    free(iodata);
    free(sets);
    free(set_elems);
    free(coefs);
    free(vars);
    free(var_inter);
    free(ele_inter);
    free(closure_vals);
    free(elem_vals);
    free(teems_cl_flags);
    free(teems_cl_shock);
    free(teems_set_isprod);
    free(teems_set_prod1);
    free(teems_set_prod2);
    free(teems_set_itstem);
    free(teems_set_itfirst);
    MPI_Comm_free(&node_comm);
    MPI_Comm_free(&node_tail_comm);
    PetscFinalize();
    return teems_error_count>0?1:0;
  }
  gettimeofday(&endtime, NULL);
  if(rank==0)logmsg(1,"Variable calculation time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
  teems_rss_probe("variable calculation");
  if(nohsl) { //Overcome MPI_Bcast limit
    if((nvarele+ncofele)*sizeof(elem_value)>1500000000) {
      j1=1500000000/sizeof(elem_value);
      i=(nvarele+ncofele)/j1;
      for(j=0; j<i; j++) {
        MPI_Bcast(elem_vals+j*j1,j1*sizeof(elem_value), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
      i=nvarele+ncofele-j*j1;
      MPI_Bcast(elem_vals+j*j1,i*sizeof(elem_value), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(elem_vals,(nvarele+ncofele)*sizeof(elem_value), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    ord_coverage_bcast(tabfile,rank);
  }
  //**************************************************************************************
  //****************** END CALCULATE VARIABLE, COEFFICIENT VALUE FROM FORMULA*************
  //**************************************************************************************
  gettimeofday(&begintime, NULL);
  if(rank==0)logmsg(1,"Variable broadcast time %.2f s\n",(begintime.tv_sec - endtime.tv_sec)+((double)(begintime.tv_usec - endtime.tv_usec))/ 1000000);
  teems_rss_probe("variable broadcast");
  //**************************************************************************************
  //****************************** MATRIX FROM FORMULA************************************
  //**************************************************************************************
  VecSize = (PetscInt) (nvarele-nexo-nbselems);
  if(rank==0) {
    if(nbselems>0)logmsg(1,"System size %d equations (%ld exogenous, %ld backsolved)\n",VecSize,nexo,nbselems);
    else logmsg(1,"System size %d equations (%ld exogenous)\n",VecSize,nexo);
  }
  strcpy(commsyntax,"equation");
  offset_t neq=0,neq1;
  if(rank==0) {
    neq=tab_count_statements(tabfile,commsyntax);
    neq1=neq;
  }
  if(nohsl) {
    MPI_Bcast(&neq1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  neq=neq1;
  if(rank==rank_hsl) {
    logmsg(2,"neq %ld\n",neq);
  }
  /* ---------------- structural partition resolution ----------------
     Derive the ordering dimensions from the equation system itself:
     the chain dimension from the lead/lag offsets the equations
     actually use, and the diagonal-block partition from probing every
     structurally eligible set. Runs on every rank that performs the
     ordering (all ranks under nohsl, rank 0 under HSL) as pure
     integer analysis of identical data, so all ranks reach the same
     result without communication. */
  const char *chain_source="none",*partition_source="none";
  char *partition_auto_json=NULL;
  if((structural_time||structural_reg)&&rank==rank_hsl) {
    if(structural_time) {
      offset_t *chainrefs=(offset_t *) calloc (nset,sizeof(offset_t));
      chain_refs_scan(tabfile,sets,nset,coefs,ncof,vars,nvar,elem_vals,chainrefs);
      alltimeset=chain_set_select(sets,nset,chainrefs,rank);
      free(chainrefs);
      if(alltimeset>=0) {
        chain_flags_apply(sets,nset,alltimeset);
        chain_source="structural";
      }
    }
    if(alltimeset>=0)ntime=sets[alltimeset].size;
    /* the probe measures the nested (NDBBD) partition geometry when a
       chain exists, the flat (DBBD) one otherwise; the -matsol it was
       launched with does not decide */
    if(solmethod==SM_PROBE)nesteddbbd=(alltimeset>=0)?1:0;
    if(structural_reg) {
      allregset=partition_auto_select(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,nvarele,closure_vals,neq,(offset_t)VecSize,alltimeset,ntime,(alltimeset>=0)?nesteddbbd:0,(long)mpisize,rank,&partition_auto_json);
      if(allregset>=0) {
        partition_flags_apply(sets,nset,allregset);
        partition_source="structural";
      }
    }
    if(allregset>=0)nreg=sets[allregset].size;
    /* rebuild the block geometry from the resolved dimensions */
    if(alltimeset>=0&&allregset>=0) {
      if(nesteddbbd==1)ndblock=ntime*(nreg+1);
      else ndblock=ntime*nreg;
    }
    else if(alltimeset>=0)ndblock=ntime;
    else if(allregset>=0)ndblock=nreg;
    ndblock1=ndblock;
  }
  if(solmethod!=SM_PROBE&&matsol==MM_NDBBD&&(alltimeset<0||allregset<0)) {
    if(rank==0) {
      if(alltimeset<0)errmsg("Error: NDBBD (-matsol 3) needs a chain dimension, but no equation couples set elements through lead/lag offsets (or the set that carries them has one element); use -matsol 2 (DBBD) or -matsol 0 (LU).\n");
      if(allregset<0)errmsg("Error: NDBBD (-matsol 3) needs a diagonal-block partition and no viable set was detected (see the candidate table above); choose another -matsol.\n");
    }
    PetscFinalize();
    return 1;
  }
  if(solmethod!=SM_PROBE&&matsol==MM_DBBD&&alltimeset<0&&allregset<0) {
    if(rank==0)errmsg("Error: DBBD (-matsol 2) needs a diagonal-block partition and no viable set was detected (see the candidate table above); use -matsol 0 (LU).\n");
    PetscFinalize();
    return 1;
  }
  /* method-vs-structure check: SBBD without a chain dimension hands
     HSL_MP48 zero blocks — the factorization errors out and the run
     used to finish with exit 0 and no solution. Abort cleanly instead.
     Under HSL only rank 0 resolves alltimeset, so the verdict is
     broadcast to keep the exit collective. */
  if(solmethod!=SM_PROBE&&matsol==MM_SBBD) {
    int sbbd_nochain=(rank==0&&alltimeset<0)?1:0;
    MPI_Bcast(&sbbd_nochain,1,MPI_INT,0,PETSC_COMM_WORLD);
    if(sbbd_nochain) {
      if(rank==0)errmsg("Error: SBBD (-matsol 1) requires a chain dimension, but the equations couple no set through lead/lag offsets (or the set that carries them has one element); use -matsol 0 (LU) or -matsol 2 (DBBD) for static and one-period models.\n");
      PetscFinalize();
      return 1;
    }
  }
  /* C3 accurate-run re-entry (design doc section 8 tail): after the
     approximate pass the 51.7.1 closure modification changes WHICH
     components are exogenous (never how many), so everything downstream
     of the closure -- border classification, exo_index numbering,
     block/equation addressing, preallocation, the shock vector -- is
     rebuilt by re-running this section. All ranks take the jump
     together (comp_acc_phase is broadcast-derived). */
comp_accurate_reentry:
  /* the approximate driver leaves commsyntax at "formula"; the
     ordering below must parse equations */
  strcpy(commsyntax,"equation");
  if(nesteddbbd==1)ndbbddrank1=(PetscInt *) calloc(ntime,sizeof(PetscInt));
  offset_t *countvarintra1= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  array_def *eq_defs= (array_def *) calloc (neq,sizeof(array_def));//recycle ha_cgeset
  bool *eq_intertemp= (bool *) calloc (neq,sizeof(bool));//recycle ha_cgeset
  dim_t *orderintra= (dim_t *) malloc (nvar*sizeof(dim_t));
  dim_t *orderreg= (dim_t *) malloc (nvar*sizeof(dim_t));
  for(i=0; i<nvar; i++) {
    orderintra[i]=-1;
    orderreg[i]=-1;
  }
  dim_t *eq_time= (dim_t *) malloc (neq*sizeof(dim_t));//recycle ha_cgeset
  dim_t *eq_reg= (dim_t *) malloc (neq*sizeof(dim_t));//recycle ha_cgeset
  for(i=0; i<neq; i++) {
    eq_time[i]=-1;
    eq_reg[i]=-1;
  }
  offset_t nintraendovar;
  if(rank==rank_hsl) {
    if(nesteddbbd==1) {
      if(!equation_order_read_nested(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,var_inter,ele_inter,eq_defs,eq_intertemp,eq_time,eq_reg,allregset,alltimeset,orderintra,orderreg))MPI_Abort(PETSC_COMM_WORLD,1);
    }
    else if(!equation_order_read(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,var_inter,ele_inter,eq_defs,eq_intertemp,eq_time,eq_reg,allregset,alltimeset,orderintra,orderreg))MPI_Abort(PETSC_COMM_WORLD,1);
    if(alltimeset>=0||allregset>=0)for(i=0; i<neq; i++)eq_intertemp[i]=!eq_intertemp[i];
    /* initial closure check (manual 23.2.7): the number of endogenous
       components must equal the number of equation rows, checked
       before any matrix is sized from it (a non-square closure used to
       surface as an unnamed MA48/PETSc failure or an overrun) */
    {
      offset_t neqrows=0;
      for(i=0; i<neq; i++)neqrows+=eq_defs[i].nelem;
      if(neqrows!=(offset_t)VecSize) {
        errmsg("Error: initial closure check: %ld endogenous components is not equal to the number of equation rows (%ld); %ld variable components, %ld exogenous, %ld backsolved -- make %ld more component(s) %s (GEMPACK manual 23.2.7)\n",(long)VecSize,(long)neqrows,(long)nvarele,(long)nexo,(long)nbselems,(long)(neqrows>(offset_t)VecSize?neqrows-(offset_t)VecSize:(offset_t)VecSize-neqrows),neqrows>(offset_t)VecSize?"endogenous":"exogenous");
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
    }
    j0=0;
    for(i=0; i<nvarele; i++)if(ele_inter[i])j0++;
    if(j0>0&&rank==0)logmsg(1,"%sElement-level border classification: %ld elements bordered by sliced lead/lag references\n",PROBE_PFX,j0);
  }
  switch (nesteddbbd) {
  case 1 : ;/* missing partition sets already abort above (NDBBD requires both) */
    offset_t *countvarintra= (offset_t *) calloc (ndblock,sizeof(offset_t));
    block_var_count(vars,nvar,sets,set_elems,closure_vals,var_inter,ele_inter,orderintra,orderreg,alltimeset,allregset,nreg,nesteddbbd,countvarintra);
    if(rank==rank_hsl) {
      j2=countvarintra[0];
      countvarintra[0]=0;
      nintraendovar=j2;
    }
    for (i=1; i<ndblock; i++) {
      nintraendovar+=countvarintra[i];
      j3=countvarintra[i];
      countvarintra[i]=countvarintra[i-1]+j2;
      countvarintra1[i]=countvarintra[i];
      j2=j3;
    }
    j0=0;
    j1=0;
    j2=0;
    j3=0;
    j4=nintraendovar;
    for (i=0; i<nvar; i++) {
      for (j=0; j<vars[i].nelem; j++) {
        j5=j3+j;
        if(!CL_EXO(j5)&&!CL_BS(j5)) {
          if(!var_inter[i]&&!ele_inter[j5]) {
            j0=j;
            j2=-1;
            for(j1=0; j1<orderintra[i]+1; j1++) {
              j2=j0/vars[i].strides[j1];
              j0-=j2*vars[i].strides[j1];
            }
            j0=j;
            j6=-1;
            for(j1=0; j1<orderreg[i]+1; j1++) {
              j6=j0/vars[i].strides[j1];
              j0-=j6*vars[i].strides[j1];
            }
            if(sets[vars[i].setid[orderintra[i]]].intsup>0)j2=set_elems[sets[vars[i].setid[orderintra[i]]].offset+j2].superset_pos[sets[vars[i].setid[orderintra[i]]].intsup];
            if(j6>-1)if(sets[vars[i].setid[orderreg[i]]].regsup>0)j6=set_elems[sets[vars[i].setid[orderreg[i]]].offset+j6].superset_pos[sets[vars[i].setid[orderreg[i]]].regsup];
            if(orderreg[i]>-1) {
              closure_vals[j5].exo_index=countvarintra[j2*(nreg+1)+j6];
              countvarintra[j2*(nreg+1)+j6]++;
            }
            else {
              closure_vals[j5].exo_index=countvarintra[j2*(nreg+1)+nreg];
              countvarintra[j2*(nreg+1)+nreg]++;
            }
          }
          else {
            closure_vals[j5].exo_index=j4;
            j4++;
          }
        }
      }
      j3+=vars[i].nelem;
    }
    if(rank==rank_hsl) {
      countvarintra1[ndblock]=countvarintra[ndblock-1];
    }
    j1=0;
    for (i=0; i<nvarele; i++) {
      if (CL_EXO(i)) {
        closure_vals[i].exo_index+=j1;
        j1++;
      }
    }
    free(countvarintra);

    break;
  default :
    if(alltimeset>=0) {
      offset_t *countvarintra= (offset_t *) calloc (ndblock,sizeof(offset_t));
      block_var_count(vars,nvar,sets,set_elems,closure_vals,var_inter,ele_inter,orderintra,orderreg,alltimeset,allregset,nreg,nesteddbbd,countvarintra);
      if(rank==rank_hsl) {
        j2=countvarintra[0];
        countvarintra[0]=0;
        nintraendovar=j2;
      }
      for (i=1; i<ndblock; i++) {
        nintraendovar+=countvarintra[i];
        j3=countvarintra[i];
        countvarintra[i]=countvarintra[i-1]+j2;
        countvarintra1[i]=countvarintra[i];
        j2=j3;
      }
      j0=0;
      j1=0;
      j2=0;
      j3=0;
      j4=nintraendovar;
      for (i=0; i<nvar; i++) {
        for (j=0; j<vars[i].nelem; j++) {
          j5=j3+j;
          if(!CL_EXO(j5)&&!CL_BS(j5)) {
            if(!var_inter[i]&&!ele_inter[j5]) {
              j0=j;
              for(j1=0; j1<orderintra[i]+1; j1++) {
                j2=j0/vars[i].strides[j1];
                j0-=j2*vars[i].strides[j1];
              }
              if(allregset>=0) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j6=j0/vars[i].strides[j1];
                  j0-=j6*vars[i].strides[j1];
                }
                if(sets[vars[i].setid[orderreg[i]]].regsup>0)j6=set_elems[sets[vars[i].setid[orderreg[i]]].offset+j6].superset_pos[sets[vars[i].setid[orderreg[i]]].regsup];
                closure_vals[j5].exo_index=countvarintra[j2*sets[allregset].size+j6];
                countvarintra[j2*sets[allregset].size+j6]++;
              }
              else {
                closure_vals[j5].exo_index=countvarintra[j2];
                countvarintra[j2]++;
              }
            }
            else {
              closure_vals[j5].exo_index=j4;
              j4++;
            }
          }
        }
        j3+=vars[i].nelem;
      }
      if(rank==rank_hsl) {
        countvarintra1[ndblock]=countvarintra[ndblock-1];
      }
      j1=0;
      for (i=0; i<nvarele; i++) {
        if (CL_EXO(i)) {
          closure_vals[i].exo_index+=j1;
          j1++;
        }
      }
      free(countvarintra);
    }
    else {
      if(allregset>=0) {
        offset_t *countvarintra= (offset_t *) calloc (ndblock,sizeof(offset_t));
        block_var_count(vars,nvar,sets,set_elems,closure_vals,var_inter,ele_inter,orderintra,orderreg,alltimeset,allregset,nreg,nesteddbbd,countvarintra);
        if(rank==rank_hsl) {
          j2=countvarintra[0];
          countvarintra[0]=0;
          nintraendovar=j2;
        }
        for (i=1; i<ndblock; i++) {
          nintraendovar+=countvarintra[i];
          j3=countvarintra[i];
          countvarintra[i]=countvarintra[i-1]+j2;
          countvarintra1[i]=countvarintra[i];
          j2=j3;
        }
        j0=0;
        j1=0;
        j2=0;
        j3=0;
        j4=nintraendovar;
        j6=0; /* flow analysis only: the stride loop assigns it, orderreg[i]>=0 for every intra variable */
        for (i=0; i<nvar; i++) {
          for (j=0; j<vars[i].nelem; j++) {
            j5=j3+j;
            if(!CL_EXO(j5)&&!CL_BS(j5)) {
              if(!var_inter[i]&&!ele_inter[j5]) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j6=j0/vars[i].strides[j1];
                  j0-=j6*vars[i].strides[j1];
                }
                if(sets[vars[i].setid[orderreg[i]]].regsup>0)j6=set_elems[sets[vars[i].setid[orderreg[i]]].offset+j6].superset_pos[sets[vars[i].setid[orderreg[i]]].regsup];
                closure_vals[j5].exo_index=countvarintra[j6];
                countvarintra[j6]++;
              }
              else {
                closure_vals[j5].exo_index=j4;
                j4++;
              }
            }
          }
          j3+=vars[i].nelem;
        }
        if(rank==rank_hsl) {
          countvarintra1[ndblock]=countvarintra[ndblock-1];
        }
        j1=0;
        for (i=0; i<nvarele; i++) {
          if (CL_EXO(i)) {
            closure_vals[i].exo_index+=j1;
            j1++;
          }
        }
        free(countvarintra);
      }
      else {
        j1=0;
        j2=0;
        for (i=0; i<nvar; i++) {
          for (j=0; j<vars[i].nelem; j++) {
            j3=j0+j;

            if (!CL_EXO(j3)&&!CL_BS(j3)) {
              closure_vals[j3].exo_index+=j2;
              j2++;
            }
            if (CL_EXO(j3)) {
              closure_vals[j3].exo_index+=j1;
              j1++;
            }
          }
          j0+=vars[i].nelem;
        }
      }
    }
  }
  free(orderintra);
  free(orderreg);

  free(var_inter);
  free(ele_inter);
  strcpy(commsyntax,"equation");
  PetscInt *eq_addr= (PetscInt *) calloc (VecSize,sizeof(PetscInt)); /* row index per equation position: PetscInt by meaning, half the offset_t footprint (6.16(b)) */
  offset_t *eq_time_offsets= (offset_t *) calloc (neq,sizeof(offset_t));
  offset_t *eq_reg_offsets= (offset_t *) calloc (neq,sizeof(offset_t));
  offset_t *counteq= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  offset_t *counteqnoadd= (offset_t *) calloc (ndblock,sizeof(offset_t));
  offset_t nintraeq=0;
  if(alltimeset>=0&&allregset<0) {
    for(i=0; i<neq; i++) {
      j3=1;
      if(eq_time[i]>-1)eq_time_offsets[i]=sets[eq_defs[i].setid[eq_time[i]]].offset;
      /* scalar equations (size 0) wrote strides[-1] = setid[9] */
      if(eq_defs[i].size>0)eq_defs[i].strides[eq_defs[i].size-1]=1;
      if(eq_defs[i].size>1) {
        for (j2=eq_defs[i].size-2; j2>-1; j2--) {
          eq_defs[i].strides[j2]=eq_defs[i].strides[j2+1]*sets[eq_defs[i].setid[j2+1]].size;
        }
      }
    }
    j3=0;
    for (i=0; i<neq; i++) {
      if(eq_intertemp[i]) {
        for (j=0; j<eq_defs[i].nelem; j++) {
          j0=j;
          for(j1=0; j1<eq_time[i]+1; j1++) {
            j2=j0/eq_defs[i].strides[j1];
            j0-=j2*eq_defs[i].strides[j1];
          }
          if(eq_time[i]>-1){ if(eq_defs[i].setid[eq_time[i]]==alltimeset)counteq[set_elems[eq_time_offsets[i]+j2].superset_pos[0]]++;
            else {
              for(j3=1; j3<MAXSUPSET; j3++)if(sets[eq_defs[i].setid[eq_time[i]]].subsetid[j3]==alltimeset)break;
              counteq[set_elems[eq_time_offsets[i]+j2].superset_pos[j3]]++;
            } }
        }
      }
    }
    if(rank==rank_hsl) {
      j2=counteq[0];
      counteqnoadd[0]=counteq[0];
      counteq[0]=0;
      nintraeq=j2;
    }
    for (i=1; i<ndblock; i++) {
      counteqnoadd[i]=counteq[i];
      nintraeq+=counteq[i];
      j3=counteq[i];
      counteq[i]=counteq[i-1]+j2;
      j2=j3;
    }
    if(rank==rank_hsl) {
      counteq[ndblock]=VecSize;
      netcut=VecSize-countvarintra1[ndblock];
    }
  }
  if(alltimeset<0&&allregset>=0) {
    for(i=0; i<neq; i++) {
      j3=1;
      if(eq_reg[i]>-1)eq_time_offsets[i]=sets[eq_defs[i].setid[eq_reg[i]]].offset;
      /* scalar equations (size 0) wrote strides[-1] = setid[9] */
      if(eq_defs[i].size>0)eq_defs[i].strides[eq_defs[i].size-1]=1;
      if(eq_defs[i].size>1) {
        for (j2=eq_defs[i].size-2; j2>-1; j2--) {
          eq_defs[i].strides[j2]=eq_defs[i].strides[j2+1]*sets[eq_defs[i].setid[j2+1]].size;
        }
      }
    }
    j3=0;
    for (i=0; i<neq; i++) {
      if(eq_intertemp[i]) { //for (j=0; j<ha_set[allregset].size; j++)counteq[j]+=(uvadd)ha_eq[i].matsize/ha_set[allregset].size;
        for (j=0; j<eq_defs[i].nelem; j++) {
          j0=j;
          for(j1=0; j1<eq_reg[i]+1; j1++) {
            j2=j0/eq_defs[i].strides[j1];
            j0-=j2*eq_defs[i].strides[j1];
          }
          counteq[set_elems[eq_time_offsets[i]+j2].superset_pos[0]]++;
        }
      }
    }
    if(rank==rank_hsl) {
      counteqnoadd[0]=counteq[0];
      j2=counteq[0];
      counteq[0]=0;
      nintraeq=j2;
    }
    for (i=1; i<ndblock; i++) {
      counteqnoadd[i]=counteq[i];
      nintraeq+=counteq[i];
      j3=counteq[i];
      counteq[i]=counteq[i-1]+j2;
      j2=j3;
    }
    if(rank==rank_hsl) {
      counteq[ndblock]=VecSize;
      netcut=VecSize-countvarintra1[ndblock];
    }
  }

  if(alltimeset>=0&&allregset>=0) {
    switch (nesteddbbd) {
    case 1 :
      for (i=0; i<neq; i++) {
        if(eq_time[i]>-1)eq_time_offsets[i]=sets[eq_defs[i].setid[eq_time[i]]].offset;
        if(eq_reg[i]>-1)eq_reg_offsets[i]=sets[eq_defs[i].setid[eq_reg[i]]].offset;
        if(eq_intertemp[i]) {
          for (j=0; j<sets[eq_defs[i].setid[eq_time[i]]].size; j++)
            if(eq_reg[i]>-1)for(j1=0; j1<sets[eq_defs[i].setid[eq_reg[i]]].size; j1++)
                counteq[set_elems[eq_time_offsets[i]+j].superset_pos[sets[eq_defs[i].setid[eq_time[i]]].intsup]*(nreg+1)+set_elems[eq_reg_offsets[i]+j1].superset_pos[sets[eq_defs[i].setid[eq_reg[i]]].regsup]]+=eq_defs[i].nelem/sets[eq_defs[i].setid[eq_time[i]]].size/sets[eq_defs[i].setid[eq_reg[i]]].size;
            else counteq[set_elems[eq_time_offsets[i]+j].superset_pos[sets[eq_defs[i].setid[eq_time[i]]].intsup]*(nreg+1)+nreg]+=eq_defs[i].nelem/sets[eq_defs[i].setid[eq_time[i]]].size;
        }
      }
      if(rank==rank_hsl) {
        counteqnoadd[0]=counteq[0];
        j2=counteq[0];
        counteq[0]=0;
        nintraeq=j2;
      }
      for (i=1; i<ndblock; i++) {
        counteqnoadd[i]=counteq[i];
        nintraeq+=counteq[i];
        j3=counteq[i];
        counteq[i]=counteq[i-1]+j2;
        j2=j3;
      }
      if(rank==rank_hsl) {
        counteq[ndblock]=VecSize;//Attention!!!!!!!!!!!! Different from countvarintra1. Unchanged for not affecting previous method
        netcut=VecSize-countvarintra1[ndblock];
      }
      break;

    default :
      for(i=0; i<neq; i++) {
        if(eq_time[i]>-1)eq_time_offsets[i]=sets[eq_defs[i].setid[eq_time[i]]].offset;
        if(eq_reg[i]>-1)eq_reg_offsets[i]=sets[eq_defs[i].setid[eq_reg[i]]].offset;
      }
      for (i=0; i<neq; i++) {
        if(eq_intertemp[i]) {
          for (j=0; j<sets[eq_defs[i].setid[eq_time[i]]].size; j++)
            if(eq_defs[i].setid[eq_time[i]]==alltimeset)
              for(j1=0; j1<sets[eq_defs[i].setid[eq_reg[i]]].size; j1++)
                counteq[set_elems[eq_time_offsets[i]+j].superset_pos[0]*sets[eq_defs[i].setid[eq_reg[i]]].size+set_elems[eq_reg_offsets[i]+j1].superset_pos[0]]+=eq_defs[i].nelem/sets[eq_defs[i].setid[eq_time[i]]].size/sets[eq_defs[i].setid[eq_reg[i]]].size;
            else {
              for(j3=1; j3<MAXSUPSET; j3++)if(sets[eq_defs[i].setid[eq_time[i]]].subsetid[j3]==alltimeset)break;
              for(j1=0; j1<sets[eq_defs[i].setid[eq_reg[i]]].size; j1++)
                counteq[set_elems[eq_time_offsets[i]+j].superset_pos[j3]*sets[eq_defs[i].setid[eq_reg[i]]].size+set_elems[eq_reg_offsets[i]+j1].superset_pos[0]]+=eq_defs[i].nelem/sets[eq_defs[i].setid[eq_time[i]]].size/sets[eq_defs[i].setid[eq_reg[i]]].size;
            }
        }
      }
      if(rank==rank_hsl) {
        counteqnoadd[0]=counteq[0];
        j2=counteq[0];
        counteq[0]=0;
        nintraeq=j2;
      }
      for (i=1; i<ndblock; i++) {
        counteqnoadd[i]=counteq[i];
        nintraeq+=counteq[i];
        j3=counteq[i];
        counteq[i]=counteq[i-1]+j2;
        j2=j3;
      }
      if(rank==rank_hsl) {
        counteq[ndblock]=VecSize;
        netcut=VecSize-countvarintra1[ndblock];
      }
    }
  }
  if(rank==rank_hsl) {
    logmsg(1,"%sBorder netcut %ld, intra-block equations %ld\n",PROBE_PFX,netcut,nintraeq);
  }
  /* DBBD on the chain alone (no partition set): when the chain blocks
     hold under half the system the border solve is the whole problem,
     and the degenerate split crashed in the block extraction */
  if(solmethod!=SM_PROBE&&matsol==MM_DBBD&&allregset<0&&alltimeset>=0&&2*netcut>VecSize) {
    if(rank==0)errmsg("Error: DBBD (-matsol 2) found no partition set and its chain blocks over %s hold only %ld of %ld equations (border %.1f%%); use -matsol 1 (SBBD) or -matsol 0 (LU)\n",sets[alltimeset].setname,(long)(VecSize-netcut),(long)VecSize,VecSize>0?100.0*netcut/VecSize:0.0);
    PetscFinalize();
    return 1;
  }
  /* rank 0 always holds valid ordering data: rank_hsl==0 under HSL, and
     under nohsl every rank computes the full ordering */
  if(rank==0) {
    stats_run_options ropt;
    ropt.postsim_on=postsim_on;
    ropt.cofdump=cofdump;
    ropt.subints=(long)subints;
    ropt.adaptive=(int)adaptive;
    ropt.rk_chart=rko.chart;
    ropt.rk_norm=rko.norm;
    ropt.rk_ctrl=rko.ctrl;
    ropt.rk_scope=rko.scope;
    ropt.rk_h0=rko.h0;
    ropt.epstol=(double)epstol;
    ropt.maxretries=maxretries;
    ropt.retryadj=(double)retryadj;
    ropt.laA=(long)laA;
    ropt.laDi=(long)laDi;
    ropt.laD=(long)laD;
    ropt.comp_steps=comp_steps;
    ropt.comp_redo=comp_redo;
    ropt.comp_do_approx=comp_do_approx;
    ropt.comp_do_acc=comp_do_acc;
    ropt.comp_sberr_warn=comp_sberr_warn;
    ropt.comp_pass=comp_acc_phase;
    ropt.comp_minfrac=comp_minfrac;
    ordering_stats_write(iodata,niodata,noutdata,nsoldata,(long)VecSize,nvarele,nexo,matsol,solmed,nesteddbbd,(long)mpisize,mc66,alltimeset,allregset,sets,ntime,nreg,ndblock,netcut,nintraeq,countvarintra1,counteqnoadd,chain_source,partition_source,partition_auto_json,&ropt);
  }
  if(partition_auto_json!=NULL) {
    free(partition_auto_json);
    partition_auto_json=NULL;
  }
  free(eq_defs);
  free(eq_time_offsets);
  free(eq_reg_offsets);
  if(nohsl) {
    VecCreate(PETSC_COMM_WORLD,&vece);
  }
  else {
    VecCreate(PETSC_COMM_SELF,&vece);
  }
  if(nohsl) {
    VecSetType(vece,VECMPI);
  }
  else {
    VecSetType(vece,VECSEQ);
  }
  int localsize=0,nmatint,localbeg,localend;
  int *locals= (int *) calloc (mpisize,sizeof(int));
  if(nesteddbbd==1) {
    nmatint=ntime/mpisize;
    for(i=0; i<mpisize; i++)locals[i]=nmatint;
    for(i=0; i<mpisize; i++)if(i<ntime-mpisize*nmatint)locals[i]++;
    localbeg=0;
    for(i=0; i<mpisize; i++)if(i<rank)localbeg+=locals[i]*(nreg+1);
    localend=0;
    for(i=0; i<mpisize; i++)if(i<rank+1)localend+=locals[i]*(nreg+1);
    if(rank==mpisize-1)localend=ndblock;
    logmsg(2,"rank %d localbeg %d localend %d\n",rank,localbeg,localend);
    localsize=0;
    for (i=1; i<ndblock+1; i++)if(i>localbeg&&i<=localend)localsize+=counteq[i]-counteq[i-1];
    logmsg(2,"rank %d localsize %d\n",rank,localsize);
  }
  /* exogenous columns run 0..nexo-1; condensation can push nexo past
     VecSize, so the shock vector and B's columns span BSize */
  PetscInt BSize=(VecSize>(PetscInt)nexo)?VecSize:(PetscInt)nexo;
  shock_vec_set_sizes(vece,nesteddbbd,localsize,VecSize,BSize);
  VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
  free(locals);
  /* Preallocation is indexed by this rank's rows of A, while B's diagonal
     block is delimited by the exogenous column range vece carries.  The
     two coincide for a square B -- the layout vece used to define on its
     own -- but part company once nexo > VecSize, so A's row split is
     taken from PETSc directly (the same PetscLayout arithmetic MatSetUp
     will apply to A) instead of being read off the shock vector. */
  PetscInt Cstart,Cend;
  ierr = VecGetOwnershipRange(vece,&Cstart,&Cend);
  CHKERRQ(ierr);
  if(nohsl) {
    PetscInt nrowlocal=(nesteddbbd==1)?(PetscInt)localsize:PETSC_DECIDE,nrows=VecSize;
    ierr = PetscSplitOwnership(PETSC_COMM_WORLD,&nrowlocal,&nrows);
    CHKERRQ(ierr);
    MPI_Scan(&nrowlocal,&Iend,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);
    Istart=Iend-nrowlocal;
  }
  else {
    Istart=0;
    Iend=VecSize;
  }
  if(rank==rank_hsl) {
    ierr = PetscMalloc((Iend-Istart)*sizeof(PetscInt),&dnnz);
    CHKERRQ(ierr);
    ierr = PetscMalloc((Iend-Istart)*sizeof(PetscInt),&onnz);
    CHKERRQ(ierr);
    ierr = PetscMalloc((Iend-Istart)*sizeof(PetscInt),&dnnzB);
    CHKERRQ(ierr);
    ierr = PetscMalloc((Iend-Istart)*sizeof(PetscInt),&onnzB);
    CHKERRQ(ierr);
    for (i=Istart; i<Iend; i++) {
      dnnz[i-Istart]=1;
      onnz[i-Istart]=0;
      dnnzB[i-Istart]=1;
      onnzB[i-Istart]=0;
    }
  }
  logmsg(2,"rank11 %d Istart %d I end %d\n",rank, Istart,Iend);
  /* capture per-statement row-addressing metadata (side table, a few
     hundred KB) so the probe AND the on-failure diagnosis can name
     defective equation elements */
  eq_probe_meta *eqmeta=NULL;
  offset_t neqmeta=0;
  if(rank==rank_hsl)eqmeta= (eq_probe_meta *) calloc (neq+1,sizeof(eq_probe_meta));
  if(rank==rank_hsl) {
    jacobian_preallocate(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,nexo,closure_vals,ndblock,alltimeset,allregset,eq_intertemp,eq_addr,eq_time,eq_reg,counteq,nintraeq,&sbbd_overrid,VecSize,Istart,Iend,Cstart,Cend,&dnz,dnnz,&onz,onnz,&dnzB,dnnzB,&onzB,onnzB,nesteddbbd,eqmeta,&neqmeta);
    probe_onfail_context(sets,set_elems,vars,nvar,closure_vals,nvarele,eq_addr,eqmeta,neqmeta,VecSize);
  }
  /* -jacdump (Tier B5): the base point is now -- after the Reads and
     Formulas, before the first step. The complementarity rows take their
     pre-simulation state branch, as in the probe. Reported, not fatal. */
  if(jacdump&&comp_acc_phase==0&&rank==0) {
    elem_value *cvsave=NULL;
    if(teems_ncomp>0) {
      /* the state weights are written into coefficient slots; the solve
         must find the slots as the formulas left them */
      cvsave=(elem_value *) malloc ((ncofele>0?ncofele:1)*sizeof(elem_value));
      if(cvsave!=NULL) {
        memcpy(cvsave,elem_vals,ncofele*sizeof(elem_value));
        comp_states_set(sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals);
      }
      else printf("Warning: -jacdump: no memory to protect the complementarity state coefficients; the E_$comp rows are exported with their formula values\n");
    }
    jacobian_dump(teems_sol_stem,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,nvarele,closure_vals,eq_addr,VecSize,eqmeta,neqmeta);
    if(cvsave!=NULL) {
      comp_states_free();
      memcpy(elem_vals,cvsave,ncofele*sizeof(elem_value));
      free(cvsave);
    }
  }
  if(jacdump==2&&comp_acc_phase==0) {
    if(rank==0)logmsg(1,"Jacobian export only (-jacdump 2): the solve is skipped\n");
    jac_only=true;
    solmethod=SM_PROBE;
  }
  if(rank==0&&sbbd_overrid&&alltimeset<0&&(set_find_alltime(sets,nset)<0||sets[set_find_alltime(sets,nset)].size>1)) {
    printf("Warning: the equations reference intertemporal sets but this run's ordering ignores that structure; a bordered matrix method (-matsol 1/2/3) would detect and exploit it\n");
  }
  free(eq_intertemp);
  free(eq_time);
  free(eq_reg);
  dnz=0;
  dnzB=0;
  onz=0;
  onzB=0;
  for (i=Istart; i<Iend; i++) {
    if (dnnzB[i-Istart]+onnzB[i-Istart]>nexo-1&&dnnzB[i-Istart]>1) {
      dnnzB[i-Istart]--;
    }
    if (dnnz[i-Istart]+onnz[i-Istart]>nvarele-nexo&&dnnz[i-Istart]>1) {
      dnnz[i-Istart]--;
    }
    if (dnnz[i-Istart]>dnz) {
      dnz=dnnz[i-Istart];
    }
    if (onnz[i-Istart]>dnz) {
      onz=dnnz[i-Istart];
    }
    if (dnnzB[i-Istart]>dnzB) {
      dnzB=dnnzB[i-Istart];
    }
    if (onnzB[i-Istart]>onzB) {
      onzB=onnzB[i-Istart];
    }
  }
  /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
  MPI_Fint fcomm;
  fcomm = MPI_Comm_c2f(PETSC_COMM_WORLD);

  fortran_int indata[5];
  indata[1]=VecSize;
  indata[2]=mpisize;
  indata[3]=ndblock;
  indata[4]=nsbbdblocks;
  fortran_int *ptx=NULL;
  ptx = indata;

  solve_real *xcf=NULL;
  solve_real *accmetric=NULL; /* embedded-RK cumulative error metrics */
  solve_real *x0=NULL;// (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, solve_real *VA, solve_real *B, solve_real *X);
  extern void spec48_single_(fortran_int *indata,int *irn, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
  extern void spec48_nomc66_(fortran_int *indata, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow, MPI_Fint *fcomm,fortran_int *rowptrin, fortran_int *colptrin);

  logmsg(2,"rank %d ncof %ld\n",rank,ncof);
  

  if(inmemory) {
    /* Residency cost of skipping the driver spills: value arrays plus the
       Gragg step state. Fall back to scratch files unless it
       fits comfortably in available memory. */
    long need=(long)(ncofele+nvarele)*sizeof(elem_value)
             +(long)nvarele*(sizeof(closure_entry)+sizeof(unsigned char)+sizeof(store_real)+6*sizeof(solve_real)+sizeof(int));
    double av=teems_mem_avail_bytes();
    long avail=av>0?(long)av:-1; /* MemAvailable capped by the container's cgroup limit */
    if(avail>0&&2*need>avail) {
      if(rank==0)printf("Warning: -inmemory needs ~%ld MB per rank but only ~%ld MB is available; using scratch files instead\n",need/1048576,avail/1048576);
      inmemory=0;
    } else if(rank==0)logmsg(1,"inmemory: keeping ~%ld MB of value arrays resident per rank\n",need/1048576);
  }
  /* C2/C3 (manual 51.1.3): endogenous complementarity-variable
     components route the first pass to the approximate (Euler) run
     with the E_$comp state machinery live; the 51.7.1 closure/shock
     modification then re-enters the pipeline and the REQUESTED method
     solves the accurate run (comp_acc_phase 1). */
  if(has_sub&&solmethod!=SM_PROBE) {
    int r;
    MPI_Bcast(&teems_nsub,1,MPI_INT,0,PETSC_COMM_WORLD);
    if(rank!=0)teems_subs=(sub_group *) calloc (teems_nsub,sizeof(sub_group));
    for(r=0; r<teems_nsub; r++) {
      MPI_Bcast(&teems_subs[r].nmem,sizeof(offset_t),MPI_BYTE,0,PETSC_COMM_WORLD);
      if(rank!=0)teems_subs[r].mem=(offset_t *) malloc ((teems_subs[r].nmem>0?teems_subs[r].nmem:1)*sizeof(offset_t));
      MPI_Bcast(teems_subs[r].mem,teems_subs[r].nmem*sizeof(offset_t),MPI_BYTE,0,PETSC_COMM_WORLD);
      free(teems_subs[r].exo);
      teems_subs[r].exo=(exo_idx_t *) calloc (teems_subs[r].nmem>0?teems_subs[r].nmem:1,sizeof(exo_idx_t));
      if(rank==rank_hsl) {
        offset_t k;
        for(k=0; k<teems_subs[r].nmem; k++)teems_subs[r].exo[k]=closure_vals[teems_subs[r].mem[k]].exo_index;
      }
    }
    teems_sub_active=1;
  }
  bool comp_dispatch=(teems_comp_active>0&&solmethod!=SM_PROBE&&comp_acc_phase==0);
  if(comp_dispatch) {
    if(comp_flags0==NULL&&comp_sub==0&&subints>1&&comp_do_acc) {
      teems_comp_nsub=(int)subints;
      if(rank==rank_hsl) {
        comp_flags0=(unsigned char *) malloc (nvarele>0?nvarele:1);
        comp_shock0=(store_real *) malloc ((nvarele>0?nvarele:1)*sizeof(store_real));
        memcpy(comp_flags0,teems_cl_flags,nvarele);
        memcpy(comp_shock0,teems_cl_shock,nvarele*sizeof(store_real));
      }
    }
    if(teems_comp_nsub>1) {
      if(rank==0)printf("Complementarity: subinterval %d of %d -- an approximate and an accurate run (GEMPACK manual 51.7.4)\n",comp_sub+1,teems_comp_nsub);
      /* this subinterval's shocks relative to its start: the shares are
         linear in the pre-simulation levels, as for any subinterval run */
      if(rank==rank_hsl)for(i=0; i<nvar; i++) {
        offset_t e,x;
        for(e=0; e<vars[i].nelem; e++) {
          x=vars[i].offset+e;
          if(!(comp_flags0[x]&CL_F_EXO))continue;
          if(vars[i].change_real)CL_SHOCK(x)=comp_shock0[x];
          else CL_SHOCK(x)=(store_real)(100.0*((100.0+(comp_sub+1)*(double)comp_shock0[x])/(100.0+comp_sub*(double)comp_shock0[x])-1.0));
        }
      }
    }
    if(comp_do_approx) {
      if(rank==0)printf("Complementarity: %ld active (endogenous) component(s); approximate simulation as forward Euler with %d steps (GEMPACK manual 51.1.2)\n",(long)teems_comp_active,comp_steps);
      if(rank==0&&subints>1&&!comp_do_acc)printf("Warning: without an accurate run the complementarity approximate run treats the simulation as one interval (GEMPACK manual 51.7.4 pairs need the accurate run)\n");
      solve_comp_approx(nohsl,VecSize,dnz,dnnz,onz,onnz,dnzB,dnnzB,onzB,onnzB,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,teems_ma48_cntl4,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,fcomm,comp_steps,comp_redo,comp_minfrac,&xcf);
      if(rank==0&&comp_do_acc&&xcf!=NULL) {
        free(comp_approx_col);
        comp_approx_col=(solve_real *) malloc (nvarele*sizeof(solve_real));
        if(comp_approx_col!=NULL)memcpy(comp_approx_col,xcf,nvarele*sizeof(solve_real));
        if(comp_approx_col!=NULL&&teems_comp_nsub>1) {
          if(comp_apx_tot==NULL)comp_apx_tot=(solve_real *) calloc (nvarele>0?nvarele:1,sizeof(solve_real));
          for(i=0; i<nvar; i++) {
            offset_t e;
            for(e=vars[i].offset; e<vars[i].offset+vars[i].nelem; e++) {
              if(vars[i].change_real)comp_apx_tot[e]+=comp_approx_col[e];
              else comp_apx_tot[e]+=comp_approx_col[e]*(1+comp_apx_tot[e]/100);
            }
          }
        }
      }
    }
    else {
      if(rank==0)printf("Complementarity: -comp_do_approx 0; taking the pre-simulation states as the accurate run's targets (GEMPACK manual 51.6)\n");
      VecDestroy(&vece); /* the skipped approximate driver would have consumed it */
    }
    if(comp_do_acc) {
      /* capture the target states from the current values (post-
         approximate, or pre-sim when the approximate run was skipped) */
      if(rank==rank_hsl) {
        if(comp_accurate_prepare(sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals)<0)MPI_Abort(PETSC_COMM_WORLD,1);
      }
      /* restore the pre-simulation values: the accurate run solves the
         original data with the modified closure (drivers re-snapshot
         .initial on entry; formulas/reads state is exactly the .initial
         snapshot taken before the first approximate step) */
      if(comp_do_approx) {
        for(i=0; i<ncofele+nvarele; i++) {
          elem_vals[i].value=elem_vals[i].initial;
          elem_vals[i].substep_base=0;
        }
      }
      /* 51.7.1 closure/shock modification on the pre-sim values (the
         shocks are pre-sim-to-target differences); rank coverage
         mirrors the update passes -- rank 0 under HSL, every rank
         under nohsl -- so the numbering below reads the same closure
         it read on pass 1 */
      if(rank==rank_hsl) {
        if(comp_accurate_closure(closure_vals,vars,nvar,coefs,ncof,sets,nset,set_elems,elem_vals)<0)MPI_Abort(PETSC_COMM_WORLD,1);
        /* a percentage-change complementarity variable whose bound is 0
           gets a -100 percent shock in the accurate run (51.7.1); under
           Gragg its update left the reals with no named cause, and the
           Runge-Kutta drivers halve their step towards it without end
           (the log chart then ends wrong with no message) */
        if(solmethod==SM_GRAGG||isrk)for(i=0; i<nvar; i++) {
          offset_t e;
          if(vars[i].change_real)continue;
          for(e=0; e<vars[i].nelem; e++)if(CL_EXO(vars[i].offset+e)&&CL_SHOCK(vars[i].offset+e)<=-100+1e-4) {
            char lab[4*NAMESIZE];
            array_element_label(&vars[i],e,lab,sizeof(lab));
            errmsg("Error: the complementarity accurate run (GEMPACK manual 51.7.1) shocks %s by -100 percent to reach its bound of zero, which the %s method cannot take (GEMPACK manual 30.2); declare %s (change,levels), or use the midpoint or Euler method\n",lab,solmed,vars[i].cofname);
            MPI_Abort(PETSC_COMM_WORLD,1);
          }
        }
      }
      if(rank==0)printf("Complementarity: accurate simulation with the %s method (closure/shocks modified per GEMPACK manual 51.7.1)\n",solmed);
      teems_comp_xac_base=comp_tot;
      comp_next_phase=1;
      /* tear down pass-1 state and re-enter the closure-dependent
         pipeline: exo_index numbering restarts from zero (backsolved
         ordinals are assigned by backsolve_read and keep theirs); the
         next subinterval's approximate run re-enters here too */
comp_teardown:
      for(i=0; i<nvarele; i++)if(!CL_BS(i))closure_vals[i].exo_index=0;
      jacobian_cache_free();
      backsolve_cache_free();
      free(xcf);
      xcf=NULL;
      free(eqmeta);
      free(eq_addr);
      free(counteq);
      free(counteqnoadd);
      free(countvarintra1);
      if(nesteddbbd==1) {
        free(ndbbddrank1);
        ndbbddrank1=NULL;
        ndbbd_cut_cache_free();
      }
      if(rank==rank_hsl) {
        PetscFree(dnnz);
        PetscFree(onnz);
        PetscFree(dnnzB);
        PetscFree(onnzB);
      }
      var_inter=(bool *) calloc (nvar,sizeof(bool));
      ele_inter=(bool *) calloc (nvarele,sizeof(bool));
      comp_acc_phase=comp_next_phase;
      goto comp_accurate_reentry;
    }
    else {
      if(rank==0)printf("Complementarity: -comp_do_acc 0; the approximate run's solution is the simulation result (GEMPACK manual 51.6)\n");
      comp_states_free();
    }
  }

  if(!comp_dispatch&&solmethod==SM_JOHANSEN)solve_johansen(nohsl,VecSize,A,dnz,dnnz,onz,onnz,B,dnzB,dnnzB,onzB,onnzB,vecb,vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele+nvarele,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,teems_ma48_cntl4,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,&xcf);

  FILE* solution;

  if(!comp_dispatch&&(solmethod==SM_GRAGG||solmethod==SM_MIDPOINT||solmethod==SM_EULER))solve_gragg(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele+nvarele,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,teems_ma48_cntl4,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,(teems_comp_nsub>1)?1:subints,fcomm,solmethod,&xcf);

  if(!comp_dispatch&&isrk)solve_rk(nohsl,VecSize,dnz,dnnz,onz,onnz,dnzB,dnnzB,onzB,onnzB,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,teems_ma48_cntl4,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,fcomm,solmethod,adaptive,(double)epstol,(double)retryadj,maxretries,&rko,&xcf,&accmetric);

  /* C3: after the accurate run, every component must sit in its
     approximate-run state with the variable inside its bounds
     (manual 51.5.4/51.7.5); fatal unless the CMF downgrades with
     -comp_sberr_warn 1 */
  if(comp_acc_phase==1) {
    offset_t comp_nbad=0;
    if(rank==rank_hsl)comp_nbad=comp_verify_states(sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,comp_sberr_warn);
    if(rank==rank_hsl&&comp_nbad>0) {
      if(comp_sberr_warn)printf("Warning: %ld complementarity state/bound error(s) after the accurate run (treated as warnings per -comp_sberr_warn; check the log carefully, GEMPACK manual 51.6)\n",(long)comp_nbad);
      else {
        errmsg("Error: %ld complementarity state/bound error(s) after the accurate run; rerun with more Euler steps (-comp_steps) or smaller shocks, or downgrade with -comp_sberr_warn 1 (GEMPACK manual 51.5.4/51.6)\n",(long)comp_nbad);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
    }
    comp_states_free();
    if(teems_comp_nsub>1) {
      if(xcf!=NULL) {
        if(comp_tot==NULL)comp_tot=(solve_real *) calloc (nvarele>0?nvarele:1,sizeof(solve_real));
        for(i=0; i<nvar; i++) {
          offset_t e;
          for(e=vars[i].offset; e<vars[i].offset+vars[i].nelem; e++) {
            if(vars[i].change_real)comp_tot[e]+=xcf[e];
            else comp_tot[e]+=xcf[e]*(1+comp_tot[e]/100);
          }
        }
      }
      /* the embedded Runge-Kutta estimate accumulates over steps, and so
         over subintervals */
      if(accmetric!=NULL) {
        if(comp_est_tot==NULL)comp_est_tot=(solve_real *) calloc (nvarele>0?nvarele:1,sizeof(solve_real));
        for(i=0; i<nvarele; i++)comp_est_tot[i]+=accmetric[i];
      }
      teems_comp_xac_base=NULL;
      comp_sub++;
      if(comp_sub<teems_comp_nsub) {
        /* back to the user's closure for the next approximate run; its
           shocks are set from the shares at the dispatch */
        if(rank==rank_hsl) {
          memcpy(teems_cl_flags,comp_flags0,nvarele);
          memcpy(teems_cl_shock,comp_shock0,nvarele*sizeof(store_real));
        }
        free(accmetric); /* re-declared NULL after the re-entry */
        accmetric=NULL;
        comp_next_phase=0;
        goto comp_teardown;
      }
      if(xcf!=NULL&&comp_tot!=NULL)memcpy(xcf,comp_tot,nvarele*sizeof(solve_real));
      if(comp_approx_col!=NULL&&comp_apx_tot!=NULL)memcpy(comp_approx_col,comp_apx_tot,nvarele*sizeof(solve_real));
      if(accmetric!=NULL&&comp_est_tot!=NULL)memcpy(accmetric,comp_est_tot,nvarele*sizeof(solve_real));
      free(comp_tot); free(comp_apx_tot); free(comp_est_tot); free(comp_flags0); free(comp_shock0);
      comp_tot=comp_apx_tot=comp_est_tot=NULL; comp_flags0=NULL; comp_shock0=NULL;
    }
  }

  /* -solmed probe: MC79 structural diagnosis rides the probe after the
     (skipped) solve dispatch — assemble the Jacobian and run the
     matching / Dulmage-Mendelsohn analysis with named defects */
  if(solmethod==SM_PROBE) {
    /* realize the pre-simulation complementarity states so the
       E_$comp rows probe with their genuine (state-branch) pattern
       rather than the all-zero weights */
    if(!jac_only&&rank==rank_hsl&&teems_ncomp>0) {
      comp_states_set(sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals);
      comp_states_free();
    }
    if(!jac_only)probe_structural(VecSize,nvarele,ncofele,dnz,dnnz,dnzB,dnnzB,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,closure_vals,ndblock,alltimeset,allregset,eq_addr,counteq,nintraeq,eqmeta,neqmeta,iodata,niodata,noutdata,nsoldata,probefine,mpisize,rank);
    VecDestroy(&vece); /* the skipped solve driver would have destroyed it */
    if(rank==0) {
      gettimeofday(&endtime, NULL);
      logmsg(1,"probe: done in %.2f s\n",(endtime.tv_sec - probe_begin.tv_sec)+((double)(endtime.tv_usec - probe_begin.tv_usec))/ 1000000);
    }
  }
  free(eqmeta);

  jacobian_cache_free();
  backsolve_cache_free();
  lu_fastrefac_free();
  sbbd_fastrefac_free();
  dbbd_fastrefac_free();
  dbbd_fastextract_free();
  ndbbd_fastrefac_free();
  ndbbd_cut_cache_free();
  free(backsolves);


  if(rank==rank_hsl) {
    for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
      if (strcmp("solfiles",iodata[i].logname)==0) {
        strcpy(tempchar,iodata[i].filname);
        break;
      }
    }
    if(i==niodata+noutdata+nsoldata) {
      strcpy(tempchar,"solution");
    }
    /* probe runs have no solution: leave a previous run's sol.bin (and
       its .est) alone rather than truncating them */
    if(xcf!=NULL) {
      strcpy(solchar,tempchar);
      strcat(solchar,".bin");
      logmsg(2,"solchar %s\n",solchar);
      if ( (solution = fopen(solchar, "wb")) == NULL ) {
        errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",solchar,strerror(errno),(int)getuid());
        return 1;
      }
      fwrite(xcf, sizeof(solve_real),nvarele, solution);
      fclose(solution);
      outputs_note(solchar,"solution",0);
    }
    if(xcf!=NULL&&accmetric==NULL) {
      /* a non-embedded run after an embedded one in the same directory:
         drop the stale estimate so it cannot be read against this
         run's solution (teems-R attaches <sol>.est when present) */
      strcpy(solchar,tempchar);
      strcat(solchar,".est");
      remove(solchar);
      strcpy(solchar,tempchar);
      strcat(solchar,".acc");
      remove(solchar);
    }
    if(accmetric!=NULL) {
      /* embedded-RK estimated error metrics, one double per variable
         element in .bin order (the accumulated per-step embedded
         estimate -- an indicator, not a bound; GEMPACK writes the same
         quantity as <sol>.acc) */
      strcpy(solchar,tempchar);
      strcat(solchar,".est");
      logmsg(2,"solchar %s\n",solchar);
      if ( (solution = fopen(solchar, "wb")) == NULL ) {
        errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",solchar,strerror(errno),(int)getuid());
        return 1;
      }
      fwrite(accmetric, sizeof(solve_real),nvarele, solution);
      fclose(solution);
      outputs_note(solchar,"rk_error_estimate",0);
    }
    if(solmethod!=SM_PROBE&&structure_files_write(tempchar,vars,nvar,sets,nset,set_elems,nsetspace,nvarele))return 1;
    if(rank==0&&comp_approx_col!=NULL&&xcf!=NULL) {
      char label[96];
      const solve_real *cp[1]={comp_approx_col};
      const char *lp[1]={label};
      snprintf(label,sizeof(label),"complementarity approximate run (Euler, %d steps)",comp_steps);
      cols_write(teems_sol_stem,3,1,nvarele,cp,lp,NULL);
    }
    if(rank==0&&teems_sub_active&&sub_columns()!=NULL&&xcf!=NULL) {
      const char **lp=(const char **) malloc (teems_nsub*sizeof(char*));
      const char **sp=(const char **) malloc (teems_nsub*sizeof(char*));
      int r;
      for(r=0; r<teems_nsub; r++) {
        lp[r]=teems_subs[r].label;
        sp[r]=teems_subs[r].spec;
      }
      cols_write(teems_sol_stem,solmethod==SM_JOHANSEN?2:1,teems_nsub,nvarele,(const solve_real *const *)sub_columns(),lp,sp);
      free(lp);
      free(sp);
    }
    sub_free();
    if(teems_subs!=NULL) {
      int r;
      for(r=0; r<teems_nsub; r++) {
        free(teems_subs[r].mem);
        free(teems_subs[r].exo);
      }
      free(teems_subs);
      teems_subs=NULL;
      teems_nsub=0;
    }
  }
  free(comp_approx_col);
  MPI_Barrier(PETSC_COMM_WORLD);
  /* la* auto-sizing record: reduce the per-rank grown maxima and patch
     the effective percents into stats.json */
  {
    long laused[3]={teems_laA_used,teems_laDi_used,teems_laD_used},lamax[3];
    MPI_Allreduce(laused,lamax,3,MPI_LONG,MPI_MAX,PETSC_COMM_WORLD);
    if(rank==0)stats_la_used_patch(iodata,niodata,noutdata,nsoldata,
                                   lamax[0]>(long)laA?lamax[0]:(long)laA,
                                   lamax[1]>(long)laDi?lamax[1]:(long)laDi,
                                   lamax[2]>(long)laD?lamax[2]:(long)laD);
  }
  /* -condest record: the measures accumulate on the solving rank
     (rank_hsl); reduce and patch alongside la_used */
  if(teems_condest) {
    double cdd[3]={teems_condest_kw1max,teems_condest_kw2max,teems_condest_omegamax},cddmax[3];
    long cdl[2]={teems_condest_solves,teems_condest_skips},cdlsum[2];
    MPI_Allreduce(cdd,cddmax,3,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
    MPI_Allreduce(cdl,cdlsum,2,MPI_LONG,MPI_SUM,PETSC_COMM_WORLD);
    if(rank==0)stats_condest_patch(iodata,niodata,noutdata,nsoldata,
                                   cddmax[0],cddmax[1],cddmax[2],cdlsum[0],cdlsum[1]);
  }
  {
    double rmax=0.0;
    long rl[3]={teems_resid_solves,teems_resid_warn,teems_resid_skipped},rls[3];
    MPI_Allreduce(&teems_resid_max,&rmax,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
    MPI_Allreduce(rl,rls,3,MPI_LONG,MPI_SUM,PETSC_COMM_WORLD);
    if(rank==0&&(rls[0]>0||rls[2]>0)) {
      if(rls[0]>0)logmsg(1,"Maximum residual ratio across the whole simulation is %.8e (%ld solves checked; GEMPACK manual 30.1.5)\n",rmax,rls[0]);
      if(rls[2]>0)logmsg(1,"Residual ratios not checked for %ld solve(s): this matrix method releases the LHS matrix before the solve completes\n",rls[2]);
      stats_resid_patch(iodata,niodata,noutdata,nsoldata,rmax,rls[0],rls[1],rls[2]);
      if(rls[1]>0)printf("Warning: there have been %ld warnings about equations not being satisfied very accurately; the more there are, the more likely it is that the solution is not valid (GEMPACK manual 30.6.1)\n",rls[1]);
      /* GEMPACK finishes and saves every file, then ends with an error
         (30.6.1); errmsg sets the exit status without aborting */
      if(rls[1]>100)errmsg("Error: more than 100 equations were not satisfied very accurately (%ld warnings); all files are written, but the solution may not be valid (GEMPACK manual 30.6.1)\n",rls[1]);
    }
  }
  {
    long rfs=0;
    double rfb=0.0,rfa=0.0,rfsec=0.0;
    fh_refine_totals(&rfs,&rfb,&rfa,&rfsec);
    if(rank==0&&matsol==MM_DBBD) {
      if(rfs>0) {
        logmsg(1,"Refinement (DBBD, one step per solve): %ld solves refined; worst residual ratio %.3e before the step, %.3e after; %.2f s\n",rfs,rfb,rfa,rfsec);
        stats_refine_patch(iodata,niodata,noutdata,nsoldata,rfs,rfb,rfa,rfsec);
      }
      else if(!teems_refine)logmsg(1,"Refinement off (-refine 0): DBBD solves are not refined\n");
    }
  }
  fh_selftest_summary(rank);
  if(isrk&&!comp_dispatch&&rank==0)stats_rk_patch(iodata,niodata,noutdata,nsoldata);
  /* phase resident-memory record: a last probe for the run's high-water
     mark (collective), then patch the per-phase table */
  teems_rss_probe(solmethod==SM_PROBE?"structural diagnosis":"solution write");
  if(rank==0)stats_rss_patch(iodata,niodata,noutdata,nsoldata);
  if(rank==0)stats_ndbbd_threads_patch(iodata,niodata,noutdata,nsoldata);
  /* PostSim foundation F3 (early Tier 0): after the solve, coefficient
     slots hold post-simulation (updated) values and xcf holds the
     composed solution; expose the solution to the formula engine and
     evaluate (postsim) assertions. Zero cost when the TAB has none. */
  if(rank==0&&postsim_on&&xcf!=NULL&&elem_vals!=NULL&&(npostsim>0||tab_has_postsim_assertions(tabfile))) {
    postsim_expose_results(elem_vals,ncofele,nvarele,xcf);
    if(npostsim>0) {
      logmsg(1,"postsim: running %d statement(s)\n",npostsim);
      teems_ps_pass=1;
      /* Reads, Formulas and Assertions in file order (manual 12.2.1);
         the Read targets are checked first (12.2.3) */
      if(postsim_reads_execute(psfile,niodata,iodata,sets,nset,set_elems,coefs,ncof,ncofele,vars,nvar,nvarele,elem_vals)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
      ord_plan_build(psfile,1,coefs,ncof,vars,nvar);
      strcpy(commsyntax,"formula");
      statements_execute(psfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,true,0);
      teems_ps_pass=0;
      teems_ps_ran=1;
    }
    assertions_execute(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,true,teems_assertions_mode,1);
  }
  {
    int cofdumped=0,wr;
    long nsets_w=0,nother_w=0,nskip_w=0;
    if(cofdump&&rank==0&&elem_vals!=NULL) {
      for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++)
        if (strcmp("solfiles",iodata[i].logname)==0)break;
      if(coefficients_dump(i<niodata+noutdata+nsoldata?iodata[i].filname:"solution",coefs,ncof,ncofele,elem_vals))MPI_Abort(PETSC_COMM_WORLD,1);
      cofdumped=1;
    }
    if(nowrites==0&&rank==0)for(i=0; i<noutdata; i++){
      if(postsim_write_skipped(iodata[i+niodata].logname)) { nskip_w++; continue; }
      wr=outputs_write_csv(tabfile,iodata[i+niodata].logname,iodata[i+niodata].filname,sets,nset,set_elems,coefs,ncof,ncofele,vars,nvar,nvarele,elem_vals);
      if(wr!=-1) {
        if(wr==1)nsets_w++; else nother_w++;
        outputs_note(iodata[i+niodata].filname,"csv",0);
      }
    }
    if(rank==0)outputs_summary_log(cofdumped,coefs,ncof,ncofele,nsets_w,nother_w,nskip_w,nowrites==0,xcf!=NULL,nvar,nvarele);
  }
  /* completion marker, the run's last write: only a run that printed no
     Error line on any rank gets one */
  {
    int nerr=teems_error_count,nerrsum=0;
    MPI_Allreduce(&nerr,&nerrsum,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
    if(rank==0&&nerrsum==0)outputs_json_write(teems_sol_stem,run_id);
  }
  free(iodata);
  free(countvarintra1);
  free(eq_addr);
  free(ndbbddrank1);
  free(counteq);
  free(counteqnoadd);
  free(sets);
  free(set_elems);
  free(coefs);
  free(vars);
  free(closure_vals);
  free(elem_vals);
  free(xcf);
  free(accmetric);
  if(x0!=NULL)free(x0);
//**************************************************************************************
//**************************************END HSL*****************************************
//**************************************************************************************

  ierr = PetscFree(dnnzB);
  CHKERRQ(ierr);
  ierr = PetscFree(onnzB);
  CHKERRQ(ierr);
  ierr = PetscFree(dnnz);
  CHKERRQ(ierr);
  ierr = PetscFree(onnz);
  CHKERRQ(ierr);
  MPI_Comm_free(&node_comm);
  MPI_Comm_free(&node_tail_comm);
  ierr = PetscFinalize();
  CHKERRQ(ierr);
  /* exit-status contract: 0 only when no Error line was printed by
     this rank (mpiexec propagates any rank's nonzero status) */
  return teems_error_count>0?1:0;
}


