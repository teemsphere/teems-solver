// Pham Van Ha & Tom Kompas
// Ensure that sum over time involves only intertemporal variables
// If otherwise switch back to normal MC66
// Currently no sub-reg set, can try different varreg to get opt: comment break in nestedmatvarread
// Currently logic fdim <> only in formula
#include <teems_solver.h>

static char help[] = "Solves a CGE model in parallel with KSP.\n\
           Input parameters include:\n\
           -None at the moment\n\n";

/* Per-run ordering statistics (<solfiles>.stats.json): netcut, border
   sizes and per-block variable/equation counts. Written before the
   solve so failed runs still record their ordering; consumed by
   teems-R for matrix_method calibration. Fields are null when no
   bordered ordering was built (no chain/partition dimension).
   chain_source/partition_source say how each dimension was determined
   ("structural" detection or "none"); auto_json
   carries the partition-candidate table (NULL when nothing was
   probed). */
static void ordering_stats_write(cmf_file_entry *iodata, int niodata, int noutdata, int nsoldata,
                                 long VecSize, offset_t nvarele, offset_t nexo,
                                 dim_t matsol, char *solmed, dim_t nesteddbbd, long mpisize, dim_t mc66,
                                 offset_t alltimeset, offset_t allregset, set_def *sets,
                                 offset_t ntime, offset_t nreg, offset_t ndblock,
                                 offset_t netcut, offset_t nintraeq,
                                 offset_t *countvarintra1, offset_t *counteqnoadd,
                                 const char *chain_source, const char *partition_source, const char *auto_json) {
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
  bool bordered=alltimeset>=0||allregset>=0;
  fprintf(fp,"{\n");
  fprintf(fp,"  \"version\": 2,\n");
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
  fprintf(fp,"]\n}\n");
  fclose(fp);
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
static void chain_flags_apply(set_def *sets, dim_t nset) {
  offset_t i,j;
  for(i=0; i<nset; i++) {
    for(j=1; j<MAXSUPSET; j++) {
      if(sets[i].subsetid[j]>-1) {
        if(sets[sets[i].subsetid[j]].intertemp) {
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
        if(!closure_vals[j3+j].is_exogenous&&!closure_vals[j3+j].is_backsolved) {
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
        if(!closure_vals[j3+j].is_exogenous&&!closure_vals[j3+j].is_backsolved) {
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
        if(!closure_vals[j3+j].is_exogenous&&!closure_vals[j3+j].is_backsolved) {
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
  filehandle=fopen(fname,"r");
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
    while (str_replace_all(readitem," ", ""));
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
    if(chain>=0)logmsg(1,"Chain dimension detected structurally: set %s (size %d), %ld lead/lag references\n",sets[chain].setname,sets[chain].size,topcount[chain]);
  }
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
  if(rank==0)logmsg(1,"Partition detection: probing %ld candidate sets against the equation structure\n",ncand);
  for(s=0; s<nset; s++) {
    if(!iscand[s])continue;
    partition_probe(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,nvarele,closure_vals,neq,VecSize,alltimeset,ntime,nesteddbbd,s,&cnetcut[s],&cnblocks[s],&cmin[s],&cmax[s]);
    viable[s]=cnblocks[s]>=mpisize&&2*cnetcut[s]<VecSize&&cmin[s]>0;
    if(rank==0)logmsg(1,"  %-14s blocks %6ld  border %8ld (%5.1f%%)  block min/max %ld/%ld  %s\n",
                        sets[s].setname,cnblocks[s],cnetcut[s],VecSize>0?100.0*cnetcut[s]/VecSize:0.0,
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
    if(chosen>=0)logmsg(1,"Partition detection: selected set %s (%ld blocks, border %.1f%% of system)\n",
                          sets[chosen].setname,cnblocks[chosen],VecSize>0?100.0*cnetcut[chosen]/VecSize:0.0);
    else if(ncand>0)logmsg(1,"Partition detection: no viable candidate (need >=%ld nonempty blocks and a border below half the system)\n",mpisize);
    else logmsg(1,"Partition detection: no candidate sets to probe\n");
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
  Vec      vecb,vece,x;  /* approx solution, RHS, exact solution */
  Mat      A,B;    /* linear system matrix */
  PetscInt  rank=0,mpisize,rank_hsl=0;
  PetscInt     VecSize=0,Istart=0,Iend=0,dnz=0,onz=0,dnzB=0,onzB=0,count,*onnz,*dnnz,*onnzB,*dnnzB,its;
  PetscErrorCode ierr;
  PetscBool   flg;
  PetscScalar  value,zero=0;
  PetscLogDouble time0,time1;
  clock_t timestr,timeend,timemulti;
  struct timeval begintime, endtime,gettime_now;
  size_t freadresult;
  offset_t lasize=1;
  offset_t i,j;
  offset_t j2=0,j1=0,j0=0,j3,j4,j5,j6;
  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize);

  char processor_name[MPI_MAX_PROCESSOR_NAME+1];
  int name_len,name_len_max,name_beg,class_size,color,group_size,node_rank;
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
  bool sbbd_overrid=false;
  PetscBool nohsl=false;
  char ch='y';
  if(rank==0) {
    logmsg(2,"Notes:\n  Shock statement values follow GEMPACK ordering (first subscript varies fastest).\n  Declare intertemporal variables with minimum dimension to minimise the net cut,\n  e.g. capital(REG,TIME)=qo(\"capital\",REG,TIME) rather than shocking qo(COM,REG,TIME).\n  laA/laDi control solver workspace sizes; use the smallest that solves.\n  Beware CRLF line endings in model text files.\n");
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  //**************************************************************************************
  //****************************** READ SET ELEMENT***************************************
  //**************************************************************************************
  char tabfile[TABREADLINE],newtabfile[TABREADLINE]="_temp_tab_file",newtabfile1[TABREADLINE]="_temp_tab_new_file",closure[TABREADLINE],shock[TABREADLINE],filename[TABREADLINE],longname[TABREADLINE],vname[NAMESIZE],copyline[TABREADLINE];
  char psfile[TABREADLINE]="\0";
  int npostsim=0,postsim_on=1;
  char tempfilenam[255],tempchar[255],solmed[NAMESIZE],solchar[255];
  int niodata=0,nj,mem_fac=0,noutdata=0,nsoldata=0,nowrites=0;
  offset_t nsetspace=0,dcount,ndblock=0,netcut=0,ndblock1,nreg=0,ntime=0;
  dim_t nset=0,vsize,dim1,nlength=0,matsol=0,laA=2,laDi=2,laD=2,nsbbdblocks=2,nesteddbbd=0,mc66=0,subints=1,subindx;
  offset_t alltimeset=-1,allregset=-1;
  map_def *maps=NULL;
  dim_t nmap=0;
  PetscReal cntl6=0,cntl3;
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
  PetscOptionsGetInt(NULL,NULL,"-laA",&laA,NULL);
  if(laA==0)laA=2;
  PetscOptionsGetInt(NULL,NULL,"-laD",&laD,NULL);
  if(laD==0)laD=2;
  PetscOptionsGetInt(NULL,NULL,"-laDi",&laDi,NULL);
  if(laDi==0)laDi=2;
  PetscOptionsGetInt(NULL,NULL,"-withmc66",&mc66,NULL);
  /* -gpzerodivide 1: GEMPACK dual-class ZERODIVIDE semantics in formulas
     (manual 10.11/10.11.1; plan A1). Default 0 keeps the legacy single
     conflated default -- adoption is a re-anchor-class change. */
  PetscOptionsGetInt(NULL,NULL,"-gpzerodivide",&teems_gpzerodivide,NULL);
  PetscOptionsGetInt(NULL,NULL,"-step1",&steps1,NULL);
  if(steps1==0)steps1=2;
  PetscOptionsGetInt(NULL,NULL,"-step2",&steps2,NULL);
  if(steps2==0)steps2=4;
  PetscOptionsGetInt(NULL,NULL,"-step3",&steps3,NULL);
  if(steps3==0)steps3=8;
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
  PetscOptionsGetReal(NULL,NULL,"-cntl_6",&cntl6,NULL); /* CNTL6 in Mat Order */
  PetscOptionsGetReal(NULL,NULL,"-cntl_3",&cntl3,NULL);/*Iterative threshold */
  {
    PetscInt verb=verbosity;
    PetscOptionsGetInt(NULL,NULL,"-verbosity",&verb,NULL);/* 0 errors/results, 1 progress (default), 2 debug */
    verbosity=(int)verb;
    char verbstr[8];
    sprintf(verbstr,"%d",verbosity);
    setenv("TEEMS_VERBOSITY",verbstr,1);/* for the Fortran kernels (hsl_kernels.f90) */
  }
  inmemory=-1;
  PetscOptionsGetInt(NULL,NULL,"-inmemory",&inmemory,NULL);/* keep value arrays resident instead of spilling to scratch */
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
        if(rank==0)printf("Error: -tempdir must be 1-200 characters: %s\n",tmpopt);
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
      if(rank==0)printf("Error: -tempdir is not a writable directory: %s\n",scratch_dir);
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
    if(rank==0)printf("Warning: -solmed Mmid is deprecated; the method is Gragg's (smoothed modified midpoint) — use -solmed Gragg\n");
    strcpy(solmed,"Gragg");
  }
  if(strcmp(solmed,"NoSol")==0) {/* transitional alias: renamed — it is a structure probe, not a degenerate solve */
    if(rank==0)printf("Warning: -solmed NoSol is deprecated — use -solmed probe\n");
    strcpy(solmed,"probe");
  }
  int solmethod=0;
  if(strcmp(solmed,"Gragg")==0)solmethod=SM_GRAGG;
  if(strcmp(solmed,"Euler")==0)solmethod=SM_EULER;
  if(strcmp(solmed,"RK2")==0)solmethod=SM_RK2;
  if(strcmp(solmed,"RK4")==0)solmethod=SM_RK4;
  if(strcmp(solmed,"BoSha32")==0)solmethod=SM_BOSHA32;
  if(strcmp(solmed,"DoPri54")==0)solmethod=SM_DOPRI54;
  if(strcmp(solmed,"Johansen")==0)solmethod=SM_JOHANSEN;
  if(strcmp(solmed,"probe")==0)solmethod=SM_PROBE;
  if(solmethod==0) {
    if(rank==0)printf("Error: unknown -solmed %s (valid: Gragg, Euler, RK2, RK4, BoSha32, DoPri54, Johansen, probe)\n",solmed);
    PetscFinalize();
    return 1;
  }
  logmsg(2,"Sol med %d\n",solmethod);
  /* -probefine: with -solmed probe, add the MC79 fine Dulmage-Mendelsohn
     report (strongly connected components of the well-determined block) */
  dim_t probefine=0;
  PetscOptionsGetInt(NULL,NULL,"-probefine",&probefine,NULL);
  bool isrk=(solmethod==SM_RK2||solmethod==SM_RK4||solmethod==SM_BOSHA32||solmethod==SM_DOPRI54);
  bool isrk_embedded=(solmethod==SM_BOSHA32||solmethod==SM_DOPRI54);
  /* Runge-Kutta controls (GEMPACK 26.5.1/26.5.2): -adaptive
     no|yes|accuracy-only (embedded flavors only; accuracy-only skips
     the check-failure retries), -epstol the per-step error-metric
     bound, -retryadj/-maxretries the check-failure retry policy */
  int adaptive=0,maxretries=3;
  PetscReal epstol=0.1,retryadj=0.5;
  {
    char adaptbuf[NAMESIZE];
    PetscOptionsGetString(NULL,NULL,"-adaptive",adaptbuf,NAMESIZE,&flg);
    if(flg) {
      if(strcmp(adaptbuf,"no")==0)adaptive=0;
      else if(strcmp(adaptbuf,"yes")==0)adaptive=1;
      else if(strcmp(adaptbuf,"accuracy-only")==0)adaptive=2;
      else {
        if(rank==0)printf("Error: unknown -adaptive %s (valid: no, yes, accuracy-only)\n",adaptbuf);
        PetscFinalize();
        return 1;
      }
      if(adaptive&&!isrk_embedded) {
        if(rank==0)printf("Error: -adaptive requires an embedded Runge-Kutta method (-solmed BoSha32 or DoPri54)\n");
        PetscFinalize();
        return 1;
      }
    }
    PetscOptionsGetReal(NULL,NULL,"-epstol",&epstol,NULL);
    PetscOptionsGetReal(NULL,NULL,"-retryadj",&retryadj,NULL);
    PetscOptionsGetInt(NULL,NULL,"-maxretries",&maxretries,NULL);
    if(isrk) {
      if(steps1<1) {
        if(rank==0)printf("Error: -step1 must be at least 1 for Runge-Kutta methods (got %d)\n",steps1);
        PetscFinalize();
        return 1;
      }
      if(subints>1) {
        if(rank==0)printf("Error: subintervals are not available with Runge-Kutta methods (got -nsubints %d); increase -step1 instead\n",subints);
        PetscFinalize();
        return 1;
      }
      if(epstol<=0||retryadj<=0||retryadj>=1||maxretries<1) {
        if(rank==0)printf("Error: -epstol and -retryadj must be positive (-retryadj below 1) and -maxretries at least 1\n");
        PetscFinalize();
        return 1;
      }
      if(adaptive&&epstol<0.005) {
        if(rank==0)printf("Warning: -epstol %g is below 0.005; tolerances this tight are hard to achieve numerically and may reject many steps\n",(double)epstol);
      }
    }
  }
  if(solmethod==SM_GRAGG) {
    /* Gragg's h^2 error expansion (Pearson 1991 Thm 6.1) holds for even
       step counts only; mixed parity also breaks the shared-power
       extrapolation */
    step_ratio2=steps1/(double)2;
    i=(offset_t)steps1/2;
    if(step_ratio2!=i){
      //odd
      step_ratio2=steps2/(double)2;
      i=(offset_t)steps2/2;
      if(step_ratio2==i){
        printf("Error: -step1/-step2/-step3 must be all odd or all even (got %d %d %d)\n",steps1,steps2,steps3);
        return 0;
      }
      step_ratio2=steps3/(double)2;
      i=(offset_t)steps3/2;
      if(step_ratio2==i){
        printf("Error: -step1/-step2/-step3 must be all odd or all even (got %d %d %d)\n",steps1,steps2,steps3);
        return 0;
      }
    }else{
      //even
      step_ratio2=steps2/(double)2;
      i=(offset_t)steps2/2;
      if(step_ratio2!=i){
        printf("Error: -step1/-step2/-step3 must be all odd or all even (got %d %d %d)\n",steps1,steps2,steps3);
        return 0;
      }
      step_ratio2=steps3/(double)2;
      i=(offset_t)steps3/2;
      if(step_ratio2!=i){
        printf("Error: -step1/-step2/-step3 must be all odd or all even (got %d %d %d)\n",steps1,steps2,steps3);
        return 0;
      }
    }
  }
  if((solmethod==SM_GRAGG||solmethod==SM_EULER)&&!(steps1<steps2&&steps2<steps3)) {
    /* extrapolation needs three distinct step sizes (GEMPACK: i<j<k) */
    printf("Error: -step1/-step2/-step3 must be strictly increasing (got %d %d %d)\n",steps1,steps2,steps3);
    return 0;
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
  if(rank==rank_hsl) {
    mem_fac=1;
  }
  char *readitem=NULL;
  if(rank==0) {
    niodata=cmf_count_files(filename,"iodata");
    if(niodata==-1)return 0;
    noutdata=cmf_count_files(filename,"outdata");
    nsoldata=cmf_count_files(filename,"soldata");
  }
  if(nohsl) {
    MPI_Bcast(&niodata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&noutdata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nsoldata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  cmf_file_entry *iodata= (cmf_file_entry *) calloc (niodata+noutdata+nsoldata,sizeof(cmf_file_entry));
  if(rank==rank_hsl) {
    cmf_read(filename,niodata,iodata,tabfile,closure,shock);
    teems_assertions_mode=cmf_assertions_mode(filename);
    cmf_range_test_modes(filename,&teems_range_test_initial,&teems_range_test_updated);
    postsim_on=cmf_postsim_on(filename);
    for (nj=0; nj<niodata+noutdata+nsoldata; nj++) logmsg(2,"rank %d logname %s fname %s\n",rank,iodata[nj].logname,iodata[nj].filname);
    if(tab_preprocess(tabfile,newtabfile)==-1)return 0;
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
  }

  strcpy(tabfile,newtabfile);
  if(rank==0)nset=sets_count(tabfile);
  if(nohsl) {
    MPI_Bcast(iodata,niodata*sizeof(cmf_file_entry), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(closure,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(shock,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nset,sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  set_def *sets= (set_def *) calloc (nset,sizeof(set_def));
  for(i=0; i<nset; i++) {
    sets[i].subsetid[0]=i;
    for(j=0; j<MAXSUPSET; j++)sets[i].subsetid[j]=-1;
  }
  if(rank==0) {
    if(sets_read(tabfile,niodata,iodata, sets,nset)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    if(sets_read_intertemporal(tabfile,niodata,iodata, sets,nset)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
    for (i=0; i<nset; i++) {
      if(sets[i].size<0){
        printf("Error: set %s has a negative size in TAB file\n",sets[i].setname);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      sets[i].offset=nsetspace;
      nsetspace=nsetspace+sets[i].size;
    }
  }
  if(nohsl) {
    MPI_Bcast(sets,nset*sizeof(set_def), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nsetspace,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  set_element *set_elems= (set_element *) calloc (nsetspace,sizeof(set_element));
  logmsg(2,"nset %d nsetspace %ld\n",nset,nsetspace);
  for (i=0; i<nsetspace; i++)for (j=0; j<MAXSUPSET; j++)set_elems[i].superset_pos[j]=-1;
  if(rank==0) {
    for (i=0; i<nset; i++) {
      strcpy(vname,sets[i].header);
      nlength=0;
      while (vname[nlength] != '\0') {
        nlength++;
      }
      if (nlength>0) {
        datafile_read_header_info(vname,iodata[sets[i].fileid].filname,&vsize,longname,&dim1);
        datafile_labels *matvar1= (datafile_labels *) calloc (dim1,sizeof(datafile_labels));
        datafile_read_labels(vname,iodata[sets[i].fileid].filname,dim1,matvar1);
        for (j=0; j<dim1; j++) {
          nj=0;
          while(matvar1[j].ch[nj]!='\0') {
            matvar1[j].ch[nj]=tolower((int) matvar1[j].ch[nj]);
            nj++;
          }
          strncpy(set_elems[j+sets[i].offset].setele,matvar1[j].ch,strlen(matvar1[j].ch));
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
                while (str_replace_all(copyline," ", ""));
                readitem = strtok(copyline,",");
                if (readitem==NULL||strlen(readitem)>=NAMESIZE) {
                  printf("Error: malformed element list for set %s\n",sets[i].setname);
                  MPI_Abort(PETSC_COMM_WORLD,1);
                }
                strcpy(set_elems[sets[i].offset].setele,readitem);
                set_elems[sets[i].offset].superset_pos[0]=0;
                for (j=1; j<dim1; j++) {
                  readitem = strtok(NULL,",");
                  if (readitem==NULL||strlen(readitem)>=NAMESIZE) {
                    printf("Error: malformed element list for set %s\n",sets[i].setname);
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
    ndblock1=ndblock;
    /* MAPPING statements (manual 11.9): declarations + by_elements
       values, resolved against the set elements built above so they
       exist before any coefficient/equation machinery runs */
    nmap=tab_count_statements(tabfile,"mapping");
    if(nmap>0) {
      maps= (map_def *) calloc (nmap,sizeof(map_def));
      mappings_read(tabfile,maps,nmap,sets,nset);
      mapping_values_read(tabfile,niodata,iodata,maps,nmap,sets,set_elems);
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
    MPI_Bcast(&nmap,sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    if(nmap>0) {
      if(rank!=0) maps= (map_def *) calloc (nmap,sizeof(map_def));
      MPI_Bcast(maps,nmap*sizeof(map_def), MPI_BYTE,0, PETSC_COMM_WORLD);
      for(i=0; i<nmap; i++) {
        if(rank!=0) maps[i].values= (dim_t *) calloc (sets[maps[i].fromset].size>0?sets[maps[i].fromset].size:1,sizeof(dim_t));
        MPI_Bcast(maps[i].values,sets[maps[i].fromset].size*sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
    }
    /* C2: comp records leave rank 0 (fixed-size PODs; the state
       machinery and the closure/driver dispatch read them) */
    MPI_Bcast(&teems_ncomp,sizeof(dim_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    if(teems_ncomp>0) {
      if(rank!=0) teems_comps= (comp_def *) calloc (teems_ncomp,sizeof(comp_def));
      MPI_Bcast(teems_comps,teems_ncomp*sizeof(comp_def), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
  }
  teems_maps=maps;
  teems_nmap=nmap;
  /* The chain dimension and the diagonal-block partition are derived
     structurally just before the ordering, once the equations are
     readable. Only the bordered methods consume these dimensions. */
  bool structural_time=(matsol==MM_SBBD||matsol==MM_DBBD||matsol==MM_NDBBD);
  bool structural_reg=(matsol==MM_DBBD||matsol==MM_NDBBD);
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
    /* interim M2b guard: mapped equation columns are invisible to the
       bordered methods' block classification, so they run only under LU */
    if(matsol!=MM_LU)mapping_eq_matsol_guard(newtabfile1);
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
  elem_value *elem_vals1;
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
  /* element-level border marks (6.5 E3), populated by the ordering scan
     alongside var_inter on the same ranks */
  bool *ele_inter= (bool *) calloc (nvarele,sizeof(bool));
  if(rank==0) {
    strcpy(commsyntax,"exogenous");
    nexo=closure_read(closure,commsyntax,closure_vals,vars,nvar,sets,nset,set_elems);
    if(nexo==-1)return 0;
    nexo1=nexo;
    strcpy(commsyntax,"shock");
    if(shocks_read(shock,commsyntax,closure_vals,nvarele,vars,nvar,sets,nset,set_elems,subints)==-1)return 0;
    /* backsolve statements: mark the eliminated elements (the flags ride
       the closure broadcast) and check the condensed system's references
       before any equation scan runs with the filter active */
    if(backsolve_read(tabfile,vars,nvar,closure_vals)==-1)return 0;
    if(backsolve_validate_refs(tabfile,vars)==-1)return 0;
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
  /* approximate/accurate-run controls (manual 51.6 semantics, TEEMS
     command-line flags -- teems-R passes them from
     ems_complementarity(); no CMF statements). Defaults: step count
     = the accurate method's step sum, redo on, both runs on, errors
     fatal. PetscOptions are identical on every rank, so no
     broadcasts. comp_acc_phase: 0 = approximate pass (or first
     pass), 1 = the accurate pass after the 51.7.1 closure
     modification (C3). */
  int comp_steps=0,comp_redo=1,comp_do_approx=1,comp_do_acc=1,comp_sberr_warn=0,comp_acc_phase=0;
  double comp_minfrac=0.005;
  if(teems_comp_active>0) {
    PetscReal minfrac_opt=comp_minfrac;
    dim_t iopt;
    /* 51.6 default = the accurate run's step sum; steps2/steps3 were
       folded into ratios above, so rebuild the three counts the
       multistep driver will use */
    comp_steps=(solmethod==SM_GRAGG||solmethod==SM_EULER)
      ?steps1+(int)llround(steps1*step_ratio2)+(int)llround(steps1*step_ratio3)
      :steps1;
    if(comp_steps<1)comp_steps=10;
    iopt=comp_steps;
    PetscOptionsGetInt(NULL,NULL,"-comp_steps",&iopt,NULL);
    if(iopt<1) {
      if(rank==0)printf("Error: -comp_steps must be a positive Euler step count (manual 51.6)\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    comp_steps=(int)iopt;
    iopt=comp_redo;
    PetscOptionsGetInt(NULL,NULL,"-comp_redo",&iopt,NULL);
    comp_redo=iopt?1:0;
    PetscOptionsGetReal(NULL,NULL,"-comp_redo_min_frac",&minfrac_opt,NULL);
    if(minfrac_opt<=0||minfrac_opt>1) {
      if(rank==0)printf("Error: -comp_redo_min_frac must lie in (0,1] (manual 51.6)\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    comp_minfrac=(double)minfrac_opt;
    iopt=comp_do_approx;
    PetscOptionsGetInt(NULL,NULL,"-comp_do_approx",&iopt,NULL);
    comp_do_approx=iopt?1:0;
    iopt=comp_do_acc;
    PetscOptionsGetInt(NULL,NULL,"-comp_do_acc",&iopt,NULL);
    comp_do_acc=iopt?1:0;
    iopt=comp_sberr_warn;
    PetscOptionsGetInt(NULL,NULL,"-comp_sberr_warn",&iopt,NULL);
    comp_sberr_warn=iopt?1:0;
    if(!comp_do_approx&&!comp_do_acc) {
      if(rank==0)printf("Error: -comp_do_approx and -comp_do_acc cannot both be 0\n");
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
    formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
assertions_execute(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,teems_assertions_mode,0);
  }
  gettimeofday(&endtime, NULL);
  if(rank==0)logmsg(1,"Variable calculation time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
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
  }
  //**************************************************************************************
  //****************** END CALCULATE VARIABLE, COEFFICIENT VALUE FROM FORMULA*************
  //**************************************************************************************
  gettimeofday(&begintime, NULL);
  if(rank==0)logmsg(1,"Variable broadcast time %.2f s\n",(begintime.tv_sec - endtime.tv_sec)+((double)(begintime.tv_usec - endtime.tv_usec))/ 1000000);
  //**************************************************************************************
  //****************************** MATRIX FROM FORMULA************************************
  //**************************************************************************************
  VecSize = (PetscInt) (nvarele-nexo-nbselems);
  if(rank==0) {
    if(nbselems>0)logmsg(1,"System size %d equations (%ld exogenous, %ld backsolved)\n",VecSize,nexo,nbselems);
    else logmsg(1,"System size %d equations (%ld exogenous)\n",VecSize,nexo);
  }
  /* Heavy condensation can leave fewer unknowns than exogenous elements.
     The sequential methods handle that with a widened shock vector (BSize
     columns in B), but the parallel-layout preallocation still assumes a
     square B, so the MPI methods abort cleanly instead of corrupting. */
  if(nohsl&&nexo>VecSize) {
    if(rank==0)printf("Error: the condensed system has %d unknowns but %ld exogenous elements; DBBD/NDBBD (-matsol 2/3) do not yet support nexo > system size -- use -matsol 0 (LU) or 1 (SBBD), or condense fewer variables\n",VecSize,nexo);
    PetscFinalize();
    return 1;
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
        chain_flags_apply(sets,nset);
        chain_source="structural";
      }
    }
    if(alltimeset>=0)ntime=sets[alltimeset].size;
    if(structural_reg) {
      allregset=partition_auto_select(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,nvarele,closure_vals,neq,(offset_t)VecSize,alltimeset,ntime,nesteddbbd,(long)mpisize,rank,&partition_auto_json);
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
  if((structural_time||structural_reg)&&matsol==MM_NDBBD&&(alltimeset<0||allregset<0)) {
    if(rank==0) {
      if(alltimeset<0)printf("Error: NDBBD (-matsol 3) needs a chain dimension, but no equation couples set elements through lead/lag offsets.\n");
      if(allregset<0)printf("Error: NDBBD (-matsol 3) needs a diagonal-block partition and no viable set was detected (see the candidate table above); choose another -matsol.\n");
    }
    PetscFinalize();
    return 1;
  }
  if(structural_reg&&matsol==MM_DBBD&&alltimeset<0&&allregset<0) {
    if(rank==0)printf("Error: DBBD (-matsol 2) needs a diagonal-block partition and no viable set was detected (see the candidate table above); use -matsol 0 (LU).\n");
    PetscFinalize();
    return 1;
  }
  /* method-vs-structure check: SBBD without a chain dimension hands
     HSL_MP48 zero blocks — the factorization errors out and the run
     used to finish with exit 0 and no solution. Abort cleanly instead.
     Under HSL only rank 0 resolves alltimeset, so the verdict is
     broadcast to keep the exit collective. */
  if(matsol==MM_SBBD) {
    int sbbd_nochain=(rank==0&&alltimeset<0)?1:0;
    MPI_Bcast(&sbbd_nochain,1,MPI_INT,0,PETSC_COMM_WORLD);
    if(sbbd_nochain) {
      if(rank==0)printf("Error: SBBD (-matsol 1) requires a chain dimension, but the equations couple no set through lead/lag offsets; use -matsol 0 (LU) or -matsol 2 (DBBD) for static models.\n");
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
  offset_t nintraendovar,summat;
  if(rank==rank_hsl) {
    if(nesteddbbd==1) {
      if(!equation_order_read_nested(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,var_inter,ele_inter,eq_defs,eq_intertemp,eq_time,eq_reg,allregset,alltimeset,orderintra,orderreg))MPI_Abort(PETSC_COMM_WORLD,1);
    }
    else if(!equation_order_read(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,var_inter,ele_inter,eq_defs,eq_intertemp,eq_time,eq_reg,allregset,alltimeset,orderintra,orderreg))MPI_Abort(PETSC_COMM_WORLD,1);
    if(alltimeset>=0||allregset>=0)for(i=0; i<neq; i++)eq_intertemp[i]=!eq_intertemp[i];
    j0=0;
    for(i=0; i<nvarele; i++)if(ele_inter[i])j0++;
    if(j0>0&&rank==0)logmsg(1,"Element-level border classification: %ld elements bordered by sliced lead/lag references\n",j0);
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
        if(!closure_vals[j5].is_exogenous&&!closure_vals[j5].is_backsolved) {
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
      if (closure_vals[i].is_exogenous) {
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
          if(!closure_vals[j5].is_exogenous&&!closure_vals[j5].is_backsolved) {
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
        if (closure_vals[i].is_exogenous) {
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
        for (i=0; i<nvar; i++) {
          for (j=0; j<vars[i].nelem; j++) {
            j5=j3+j;
            if(!closure_vals[j5].is_exogenous&&!closure_vals[j5].is_backsolved) {
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
          if (closure_vals[i].is_exogenous) {
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

            if (!closure_vals[j3].is_exogenous&&!closure_vals[j3].is_backsolved) {
              closure_vals[j3].exo_index+=j2;
              j2++;
            }
            if (closure_vals[j3].is_exogenous) {
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
  offset_t *eq_addr= (offset_t *) calloc (VecSize,sizeof(offset_t));//recycle ha_cgeset
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
          if(eq_time[i]>-1)if(eq_defs[i].setid[eq_time[i]]==alltimeset)counteq[set_elems[eq_time_offsets[i]+j2].superset_pos[0]]++;
            else {
              for(j3=1; j3<MAXSUPSET; j3++)if(sets[eq_defs[i].setid[eq_time[i]]].subsetid[j3]==alltimeset)break;
              counteq[set_elems[eq_time_offsets[i]+j2].superset_pos[j3]]++;
            }
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
    logmsg(1,"Border netcut %ld, intra-block equations %ld\n",netcut,nintraeq);
  }
  /* rank 0 always holds valid ordering data: rank_hsl==0 under HSL, and
     under nohsl every rank computes the full ordering */
  if(rank==0)ordering_stats_write(iodata,niodata,noutdata,nsoldata,(long)VecSize,nvarele,nexo,matsol,solmed,nesteddbbd,(long)mpisize,mc66,alltimeset,allregset,sets,ntime,nreg,ndblock,netcut,nintraeq,countvarintra1,counteqnoadd,chain_source,partition_source,partition_auto_json);
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
    VecSetSizes(vece,localsize,VecSize);
  }
  else {
    /* exogenous columns run 0..nexo-1; condensation can push nexo past
       VecSize, so the shock vector and B's columns span BSize */
    VecSetSizes(vece,PETSC_DECIDE,(VecSize>(PetscInt)nexo)?VecSize:(PetscInt)nexo);
  }
  VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
  free(locals);
  if(rank==rank_hsl) {
    ierr = VecGetOwnershipRange(vece,&Istart,&Iend);
    CHKERRQ(ierr);
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
    jacobian_preallocate(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,nexo,closure_vals,ndblock,alltimeset,allregset,eq_intertemp,eq_addr,eq_time,eq_reg,counteq,nintraeq,&sbbd_overrid,Istart,Iend,&dnz,dnnz,&onz,onnz,&dnzB,dnnzB,&onzB,onnzB,nesteddbbd,eqmeta,&neqmeta);
    probe_onfail_context(sets,set_elems,vars,nvar,closure_vals,nvarele,eq_addr,eqmeta,neqmeta,VecSize);
  }
  if(rank==0&&sbbd_overrid&&alltimeset<0) {
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
  FILE* tempvar;
  MPI_Fint fcomm;
  fcomm = MPI_Comm_c2f(PETSC_COMM_WORLD);
  PetscInt nz01,*ai,*aj;
  PetscScalar *vals;
  FILE *ofp;
  int sol;

  fortran_int indata[5];
  indata[1]=VecSize;
  indata[2]=mpisize;
  indata[3]=ndblock;
  indata[4]=nsbbdblocks;
  fortran_int *ptx=NULL;
  ptx = indata;

  solve_real *x1=NULL;//= (ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));
  solve_real *xcf=NULL;
  solve_real *accmetric=NULL; /* embedded-RK cumulative error metrics */
  solve_real *x0=NULL;// (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  solve_real *b1=NULL;//= (ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));
  extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, solve_real *VA, solve_real *B, solve_real *X);
  extern void spec48_single_(fortran_int *indata,int *irn, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
  extern void spec48_nomc66_(fortran_int *indata, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow, MPI_Fint *fcomm,fortran_int *rowptrin, fortran_int *colptrin);

  fortran_int k=0,m=1;
  solve_real temp1,temp2;
  fortran_int tindx1;
  logmsg(2,"rank %d ncof %ld\n",rank,ncof);
  

  if(inmemory) {
    /* Residency cost of skipping the driver spills: value arrays plus the
       Gragg step state. Fall back to scratch files unless it
       fits comfortably in available memory. */
    long need=(long)(ncofele+nvarele)*sizeof(elem_value)
             +(long)nvarele*(sizeof(closure_entry)+6*sizeof(solve_real)+sizeof(int));
    long avail=-1;
    FILE *mi=fopen("/proc/meminfo","r");
    if(mi) {
      char mline[256];
      while(fgets(mline,sizeof(mline),mi))if(sscanf(mline,"MemAvailable: %ld kB",&avail)==1)break;
      fclose(mi);
      if(avail>0)avail*=1024;
    }
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
  bool comp_dispatch=(teems_comp_active>0&&solmethod!=SM_PROBE&&comp_acc_phase==0);
  if(comp_dispatch) {
    if(comp_do_approx) {
      if(rank==0)printf("Complementarity: %ld active (endogenous) component(s); approximate simulation as forward Euler with %d steps (manual 51.1.2)\n",(long)teems_comp_active,comp_steps);
      if(rank==0&&subints>1)printf("Warning: the complementarity approximate run treats the simulation as one interval (per-subinterval approximate+accurate pairs, manual 51.7.4, are not implemented)\n");
      solve_comp_approx(nohsl,VecSize,dnz,dnnz,onz,onnz,dnzB,dnnzB,onzB,onnzB,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,fcomm,comp_steps,comp_redo,comp_minfrac,&xcf);
    }
    else {
      if(rank==0)printf("Complementarity: -comp_do_approx 0; taking the pre-simulation states as the accurate run's targets (manual 51.6)\n");
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
      }
      if(rank==0)printf("Complementarity: accurate simulation with the %s method (closure/shocks modified per manual 51.7.1)\n",solmed);
      /* tear down pass-1 state and re-enter the closure-dependent
         pipeline: exo_index numbering restarts from zero (backsolved
         ordinals are assigned by backsolve_read and keep theirs) */
      for(i=0; i<nvarele; i++)if(!closure_vals[i].is_backsolved)closure_vals[i].exo_index=0;
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
      }
      if(rank==rank_hsl) {
        PetscFree(dnnz);
        PetscFree(onnz);
        PetscFree(dnnzB);
        PetscFree(onnzB);
      }
      var_inter=(bool *) calloc (nvar,sizeof(bool));
      ele_inter=(bool *) calloc (nvarele,sizeof(bool));
      comp_acc_phase=1;
      goto comp_accurate_reentry;
    }
    else {
      if(rank==0)printf("Complementarity: -comp_do_acc 0; the approximate run's solution is the simulation result (manual 51.6)\n");
      comp_states_free();
    }
  }

  if(!comp_dispatch&&solmethod==SM_JOHANSEN)solve_johansen(nohsl,VecSize,A,dnz,dnnz,onz,onnz,B,dnzB,dnnzB,onzB,onnzB,vecb,vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele+nvarele,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,&xcf);

  FILE* solution;

  if(!comp_dispatch&&(solmethod==SM_GRAGG||solmethod==SM_EULER))solve_gragg(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele+nvarele,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,solmethod,&xcf);

  if(!comp_dispatch&&isrk)solve_rk(nohsl,VecSize,dnz,dnnz,onz,onnz,dnzB,dnnzB,onzB,onnzB,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,fcomm,solmethod,adaptive,(double)epstol,(double)retryadj,maxretries,&xcf,&accmetric);

  /* C3: after the accurate run, every component must sit in its
     approximate-run state with the variable inside its bounds
     (manual 51.5.4/51.7.5); fatal unless the CMF downgrades with
     -comp_sberr_warn 1 */
  if(comp_acc_phase==1) {
    offset_t comp_nbad=0;
    if(rank==rank_hsl)comp_nbad=comp_verify_states(sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals);
    if(rank==rank_hsl&&comp_nbad>0) {
      if(comp_sberr_warn)printf("Warning: %ld complementarity state/bound error(s) after the accurate run (treated as warnings per -comp_sberr_warn; check the log carefully, manual 51.6)\n",(long)comp_nbad);
      else {
        printf("Error: %ld complementarity state/bound error(s) after the accurate run; rerun with more Euler steps (-comp_steps) or smaller shocks, or downgrade with -comp_sberr_warn 1 (manual 51.5.4/51.6)\n",(long)comp_nbad);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
    }
    comp_states_free();
  }

  /* -solmed probe: MC79 structural diagnosis rides the probe after the
     (skipped) solve dispatch — assemble the Jacobian and run the
     matching / Dulmage-Mendelsohn analysis with named defects */
  if(solmethod==SM_PROBE) {
    /* realize the pre-simulation complementarity states so the
       E_$comp rows probe with their genuine (state-branch) pattern
       rather than the all-zero weights */
    if(rank==rank_hsl&&teems_ncomp>0) {
      comp_states_set(sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals);
      comp_states_free();
    }
    probe_structural(VecSize,nvarele,ncofele,dnz,dnnz,dnzB,dnnzB,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,closure_vals,ndblock,alltimeset,allregset,eq_addr,counteq,nintraeq,eqmeta,neqmeta,iodata,niodata,noutdata,nsoldata,probefine,mpisize,rank);
    VecDestroy(&vece); /* the skipped solve driver would have destroyed it */
  }
  free(eqmeta);

  jacobian_cache_free();
  backsolve_cache_free();
  lu_fastrefac_free();
  sbbd_fastrefac_free();
  dbbd_fastrefac_free();
  dbbd_fastextract_free();
  ndbbd_fastrefac_free();
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
    strcpy(solchar,tempchar);
    strcat(solchar,".bin");
    logmsg(2,"solchar %s\n",solchar);
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error: cannot open %s for writing\n",solchar);
      return 1;
    }
    if(xcf!=NULL)fwrite(xcf, sizeof(solve_real),nvarele, solution); /* probe runs have no solution; sol.bin stays empty */
    fclose(solution);
    if(accmetric!=NULL) {
      /* embedded-RK cumulative error metrics, one double per variable
         element in .bin order (GEMPACK's <sol>.acc equivalent) */
      strcpy(solchar,tempchar);
      strcat(solchar,".acc");
      logmsg(2,"solchar %s\n",solchar);
      if ( (solution = fopen(solchar, "wb")) == NULL ) {
        printf("Error: cannot open %s for writing\n",solchar);
        return 1;
      }
      fwrite(accmetric, sizeof(solve_real),nvarele, solution);
      fclose(solution);
    }
    strcpy(solchar,tempchar);
    strcat(solchar,".var");
    logmsg(2,"solchar %s\n",solchar);
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error: cannot open %s for writing\n",solchar);
      return 1;
    }
    fwrite(vars, sizeof(array_def),nvar, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".set");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error: cannot open %s for writing\n",solchar);
      return 1;
    }
    fwrite(sets, sizeof(set_def),nset, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".sel");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error: cannot open %s for writing\n",solchar);
      return 1;
    }
    fwrite(set_elems, sizeof(set_element),nsetspace, solution);
    fclose(solution);
    offset_t modeldes[4];
    modeldes[0]=nsetspace;
    modeldes[1]=nvar;
    modeldes[2]=nvarele;
    modeldes[3]=(offset_t)nset;
    strcpy(solchar,tempchar);
    strcat(solchar,".mds");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error: cannot open %s for writing\n",solchar);
      return 1;
    }
    fwrite(modeldes, sizeof(offset_t),4, solution);
    fclose(solution);
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  /* PostSim foundation F3 (early Tier 0): after the solve, coefficient
     slots hold post-simulation (updated) values and xcf holds the
     composed solution; expose the solution to the formula engine and
     evaluate (postsim) assertions. Zero cost when the TAB has none. */
  if(rank==0&&postsim_on&&xcf!=NULL&&elem_vals!=NULL&&(npostsim>0||tab_has_postsim_assertions(tabfile))) {
    postsim_expose_results(elem_vals,ncofele,nvarele,xcf);
    if(npostsim>0) {
      logmsg(1,"postsim: running %d statement(s)\n",npostsim);
      teems_ps_pass=1;
      /* reads first (each statement kind runs as its own scan; within
         a kind, file order -- manual 12.2.1) */
      if(postsim_reads_execute(psfile,niodata,iodata,sets,nset,set_elems,coefs,ncof,ncofele,vars,nvar,nvarele,elem_vals)==-1)MPI_Abort(PETSC_COMM_WORLD,1);
      strcpy(commsyntax,"formula");
      formulas_execute(psfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,true);
      assertions_execute(psfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,true,teems_assertions_mode,0);
      teems_ps_pass=0;
    }
    assertions_execute(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,true,teems_assertions_mode,1);
  }
  if(nowrites==0&&rank==0)for(i=0; i<noutdata; i++){
    outputs_write_csv(tabfile,iodata[i+niodata].logname,iodata[i+niodata].filname,sets,nset,set_elems,coefs,ncof,ncofele,vars,nvar,nvarele,elem_vals);
    logmsg(1,"Wrote %s\n",iodata[i+niodata].logname);
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
  return 0;
}


