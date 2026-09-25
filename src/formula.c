#include <teems_solver.h>

/* Formulas are recompiled every step, so a per-bind warning would repeat
   thousands of times over a solve; warn once per array. */
static char *novalue_warned[2]={NULL,NULL};   /* [0] coefficients, [1] variables */
static offset_t novalue_warned_n[2]={0,0};
static void warn_no_values(const char *name, offset_t idx, int which) {
  /* complementarity-derived names (name_is_comp_derived) take their
     values from C passes, not from reads/formulas */
  if(name_is_comp_derived(name)) return;
  if(idx>=novalue_warned_n[which]) {
    novalue_warned[which]=realloc(novalue_warned[which],idx+1);
    memset(novalue_warned[which]+novalue_warned_n[which],0,idx+1-novalue_warned_n[which]);
    novalue_warned_n[which]=idx+1;
  }
  if(!novalue_warned[which][idx]) {
    novalue_warned[which][idx]=1;
    printf("Warning: %s has no supplied values\n",name);
  }
}

solve_real formula_subst_scalar(char *var2, elem_value *record, array_def *coefs,offset_t ncof) {
  offset_t index;
  if (var2==NULL) return 0;
  str_delete_char(var2,' ');
  solve_real eval=0;
  if (var2[0]>='0'&&var2[0]<='9') {
    eval=atof(var2);
    return eval;
  }
  if (ncof<=0) return eval; /* no coefficient to resolve against: the do-while below started at coefs[-1] */
  index=ncof-1;
  do {
    if (strcmp(coefs[index].cofname,var2)==0) {
      eval=record[coefs[index].offset].value;
      break;
    }
  } while (index--);
  return eval;
}


/* split a lowered mapping token "map~idx" (mapping_lower_calls):
   returns the index part and sets *mp to mapping id + 1, or returns
   the token unchanged with *mp = 0.  Shared with the jacobian-side
   linear-variable index walk (design doc M2b). */
char *mapping_token_split(char *p, int *mp) {
  char *at;
  dim_t m;
  *mp=0;
  if (p==NULL) return p;
  at=strchr(p,MAPMARK);
  if (at==NULL) return p;
  if (strchr(at+1,MAPMARK)!=NULL) {
    errmsg("Error: composition of set mappings is not supported (manual 11.9.6)\n");
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  *at='\0';
  for (m=0; m<teems_nmap; m++) if (strcmp(p,teems_maps[m].mapname)==0) {
      *mp=(int)m+1;
      return at+1;
    }
  errmsg("Error: unknown mapping %s in an index expression\n",p);
  MPI_Abort(PETSC_COMM_WORLD,1);
  return at+1;
}

/* bind one mapped dimension (manual 11.9.4/11.9.7): the index must
   range over the mapping's domain set exactly and the argument
   position must be the codomain set exactly (subset routing around a
   mapped argument is deferred -- named fatal, not a mis-bind) */
static void map_dim_bind(dim_addr *Dm, int mp, dim_t frame_setid, offset_t arg_setid, offset_t stride, int leadlag, const char *symname, set_def *sets) {
  map_def *md=&teems_maps[mp-1];
  dim_t dss=0,css=0;
  /* the index may range over a declared subset of the mapping's domain
     (its elements route to their domain positions first), and the
     codomain may be a declared subset of the argument set (the mapped
     position routes on into the argument set), manual 11.9.7 --
     potential-models vetting G9/G17/G-S5. Anything else stays fatal. */
  if ((dim_t)md->fromset!=frame_setid) {
    dss=set_supset_slot(sets,frame_setid,(dim_t)md->fromset);
    if (dss<0) {
      errmsg("Error: the index of mapping %s does not range over its domain set %s or a declared subset of it (it ranges over %s, in %s)\n",md->mapname,sets[md->fromset].setname,sets[frame_setid].setname,symname);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
  if ((offset_t)md->toset!=arg_setid) {
    css=set_supset_slot(sets,(dim_t)md->toset,(dim_t)arg_setid);
    if (css<0) {
      errmsg("Error: mapping %s does not map into the argument set at that position of %s (its codomain %s is neither %s nor a declared subset of it; manual 11.9.7)\n",md->mapname,symname,sets[md->toset].setname,sets[arg_setid].setname);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
  if (!md->has_values) {
    errmsg("Error: mapping %s is used (in %s) before a Formula has assigned all of its values (manual 10.13.1/11.9.1)\n",md->mapname,symname);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  Dm->ADims=stride;
  Dm->leadlag=leadlag;
  Dm->SupSet=(css>0)?1:0;
  Dm->SSIndx=(int)css;
  Dm->MapId=mp;
  Dm->MapDomSS=(int)dss;
  md->used=true;
}

/* bind one argument position of a coefficient/variable reference to
   the frame dimension its index ranges over. A repeated index
   (T(d,d), TRADE(c,"dom",r,r)) binds the same frame dimension at
   several positions: each position contributes its own stride, so the
   contributions add -- the stride used to be assigned, and T(d,d) read
   T(first,d) (potential-models vetting SW2). Positions routed alike
   (same superset slot, mapping, lead/lag) share one accumulated
   stride; a second, differently routed group takes the Rep slot. */
static void dim_bind_pos(dim_addr *D, bool *seen, char *p, int mp, int leadlag, dim_t frame_setid, offset_t arg_setid, offset_t stride, const char *symname, set_def *sets) {
  dim_addr b;
  memset(&b,0,sizeof(b));
  if (mp>0) map_dim_bind(&b,mp,frame_setid,arg_setid,stride,leadlag,symname,sets);
  else {
    dim_t ss=set_supset_slot(sets,frame_setid,(dim_t)arg_setid);
    if (ss<0) set_supset_fatal(p,symname,NULL,sets,frame_setid,(dim_t)arg_setid);
    if (ss>0) { b.SupSet=1; b.SSIndx=(int)ss; }
    b.ADims=stride;
    b.leadlag=leadlag;
  }
  if (!*seen) {
    *D=b;
    *seen=true;
    return;
  }
  if (D->MapId==b.MapId&&D->MapDomSS==b.MapDomSS&&D->SupSet==b.SupSet&&D->SSIndx==b.SSIndx&&D->leadlag==b.leadlag) {
    D->ADims+=b.ADims;
    return;
  }
  if (!D->Rep) {
    D->Rep=1;
    D->ADims2=b.ADims;
    D->SupSet2=b.SupSet;
    D->SSIndx2=b.SSIndx;
    D->leadlag2=b.leadlag;
    D->MapId2=b.MapId;
    D->MapDomSS2=b.MapDomSS;
    return;
  }
  if (D->MapId2==b.MapId&&D->MapDomSS2==b.MapDomSS&&D->SupSet2==b.SupSet&&D->SSIndx2==b.SSIndx&&D->leadlag2==b.leadlag) {
    D->ADims2+=b.ADims;
    return;
  }
  errmsg("Error: index %s is repeated in %s at argument positions needing three different set routings; this is not supported\n",p,symname);
  MPI_Abort(PETSC_COMM_WORLD,1);
}

/* $POS(...) (manual 11.5.6) inside a formula text (spaces stripped,
   lowercased, mapping calls already lowered to map~idx): each call is
   compiled to an OP_LOAD/OT_POS op emitted ahead of the expression and
   replaced in the text by that op's temp name, which the operand binder
   resolves like any other generated temp.  Forms:
     $pos(i)        position of quantifier index i in the set it ranges over
     $pos(i,S)      its position in S, a superset of that set
     $pos(map~i)    position of the mapped element in the codomain (11.9.5)
     $pos(map~i,S)  ... lifted into a superset S of the codomain
     $pos("el",S)   a constant: the position of the element in S
   Returns 0 (with a message) on a malformed or unresolvable call. */
static int pos_lower(char *f, set_def *sets, quantifier *arSet, dim_t fdim, formula_op *ops, dim_t *nops) {
  char *q;
  /* index lists ride inside {} once formula_normalize has run, so the
     call may read $pos{...} as well as $pos(...) */
  while ((q=strstr(f,"$pos("))!=NULL||(q=strstr(f,"$pos{"))!=NULL) {
    char inner[NAMESIZE],arg1[NAMESIZE],arg2[NAMESIZE],repl[NAMESIZE];
    char *p=q+5,*close;
    int depth=1,k=0;
    for (close=p; *close!='\0'; close++) {
      if (*close=='('||*close=='{') depth++;
      else if (*close==')'||*close=='}') { depth--; if (depth==0) break; }
    }
    if (*close=='\0') { errmsg("Error: unbalanced $POS( in formula: %s\n",f); return 0; }
    if (close-p>=NAMESIZE-1) { errmsg("Error: $POS argument too long in formula: %s\n",f); return 0; }
    strncpy(inner,p,close-p);
    inner[close-p]='\0';
    /* split on the top-level comma */
    arg2[0]='\0';
    {
      char *c=inner;
      int d=0;
      for (; *c!='\0'; c++) {
        if (*c=='"') { c++; while (*c!='"'&&*c!='\0') c++; if (*c=='\0') break; continue; }
        if (*c=='('||*c=='{') d++;
        else if (*c==')'||*c=='}') d--;
        else if (*c==','&&d==0) break;
      }
      if (*c==',') { *c='\0'; strcpy(arg2,c+1); }
      strcpy(arg1,inner);
    }
    if (arg1[0]=='\0') { errmsg("Error: $POS needs an index, element or index expression argument (manual 11.5.6): %s\n",f); return 0; }
    dim_t setS=-1;
    if (arg2[0]!='\0') {
      for (setS=0; setS<teems_nset; setS++) if (strcmp(arg2,teems_sets[setS].setname)==0) break;
      if (setS==teems_nset) { errmsg("Error: $POS: %s is not a declared set (manual 11.5.6)\n",arg2); return 0; }
    }
    if (arg1[0]=='"') {
      /* $POS("el",S): a constant */
      char el[NAMESIZE];
      dim_t e;
      if (setS<0) { errmsg("Error: $POS(\"%s\") needs the set as second argument (manual 11.5.6)\n",arg1); return 0; }
      k=0;
      for (p=arg1+1; *p!='"'&&*p!='\0'&&k<NAMESIZE-1; p++) el[k++]=tolower((int)*p);
      el[k]='\0';
      for (e=0; e<teems_sets[setS].size; e++) if (strcmp(el,teems_set_elems[teems_sets[setS].offset+e].setele)==0) break;
      if (e==teems_sets[setS].size) { errmsg("Error: $POS: %s is not an element of set %s (manual 11.5.6)\n",el,arg2); return 0; }
      snprintf(repl,NAMESIZE,"%d",(int)e+1);
    } else {
      char *at=strchr(arg1,MAPMARK),*idx=arg1;
      int mp=0;
      dim_t l,base,sup=0;
      if (at!=NULL) {
        dim_t m;
        *at='\0';
        idx=at+1;
        for (m=0; m<teems_nmap; m++) if (strcmp(arg1,teems_maps[m].mapname)==0) { mp=(int)m+1; break; }
        if (mp==0) { errmsg("Error: $POS: unknown mapping %s in an index expression\n",arg1); return 0; }
      }
      for (l=0; l<fdim; l++) if (strcmp(idx,arSet[l].index_name)==0) break;
      if (l==fdim) { errmsg("Error: $POS: %s is not an index of the statement's quantifiers (manual 11.5.6): %s\n",idx,f); return 0; }
      base=(dim_t)arSet[l].setid;
      if (mp>0) {
        if ((dim_t)teems_maps[mp-1].fromset!=base) { errmsg("Error: $POS: the index of mapping %s does not range over its domain set\n",teems_maps[mp-1].mapname); return 0; }
        if (!teems_maps[mp-1].has_values) { errmsg("Error: mapping %s is used (in $POS) before a Formula has assigned all of its values (manual 10.13.1/11.9.1)\n",teems_maps[mp-1].mapname); return 0; }
        teems_maps[mp-1].used=true;
        base=(dim_t)teems_maps[mp-1].toset;
      }
      if (setS>=0&&setS!=base) {
        for (sup=1; sup<MAXSUPSET; sup++) if (sets[base].subsetid[sup]==setS) break;
        if (sup==MAXSUPSET) { errmsg("Error: $POS(%s,%s): %s does not range over a subset of %s (manual 11.5.6)\n",idx,arg2,idx,arg2); return 0; }
      }
      ops[*nops].Oper=OP_LOAD;
      ops[*nops].Var1Type=OT_POS;
      ops[*nops].Var1BegAdd=(offset_t)l;
      ops[*nops].Var1Dims[0].MapId=mp;
      ops[*nops].Var1Dims[0].SupSet=(setS>=0&&setS!=base)?1:0;
      ops[*nops].Var1Dims[0].SSIndx=(int)sup;
      ops[*nops].Var1Dims[0].ADims=(offset_t)base;
      ops[*nops].Var1Dims[0].leadlag=0;
      snprintf(ops[*nops].TmpVarName,sizeof(ops[*nops].TmpVarName),"gen_pos%03d",(int)*nops);
      strcpy(repl,ops[*nops].TmpVarName);
      *nops=*nops+1;
    }
    /* splice: text before $pos( + replacement + text after the ')' */
    {
      char tail[TABREADLINE];
      strcpy(tail,close+1);
      if ((q-f)+strlen(repl)+strlen(tail)>=TABREADLINE) { errmsg("Error: formula too long after $POS expansion: %s\n",f); return 0; }
      *q='\0';
      strcat(f,repl);
      strcat(f,tail);
    }
  }
  return 1;
}

/* next index token of a reference's argument list: an empty index
   (x( ), x(c,,t) -- normalize strips the blanks) left strtok NULL and the
   binder compared it (SEGV, fuzz batch 13 with the real declaration
   context); abort by name instead */
static char *bind_next_index(const char *delim, const char *tokcopy) {
  char *p=strtok(NULL,delim);
  if (p==NULL) {
    errmsg("Error: %s has an empty index in its argument list; a reference must carry exactly the declared indices (manual 10.3, 11.4.10)\n",tokcopy);
    fflush(stdout);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  return p;
}

int formula_bind_operand(char *var2, set_def *sets,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int nops,quantifier *arSet,dim_t fdim,int varindex) {
  offset_t index;
  char *p=NULL;//,copyvar[TABREADLINE];//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  dim_t l1,l;//,svar2;//=0,i2=0,i3=0,i4=0,svar1,svar2,checkvar20=0,checkvar10=0,checkvar11=0,checkvar12=0,checkvar16=0,checkvar14=0,l;
  int leadlag;
  int mp=0;
  bool IsChange=false;
  bool seen[MAXVARDIM]={false};
  /* index-count guard: a reference must carry exactly the declared
     number of indices (manual 10.3/11.4.10). The binders below strtok
     the argument list positionally, so a bare reference to an indexed
     array (esun for esun(t)) handed a NULL index to parse_index_leadlag
     (SEGV in an equation term) or bound element 0 in a formula, and a
     surplus index was silently dropped. Counted before strtok cuts the
     token; the copy keeps the reference text for the message. */
  char tokcopy[NAMESIZE];
  int nargs_tok=0;
  {
    const char *b=strchr(var2,'{');
    strncpy(tokcopy,var2,NAMESIZE-1);
    tokcopy[NAMESIZE-1]='\0';
    if(b!=NULL) {
      const char *q;
      nargs_tok=1;
      for(q=b+1; *q!='\0'&&*q!='}'; q++) if(*q==',') nargs_tok++;
    }
  }
  {
    /* ops slots are reused across compiles: the repeat/subset routing
       fields are set only by the binders below, so clear them first */
    dim_addr *D=(varindex==2)?ops[nops].Var2Dims:ops[nops].Var1Dims;
    for (l=0; l<MAXVARDIM; l++) { D[l].Rep=0; D[l].MapDomSS=0; }
  }
  p= strtok(var2,"{");
  if (p==NULL) {
    p=&var2[0];
  }
  if(var2[0]=='p'&&var2[1]=='_') {
    IsChange=true;
    p++;
    p++;
  }
  /* p_-prefixed tokens name LINEAR variables: when the token resolves
     to one -- bare TEEMS convention, the C0 levels auto-pair, or a
     declared p_-/c_-leading name (design doc sections 5/6) -- the
     variable wins and the coefficient search is skipped; bare tokens
     keep coefficient-first (value reference). Bit-neutral otherwise:
     11.2.1 uniqueness plus the names_validate ambiguity fatal make
     the resolution unique. */
  bool PairSkipCoefs=false;
  if(IsChange) {
    offset_t lr=linvar_resolve(p,vars,nvar);
    if(lr>=0) {
      PairSkipCoefs=true;
      /* redirect the search token to the declared name so the
         variable loop below binds the resolved record */
      p=vars[lr].cofname;
    }
  }
  index=ncof-1;
  if(!PairSkipCoefs) do {
    if (strcmp(coefs[index].cofname,p)==0) {
      if((int)coefs[index].size!=nargs_tok) {
        errmsg("Error: coefficient %s is declared with %d %s but is referenced with %d in %s; a reference must carry exactly the declared indices (manual 10.3, 11.4.10)\n",coefs[index].cofname,(int)coefs[index].size,coefs[index].size==1?"index":"indices",nargs_tok,tokcopy);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      if(!coefs[index].suplval)warn_no_values(coefs[index].cofname,index,0);
      if(varindex==2) {
        ops[nops].Var2BegAdd=coefs[index].offset;
        ops[nops].Var2Type=OT_ARRAY;
      } else {
        ops[nops].Var1BegAdd=coefs[index].offset;
        ops[nops].Var1Type=OT_ARRAY;
      }
      switch(coefs[index].size) {
      case 0:
        for (l=0; l<fdim; l++) {
          if(varindex==2) {
            ops[nops].Var2Dims[l].ADims=0;
            ops[nops].Var2Dims[l].MapId=0;
          } else {
            ops[nops].Var1Dims[l].ADims=0;
            ops[nops].Var1Dims[l].MapId=0;
          }
        }
        break;
      case 1:
        p=bind_next_index("}",tokcopy);
        leadlag=0;
        parse_index_leadlag(p,&leadlag);
          p=mapping_token_split(p,&mp);
       for (l=0; l<fdim; l++) {
             if(varindex==2) {
              ops[nops].Var2Dims[l].ADims=0;
              ops[nops].Var2Dims[l].leadlag=0;
              ops[nops].Var2Dims[l].SupSet=0;
              ops[nops].Var2Dims[l].MapId=0;
            } else {
              ops[nops].Var1Dims[l].ADims=0;
              ops[nops].Var1Dims[l].leadlag=0;
              ops[nops].Var1Dims[l].SupSet=0;
              ops[nops].Var1Dims[l].MapId=0;
            }
          if (strcmp(p,arSet[l].index_name)==0) {
            if(mp>0) {
              map_dim_bind(varindex==2?&ops[nops].Var2Dims[l]:&ops[nops].Var1Dims[l],mp,arSet[l].setid,coefs[index].setid[0],coefs[index].strides[0],leadlag,coefs[index].cofname,sets);
            } else
            if(varindex==2) {
              { dim_t ss=set_supset_slot(sets,arSet[l].setid,coefs[index].setid[0]); if(ss<0)set_supset_fatal(p,coefs[index].cofname,NULL,sets,arSet[l].setid,coefs[index].setid[0]); if(ss>0){ops[nops].Var2Dims[l].SupSet=1; ops[nops].Var2Dims[l].SSIndx=(int)ss;} }
              ops[nops].Var2Dims[l].ADims=coefs[index].strides[0];
              ops[nops].Var2Dims[l].leadlag=leadlag;
            } else {
              { dim_t ss=set_supset_slot(sets,arSet[l].setid,coefs[index].setid[0]); if(ss<0)set_supset_fatal(p,coefs[index].cofname,NULL,sets,arSet[l].setid,coefs[index].setid[0]); if(ss>0){ops[nops].Var1Dims[l].SupSet=1; ops[nops].Var1Dims[l].SSIndx=(int)ss;} }
              ops[nops].Var1Dims[l].ADims=coefs[index].strides[0];
              ops[nops].Var1Dims[l].leadlag=leadlag;
            }
          }
        }
        break;
      default:
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) {
            ops[nops].Var2Dims[l1].ADims=0;
            ops[nops].Var2Dims[l1].leadlag=0;
            ops[nops].Var2Dims[l1].SupSet=0;
            ops[nops].Var2Dims[l1].MapId=0;
          } else {
            ops[nops].Var1Dims[l1].ADims=0;
            ops[nops].Var1Dims[l1].leadlag=0;
            ops[nops].Var1Dims[l1].SupSet=0;
            ops[nops].Var1Dims[l1].MapId=0;
          }
        }
        for (l=0; l<coefs[index].size-1; l++) {
          p=bind_next_index(",",tokcopy);
          leadlag=0;
          parse_index_leadlag(p,&leadlag);
          p=mapping_token_split(p,&mp);
          for (l1=0; l1<fdim; l1++) {
            if (strcmp(p,arSet[l1].index_name)==0) {
              dim_bind_pos(varindex==2?&ops[nops].Var2Dims[l1]:&ops[nops].Var1Dims[l1],&seen[l1],p,mp,leadlag,arSet[l1].setid,coefs[index].setid[l],coefs[index].strides[l],coefs[index].cofname,sets);
              break;
            }
          }
        }
        p=bind_next_index("}",tokcopy);
        leadlag=0;
        parse_index_leadlag(p,&leadlag);
          p=mapping_token_split(p,&mp);
        for (l1=0; l1<fdim; l1++) {
          if (strcmp(p,arSet[l1].index_name)==0) {
            dim_bind_pos(varindex==2?&ops[nops].Var2Dims[l1]:&ops[nops].Var1Dims[l1],&seen[l1],p,mp,leadlag,arSet[l1].setid,coefs[index].setid[l],coefs[index].strides[l],coefs[index].cofname,sets);
            break;
          }
        }
        break;
      }
      return 1;
    }
  } while (index--);

  index=nvar-1;
  do {
    if (strcmp(vars[index].cofname,p)==0) {
      if((int)vars[index].size!=nargs_tok) {
        errmsg("Error: variable %s is declared with %d %s but is referenced with %d in %s; a reference must carry exactly the declared indices (manual 10.3, 11.4.10)\n",vars[index].cofname,(int)vars[index].size,vars[index].size==1?"index":"indices",nargs_tok,tokcopy);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      if(varindex==2) {
        ops[nops].Var2BegAdd=vars[index].offset+ncofele;
        if(IsChange) ops[nops].Var2Type=OT_CHANGE;
        else ops[nops].Var2Type=OT_ARRAY;
      } else {
        ops[nops].Var1BegAdd=vars[index].offset+ncofele;
        if(IsChange) ops[nops].Var1Type=OT_CHANGE;
        else ops[nops].Var1Type=OT_ARRAY;
      }
      switch(vars[index].size) {
      case 0:
        for (l=0; l<fdim; l++) {
          if(varindex==2) {
            ops[nops].Var2Dims[l].ADims=0;
            ops[nops].Var2Dims[l].leadlag=0;
            ops[nops].Var2Dims[l].MapId=0;
          } else {
            ops[nops].Var1Dims[l].ADims=0;
            ops[nops].Var1Dims[l].leadlag=0;
            ops[nops].Var1Dims[l].MapId=0;
          }
        }
        break;
      case 1:
        p=bind_next_index("}",tokcopy);
        leadlag=0;
        parse_index_leadlag(p,&leadlag);
          p=mapping_token_split(p,&mp);
        for (l=0; l<fdim; l++) {
             if(varindex==2) {
              ops[nops].Var2Dims[l].ADims=0;
              ops[nops].Var2Dims[l].leadlag=0;
              ops[nops].Var2Dims[l].SupSet=0;
              ops[nops].Var2Dims[l].MapId=0;
            } else {
              ops[nops].Var1Dims[l].ADims=0;
              ops[nops].Var1Dims[l].leadlag=0;
              ops[nops].Var1Dims[l].SupSet=0;
              ops[nops].Var1Dims[l].MapId=0;
            }
          if (strcmp(p,arSet[l].index_name)==0) {
            if(mp>0) {
              map_dim_bind(varindex==2?&ops[nops].Var2Dims[l]:&ops[nops].Var1Dims[l],mp,arSet[l].setid,vars[index].setid[0],vars[index].strides[0],leadlag,vars[index].cofname,sets);
            } else
            if(varindex==2) {
              { dim_t ss=set_supset_slot(sets,arSet[l].setid,vars[index].setid[0]); if(ss<0)set_supset_fatal(p,vars[index].cofname,NULL,sets,arSet[l].setid,vars[index].setid[0]); if(ss>0){ops[nops].Var2Dims[l].SupSet=1; ops[nops].Var2Dims[l].SSIndx=(int)ss;} }
              ops[nops].Var2Dims[l].ADims=vars[index].strides[0];
              ops[nops].Var2Dims[l].leadlag=leadlag;
            } else {
              { dim_t ss=set_supset_slot(sets,arSet[l].setid,vars[index].setid[0]); if(ss<0)set_supset_fatal(p,vars[index].cofname,NULL,sets,arSet[l].setid,vars[index].setid[0]); if(ss>0){ops[nops].Var1Dims[l].SupSet=1; ops[nops].Var1Dims[l].SSIndx=(int)ss;} }
              ops[nops].Var1Dims[l].ADims=vars[index].strides[0];
              ops[nops].Var1Dims[l].leadlag=leadlag;
            }
          }
        }
        break;
      default:
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) {
            ops[nops].Var2Dims[l1].ADims=0;
            ops[nops].Var2Dims[l1].leadlag=0;
            ops[nops].Var2Dims[l1].SupSet=0;
            ops[nops].Var2Dims[l1].MapId=0;
          } else {
            ops[nops].Var1Dims[l1].ADims=0;
            ops[nops].Var1Dims[l1].leadlag=0;
            ops[nops].Var1Dims[l1].SupSet=0;
            ops[nops].Var1Dims[l1].MapId=0;
          }
        }
        for (l=0; l<vars[index].size-1; l++) {
          p=bind_next_index(",",tokcopy);
          leadlag=0;
          parse_index_leadlag(p,&leadlag);
          p=mapping_token_split(p,&mp);
          for (l1=0; l1<fdim; l1++) {
            if (strcmp(p,arSet[l1].index_name)==0) {
              dim_bind_pos(varindex==2?&ops[nops].Var2Dims[l1]:&ops[nops].Var1Dims[l1],&seen[l1],p,mp,leadlag,arSet[l1].setid,vars[index].setid[l],vars[index].strides[l],vars[index].cofname,sets);
              break;
            }
          }
        }
        p=bind_next_index("}",tokcopy);
        leadlag=0;
        parse_index_leadlag(p,&leadlag);
          p=mapping_token_split(p,&mp);
        for (l1=0; l1<fdim; l1++) {
          if (strcmp(p,arSet[l1].index_name)==0) {
            dim_bind_pos(varindex==2?&ops[nops].Var2Dims[l1]:&ops[nops].Var1Dims[l1],&seen[l1],p,mp,leadlag,arSet[l1].setid,vars[index].setid[l],vars[index].strides[l],vars[index].cofname,sets);
            break;
          }
        }
        break;
      }
      return 1;
    }
  } while (index--);

  for (index=totalsum-1; index>-1; index--) {
    if (strcmp(sum_cof[index].sumname,p)==0) {
      if(varindex==2) {
        ops[nops].Var2BegAdd=sum_cof[index].offset;
        ops[nops].Var2Type=OT_SUM;
        for (l1=0; l1<sum_cof[index].size; l1++)ops[nops].Var2Dims[l1].SSIndx=0;
      } else {
        ops[nops].Var1BegAdd=sum_cof[index].offset;
        ops[nops].Var1Type=OT_SUM;
        for (l1=0; l1<sum_cof[index].size; l1++)ops[nops].Var1Dims[l1].SSIndx=0;
      }
      /* ops slots are reused across the sum-body compiles and the outer
         statement compile of one statement, and this block only ever
         wrote ADims: a generated-sum dim inheriting a stale MapId from
         a mapped body operand routed the OUTER read through the
         mapping's value table (M2c kit, carried leg) -- clear the
         routing fields for every frame dim */
      for (l1=0; l1<fdim; l1++) {
        if(varindex==2) {
          ops[nops].Var2Dims[l1].leadlag=0;
          ops[nops].Var2Dims[l1].SupSet=0;
          ops[nops].Var2Dims[l1].MapId=0;
        } else {
          ops[nops].Var1Dims[l1].leadlag=0;
          ops[nops].Var1Dims[l1].SupSet=0;
          ops[nops].Var1Dims[l1].MapId=0;
        }
      }
      switch(sum_cof[index].size) {
      case 0:
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) ops[nops].Var2Dims[l1].ADims=0;
          else ops[nops].Var1Dims[l1].ADims=0;
        }
        break;
      case 1:
        p=bind_next_index("}",tokcopy);
        for (l=0; l<fdim; l++) {
          if (strcmp(p,arSet[l].index_name)==0) {
            if(varindex==2) {
              ops[nops].Var2Dims[l].ADims=sum_cof[index].strides[0];
            } else {
              ops[nops].Var1Dims[l].ADims=sum_cof[index].strides[0];
            }
          } else {
            if(varindex==2) {
              ops[nops].Var2Dims[l].ADims=0;
            } else {
              ops[nops].Var1Dims[l].ADims=0;
            }
          }
        }
        break;
      default:
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) ops[nops].Var2Dims[l1].ADims=0;
          else ops[nops].Var1Dims[l1].ADims=0;
        }
        for (l=0; l<sum_cof[index].size-1; l++) {
          p=bind_next_index(",",tokcopy);
          for (l1=0; l1<fdim; l1++) {
            if (strcmp(p,arSet[l1].index_name)==0) {
              if(varindex==2) {
                ops[nops].Var2Dims[l1].ADims=sum_cof[index].strides[l];
              } else {
                ops[nops].Var1Dims[l1].ADims=sum_cof[index].strides[l];
              }
              break;
            }
          }
        }
        p=bind_next_index("}",tokcopy);
        for (l1=0; l1<fdim; l1++) {
          if (strcmp(p,arSet[l1].index_name)==0) {
            if(varindex==2) {
              ops[nops].Var2Dims[l1].ADims=sum_cof[index].strides[l];
            } else {
              ops[nops].Var1Dims[l1].ADims=sum_cof[index].strides[l];
            }
            break;
          }
        }
        break;
      }
      return 1;
    }
  }

  for (l1=0; l1<nops; l1++) {
    if (strcmp(var2,ops[l1].TmpVarName)==0) {
      if(varindex==2) {
        ops[nops].Var2Type=OT_TEMP;
        ops[nops].Var2BegAdd=l1;
      } else {
        ops[nops].Var1Type=OT_TEMP;
        ops[nops].Var1BegAdd=l1;
      }
      return 1;
    }
  }
  if (var2[0]>='0'&&var2[0]<='9') {
    if(varindex==2) {
      ops[nops].Var2Type=OT_CONST;
      ops[nops].Var2Val=atof(var2);
    } else {
      ops[nops].Var1Type=OT_CONST;
      ops[nops].Var1Val=atof(var2);
    }
    return 1;
  }
  {
    /* a signed or leading-dot number is a constant only when the whole
       token reads as one */
    char *e;
    double v=strtod(var2,&e);
    if (e!=var2&&*e=='\0') {
      if(varindex==2) {
        ops[nops].Var2Type=OT_CONST;
        ops[nops].Var2Val=v;
      } else {
        ops[nops].Var1Type=OT_CONST;
        ops[nops].Var1Val=v;
      }
      return 1;
    }
  }
  /* anything else -- a quantifier index, a quoted element, an undeclared
     name -- used to fall through with the operand left unset and
     evaluate as 0 without a message (the IF rewrite's `[r] - ["usa"]`
     helper solved to 0 everywhere, so every such condition held).
     Callers do not test the return value, so abort here by name. */
  errmsg("Error: %s is not a coefficient, variable or number and cannot be an arithmetic operand (an index or quoted element compares through $POS, manual 11.5.6/11.4.11)\n",tokcopy);
  fflush(stdout);
  MPI_Abort(PETSC_COMM_WORLD,1);
  return 0;
}

int leadlag_encode(char *line) {
  char *p1,*p2,*p3;
  p1=strchr(line,'{');
  while (p1!=NULL) {
    p3=strchr(p1,'}');
    if(p3==NULL)break;
    p2=strchr(p1,'+');
    if(p2!=NULL) {
      while(p3-p2>0) {
        *p2='#';
        p2=strchr(p1,'+');
        if(p2==NULL)p2=p3;
      }
    }
    p2=strchr(p1,'-');
    if(p2!=NULL) {
      while(p3-p2>0) {
        *p2='!';
        p2=strchr(p1,'-');
        if(p2==NULL)p2=p3;
      }
    }
    p1=strchr(p3,'{');
  }
  return 1;
}

int parse_index_leadlag(char *p,int *leadlag) {
  char *plussign,*minsign;
  if(p==NULL) return 0;
  plussign=strchr(p,'#');
  minsign=strchr(p,'!');
  if(plussign!=NULL) {
    *leadlag=atoi(plussign+1);
    *plussign='\0';
  }
  if(minsign!=NULL) {
    *leadlag=-atoi(minsign+1);
    *minsign='\0';
  }
  return 1;
}

int formula_compile(char *fomulain, set_def *sets,array_def *coefs, offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,dim_t totalsum,formula_op *ops,dim_t *nops,quantifier *arSet,dim_t fdim) {
  int npar=0,npow=0,nmul=0,ndiv=0,nplu=0,nmin=0,j;
  *nops=0;
  /* $POS (manual 11.5.6): compiled to position ops ahead of the
     expression, the calls replaced by their temp names */
  if(strstr(fomulain,"$pos(")!=NULL||strstr(fomulain,"$pos{")!=NULL) {
    if(!pos_lower(fomulain,sets,arSet,fdim,ops,nops)) return 0;
  } else if(strstr(fomulain,"$pos")!=NULL) {
    errmsg("Error: malformed $POS call in formula: %s\n",fomulain);
    return 0;
  }
  npar=str_count_char(fomulain, ')');
  if (npar==0) {
    npow=str_count_char(fomulain, '^');
    if (npow>0) {
      if(!formula_compile_pow(fomulain,sets,npow,0,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
    }
    nmul=str_count_char(fomulain, '*');
    ndiv=str_count_char(fomulain, '/');
    nmul=nmul+ndiv;
    if (nmul>0) {
      if(!formula_compile_muldiv(fomulain,sets,nmul,0,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
    }
    nplu=str_count_char(fomulain, '+');
    nmin=str_count_char(fomulain, '-');
    nplu=nplu+nmin;
    if (nplu>0) {
      if(!formula_compile_addsub(fomulain,sets,nplu,0,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
    }
      formula_bind_operand(fomulain,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
      ops[*nops].Oper=OP_LOAD;
      ops[*nops].TmpVarName[0]='\0';
      *nops=*nops+1;
    return 1;
  }
  int i;
  char *p=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE];
  char interchar[TABREADLINE],interchar1[TABREADLINE];
  /* fresh group numbers for multi-arg function arguments: continue
     past the paren-group range so gen_* temp names never collide */
  int arggrp=npar+2;

  for (i=1; i<npar+2; i++) {
    p=strchr(fomulain,')');
    if (p!=NULL) {
      strncpy(fpart1, fomulain, p-fomulain);
      fpart1[p-fomulain] = '\0';
      strcpy(fpart3,p+1);
      p=strrchr(fpart1,'(');
      if (p!=NULL) {
        strcpy(fpart2, p+1);
        strncpy(fpart1,fomulain,p-fpart1);
        fpart1[p-fpart1] = '\0';
      } else {
        errmsg("Error: malformed formula (unbalanced parentheses)\n");
        return 0;
      }
    } else {
      fpart1[0]='\0';
      strcpy(fpart2,fomulain);
      fpart3[0]='\0';
    }
    /* multi-arg intrinsics MAX/MIN/ID0V (manual 11.5/11.5.1): split
       the group on brace-depth-0 commas (indices ride inside {}),
       compile each argument to its own temp under a fresh group
       number, then fold pairwise; ID0V(x,v) = x unless x is 0 */
    int fnmulti=0;
    j=strlen(fpart1);
    if (j>=3&&strncmp(fpart1+j-3,"max",3)==0&&(j==3||fpart1[j-4]==' '||fpart1[j-4]=='('||fpart1[j-4]=='+'||fpart1[j-4]=='-'||fpart1[j-4]=='*'||fpart1[j-4]=='/'||fpart1[j-4]=='^'||fpart1[j-4]==','))fnmulti=1;
    if (j>=3&&strncmp(fpart1+j-3,"min",3)==0&&(j==3||fpart1[j-4]==' '||fpart1[j-4]=='('||fpart1[j-4]=='+'||fpart1[j-4]=='-'||fpart1[j-4]=='*'||fpart1[j-4]=='/'||fpart1[j-4]=='^'||fpart1[j-4]==','))fnmulti=2;
    if (j>=4&&strncmp(fpart1+j-4,"id0v",4)==0&&(j==4||fpart1[j-5]==' '||fpart1[j-5]=='('||fpart1[j-5]=='+'||fpart1[j-5]=='-'||fpart1[j-5]=='*'||fpart1[j-5]=='/'||fpart1[j-5]=='^'||fpart1[j-5]==','))fnmulti=3;
    if (fnmulti>0) {
      char abuf[TABREADLINE];
      offset_t argres[16];
      int nargs=0,k,adepth=0,a0=0,fl,anpow,anmul,anplu;
      fpart1[j-((fnmulti==3)?4:3)]='\0';
      for (k=0;; k++) {
        if(fpart2[k]=='{')adepth++;
        if(fpart2[k]=='}')adepth--;
        if((fpart2[k]==','&&adepth==0)||fpart2[k]=='\0') {
          if(nargs>=16) {
            errmsg("Error: too many arguments in an intrinsic function call: %s\n",fomulain);
            return 0;
          }
          fl=k-a0;
          if(fl<=0) {
            errmsg("Error: empty argument in an intrinsic function call: %s\n",fomulain);
            return 0;
          }
          strncpy(abuf,fpart2+a0,fl);
          abuf[fl]='\0';
          anpow=str_count_char(abuf, '^');
          if (anpow>0) {
            if(!formula_compile_pow(abuf,sets,anpow,arggrp,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
          }
          anmul=str_count_char(abuf, '*')+str_count_char(abuf, '/');
          if (anmul>0) {
            if(!formula_compile_muldiv(abuf,sets,anmul,arggrp,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
          }
          anplu=str_count_char(abuf, '+')+str_count_char(abuf, '-');
          if (anplu>0) {
            if(!formula_compile_addsub(abuf,sets,anplu,arggrp,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
          }
          arggrp++;
          formula_bind_operand(abuf,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
          ops[*nops].Oper=OP_LOAD;
          ops[*nops].TmpVarName[0]='\0';
          *nops=*nops+1;
          argres[nargs++]=*nops-1;
          a0=k+1;
          if(fpart2[k]=='\0')break;
        }
      }
      if(fnmulti==3&&nargs!=2) {
        errmsg("Error: ID0V takes exactly 2 arguments: %s\n",fomulain);
        return 0;
      }
      if(nargs<2) {
        errmsg("Error: %s takes at least 2 arguments: %s\n",(fnmulti==1)?"MAX":"MIN",fomulain);
        return 0;
      }
      for (k=1; k<nargs; k++) {
        ops[*nops].Oper=(fnmulti==1)?OP_MAXF:(fnmulti==2)?OP_MINF:OP_ID0VF;
        ops[*nops].Var1Type=OT_TEMP;
        ops[*nops].Var1BegAdd=(k==1)?argres[0]:(*nops-1);
        ops[*nops].Var2Type=OT_TEMP;
        ops[*nops].Var2BegAdd=argres[k];
        ops[*nops].TmpVarName[0]='\0';
        *nops=*nops+1;
      }
    } else {
    npow=str_count_char(fpart2, '^');
    if (npow>0) {
      if(!formula_compile_pow(fpart2,sets,npow,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
    }
    nmul=str_count_char(fpart2, '*');
    ndiv=str_count_char(fpart2, '/');
    nmul=nmul+ndiv;
    if (nmul>0) {
      if(!formula_compile_muldiv(fpart2,sets,nmul,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
    }
    nplu=str_count_char(fpart2, '+');
    nmin=str_count_char(fpart2, '-');
    nplu=nplu+nmin;
    if (nplu>0) {
      if(!formula_compile_addsub(fpart2,sets,nplu,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
    }
      if(strpbrk(fpart2,"=<>")==NULL){
      formula_bind_operand(fpart2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
      ops[*nops].Oper=OP_LOAD;
      ops[*nops].TmpVarName[0]='\0';
      *nops=*nops+1;
      }
    }
    sprintf(interchar1, "%d", i);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"gen_par000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"gen_par00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"gen_par0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"gen_par");
    }
    if (i>10000) {
      strcat(interchar,"gen_par");
      printf("Warning: too many nested parentheses in formula\n");
    }
    strcat(interchar,interchar1);
    strcpy(ops[*nops].TmpVarName,interchar);
    ops[*nops].Oper=OP_LOAD;
    ops[*nops].Var1BegAdd=*nops-1;
    j=strlen(fpart1);
    ops[*nops].Var1Type=OT_TEMP;
    if (j>3) if (fpart1[j-1]=='1'&&fpart1[j-2]=='0'&&tolower((int)fpart1[j-3])=='d'&&tolower((int)fpart1[j-4])=='i') {
        ops[*nops].Var1Type=OT_TEMP_ID01;
        fpart1[j-4]='\0';
      }
    if (j>3) if (fpart1[j-1]=='e'&&fpart1[j-2]=='g'&&fpart1[j-3]=='o'&&fpart1[j-4]=='l') {
        ops[*nops].Var1Type=OT_TEMP_LOG;
        fpart1[j-4]='\0';
      }
    if (j==3) if (fpart1[j-1]=='s'&&fpart1[j-2]=='b'&&fpart1[j-3]=='a') {
        ops[*nops].Var1Type=OT_TEMP_ABS;
        fpart1[j-3]='\0';
      }
    if (j>3) if (fpart1[j-1]=='s'&&fpart1[j-2]=='b'&&fpart1[j-3]=='a') if(j==3||fpart1[j-4]==' '||fpart1[j-4]=='('||fpart1[j-4]=='+'||fpart1[j-4]=='-'||fpart1[j-4]=='*'||fpart1[j-4]=='/'||fpart1[j-4]=='^'||fpart1[j-4]==',') {
          ops[*nops].Var1Type=OT_TEMP_ABS;
          fpart1[j-3]='\0';
        }
    /* intrinsic functions EXP/SQRT/LOG10/ROUND (manual 11.5, parity
       plan 3.1) -- same reversed-name recognition as loge/abs/id01,
       with the delimiter boundary check so names merely ending in a
       function word cannot trigger */
    if (j>=3) if (fpart1[j-1]=='p'&&fpart1[j-2]=='x'&&fpart1[j-3]=='e') if(j==3||fpart1[j-4]==' '||fpart1[j-4]=='('||fpart1[j-4]=='+'||fpart1[j-4]=='-'||fpart1[j-4]=='*'||fpart1[j-4]=='/'||fpart1[j-4]=='^'||fpart1[j-4]==',') {
          ops[*nops].Var1Type=OT_TEMP_EXP;
          fpart1[j-3]='\0';
        }
    if (j>=4) if (fpart1[j-1]=='t'&&fpart1[j-2]=='r'&&fpart1[j-3]=='q'&&fpart1[j-4]=='s') if(j==4||fpart1[j-5]==' '||fpart1[j-5]=='('||fpart1[j-5]=='+'||fpart1[j-5]=='-'||fpart1[j-5]=='*'||fpart1[j-5]=='/'||fpart1[j-5]=='^'||fpart1[j-5]==',') {
          ops[*nops].Var1Type=OT_TEMP_SQRT;
          fpart1[j-4]='\0';
        }
    if (j>=5) if (fpart1[j-1]=='0'&&fpart1[j-2]=='1'&&fpart1[j-3]=='g'&&fpart1[j-4]=='o'&&fpart1[j-5]=='l') if(j==5||fpart1[j-6]==' '||fpart1[j-6]=='('||fpart1[j-6]=='+'||fpart1[j-6]=='-'||fpart1[j-6]=='*'||fpart1[j-6]=='/'||fpart1[j-6]=='^'||fpart1[j-6]==',') {
          ops[*nops].Var1Type=OT_TEMP_LOG10;
          fpart1[j-5]='\0';
        }
    if (j>=5) if (fpart1[j-1]=='d'&&fpart1[j-2]=='n'&&fpart1[j-3]=='u'&&fpart1[j-4]=='o'&&fpart1[j-5]=='r') if(j==5||fpart1[j-6]==' '||fpart1[j-6]=='('||fpart1[j-6]=='+'||fpart1[j-6]=='-'||fpart1[j-6]=='*'||fpart1[j-6]=='/'||fpart1[j-6]=='^'||fpart1[j-6]==',') {
          ops[*nops].Var1Type=OT_TEMP_ROUND;
          fpart1[j-5]='\0';
        }
    /* trunc0 (toward zero) / truncb (toward -inf) -- plan 3.1 */
    if (j>=6) if (strncmp(fpart1+j-6,"trunc0",6)==0) if(j==6||fpart1[j-7]==' '||fpart1[j-7]=='('||fpart1[j-7]=='+'||fpart1[j-7]=='-'||fpart1[j-7]=='*'||fpart1[j-7]=='/'||fpart1[j-7]=='^'||fpart1[j-7]==',') {
          ops[*nops].Var1Type=OT_TEMP_TRUNC0;
          fpart1[j-6]='\0';
        }
    if (j>=6) if (strncmp(fpart1+j-6,"truncb",6)==0) if(j==6||fpart1[j-7]==' '||fpart1[j-7]=='('||fpart1[j-7]=='+'||fpart1[j-7]=='-'||fpart1[j-7]=='*'||fpart1[j-7]=='/'||fpart1[j-7]=='^'||fpart1[j-7]==',') {
          ops[*nops].Var1Type=OT_TEMP_TRUNCB;
          fpart1[j-6]='\0';
        }
    if (j==2) if (fpart1[j-1]=='f'&&fpart1[j-2]=='i') {
        if(!formula_compile_if(fpart2,sets,2,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
        fpart1[j-2]='\0';
      }
    if (j>2) if (fpart1[j-1]=='f'&&fpart1[j-2]=='i') if(j==2||fpart1[j-3]==' '||fpart1[j-3]=='('||fpart1[j-3]=='+'||fpart1[j-3]=='-'||fpart1[j-3]=='*'||fpart1[j-3]=='/'||fpart1[j-3]=='^'||fpart1[j-3]==',') {
          if(!formula_compile_if(fpart2,sets,2,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim))return 0;
          fpart1[j-2]='\0';
        }
    *nops=*nops+1;
    if(strlen(fpart1)+strlen(interchar)+strlen(fpart3)>=TABREADLINE){
      errmsg("Error: formula too long to compile (exceeds %d chars): %s\n",TABREADLINE,fomulain);
      return 0;
    }
    strcat(fpart1, interchar);
    strcat(fpart1, fpart3);
    strcpy(fomulain,fpart1);
  }
  return 1;
}

/* one declared-range check over a coefficient's elements (manual
   25.4/11.6.7): prints per the CMF mode (Warning/Error) on the first
   violating element; returns 1 on violation. Called once per bound
   slot -- one lower + one upper per declaration (audit A9). */
static int coef_range_check(array_def *coefs,offset_t index,offset_t offset,offset_t varsize,elem_value *elem_vals,int glmode,int gltype,store_real glval,int updated) {
  offset_t l;
  int glviol=0;
  const char *msg="";
  /* "updated value" marks a post-update check so the front end can tell
     a step overshoot (violations, then a singular LU) from a
     structurally deficient closure (singular at the first LU) */
  const char *phase=updated?"an updated value":"a value";
  if(gltype==BT_GE)msg="below its declared lower bound";
  if(gltype==BT_GT)msg="at or below its declared strict lower bound";
  if(gltype==BT_LE)msg="above its declared upper bound";
  if(gltype==BT_LT)msg="at or above its declared strict upper bound";
  #pragma omp parallel for private(l) shared(elem_vals,coefs,offset,index,glviol,gltype,glval,glmode,msg,phase,varsize)
  for (l=0; l<varsize; l++) {
    int bad=0;
    if(gltype==BT_GE&&elem_vals[offset+l].value<glval)bad=1;
    if(gltype==BT_GT&&elem_vals[offset+l].value<=glval)bad=1;
    if(gltype==BT_LE&&elem_vals[offset+l].value>glval)bad=1;
    if(gltype==BT_LT&&elem_vals[offset+l].value>=glval)bad=1;
    if(bad&&!glviol) {
      printf("%s: coefficient %s has %s %s %f\n",glmode==2?"Error":"Warning",coefs[index].cofname,phase,msg,(double)glval);
      glviol=1;
    }
  }
  return glviol;
}

/* active dual-class zerodivide state for the CURRENT statement
   (plan A1): captured from the scanner position by the formula/
   assertion executors, disabled during update and equation
   evaluation (manual 10.11.1: zerodivide never applies there --
   those keep the legacy single default until their own pass) */
static zdiv_state zdiv_active = { 0, 0, 1, 0 };
static int zdiv_enabled = 0;

void zdiv_capture(void) {
  zdiv_active=teems_zdiv_scan;
  zdiv_enabled=teems_gpzerodivide;
}

void zdiv_disable(void) {
  zdiv_enabled=0;
}

/* per-tuple element offset of one operand: each quantifier position
   contributes stride * (element index -- routed through superset_pos
   when the operand's set is wider than the frame's -- plus leadlag).
   The single authority for operand addressing; previously ~30 verbatim
   copies across the eval switch (refactored for the mapping mode,
   docs/mapping_complementarity_design.md M2). Sum operands arrive with
   SupSet/leadlag zeroed, so the full form is exact for them too. */
static inline offset_t dims_route(offset_t ADims, int SupSet, int SSIndx, int leadlag, int MapId, int MapDomSS, const quantifier *q, const set_def *sets, const set_element *set_elems) {
  if(MapId>0) {
    /* mapped argument (manual 11.9.4): the domain element routes
       through the mapping's codomain-position table -- from a subset
       of the domain via its domain position, and on into a superset
       argument set via the codomain element's position there */
    const map_def *md=&teems_maps[MapId-1];
    offset_t di=q->indx;
    offset_t v;
    if(MapDomSS>0) di=set_elems[sets[q->setid].offset+di].superset_pos[MapDomSS];
    v=md->values[di+leadlag];
    if(SupSet==1) v=set_elems[sets[md->toset].offset+v].superset_pos[SSIndx];
    return ADims*v;
  } else if(SupSet==1) {
    return ADims*(set_elems[sets[q->setid].offset+q->indx].superset_pos[SSIndx]+leadlag);
  }
  return ADims*(q->indx+leadlag);
}

static inline offset_t dims_offset(const dim_addr *D, const quantifier *arSet, dim_t fdim, const set_def *sets, const set_element *set_elems) {
  offset_t l=0;
  dim_t j;
  for (j=0; j<fdim; j++) {
    if(D[j].MapId>0) {
      l+=dims_route(D[j].ADims,D[j].SupSet,D[j].SSIndx,D[j].leadlag,D[j].MapId,D[j].MapDomSS,&arSet[j],sets,set_elems);
    } else if(D[j].SupSet==1) {
      l+=D[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[D[j].SSIndx]+D[j].leadlag);
    } else {
      l+=D[j].ADims*(arSet[j].indx+D[j].leadlag);
    }
    /* a repeated index whose positions route differently (dim_bind_pos) */
    if(D[j].Rep) l+=dims_route(D[j].ADims2,D[j].SupSet2,D[j].SSIndx2,D[j].leadlag2,D[j].MapId2,D[j].MapDomSS2,&arSet[j],sets,set_elems);
  }
  return l;
}

solve_real formula_eval(elem_value *record,set_def *sets,set_element *set_elems,sum_value *sum_vals,formula_op *ops,int nops,quantifier *arSet,dim_t fdim, solve_real zerodivide) {
  int i;
  offset_t l=0,l1=0;
  solve_real eval1=0,eval2=0,eval3=0;
  for (i=0; i<nops; i++) {
    switch(ops[i].Oper) {
    case OP_LOAD:
      if (ops[i].Var1Type==OT_ARRAY) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        ops[i].TmpVarVal=record[ops[i].Var1BegAdd+l].value;
        break;
      }
      if (ops[i].Var1Type==OT_SUM) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        ops[i].TmpVarVal=sum_vals[ops[i].Var1BegAdd+l].value;
        break;
      }
      if (ops[i].Var1Type==OT_TEMP) {
        ops[i].TmpVarVal=ops[ops[i].Var1BegAdd].TmpVarVal;
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_ID01) {
        if(ops[ops[i].Var1BegAdd].TmpVarVal==0)ops[i].TmpVarVal=1;
        else ops[i].TmpVarVal=ops[ops[i].Var1BegAdd].TmpVarVal;
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_ABS) {
        if(ops[ops[i].Var1BegAdd].TmpVarVal<0)ops[i].TmpVarVal=-ops[ops[i].Var1BegAdd].TmpVarVal;
        else ops[i].TmpVarVal=ops[ops[i].Var1BegAdd].TmpVarVal;
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_LOG) {
        ops[i].TmpVarVal=log(ops[ops[i].Var1BegAdd].TmpVarVal);
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_EXP) {
        ops[i].TmpVarVal=exp(ops[ops[i].Var1BegAdd].TmpVarVal);
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_SQRT) {
        ops[i].TmpVarVal=sqrt(ops[ops[i].Var1BegAdd].TmpVarVal);
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_LOG10) {
        ops[i].TmpVarVal=log10(ops[ops[i].Var1BegAdd].TmpVarVal);
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_ROUND) {
        ops[i].TmpVarVal=round(ops[ops[i].Var1BegAdd].TmpVarVal);
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_TRUNC0) {
        ops[i].TmpVarVal=trunc(ops[ops[i].Var1BegAdd].TmpVarVal);
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_TRUNCB) {
        ops[i].TmpVarVal=floor(ops[ops[i].Var1BegAdd].TmpVarVal);
        break;
      }
      if (ops[i].Var1Type==OT_CONST) {
        ops[i].TmpVarVal=ops[i].Var1Val;
        break;
      }
      if (ops[i].Var1Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        ops[i].TmpVarVal=record[ops[i].Var1BegAdd+l].substep_base;
        break;
      }
      if (ops[i].Var1Type==OT_POS) {
        /* $POS: 1-based position of a quantifier index, optionally
           routed through a mapping and lifted into a superset */
        offset_t pidx=arSet[ops[i].Var1BegAdd].indx;
        if (ops[i].Var1Dims[0].MapId>0) pidx=teems_maps[ops[i].Var1Dims[0].MapId-1].values[pidx];
        if (ops[i].Var1Dims[0].SupSet==1) pidx=set_elems[sets[ops[i].Var1Dims[0].ADims].offset+pidx].superset_pos[ops[i].Var1Dims[0].SSIndx];
        ops[i].TmpVarVal=(store_real)(pidx+1);
        break;
      }
      break;
    case OP_MUL:
      if(ops[i].Var1Type<3||ops[i].Var1Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
      }
      if(ops[i].Var2Type<3||ops[i].Var2Type==OT_CHANGE) {
        l1=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
      }
      if(ops[i].Var1Type==OT_ARRAY) eval1=record[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_SUM) eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_TEMP) eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      if(ops[i].Var1Type==OT_CONST) eval1=ops[i].Var1Val;
      if(ops[i].Var1Type==OT_CHANGE) eval1=record[ops[i].Var1BegAdd+l].substep_base;
      if(ops[i].Var2Type==OT_ARRAY) eval2=record[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_SUM) eval2=sum_vals[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_TEMP) eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      if(ops[i].Var2Type==OT_CONST) eval2=ops[i].Var2Val;
      if(ops[i].Var2Type==OT_CHANGE) eval2=record[ops[i].Var2BegAdd+l1].substep_base;
      ops[i].TmpVarVal=eval1*eval2;
      break;
    case OP_DIV:
      if(ops[i].Var1Type<3||ops[i].Var1Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
      }
      if(ops[i].Var2Type<3||ops[i].Var2Type==OT_CHANGE) {
        l1=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
      }
      if(ops[i].Var1Type==OT_ARRAY) eval1=record[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_SUM) eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_TEMP) eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      if(ops[i].Var1Type==OT_CONST) eval1=ops[i].Var1Val;
      if(ops[i].Var1Type==OT_CHANGE) eval1=record[ops[i].Var1BegAdd+l].substep_base;
      if(ops[i].Var2Type==OT_ARRAY) eval2=record[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_SUM) eval2=sum_vals[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_TEMP) eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      if(ops[i].Var2Type==OT_CONST) eval2=ops[i].Var2Val;
      if(ops[i].Var2Type==OT_CHANGE) eval2=record[ops[i].Var2BegAdd+l1].substep_base;
      if(eval2==0) {
        if(zdiv_enabled) {
          if(eval1==0) {
            if(zdiv_active.zbz_on)ops[i].TmpVarVal=zdiv_active.zbz_val;
            else {
              errmsg("Error: zero divided by zero in a formula while Zerodivide (zero_by_zero) is off\n");
              MPI_Abort(PETSC_COMM_WORLD,1);
            }
          } else {
            if(zdiv_active.nbz_on)ops[i].TmpVarVal=zdiv_active.nbz_val;
            else {
              errmsg("Error: division by zero in a formula; Zerodivide (nonzero_by_zero) is off (GEMPACK default) -- set a default or guard with ID01\n");
              MPI_Abort(PETSC_COMM_WORLD,1);
            }
          }
        } else {
          ops[i].TmpVarVal=zerodivide;
        }
      } else {
        ops[i].TmpVarVal=eval1/eval2;
      }
      break;
    case OP_ADD:
      if(ops[i].Var1Type==OT_ARRAY) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        eval1=record[ops[i].Var1BegAdd+l].value;
      }
      if(ops[i].Var1Type==OT_SUM) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      }
      if (ops[i].Var1Type==OT_TEMP) {
        eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      }
      if (ops[i].Var1Type==OT_CONST) {
        eval1=ops[i].Var1Val;
      }
      if(ops[i].Var1Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        eval1=record[ops[i].Var1BegAdd+l].substep_base;
      }
      if(ops[i].Var2Type==OT_ARRAY) {
        l=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
        eval2=record[ops[i].Var2BegAdd+l].value;
      }
      if(ops[i].Var2Type==OT_SUM) {
        l=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
        eval2=sum_vals[ops[i].Var2BegAdd+l].value;
      }
      if (ops[i].Var2Type==OT_TEMP) {
        eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      }
      if (ops[i].Var2Type==OT_CONST) {
        eval2=ops[i].Var2Val;
      }
      if(ops[i].Var2Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
        eval2=record[ops[i].Var2BegAdd+l].substep_base;
      }
      ops[i].TmpVarVal=eval1+eval2;
      break;
    case OP_SUB:
      if(ops[i].Var1Type==OT_ARRAY) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        eval1=record[ops[i].Var1BegAdd+l].value;
      }
      if(ops[i].Var1Type==OT_SUM) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      }
      if (ops[i].Var1Type==OT_TEMP) {
        eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      }
      if (ops[i].Var1Type==OT_CONST) {
        eval1=ops[i].Var1Val;
      }
      if(ops[i].Var1Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        eval1=record[ops[i].Var1BegAdd+l].substep_base;
      }
      if(ops[i].Var2Type==OT_ARRAY) {
        l=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
        eval2=record[ops[i].Var2BegAdd+l].value;
      }
      if(ops[i].Var2Type==OT_SUM) {
        l=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
        eval2=sum_vals[ops[i].Var2BegAdd+l].value;
      }
      if (ops[i].Var2Type==OT_TEMP) {
        eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      }
      if (ops[i].Var2Type==OT_CONST) {
        eval2=ops[i].Var2Val;
      }
      if(ops[i].Var2Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
        eval2=record[ops[i].Var2BegAdd+l].substep_base;
      }
      ops[i].TmpVarVal=eval1-eval2;
      break;
    case OP_POW:
      if(ops[i].Var1Type<3||ops[i].Var1Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
      }
      if(ops[i].Var2Type<3||ops[i].Var2Type==OT_CHANGE) {
        l1=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
      }
      if(ops[i].Var1Type==OT_ARRAY) eval1=record[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_SUM) eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_TEMP) eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      if(ops[i].Var1Type==OT_CONST) eval1=ops[i].Var1Val;
      if(ops[i].Var1Type==OT_CHANGE) eval1=record[ops[i].Var1BegAdd+l].substep_base;
      if(ops[i].Var2Type==OT_ARRAY) eval2=record[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_SUM) eval2=sum_vals[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_TEMP) eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      if(ops[i].Var2Type==OT_CONST) eval2=ops[i].Var2Val;
      if(ops[i].Var2Type==OT_CHANGE) eval2=record[ops[i].Var2BegAdd+l1].substep_base;
      if(eval1==0&&eval2<0) {
        ops[i].TmpVarVal=zerodivide;
      } else {
        if(eval1<0&&eval2-floor(eval2)!=0)errmsg("Error: fractional power of a negative number in formula evaluation\n");
        ops[i].TmpVarVal=pow(eval1,eval2);
      }
      break;
    case OP_MAXF:
    case OP_MINF:
    case OP_ID0VF:
      /* multi-arg intrinsic folds: operands are always compiled temps */
      eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      if(ops[i].Oper==OP_MAXF)ops[i].TmpVarVal=(eval1>eval2)?eval1:eval2;
      else if(ops[i].Oper==OP_MINF)ops[i].TmpVarVal=(eval1<eval2)?eval1:eval2;
      else ops[i].TmpVarVal=(eval1!=0)?eval1:eval2;
      break;
    default:
      if(ops[i].Var1Type==OT_ARRAY) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        eval1=record[ops[i].Var1BegAdd+l].value;
      }
      if(ops[i].Var1Type==OT_SUM) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      }
      if (ops[i].Var1Type==OT_TEMP) {
        eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      }
      if (ops[i].Var1Type==OT_CONST) {
        eval1=ops[i].Var1Val;
      }
      if(ops[i].Var1Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var1Dims,arSet,fdim,sets,set_elems);
        eval1=record[ops[i].Var1BegAdd+l].substep_base;
      }
      if(ops[i].Var2Type==OT_ARRAY) {
        l=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
        eval2=record[ops[i].Var2BegAdd+l].value;
      }
      if(ops[i].Var2Type==OT_SUM) {
        l=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
        eval2=sum_vals[ops[i].Var2BegAdd+l].value;
      }
      if (ops[i].Var2Type==OT_TEMP) {
        eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      }
      if (ops[i].Var2Type==OT_CONST) {
        eval2=ops[i].Var2Val;
      }
      if(ops[i].Var2Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var2Dims,arSet,fdim,sets,set_elems);
        eval2=record[ops[i].Var2BegAdd+l].substep_base;
      }

      if(ops[i].Var3Type==OT_ARRAY) {
        l=dims_offset(ops[i].Var3Dims,arSet,fdim,sets,set_elems);
        eval3=record[ops[i].Var3BegAdd+l].value;
      }
      if(ops[i].Var3Type==OT_SUM) {
        l=dims_offset(ops[i].Var3Dims,arSet,fdim,sets,set_elems);
        eval3=sum_vals[ops[i].Var3BegAdd+l].value;
      }
      if (ops[i].Var3Type==OT_TEMP) {
        eval3=ops[ops[i].Var3BegAdd].TmpVarVal;
      }
      if (ops[i].Var3Type==OT_CONST) {
        eval3=ops[i].Var3Val;
      }
      if(ops[i].Var3Type==OT_CHANGE) {
        l=dims_offset(ops[i].Var3Dims,arSet,fdim,sets,set_elems);
        eval3=record[ops[i].Var3BegAdd+l].substep_base;
      }
      if(ops[i].Oper==OP_IF_EQ){ if(eval1==eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0; }
      if(ops[i].Oper==OP_IF_GT){ if(eval1>eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0; }
      if(ops[i].Oper==OP_IF_LT){ if(eval1<eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0; }
      if(ops[i].Oper==OP_IF_NE){ if(eval1!=eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0; }
      if(ops[i].Oper==OP_IF_LE){ if(eval1<=eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0; }
      if(ops[i].Oper==OP_IF_GE){ if(eval1>=eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0; }
      break;
    }
  }
  return ops[i-1].TmpVarVal;
}




int formula_compile_pow(char *fomulain, set_def *sets,int npow,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim) {
  int i,i1,ibar=0,index,j,j1,i5,p1;
  char *p=NULL;//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE],var1[TABREADLINE],var2[TABREADLINE],interchar[TABREADLINE],interchar1[TABREADLINE];

  for (i=1; i<npow+1; i++) {
    p=strchr(fomulain,'^');
    index=p-fomulain;
    strncpy(fpart1, fomulain, index);
    fpart1[index] = '\0';
    strcpy(fpart3,fomulain+index+1);

    i1=-1;
    p=str_rfind_any(fpart1,"^*/+-=<>");
    if (p!=NULL) {
      i1=p-fpart1;
    }
    i5=-1;
    j=strlen(fpart1)-1;
    while(j>-1){
      if(fpart1[j]==','){
        p1=-1;
        for(j1=j-1;j1>-1;j1--){
          if(fpart1[j1]=='{')p1=j1;
          if(fpart1[j1]=='}')break;
        }
        if(p1==-1)break;
        else j=p1;
      }
      j--;
    }
    if(j>0)i5=j;
    index=i1;
    if (index<i5) {
      index=i5;
    }
    strcpy(var1, fpart1+index+1);
    strncpy(fpart1,fomulain,index+1);
    fpart1[index+1] = '\0';

    ibar=0;
    while (fpart3[ibar] != '\0') {
      ibar++;
    }
    p=strpbrk(fpart3,"^*/+-=<>");
    if(p==NULL)index=ibar;
    else index=p-fpart3;
    strncpy(var2, fpart3, index);
    var2[index] = '\0';
    strcpy(fpart2,fpart3+index);

    formula_bind_operand(var1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
    formula_bind_operand(var2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);

    sprintf(interchar1, "%d", ipar);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"gen_pow");
      strcat(interchar,interchar1);
      strcat(interchar,"000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"gen_pow");
      strcat(interchar,interchar1);
      strcat(interchar,"00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"gen_pow");
      strcat(interchar,interchar1);
      strcat(interchar,"0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"gen_pow");
      strcat(interchar,interchar1);
    }
    sprintf(interchar1, "%d", i);
    strcat(interchar,interchar1);
    strcpy(ops[*nops].TmpVarName,interchar);
    ops[*nops].Oper=OP_POW;
    *nops=*nops+1;
    if(strlen(fpart1)+strlen(interchar)+strlen(fpart2)>=TABREADLINE){
      errmsg("Error: formula too long to compile (exceeds %d chars): %s\n",TABREADLINE,fomulain);
      return 0;
    }
    strcat(fpart1, interchar);
    strcat(fpart1, fpart2);
    strcpy(fomulain,fpart1);
  }
  return 1;
}
int formula_compile_muldiv(char *fomulain, set_def *sets,int nmul,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim) {
  int i,i1,i5,ibar=0,index,j,j1,p1;
  char *p=NULL;//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE],var1[TABREADLINE],var2[TABREADLINE],interchar[TABREADLINE],interchar1[TABREADLINE];

  for (i=1; i<nmul+1; i++) {
    index=0;
    p=strpbrk(fomulain,"*/");
    if (*p=='/') {
      ops[*nops].Oper=OP_DIV;
    } else {
      ops[*nops].Oper=OP_MUL;
    }
    index=p-fomulain;

    strncpy(fpart1, fomulain, index);
    fpart1[index] = '\0';
    strcpy(fpart3,fomulain+index+1);
    i1=-1;
    p=str_rfind_any(fpart1,"*/+-=<>");
    if (p!=NULL) {
      i1=p-fpart1;
    }
    i5=-1;
    j=strlen(fpart1)-1;
    while(j>-1){
      if(fpart1[j]==','){
        p1=-1;
        for(j1=j-1;j1>-1;j1--){
          if(fpart1[j1]=='{')p1=j1;
          if(fpart1[j1]=='}')break;
        }
        if(p1==-1)break;
        else j=p1;
      }
      j--;
    }
    if(j>0)i5=j;
    index=i1;
    if (index<i5) {
      index=i5;
    }
    strcpy(var1, fpart1+index+1);
    strncpy(fpart1,fomulain,index+1);
    fpart1[index+1] = '\0';

    ibar=0;
    while (fpart3[ibar] != '\0') {
      ibar++;
    }
    p=strpbrk(fpart3,"*/+-=<>");
    if(p==NULL)index=ibar;
    else index=p-fpart3;
    strncpy(var2, fpart3, index);
    var2[index] = '\0';
    strcpy(fpart2,fpart3+index);

    formula_bind_operand(var1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
    formula_bind_operand(var2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);

    sprintf(interchar1, "%d", ipar);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"gen_mul");
      strcat(interchar,interchar1);
      strcat(interchar,"000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"gen_mul");
      strcat(interchar,interchar1);
      strcat(interchar,"00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"gen_mul");
      strcat(interchar,interchar1);
      strcat(interchar,"0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"gen_mul");
      strcat(interchar,interchar1);
    }
    sprintf(interchar1, "%d", i);
    strcat(interchar,interchar1);
    strcpy(ops[*nops].TmpVarName,interchar);
    *nops=*nops+1;
    if(strlen(fpart1)+strlen(interchar)+strlen(fpart2)>=TABREADLINE){
      errmsg("Error: formula too long to compile (exceeds %d chars): %s\n",TABREADLINE,fomulain);
      return 0;
    }
    strcat(fpart1, interchar);
    strcat(fpart1, fpart2);
    strcpy(fomulain,fpart1);
  }
  return 1;
}

int formula_compile_addsub(char *fomulain, set_def *sets,int nplu,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim) {
  int i,i1,i5,ibar=0,index,j,j1,p1;
  char *p=NULL;//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE],var1[TABREADLINE],var2[TABREADLINE],interchar[TABREADLINE],interchar1[TABREADLINE];

  for (i=1; i<nplu+1; i++) {
    index=0;
    p=strpbrk(fomulain,"+-");
    if (*p=='+') {
      ops[*nops].Oper=OP_ADD;
    } else {
      ops[*nops].Oper=OP_SUB;
    }
    index=p-fomulain;

    strncpy(fpart1, fomulain, index);
    fpart1[index] = '\0';
    strcpy(fpart3,fomulain+index+1);
    i1=-1;
    p=str_rfind_any(fpart1,"+-=<>");
    if (p!=NULL) {
      i1=p-fpart1;
    }
    i5=-1;
    j=strlen(fpart1)-1;
    while(j>-1){
      if(fpart1[j]==','){
        p1=-1;
        for(j1=j-1;j1>-1;j1--){
          if(fpart1[j1]=='{')p1=j1;
          if(fpart1[j1]=='}')break;
        }
        if(p1==-1)break;
        else j=p1;
      }
      j--;
    }
    if(j>0)i5=j;
    index=i1;
    if (index<i5) {
      index=i5;
    }
    strcpy(var1, fpart1+index+1);
    strncpy(fpart1,fomulain,index+1);
    fpart1[index+1] = '\0';

    ibar=0;
    while (fpart3[ibar] != '\0') {
      ibar++;
    }
    p=strpbrk(fpart3,"+-=<>");
    if(p==NULL)index=ibar;
    else index=p-fpart3;
    strncpy(var2, fpart3, index);
    var2[index] = '\0';
    strcpy(fpart2,fpart3+index);

    if(i==1&&var1[0]=='\0'){
        ops[*nops].Var1Type=OT_CONST;
        ops[*nops].Var1Val=0;
    }
    else formula_bind_operand(var1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
    if(i==1&&var2[0]=='\0'){
        ops[*nops].Var2Type=OT_CONST;
        ops[*nops].Var2Val=0;
    }
    else formula_bind_operand(var2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);

    sprintf(interchar1, "%d", ipar);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"gen_add");
      strcat(interchar,interchar1);
      strcat(interchar,"000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"gen_add");
      strcat(interchar,interchar1);
      strcat(interchar,"00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"gen_add");
      strcat(interchar,interchar1);
      strcat(interchar,"0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"gen_add");
      strcat(interchar,interchar1);
    }
    sprintf(interchar1, "%d", i);
    strcat(interchar,interchar1);
    strcpy(ops[*nops].TmpVarName,interchar);
    *nops=*nops+1;
    if(strlen(fpart1)+strlen(interchar)+strlen(fpart2)>=TABREADLINE){
      errmsg("Error: formula too long to compile (exceeds %d chars): %s\n",TABREADLINE,fomulain);
      return 0;
    }
    strcat(fpart1, interchar);
    strcat(fpart1, fpart2);
    strcpy(fomulain,fpart1);
  }
  return 1;
}


int formula_compile_if(char *fomulain, set_def *sets,int nif,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim) {
  char *p=NULL,*p1,*p3,var1[NAMESIZE],var2[NAMESIZE],var3[NAMESIZE];
  int i,j1,j2,j3,l;//,varindex;
  p1=fomulain;
  p=strpbrk(p1,"=<>");
  if(p==NULL){
    errmsg("Error: malformed if() in formula (no comparison operator): %s\n",fomulain);
    return 0;
  }
  if(p-p1>=NAMESIZE){
    errmsg("Error: if() condition operand too long in formula: %s\n",fomulain);
    return 0;
  }
  strncpy(var1,p1,p-p1);
  var1[p-p1]='\0';
  p3=p+1;
  if(*p3=='='){
    if(*p=='<')ops[*nops].Oper=OP_IF_LE;
    if(*p=='>')ops[*nops].Oper=OP_IF_GE;
    p++;
    p++;
  }else if(*p3=='>'){
          ops[*nops].Oper=OP_IF_NE;
          p++;
          p++;
          }else{
          if(*p=='=')ops[*nops].Oper=OP_IF_EQ;
          if(*p=='>')ops[*nops].Oper=OP_IF_GT;
          if(*p=='<')ops[*nops].Oper=OP_IF_LT;
          p++;
        }
  l=strlen(p);
  j1=-1;j2=-1;j3=-1;
  for(i=0;i<l;i++){
    if(*(p+i)==',')j1=i;
    if(*(p+i)=='}')j2=i;
    if(*(p+i)=='{')j3=i;
    if(j3>-1){ if(j1>-1&&j1>j2)break;
    else if(j1>-1)break; }
  }
  if(j1<0){
    errmsg("Error: malformed if() in formula (missing comma before value): %s\n",fomulain);
    return 0;
  }
  if(j1>=NAMESIZE||strlen(p+j1+1)>=NAMESIZE){
    errmsg("Error: if() value operand too long in formula: %s\n",fomulain);
    return 0;
  }
  strncpy(var2,p,j1);
  var2[j1]='\0';
  strcpy(var3,p+j1+1);
  formula_bind_operand(var1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
  formula_bind_operand(var3,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);
  ops[*nops].Var3Type=ops[*nops].Var2Type;
  ops[*nops].Var3BegAdd=ops[*nops].Var2BegAdd;
  /* the whole address, repeat/mapping routing included */
  for(i=0;i<fdim;i++) ops[*nops].Var3Dims[i]=ops[*nops].Var2Dims[i];
  ops[*nops].Var3Val=ops[*nops].Var2Val;
  formula_bind_operand(var2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);
  return 1;
}

/* argument count of a condition reference's "a,b," list (the caller
   appends the trailing comma; a scalar reference leaves just ",") */
static int cond_arg_count(const char *argu) {
  int n=0;
  if(strcmp(argu,",")==0)return 0;
  for(; *argu!='\0'; argu++) if(*argu==',') n++;
  return n;
}

/* Store one mapping value from a Formula (manual 10.13.1/11.9.9):
   codomain position cod (0-based) for domain position dom.  A value
   may be re-assigned freely until the mapping has been used; after
   that a change is an error (11.9.9), a repeat is harmless (multi-step
   formula passes re-run the statement). */
static void mapping_store_value(dim_t mm, dim_t dom, dim_t cod, set_def *sets, set_element *set_elems) {
  map_def *md=&teems_maps[mm];
  dim_t n=sets[md->fromset].size;
  if (cod<0||cod>=sets[md->toset].size) {
    errmsg("Error: Formula assigns mapping %s a position outside its codomain set %s (%d elements): %d (manual 11.9.2)\n",md->mapname,sets[md->toset].setname,(int)sets[md->toset].size,(int)cod+1);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  if (md->assigned[dom]&&md->values[dom]!=cod&&md->used) {
    errmsg("Error: Formula changes the value of mapping %s for element %s after the mapping has been used (manual 11.9.9)\n",md->mapname,set_elems[sets[md->fromset].offset+dom].setele);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  md->values[dom]=cod;
  if (!md->assigned[dom]) { md->assigned[dom]=1; md->nassigned++; }
  if (md->nassigned==n&&!md->has_values) {
    md->has_values=true;
    if (mapping_check_onto(teems_maps,mm,sets,set_elems)==-1) MPI_Abort(PETSC_COMM_WORLD,1);
  }
}

/* element name -> position in set s, or -1 */
static dim_t set_element_pos(const char *el, dim_t s, set_def *sets, set_element *set_elems) {
  dim_t e;
  for (e=0; e<sets[s].size; e++) if (strcmp(el,set_elems[sets[s].offset+e].setele)==0) return e;
  return -1;
}

/* MAP("dom") = "cod" | <position> (manual 10.13.1): lhs is the
   space-stripped LHS text, rhs the RHS text (';' may trail) */
static void mapping_assign_literal(dim_t mm, char *lhs, char *rhs, set_def *sets, set_element *set_elems, int byele) {
  map_def *md=&teems_maps[mm];
  char el[NAMESIZE];
  char *p=strchr(lhs,'('),*q;
  int k=0;
  dim_t dom,cod;
  /* mapping_lower_calls has rewritten map("el") to map~"el" */
  if (p==NULL) p=strchr(lhs,MAPMARK);
  if (p==NULL||p[1]!='"') {
    errmsg("Error: unquantified Formula for mapping %s must name a domain element in quotes, e.g. %s(\"food\") (manual 10.13.1)\n",md->mapname,md->mapname);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  for (p+=2; *p!='"'&&*p!='\0'&&k<NAMESIZE-1; p++) el[k++]=tolower((int)*p);
  el[k]='\0';
  dom=set_element_pos(el,md->fromset,sets,set_elems);
  if (dom<0) {
    errmsg("Error: %s is not an element of set %s, the domain of mapping %s (manual 11.9.2)\n",el,sets[md->fromset].setname,md->mapname);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  while (*rhs==' ') rhs++;
  if ((q=strchr(rhs,';'))!=NULL) *q='\0';
  if (rhs[0]=='"') {
    k=0;
    for (p=rhs+1; *p!='"'&&*p!='\0'&&k<NAMESIZE-1; p++) el[k++]=tolower((int)*p);
    el[k]='\0';
    cod=set_element_pos(el,md->toset,sets,set_elems);
    if (cod<0) {
      errmsg("Error: %s is not an element of set %s, the codomain of mapping %s (manual 11.9.2)\n",el,sets[md->toset].setname,md->mapname);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  } else {
    if (byele) {
      errmsg("Error: Formula (by_elements) for mapping %s needs a quoted codomain element on the right-hand side (manual 10.13.1)\n",md->mapname);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    cod=(dim_t)atoi(rhs)-1;
  }
  mapping_store_value(mm,dom,cod,sets,set_elems);
}

/* (all,i,DOM) MAP(i) = <expr> (manual 10.13.1): the single LHS index
   ranges over the domain set or a subset of it; the RHS is an integer
   codomain position (typically a $POS shape) or, under (by_elements),
   a quoted codomain element.  Sums are already evaluated (sum_vals). */
static void mapping_assign_formula(dim_t mm, char *vname, char *rhs, int byele, set_def *sets, set_element *set_elems, sum_value *sum_vals, elem_value *elem_vals, quantifier *arSet, dim_t fdim, offset_t *dcountdim1, offset_t nloops, formula_op *ops, array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar, offset_t ncofele, sum_def *sum_cof, dim_t totalsum, solve_real zerodivide) {
  map_def *md=&teems_maps[mm];
  char arg[NAMESIZE],el[NAMESIZE];
  char *p=strchr(vname,'(');
  int k=0,sup=0,litdom=0,idxmode=0;
  dim_t l=-1,dom=-1,codlit=-1,nops=0,idxl=0,idxss=0,cpdss=0,cpmap=0;
  offset_t t,i4,i3;
  dim_t dcount;
  /* mapping_lower_calls has rewritten map(i) to map~i */
  if (p==NULL) p=strchr(vname,MAPMARK);
  if (p==NULL) {
    errmsg("Error: Formula for mapping %s has no argument (manual 10.13.1)\n",md->mapname);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  for (p++; *p!=')'&&*p!='}'&&*p!='\0'&&k<NAMESIZE-1; p++) arg[k++]=*p;
  arg[k]='\0';
  if (strchr(arg,',')!=NULL) {
    errmsg("Error: mapping %s takes one argument in a Formula (manual 10.13.1): %s\n",md->mapname,vname);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  if (arg[0]=='"') {
    k=0;
    for (p=arg+1; *p!='"'&&*p!='\0'&&k<NAMESIZE-1; p++) el[k++]=tolower((int)*p);
    el[k]='\0';
    dom=set_element_pos(el,md->fromset,sets,set_elems);
    if (dom<0) {
      errmsg("Error: %s is not an element of set %s, the domain of mapping %s (manual 11.9.2)\n",el,sets[md->fromset].setname,md->mapname);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    litdom=1;
  } else {
    for (l=0; l<fdim-1; l++) if (strcmp(arg,arSet[l].index_name)==0) break;
    if (l==fdim-1) {
      errmsg("Error: the argument %s of mapping %s is not one of the Formula's quantifier indices\n",arg,md->mapname);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    if ((dim_t)arSet[l].setid!=md->fromset) {
      for (sup=1; sup<MAXSUPSET; sup++) if (sets[arSet[l].setid].subsetid[sup]==md->fromset) break;
      if (sup==MAXSUPSET) {
        errmsg("Error: index %s of the Formula for mapping %s does not range over its domain set %s or a subset of it (manual 11.9.2)\n",arg,md->mapname,sets[md->fromset].setname);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
    }
  }
  while (*rhs==' ') rhs++;
  if (rhs[0]=='"') {
    k=0;
    for (p=rhs+1; *p!='"'&&*p!='\0'&&k<NAMESIZE-1; p++) el[k++]=tolower((int)*p);
    el[k]='\0';
    codlit=set_element_pos(el,md->toset,sets,set_elems);
    if (codlit<0) {
      errmsg("Error: %s is not an element of set %s, the codomain of mapping %s (manual 11.9.2)\n",el,sets[md->toset].setname,md->mapname);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  } else {
    /* identity and copy forms (manual 10.13.1; potential-models vetting
       G-S9/ORANI G17): the RHS is a quantifier index whose set is the
       codomain or a declared subset of it (MAP(c) = c, AGG(s) = s), or
       another mapping of an index (A(z) = B(z)); its value is the
       element's codomain position. Both used to be rejected. */
    char rr[NAMESIZE];
    k=0;
    for (p=rhs; *p!='\0'&&*p!=';'&&k<NAMESIZE-1; p++) if (*p!=' ') rr[k++]=*p;
    rr[k]='\0';
    for (idxl=0; idxl<fdim-1; idxl++) if (strcmp(rr,arSet[idxl].index_name)==0) break;
    if (idxl<fdim-1) {
      idxss=set_supset_slot(sets,(dim_t)arSet[idxl].setid,md->toset);
      if (idxss<0) {
        errmsg("Error: Formula for mapping %s: index %s ranges over set %s, which is neither the codomain %s nor a declared subset of it (manual 10.13.1)\n",md->mapname,rr,sets[arSet[idxl].setid].setname,sets[md->toset].setname);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      idxmode=1;
    } else if ((p=strchr(rr,MAPMARK))!=NULL&&strpbrk(rr,"+-*/^(){}")==NULL) {
      char bname[NAMESIZE];
      dim_t mb;
      memcpy(bname,rr,p-rr);
      bname[p-rr]='\0';
      for (mb=0; mb<teems_nmap; mb++) if (strcmp(bname,teems_maps[mb].mapname)==0) break;
      for (idxl=0; idxl<fdim-1; idxl++) if (strcmp(p+1,arSet[idxl].index_name)==0) break;
      if (mb==teems_nmap||idxl==fdim-1) {
        errmsg("Error: Formula for mapping %s: %s is not a mapping of one of the Formula's quantifier indices (manual 10.13.1)\n",md->mapname,rr);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      if (!teems_maps[mb].has_values) {
        errmsg("Error: mapping %s is used (in the Formula for mapping %s) before all of its values are assigned (manual 10.13.1/11.9.1)\n",teems_maps[mb].mapname,md->mapname);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      cpdss=set_supset_slot(sets,(dim_t)arSet[idxl].setid,teems_maps[mb].fromset);
      idxss=set_supset_slot(sets,teems_maps[mb].toset,md->toset);
      if (cpdss<0||idxss<0) {
        errmsg("Error: Formula for mapping %s: %s needs its index over the domain of %s (or a subset) and a codomain equal to or a subset of %s (manual 10.13.1)\n",md->mapname,rr,teems_maps[mb].mapname,sets[md->toset].setname);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      teems_maps[mb].used=true;
      cpmap=mb;
      idxmode=2;
    } else {
      if (byele) {
        errmsg("Error: Formula (by_elements) for mapping %s needs a quoted codomain element, an index of the codomain or another mapping on the right-hand side (manual 10.13.1)\n",md->mapname);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      if(!formula_compile(rhs,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdim-1))MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
  for (t=0; t<nloops; t++) {
    dim_t cod,d;
    i4=t;
    for (dcount=0; dcount<fdim-1; dcount++) {
      i3=(offset_t) i4/dcountdim1[dcount];
      arSet[dcount].indx=i3;
      i4=i4-i3*dcountdim1[dcount];
    }
    if (litdom) d=dom;
    else {
      d=arSet[l].indx;
      if (sup>0) d=set_elems[sets[arSet[l].setid].offset+d].superset_pos[sup];
    }
    if (codlit>=0) cod=codlit;
    else if (idxmode==1) {
      cod=(dim_t)arSet[idxl].indx;
      if (idxss>0) cod=set_elems[sets[arSet[idxl].setid].offset+cod].superset_pos[idxss];
    } else if (idxmode==2) {
      dim_t bd=(dim_t)arSet[idxl].indx;
      if (cpdss>0) bd=set_elems[sets[arSet[idxl].setid].offset+bd].superset_pos[cpdss];
      cod=teems_maps[cpmap].values[bd];
      if (idxss>0) cod=set_elems[sets[teems_maps[cpmap].toset].offset+cod].superset_pos[idxss];
    } else {
      solve_real v=formula_eval(elem_vals,sets,set_elems,sum_vals,ops,nops,arSet,fdim-1,zerodivide);
      cod=(dim_t)lround((double)v)-1;
    }
    mapping_store_value(mm,d,cod,sets,set_elems);
  }
}

/* The statement's top-level '=' , or a named abort. str_rfind_toplevel
   returns NULL when no '=' has balanced parentheses to its right: a
   statement with no '=' at all, or -- the case that reached users --
   an unclosed '(' in the right-hand side, as in
   `Formula (all,r,REG) Z(r) = sum(c,COMM, 1.0;`. All three callers
   below dereferenced the result unchecked, so such a statement faulted
   inside formulas_execute instead of being reported (fuzz batch 13; the
   sum_extract guard does not cover this path). teems-R's
   .chk_tab_parens rejects it before the solver runs, so this names it
   for the standalone solver. */
static char *formula_toplevel_eq(char *line, const char *commsyntax) {
  char *p=str_rfind_toplevel(line,'=');
  if (p==NULL) {
    errmsg("Error: %s statement has no top-level '=' (a missing '=' or unbalanced parentheses): %s\n",commsyntax,line);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  return p;
}

/* '(' inside quantifier conditions of a Formula LHS text (spaces
   stripped): they are not quantifiers or argument lists, so the fdim
   count by '(' must not include them */
static int quant_cond_paren_count(const char *s) {
  int n=0;
  const char *p=s;
  while (str_find_ci((char *)p,"(all,")==0) {
    int d=0,incond=0;
    for (; *p!='\0'; p++) {
      if (*p=='(') { if (incond) n++; d++; }
      else if (*p==')') { d--; if (d==0) { p++; break; } }
      else if (*p==':'&&d==1) incond=1;
    }
    if (d!=0) break;
  }
  return n;
}

static void cond_programs_free(char **condL, char **condR, formula_op *condops[][2]) {
  int i;
  for (i=0; i<MAXVARDIM; i++) {
    free(condL[i]); condL[i]=NULL;
    free(condR[i]); condR[i]=NULL;
    free(condops[i][0]); condops[i][0]=NULL;
    free(condops[i][1]); condops[i][1]=NULL;
  }
}

offset_t formulas_execute(char *fname, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,bool IsIni) {
  FILE * filehandle;
  char line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE],condvar[MAXVARDIM][NAMESIZE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE],tempset[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  offset_t i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],ncond,nloops,logioper[MAXVARDIM],logi,logiantidim[MAXVARDIM][MAXVARDIM],logisup[MAXVARDIM][MAXVARDIM],logivarindx[MAXVARDIM],logivartype[MAXVARDIM];//m,
  dim_t fdim,dcount,neqsign=0,varsupsetid[MAXVARDIM];
  int nops=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx,b=0;
  offset_t varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  solve_real zerodivide=0,cond[MAXVARDIM],eval;
  bool IsFomIni=false,IsDefFomIni=false;
  char fdefval[NAMESIZE];
  quantifier *arSet1=NULL;
  formula_op *ops1= NULL;
  /* general quantifier conditions (anything but "COEF(idx) op number"):
     both sides compiled as expressions over the statement frame */
  int condgen[MAXVARDIM];
  char *condL[MAXVARDIM],*condR[MAXVARDIM];
  formula_op *condops[MAXVARDIM][2];
  dim_t condnops[MAXVARDIM][2];
  int condcap[MAXVARDIM][2];
  for (i=0; i<MAXVARDIM; i++) { condL[i]=NULL; condR[i]=NULL; condops[i][0]=condops[i][1]=NULL; }

  filehandle = fopen(fname,"r");
  /* each pass rescans from the top, so the positional zerodivide state
     replays in file order (manual 10.11.1); this also gives the PostSim
     pass its fresh initial state (manual 12.2.4) */
  zdiv_scan_reset();
  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    /* mapping calls lower to flat map~idx tokens before any brace
       tokenizer runs (manual 11.9.4; design doc M2); sum carried-dim
       discovery identifies them by their domain index (M2c) */
    if (teems_nmap>0) mapping_lower_calls(line);
    /* audit A7 / plan 2.2: FORMULA & EQUATION expands to
       Formula (initial) + Equation (levels) (manual 10.9.1) -- the
       levels-equation half is outside this solver's linearized-only
       scope, so the combined form stays fatal by design rather than
       running the formula half while the equation silently vanishes */
    if (strstr(line,"& equation")!=NULL||strstr(line,"&equation")!=NULL) {
      errmsg("Error: 'Formula & Equation' is not supported: its expansion needs a levels equation (manual 10.9.1); linearize the equation and set the base value with Formula (initial)\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    /* positional INITIAL/ALWAYS default (manual 10.19; audit A6):
       both directions and spaced forms; values validated up front by
       tab_defaults_validate. Integer-LHS formulas nominally keep
       INITIAL regardless (10.19) -- not distinguished here. */
    if (tab_default_value(line,fdefval)) {
      if(strcmp(fdefval,"initial")==0)IsDefFomIni=true;
      else if(strcmp(fdefval,"always")==0)IsDefFomIni=false;
      continue;
    }
    {
      IsFomIni=IsDefFomIni;
      int byele=0;
      if(strstr(line, "(initial)")!=NULL) {
        str_replace_first(line, "(initial)", "");
        IsFomIni=true;
      }
      if(strstr(line, "(always)")!=NULL) {
        str_replace_first(line, "(always)", "");
        IsFomIni=false;
      }
      /* Formula (by_elements): a mapping assigned by codomain element
         name rather than position (manual 10.13.1/11.9.12) */
      if(strstr(line, "(by_elements)")!=NULL) {
        str_replace_first(line, "(by_elements)", "");
        byele=1;
      }
      if(!(IsFomIni&&!IsIni)) {
        zdiv_capture();
        ncond=0;
        for (i=0; i<MAXVARDIM; i++)logioper[i]=0;
        str_replace_first(line, commsyntax, "");
        /* if() conditions may spell the comparison as a word
           (GTAP-AEZ: IF[X LE 0, ..]); rewritten before the spaces go */
        tab_wordops_normalize(line);
        str_delete_char(line,' ');
        while (str_replace_char(line, '[', '('));
        while (str_replace_char(line, ']', ')'));
        while (str_replace_char(line, '{', '('));
        while (str_replace_char(line, '}', ')'));
        strcpy(linecopy,line);
        totalsum=str_count_ci(line, "sum(");
        neqsign=str_count_char(line, '=');
        readitem = strtok(line,"=");
        for(i=1; i<neqsign; i++)readitem = strtok(NULL,"=");
        readitem = strtok(NULL,";");
        npow=str_count_char(readitem, '^');
        nmul=str_count_char(readitem, '*');
        ndiv=str_count_char(readitem, '/');
        nmul=nmul+ndiv;
        nplu=str_count_char(readitem, '+');
        nmin=str_count_char(readitem, '-');
        nplu=nplu+nmin;

        strcpy(line,linecopy);
        readitem =formula_toplevel_eq(line,commsyntax);
        line[readitem-line]='\0';
        fdim=str_count_char(line, '(')-quant_cond_paren_count(line);
        /* unquantified mapping assignment MAP("dom") = "cod" | <pos>
           (manual 10.13.1): no tuple machinery involved */
        if (teems_nmap>0&&strstr(line,"(all,")==NULL&&(strchr(line,'(')!=NULL||strchr(line,MAPMARK)!=NULL)) {
          char mname[NAMESIZE];
          dim_t mm;
          int mk=0;
          for (p=line; *p!='\0'&&*p!='('&&*p!=MAPMARK&&mk<NAMESIZE-1; p++) mname[mk++]=*p;
          mname[mk]='\0';
          for (mm=0; mm<teems_nmap; mm++) if (strcmp(mname,teems_maps[mm].mapname)==0) break;
          if (mm<teems_nmap) {
            char *rhs;
            /* values are fixed once the mapping has been used (manual
               11.9.9): later formula passes re-run the statement
               sequence, whose final state is already in the table */
            if (!teems_maps[mm].used) {
              strcpy(line,linecopy);
              rhs=formula_toplevel_eq(line,commsyntax);
              *rhs='\0';
              rhs++;
              mapping_assign_literal(mm,line,rhs,sets,set_elems,byele);
            }
            j++;
            continue;
          }
        }
        if (fdim==1) {
          fdim=fdim+1;
        }
        quantifier *arSet= (quantifier *) calloc (fdim+1,sizeof(quantifier));

        nloops=1;
        if (fdim==0) {
          readitem = line;
        } else {
          /* quantifier groups parsed depth-aware: a condition may carry
             several parenthesised groups ((all,r,REG:$POS(c)=$POS(g)),
             :A(d)>A(d)), which the strtok walk used to count as
             quantifiers until it read past the statement (tpmh0147
             SEGV, potential-models segfaults.md #3) */
          char *qp=line;
          for (i=0; i<fdim-1; i++) {
            char grp[TABREADLINE],*qe,*c1,*c2,*cc;
            int qd=0;
            condgen[i]=0;
            if (str_find_ci(qp,"(all,")!=0) {
              errmsg("Error: malformed quantifier list in Formula: %s\n",linecopy);
              MPI_Abort(PETSC_COMM_WORLD,1);
            }
            for (qe=qp; *qe!='\0'; qe++) {
              if (*qe=='(') qd++;
              else if (*qe==')') { qd--; if (qd==0) break; }
            }
            if (*qe!=')'||qe-qp-1>=(long)sizeof(grp)) {
              errmsg("Error: unbalanced quantifier in Formula: %s\n",linecopy);
              MPI_Abort(PETSC_COMM_WORLD,1);
            }
            memcpy(grp,qp+1,qe-qp-1);
            grp[qe-qp-1]='\0';
            c1=strchr(grp,',');
            c2=(c1!=NULL)?strchr(c1+1,','):NULL;
            if (c1==NULL||c2==NULL||c2-c1-1<=0||c2-c1-1>=NAMESIZE) {
              errmsg("Error: malformed quantifier (%s) in Formula: %s\n",grp,linecopy);
              MPI_Abort(PETSC_COMM_WORLD,1);
            }
            memcpy(arSet[i].index_name,c1+1,c2-c1-1);
            arSet[i].index_name[c2-c1-1]='\0';
            cc=strchr(c2+1,':');
            if (cc==NULL) {
              if (strlen(c2+1)>=NAMESIZE) { errmsg("Error: malformed quantifier (%s) in Formula: %s\n",grp,linecopy); MPI_Abort(PETSC_COMM_WORLD,1); }
              strcpy(tempset,c2+1);
            } else {
              char *cond_s=cc+1,*op=NULL;
              int oplen=0,od=0;
              if (cc-c2-1>=NAMESIZE) { errmsg("Error: malformed quantifier (%s) in Formula: %s\n",grp,linecopy); MPI_Abort(PETSC_COMM_WORLD,1); }
              memcpy(tempset,c2+1,cc-c2-1);
              tempset[cc-c2-1]='\0';
              /* the top-level comparison of the condition (manual 11.4.11) */
              for (p=cond_s; *p!='\0'; p++) {
                if (*p=='('||*p=='{'||*p=='[') od++;
                else if (*p==')'||*p=='}'||*p==']') od--;
                else if (od==0&&(*p=='<'||*p=='>'||*p=='=')) {
                  op=p;
                  oplen=((p[0]=='<'&&(p[1]=='='||p[1]=='>'))||(p[0]=='>'&&p[1]=='='))?2:1;
                  break;
                }
              }
              if (op==NULL) {
                errmsg("Error: Formula quantifier condition %s has no comparison (=, <>, <, >, <=, >=; manual 11.4.11): %s\n",cond_s,linecopy);
                MPI_Abort(PETSC_COMM_WORLD,1);
              }
              if (op[0]=='=') logioper[i]=1;
              else if (op[0]=='<'&&op[1]=='=') logioper[i]=5;
              else if (op[0]=='>'&&op[1]=='=') logioper[i]=6;
              else if (op[0]=='<'&&op[1]=='>') logioper[i]=4;
              else if (op[0]=='>') logioper[i]=2;
              else logioper[i]=3;
              {
                char lhs_s[TABREADLINE],*rhs_s=op+oplen,*endn=NULL;
                double num;
                int simple=0;
                memcpy(lhs_s,cond_s,op-cond_s);
                lhs_s[op-cond_s]='\0';
                /* a mapping compared with an element or another mapping
                   (MAP(i)=r) has no numeric value to compare: named
                   fatal (M3). A mapped ARGUMENT inside a numeric
                   comparison evaluates like any expression. */
                {
                  /* a bare mapping side (MAP(z) = r, MAP(z) = "ele",
                     A(z) = B(y)) compares codomain positions: it becomes
                     $POS(MAP(z)) and the other side the position of the
                     index or element in that codomain (vetting G-S7) */
                  char rhs_buf[TABREADLINE];
                  int lm=(strchr(lhs_s,MAPMARK)!=NULL&&strpbrk(lhs_s,"(){}+-*/^$")==NULL);
                  int rm=(strchr(rhs_s,MAPMARK)!=NULL&&strpbrk(rhs_s,"(){}+-*/^$")==NULL);
                  if (lm||rm) {
                    char *ms=lm?lhs_s:rhs_s,*os=lm?rhs_s:lhs_s,mname[NAMESIZE],other[TABREADLINE],*mk=strchr(ms,MAPMARK);
                    dim_t mm2,cset,q2;
                    if (logioper[i]!=1&&logioper[i]!=4) {
                      errmsg("Error: a mapping in a Formula quantifier condition compares by = or <> only (manual 11.4.11): %s\n",linecopy);
                      MPI_Abort(PETSC_COMM_WORLD,1);
                    }
                    memcpy(mname,ms,mk-ms); mname[mk-ms]='\0';
                    for (mm2=0; mm2<teems_nmap; mm2++) if (strcmp(mname,teems_maps[mm2].mapname)==0) break;
                    if (mm2==teems_nmap) { errmsg("Error: unknown mapping %s in a Formula quantifier condition: %s\n",mname,linecopy); MPI_Abort(PETSC_COMM_WORLD,1); }
                    cset=teems_maps[mm2].toset;
                    if (strchr(os,MAPMARK)!=NULL&&strpbrk(os,"(){}+-*/^$")==NULL) { if (snprintf(other,sizeof(other),"$pos(%s)",os)>=(int)sizeof(other)) { errmsg("Error: Formula quantifier condition too long: %s\n",linecopy); MPI_Abort(PETSC_COMM_WORLD,1); } }
                    else {
                      for (q2=0; q2<=i; q2++) if (strcmp(os,arSet[q2].index_name)==0) break;
                      if (q2<=i) { if (snprintf(other,sizeof(other),"$pos(%s,%s)",os,sets[cset].setname)>=(int)sizeof(other)) { errmsg("Error: Formula quantifier condition too long: %s\n",linecopy); MPI_Abort(PETSC_COMM_WORLD,1); } }
                      else {
                        char el[NAMESIZE];
                        int ek=0;
                        const char *oq=os;
                        dim_t e;
                        for (; *oq!='\0'&&ek<NAMESIZE-1; oq++) if (*oq!='\"') el[ek++]=(char)tolower((int)*oq);
                        el[ek]='\0';
                        for (e=0; e<sets[cset].size; e++) if (strcmp(el,set_elems[sets[cset].offset+e].setele)==0) break;
                        if (e==sets[cset].size) {
                          errmsg("Error: %s in a Formula quantifier condition is neither a quantifier index nor an element of %s, the codomain of mapping %s: %s\n",os,sets[cset].setname,mname,linecopy);
                          MPI_Abort(PETSC_COMM_WORLD,1);
                        }
                        snprintf(other,sizeof(other),"%d",(int)e+1);
                      }
                    }
                    if (snprintf(rhs_buf,sizeof(rhs_buf),"$pos(%s)",ms)>=(int)sizeof(rhs_buf)) { errmsg("Error: Formula quantifier condition too long: %s\n",linecopy); MPI_Abort(PETSC_COMM_WORLD,1); }
                    strcpy(lhs_s,rhs_buf);
                    strcpy(rhs_buf,other);
                    rhs_s=rhs_buf;
                  }
                num=strtod(rhs_s,&endn);
                /* the historical fast form: COEF or COEF(indices) against
                   a number */
                if (endn!=rhs_s&&*endn=='\0'&&strlen(lhs_s)<NAMESIZE&&strchr(lhs_s,MAPMARK)==NULL) {
                  char *q=lhs_s;
                  if (isalpha((unsigned char)*q)) {
                    while (isalnum((unsigned char)*q)||*q=='_'||*q=='@') q++;
                    if (*q=='\0') simple=1;
                    else if (*q=='('&&strchr(q+1,'(')==NULL&&strchr(q,')')==q+strlen(q)-1) simple=1;
                  }
                }
                if (simple) {
                  strcpy(condvar[i],lhs_s);
                  cond[i]=num;
                } else {
                  /* anything else -- $POS, arithmetic, a coefficient on
                     the right (it was atof'd to 0) -- is evaluated as two
                     expressions over the frame */
                  if (str_count_ci(cond_s,"sum(")>0) {
                    errmsg("Error: sums in Formula quantifier conditions are not supported; compute the sum in a coefficient first: %s\n",linecopy);
                    MPI_Abort(PETSC_COMM_WORLD,1);
                  }
                  condgen[i]=1;
                  condL[i]=strdup(lhs_s);
                  condR[i]=strdup(rhs_s);
                  condvar[i][0]='\0';
                  cond[i]=0;
                }
                }
              }
              ncond++;
            }
            for (i4=0; i4<nset; i4++) if(strcmp(tempset,sets[i4].setname)==0) {
                arSet[i].setid=i4;
                break;
              }
            nloops=nloops*sets[arSet[i].setid].size;
            qp=qe+1;
          }
          readitem = qp;
          if (*readitem=='\0') {
            errmsg("Error: Formula has no left-hand side after its quantifiers: %s\n",linecopy);
            MPI_Abort(PETSC_COMM_WORLD,1);
          }
          dcountdim1[fdim-2]=1;
          for (i=fdim-3; i>-1; i--) {
            dcountdim1[i]=sets[arSet[i+1].setid].size*dcountdim1[i+1];
          }
        }
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem =formula_toplevel_eq(line,commsyntax);
        readitem++;
        readitem = strtok(readitem,";");
        while (formula_normalize(readitem)==1);
        leadlag_encode(readitem);
        npar=str_count_char(readitem,'(')+str_count_char(readitem,',')+str_count_ci(readitem,"$pos");/* comma slack: multi-arg intrinsics emit per-arg loads + folds (plan 3.1); $POS calls: one position op each, their index lists ride in {} after normalize (fuzz batch 13) */
        strcpy(sumsyntax,"sum(");
        totalsum=sum_count(readitem,sumsyntax);
        sum_def *sum_cof= (sum_def *) calloc (totalsum,sizeof(sum_def));
        sumcount=0;
        strcpy(line1,readitem);
        strcpy(line2,line1);
        readitem=line2;
        while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        i3=0;
        for (i=0; i<totalsum; i++) {
          i1=1;
          for(j=0; j<sum_cof[i].size; j++) {
            i1=i1*sets[sum_cof[i].setid[j]].size;
          }
          sum_cof[i].offset=i3;
          sum_cof[i].summatsize=i1;
          i3=i3+i1;
        }
        nsumele=i3;
        for (i=0; i<totalsum; i++) {
          i1=1;
          /* a scalar sum store (size 0) has no strides; the unguarded
             write was strides[-1] = setid[MAXVARDIM-1], benign only
             because that slot is never read at size 0 */
          if (sum_cof[i].size==0) continue;
          sum_cof[i].strides[sum_cof[i].size-1]=1;
          for(j=sum_cof[i].size-2; j>-1; j--) {
            sum_cof[i].strides[j]=sum_cof[i].strides[j+1]*sets[sum_cof[i].setid[j+1]].size;
          }
        }
        formula_op *ops= (formula_op *) calloc (npow+nmul+nplu+2*(npar+2),sizeof(formula_op));
        sum_value *sum_vals= (sum_value *) calloc (nsumele,sizeof(sum_value));
        sumcount=0;
        strcpy(line2,line1);
        readitem=line2;
        sumindx=0;
        while (sum_eval(readitem,sumsyntax,sets,nset,set_elems,elem_vals,ncofvar,ncofele,coefs,ncof,vars,nvar,sum_cof,totalsum,sum_vals,nsumele,ops,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
          sumcount++;
        }
        strcpy(line1,readitem);
        /* quantified mapping assignment (manual 10.13.1): the LHS is a
           declared mapping, the RHS a codomain position ($POS shapes) or,
           under (by_elements), a quoted codomain element */
        if (teems_nmap>0) {
          char mname[NAMESIZE];
          dim_t mm;
          int mk=0;
          for (p=vname; *p!='\0'&&*p!='('&&*p!=MAPMARK&&mk<NAMESIZE-1; p++) mname[mk++]=*p;
          mname[mk]='\0';
          for (mm=0; mm<teems_nmap; mm++) if (strcmp(mname,teems_maps[mm].mapname)==0) break;
          if (mm<teems_nmap) {
            if (!teems_maps[mm].used) mapping_assign_formula(mm,vname,line1,byele,sets,set_elems,sum_vals,elem_vals,arSet,fdim,dcountdim1,nloops,ops,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,zerodivide);
            free(sum_cof);
            free(sum_vals);
            free(arSet);
            free(ops);
            cond_programs_free(condL,condR,condops);
            j++;
            continue;
          }
        }
        offset_t index=ncof-1, offset=0;//,simpl=0;
        bool check10=true;
        offset_t varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            /* F2: parameters may only be set by Read or Formula
               (Initial) -- warn once (first pass) rather than abort,
               since TEEMS never enforced the distinction */
            if(!teems_ps_pass&&IsIni&&!IsFomIni&&teems_coef_is_param!=NULL&&teems_coef_is_param[index])printf("Warning: Formula (always) assigns (parameter) coefficient %s; GEMPACK only allows Read or Formula (Initial) for parameters\n",coefs[index].cofname);
            coefs[index].suplval=true;
            offset=coefs[index].offset;
            varsize=coefs[index].size;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(vars[index].cofname,p)==0) {
              offset=ncofele+vars[index].offset;
              varsize=vars[index].size;
              vars[index].suplval=true;
              if(vars[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        /* a satisfied Read (IfHeaderExists) supersedes formulas
           assigning the read coefficient (manual 11.11.8 idiom under
           the reads-then-formulas pass order, plan 3.9): keep all the
           statement plumbing, evaluate nothing */
        if (!check10&&teems_coef_ifhdr!=NULL&&teems_coef_ifhdr[index]) {
          nloops=0;
        }
        /* PostSim Formula LHS must be a PostSim Coefficient (manual
           12.2.2): never a Variable, never an ordinary Coefficient */
        if (teems_ps_pass) {
          if (check10) {
            errmsg("Error: PostSim Formula assigns variable %s; simulation results cannot be changed (manual 12.2.2)\n",vars[index].cofname);
            MPI_Abort(PETSC_COMM_WORLD,1);
          }
          if (teems_coef_is_ps==NULL||!teems_coef_is_ps[index]) {
            errmsg("Error: PostSim Formula assigns ordinary coefficient %s; the LHS must be a PostSim Coefficient (manual 12.2.2)\n",coefs[index].cofname);
            MPI_Abort(PETSC_COMM_WORLD,1);
          }
        }
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<vars[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=vars[index].strides[dcount];
                  vararset[dcount]=l+1;
                  { dim_t ss=set_supset_slot(sets,arSet[l].setid,vars[index].setid[dcount]); if(ss<0)set_supset_fatal(p,vars[index].cofname,NULL,sets,arSet[l].setid,vars[index].setid[dcount]); if(ss>0){varsubset[dcount]=1; varsupsetid[dcount]=ss;} }
                  break;
                }
            }
        } else {
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=coefs[index].strides[dcount];
                  vararset[dcount]=l+1;
                  { dim_t ss=set_supset_slot(sets,arSet[l].setid,coefs[index].setid[dcount]); if(ss<0)set_supset_fatal(p,coefs[index].cofname,NULL,sets,arSet[l].setid,coefs[index].setid[dcount]); if(ss>0){varsubset[dcount]=1; varsupsetid[dcount]=ss;} }
                  break;
                }
            }
        }
        if(!formula_compile(line1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdim-1))MPI_Abort(PETSC_COMM_WORLD,1);
        if(ncond>0) {
          for(i=0; i<MAXVARDIM; i++)for(j=0; j<MAXVARDIM; j++){
            logiantidim[i][j]=0;
            logisup[i][j]=0;
          }
          for(i1=0; i1<fdim-1; i1++) if(logioper[i1]>0&&condgen[i1]) {
            int sd;
            for (sd=0; sd<2; sd++) {
              char ctext[TABREADLINE];
              int cn;
              strcpy(ctext,sd==0?condL[i1]:condR[i1]);
              while (formula_normalize(ctext)==1);
              leadlag_encode(ctext);
              cn=str_count_char(ctext,'^')+str_count_char(ctext,'*')+str_count_char(ctext,'/')+str_count_char(ctext,'+')+str_count_char(ctext,'-')
                 +2*(str_count_char(ctext,'(')+str_count_char(ctext,',')+str_count_char(ctext,'{')+str_count_ci(ctext,"$pos")+2)+4;
              condops[i1][sd]=(formula_op *)calloc(cn,sizeof(formula_op));
              condcap[i1][sd]=cn;
              condnops[i1][sd]=0;
              if(!formula_compile(ctext,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,condops[i1][sd],&condnops[i1][sd],arSet,fdim-1)) {
                errmsg("Error: cannot evaluate the Formula quantifier condition %s: %s\n",sd==0?condL[i1]:condR[i1],linecopy);
                MPI_Abort(PETSC_COMM_WORLD,1);
              }
            }
          }
          for(i1=0; i1<fdim-1; i1++) if(logioper[i1]>0&&!condgen[i1]){
            index=ncof-1;
            p=strtok(condvar[i1],"(");
            b=0;
            do {
              if (strcmp(coefs[index].cofname,condvar[i1])==0) {
                if(!coefs[index].suplval)warn_no_values(coefs[index].cofname,index,0);
                logivarindx[i1]=index;
                logivartype[i1]=0;
                b++;
                /* a scalar condition coefficient has no argument list:
                   strtok returns NULL and the copy used to SEGV */
                if(p!=NULL)p=strtok(NULL,")");
                if(p!=NULL)strcpy(argu,p);
                else argu[0]='\0';
                strcat(argu,",");
                if(cond_arg_count(argu)!=(int)coefs[index].size) {
                  errmsg("Error: quantifier condition %s carries %d argument(s); %s is declared with %d (manual 11.4.11)\n",coefs[index].cofname,cond_arg_count(argu),coefs[index].cofname,(int)coefs[index].size);
                  MPI_Abort(PETSC_COMM_WORLD,1);
                }
                for(i=0; i<coefs[index].size; i++) {
                  if(i==0)p=strtok(argu,",");
                  else p=strtok(NULL,",");
                  for(j=0; j<fdim; j++) {
                    if(strcmp(arSet[j].index_name,p)==0) {
                      logiantidim[i1][j]=coefs[index].strides[i];
                      { dim_t ss=set_supset_slot(sets,arSet[j].setid,coefs[index].setid[i]); if(ss<0)set_supset_fatal(p,coefs[index].cofname,NULL,sets,arSet[j].setid,coefs[index].setid[i]); logisup[i1][j]=(ss>0)?ss:0; }
                      break;
                    }
                  }
                }
                break;
              }
            } while (index--);
            if(b==0) {
              index=nvar-1;
              p=strtok(condvar[i1],"(");
              do {
                if (strcmp(vars[index].cofname,condvar[i1])==0) {
                if(!vars[index].suplval)warn_no_values(vars[index].cofname,index,1);
                  logivarindx[i1]=index;
                  logivartype[i1]=1;
                  b++;
                  if(p!=NULL)p=strtok(NULL,")");
                  if(p!=NULL)strcpy(argu,p);
                  else argu[0]='\0';
                  strcat(argu,",");
                  if(cond_arg_count(argu)!=(int)vars[index].size) {
                    errmsg("Error: quantifier condition %s carries %d argument(s); %s is declared with %d (manual 11.4.11)\n",vars[index].cofname,cond_arg_count(argu),vars[index].cofname,(int)vars[index].size);
                    MPI_Abort(PETSC_COMM_WORLD,1);
                  }
                  for(i=0; i<vars[index].size; i++) {
                    if(i==0)p=strtok(argu,",");
                    else p=strtok(NULL,",");
                    for(j=0; j<fdim-1; j++) {
                      if(strcmp(arSet[j].index_name,p)==0) {
                        logiantidim[i1][j]=vars[index].strides[i];
                        { dim_t ss=set_supset_slot(sets,arSet[j].setid,vars[index].setid[i]); if(ss<0)set_supset_fatal(p,vars[index].cofname,NULL,sets,arSet[j].setid,vars[index].setid[i]); logisup[i1][j]=(ss>0)?ss:0; }
                        break;
                      }
                    }
                  }
                  break;
                }
              } while (index--);
            }
            if(b==0) {
              errmsg("Error: quantifier condition refers to %s, which is not a declared coefficient or variable (manual 11.4.11)\n",condvar[i1]);
              MPI_Abort(PETSC_COMM_WORLD,1);
            }
          }
        /* i: inner condition loop's counter — was shared-by-omission,
           a result-corrupting write-write race at maxthreads>1 */
        #pragma omp parallel private(l,l2,i,i4,dcount,i3,i1,arSet1,logi,index,eval,ops1) shared(elem_vals,arSet)
        {
        formula_op *cops[MAXVARDIM][2];
        {
          int ci,sd;
          for (ci=0; ci<MAXVARDIM; ci++) for (sd=0; sd<2; sd++) {
            cops[ci][sd]=NULL;
            if (ci<fdim-1&&condgen[ci]&&logioper[ci]>0) {
              if (omp_get_thread_num()!=0) {
                cops[ci][sd]=malloc(condcap[ci][sd]*sizeof(formula_op));
                memcpy(cops[ci][sd],condops[ci][sd],condcap[ci][sd]*sizeof(formula_op));
              } else cops[ci][sd]=condops[ci][sd];
            }
          }
        }
        if(omp_get_thread_num()!=0){
          arSet1=malloc((fdim+1)*sizeof(quantifier));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(quantifier));
          ops1=malloc((npow+nmul+nplu+2*(npar+2))*sizeof(formula_op));
          memcpy (ops1,ops,(npow+nmul+nplu+2*(npar+2))*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet1=arSet;
        }
    /* workers must finish copying the shared scratch (master aliases it)
       before the master starts mutating per-iteration state */
    #pragma omp barrier
        #pragma omp for
          for (l=0; l<nloops; l++) {
            l2=0;
            i4=l;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(offset_t) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<fdim-1||varsize==fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
            logi=0;
            index=0;
            for(i1=0; i1<fdim-1; i1++) {
              if(logioper[i1]>0){
              solve_real cv=cond[i1];
              if(condgen[i1]) {
                eval=formula_eval(elem_vals,sets,set_elems,sum_vals,cops[i1][0],condnops[i1][0],arSet1,fdim-1,zerodivide);
                cv=formula_eval(elem_vals,sets,set_elems,sum_vals,cops[i1][1],condnops[i1][1],arSet1,fdim-1,zerodivide);
              } else {
              index=0;
              for(i=0; i<fdim-1; i++){
                index+=set_elems[sets[arSet1[i].setid].offset+arSet1[i].indx].superset_pos[logisup[i1][i]]*logiantidim[i1][i];
              }
              if(logivartype[i1]==0)eval=elem_vals[coefs[logivarindx[i1]].offset+index].value;
              else eval=elem_vals[ncofele+vars[logivarindx[i1]].offset+index].value;
              }
              if(logioper[i1]==1)if(eval==cv)logi++;
              if(logioper[i1]==2)if(eval>cv)logi++;
              if(logioper[i1]==3)if(eval<cv)logi++;
              if(logioper[i1]==4)if(eval!=cv)logi++;
              if(logioper[i1]==5)if(eval<=cv)logi++;
              if(logioper[i1]==6)if(eval>=cv)logi++;
              }
            }
            if(logi==ncond)elem_vals[offset+l2].value=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
          }
        if(omp_get_thread_num()!=0){
          int ci;
          for (ci=0; ci<MAXVARDIM; ci++) { free(cops[ci][0]); free(cops[ci][1]); }
          free(arSet1);
          arSet1=NULL;
          free(ops1);
          ops1=NULL;
        }else{
          ops1=NULL;
          arSet1=NULL;
        }
        }
        } else {
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,arSet1,ops1) shared(elem_vals,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=malloc((fdim+1)*sizeof(quantifier));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(quantifier));
          ops1=malloc((npow+nmul+nplu+2*(npar+2))*sizeof(formula_op));
          memcpy (ops1,ops,(npow+nmul+nplu+2*(npar+2))*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet1=arSet;
        }
    /* workers must finish copying the shared scratch (master aliases it)
       before the master starts mutating per-iteration state */
    #pragma omp barrier
        #pragma omp for
          for (l=0; l<nloops; l++) {
            l2=0;
            i4=l;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(offset_t) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
            elem_vals[offset+l2].value=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
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

        int glmode=IsIni?teems_range_test_initial:teems_range_test_updated;
        if(teems_rk_stage_checks&&glmode==2)glmode=1;   /* RK stage state: warn and let the driver retry */
        if(glmode>0){
          int glviol=0;
          if(coefs[index].gltype>0)glviol|=coef_range_check(coefs,index,offset,varsize,elem_vals,glmode,coefs[index].gltype,coefs[index].glval,!IsIni);
          if(teems_coef_gltype2!=NULL&&teems_coef_gltype2[index]>0)glviol|=coef_range_check(coefs,index,offset,varsize,elem_vals,glmode,teems_coef_gltype2[index],teems_coef_glval2[index],!IsIni);
          /* fatal only when requested (-range_test_* 2; manual
             25.4.4); the GEMPACK default is warn */
          if(glviol) {
            if(teems_rk_stage_checks)teems_check_viol_range++;
            else if(glmode==2)MPI_Abort(PETSC_COMM_WORLD,1);
          }
        }
        
        
        
        free(sum_cof);
        free(sum_vals);
        free(arSet);
        free(ops);
        cond_programs_free(condL,condR,condops);
      }
    }
  }
  fclose(filehandle);
  /* zerodivide never applies outside formulas (manual 10.11.1) --
     downstream update/equation evaluation must see it off */
  zdiv_disable();
  return j;
}

/* Applies the TAB update statements to coefficient values. With
   midpoint!=0 the modified-midpoint correction is used: the value is
   advanced from the sub-step base (csolpupd) by twice the computed
   change, and csolpupd retains the pre-update value. */
/* ---- one Update pass (manual 11.12, 4.4.4): the right-hand side of
   every Update statement uses the values at the START of the step, and
   several Updates of one coefficient apply in order, the later one
   overriding. Statements used to write in place, so a later statement
   read an earlier one's new value (Update V2 = x; Update (change) V =
   0.01*V2*x gave 30 instead of 20; tpmh0114 LAGY 8.24 instead of 4.04 --
   potential-models vetting SW5). A statement whose results another
   statement of the pass could read, or which shares its target with
   another statement, holds its results until the pass ends; the rest
   (product updates of a coefficient nothing else touches -- they read
   only their own element) still write in place. ---- */
static int upd_nstmt=-1;
static unsigned char *upd_defer=NULL;
static offset_t upd_pend_n=0,upd_pend_cap=0;
static offset_t *upd_pend_off=NULL,*upd_pend_sidx=NULL;
static store_real *upd_pend_val=NULL,*upd_pend_base=NULL;
static double *upd_pend_sv=NULL,*upd_pend_sb=NULL;
typedef struct { offset_t index,offset,varsize; } upd_check;
static upd_check *upd_checks=NULL;
static int upd_nchecks=0,upd_capchecks=0;

/* (change)/(explicit) targets whose final value is extrapolated from the
   per-pass path values like the variables (manual 26.2: updated data is
   extrapolated from the separate multi-step calculations) */
static offset_t *upd_pathbase=NULL;     /* per coefficient: base in upd_pathacc, -1 = not a target */
static double *upd_pathacc=NULL;
static offset_t upd_pathn=0;
int teems_upd_pathuse=0;
/* double-precision shadow of those targets' path (value, previous-step
   value): single-precision storage rounds each step's increment onto
   the value, and a change update adding ~1e-6 per step onto ~1 loses a
   tenth of it; the shadow accumulates the increments exactly and only
   feeds the extrapolation -- the stored values the steps use are
   unchanged */
static double *upd_shv=NULL,*upd_shb=NULL;
static int upd_shvalid=0;

static offset_t upd_coef_lookup(const char *name, array_def *coefs, offset_t ncof) {
  offset_t k;
  for (k=0; k<ncof; k++) if (strcmp(coefs[k].cofname,name)==0) return k;
  return -1;
}

static void upd_analyse(char *fname, elem_value *elem_vals, array_def *coefs, offset_t ncof) {
  FILE *f;
  char line[TABREADLINE],cs[NAMESIZE],name[NAMESIZE];
  solve_real zd=0;
  int n=0,cap=16,k;
  offset_t *tgt=NULL;
  unsigned char *ce=NULL,*cnt=NULL,*readset=NULL;
  if (upd_nstmt>=0) return;
  tgt=malloc(cap*sizeof(offset_t));
  ce=malloc(cap);
  cnt=calloc(ncof+1,1);
  readset=calloc(ncof+1,1);
  upd_pathbase=malloc((ncof+1)*sizeof(offset_t));
  for (k=0; k<ncof; k++) upd_pathbase[k]=-1;
  strcpy(cs,"update");
  f=fopen(fname,"r");
  while (f!=NULL&&tab_next_statement_resolved(cs,f,line,elem_vals,coefs,ncof,&zd,TABREADLINE)) {
    char *p=line,*eq;
    int isce=0,j;
    if (n==cap) { cap*=2; tgt=realloc(tgt,cap*sizeof(offset_t)); ce=realloc(ce,cap); }
    isce=(str_find_ci(line,"(change)")>=0||str_find_ci(line,"(explicit)")>=0);
    while (*p==' ') p++;
    if (str_find_ci(p,"update")==0) p+=6;
    for (;;) {
      int d=0;
      while (*p==' ') p++;
      if (*p!='(') break;
      for (; *p!='\0'; p++) { if (*p=='(') d++; else if (*p==')') { d--; if (d==0) { p++; break; } } }
    }
    j=0;
    while ((isalnum((unsigned char)*p)||*p=='_'||*p=='@')&&j<NAMESIZE-1) name[j++]=*p++;
    name[j]='\0';
    tgt[n]=upd_coef_lookup(name,coefs,ncof);
    ce[n]=(unsigned char)isce;
    if (tgt[n]>=0&&cnt[tgt[n]]<2) cnt[tgt[n]]++;
    if (isce&&tgt[n]>=0&&upd_pathbase[tgt[n]]<0) { upd_pathbase[tgt[n]]=upd_pathn; upd_pathn+=coefs[tgt[n]].nelem; }
    eq=strchr(line,'=');
    if (isce&&eq!=NULL) {
      char *q=eq+1;
      while (*q!='\0') {
        if (isalpha((unsigned char)*q)) {
          j=0;
          while ((isalnum((unsigned char)*q)||*q=='_'||*q=='@')&&j<NAMESIZE-1) name[j++]=*q++;
          name[j]='\0';
          while (isalnum((unsigned char)*q)||*q=='_'||*q=='@') q++;
          { offset_t r=upd_coef_lookup(name,coefs,ncof); if (r>=0) readset[r]=1; }
        } else q++;
      }
    }
    n++;
  }
  if (f!=NULL) fclose(f);
  upd_defer=malloc(n>0?n:1);
  for (k=0; k<n; k++) upd_defer[k]=(unsigned char)(ce[k]||tgt[k]<0||cnt[tgt[k]]>1||readset[tgt[k]]);
  upd_nstmt=n;
  free(tgt); free(ce); free(cnt); free(readset);
}

static int upd_stmt_defer(int s) {
  return (s<upd_nstmt)?upd_defer[s]:1;
}

static offset_t upd_pend_reserve(offset_t n) {
  offset_t b=upd_pend_n;
  if (upd_pend_n+n>upd_pend_cap) {
    offset_t c=upd_pend_cap?upd_pend_cap:1024;
    while (c<upd_pend_n+n) c*=2;
    upd_pend_off=realloc(upd_pend_off,c*sizeof(offset_t));
    upd_pend_sidx=realloc(upd_pend_sidx,c*sizeof(offset_t));
    upd_pend_val=realloc(upd_pend_val,c*sizeof(store_real));
    upd_pend_base=realloc(upd_pend_base,c*sizeof(store_real));
    upd_pend_sv=realloc(upd_pend_sv,c*sizeof(double));
    upd_pend_sb=realloc(upd_pend_sb,c*sizeof(double));
    if (upd_pend_off==NULL||upd_pend_sidx==NULL||upd_pend_val==NULL||upd_pend_base==NULL||upd_pend_sv==NULL||upd_pend_sb==NULL) {
      errmsg("Error: out of memory holding Update results\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    upd_pend_cap=c;
  }
  upd_pend_n+=n;
  return b;
}

static void upd_check_later(offset_t index, offset_t offset, offset_t varsize) {
  if (upd_nchecks==upd_capchecks) {
    upd_capchecks=upd_capchecks?2*upd_capchecks:16;
    upd_checks=realloc(upd_checks,upd_capchecks*sizeof(upd_check));
  }
  upd_checks[upd_nchecks].index=index;
  upd_checks[upd_nchecks].offset=offset;
  upd_checks[upd_nchecks].varsize=varsize;
  upd_nchecks++;
}

static void upd_range_check(array_def *coefs, offset_t index, offset_t offset, offset_t varsize, elem_value *elem_vals) {
  int glmode=teems_range_test_updated;
  if(teems_rk_stage_checks&&glmode==2)glmode=1;   /* RK stage state: warn and let the driver retry */
  if(glmode>0){
    int glviol=0;
    if(coefs[index].gltype>0)glviol|=coef_range_check(coefs,index,offset,varsize,elem_vals,glmode,coefs[index].gltype,coefs[index].glval,1);
    if(teems_coef_gltype2!=NULL&&teems_coef_gltype2[index]>0)glviol|=coef_range_check(coefs,index,offset,varsize,elem_vals,glmode,teems_coef_gltype2[index],teems_coef_glval2[index],1);
    /* fatal only when requested (-range_test_* 2; manual 25.4.4); the
       GEMPACK default is warn */
    if(glviol) {
      if(teems_rk_stage_checks)teems_check_viol_range++;
      else if(glmode==2)MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
}

/* end of a pass: held results land in statement order (a later Update
   of the same element overrides), then their range checks run */
static void upd_pass_flush(array_def *coefs, elem_value *elem_vals) {
  offset_t k;
  int c;
  for (k=0; k<upd_pend_n; k++) {
    elem_vals[upd_pend_off[k]].value=upd_pend_val[k];
    elem_vals[upd_pend_off[k]].substep_base=upd_pend_base[k];
    if (upd_pend_sidx[k]>=0) {
      upd_shv[upd_pend_sidx[k]]=upd_pend_sv[k];
      upd_shb[upd_pend_sidx[k]]=upd_pend_sb[k];
    }
  }
  upd_pend_n=0;
  for (c=0; c<upd_nchecks; c++) upd_range_check(coefs,upd_checks[c].index,upd_checks[c].offset,upd_checks[c].varsize,elem_vals);
  upd_nchecks=0;
}

/* after the last step of multi-step pass `pass` (0..2): fold the path
   values of the (change)/(explicit) targets into the extrapolation with
   that pass's Richardson weight (the weights the variables use) */
int updates_path_active(void) {
  return upd_pathn>0;
}

/* a multi-step pass restarts from the initial data: the shadow re-syncs
   from the stored values at the next Update pass */
void updates_path_restart(void) {
  upd_shvalid=0;
}

static void upd_shadow_sync(array_def *coefs, offset_t ncof, elem_value *elem_vals) {
  offset_t k,e;
  if (upd_pathn==0||upd_shvalid) return;
  if (upd_shv==NULL) {
    upd_shv=malloc(upd_pathn*sizeof(double));
    upd_shb=malloc(upd_pathn*sizeof(double));
    if (upd_shv==NULL||upd_shb==NULL) {
      errmsg("Error: out of memory for the updated-data extrapolation\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
  for (k=0; k<ncof; k++) {
    offset_t b=upd_pathbase[k];
    if (b<0) continue;
    for (e=0; e<coefs[k].nelem; e++) {
      upd_shv[b+e]=elem_vals[coefs[k].offset+e].value;
      upd_shb[b+e]=elem_vals[coefs[k].offset+e].substep_base;
    }
  }
  upd_shvalid=1;
}

void updates_path_accumulate(array_def *coefs, offset_t ncof, elem_value *elem_vals, double w, int first) {
  offset_t k,e;
  if (upd_pathbase==NULL||upd_pathn==0) return;
  if (upd_pathacc==NULL) upd_pathacc=calloc(upd_pathn,sizeof(double));
  for (k=0; k<ncof; k++) {
    offset_t b=upd_pathbase[k];
    if (b<0) continue;
    for (e=0; e<coefs[k].nelem; e++) {
      double v=w*(upd_shvalid?upd_shv[b+e]:(double)elem_vals[coefs[k].offset+e].value);
      upd_pathacc[b+e]=first?v:upd_pathacc[b+e]+v;
    }
  }
}

offset_t updates_apply(char *fname,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,int midpoint) {
  FILE * filehandle;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  offset_t i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops;//m,
  dim_t fdim,dcount,varsupsetid[MAXVARDIM];
  int nops=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx;
  bool IsChange=false,IsExplicit=false;
  solve_real zerodivide=0,temp1,temp2;
  offset_t varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  quantifier *arSet1=NULL;
  formula_op *ops1= NULL;
  strcpy(commsyntax,"update");
  /* division by zero is never allowed in UPDATEs (manual 10.11.1):
     the dual-class state must not leak in; legacy single default only */
  zdiv_disable();
  upd_analyse(fname,elem_vals,coefs,ncof);
  upd_shadow_sync(coefs,ncof,elem_vals);
  int stmt=0;
  filehandle = fopen(fname,"r");
  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    /* a mapped argument on the RHS lowers to map~idx and binds through
       the formula operand binder like any expression (manual 11.9.4;
       GTAP-E NCTAXLEV(r) = del_nctaxb(REGTOBLOC(r))); on the LHS it
       stays a named fatal (11.9.9) */
    mapping_reject_lhs(line,"Update");
    if (teems_nmap>0) mapping_lower_calls(line);
    /* update statements have no condition machinery: a ':' used to make
       the set lookup miss and expand over sets[0] in silence (M3) */
    if (strchr(line,':')!=NULL) {
      errmsg("Error: conditions in Update statements are not supported\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    IsChange=false;
    IsExplicit=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      str_replace_first(line, "(change)", "");
    }
    if(strstr(line, "(explicit)")!=NULL) {
      IsExplicit=true;
      str_replace_first(line, "(explicit)", "");
    }
    /* (product) is GEMPACK's explicit spelling of the default update
       form (parity plan 3.4) -- accept and ignore */
    if(strstr(line, "(product)")!=NULL)str_replace_first(line, "(product)", "");
    str_replace_first(line, commsyntax, "");
    str_delete_char(line,' ');
    while (str_replace_char(line, '[', '('));
    while (str_replace_char(line, ']', ')'));
    while (str_replace_char(line, '{', '('));
    while (str_replace_char(line, '}', ')'));
    strcpy(linecopy,line);
    readitem = strtok(line,"=");
    fdim=str_count_char(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    quantifier *arSet= (quantifier *) calloc (fdim+1,sizeof(quantifier));

    nloops=1;
    if (fdim==0) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
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
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-2]=1;
      for (i=fdim-3; i>-1; i--) {
        dcountdim1[i]=sets[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    strcpy(vname,readitem);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    if(IsChange==false&&IsExplicit==false) {
      while (str_replace_all(readitem,"*", "+"));
      strcpy(line1,vname);
      strcat(line1,"*(1+(");
      strcat(line1,readitem);
      strcat(line1,")/100)");
      readitem=line1;
    }
    if(IsChange==true) {
      strcpy(line1,vname);
      if(readitem[0]!='+'||readitem[0]!='-') strcat(line1,"+");
      strcat(line1,readitem);
      readitem=line1;
    }
    npow=str_count_char(readitem, '^');
    nmul=str_count_char(readitem, '*');
    ndiv=str_count_char(readitem, '/');
    nmul=nmul+ndiv;
    nplu=str_count_char(readitem, '+');
    nmin=str_count_char(readitem, '-');
    nplu=nplu+nmin;

    while (formula_normalize(readitem)==1);
    leadlag_encode(readitem);
    npar=str_count_char(readitem,'(')+str_count_char(readitem,',')+str_count_ci(readitem,"$pos");/* comma slack: multi-arg intrinsics emit per-arg loads + folds (plan 3.1); $POS calls: one position op each, their index lists ride in {} after normalize (fuzz batch 13) */
    strcpy(sumsyntax,"sum(");
    totalsum=sum_count(readitem,sumsyntax);
    sum_def *sum_cof= (sum_def *) calloc (totalsum,sizeof(sum_def));
    sumcount=0;
    /* readitem aliases line1 on the non-explicit branches above */
    if (readitem!=line1) strcpy(line1,readitem);
    strcpy(line2,line1);
    readitem=line2;
    while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim,sumcount)==1) {
      sumcount++;
    }
    totalsum=sumcount;
    i3=0;
    for (i=0; i<totalsum; i++) {
      i1=1;
      for(j=0; j<sum_cof[i].size; j++) {
        i1=i1*sets[sum_cof[i].setid[j]].size;
      }
      sum_cof[i].offset=i3;
      i3=i3+i1;
    }
    nsumele=i3;
    for (i=0; i<totalsum; i++) {
      i1=1;
      /* scalar sum store (size 0): no strides (formulas_execute
         guard; strides[-1] hit setid[MAXVARDIM-1], unread at
         size 0) */
      if (sum_cof[i].size==0) continue;
      sum_cof[i].strides[sum_cof[i].size-1]=1;
      for(j=sum_cof[i].size-2; j>-1; j--) {
        sum_cof[i].strides[j]=sum_cof[i].strides[j+1]*sets[sum_cof[i].setid[j+1]].size;
      }
    }
    formula_op *ops= (formula_op *) calloc (npow+nmul+nplu+2*npar+2,sizeof(formula_op));
    sum_value *sum_vals= (sum_value *) calloc (nsumele,sizeof(sum_value));
    sumcount=0;
    strcpy(line2,line1);
    readitem=line2;
    sumindx=0;
    while (sum_eval(readitem,sumsyntax,sets,nset,set_elems,elem_vals,ncofvar,ncofele,coefs,ncof,vars,nvar,sum_cof,totalsum,sum_vals,nsumele,ops,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
      sumcount++;
    }
    strcpy(line1,readitem);
    offset_t index=ncof-1, offset=0;//,simpl=0;
    bool check10=true;
        offset_t varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            offset=coefs[index].offset;
            varsize=coefs[index].size;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(vars[index].cofname,p)==0) {
              offset=ncofele+vars[index].offset;
              varsize=vars[index].size;
              if(vars[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<vars[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=vars[index].strides[dcount];
                  vararset[dcount]=l+1;
                  { dim_t ss=set_supset_slot(sets,arSet[l].setid,vars[index].setid[dcount]); if(ss<0)set_supset_fatal(p,vars[index].cofname,NULL,sets,arSet[l].setid,vars[index].setid[dcount]); if(ss>0){varsubset[dcount]=1; varsupsetid[dcount]=ss;} }
                  break;
                }
            }
        } else {
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=coefs[index].strides[dcount];
                  vararset[dcount]=l+1;
                  { dim_t ss=set_supset_slot(sets,arSet[l].setid,coefs[index].setid[dcount]); if(ss<0)set_supset_fatal(p,coefs[index].cofname,NULL,sets,arSet[l].setid,coefs[index].setid[dcount]); if(ss>0){varsubset[dcount]=1; varsupsetid[dcount]=ss;} }
                  break;
                }
            }
        }
    if(!formula_compile(line1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdim-1))MPI_Abort(PETSC_COMM_WORLD,1);
    int dfr=upd_stmt_defer(stmt);
    offset_t pb=dfr?upd_pend_reserve(nloops):0;
    offset_t shb0=((IsChange||IsExplicit)&&!check10&&upd_shvalid&&upd_pathbase[index]>=0)?upd_pathbase[index]:-1;
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,temp1,temp2,arSet1,ops1) shared(elem_vals,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=malloc((fdim+1)*sizeof(quantifier));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(quantifier));
          ops1=malloc(nops*sizeof(formula_op));
          memcpy (ops1,ops,nops*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet1=arSet;
        }
    /* workers must finish copying the shared scratch (master aliases it)
       before the master starts mutating per-iteration state */
    #pragma omp barrier
        #pragma omp for
    for (l=0; l<nloops; l++) {
      l2=0;
      i4=l;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(offset_t) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
      if(dfr){
        store_real v0=elem_vals[offset+l2].value,nv=v0;
        temp1=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
        if(midpoint==2)nv=0.5*(elem_vals[offset+l2].substep_base+temp1);
        else if(temp1-v0>0.000000001||temp1-v0<-0.000000001)nv=midpoint?elem_vals[offset+l2].substep_base+2*(temp1-v0):temp1;
        upd_pend_off[pb+l]=offset+l2;
        upd_pend_val[pb+l]=nv;
        upd_pend_base[pb+l]=v0;
        upd_pend_sidx[pb+l]=-1;
        if(shb0>=0){
          offset_t si=shb0+l2;
          double s0=upd_shv[si],sb=upd_shb[si],d=temp1-(double)v0,sv;
          if(IsExplicit) sv=(midpoint==2)?0.5*(sb+temp1):(midpoint?sb+2*(temp1-s0):temp1);
          else sv=(midpoint==2)?0.5*(sb+s0+d):(midpoint?sb+2*d:s0+d);
          upd_pend_sidx[pb+l]=si;
          upd_pend_sv[pb+l]=sv;
          upd_pend_sb[pb+l]=s0;
        }
      }else if(midpoint==2){
        /* Gragg terminal smoothing of the data: (C[n-1] + C[n] + dC)/2 */
        temp2=elem_vals[offset+l2].value;
        temp1=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
        elem_vals[offset+l2].value=0.5*(elem_vals[offset+l2].substep_base+temp1);
        elem_vals[offset+l2].substep_base=temp2;
      }else if(midpoint){
        temp2=elem_vals[offset+l2].value;
        temp1=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
        if(temp1-elem_vals[offset+l2].value>0.000000001||temp1-elem_vals[offset+l2].value<-0.000000001)elem_vals[offset+l2].value=elem_vals[offset+l2].substep_base+2*(temp1-elem_vals[offset+l2].value);
        elem_vals[offset+l2].substep_base=temp2;
      }else{
        elem_vals[offset+l2].substep_base=elem_vals[offset+l2].value;
        temp1=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
        if(temp1-elem_vals[offset+l2].value>0.000000001||temp1-elem_vals[offset+l2].value<-0.000000001)elem_vals[offset+l2].value=temp1;
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
    free(sum_cof);
    free(sum_vals);
    free(arSet);
    free(ops);

        if(dfr) upd_check_later(index,offset,varsize);
        else upd_range_check(coefs,index,offset,varsize,elem_vals);
        stmt++;
    
  }
  fclose(filehandle);
  upd_pass_flush(coefs,elem_vals);
  return j;
}



offset_t updates_apply_product(char *fname,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele) {
  FILE * filehandle;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  offset_t i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops;//m,
  dim_t fdim,dcount,varsupsetid[MAXVARDIM];
  int nops=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx;
  bool IsChange=false,IsExplicit=false;
  solve_real zerodivide=0,temp1;
  offset_t varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  quantifier *arSet1=NULL;
  formula_op *ops1= NULL;
  strcpy(commsyntax,"update");
  /* division by zero is never allowed in UPDATEs (manual 10.11.1):
     the dual-class state must not leak in; legacy single default only */
  zdiv_disable();
  upd_analyse(fname,elem_vals,coefs,ncof);
  int stmt=0;
  filehandle = fopen(fname,"r");
  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    /* a mapped argument on the RHS lowers to map~idx and binds through
       the formula operand binder like any expression (manual 11.9.4;
       GTAP-E NCTAXLEV(r) = del_nctaxb(REGTOBLOC(r))); on the LHS it
       stays a named fatal (11.9.9) */
    mapping_reject_lhs(line,"Update");
    if (teems_nmap>0) mapping_lower_calls(line);
    /* update statements have no condition machinery: a ':' used to make
       the set lookup miss and expand over sets[0] in silence (M3) */
    if (strchr(line,':')!=NULL) {
      errmsg("Error: conditions in Update statements are not supported\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    IsChange=false;
    IsExplicit=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      str_replace_first(line, "(change)", "");
    }
    if(strstr(line, "(explicit)")!=NULL) {
      IsExplicit=true;
      str_replace_first(line, "(explicit)", "");
    }
    /* (product) is GEMPACK's explicit spelling of the default update
       form (parity plan 3.4) -- accept and ignore */
    if(strstr(line, "(product)")!=NULL)str_replace_first(line, "(product)", "");
    str_replace_first(line, commsyntax, "");
    str_delete_char(line,' ');
    while (str_replace_char(line, '[', '('));
    while (str_replace_char(line, ']', ')'));
    while (str_replace_char(line, '{', '('));
    while (str_replace_char(line, '}', ')'));
    strcpy(linecopy,line);
    readitem = strtok(line,"=");
    fdim=str_count_char(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    quantifier *arSet= (quantifier *) calloc (fdim+1,sizeof(quantifier));

    nloops=1;
    if (fdim==0) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
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
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-2]=1;
      for (i=fdim-3; i>-1; i--) {
        dcountdim1[i]=sets[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    strcpy(vname,readitem);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    if(IsChange==false&&IsExplicit==false) {
      while (str_replace_all(readitem,"*", "/100)!(1+"));
      while (str_replace_all(readitem,"!", "*"));
      strcpy(line1,vname);
      strcat(line1,"*(1+");
      strcat(line1,readitem);
      strcat(line1,"/100)");
      readitem=line1;
    }
    if(IsChange==true) {
      strcpy(line1,vname);
      if(readitem[0]!='+'||readitem[0]!='-') strcat(line1,"+");
      strcat(line1,readitem);
      readitem=line1;
    }
    npow=str_count_char(readitem, '^');
    nmul=str_count_char(readitem, '*');
    ndiv=str_count_char(readitem, '/');
    nmul=nmul+ndiv;
    nplu=str_count_char(readitem, '+');
    nmin=str_count_char(readitem, '-');
    nplu=nplu+nmin;

    while (formula_normalize(readitem)==1);
    leadlag_encode(readitem);
    npar=str_count_char(readitem,'(')+str_count_char(readitem,',')+str_count_ci(readitem,"$pos");/* comma slack: multi-arg intrinsics emit per-arg loads + folds (plan 3.1); $POS calls: one position op each, their index lists ride in {} after normalize (fuzz batch 13) */
    strcpy(sumsyntax,"sum(");
    totalsum=sum_count(readitem,sumsyntax);
    sum_def *sum_cof= (sum_def *) calloc (totalsum,sizeof(sum_def));
    sumcount=0;
    /* readitem aliases line1 on the non-explicit branches above */
    if (readitem!=line1) strcpy(line1,readitem);
    strcpy(line2,line1);
    readitem=line2;
    while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim,sumcount)==1) {
      sumcount++;
    }
    totalsum=sumcount;
    i3=0;
    for (i=0; i<totalsum; i++) {
      i1=1;
      for(j=0; j<sum_cof[i].size; j++) {
        i1=i1*sets[sum_cof[i].setid[j]].size;
      }
      sum_cof[i].offset=i3;
      i3=i3+i1;
    }
    nsumele=i3;
    for (i=0; i<totalsum; i++) {
      i1=1;
      /* scalar sum store (size 0): no strides (formulas_execute
         guard; strides[-1] hit setid[MAXVARDIM-1], unread at
         size 0) */
      if (sum_cof[i].size==0) continue;
      sum_cof[i].strides[sum_cof[i].size-1]=1;
      for(j=sum_cof[i].size-2; j>-1; j--) {
        sum_cof[i].strides[j]=sum_cof[i].strides[j+1]*sets[sum_cof[i].setid[j+1]].size;
      }
    }
    formula_op *ops= (formula_op *) calloc (npow+nmul+nplu+2*npar+2,sizeof(formula_op));
    sum_value *sum_vals= (sum_value *) calloc (nsumele,sizeof(sum_value));
    sumcount=0;
    strcpy(line2,line1);
    readitem=line2;
    sumindx=0;
    while (sum_eval(readitem,sumsyntax,sets,nset,set_elems,elem_vals,ncofvar,ncofele,coefs,ncof,vars,nvar,sum_cof,totalsum,sum_vals,nsumele,ops,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
      sumcount++;
    }
    strcpy(line1,readitem);
    offset_t index=ncof-1, offset=0;//,simpl=0;
    bool check10=true;
        offset_t varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            offset=coefs[index].offset;
            varsize=coefs[index].size;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(vars[index].cofname,p)==0) {
              offset=ncofele+vars[index].offset;
              varsize=vars[index].size;
              if(vars[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<vars[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=vars[index].strides[dcount];
                  vararset[dcount]=l+1;
                  { dim_t ss=set_supset_slot(sets,arSet[l].setid,vars[index].setid[dcount]); if(ss<0)set_supset_fatal(p,vars[index].cofname,NULL,sets,arSet[l].setid,vars[index].setid[dcount]); if(ss>0){varsubset[dcount]=1; varsupsetid[dcount]=ss;} }
                  break;
                }
            }
        } else {
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=coefs[index].strides[dcount];
                  vararset[dcount]=l+1;
                  { dim_t ss=set_supset_slot(sets,arSet[l].setid,coefs[index].setid[dcount]); if(ss<0)set_supset_fatal(p,coefs[index].cofname,NULL,sets,arSet[l].setid,coefs[index].setid[dcount]); if(ss>0){varsubset[dcount]=1; varsupsetid[dcount]=ss;} }
                  break;
                }
            }
        }
    if(!formula_compile(line1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdim-1))MPI_Abort(PETSC_COMM_WORLD,1);
    int dfr=upd_stmt_defer(stmt);
    offset_t pb=dfr?upd_pend_reserve(nloops):0;
    /* the post-simulation pass of an extrapolated solve: a (change) or
       (explicit) target takes its extrapolated path value; product
       updates stay the exact one-shot form */
    const double *pathv=(teems_upd_pathuse&&(IsChange||IsExplicit)&&!check10&&upd_pathacc!=NULL&&upd_pathbase[index]>=0)?upd_pathacc+upd_pathbase[index]:NULL;
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,temp1,arSet1,ops1) shared(elem_vals,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=malloc((fdim+1)*sizeof(quantifier));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(quantifier));
          ops1=malloc(nops*sizeof(formula_op));
          memcpy (ops1,ops,nops*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet1=arSet;
        }
    /* workers must finish copying the shared scratch (master aliases it)
       before the master starts mutating per-iteration state */
    #pragma omp barrier
        #pragma omp for
    for (l=0; l<nloops; l++) {
      l2=0;
      i4=l;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(offset_t) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
      temp1=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
      if(pathv!=NULL)temp1=pathv[l2];
      if(dfr){
        store_real v0=elem_vals[offset+l2].value,nv=v0;
        if(temp1-v0>0.000000001||temp1-v0<-0.000000001)nv=temp1;
        upd_pend_off[pb+l]=offset+l2;
        upd_pend_val[pb+l]=nv;
        upd_pend_base[pb+l]=elem_vals[offset+l2].substep_base;
        upd_pend_sidx[pb+l]=-1;
      }else
      if(temp1-elem_vals[offset+l2].value>0.000000001||temp1-elem_vals[offset+l2].value<-0.000000001)elem_vals[offset+l2].value=temp1;
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
    free(sum_cof);
    free(sum_vals);
    free(arSet);
    free(ops);

        if(dfr) upd_check_later(index,offset,varsize);
        else upd_range_check(coefs,index,offset,varsize,elem_vals);
        stmt++;
    
  }
  fclose(filehandle);
  upd_pass_flush(coefs,elem_vals);
  return j;
}

int sum_eval(char *formulain, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele, array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,sum_def *sum_cof,int totalsum,sum_value *sum_vals,offset_t nsumele,formula_op *ops,quantifier *arSet1,dim_t fdim,int *sumindx,int j, solve_real zerodivide) {
  char *readitem,*p;//,*p1,interchar2[NAMESIZE],line5[TABREADLINE];
  char interchar[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE];//,line3[TABREADLINE],line4[TABREADLINE];//,interchar1[NAMESIZE]
  int nops,length,k=0,k1=0,i=0;
  dim_t fdimsumcof,superset_pos,dcount;
  offset_t l,l1,l2,nloops,dcountdim1[4*MAXVARDIM];
  solve_real vval;
  quantifier *arSet2=NULL;
  formula_op *ops1= NULL;
  offset_t arsetsize;
  int condpos=-1;
  sum_cofcond cofcond;
  cofcond.cofid=-1;
  offset_t condfix=-1;
  dim_t condss=0;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=str_find_token_ci(formulain,readitem,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      sum_extract(line);
      k1=str_find_token_ci(line,line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
        strcpy(line1,line);
        /* depth-aware body extraction: a coefficient condition may
           carry commas the old comma-splits broke on (gap 2) */
        p=sum_body_extract(line1);
        if (p==NULL) {
          errmsg("Error: malformed sum statement '%s'\n",line1);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
        strcpy(line2,p);
        arsetsize=sum_cof[j].size+1;
        quantifier *arSet= (quantifier *) calloc (arsetsize,sizeof(quantifier));
        for (l=0; l<sum_cof[j].size; l++) {
          arSet[l].setid=sum_cof[j].setid[l];
          strcpy(arSet[l].index_name,sum_cof[j].dimnames[l]);
        }
        nloops=1;
        for (l=0; l<sum_cof[j].size; l++) {
          nloops=nloops*sets[arSet[l].setid].size;
          dcount=sum_cof[j].size-l;
          if(dcount==sum_cof[j].size) {
            dcountdim1[dcount-1]=1;
          } else {
            dcountdim1[dcount-1]=dcountdim1[dcount]*sets[arSet[dcount].setid].size;
          }
        }
        arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
        strcpy(arSet[sum_cof[j].size].index_name,sum_cof[j].sumindx);
        fdimsumcof=sum_cof[j].size+1;
        sum_cond_rhs_resolve(sum_cof[j].cond_mapid,sum_cof[j].cond_rhs,arSet,fdimsumcof,sets,set_elems,&condpos,&condfix,&condss);
        sum_cond_coef_resolve(&sum_cof[j],arSet,fdimsumcof,sets,set_elems,coefs,ncof,&cofcond);
        nops=0;
        if(!formula_compile(p,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdimsumcof))MPI_Abort(PETSC_COMM_WORLD,1);
        #pragma omp parallel private(l,l1,l2,dcount,superset_pos,vval,arSet2,ops1) shared(elem_vals,arSet,sum_vals)
        {
        if(omp_get_thread_num()!=0){
          arSet2=malloc(arsetsize*sizeof(quantifier));
          memcpy(arSet2,arSet,arsetsize*sizeof(quantifier));
          ops1=malloc(nops*sizeof(formula_op));
          memcpy(ops1,ops,nops*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet2=arSet;
        }
    /* workers must finish copying the shared scratch (master aliases it)
       before the master starts mutating per-iteration state */
    #pragma omp barrier
        #pragma omp for
        for (l=0; l<nloops; l++) {
          l2=l;
          for (dcount=0; dcount<sum_cof[j].size; dcount++) {
            superset_pos=(dim_t) l2/dcountdim1[dcount];
            arSet2[dcount].indx=superset_pos;
            l2=l2-superset_pos*dcountdim1[dcount];
          }
          vval=0;
          for (l1=0; l1<sets[sum_cof[j].sumsetid].size; l1++) {
            /* mapping-equality condition (M3): only domain elements
               mapping to the target codomain position contribute */
            if (sum_cof[j].cond_mapid>0&&(offset_t)teems_maps[sum_cof[j].cond_mapid-1].values[l1]!=sum_cond_target(condpos,condss,condfix,arSet2,sets,set_elems)) continue;
            /* coefficient-comparison condition (11.4.11; IF-survey gap 2) */
            if (cofcond.cofid>=0&&!sum_cofcond_test(&cofcond,elem_vals,arSet2,l1)) continue;
            arSet2[sum_cof[j].size].indx=l1;
            vval+=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet2,fdimsumcof,zerodivide);
          }
          sum_vals[*sumindx+l].value=vval;
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
        *sumindx=*sumindx+nloops;
        strcpy(interchar,sum_cof[j].sumname);
        strcat(interchar,"{");
        for (l=0; l<sum_cof[j].size; l++) {
          strcat(interchar,sum_cof[j].dimnames[l]);
          strcat(interchar,",");
        }
        if (interchar[strlen(interchar)-1]==',') {
          interchar[strlen(interchar)-1]='}';
        } else {
          if (interchar[strlen(interchar)-1]=='{') {
            interchar[strlen(interchar)-1]='\0';
          } else {
            strcat(interchar,"}");
          }
        }
        while(str_replace_all(formulain,line1,interchar)!=NULL);
        free(arSet);
        return 1;
      }
    } else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      sum_extract(line);
      k1=str_find_token_ci(line,line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
        strcpy(line1,line);
        /* depth-aware body extraction: a coefficient condition may
           carry commas the old comma-splits broke on (gap 2) */
        p=sum_body_extract(line1);
        if (p==NULL) {
          errmsg("Error: malformed sum statement '%s'\n",line1);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
        arsetsize=sum_cof[j].size+1;
        quantifier *arSet= (quantifier *) calloc (arsetsize,sizeof(quantifier));
        for (l=0; l<sum_cof[j].size; l++) {
          arSet[l].setid=sum_cof[j].setid[l];
          strcpy(arSet[l].index_name,sum_cof[j].dimnames[l]);
        }
        nloops=1;
        for (l=0; l<sum_cof[j].size; l++) {
          nloops=nloops*sets[arSet[l].setid].size;//sum_cof[j].dims[l];
          dcount=sum_cof[j].size-l;
          if(dcount==sum_cof[j].size) {
            dcountdim1[dcount-1]=1;
          } else {
            dcountdim1[dcount-1]=dcountdim1[dcount]*sets[arSet[dcount].setid].size;
          }
        }
        arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
        strcpy(arSet[sum_cof[j].size].index_name,sum_cof[j].sumindx);
        fdimsumcof=sum_cof[j].size+1;
        sum_cond_rhs_resolve(sum_cof[j].cond_mapid,sum_cof[j].cond_rhs,arSet,fdimsumcof,sets,set_elems,&condpos,&condfix,&condss);
        sum_cond_coef_resolve(&sum_cof[j],arSet,fdimsumcof,sets,set_elems,coefs,ncof,&cofcond);
        nops=0;
        if(!formula_compile(p,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdimsumcof))MPI_Abort(PETSC_COMM_WORLD,1);
        #pragma omp parallel private(l,l1,l2,dcount,superset_pos,vval,arSet2,ops1) shared(elem_vals,arSet,sum_vals)
        {
        if(omp_get_thread_num()!=0){
          arSet2=malloc(arsetsize*sizeof(quantifier));
          memcpy(arSet2,arSet,arsetsize*sizeof(quantifier));
          ops1=malloc(nops*sizeof(formula_op));
          memcpy(ops1,ops,nops*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet2=arSet;
        }
    /* workers must finish copying the shared scratch (master aliases it)
       before the master starts mutating per-iteration state */
    #pragma omp barrier
        #pragma omp for
        for (l=0; l<nloops; l++) {
          l2=l;
          for (dcount=0; dcount<sum_cof[j].size; dcount++) {
            superset_pos=(dim_t) l2/dcountdim1[dcount];
            arSet2[dcount].indx=superset_pos;
            l2=l2-superset_pos*dcountdim1[dcount];
          }
          vval=0;
          for (l1=0; l1<sets[sum_cof[j].sumsetid].size; l1++) {
            /* mapping-equality condition (M3): only domain elements
               mapping to the target codomain position contribute */
            if (sum_cof[j].cond_mapid>0&&(offset_t)teems_maps[sum_cof[j].cond_mapid-1].values[l1]!=sum_cond_target(condpos,condss,condfix,arSet2,sets,set_elems)) continue;
            /* coefficient-comparison condition (11.4.11; IF-survey gap 2) */
            if (cofcond.cofid>=0&&!sum_cofcond_test(&cofcond,elem_vals,arSet2,l1)) continue;
            arSet2[sum_cof[j].size].indx=l1;
            vval+=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet2,fdimsumcof,zerodivide);
          }
          sum_vals[*sumindx+l].value=vval;//ha_sumele[*sumindx+l2].varval=vval;
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
        *sumindx=*sumindx+nloops;
        strcpy(interchar,sum_cof[j].sumname);
        strcat(interchar,"{");
        for (l=0; l<sum_cof[j].size; l++) {
          strcat(interchar,sum_cof[j].dimnames[l]);
          strcat(interchar,",");
        }
        if (interchar[strlen(interchar)-1]==',') {
          interchar[strlen(interchar)-1]='}';
        } else {
          if (interchar[strlen(interchar)-1]=='{') {
            interchar[strlen(interchar)-1]='\0';
          } else {
            strcat(interchar,"}");
          }
        }
        while(str_replace_all(formulain,line1,interchar));
        free(arSet);
        return 1;
      }
    } else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

/* ASSERTION evaluation (GEMPACK manual 10.14 / 25.3; parity plan 1.1).
   The condition "lhs RELOP rhs" is evaluated as the residual lhs-(rhs)
   through the ordinary formula engine, per quantifier tuple, against
   the current coefficient values -- so an (always) assertion rides
   every formulas_execute pass and an (initial) one only the first.
   mode (-assertions run switch, default fatal): 0 = skip,
   1 = warn and continue, 2 = report failing elements and abort.
   Unsupported condition forms (conditional quantifiers, functions the
   compiler lacks) warn and skip the assertion rather than turning a
   previously-running model into an abort. */
offset_t assertions_execute(char *fname,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,bool IsIni,int mode,int postsim_pass) {
  FILE *filehandle;
  char line[TABREADLINE],linecopy[TABREADLINE],resid[TABREADLINE],msg[TABREADLINE];
  char sumsyntax[NAMESIZE],tempset[NAMESIZE];
  char *p=NULL,*q=NULL,*right=NULL;
  static offset_t nassert=-1;
  static char *nassert_fname=NULL;
  offset_t i,l,i3,i4,nloops,len,dcountdim1[4*MAXVARDIM],total_fail=0;
  dim_t nq,dcount;
  int nops,npow,nmul,nplu,npar,relop,fails,shown,totalsum,sumcount,sumindx,depth;
  offset_t nsumele;
  solve_real zerodivide=0,r;
  bool IsAssIni,ok,skip;
  /* conditional quantifiers "(all,i,S: expr RELOP expr)" (manual
     10.14 / 11.4.11): each condition is compiled as the residual
     lhs-(rhs) through the formula engine, like the assertion itself,
     and evaluated per tuple; a false condition excludes the tuple
     (ORANI-G's zero-flow margin and tax checks, 2026-09-24) */
  formula_op *cond_ops[4*MAXVARDIM];
  dim_t cond_nops[4*MAXVARDIM];
  int cond_relop[4*MAXVARDIM],ncond;
  char condtxt[TABREADLINE];
  if(mode==0)return 0;
  /* count cached per file buffer: the PostSim pass runs this on the
     _ps companion with its own count */
  if(nassert<0||fname!=nassert_fname) {
    nassert=tab_count_statements(fname,"assertion");
    nassert_fname=fname;
  }
  if(nassert==0)return 0;
  filehandle=fopen(fname,"r");
  if(filehandle==NULL)return 0;
  /* fresh positional zerodivide state for this rescan (manual 10.11.1) */
  zdiv_scan_reset();
  strcpy(sumsyntax,"sum(");
  while (tab_next_statement_resolved("assertion",filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    mapping_reject_in(line,"Assertion");
    /* (postsim) assertions run only in the post-solve pass, where the
       initial/always qualifiers are ignored (manual 12.2.4) */
    if(strstr(line,"(postsim)")!=NULL) {
      if(postsim_pass==0)continue;
      str_replace_first(line,"(postsim)","");
    }
    else if(postsim_pass==1)continue;
    IsAssIni=false;
    if(strstr(line,"(initial)")!=NULL) {
      str_replace_first(line,"(initial)","");
      IsAssIni=true;
    }
    if(strstr(line,"(always)")!=NULL)str_replace_first(line,"(always)","");
    if(postsim_pass==0&&IsAssIni&&!IsIni)continue;
    /* assertions are formula-class statements: the dual-class
       zerodivide state as of this file position applies (plan A1) */
    zdiv_capture();
    /* optional # message #, captured before whitespace stripping */
    msg[0]='\0';
    p=strchr(line,'#');
    if(p!=NULL) {
      q=strchr(p+1,'#');
      if(q!=NULL) {
        len=q-p-1;
        if(len>=TABREADLINE)len=TABREADLINE-1;
        strncpy(msg,p+1,len);
        msg[len]='\0';
        while(len>0&&msg[len-1]==' ') msg[--len]='\0';
        while(msg[0]==' ') memmove(msg,msg+1,strlen(msg));
        memmove(p,q+1,strlen(q+1)+1);
      }
    }
    /* word comparison operators need their delimiting spaces */
    tab_wordops_normalize(line);
    str_replace_first(line,"assertion","");
    str_delete_char(line,' ');
    while (str_replace_char(line, '[', '('));
    while (str_replace_char(line, ']', ')'));
    while (str_replace_char(line, '{', '('));
    while (str_replace_char(line, '}', ')'));
    strcpy(linecopy,line);
    /* leading (all,index,SET) quantifiers */
    quantifier *arSet= (quantifier *) calloc (4*MAXVARDIM+1,sizeof(quantifier));
    nq=0;
    nloops=1;
    skip=false;
    ncond=0;
    p=linecopy;
    while(strncmp(p,"(all,",5)==0&&nq<4*MAXVARDIM) {
      p+=5;
      q=strchr(p,',');
      if(q==NULL) {
        skip=true;
        break;
      }
      len=q-p;
      if(len>=NAMESIZE)len=NAMESIZE-1;
      strncpy(arSet[nq].index_name,p,len);
      arSet[nq].index_name[len]='\0';
      p=q+1;
      /* the quantifier's own ')' : a condition may carry argument
         lists, "(all,i,IND:V1BAS(c,s,i)=0)" */
      depth=1;
      for(q=p; *q!='\0'; q++) { if(*q=='(')depth++; else if(*q==')'&&--depth==0)break; }
      if(*q=='\0') {
        skip=true;
        break;
      }
      len=q-p;
      if(len>=NAMESIZE)len=NAMESIZE-1;
      strncpy(tempset,p,len);
      tempset[len]='\0';
      cond_ops[nq]=NULL;
      if(strchr(tempset,':')!=NULL) {
        char *c=strchr(tempset,':');
        *c='\0';
        c++;
        if(strlen(c)+1>=sizeof(condtxt)||strstr(c,"sum(")!=NULL||strstr(c,"$pos")!=NULL) {
          printf("Warning: assertion quantifier condition not supported (sums and $POS) -- assertion skipped: %s\n",linecopy);
          skip=true;
          break;
        }
        strcpy(condtxt,c);
        /* top-level comparison operator of the condition */
        {
          char *cq,*cright;
          int cdepth=0,crel=0;
          for(cq=condtxt; *cq!='\0'; cq++) {
            if(*cq=='(')cdepth++;
            else if(*cq==')')cdepth--;
            else if(cdepth==0) {
              if(*cq=='<'&&*(cq+1)=='>') { crel=2; break; }
              if(*cq=='>'&&*(cq+1)=='=') { crel=3; break; }
              if(*cq=='<'&&*(cq+1)=='=') { crel=4; break; }
              if(*cq=='>') { crel=5; break; }
              if(*cq=='<') { crel=6; break; }
              if(*cq=='=') { crel=1; break; }
            }
          }
          if(crel==0) {
            printf("Warning: assertion quantifier condition has no comparison operator -- assertion skipped: %s\n",linecopy);
            skip=true;
            break;
          }
          cright=cq+((crel==1||crel==5||crel==6)?1:2);
          *cq='\0';
          if(snprintf(resid,TABREADLINE,"%s-(%s)",condtxt,cright)>=TABREADLINE-1) {
            printf("Warning: assertion quantifier condition too long -- skipped\n");
            skip=true;
            break;
          }
          while (formula_normalize(resid)==1);
          leadlag_encode(resid);
          {
            int cnpow=str_count_char(resid,'^'),cnmul=str_count_char(resid,'*')+str_count_char(resid,'/');
            int cnplu=str_count_char(resid,'+')+str_count_char(resid,'-'),cnpar=str_count_char(resid,'(')+str_count_char(resid,',');
            cond_ops[nq]=(formula_op *) calloc (cnpow+cnmul+cnplu+2*(cnpar+2),sizeof(formula_op));
          }
          cond_nops[nq]=0;
          cond_relop[nq]=crel;
          /* the condition may only use quantifiers declared so far
             (this one included): compile against the frame after this
             quantifier's set is resolved, below */
          ncond++;
        }
      }
      for(i=0; i<nset; i++)if(strcmp(sets[i].setname,tempset)==0)break;
      if(i>=nset) {
        printf("Warning: assertion references unknown set '%s' -- assertion skipped\n",tempset);
        skip=true;
        break;
      }
      arSet[nq].setid=i;
      nloops=nloops*sets[i].size;
      nq++;
      if(cond_ops[nq-1]!=NULL) {
        if(!formula_compile(resid,sets,coefs,ncof,vars,nvar,ncofele,NULL,0,cond_ops[nq-1],&cond_nops[nq-1],arSet,nq)) {
          printf("Warning: assertion quantifier condition could not be compiled -- assertion skipped: %s\n",linecopy);
          skip=true;
          break;
        }
      }
      p=q+1;
    }
    if(skip) {
      for(dcount=0; dcount<nq; dcount++)if(cond_ops[dcount]!=NULL)free(cond_ops[dcount]);
      if(nq<4*MAXVARDIM&&cond_ops[nq]!=NULL)free(cond_ops[nq]);
      free(arSet);
      continue;
    }
    q=strchr(p,';');
    if(q!=NULL)*q='\0';
    /* top-level comparison operator: 1 = 2 <> 3 >= 4 <= 5 > 6 < */
    relop=0;
    depth=0;
    for(q=p; *q!='\0'; q++) {
      if(*q=='(')depth++;
      else if(*q==')')depth--;
      else if(depth==0) {
        if(*q=='<'&&*(q+1)=='>') {
          relop=2;
          break;
        }
        if(*q=='>'&&*(q+1)=='=') {
          relop=3;
          break;
        }
        if(*q=='<'&&*(q+1)=='=') {
          relop=4;
          break;
        }
        if(*q=='>') {
          relop=5;
          break;
        }
        if(*q=='<') {
          relop=6;
          break;
        }
        if(*q=='=') {
          relop=1;
          break;
        }
      }
    }
    if(relop==0) {
      printf("Warning: assertion has no comparison operator -- skipped: %s\n",linecopy);
      free(arSet);
      continue;
    }
    right=q+((relop==1||relop==5||relop==6)?1:2);
    *q='\0';
    /* strlen(p)+strlen(right)+3 chars plus the terminator must fit TABREADLINE; the
       return value carries the untruncated length, so a skip on >= is the old guard */
    if(snprintf(resid,TABREADLINE,"%s-(%s)",p,right)>=TABREADLINE-1) {
      printf("Warning: assertion condition too long -- skipped\n");
      free(arSet);
      continue;
    }
    while (formula_normalize(resid)==1);
    leadlag_encode(resid);
    npow=str_count_char(resid,'^');
    nmul=str_count_char(resid,'*')+str_count_char(resid,'/');
    nplu=str_count_char(resid,'+')+str_count_char(resid,'-');
    npar=str_count_char(resid,'(')+str_count_char(resid,',')+str_count_ci(resid,"$pos");/* comma slack for multi-arg intrinsics (plan 3.1); $POS calls: one position op each, their index lists ride in {} after normalize (fuzz batch 13) */
    totalsum=sum_count(resid,sumsyntax);
    sum_def *sum_cof= (sum_def *) calloc (totalsum+1,sizeof(sum_def));
    sumcount=0;
    strcpy(line,resid);
    p=line;
    while (sum_parse(p,sumsyntax,sum_cof,arSet,sets,nset,nq+1,sumcount)==1)sumcount++;
    totalsum=sumcount;
    i3=0;
    for(i=0; i<totalsum; i++) {
      i4=1;
      for(l=0; l<sum_cof[i].size; l++)i4=i4*sets[sum_cof[i].setid[l]].size;
      sum_cof[i].offset=i3;
      sum_cof[i].summatsize=i4;
      i3=i3+i4;
      /* a scalar sum store (size 0) has no strides; the unguarded
         write was strides[-1] = setid[MAXVARDIM-1] (the class fixed
         in formulas_execute and the equation builder; this is the
         assertion-executor twin, surfaced by the C1 kit's scalar
         sum assertions) */
      if (sum_cof[i].size==0) continue;
      sum_cof[i].strides[sum_cof[i].size-1]=1;
      for(l=sum_cof[i].size-2; l>-1; l--)sum_cof[i].strides[l]=sum_cof[i].strides[l+1]*sets[sum_cof[i].setid[l+1]].size;
    }
    nsumele=i3;
    formula_op *ops= (formula_op *) calloc (npow+nmul+nplu+2*(npar+2),sizeof(formula_op));
    sum_value *sum_vals= (sum_value *) calloc (nsumele+1,sizeof(sum_value));
    sumcount=0;
    sumindx=0;
    strcpy(line,resid);
    p=line;
    while (sum_eval(p,sumsyntax,sets,nset,set_elems,elem_vals,ncofvar,ncofele,coefs,ncof,vars,nvar,sum_cof,totalsum,sum_vals,nsumele,ops,arSet,nq+1,&sumindx,sumcount,zerodivide)==1)sumcount++;
    nops=0;
    if(!formula_compile(p,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,nq)) {
      printf("Warning: assertion condition could not be compiled -- assertion skipped: %s\n",linecopy);
      for(dcount=0; dcount<nq; dcount++)if(cond_ops[dcount]!=NULL)free(cond_ops[dcount]);
      free(arSet);
      free(sum_cof);
      free(sum_vals);
      free(ops);
      continue;
    }
    if(nq>0) {
      dcountdim1[nq-1]=1;
      for(i=nq-2; i>-1; i--)dcountdim1[i]=sets[arSet[i+1].setid].size*dcountdim1[i+1];
    }
    fails=0;
    shown=0;
    for(l=0; l<nloops; l++) {
      i4=l;
      for(dcount=0; dcount<nq; dcount++) {
        i3=(offset_t) i4/dcountdim1[dcount];
        arSet[dcount].indx=i3;
        i4=i4-i3*dcountdim1[dcount];
      }
      if(ncond>0) {
        bool cond_ok=true;
        for(dcount=0; dcount<nq&&cond_ok; dcount++) {
          solve_real cr;
          if(cond_ops[dcount]==NULL)continue;
          cr=formula_eval(elem_vals,sets,set_elems,NULL,cond_ops[dcount],cond_nops[dcount],arSet,dcount+1,zerodivide);
          if(cond_relop[dcount]==1)cond_ok=(cr==0);
          if(cond_relop[dcount]==2)cond_ok=(cr!=0);
          if(cond_relop[dcount]==3)cond_ok=(cr>=0);
          if(cond_relop[dcount]==4)cond_ok=(cr<=0);
          if(cond_relop[dcount]==5)cond_ok=(cr>0);
          if(cond_relop[dcount]==6)cond_ok=(cr<0);
        }
        if(!cond_ok)continue;
      }
      r=formula_eval(elem_vals,sets,set_elems,sum_vals,ops,nops,arSet,nq,zerodivide);
      ok=true;
      if(relop==1)ok=(r==0);
      if(relop==2)ok=(r!=0);
      if(relop==3)ok=(r>=0);
      if(relop==4)ok=(r<=0);
      if(relop==5)ok=(r>0);
      if(relop==6)ok=(r<0);
      if(!ok) {
        fails++;
        if(shown<10) {
          printf("%%%% Assertion '%s' does not hold",(msg[0]!='\0')?msg:linecopy);
          for(dcount=0; dcount<nq; dcount++)printf(" (quantifier number %d is '%s')",(int)(dcount+1),set_elems[sets[arSet[dcount].setid].offset+arSet[dcount].indx].setele);
          printf("\n");
          shown++;
        }
      }
    }
    if(fails>0) {
      if(fails>shown)printf("%%%% (%d further failing elements not listed)\n",fails-shown);
      printf("Assertion '%s' does not hold.\n",(msg[0]!='\0')?msg:linecopy);
      total_fail+=fails;
      if(teems_rk_stage_checks) {
        teems_check_viol_assert++;
        printf("Warning: assertion failed at a Runge-Kutta stage state; the step will be retried\n");
      }
      else if(mode==2) {
        errmsg("Error: assertion failed (-assertions 1 downgrades this abort to a warning, -assertions 0 skips the checks)\n");
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
    }
    for(dcount=0; dcount<nq; dcount++)if(cond_ops[dcount]!=NULL)free(cond_ops[dcount]);
    free(arSet);
    free(sum_cof);
    free(sum_vals);
    free(ops);
  }
  fclose(filehandle);
  zdiv_disable();
  if(mode==1&&total_fail>0)printf("Warning: %ld assertion failure%s in this pass (-assertions 1)\n",(long)total_fail,(total_fail==1)?"":"s");
  return total_fail;
}

/* PostSim foundation F3: expose the composed solution to the formula
   engine -- copy xcf (per variable element: the %-change / change
   results, exogenous shocks included) into the variables' elem_vals
   slots, so post-solve statements read simulation results as if the
   variables were coefficients (manual ch.12). */
void postsim_expose_results(elem_value *elem_vals,offset_t ncofele,offset_t nvarele,solve_real *xcf) {
  offset_t j5;
  /* both slots: a bare name binds .value, a p_/c_ name (the linear
     variable of a levels variable, p_V1BAS) binds the change slot
     .substep_base, which the final update pass had zeroed -- p_X read 0
     in PostSim (potential-models vetting, ADJUST) */
  for(j5=0; j5<nvarele; j5++) {
    elem_vals[ncofele+j5].value=(store_real)xcf[j5];
    elem_vals[ncofele+j5].substep_base=(store_real)xcf[j5];
  }
}

/* does the TAB carry any (postsim) assertions? (one scan, post-solve) */
int tab_has_postsim_assertions(char *fname) {
  FILE *filehandle;
  char line[TABREADLINE];
  int found=0;
  filehandle=fopen(fname,"r");
  if(filehandle==NULL)return 0;
  while (tab_next_statement("assertion",filehandle,line,TABREADLINE)) {
    if(strstr(line,"(postsim)")!=NULL) {
      found=1;
      break;
    }
  }
  fclose(filehandle);
  return found;
}
