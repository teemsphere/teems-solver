// Pham Van Ha & Tom Kompas
// Ensure that sum over time involves only intertemporal variables
// If otherwise switch back to normal MC66
// Currently no sub-reg set, can try different varreg to get opt: comment break in nestedmatvarread
// Currently logic fdim <> only in formula
#include <teems_solver.h>

static char help[] = "Solves a CGE model in parallel with KSP.\n\
           Input parameters include:\n\
           -None at the moment\n\n";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args) {
  Vec      vecb,vece,x;  /* approx solution, RHS, exact solution */
  Mat      A,B;    /* linear system matrix */
  KSP      ksp;   /* linear solver context */
  PC       pc;   /* preconditionercontext */
  PetscRandom  rctx;   /* random number generator context */
  PetscInt  rank=0,mpisize,rank_hsl=0;
  PetscInt     VecSize=0,Istart=0,Iend=0,dnz=0,onz=0,dnzB=0,onzB=0,count,*onnz,*dnnz,*onnzB,*dnnzB,its,ndbbdrank=0;
  PetscErrorCode ierr;
  PetscBool   flg,presol;
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
  printf("rank %d name len %d proc name %s\n",rank,name_len,processor_name);
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
  PetscBool sbbd_overuser=false,nohsl=false;
  char ch='y';
  if(rank==0) {
    printf("***********\nNotes:\n***********\nShock statement values follow GEMPACK ordering (first subscript varies fastest).\nIntertemporal variables (equations) should be declared with minimum dimension to minimise\nthe net cut, e.g. capital(REG,TIME)=qo(\"capital\",REG,TIME) rather than shocking qo(COM,REG,TIME).\nlaA/laDi control solver workspace and temp sizes; use the smallest that solves.\nBeware CRLF line endings in model text files.\n");
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  //**************************************************************************************
  //****************************** READ SET ELEMENT***************************************
  //**************************************************************************************
  char tabfile[TABREADLINE],newtabfile[TABREADLINE]="_temp_tab_file",newtabfile1[TABREADLINE]="_temp_tab_new_file",closure[TABREADLINE],shock[TABREADLINE],filename[TABREADLINE],longname[TABREADLINE],vname[NAMESIZE],copyline[TABREADLINE],regset[NAMESIZE];
  char tempfilenam[255],tempchar[255],solmed[NAMESIZE],solchar[255];
  int niodata=0,nj,mem_fac=0,noutdata=0,nsoldata=0,nowrites=0;
  offset_t nsetspace=0,dcount,ndblock=0,netcut=0,ndblock1,nreg=0,ntime=0;
  dim_t nset=0,vsize,dim1,nlength=0,matsol=0,laA=2,laDi=2,laD=2,nsbbdblocks=2,nesteddbbd=0,mc66=0,subints=1,subindx,StoIter=1;
  offset_t alltimeset=-1,allregset=-1;
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
  isLinux=0;
  PetscOptionsGetInt(NULL,NULL,"-laA",&laA,NULL);
  if(laA==0)laA=2;
  PetscOptionsGetInt(NULL,NULL,"-laD",&laD,NULL);
  if(laD==0)laD=2;
  PetscOptionsGetInt(NULL,NULL,"-laDi",&laDi,NULL);
  if(laDi==0)laDi=2;
  PetscOptionsGetInt(NULL,NULL,"-withmc66",&mc66,NULL);
  PetscOptionsGetInt(NULL,NULL,"-step1",&steps1,NULL);
  if(steps1==0)steps1=2;
  PetscOptionsGetInt(NULL,NULL,"-step2",&steps2,NULL);
  if(steps2==0)steps2=4;
  PetscOptionsGetInt(NULL,NULL,"-step3",&steps3,NULL);
  if(steps3==0)steps3=8;
  step_ratio2=steps1/(double)2;
  i=(offset_t)steps1/2;
  if(step_ratio2!=i){
    //odd
    step_ratio2=steps2/(double)2;
    i=(offset_t)steps2/2;
    if(step_ratio2==i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
    step_ratio2=steps3/(double)2;
    i=(offset_t)steps3/2;
    if(step_ratio2==i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
  }else{
    //even
    step_ratio2=steps2/(double)2;
    i=(offset_t)steps2/2;
    if(step_ratio2!=i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
    step_ratio2=steps3/(double)2;
    i=(offset_t)steps3/2;
    if(step_ratio2!=i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
  }
  
  step_ratio2=steps2/(double)steps1;
  steps2=(PetscInt)steps2/steps1;
  step_ratio3=steps3/(double)steps1;
  steps3=(PetscInt)steps3/steps1;
  section_threads=0;
  max_threads=1;
  PetscOptionsGetInt(NULL,NULL,"-maxthreads",&max_threads,NULL);
  if(max_threads>1&&max_threads<=omp_get_max_threads( )){
    omp_set_num_threads(max_threads);
  }else{
    printf("Max Threads Num = %d\nPlease set OMP_NUM_THREADS <=maxthreads!!! I am setting it to 1!\n",max_threads);
    max_threads=1;
    omp_set_num_threads(max_threads);
  }
  PetscOptionsGetInt(NULL,NULL,"-smllthreads",&section_threads,NULL);
  if(section_threads==0)section_threads=max_threads;
  PetscOptionsGetInt(NULL,NULL,"-nsubints",&subints,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nsbbdblocks",&nsbbdblocks,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nesteddbbd",&nesteddbbd,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nowrites",&nowrites,NULL);
  PetscOptionsGetBool(NULL,NULL,"-presol",&presol,NULL); /* preparation for next solution */
  PetscOptionsGetReal(NULL,NULL,"-cntl_6",&cntl6,NULL); /* CNTL6 in Mat Order */
  PetscOptionsGetReal(NULL,NULL,"-cntl_3",&cntl3,NULL);/*Iterative threshold */
  PetscOptionsGetInt(NULL,NULL,"-ndbbd_bl_rank",&ndbbdrank,NULL);/*Override default rank for last block in NDBBD method. Read text file. >0 text column. Use with care*/
  PetscOptionsGetInt(NULL,NULL,"-stoiter",&StoIter,NULL);
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
        if(rank==0)printf("Error!!! tempdir must be 1-200 characters: %s\n",tmpopt);
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
      if(rank==0)printf("inmemory: scratch on tmpfs (%s)\n",scratch_dir);
    }
    /* export for the Fortran kernels, which build their factor-file
       paths themselves (hsl_kernels.f90) */
    setenv("TEEMS_SCRATCH",scratch_dir,1);
    if(access(scratch_dir,W_OK)!=0) {
      if(rank==0)printf("Error!!! tempdir is not a writable directory: %s\n",scratch_dir);
      PetscFinalize();
      return 1;
    }
  }
  PetscInt *ndbbddrank1= (PetscInt *) calloc (ndbbdrank,sizeof(PetscInt));
  PetscOptionsGetString(NULL,NULL,"-nestfile",filename,TABREADLINE,&flg);
  if (!flg) {
    strcpy(filename,"./ndbbd_drank.csv");//orani03.cmf");
  }
  if (flg)csv_read_ints(filename,ndbbddrank1,ndbbdrank);//
  printf("matsol %d\n",matsol);
  PetscOptionsGetString(NULL,NULL,"-cmdfile",filename,TABREADLINE,&flg);
  if (!flg) {
    strcpy(filename,"./reg.cmf");//orani03.cmf");
  }
  PetscOptionsGetBool(PETSC_NULLPTR,NULL, "-enable_time", &sbbd_overuser,PETSC_NULLPTR);/* Overrid MC66 ordering */
  regset[0]='\0';
  PetscOptionsGetString(NULL,NULL,"-regset",regset,NAMESIZE,&flg);
  if(regset[0]!='\0')for(i=0; i<NAMESIZE; i++) {
      regset[i]=tolower((int)regset[i]);
    }
  PetscOptionsGetString(NULL,NULL,"-solmed",solmed,NAMESIZE,&flg);
  if (!flg) {
    strcpy(solmed,"Mmid");//orani03.cmf");
  }
  int solmethod;
  if(strcmp(solmed,"Mmid")==0)solmethod=SM_MODIFIED_MIDPOINT;
  if(strcmp(solmed,"Johansen")==0)solmethod=SM_JOHANSEN;
  if(strcmp(solmed,"Stochastic")==0)solmethod=SM_STOCHASTIC;
  if(strcmp(solmed,"StoSim")==0)solmethod=SM_STOSIM;
  if(strcmp(solmed,"NoSol")==0)solmethod=SM_NOSOLVE;
  printf("Sol med %d regset %s\n",solmethod,regset);

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
    printf("OK!!!\n");
    niodata=cmf_count_files(filename,"iodata");
    if(niodata==-1)return 0;
    printf("OK!!!\n");
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
    for (nj=0; nj<niodata+noutdata+nsoldata; nj++) printf("rank %d logname %s fname %s\n",rank,iodata[nj].logname,iodata[nj].filname);
    if(tab_preprocess(tabfile,newtabfile)==-1)return 0;
    printf("OK1!\n");
  }
  printf("rank %d OKK1\n",rank);

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
    sets_read(tabfile,niodata,iodata, sets,nset);
    sets_read_intertemporal(tabfile,niodata,iodata, sets,nset);
    for (i=0; i<nset; i++) {
      sets[i].offset=nsetspace;
      nsetspace=nsetspace+sets[i].size;
    }
    if(regset[0]!='\0')for (i=0; i<nset; i++)if(strcmp(regset,sets[i].setname)==0) {
          allregset=i;
          break;
        }
  }
  if(nohsl) {
    MPI_Bcast(sets,nset*sizeof(set_def), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&allregset,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nsetspace,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  printf("rank %d regset %s indx %ld\n",rank,regset,allregset);
  set_element *set_elems= (set_element *) calloc (nsetspace,sizeof(set_element));
  printf("nset %d nsetspace %ld\n",nset,nsetspace);
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
        if (sets[i].readele[0]=='-'&&sets[i].readele[1]==',') {
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
                dim1=sets[i].size;
                strcpy(copyline,sets[i].readele);
                while (str_replace_all(copyline," ", ""));
                readitem = strtok(copyline,"=");
                readitem = strtok(readitem,"=");
                j1=atoi(readitem);
                for (j=0; j<dim1; j++) {
                  strcpy(set_elems[j+sets[i].offset].setele,set_elems[j+sets[j1].offset].setele);
                  set_elems[j+sets[i].offset].superset_pos[0]=j;
                }
              }
              else {
                dim1=sets[i].size;
                strcpy(copyline,sets[i].readele);
                strcat(copyline,",");
                while (str_replace_all(copyline," ", ""));
                readitem = strtok(copyline,",");
                strcpy(set_elems[sets[i].offset].setele,readitem);
                set_elems[sets[i].offset].superset_pos[0]=0;
                for (j=1; j<dim1; j++) {
                  readitem = strtok(NULL,",");
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
    if(sbbd_overuser) {
      alltimeset=set_find_alltime(sets,nset);
    }
    if(alltimeset>=0||allregset>=0) {
      if(alltimeset>=0&&allregset>=0) {
        ndblock=sets[alltimeset].size*sets[allregset].size;
        nreg=sets[allregset].size;
        ntime=sets[alltimeset].size;
      }
      if(alltimeset>=0&&allregset>=0&&nesteddbbd==1) {
        ndblock=sets[alltimeset].size*(sets[allregset].size+1);
        nreg=sets[allregset].size;
        ntime=sets[alltimeset].size;
      }
      if(alltimeset<0&&allregset>=0) {
        ndblock=sets[allregset].size;
        nreg=ndblock;
      }
      if(alltimeset>=0&&allregset<0) {
        ndblock=sets[alltimeset].size;
        ntime=ndblock;
      }
    }
    if(allregset>=0) {
      sets[allregset].regional=true;
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
    if(alltimeset>=0) {
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
    ndblock1=ndblock;
    printf("rank %d alltime set %ld allreg %ld size %ld\n",rank,alltimeset,allregset,ndblock);
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  if(nohsl) {
    MPI_Bcast(&alltimeset,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&allregset,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(sets,nset*sizeof(set_def), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(set_elems,nsetspace*sizeof(set_element), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&ndblock1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&ntime,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nreg,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  if(nesteddbbd==1&&ntime!=ndbbdrank) {
    ndbbdrank=ntime;
    free(ndbbddrank1);
    ndbbddrank1=(PetscInt *) calloc(ndbbdrank,sizeof(PetscInt));
  }
  ndblock=ndblock1;
  printf("rank %d ndblock %ld allreg %ld\n",rank,ndblock,allregset);

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
    printf("tabfile %s\n",tabfile);
    ncof=tab_count_statements(tabfile,commsyntax);
    printf("tabfile %s ncof %ld\n",tabfile,ncof);
    ncof1=ncof;
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  printf("rank %d ncof %ld\n",rank,ncof);
  if(nohsl)MPI_Bcast(&ncof1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
  ncof=ncof1;
  array_def *coefs= (array_def *) calloc (ncof,sizeof(array_def));//recycle ha_cgeset
  if(rank==0) {
    ncofele=coefficients_read(tabfile,commsyntax,coefs,ncof,sets,nset);
    if(ncofele==-1)return 0;
    ncofele1=ncofele;
  }
  if(nohsl) {
    MPI_Bcast(&ncofele1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(coefs,ncof*sizeof(array_def), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  ncofele=ncofele1;
  printf("rank %d ncofele %ld\n",rank,ncofele);
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
  printf("nvarele %ld\n",nvarele);
  if(rank==0) {
    variables_read_defaults(tabfile,vars,nvar);
    nvarele=variables_read(tabfile,commsyntax,vars,nvar,sets,nset);
    nvarele1=nvarele;
  }
  printf("nvarele %ld\n",nvarele);
  if(nohsl) {
    MPI_Bcast(&nvarele1,sizeof(offset_t), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(vars,nvar*sizeof(array_def), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  nvarele=nvarele1;
  if(rank==0){
  for (i=0;i<nvar;i++){
    if(vars[i].cofname[1]=='_'){
    for (j=0;j<ncof;j++)if (strcmp(coefs[j].cofname,vars[i].cofname+2)==0){
      printf("Error!!! Same variable and coefficient names are not supported in this version, even with p_ or c_\nPlease change the name of coefficient %s or variable %s\n",coefs[j].cofname,vars[i].cofname);
      return 0;
    }
    }
  }
  }
  printf("nvarele %ld\n",nvarele);
  elem_value *elem_vals= (elem_value *) calloc ((ncofele+nvarele),sizeof(elem_value));
  elem_store *coef_store= (elem_store *) calloc (ncofele,sizeof(elem_store));
  elem_store *var_store= (elem_store *) calloc (nvarele,sizeof(elem_store));
  printf("rankasd %d nvar %ld\n",rank,nvar);
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
  printf("rank %d ncofele %ld\n",rank,ncofele);

  if(rank==0)coef_resolve_sets(vars,nvar,sets,nset,var_store);
  printf("rank %d OK!!!\n",rank);

  if(rank==rank_hsl)tab_write_variables(tabfile,newtabfile1,vars,nvar);
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
    if(data_read_files(tabfile,niodata,iodata,commsyntax,sets,nset,set_elems,coefs,ncof,coef_store,ncofele,vars,nvar,var_store,nvarele)==-1)return 0;
  }
  printf("rank %d OK???\n",rank);
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
  printf("rank %d ncofvar %ld\n",rank,ncofele+nvarele);
  free(coef_store);
  free(var_store);
  closure_entry *closure_vals= (closure_entry *) calloc (nvarele,sizeof(closure_entry));
  if(rank==0) {
    strcpy(commsyntax,"exogenous");
    nexo=closure_read(closure,commsyntax,closure_vals,vars,nvar,sets,nset,set_elems);
    nexo1=nexo;
    strcpy(commsyntax,"shock");
    if(shocks_read(shock,commsyntax,closure_vals,nvarele,vars,nvar,sets,nset,set_elems,subints)==-1)return 0;
  }
  printf("rank %d OK???\n",rank);
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
  }
  nexo=nexo1;
  strcpy(commsyntax,"formula");
  bool IsIni=true;
  printf("OK???\n");
  if(rank==0) {
    formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
  }
  gettimeofday(&endtime, NULL);
  if(rank==0)printf("Calculation of variables time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
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
  if(rank==0)printf("Broadcast of variables time %f\n",(begintime.tv_sec - endtime.tv_sec)+((double)(begintime.tv_usec - endtime.tv_usec))/ 1000000);
  //**************************************************************************************
  //****************************** MATRIX FROM FORMULA************************************
  //**************************************************************************************
  VecSize = (PetscInt) nvarele-nexo;
  PetscPrintf(PETSC_COMM_SELF,"VecSize %d exo %ld\n",VecSize,nexo);
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
    printf("neq %ld\n",neq);
  }
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
  printf("rank %d\n",rank);
  if(rank==rank_hsl) {
    if(nesteddbbd==1)equation_order_read_nested(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,var_inter,eq_defs,eq_intertemp,eq_time,eq_reg,allregset,alltimeset,orderintra,orderreg);
    else equation_order_read(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,var_inter,eq_defs,eq_intertemp,eq_time,eq_reg,allregset,alltimeset,orderintra,orderreg);
    if(alltimeset>=0||allregset>=0)for(i=0; i<neq; i++)eq_intertemp[i]=!eq_intertemp[i];
  }
  switch (nesteddbbd) {
  case 1 :
    if(!(alltimeset>=0&&allregset>=0))printf("Not a intertemporal regional CGE model!\n");
    offset_t *countvarintra= (offset_t *) calloc (ndblock,sizeof(offset_t));
    j3=0;
    for (i=0; i<nvar; i++) {
      for (j=0; j<vars[i].nelem; j++) {
        if(!closure_vals[j3+j].is_exogenous) {
          if(!var_inter[i]) {
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
            if(orderreg[i]>-1)countvarintra[j2*(nreg+1)+j4]++;
            else countvarintra[j2*(nreg+1)+nreg]++;
          }
        }
      }
      j3+=vars[i].nelem;
    }
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
        if(!closure_vals[j5].is_exogenous) {
          if(!var_inter[i]) {
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
      j3=0;
      for (i=0; i<nvar; i++) {
        for (j=0; j<vars[i].nelem; j++) {
          if(!closure_vals[j3+j].is_exogenous) {
            if(!var_inter[i]) {
              j0=j;
              for(j1=0; j1<orderintra[i]+1; j1++) {
                j2=j0/vars[i].strides[j1];
                j0-=j2*vars[i].strides[j1];
              }
              if(allregset>=0) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j4=j0/vars[i].strides[j1];
                  j0-=j4*vars[i].strides[j1];
                }
                if(sets[vars[i].setid[orderreg[i]]].regsup>0)j4=set_elems[sets[vars[i].setid[orderreg[i]]].offset+j4].superset_pos[sets[vars[i].setid[orderreg[i]]].regsup];
                countvarintra[j2*nreg+j4]++;
              }
              else {
                countvarintra[j2]++;
              }
            }
          }
        }
        j3+=vars[i].nelem;
      }
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
          if(!closure_vals[j5].is_exogenous) {
            if(!var_inter[i]) {
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
        j3=0;
        for (i=0; i<nvar; i++) {
          for (j=0; j<vars[i].nelem; j++) {
            if(!closure_vals[j3+j].is_exogenous) {
              if(!var_inter[i]) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j4=j0/vars[i].strides[j1];
                  j0-=j4*vars[i].strides[j1];
                }
                if(sets[vars[i].setid[orderreg[i]]].regsup>0)j4=set_elems[sets[vars[i].setid[orderreg[i]]].offset+j4].superset_pos[sets[vars[i].setid[orderreg[i]]].regsup];
                countvarintra[j4]++;
              }
            }
          }
          j3+=vars[i].nelem;
        }
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
            if(!closure_vals[j5].is_exogenous) {
              if(!var_inter[i]) {
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

            if (!closure_vals[j3].is_exogenous) {
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
      eq_defs[i].strides[eq_defs[i].size-1]=1;
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
              for(j3=1; j3<MAXSUPSET; j3++)if(sets[eq_defs[i].setid[eq_time[i]]].subsetid[j3]=alltimeset)break;
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
      eq_defs[i].strides[eq_defs[i].size-1]=1;
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
  printf("rank1 %d\n",rank);

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
              for(j3=1; j3<MAXSUPSET; j3++)if(sets[eq_defs[i].setid[eq_time[i]]].subsetid[j3]=alltimeset)break;
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
    printf("netcut %ld nintraeq %ld\n",netcut,nintraeq);
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
    printf("rank %d localbeg %d localend %d\n",rank,localbeg,localend);
    localsize=0;
    for (i=1; i<ndblock+1; i++)if(i>localbeg&&i<=localend)localsize+=counteq[i]-counteq[i-1];
    printf("rank %d localsize %d\n",rank,localsize);
    VecSetSizes(vece,localsize,VecSize);
  }
  else {
    VecSetSizes(vece,PETSC_DECIDE,VecSize);
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
  printf("rank11 %d Istart %d I end %d\n",rank, Istart,Iend);
  if(rank==rank_hsl) {
    jacobian_preallocate(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,nexo,closure_vals,ndblock,alltimeset,allregset,eq_intertemp,eq_addr,eq_time,eq_reg,counteq,nintraeq,&sbbd_overrid,Istart,Iend,&dnz,dnnz,&onz,onnz,&dnzB,dnnzB,&onzB,onnzB,nesteddbbd);
  }
  printf("OKla!!!\n");
  if(sbbd_overrid&&!sbbd_overuser) {
    printf(" It looks like you have an intertemporal model, \n please do -enable_time for more accurate results if \n you know what you are doing!\n");
  }
  if(sbbd_overuser) {
    sbbd_overrid=false;
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
  printf("OK1!!! rank %d\n",rank);
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
  solve_real *x0=NULL;// (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  solve_real *b1=NULL;//= (ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));
  extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, solve_real *VA, solve_real *B, solve_real *X);
  extern void spec48_single_(fortran_int *indata,int *irn, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
  extern void spec48_nomc66_(fortran_int *indata, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow, MPI_Fint *fcomm,fortran_int *rowptrin, fortran_int *colptrin);

  fortran_int k=0,m=1;
  solve_real temp1,temp2;
  fortran_int tindx1;
  printf("rank %d ncof %ld\n",rank,ncof);
  

  if(inmemory) {
    /* Residency cost of skipping the driver spills: value arrays plus the
       modified-midpoint step state. Fall back to scratch files unless it
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
      if(rank==0)printf("Warning!!! -inmemory needs ~%ld MB per rank but only ~%ld MB available; using scratch files instead\n",need/1048576,avail/1048576);
      inmemory=0;
    } else if(rank==0)printf("inmemory: keeping ~%ld MB of value arrays resident per rank\n",need/1048576);
  }
  if(solmethod==SM_JOHANSEN)solve_johansen(nohsl,VecSize,A,dnz,dnnz,onz,onnz,B,dnzB,dnnzB,onzB,onnzB,vecb,vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele+nvarele,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,&xcf);
  
  int stepcount;
  int nsteps=3;
  solve_real vpercents=1.0,perprecis=0;
  FILE* solution;
  int maxsol=3;

    if(solmethod==SM_MODIFIED_MIDPOINT)solve_modified_midpoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele+nvarele,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,0);
    
  if(solmethod==SM_STOCHASTIC){
    MPI_Barrier(PETSC_COMM_WORLD);
  for(j=0;j<StoIter;j++){
    if(j>0){
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
              if(nesteddbbd==1)VecSetSizes(vece,localsize,VecSize);
              else VecSetSizes(vece,PETSC_DECIDE,VecSize);
              VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    }
      if(rank==rank_hsl)subinterval_update(rank,tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,nvarele,10*laA,subints,1,1,0);
      
      solve_modified_midpoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele+nvarele,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,0);
  printf("ncof %ld rank %d\n",ncof,rank);
  }
  }
  if(solmethod==SM_STOSIM){
    MPI_Barrier(PETSC_COMM_WORLD);
    printf("Heere rank %d\n",rank);
      solve_modified_midpoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,&elem_vals,ncofele+nvarele,ncofele,nvarele,&closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,2);
  }
    
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
    printf("solchar %s\n",solchar);
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(xcf, sizeof(solve_real),nvarele, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".var");
    printf("solchar %s\n",solchar);
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(vars, sizeof(array_def),nvar, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".set");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(sets, sizeof(set_def),nset, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".sel");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
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
      printf("Error opening file\n");
      return 1;
    }
    fwrite(modeldes, sizeof(offset_t),4, solution);
    fclose(solution);
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  if(nowrites==0&&rank==0)for(i=0; i<noutdata; i++){
    outputs_write_csv(tabfile,iodata[i+niodata].logname,iodata[i+niodata].filname,sets,nset,set_elems,coefs,ncof,ncofele,vars,nvar,nvarele,elem_vals);
    printf("outfile %s\n",iodata[i+niodata].logname);
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
  printf("Hello world2! sof int %ld sof PetscInt %ld rank %d rankhsl %d\n",sizeof(int),sizeof(PetscInt),rank,rank_hsl);
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


