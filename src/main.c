// Pham Van Ha & Tom Kompas
// Ensure that sum over time involves only intertemporal variables
// If otherwise switch back to normal MC66
// Currently no sub-reg set, can try different varreg to get opt: comment break in nestedmatvarread
// Currently logic fdim <> only in formula
#include <ha_cgeglobal.h>

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
  uvadd lasize=1;
  uvadd i,j;
  uvadd j2=0,j1=0,j0=0,j3,j4,j5,j6;
  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize);

  char processor_name[MPI_MAX_PROCESSOR_NAME+1];
  int name_len,name_len_max,name_beg,class_size,color,group_size,ha_id;
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
  MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&HA_COMM);
  MPI_Comm_rank( HA_COMM, &ha_id);
  MPI_Comm_size(HA_COMM,&group_size);
  if(ha_id==group_size-1)color=1;
  else color=0;
  MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&HA1_COMM);

  gettimeofday(&begintime, NULL);
  bool sbbd_overrid=false;
  PetscBool sbbd_overuser=false,nohsl=false;
  char ch='y';
  if(rank==0) {
    printf("***********\nWarnings:\n***********\nA. This software was written by Pham Van Ha and Tom Kompas for their own research. Use it at your own risk!\nB. The shocks to only 2 dimensions variable are different from GEMPACK. In fact the matrix needs to be transposed.\nIntertemporal variables (equation) must be declared with minimum dimension to minimise\nnet cut (e.g. qo(COM,REG,TIME) must be capital(REG,TIME)=qo(\"capital\",REG,TIME)\nFor NDBBD, LA in Di is not impotant! The procedure automatically determine rank!\nTake care with text file created in windows, dos environment when using gedit. Gedit will use cr for new line! Bug!\nLaA, LaDi are very important. They control the size of temporary files. Enter smallest posible!!!\nDo you agree with the term (y/n)?\n");
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  //**************************************************************************************
  //****************************** READ SET ELEMENT***************************************
  //**************************************************************************************
  char tabfile[TABREADLINE],newtabfile[TABREADLINE]="_temp_tab_file",newtabfile1[TABREADLINE]="_temp_tab_new_file",closure[TABREADLINE],shock[TABREADLINE],filename[TABREADLINE],longname[TABREADLINE],vname[NAMESIZE],copyline[TABREADLINE],regset[NAMESIZE];
  char tempfilenam[255],tempchar[255],solmed[NAMESIZE],solchar[255];
  int niodata=0,nj,mem_fac=0,noutdata=0,nsoldata=0,nowrites=0;
  uvadd nsetspace=0,dcount,ndblock=0,netcut=0,ndblock1,nreg=0,ntime=0;
  uvdim nset=0,vsize,dim1,nlength=0,matsol=0,laA=2,laDi=2,laD=2,nsbbdblocks=2,nesteddbbd=0,mc66=0,subints=1,subindx,StoIter=1;
  uvadd alltimeset=-1,allregset=-1;
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
  PetscOptionsGetInt(NULL,NULL,"-matsol",&matsol,NULL);//0 MA48 (nproc must be 1) 1 SBBD 2 DBBD, if >=1 should enable reg or time. First reg set will be regorder in var. Note INCL(4) in hsl_mp48ss
  if(matsol==2)nohsl=true;
  if(matsol==3)nohsl=true;
  isLinux=0;
  PetscOptionsGetInt(NULL,NULL,"-laA",&laA,NULL);
  if(laA==0)laA=2;
  PetscOptionsGetInt(NULL,NULL,"-laD",&laD,NULL);
  if(laD==0)laD=2;
  PetscOptionsGetInt(NULL,NULL,"-laDi",&laDi,NULL);
  if(laDi==0)laDi=2;
  PetscOptionsGetInt(NULL,NULL,"-withmc66",&mc66,NULL);
  PetscOptionsGetInt(NULL,NULL,"-step1",&step1,NULL);
  if(step1==0)step1=2;
  PetscOptionsGetInt(NULL,NULL,"-step2",&step2,NULL);
  if(step2==0)step2=4;
  PetscOptionsGetInt(NULL,NULL,"-step3",&step3,NULL);
  if(step3==0)step3=8;
  kindx1=step1/(double)2;
  i=(uvadd)step1/2;
  if(kindx1!=i){
    //odd
    kindx1=step2/(double)2;
    i=(uvadd)step2/2;
    if(kindx1==i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
    kindx1=step3/(double)2;
    i=(uvadd)step3/2;
    if(kindx1==i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
  }else{
    //even
    kindx1=step2/(double)2;
    i=(uvadd)step2/2;
    if(kindx1!=i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
    kindx1=step3/(double)2;
    i=(uvadd)step3/2;
    if(kindx1!=i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
  }
  
  kindx1=step2/(double)step1;
  step2=(PetscInt)step2/step1;
  kindx2=step3/(double)step1;
  step3=(PetscInt)step3/step1;
  smallthreads=0;
  mymaxnumthrd=1;
  PetscOptionsGetInt(NULL,NULL,"-maxthreads",&mymaxnumthrd,NULL);
  if(mymaxnumthrd>1&&mymaxnumthrd<=omp_get_max_threads( )){
    omp_set_num_threads(mymaxnumthrd);
  }else{
    printf("Max Threads Num = %d\nPlease set OMP_NUM_THREADS <=maxthreads!!! I am setting it to 1!\n",mymaxnumthrd);
    mymaxnumthrd=1;
    omp_set_num_threads(mymaxnumthrd);
  }
  PetscOptionsGetInt(NULL,NULL,"-smllthreads",&smallthreads,NULL);
  if(smallthreads==0)smallthreads=mymaxnumthrd;
  PetscOptionsGetInt(NULL,NULL,"-nsubints",&subints,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nsbbdblocks",&nsbbdblocks,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nesteddbbd",&nesteddbbd,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nowrites",&nowrites,NULL);
  PetscOptionsGetBool(NULL,NULL,"-presol",&presol,NULL); /* preparation for next solution */
  PetscOptionsGetReal(NULL,NULL,"-cntl_6",&cntl6,NULL); /* CNTL6 in Mat Order */
  PetscOptionsGetReal(NULL,NULL,"-cntl_3",&cntl3,NULL);/*Iterative threshold */
  PetscOptionsGetInt(NULL,NULL,"-ndbbd_bl_rank",&ndbbdrank,NULL);/*Override default rank for last block in NDBBD method. Read text file. >0 text column. Use with care*/
  PetscOptionsGetInt(NULL,NULL,"-stoiter",&StoIter,NULL);
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
      strcpy(temdir,tmpopt);
    }
    if(access(temdir,W_OK)!=0) {
      if(rank==0)printf("Error!!! tempdir is not a writable directory: %s\n",temdir);
      PetscFinalize();
      return 1;
    }
  }
  PetscInt *ndbbddrank1= (PetscInt *) calloc (ndbbdrank,sizeof(PetscInt));
  PetscOptionsGetString(NULL,NULL,"-nestfile",filename,TABREADLINE,&flg);
  if (!flg) {
    strcpy(filename,"./ndbbd_drank.csv");//orani03.cmf");
  }
  if (flg)intreadCSV(filename,ndbbddrank1,ndbbdrank);//
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
  if(strcmp(solmed,"Mmid")==0)solmethod=1;
  if(strcmp(solmed,"Johansen")==0)solmethod=10;
  if(strcmp(solmed,"Stochastic")==0)solmethod=20;
  if(strcmp(solmed,"StoSim")==0)solmethod=21;
  if(strcmp(solmed,"NoSol")==0)solmethod=100;
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
    niodata=hcge_niodata(filename,"iodata");
    if(niodata==-1)return 0;
    printf("OK!!!\n");
    noutdata=hcge_niodata(filename,"outdata");
    nsoldata=hcge_niodata(filename,"soldata");
  }
  if(nohsl) {
    MPI_Bcast(&niodata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&noutdata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nsoldata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  hcge_iodata *iodata= (hcge_iodata *) calloc (niodata+noutdata+nsoldata,sizeof(hcge_iodata));
  if(rank==rank_hsl) {
    hcge_rcmd(filename,niodata,iodata,tabfile,closure,shock);
    for (nj=0; nj<niodata+noutdata+nsoldata; nj++) printf("rank %d logname %s fname %s\n",rank,iodata[nj].logname,iodata[nj].filname);
    if(hcge_wtab(tabfile,newtabfile)==-1)return 0;
    printf("OK1!\n");
  }
  printf("rank %d OKK1\n",rank);

  strcpy(tabfile,newtabfile);
  if(rank==0)nset=ha_cgenset(tabfile);
  if(nohsl) {
    MPI_Bcast(iodata,niodata*sizeof(hcge_iodata), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(closure,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(shock,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nset,sizeof(uvdim), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  ha_cgeset *ha_set= (ha_cgeset *) calloc (nset,sizeof(ha_cgeset));
  for(i=0; i<nset; i++) {
    ha_set[i].subsetid[0]=i;
    for(j=0; j<MAXSUPSET; j++)ha_set[i].subsetid[j]=-1;
  }
  if(rank==0) {
    ha_cgerset(tabfile,niodata,iodata, ha_set,nset);
    hcge_rinterset(tabfile,niodata,iodata, ha_set,nset);
    for (i=0; i<nset; i++) {
      ha_set[i].begadd=nsetspace;
      nsetspace=nsetspace+ha_set[i].size;
    }
    if(regset[0]!='\0')for (i=0; i<nset; i++)if(strcmp(regset,ha_set[i].setname)==0) {
          allregset=i;
          break;
        }
  }
  if(nohsl) {
    MPI_Bcast(ha_set,nset*sizeof(ha_cgeset), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&allregset,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nsetspace,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  printf("rank %d regset %s indx %ld\n",rank,regset,allregset);
  ha_cgesetele *ha_setele= (ha_cgesetele *) calloc (nsetspace,sizeof(ha_cgesetele));
  printf("nset %d nsetspace %ld\n",nset,nsetspace);
  for (i=0; i<nsetspace; i++)for (j=0; j<MAXSUPSET; j++)ha_setele[i].setsh[j]=-1;
  if(rank==0) {
    for (i=0; i<nset; i++) {
      strcpy(vname,ha_set[i].header);
      nlength=0;
      while (vname[nlength] != '\0') {
        nlength++;
      }
      if (nlength>0) {
        ha_cgerdvar1(vname,iodata[ha_set[i].fileid].filname,&vsize,longname,&dim1);
        ha_cgemvar1 *matvar1= (ha_cgemvar1 *) calloc (dim1,sizeof(ha_cgemvar1));
        ha_cgermvar1(vname,iodata[ha_set[i].fileid].filname,dim1,matvar1);
        for (j=0; j<dim1; j++) {
          nj=0;
          while(matvar1[j].ch[nj]!='\0') {
            matvar1[j].ch[nj]=tolower((int) matvar1[j].ch[nj]);
            nj++;
          }
          strncpy(ha_setele[j+ha_set[i].begadd].setele,matvar1[j].ch,strlen(matvar1[j].ch));
          ha_setele[j+ha_set[i].begadd].setsh[0]=j;
        }
        free(matvar1);
      }
      else {
        if (ha_set[i].readele[0]=='-'&&ha_set[i].readele[1]==',') {
          ha_setminus(ha_setele, ha_set,nset,i);
        }
        else {
          if (ha_set[i].readele[0]=='+'&&ha_set[i].readele[1]==',') {
            ha_setplus(ha_setele, ha_set,nset,i);
          }
          else {
            if (ha_set[i].readele[0]=='^'&&ha_set[i].readele[1]==',') {
              ha_setunion(ha_setele, ha_set,nset,i);
            }
            else {
              if(ha_set[i].readele[0]=='=') {
                dim1=ha_set[i].size;
                strcpy(copyline,ha_set[i].readele);
                while (ha_cgefrstr(copyline," ", ""));
                readitem = strtok(copyline,"=");
                readitem = strtok(readitem,"=");
                j1=atoi(readitem);
                for (j=0; j<dim1; j++) {
                  strcpy(ha_setele[j+ha_set[i].begadd].setele,ha_setele[j+ha_set[j1].begadd].setele);
                  ha_setele[j+ha_set[i].begadd].setsh[0]=j;
                }
              }
              else {
                dim1=ha_set[i].size;
                strcpy(copyline,ha_set[i].readele);
                strcat(copyline,",");
                while (ha_cgefrstr(copyline," ", ""));
                readitem = strtok(copyline,",");
                strcpy(ha_setele[ha_set[i].begadd].setele,readitem);
                ha_setele[ha_set[i].begadd].setsh[0]=0;
                for (j=1; j<dim1; j++) {
                  readitem = strtok(NULL,",");
                  strcpy(ha_setele[j+ha_set[i].begadd].setele,readitem);
                  ha_setele[j+ha_set[i].begadd].setsh[0]=j;
                }
              }
            }
          }
        }
      }

    }
    ha_cgersubset(tabfile, ha_setele, ha_set,nset);
    j2=1;
    while(j2==1)for(i=1; i<MAXSUPSET; i++)ha_cgesubsetchck(ha_setele,ha_set,nset,&j2); //printf("check %d\n",i);}
    if(sbbd_overuser) {
      alltimeset=ha_cgeralltime(ha_set,nset);
    }
    if(alltimeset>=0||allregset>=0) {
      if(alltimeset>=0&&allregset>=0) {
        ndblock=ha_set[alltimeset].size*ha_set[allregset].size;
        nreg=ha_set[allregset].size;
        ntime=ha_set[alltimeset].size;
      }
      if(alltimeset>=0&&allregset>=0&&nesteddbbd==1) {
        ndblock=ha_set[alltimeset].size*(ha_set[allregset].size+1);
        nreg=ha_set[allregset].size;
        ntime=ha_set[alltimeset].size;
      }
      if(alltimeset<0&&allregset>=0) {
        ndblock=ha_set[allregset].size;
        nreg=ndblock;
      }
      if(alltimeset>=0&&allregset<0) {
        ndblock=ha_set[alltimeset].size;
        ntime=ndblock;
      }
    }
    if(allregset>=0) {
      ha_set[allregset].regional=true;
      for(i=0; i<nset; i++) {
        for(j=1; j<MAXSUPSET; j++) {
          if(ha_set[i].subsetid[j]>-1) {
            if(ha_set[ha_set[i].subsetid[j]].regional) {
              ha_set[i].regional=true;
              ha_set[i].regsup=j;
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
          if(ha_set[i].subsetid[j]>-1) {
            if(ha_set[ha_set[i].subsetid[j]].intertemp) {
              ha_set[i].intsup=j;
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
    MPI_Bcast(&alltimeset,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&allregset,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(ha_set,nset*sizeof(ha_cgeset), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(ha_setele,nsetspace*sizeof(ha_cgesetele), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&ndblock1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&ntime,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nreg,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
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
  uvadd ncof=0,ncofele=0,ncof1,ncofele1;
  if(rank==0) {
    printf("tabfile %s\n",tabfile);
    ncof=ha_cgencof(tabfile,commsyntax);
    printf("tabfile %s ncof %ld\n",tabfile,ncof);
    ncof1=ncof;
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  printf("rank %d ncof %ld\n",rank,ncof);
  if(nohsl)MPI_Bcast(&ncof1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  ncof=ncof1;
  hcge_cof *ha_cof= (hcge_cof *) calloc (ncof,sizeof(hcge_cof));//recycle ha_cgeset
  if(rank==0) {
    ncofele=hcge_rcof(tabfile,commsyntax,ha_cof,ncof,ha_set,nset);
    if(ncofele==-1)return 0;
    ncofele1=ncofele;
  }
  if(nohsl) {
    MPI_Bcast(&ncofele1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(ha_cof,ncof*sizeof(hcge_cof), MPI_BYTE,0, PETSC_COMM_WORLD);
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
  uvadd nvar=0,nvarele=0,nvar1,nvarele1;
  if(rank==0) {
    nvar=ha_cgencof(tabfile,commsyntax);
    nvar1=nvar;
  }
  if(nohsl)MPI_Bcast(&nvar1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  nvar=nvar1;
  hcge_cof *ha_var= (hcge_cof *) calloc (nvar,sizeof(hcge_cof));//recycle ha_cgeset
  bool *var_inter= (bool *) calloc (nvar,sizeof(bool));//recycle ha_cgeset
  printf("nvarele %ld\n",nvarele);
  if(rank==0) {
    hcge_defvar(tabfile,ha_var,nvar);
    nvarele=hcge_rvar(tabfile,commsyntax,ha_var,nvar,ha_set,nset);
    nvarele1=nvarele;
  }
  printf("nvarele %ld\n",nvarele);
  if(nohsl) {
    MPI_Bcast(&nvarele1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(ha_var,nvar*sizeof(hcge_cof), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  nvarele=nvarele1;
  if(rank==0){
  for (i=0;i<nvar;i++){
    if(ha_var[i].cofname[1]=='_'){
    for (j=0;j<ncof;j++)if (strcmp(ha_cof[j].cofname,ha_var[i].cofname+2)==0){
      printf("Error!!! Same variable and coefficient names are not supported in this version, even with p_ or c_\nPlease change the name of coefficient %s or variable %s\n",ha_cof[j].cofname,ha_var[i].cofname);
      return 0;
    }
    }
  }
  }
  printf("nvarele %ld\n",nvarele);
  ha_cgevar *ha_cofvar= (ha_cgevar *) calloc ((ncofele+nvarele),sizeof(ha_cgevar));
  ha_cgecofele *ha_cofele= (ha_cgecofele *) calloc (ncofele,sizeof(ha_cgecofele));
  ha_cgecofele *ha_varele= (ha_cgecofele *) calloc (nvarele,sizeof(ha_cgecofele));
  printf("rankasd %d nvar %ld\n",rank,nvar);
  if(rank==0) {
    hcge_rcofele(ha_cof,ncof,ha_set,nset,ha_cofele);
  }
  if(nohsl) {
    if(ncofele*sizeof(ha_cgecofele)>1500000000) {
      j1=1500000000/sizeof(ha_cgecofele);
      i=ncofele/j1;
      for(j=0; j<i; j++) {
        MPI_Bcast(ha_cofele+j*j1,j1*sizeof(ha_cgecofele), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
      i=ncofele-j*j1;
      MPI_Bcast(ha_cofele+j*j1,i*sizeof(ha_cgecofele), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(ha_cofele,ncofele*sizeof(ha_cgecofele), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
  }
  printf("rank %d ncofele %ld\n",rank,ncofele);

  if(rank==0)hcge_rcofele(ha_var,nvar,ha_set,nset,ha_varele);
  printf("rank %d OK!!!\n",rank);

  if(rank==rank_hsl)hcge_wvar(tabfile,newtabfile1,ha_var,nvar);
  strcpy(tabfile,newtabfile1);
  if(nohsl) {
    MPI_Bcast(ha_set,nset*sizeof(ha_cgeset), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  //**************************************************************************************
  //****************************** END READ VARIABLE NAME*********************************
  //**************************************************************************************

  //**************************************************************************************
  //********************* READ VARIABLE, COEFFICIENT VALUE FROM FILE**********************
  //**************************************************************************************
  if(rank==0) {
    strcpy(commsyntax,"read");
    if(hcge_readff(tabfile,niodata,iodata,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_cofele,ncofele,ha_var,nvar,ha_varele,nvarele)==-1)return 0;
  }
  printf("rank %d OK???\n",rank);
  //**************************************************************************************
  //********************* END READ VARIABLE, COEFFICIENT VALUE FROM FILE******************
  //**************************************************************************************

  //**************************************************************************************
  //***************** CALCULATE VARIABLE, COEFFICIENT VALUE FROM FORMULA******************
  //**************************************************************************************

  uvadd nexo=0,nexo1;
  ha_cgevar *ha_cofvar1;
  if(rank==0) {
    for (i=0; i<ncofele; i++) {
      ha_cofvar[i].varval=ha_cofele[i].cofval;
    }
    for (i=ncofele; i<nvarele+ncofele; i++) {
      ha_cofvar[i].varval=ha_varele[i-ncofele].cofval;
    }
  }
  printf("rank %d ncofvar %ld\n",rank,ncofele+nvarele);
  free(ha_cofele);
  free(ha_varele);
  ha_cgeexovar *ha_cgeshock= (ha_cgeexovar *) calloc (nvarele,sizeof(ha_cgeexovar));
  if(rank==0) {
    strcpy(commsyntax,"exogenous");
    nexo=hcge_rexo(closure,commsyntax,ha_cgeshock,ha_var,nvar,ha_set,nset,ha_setele);
    nexo1=nexo;
    strcpy(commsyntax,"shock");
    if(hcge_rshock(shock,commsyntax,ha_cgeshock,nvarele,ha_var,nvar,ha_set,nset,ha_setele,subints)==-1)return 0;
  }
  printf("rank %d OK???\n",rank);
  if(nohsl) {
    if(nvarele*sizeof(ha_cgeexovar)>1500000000) {
      j1=1500000000/sizeof(ha_cgeexovar);
      i=nvarele/j1;
      for(j=0; j<i; j++) {
        MPI_Bcast(ha_cgeshock+j*j1,j1*sizeof(ha_cgeexovar), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
      i=nvarele-j*j1;
      MPI_Bcast(ha_cgeshock+j*j1,i*sizeof(ha_cgeexovar), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(ha_cgeshock,nvarele*sizeof(ha_cgeexovar), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    MPI_Bcast(&nexo1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  nexo=nexo1;
  strcpy(commsyntax,"formula");
  bool IsIni=true;
  printf("OK???\n");
  if(rank==0) {
    hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
  }
  gettimeofday(&endtime, NULL);
  if(rank==0)printf("Calculation of variables time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
  if(nohsl) { //Overcome MPI_Bcast limit
    if((nvarele+ncofele)*sizeof(ha_cgevar)>1500000000) {
      j1=1500000000/sizeof(ha_cgevar);
      i=(nvarele+ncofele)/j1;
      for(j=0; j<i; j++) {
        MPI_Bcast(ha_cofvar+j*j1,j1*sizeof(ha_cgevar), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
      i=nvarele+ncofele-j*j1;
      MPI_Bcast(ha_cofvar+j*j1,i*sizeof(ha_cgevar), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(ha_cofvar,(nvarele+ncofele)*sizeof(ha_cgevar), MPI_BYTE,0, PETSC_COMM_WORLD);
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
  uvadd neq=0,neq1;
  if(rank==0) {
    neq=ha_cgencof(tabfile,commsyntax);
    neq1=neq;
  }
  if(nohsl) {
    MPI_Bcast(&neq1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  neq=neq1;
  if(rank==rank_hsl) {
    printf("neq %ld\n",neq);
  }
  uvadd *countvarintra1= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  hcge_cof *ha_eq= (hcge_cof *) calloc (neq,sizeof(hcge_cof));//recycle ha_cgeset
  bool *ha_eqint= (bool *) calloc (neq,sizeof(bool));//recycle ha_cgeset
  uvdim *orderintra= (uvdim *) malloc (nvar*sizeof(uvdim));
  uvdim *orderreg= (uvdim *) malloc (nvar*sizeof(uvdim));
  for(i=0; i<nvar; i++) {
    orderintra[i]=-1;
    orderreg[i]=-1;
  }
  uvdim *ha_eqtime= (uvdim *) malloc (neq*sizeof(uvdim));//recycle ha_cgeset
  uvdim *ha_eqreg= (uvdim *) malloc (neq*sizeof(uvdim));//recycle ha_cgeset
  for(i=0; i<neq; i++) {
    ha_eqtime[i]=-1;
    ha_eqreg[i]=-1;
  }
  uvadd nintraendovar,summat;
  printf("rank %d\n",rank);
  if(rank==rank_hsl) {
    if(nesteddbbd==1)NestedMatvarRead(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,var_inter,ha_eq,ha_eqint,ha_eqtime,ha_eqreg,allregset,alltimeset,orderintra,orderreg);
    else NewMatvarRead(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,var_inter,ha_eq,ha_eqint,ha_eqtime,ha_eqreg,allregset,alltimeset,orderintra,orderreg);
    if(alltimeset>=0||allregset>=0)for(i=0; i<neq; i++)ha_eqint[i]=!ha_eqint[i];
  }
  switch (nesteddbbd) {
  case 1 :
    if(!(alltimeset>=0&&allregset>=0))printf("Not a intertemporal regional CGE model!\n");
    uvadd *countvarintra= (uvadd *) calloc (ndblock,sizeof(uvadd));
    j3=0;
    for (i=0; i<nvar; i++) {
      for (j=0; j<ha_var[i].matsize; j++) {
        if(!ha_cgeshock[j3+j].ShockId) {
          if(!var_inter[i]) {
            j0=j;
            j2=-1;
            for(j1=0; j1<orderintra[i]+1; j1++) {
              j2=j0/ha_var[i].antidims[j1];
              j0-=j2*ha_var[i].antidims[j1];
            }
            j0=j;
            j4=-1;
            for(j1=0; j1<orderreg[i]+1; j1++) {
              j4=j0/ha_var[i].antidims[j1];
              j0-=j4*ha_var[i].antidims[j1];
            }
            if(j4>-1)if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j4=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j4].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
            if(ha_set[ha_var[i].setid[orderintra[i]]].intsup>0)j2=ha_setele[ha_set[ha_var[i].setid[orderintra[i]]].begadd+j2].setsh[ha_set[ha_var[i].setid[orderintra[i]]].intsup];
            if(orderreg[i]>-1)countvarintra[j2*(nreg+1)+j4]++;
            else countvarintra[j2*(nreg+1)+nreg]++;
          }
        }
      }
      j3+=ha_var[i].matsize;
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
      for (j=0; j<ha_var[i].matsize; j++) {
        j5=j3+j;
        if(!ha_cgeshock[j5].ShockId) {
          if(!var_inter[i]) {
            j0=j;
            j2=-1;
            for(j1=0; j1<orderintra[i]+1; j1++) {
              j2=j0/ha_var[i].antidims[j1];
              j0-=j2*ha_var[i].antidims[j1];
            }
            j0=j;
            j6=-1;
            for(j1=0; j1<orderreg[i]+1; j1++) {
              j6=j0/ha_var[i].antidims[j1];
              j0-=j6*ha_var[i].antidims[j1];
            }
            if(ha_set[ha_var[i].setid[orderintra[i]]].intsup>0)j2=ha_setele[ha_set[ha_var[i].setid[orderintra[i]]].begadd+j2].setsh[ha_set[ha_var[i].setid[orderintra[i]]].intsup];
            if(j6>-1)if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j6=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j6].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
            if(orderreg[i]>-1) {
              ha_cgeshock[j5].ExoIndx=countvarintra[j2*(nreg+1)+j6];
              countvarintra[j2*(nreg+1)+j6]++;
            }
            else {
              ha_cgeshock[j5].ExoIndx=countvarintra[j2*(nreg+1)+nreg];
              countvarintra[j2*(nreg+1)+nreg]++;
            }
          }
          else {
            ha_cgeshock[j5].ExoIndx=j4;
            j4++;
          }
        }
      }
      j3+=ha_var[i].matsize;
    }
    if(rank==rank_hsl) {
      countvarintra1[ndblock]=countvarintra[ndblock-1];
    }
    j1=0;
    for (i=0; i<nvarele; i++) {
      if (ha_cgeshock[i].ShockId) {
        ha_cgeshock[i].ExoIndx+=j1;
        j1++;
      }
    }
    free(countvarintra);

    break;
  default :
    if(alltimeset>=0) {
      uvadd *countvarintra= (uvadd *) calloc (ndblock,sizeof(uvadd));
      j3=0;
      for (i=0; i<nvar; i++) {
        for (j=0; j<ha_var[i].matsize; j++) {
          if(!ha_cgeshock[j3+j].ShockId) {
            if(!var_inter[i]) {
              j0=j;
              for(j1=0; j1<orderintra[i]+1; j1++) {
                j2=j0/ha_var[i].antidims[j1];
                j0-=j2*ha_var[i].antidims[j1];
              }
              if(allregset>=0) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j4=j0/ha_var[i].antidims[j1];
                  j0-=j4*ha_var[i].antidims[j1];
                }
                if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j4=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j4].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
                countvarintra[j2*nreg+j4]++;
              }
              else {
                countvarintra[j2]++;
              }
            }
          }
        }
        j3+=ha_var[i].matsize;
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
        for (j=0; j<ha_var[i].matsize; j++) {
          j5=j3+j;
          if(!ha_cgeshock[j5].ShockId) {
            if(!var_inter[i]) {
              j0=j;
              for(j1=0; j1<orderintra[i]+1; j1++) {
                j2=j0/ha_var[i].antidims[j1];
                j0-=j2*ha_var[i].antidims[j1];
              }
              if(allregset>=0) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j6=j0/ha_var[i].antidims[j1];
                  j0-=j6*ha_var[i].antidims[j1];
                }
                if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j6=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j6].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
                ha_cgeshock[j5].ExoIndx=countvarintra[j2*ha_set[allregset].size+j6];
                countvarintra[j2*ha_set[allregset].size+j6]++;
              }
              else {
                ha_cgeshock[j5].ExoIndx=countvarintra[j2];
                countvarintra[j2]++;
              }
            }
            else {
              ha_cgeshock[j5].ExoIndx=j4;
              j4++;
            }
          }
        }
        j3+=ha_var[i].matsize;
      }
      if(rank==rank_hsl) {
        countvarintra1[ndblock]=countvarintra[ndblock-1];
      }
      j1=0;
      for (i=0; i<nvarele; i++) {
        if (ha_cgeshock[i].ShockId) {
          ha_cgeshock[i].ExoIndx+=j1;
          j1++;
        }
      }
      free(countvarintra);
    }
    else {
      if(allregset>=0) {
        uvadd *countvarintra= (uvadd *) calloc (ndblock,sizeof(uvadd));
        j3=0;
        for (i=0; i<nvar; i++) {
          for (j=0; j<ha_var[i].matsize; j++) {
            if(!ha_cgeshock[j3+j].ShockId) {
              if(!var_inter[i]) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j4=j0/ha_var[i].antidims[j1];
                  j0-=j4*ha_var[i].antidims[j1];
                }
                if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j4=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j4].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
                countvarintra[j4]++;
              }
            }
          }
          j3+=ha_var[i].matsize;
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
          for (j=0; j<ha_var[i].matsize; j++) {
            j5=j3+j;
            if(!ha_cgeshock[j5].ShockId) {
              if(!var_inter[i]) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j6=j0/ha_var[i].antidims[j1];
                  j0-=j6*ha_var[i].antidims[j1];
                }
                if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j6=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j6].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
                ha_cgeshock[j5].ExoIndx=countvarintra[j6];
                countvarintra[j6]++;
              }
              else {
                ha_cgeshock[j5].ExoIndx=j4;
                j4++;
              }
            }
          }
          j3+=ha_var[i].matsize;
        }
        if(rank==rank_hsl) {
          countvarintra1[ndblock]=countvarintra[ndblock-1];
        }
        j1=0;
        for (i=0; i<nvarele; i++) {
          if (ha_cgeshock[i].ShockId) {
            ha_cgeshock[i].ExoIndx+=j1;
            j1++;
          }
        }
        free(countvarintra);
      }
      else {
        j1=0;
        j2=0;
        for (i=0; i<nvar; i++) {
          for (j=0; j<ha_var[i].matsize; j++) {
            j3=j0+j;

            if (!ha_cgeshock[j3].ShockId) {
              ha_cgeshock[j3].ExoIndx+=j2;
              j2++;
            }
            if (ha_cgeshock[j3].ShockId) {
              ha_cgeshock[j3].ExoIndx+=j1;
              j1++;
            }
          }
          j0+=ha_var[i].matsize;
        }
      }
    }
  }
  free(orderintra);
  free(orderreg);

  free(var_inter);
  strcpy(commsyntax,"equation");
  uvadd *ha_eqadd= (uvadd *) calloc (VecSize,sizeof(uvadd));//recycle ha_cgeset
  uvadd *ha_eqtimesbegad= (uvadd *) calloc (neq,sizeof(uvadd));
  uvadd *ha_eqregsbegad= (uvadd *) calloc (neq,sizeof(uvadd));
  uvadd *counteq= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *counteqnoadd= (uvadd *) calloc (ndblock,sizeof(uvadd));
  uvadd nintraeq=0;
  if(alltimeset>=0&&allregset<0) {
    for(i=0; i<neq; i++) {
      j3=1;
      if(ha_eqtime[i]>-1)ha_eqtimesbegad[i]=ha_set[ha_eq[i].setid[ha_eqtime[i]]].begadd;
      ha_eq[i].antidims[ha_eq[i].size-1]=1;
      if(ha_eq[i].size>1) {
        for (j2=ha_eq[i].size-2; j2>-1; j2--) {
          ha_eq[i].antidims[j2]=ha_eq[i].antidims[j2+1]*ha_set[ha_eq[i].setid[j2+1]].size;
        }
      }
    }
    j3=0;
    for (i=0; i<neq; i++) {
      if(ha_eqint[i]) {
        for (j=0; j<ha_eq[i].matsize; j++) {
          j0=j;
          for(j1=0; j1<ha_eqtime[i]+1; j1++) {
            j2=j0/ha_eq[i].antidims[j1];
            j0-=j2*ha_eq[i].antidims[j1];
          }
          if(ha_eqtime[i]>-1)if(ha_eq[i].setid[ha_eqtime[i]]==alltimeset)counteq[ha_setele[ha_eqtimesbegad[i]+j2].setsh[0]]++;
            else {
              for(j3=1; j3<MAXSUPSET; j3++)if(ha_set[ha_eq[i].setid[ha_eqtime[i]]].subsetid[j3]=alltimeset)break;
              counteq[ha_setele[ha_eqtimesbegad[i]+j2].setsh[j3]]++;
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
      if(ha_eqreg[i]>-1)ha_eqtimesbegad[i]=ha_set[ha_eq[i].setid[ha_eqreg[i]]].begadd;
      ha_eq[i].antidims[ha_eq[i].size-1]=1;
      if(ha_eq[i].size>1) {
        for (j2=ha_eq[i].size-2; j2>-1; j2--) {
          ha_eq[i].antidims[j2]=ha_eq[i].antidims[j2+1]*ha_set[ha_eq[i].setid[j2+1]].size;
        }
      }
    }
    j3=0;
    for (i=0; i<neq; i++) {
      if(ha_eqint[i]) { //for (j=0; j<ha_set[allregset].size; j++)counteq[j]+=(uvadd)ha_eq[i].matsize/ha_set[allregset].size;
        for (j=0; j<ha_eq[i].matsize; j++) {
          j0=j;
          for(j1=0; j1<ha_eqreg[i]+1; j1++) {
            j2=j0/ha_eq[i].antidims[j1];
            j0-=j2*ha_eq[i].antidims[j1];
          }
          counteq[ha_setele[ha_eqtimesbegad[i]+j2].setsh[0]]++;
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
        if(ha_eqtime[i]>-1)ha_eqtimesbegad[i]=ha_set[ha_eq[i].setid[ha_eqtime[i]]].begadd;
        if(ha_eqreg[i]>-1)ha_eqregsbegad[i]=ha_set[ha_eq[i].setid[ha_eqreg[i]]].begadd;
        if(ha_eqint[i]) {
          for (j=0; j<ha_set[ha_eq[i].setid[ha_eqtime[i]]].size; j++)
            if(ha_eqreg[i]>-1)for(j1=0; j1<ha_set[ha_eq[i].setid[ha_eqreg[i]]].size; j1++)
                counteq[ha_setele[ha_eqtimesbegad[i]+j].setsh[ha_set[ha_eq[i].setid[ha_eqtime[i]]].intsup]*(nreg+1)+ha_setele[ha_eqregsbegad[i]+j1].setsh[ha_set[ha_eq[i].setid[ha_eqreg[i]]].regsup]]+=ha_eq[i].matsize/ha_set[ha_eq[i].setid[ha_eqtime[i]]].size/ha_set[ha_eq[i].setid[ha_eqreg[i]]].size;
            else counteq[ha_setele[ha_eqtimesbegad[i]+j].setsh[ha_set[ha_eq[i].setid[ha_eqtime[i]]].intsup]*(nreg+1)+nreg]+=ha_eq[i].matsize/ha_set[ha_eq[i].setid[ha_eqtime[i]]].size;
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
        if(ha_eqtime[i]>-1)ha_eqtimesbegad[i]=ha_set[ha_eq[i].setid[ha_eqtime[i]]].begadd;
        if(ha_eqreg[i]>-1)ha_eqregsbegad[i]=ha_set[ha_eq[i].setid[ha_eqreg[i]]].begadd;
      }
      for (i=0; i<neq; i++) {
        if(ha_eqint[i]) {
          for (j=0; j<ha_set[ha_eq[i].setid[ha_eqtime[i]]].size; j++)
            if(ha_eq[i].setid[ha_eqtime[i]]==alltimeset)
              for(j1=0; j1<ha_set[ha_eq[i].setid[ha_eqreg[i]]].size; j1++)
                counteq[ha_setele[ha_eqtimesbegad[i]+j].setsh[0]*ha_set[ha_eq[i].setid[ha_eqreg[i]]].size+ha_setele[ha_eqregsbegad[i]+j1].setsh[0]]+=ha_eq[i].matsize/ha_set[ha_eq[i].setid[ha_eqtime[i]]].size/ha_set[ha_eq[i].setid[ha_eqreg[i]]].size;
            else {
              for(j3=1; j3<MAXSUPSET; j3++)if(ha_set[ha_eq[i].setid[ha_eqtime[i]]].subsetid[j3]=alltimeset)break;
              for(j1=0; j1<ha_set[ha_eq[i].setid[ha_eqreg[i]]].size; j1++)
                counteq[ha_setele[ha_eqtimesbegad[i]+j].setsh[j3]*ha_set[ha_eq[i].setid[ha_eqreg[i]]].size+ha_setele[ha_eqregsbegad[i]+j1].setsh[0]]+=ha_eq[i].matsize/ha_set[ha_eq[i].setid[ha_eqtime[i]]].size/ha_set[ha_eq[i].setid[ha_eqreg[i]]].size;
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
  free(ha_eq);
  free(ha_eqtimesbegad);
  free(ha_eqregsbegad);
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
    NewMatreadele(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,nexo,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqint,ha_eqadd,ha_eqtime,ha_eqreg,counteq,nintraeq,&sbbd_overrid,Istart,Iend,&dnz,dnnz,&onz,onnz,&dnzB,dnnzB,&onzB,onnzB,nesteddbbd);
  }
  printf("OKla!!!\n");
  if(sbbd_overrid&&!sbbd_overuser) {
    printf(" It looks like you have an intertemporal model, \n please do -enable_time for more accurate results if \n you know what you are doing!\n");
  }
  if(sbbd_overuser) {
    sbbd_overrid=false;
  }
  free(ha_eqint);
  free(ha_eqtime);
  free(ha_eqreg);
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

  forint indata[5];
  indata[1]=VecSize;
  indata[2]=mpisize;
  indata[3]=ndblock;
  indata[4]=nsbbdblocks;
  forint *ptx=NULL;
  ptx = indata;

  ha_cgetype *x1=NULL;//= (ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));
  ha_cgetype *xcf=NULL;
  ha_cgetype *x0=NULL;// (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  ha_cgetype *b1=NULL;//= (ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));
  extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
  extern void spec48_single_(forint *indata,int *irn, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
  extern void spec48_nomc66_(forint *indata, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow, MPI_Fint *fcomm,forint *rowptrin, forint *colptrin);

  forint k=0,m=1;
  ha_cgetype temp1,temp2;
  forint tindx1;
  printf("rank %d ncof %ld\n",rank,ncof);
  

  if(solmethod==10)Johansen(nohsl,VecSize,A,dnz,dnnz,onz,onnz,B,dnzB,dnnzB,onzB,onnzB,vecb,vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,&xcf);
  
  int stepcount;
  int nsteps=3;
  ha_cgetype vpercents=1.0,perprecis=0;
  FILE* solution;
  int maxsol=3;

    if(solmethod==1)ModMidPoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,0);
    
  if(solmethod==20){
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
      if(rank==rank_hsl)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,10*laA,subints,1,1,0);
      
      ModMidPoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,0);
  printf("ncof %ld rank %d\n",ncof,rank);
  }
  }
  if(solmethod==21){
    MPI_Barrier(PETSC_COMM_WORLD);
    printf("Heere rank %d\n",rank);
      ModMidPoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,2);
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
    fwrite(xcf, sizeof(ha_cgetype),nvarele, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".var");
    printf("solchar %s\n",solchar);
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(ha_var, sizeof(hcge_cof),nvar, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".set");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(ha_set, sizeof(ha_cgeset),nset, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".sel");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(ha_setele, sizeof(ha_cgesetele),nsetspace, solution);
    fclose(solution);
    uvadd modeldes[4];
    modeldes[0]=nsetspace;
    modeldes[1]=nvar;
    modeldes[2]=nvarele;
    modeldes[3]=(uvadd)nset;
    strcpy(solchar,tempchar);
    strcat(solchar,".mds");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(modeldes, sizeof(uvadd),4, solution);
    fclose(solution);
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  if(nowrites==0&&rank==0)for(i=0; i<noutdata; i++){
    hcge_wdata(tabfile,iodata[i+niodata].logname,iodata[i+niodata].filname,ha_set,nset,ha_setele,ha_cof,ncof,ncofele,ha_var,nvar,nvarele,ha_cofvar);
    printf("outfile %s\n",iodata[i+niodata].logname);
  }
  free(iodata);
  free(countvarintra1);
  free(ha_eqadd);
  free(ndbbddrank1);
  free(counteq);
  free(counteqnoadd);
  free(ha_set);
  free(ha_setele);
  free(ha_cof);
  free(ha_var);
  free(ha_cgeshock);
  free(ha_cofvar);
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
  MPI_Comm_free(&HA_COMM);
  MPI_Comm_free(&HA1_COMM);
  ierr = PetscFinalize();
  CHKERRQ(ierr);
  return 0;
}


