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
  extern void spec48_ssol_(forint *INSIZE,forint *IRN, forint *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
  extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
  extern void spec48_single_(forint *indata,int *irn, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
  extern void spec48_nomc66_(forint *indata, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow, MPI_Fint *fcomm,forint *rowptrin, forint *colptrin);
  extern void my_vec_comz_(ha_cgetype vecbivi,int *biviindx,int *col, int *row, int *colsize,int *nz0,int *nz1);

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

    if(solmethod==1) { //Modified midpoint Pearson 1991
              uvadd *counteqs= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
              uvadd *counteqnoadds= (uvadd *) calloc (ndblock,sizeof(uvadd));
              uvadd *countvarintra1s= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
              memcpy(counteqs,counteq,(ndblock+1)*sizeof(uvadd));
              memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(uvadd));
              memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(uvadd));
    gettimeofday(&begintime, NULL);
    ha_cgetype *xc0=NULL;
    ha_cgetype *xc12= NULL;
    ha_cgetype *xc24= NULL;
    int *xc124= NULL;
    ha_cgetype *clag1= NULL;
    ha_cgetype *varchange= NULL;
    for(subindx=0; subindx<subints; subindx++) {
      for(sol=0; sol<maxsol; sol++) {
        if(sol==0)nsteps=step1;
        if(sol==1) nsteps=(int)step1*kindx1;
        if(sol==2) nsteps=(int)step1*kindx2;
        vpercents=(ha_cgetype)100/nsteps;
        for(stepcount=0; stepcount<nsteps; stepcount++) {
          printf("subint %d sol %d stepcount %d nsteps %d\n",subindx,sol,stepcount,nsteps);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time0);
          CHKERRQ(ierr);
          if(stepcount==0) {
            MPI_Barrier(PETSC_COMM_WORLD);
            if(!(subindx==0&&sol==0&&stepcount==0)) {
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
            if(sol==0)for(i=0; i<ncofele; i++) {
                ha_cofvar[i].var0=ha_cofvar[i].varval;
              }
            else for(i=0; i<ncofele; i++) {
                ha_cofvar[i].varval=ha_cofvar[i].var0;
              }
            ha_cofvar1=ha_cofvar+ncofele;
            for(i=0; i<nvar; i++) {
              if(ha_var[i].change_real) {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                    ha_cofvar1[tindx1].csolpupd=ha_cgeshock[tindx1].ShockVal/nsteps;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                  }
                }
              }
              else {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                    temp2=ha_cgeshock[tindx1].ShockVal;//subints;
                    ha_cofvar1[tindx1].csolpupd=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;//ha_cgeshock[ha_var[i].begadd+j].ShockVal/nsteps;//(exp(log(1+ha_cgeshock[ha_var[i].begadd+j].ShockVal/100)/nsteps)-1)*100;
                    ha_cofvar1[tindx1].csolpupd*=vpercents;//nsteps*100;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                  }
                }
              }
            }
            MPI_Barrier(PETSC_COMM_WORLD);
            ierr = VecAssemblyBegin(vece);
            CHKERRQ(ierr);
            ierr = VecAssemblyEnd(vece);
            CHKERRQ(ierr);
          }
          if(rank==rank_hsl) {

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempclag1");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(clag1, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(clag1);
            clag1=NULL;

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvarchange");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(varchange, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(varchange);
            varchange=NULL;
          }
          MPI_Barrier(PETSC_COMM_WORLD);

          strcpy(commsyntax,"equation");

          if(nohsl) {
            MatCreate(PETSC_COMM_WORLD,&A);
          }
          else {
            MatCreate(PETSC_COMM_SELF,&A);
          }
          if(nesteddbbd==1)MatSetSizes(A,localsize,localsize,VecSize,VecSize);
          else MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
          if(nohsl) {
            MatSetType(A,MATMPIAIJ);
          }
          else {
            MatSetType(A,MATSEQAIJ);
          }
          if(nohsl) {
            MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
          }
          else {
            MatSeqAIJSetPreallocation(A,dnz,dnnz);
          }

          if(nohsl) {
            MatCreate(PETSC_COMM_WORLD,&B);
          }
          else {
            MatCreate(PETSC_COMM_SELF,&B);
          }
          if(nesteddbbd==1)MatSetSizes(B,localsize,localsize,VecSize,VecSize);
          else MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
          if(nohsl) {
            MatSetType(B,MATMPIAIJ);
          }
          else {
            MatSetType(B,MATSEQAIJ);
          }
          if(nohsl) {
            MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
          }
          else {
            MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
          }

          if(rank==rank_hsl) {
            HaNewMatVal(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqadd,counteq,nintraeq,A,B);
          }
          if(rank==rank_hsl) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvar");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
            fclose(tempvar);
            free(ha_cofvar);
            ha_cofvar=NULL;

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempshock");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
            fclose(tempvar);
            free(ha_cgeshock);
            ha_cgeshock=NULL;
          }

          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = VecDuplicate(vece,&vecb);
          CHKERRQ(ierr);
          if(rank==rank_hsl) {
            ierr = MatMult(B,vece,vecb);
            CHKERRQ(ierr);
          }
          ierr = VecDestroy(&vece);
          CHKERRQ(ierr);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyBegin(vecb);
          CHKERRQ(ierr);
          ierr = VecAssemblyEnd(vecb);
          CHKERRQ(ierr);
          ierr = MatDestroy(&B);
          CHKERRQ(ierr);
          if(matsol>=2) {
            int *ha_rows= (int *) calloc (VecSize,sizeof(int));
            int *ha_cols= (int *) calloc (VecSize,sizeof(int));
            int *ha_ndblocks= (int *) calloc (ndblock,sizeof(int));
            gettimeofday(&gettime_now, NULL);

            if(matsol==2) {
              HaDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,cntl6);
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              HaDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
            }

            if(matsol==3) {
              presol=1;
              memcpy(counteq,counteqs,(ndblock+1)*sizeof(uvadd));
              memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(uvadd));
              memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(uvadd));
              HaNDBBDMatOderPre(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
              HaNDBBDParPre(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
              presol=0;
              HaNDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              HaNDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
            }

            gettimeofday(&endtime, NULL);
            MPI_Barrier(PETSC_COMM_WORLD);
            ierr = PetscGetCPUTime(&time1);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
            if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - gettime_now.tv_sec)+((double)(endtime.tv_usec - gettime_now.tv_usec))/ 1000000);
            free(ha_rows);
            free(ha_cols);
            free(ha_ndblocks);
            printf("rank %d\n",rank);
            MPI_Barrier(PETSC_COMM_WORLD);
          }
          else {
            if(matsol==1) {
              if(rank==rank_hsl) {
                Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
                ai= aa->i;
                aj= aa->j;
                vals=aa->a;
                nz01=aa->nz;
                count=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    count++;
                  }
              }
              int *irn=(int *) calloc (count,sizeof(int));
              int *irn1=(int *) calloc (nz01,sizeof(int));
              int *jcn=(int *) calloc (count,sizeof(int));
              ha_cgetype *values= (ha_cgetype *) calloc (count,sizeof(ha_cgetype));
              if(rank==rank_hsl) {
                for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                    irn1[j]=i+1;
                  }
                for(j=ai[VecSize-1]; j<nz01; j++) {
                  irn1[j]=VecSize;
                }
                j=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    irn[j]=irn1[i];
                    jcn[j]=aj[i]+1;
                    values[j]=vals[i];
                    j++;
                  }
              }
              ierr = MatDestroy(&A);
              CHKERRQ(ierr);
              free(irn1);
              b1=realloc (b1,VecSize*sizeof(ha_cgetype));
              if(rank==rank_hsl) {
                VecGetArray(vecb,&vals);
                for(i=0; i<VecSize; i++) {
                  b1[i]=vals[i];
                }
              }
              ierr = VecDestroy(&vecb);
              CHKERRQ(ierr);
              int *neleperrow= (int *) calloc (VecSize,sizeof(int));
              int *ai1= (int *) calloc (VecSize,sizeof(int));
              if(rank==rank_hsl) {
                j=1;
                k=0,m=1;
                for(i=1; i<count; i++) {
                  if(irn[i]-irn[i-1]>0) {
                    neleperrow[k]=j;
                    ai1[k]=m;
                    j=1;
                    m=i+1;
                    k++;
                  }
                  else {
                    j++;
                  }
                }
                neleperrow[k]=j;
                ai1[k]=ai1[k-1]+neleperrow[k-1];
              }
              indata[1]=VecSize;//.m
              indata[0]=count;//.nz
              ptx = indata;
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
              free(irn);
              if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              free(jcn);
              free(values);
              free(neleperrow);
              free(ai1);
              free(b1);
              b1=NULL;
            }
            else {
              if(rank==rank_hsl) {
                Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
                ai= aa->i;
                aj= aa->j;
                vals=aa->a;
                nz01=aa->nz;
                count=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    count++;
                  }

              }
              lasize=ceil((laA/100)*count);
              int *irn=(int *) calloc (lasize,sizeof(int));
              int *irn1=(int *) calloc (nz01,sizeof(int));
              int *jcn=(int *) calloc (lasize,sizeof(int));
              ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
              if(rank==rank_hsl) {
                for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                    irn1[j]=i+1;
                  }
                for(j=ai[VecSize-1]; j<nz01; j++) {
                  irn1[j]=VecSize;
                }
                j=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    irn[j]=irn1[i];
                    jcn[j]=aj[i]+1;
                    values[j]=vals[i];
                    j++;
                  }
                VecGetArray(vecb,&vals);
              }
              ierr = MatDestroy(&A);
              CHKERRQ(ierr);
              free(irn1);
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              int *insize=(int *) calloc (4,sizeof(int));
              insize[0]=VecSize;
              insize[1]=VecSize;
              insize[2]=count;
              insize[3]=laA;
              if(rank==rank_hsl)spec48_ssol2la_(insize,irn,jcn,values,vals,x1);
              ierr = VecDestroy(&vecb);
              CHKERRQ(ierr);
              free(insize);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              free(irn);
              free(jcn);
              free(values);
            }
          }
          if(rank==rank_hsl) {
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            ha_cgeshock=realloc (ha_cgeshock,(nvarele)*sizeof(ha_cgeexovar));
            freadresult=fread(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvar");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            ha_cofvar=realloc (ha_cofvar,(ncofele+nvarele)*sizeof(ha_cgevar));
            freadresult=fread(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempclag1");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            clag1=realloc (clag1,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(clag1, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvarchange");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            varchange=realloc (varchange,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(varchange, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

          }
          MPI_Barrier(PETSC_COMM_WORLD);
          printf("sol %d stepcount %d\n\n",sol,stepcount);
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
          printf("rank %d check %d\n\n",rank,1);
          if(nesteddbbd==1)VecSetSizes(vece,localsize,VecSize);
          else VecSetSizes(vece,PETSC_DECIDE,VecSize);
          VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
          printf("rank %d check %d\n\n",rank,2);
          ha_cofvar1=ha_cofvar+ncofele;
          if(stepcount==0) {
            for(i=0; i<nvar; i++) {
              if(ha_var[i].change_real) {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    ha_cofvar1[tindx1].varval+=ha_cofvar1[tindx1].csolpupd;
                    varchange[tindx1]=ha_cofvar1[tindx1].csolpupd;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    varchange[tindx1]=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].varval+=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];
                    clag1[tindx1]=0;
                  }
                }
              }
              else {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    varchange[tindx1]=ha_cofvar1[tindx1].csolpupd;
                    ha_cofvar1[tindx1].varval*=(1+ha_cofvar1[tindx1].csolpupd/100);
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd/(1+ha_cofvar1[tindx1].csolpupd/100),INSERT_VALUES);
                  }
                  else {
                    varchange[tindx1]=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].varval*=(1+ha_cofvar1[tindx1].csolpupd/100);
                    clag1[tindx1]=0;
                  }
                }
              }
            }
          }
          else {
            for(i=0; i<nvar; i++) {
              if(ha_var[i].change_real) {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    ha_cofvar1[tindx1].varval+=ha_cofvar1[tindx1].csolpupd;
                    varchange[tindx1]+=ha_cofvar1[tindx1].csolpupd;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    temp1=ha_cofvar1[tindx1].varval;//change;
                    varchange[tindx1]=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx];//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx];//
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];//ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;
                    ha_cofvar1[tindx1].varval=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx];//ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;//ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd;
                    clag1[tindx1]=temp1;
                  }
                }
              }
              else {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    temp2=ha_cgeshock[tindx1].ShockVal;//subints;
                    temp1=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;
                    temp1*=vpercents;
                    ha_cofvar1[tindx1].csolpupd=temp1/(1+varchange[tindx1]/100);
                    varchange[tindx1]+=temp1;//*(1+ha_cofvar[ncofele+ha_var[i].begadd+j].varchange/100)
                    ha_cofvar1[tindx1].varval=(1+varchange[tindx1]/100)*ha_cofvar1[tindx1].var0;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,temp1/(1+varchange[tindx1]/100),INSERT_VALUES);
                  }
                  else {
                    temp1=varchange[tindx1];
                    varchange[tindx1]=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx]*(100+temp1)/100;//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]*(1+temp1/100);//
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];//(ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1)/(1+temp1/100);
                    ha_cofvar1[tindx1].varval=(100+varchange[tindx1])/100*ha_cofvar1[tindx1].var0;
                    clag1[tindx1]=temp1;
                  }
                }
              }
            }
          }
          printf("rank %d check %d\n\n",rank,3);
          free(x1);
          x1=NULL;
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyBegin(vece);
          CHKERRQ(ierr);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyEnd(vece);
          CHKERRQ(ierr);
          printf("rank %d check %d\n\n",rank,4);
          if(rank==rank_hsl) {
            if(stepcount==0) {
              hnew_update(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
            }
            else {
              hnew_mupdate(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
            }
            strcpy(commsyntax,"formula");
            IsIni=false;
            hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
          }
          printf("rank %d check %d\n\n",rank,5);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Update time %f\n",time1-time0);
          CHKERRQ(ierr);
          ierr = PetscGetCPUTime(&time0);
          CHKERRQ(ierr);
          printf("OKKL!\n");
        }

        strcpy(commsyntax,"equation");
        if(rank==rank_hsl) {

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempclag1");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
          }
          fwrite(clag1, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          free(clag1);
          clag1=NULL;

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvarchange");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
          }
          fwrite(varchange, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          free(varchange);
          varchange=NULL;

        }

        MPI_Barrier(PETSC_COMM_WORLD);
        if(nohsl) {
          MatCreate(PETSC_COMM_WORLD,&A);
        }
        else {
          MatCreate(PETSC_COMM_SELF,&A);
        }
        if(nesteddbbd==1)MatSetSizes(A,localsize,localsize,VecSize,VecSize);
        else MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
        if(nohsl) {
          MatSetType(A,MATMPIAIJ);
        }
        else {
          MatSetType(A,MATSEQAIJ);
        }
        if(nohsl) {
          MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
        }
        else {
          MatSeqAIJSetPreallocation(A,dnz,dnnz);
        }

        if(nohsl) {
          MatCreate(PETSC_COMM_WORLD,&B);
        }
        else {
          MatCreate(PETSC_COMM_SELF,&B);
        }
        if(nesteddbbd==1)MatSetSizes(B,localsize,localsize,VecSize,VecSize);
        else MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
        if(nohsl) {
          MatSetType(B,MATMPIAIJ);
        }
        else {
          MatSetType(B,MATSEQAIJ);
        }
        if(nohsl) {
          MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
        }
        else {
          MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
        }

        if(rank==rank_hsl) {
          HaNewMatVal(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqadd,counteq,nintraeq,A,B);
        }

        if(rank==rank_hsl) {
          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvar");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
          }
          fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
          fclose(tempvar);
          free(ha_cofvar);
          ha_cofvar=NULL;

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempshock");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
          }
          fwrite(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
          fclose(tempvar);
          free(ha_cgeshock);
          ha_cgeshock=NULL;
        }

        MPI_Barrier(PETSC_COMM_WORLD);
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = VecDuplicate(vece,&vecb);
        CHKERRQ(ierr);
        if(rank==rank_hsl) {
          ierr = MatMult(B,vece,vecb);
          CHKERRQ(ierr);
        }
        ierr = VecDestroy(&vece);
        CHKERRQ(ierr);
        ierr = VecAssemblyBegin(vecb);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(vecb);
        CHKERRQ(ierr);
        ierr = MatDestroy(&B);
        CHKERRQ(ierr);
        if(matsol>=2) {
          int *ha_rows= (int *) calloc (VecSize,sizeof(int));
          int *ha_cols= (int *) calloc (VecSize,sizeof(int));
          int *ha_ndblocks= (int *) calloc (ndblock,sizeof(int));
          gettimeofday(&gettime_now, NULL);

          if(matsol==2) {
            HaDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,cntl6);
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            HaDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
          }

          if(matsol==3) {
            presol=1;
            memcpy(counteq,counteqs,(ndblock+1)*sizeof(uvadd));
            memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(uvadd));
            memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(uvadd));
            HaNDBBDMatOderPre(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
            HaNDBBDParPre(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
            presol=0;
            HaNDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            HaNDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
          }

          gettimeofday(&endtime, NULL);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
          if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - gettime_now.tv_sec)+((double)(endtime.tv_usec - gettime_now.tv_usec))/ 1000000);
          free(ha_rows);
          free(ha_cols);
          free(ha_ndblocks);
          printf("rank %d\n",rank);
          MPI_Barrier(PETSC_COMM_WORLD);
        }
        else {
          if(matsol==1) {
            if(rank==rank_hsl) {
              Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
              ai= aa->i;
              aj= aa->j;
              vals=aa->a;
              nz01=aa->nz;
              count=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  count++;
                }
            }
            int *irn=(int *) calloc (count,sizeof(int));
            int *irn1=(int *) calloc (nz01,sizeof(int));
            int *jcn=(int *) calloc (count,sizeof(int));
            ha_cgetype *values= (ha_cgetype *) calloc (count,sizeof(ha_cgetype));
            if(rank==rank_hsl) {
              for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                  irn1[j]=i+1;
                }
              for(j=ai[VecSize-1]; j<nz01; j++) {
                irn1[j]=VecSize;
              }
              j=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  irn[j]=irn1[i];
                  jcn[j]=aj[i]+1;
                  values[j]=vals[i];
                  j++;
                }
            }
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);
            free(irn1);
            b1=realloc (b1,VecSize*sizeof(ha_cgetype));
            if(rank==rank_hsl) {
              VecGetArray(vecb,&vals);
              for(i=0; i<VecSize; i++) {
                b1[i]=vals[i];
              }
            }
            ierr = VecDestroy(&vecb);
            CHKERRQ(ierr);
            int *neleperrow= (int *) calloc (VecSize,sizeof(int));
            int *ai1= (int *) calloc (VecSize,sizeof(int));
            if(rank==rank_hsl) {
              j=1;
              k=0,m=1;
              for(i=1; i<count; i++) {
                if(irn[i]-irn[i-1]>0) {
                  neleperrow[k]=j;
                  ai1[k]=m;
                  j=1;
                  m=i+1;
                  k++;
                }
                else {
                  j++;
                }
              }
              neleperrow[k]=j;
              ai1[k]=ai1[k-1]+neleperrow[k-1];
            }
            indata[1]=VecSize;//.m
            indata[0]=count;//.nz
            ptx = indata;//&
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
            free(irn);
            if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            free(jcn);
            free(values);
            free(neleperrow);
            free(ai1);
            free(b1);
            b1=NULL;
          }
          else {
            if(rank==rank_hsl) {
              Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
              ai= aa->i;
              aj= aa->j;
              vals=aa->a;
              nz01=aa->nz;
              count=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  count++;
                }
            }
            lasize=ceil((laA/100)*count);
            int *irn=(int *) calloc (lasize,sizeof(int));
            int *irn1=(int *) calloc (nz01,sizeof(int));
            int *jcn=(int *) calloc (lasize,sizeof(int));
            ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
            if(rank==rank_hsl) {
              for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                  irn1[j]=i+1;
                }
              for(j=ai[VecSize-1]; j<nz01; j++) {
                irn1[j]=VecSize;
              }
              j=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  irn[j]=irn1[i];
                  jcn[j]=aj[i]+1;
                  values[j]=vals[i];
                  j++;
                }
              VecGetArray(vecb,&vals);
            }
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);
            free(irn1);
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            int *insize=(int *) calloc (4,sizeof(int));
            insize[0]=VecSize;
            insize[1]=VecSize;
            insize[2]=count;
            insize[3]=laA;
            if(rank==rank_hsl)spec48_ssol2la_(insize,irn,jcn,values,vals,x1);
            free(insize);
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            ierr = VecDestroy(&vecb);
            CHKERRQ(ierr);
            free(irn);
            free(jcn);
            free(values);
          }
        }
        if(rank==rank_hsl) {
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          ha_cgeshock=realloc (ha_cgeshock,(nvarele)*sizeof(ha_cgeexovar));
          freadresult=fread(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvar");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          ha_cofvar=realloc (ha_cofvar,(ncofele+nvarele)*sizeof(ha_cgevar));
          freadresult=fread(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempclag1");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          clag1=realloc (clag1,(nvarele)*sizeof(ha_cgetype));
          freadresult=fread(clag1, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvarchange");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          varchange=realloc (varchange,(nvarele)*sizeof(ha_cgetype));
          freadresult=fread(varchange, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

        }
        ha_cofvar1=ha_cofvar+ncofele;
        for(i=0; i<nvar; i++) {
          if(ha_var[i].change_real) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
              if(ha_cgeshock[tindx1].ShockId) {
                ha_cofvar1[tindx1].varval=0;
              }
              else {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+x1[ha_cgeshock[tindx1].ExoIndx]);
                ha_cofvar1[tindx1].varval=0;//ha_cofvar[tindx2].var0+varchange[tindx1];//no distortion between steps
                clag1[tindx1]=0;
              }
            }
          }
          else {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
              if(ha_cgeshock[tindx1].ShockId) {
                ha_cofvar1[tindx1].varval=0;
              }
              else {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+x1[ha_cgeshock[tindx1].ExoIndx]*(1+varchange[tindx1]/100));
                ha_cofvar1[tindx1].varval=0;//ha_cofvar[tindx2].varval*varchange[tindx1]/100;
                clag1[tindx1]=0;
              }
            }
          }
        }
        if(rank==rank_hsl) {
          if(subindx!=0||sol!=0) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcf");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xcf=realloc (xcf,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(xcf, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc12");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xc12=realloc (xc12,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(xc12, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc24");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xc24=realloc (xc24,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(xc24, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            
            xc0=realloc (xc0,(nvarele)*sizeof(ha_cgetype));
            if(subindx>0&&sol>0){
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcO");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            freadresult=fread(xc0, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            if(subindx==subints-1&&sol==maxsol-1)remove(tempfilenam);
            }

          }
          if(subindx==0&&sol==0) {
            xc0=realloc (xc0,nvarele*sizeof(ha_cgetype));
            xcf=realloc (xcf,nvarele*sizeof(ha_cgetype));
            for(i=0; i<nvarele; i++)xcf[i]=0;
          }
          if(sol==0)xc12=realloc (xc12,nvarele*sizeof(ha_cgetype));
          if(sol==0)xc24=realloc (xc24,nvarele*sizeof(ha_cgetype));
          if(subindx>0) {
            if(sol==0)for(i=0; i<nvarele; i++) xc0[i]=1+xcf[i]/100;
            if(sol==0) {
              kval1=1.0/(kindx1*kindx1-1.0);
              kval2=1.0/(1-kindx1*kindx1)/(1.0-kindx2*kindx2);
              for(i=0; i<nvar; i++) {
                if(ha_var[i].change_real) {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc12[k]=xcf[k]-varchange[k]*kval1;
                    xc24[k]=xcf[k];
                    xcf[k]+=varchange[k]*kval2;
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]=xcf[k];
                    xc12[k]=xcf[k]-varchange[k]*xc0[k]*kval1;
                    xcf[k]+=varchange[k]*kval2*xc0[k];
                  }
                }
              }
            }
            if(sol==1) {
              kval1=kindx1*kindx1/(kindx1*kindx1-1.0);
              kval2=kindx1*kindx1/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx1*kindx1*kindx1*kindx1/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx1*kindx1);
              for(i=0; i<nvar; i++) {
                if(ha_var[i].change_real) {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]-=varchange[k]*kval2;
                    xc12[k]+=varchange[k]*kval1;
                    xcf[k]-=varchange[k]*kval3;
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]-=varchange[k]*xc0[k]*kval2;
                    xc12[k]+=varchange[k]*xc0[k]*kval1;
                    xcf[k]-=varchange[k]*xc0[k]*kval3;
                  }
                }
              }
            }
            if(sol==2) {
              kval2=kindx2*kindx2/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx2*kindx2*kindx2*kindx2/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx2*kindx2);
              for(i=0; i<nvar; i++) {
                if(ha_var[i].change_real) {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]+=varchange[k]*kval2;
                    xcf[k]+=varchange[k]*kval3;
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]+=varchange[k]*xc0[k]*kval2;
                    xcf[k]+=varchange[k]*xc0[k]*kval3;
                  }
                }
              }
            }
          }
          else {
            if(sol==0){
              kval1=1.0/(kindx1*kindx1-1);
              kval2=1.0/(1-kindx1*kindx1)/(1-kindx2*kindx2);
              for(i=0; i<nvarele; i++) {
                xc12[i]=-varchange[i]*kval1;
                xcf[i]+=varchange[i]*kval2;
              }
            }
            if(sol==1) {
              kval1=kindx1*kindx1/(kindx1*kindx1-1.0);
              kval2=kindx1*kindx1/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx1*kindx1*kindx1*kindx1/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx1*kindx1);
              for(i=0; i<nvarele; i++) {
                xc24[i]=-varchange[i]*kval2;
                xc12[i]+=varchange[i]*kval1;
                xcf[i]-=varchange[i]*kval3;
              }
            }
            if(sol==2) {
              kval2=kindx2*kindx2/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx2*kindx2*kindx2*kindx2/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx2*kindx2);
              for(i=0; i<nvarele; i++) {
                xc24[i]+=varchange[i]*kval2;
                xcf[i]+=varchange[i]*kval3;
              }
            }
          }

          if(sol==maxsol-1){
          if(subindx==0){
          for(i=0; i<nvar; i++) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
               ha_cofvar1[tindx1].csolpupd=xcf[tindx1];
            }
          }
          }else{
          for(i=0; i<nvar; i++) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
               ha_cofvar1[tindx1].csolpupd=(100+xcf[tindx1])/xc0[tindx1]-100;
            }
          }
          }
          for(i=0; i<ncofele; i++) ha_cofvar[i].varval=ha_cofvar[i].var0;
          hnew_gupd(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
          strcpy(commsyntax,"formula");
          IsIni=false;
          hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
          for(i=0; i<nvar; i++) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
               ha_cofvar1[tindx1].csolpupd=0;
            }
          }

          if(nohsl)MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Last Update time %f\n",time1-time0);
          CHKERRQ(ierr);
          }

          
          if(!(subindx==subints-1&&sol==maxsol-1)) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcf");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xcf, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(xcf);
            xcf=NULL;
            
            if(subindx>0&&sol==0){            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcO");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xc0, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            }
            free(xc0);
            xc0=NULL;
          }

    if(rank==rank_hsl&&sol==maxsol-1) {
            if(subindx==0){
              xc124=realloc (xc124,nvarele*sizeof(int));
              for(i=0; i<nvarele; i++)xc124[i]=6;
            }else{
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc124");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xc124=realloc (xc124,nvarele*sizeof(int));
            freadresult=fread(xc124, sizeof(int),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }

      for(i=0; i<nvarele; i++){
        j=0;
        if(xc12[i]>0){
        while (xc12[i] >= 10){
          xc12[i] /= 10;
          j++;
        }
        }else{
        while (xc12[i] <= -10){
          xc12[i] /= 10;
          j++;
        }
        }
        xc24[i]/=pow(10,j);
        j=abs(floor((xc12[i]-xc24[i])*100000));
        if(j!=0){
          if(j<10){
            if(xc124[i]>5)xc124[i]=5;
          }else{
            if(j<100){
              if(xc124[i]>4)xc124[i]=4;
            }else {
              if(j<1000){
                if(xc124[i]>3)xc124[i]=3;
              }else {
                if(j<10000){
                  if(xc124[i]>2)xc124[i]=2;
                }else{
                  if(xc124[i]>1)xc124[i]=1;
                }
              }
            }
          }
        }
      }
      if(subindx!=subints-1){
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc124");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xc124, sizeof(int),nvarele, tempvar);
            fclose(tempvar);
            free(xc124);
            xc124=NULL;
      }
    }
            
          if(!(subindx==subints-1&&sol==maxsol-1)) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc12");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xc12, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(xc12);
            xc12=NULL;
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc24");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xc24, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(xc24);
            xc24=NULL;
          }
          
          
          strcpy(solchar,temdir);
          if(subindx<10)strcat(solchar,"_tempsol0");
          else strcat(solchar,"_tempsol");
          sprintf(tempchar, "%d", subindx);
          strcat(solchar,tempchar);
          sprintf(tempchar, "%d", sol);
          strcat(solchar,tempchar);
          strcat(solchar,".bin");
          printf("solchar %s\n",solchar);
          if ( (solution = fopen(solchar, "wb")) == NULL ) {
            printf("Error opening file\n");
            return 1;
          }
          fwrite(xcf, sizeof(ha_cgetype),nvarele, solution);
          fclose(solution);
        }
        free(x1);
        x1=NULL;
      }
    }
    long int *precis= (long int *) calloc (6,sizeof(long int));
    
    if(rank==rank_hsl) {
      for(i=0; i<nvarele; i++){
   switch(xc124[i]) {
      case 6:
         precis[5]+=1;
         break;
      case 5:
         precis[4]+=1;
         break;
      case 4:
         precis[3]+=1;
         break;
      case 3:
         precis[2]+=1;
         break;
      case 2:
         precis[1]+=1;
         break;
      default :
         precis[0]+=1;
   }
    }
    if(rank==0)printf("Accurate at 6 digits        %ld\nAccurate at 5 digits        %ld\nAccurate at 4 digits        %ld\nAccurate at 3 digits        %ld\nAccurate at 2 digits        %ld\nAccurate at 1 digit or none %ld\n",precis[5],precis[4],precis[3],precis[2],precis[1],precis[0]);
    }
    free(precis);
    xc0=realloc (xc0,sizeof(ha_cgetype));
    free(xc0);
    xc0=NULL;
    free(xc12);
    xc12=NULL;
    free(xc24);
    xc24=NULL;
    free(xc124);
    xc124=NULL;
    free(clag1);
    free(varchange);
    gettimeofday(&endtime, NULL);
    if(rank==0)printf("Mmid method calculation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
              free(counteqs);
              free(counteqnoadds);
              free(countvarintra1s);
  }
    
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


