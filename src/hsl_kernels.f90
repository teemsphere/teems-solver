module constants
  implicit none
  INTEGER, PARAMETER :: DPC = KIND(0.0D0)
  INTEGER, PARAMETER :: FSORD = 1
  !switch to single by setting FSORD = 0, DPC and hsl_mp48d.f90 accordingly
contains
  ! Verbosity from the C side (main.c exports TEEMS_VERBOSITY alongside
  ! TEEMS_SCRATCH). Informational kernel output is gated at >= 2; error
  ! output is never gated. Called once per kernel entry, never in loops.
  integer function teems_verbosity()
    character(len=8) :: s
    integer :: st
    teems_verbosity=1
    call get_environment_variable("TEEMS_VERBOSITY",s,status=st)
    if (st==0) then
      read(s,*,iostat=st) teems_verbosity
      if (st/=0) teems_verbosity=1
    end if
  end function teems_verbosity
end module constants

SUBROUTINE SPEC48_SINGLE(indata,irn1,jcn1,b1,values1,x,neleperrow,ai1,fcomm)
  USE HSL_MP48_DOUBLE
  use HSL_mc66_double
  use constants
  IMPLICIT NONE

  integer (8) indata(*)
  integer (4) fcomm(*),fcomm1
  integer(4) irn1(*),jcn1(*),ai1(*),neleperrow(*)
  real (kind=DPC) b1(*),values1(*),x(*)
  ! mc66 controller:
  type (mc66_control) :: control

  ! random number seed
  type (fa14_seed) :: seed

  !nz: number of nonzeros
  !n: number of rows.
  !m: number of columns
  !irn: row indices of the matrix
  !jcn: column indices of the matrix
  integer nz,n,m

  ! rowptr: rowptr[1:nblocks] is the starting row index for
  !         each diagonal block. rowptr[nblocks+1]=m
  ! colptr: colptr[1:nblocks] is the starting column index
  !         for each diagonal block. colptr[nblocks+1]
  !         is the starting column index of the border.
  integer, pointer :: rowptr(:)
  integer, pointer :: colptr(:)

  ! column order: column_order(i) is the original column index
  !   of the i-column of the reordered matrix
  ! row ordering: row_order(i) is the original row index
  !   of the i-row of the reordered matrix
  ! row_position: row_position(i) is the new row index
  !   of the original row i.
  ! column_position: column_position(i) is the new column index
  !   of the original column i.
  integer, pointer :: row_order(:),column_order(:),column_position(:) !,row_position
  ! nblocks: number of blocks
  ! netcut: column size of size of border
  ! info: info tag
  ! kblocks: the actual number of blocks in the SBBD form
  ! rowdiff: row dimension imbalance in percentage term
  integer nblocks,netcut,info,kblocks
  DOUBLE PRECISION rowdiff
  integer i,j,k,h,o,p
  TYPE (MP48_DATA) data
  INTEGER ERCODE,ST
  INTEGER (4) ZDIAG
  LOGICAL FLAG
  integer myid, numprocs
  fcomm1=fcomm(1)

  ! Program to illustrate use of MP48.
  ! Define a communicator for the package
  data%COMM = fcomm1 !MPI_COMM_WORLD
  ! Initialize package
  data%JOB = 1
  CALL MP48AD(data)
  ! Reset control parameters (if required)
  ! Read all values on host
  data%ICNTL(7) = 3
  IF (data%RANK.EQ.0) THEN

    nz = indata(1)!indata%nz
    m = indata(2)!indata%m
    n = m
    IF (indata(5).GE.2) then !%nsbbdblocks
      nblocks = indata(5)!%nsbbdblocks
    ELSE
      nblocks = 2
    ENDIF
    allocate(row_order(m),rowptr(nblocks+1), &
      column_order(n), colptr(nblocks+1))
    control%COARSEN_SCHEME=2
    if (teems_verbosity()>=2) write(*,"(i10,i10,i10,i10)") nz,m,n,nblocks
    call mc66(m,n,nz,irn1,jcn1,nblocks,control,seed, &
      row_order,info,rowptr,column_order,colptr,&
      netcut,rowdiff,kblocks)
    if (info /= 0) then
      call mc66_print_message(info)
      if (info < 0) stop "mc66 failed"
    end if
    if (teems_verbosity()>=2) write(*,"(a,I10)") "netcut = ",netcut
    if (teems_verbosity()>=2) write(*,"(a,f10.2,'%')") "row dimension imbalance = ", rowdiff
    ! dump block sizes
    ! reorder the original matrix
    allocate(column_position(n))!, row_position(m))
    do i = 1, n
      column_position(column_order(i)) = i
    end do


    data%NEQ=m !indata%x
    data%NBLOCK=nblocks !indata%y
    data%NE=nz !indata%z
    ALLOCATE(data%NEQSB(1:data%NBLOCK),STAT=ST )
    ALLOCATE(data%EQPTR(1:data%NEQ+1),STAT=ST )
    ALLOCATE(data%EQVAR(1:data%NE),STAT=ST )
    ALLOCATE(data%VALUES(1:data%NE),STAT=ST )
    ALLOCATE(data%B(1:data%NEQ),STAT=ST )
    ! Read matrix data on host.
    do i = 1, nblocks
      data%NEQSB(i)=rowptr(i+1)-rowptr(i)!neqsb(i)
      if (teems_verbosity()>=2) write(*,"('block ',i4,' of dimension  ',i10,' X ',i10)") &
      i,rowptr(i+1)-rowptr(i),colptr(i+1)-colptr(i)
    end do
    h=1
    do j = 1, m
      data%EQPTR(j)=h
      do i = 1,neleperrow(row_order(j))
        data%EQVAR(h)=column_position(jcn1(ai1(row_order(j))+i-1))
        data%VALUES(h)=values1(ai1(row_order(j))+i-1)
        h=h+1
      end do
    end do
    data%EQPTR(m+1)=nz+1
    ! Also read right hand side
    do i = 1, m
      data%B(i)=b1(row_order(i))
    end do
  END IF
  CALL MPI_BARRIER(data%COMM,ERCODE)
  call MPI_COMM_RANK( fcomm1 , myid, ERCODE ) !MPI_COMM_WORLD
  call MPI_COMM_SIZE( fcomm1, numprocs, ERCODE ) !MPI_COMM_WORLD
  if (teems_verbosity()>=2) print *, "Process ", myid, " of ", numprocs, " is alive1"
  data%JOB = 25
  CALL MP48AD(data)
  IF (data%RANK.EQ.0) THEN
    IF (data%ERROR.LT.0 .OR. data%ERROR.EQ.2) THEN
      ! fail fast with named diagnosis for the singularity classes,
      ! mirroring SPEC48_NOMC66 (see the comment there)
      WRITE (6,'(A,I4)') 'Error STOP from MP48 with data%ERROR = ',data%ERROR
      IF (data%ERROR.EQ.-21) WRITE (6,*) ' (the SBBD interface matrix is structurally rank deficient)'
      IF (data%ERROR.EQ.2) WRITE (6,*) ' (MP48 reports the matrix is singular)'
      IF (data%ERROR.EQ.-21 .OR. data%ERROR.EQ.2) THEN
        ZDIAG=0
        CALL TEEMS_ONFAIL_DIAG(ZDIAG)
      END IF
      CALL TEEMS_ONFAIL_ABORT()
    ELSE
      do o=1,m
        x(column_order(o))=data%X(o)
      end do
    END IF
    deallocate(row_order,rowptr,column_order, colptr)
  END IF
  call MPI_COMM_RANK( fcomm1, myid, ERCODE ) !MPI_COMM_WORLD
  call MPI_COMM_SIZE( fcomm1, numprocs, ERCODE ) !MPI_COMM_WORLD
  if (teems_verbosity()>=2) print *, "Process ", myid, " of ", numprocs, " is alive2"
  data%JOB = 6
  CALL MP48AD(data)
END SUBROUTINE SPEC48_SINGLE

SUBROUTINE SPEC48_NOMC66(indata,jcn1,b1,values1,x,neleperrow,fcomm,rowptrin,colptrin)
  USE HSL_MP48_DOUBLE
  use constants
  IMPLICIT NONE

  integer myid, numprocs
  integer (8) :: i,j,o,h
  integer (8) :: nblocks,maxsbcols,ncols

  TYPE (MP48_DATA) data
  INTEGER ERCODE,ST
  INTEGER (4) ZDIAG
  LOGICAL FLAG
  integer (4) fcomm(*),fcomm1
  integer(4) neleperrow(*),jcn1(*)
  integer (8) indata(*),rowptrin(*),colptrin(*)
  real (kind=DPC) b1(*),values1(*),x(*)
  integer (8) :: nz
  integer (8) :: m
  fcomm1=fcomm(1)

  ! rowptr: rowptr[1:nblocks] is the starting row index for
  !         each diagonal block. rowptr[nblocks+1]=m
  ! colptr: colptr[1:nblocks] is the starting column index
  !         for each diagonal block. colptr[nblocks+1]
  !         is the starting column index of the border.

  data%COMM = fcomm1 !MPI_COMM_WORLD
  ! Initialize package
  data%JOB = 1
  CALL MP48AD(data)
  ! Reset control parameters (if required)
  ! Read all values on host
  data%ICNTL(7) = 3
  IF (data%RANK.EQ.0) THEN
    nz = indata(1)!%nz
    m = indata(2)!%m
    nblocks = indata(4)!%nblock
    data%NEQ=indata(2)!%m !indata%x
    data%NBLOCK=nblocks !indata%y
    data%NE=indata(1)!%nz !indata%z
    ALLOCATE(data%NEQSB(1:data%NBLOCK),STAT=ST )
    ALLOCATE(data%EQPTR(1:data%NEQ+1),STAT=ST )
    ALLOCATE(data%EQVAR(1:data%NE),STAT=ST )
    ALLOCATE(data%VALUES(1:data%NE),STAT=ST )
    ALLOCATE(data%B(1:data%NEQ),STAT=ST )
    ! Read matrix data on host.
    maxsbcols=0
    do i = 1, nblocks
      data%NEQSB(i)=rowptrin(i+1)-rowptrin(i)!neqsb(i)
      if (teems_verbosity()>=2) write(*,"('nomc block ',i4,' of dimension  ',i10,' X ',i10)") &
      i,data%NEQSB(i),colptrin(i+1)-colptrin(i)
      ncols=colptrin(i+1)-colptrin(i)
      if(maxsbcols.LT.ncols)maxsbcols=ncols
    end do
    maxsbcols=maxsbcols+m-colptrin(nblocks+1)
    if (teems_verbosity()>=2) print *, "row ", m,"nz ",nz,"maxcolsb ",maxsbcols
    h=1
    data%MAXSBCOLS=maxsbcols
    do j = 1, m
      data%EQPTR(j)=h
      do i = 1,neleperrow(j)
        data%EQVAR(h)=jcn1(h)
        data%VALUES(h)=values1(h)
        h=h+1
      end do
    end do
    data%EQPTR(m+1)=h
    ! Also read right hand side
    do i = 1, m
      data%B(i)=b1(i)
    end do
  END IF
  CALL MPI_BARRIER(data%COMM,ERCODE)
  call MPI_COMM_RANK( fcomm1 , myid, ERCODE ) !MPI_COMM_WORLD
  call MPI_COMM_SIZE( fcomm1, numprocs, ERCODE ) !MPI_COMM_WORLD
  if (teems_verbosity()>=2) print *, "Process ", myid, " of ", numprocs, " is alive1"
  data%JOB = 25
  CALL MP48AD(data)
  if (teems_verbosity()>=2) print *, "Processa ", myid, " of ", numprocs, " is alive1"
  IF (data%RANK.EQ.0) THEN
    IF (data%ERROR.LT.0 .OR. data%ERROR.EQ.2) THEN
      ! fail fast: the old path printed one line and CONTINUED with a
      ! stale x.  ERROR=-21 (interface matrix structurally rank
      ! deficient) and the ERROR=+2 warning (matrix is singular, the
      ! "solution" is not usable) are singularity-class: name the
      ! defects.  Other codes are setup/allocation/IO.
      WRITE (6,'(A,I4)') 'Error STOP from MP48 with data%ERROR = ',data%ERROR
      IF (data%ERROR.EQ.-21) WRITE (6,*) ' (the SBBD interface matrix is structurally rank deficient)'
      IF (data%ERROR.EQ.2) WRITE (6,*) ' (MP48 reports the matrix is singular)'
      IF (data%ERROR.EQ.-21 .OR. data%ERROR.EQ.2) THEN
        ZDIAG=0
        CALL TEEMS_ONFAIL_DIAG(ZDIAG)
      END IF
      CALL TEEMS_ONFAIL_ABORT()
    ELSE
      do o=1,m
        x(o)=data%X(o)
      end do
    END IF
  END IF
  if (teems_verbosity()>=2) print *, "Process ", myid, " of ", numprocs, " is alive2"
  data%JOB = 6
  CALL MP48AD(data)
END SUBROUTINE SPEC48_NOMC66


module mp48_persist
  ! State for the persistent SBBD instance (-fastrefac): the MP48_DATA
  ! structure (border lists, per-block pivot sequences, factors) lives
  ! across steps so repeat steps refill VALUES/B and refactorize with
  ! FACT_JOB=2 instead of rebuilding the whole instance.
  USE HSL_MP48_DOUBLE
  implicit none
  TYPE (MP48_DATA), save :: pdata
  logical, save :: pinit=.false.
end module mp48_persist

SUBROUTINE SPEC48_NOMC66_P(indata,jcn1,b1,values1,x,neleperrow,fcomm,rowptrin,colptrin,redo)
  ! Persistent-instance variant of SPEC48_NOMC66 (structural partition).
  ! redo(1) on entry: 0 = (re)build the instance (JOB 1,2,3 then 4 with
  ! FACT_JOB=1, then 5); 1 = values-only step (host refills VALUES/B,
  ! JOB=4 with FACT_JOB=2, JOB=5).  On exit: 0 = solved, <0 = MP48
  ! error on the fast path (caller retries with redo=0).  All ranks
  ! must call with the same redo.  Caller keeps jcn1/neleperrow
  ! unchanged between redo=1 calls and refills values1/b1 per step.
  use mp48_persist
  use constants
  IMPLICIT NONE
  integer myid, numprocs
  integer (8) :: i,j,o,h
  integer (8) :: nblocks,maxsbcols,ncols
  INTEGER ERCODE,ST
  INTEGER (4) ZDIAG
  integer (4) fcomm(*),fcomm1
  integer (4) redo(*)
  integer(4) neleperrow(*),jcn1(*)
  integer (8) indata(*),rowptrin(*),colptrin(*)
  real (kind=DPC) b1(*),values1(*),x(*)
  integer (8) :: nz
  integer (8) :: m
  fcomm1=fcomm(1)
  IF (redo(1).EQ.1 .AND. pinit) THEN
    ! values-only step: refactorize on the kept structure
    IF (pdata%RANK.EQ.0) THEN
      nz = indata(1)
      m = indata(2)
      do i = 1, nz
        pdata%VALUES(i)=values1(i)
      end do
      do i = 1, m
        pdata%B(i)=b1(i)
      end do
    END IF
    pdata%FACT_JOB=2
    pdata%JOB = 4
    CALL MP48AD(pdata)
    IF (pdata%ERROR.LT.0) THEN
      if (teems_verbosity()>=1 .AND. pdata%RANK.EQ.0) WRITE (6,'(A,I4)') &
        'Note: fast SBBD refactorize declined, MP48 code ',pdata%ERROR
      redo(1)=pdata%ERROR
      RETURN
    END IF
    IF (pdata%ERROR.EQ.2) THEN
      ! MP48 "matrix is singular" warning: the solution would be
      ! garbage — surface as -21 so the caller diagnoses and aborts
      IF (pdata%RANK.EQ.0) WRITE (6,*) ' MP48 reports the matrix is singular'
      redo(1)=-21
      RETURN
    END IF
    pdata%JOB = 5
    CALL MP48AD(pdata)
    IF (pdata%ERROR.LT.0) THEN
      redo(1)=pdata%ERROR
      RETURN
    END IF
    IF (pdata%RANK.EQ.0) THEN
      do o=1,indata(2)
        x(o)=pdata%X(o)
      end do
    END IF
    redo(1)=0
    RETURN
  END IF
  ! full (re)build
  IF (pinit) THEN
    pdata%JOB = 6
    CALL MP48AD(pdata)
    pinit=.false.
  END IF
  pdata%COMM = fcomm1
  pdata%JOB = 1
  CALL MP48AD(pdata)
  pdata%ICNTL(7) = 3
  ! the instance persists across steps: own the solution vector so
  ! repeated JOB=5 calls never re-allocate package-owned storage
  pdata%ICNTL(13) = 1
  IF (pdata%RANK.EQ.0) THEN
    nz = indata(1)
    m = indata(2)
    nblocks = indata(4)
    pdata%NEQ=indata(2)
    pdata%NBLOCK=nblocks
    pdata%NE=indata(1)
    ! JOB=6 frees only package-allocated arrays; the user-supplied
    ! components survive it and must be released before a rebuild
    IF (ALLOCATED(pdata%NEQSB)) DEALLOCATE(pdata%NEQSB)
    IF (ALLOCATED(pdata%EQPTR)) DEALLOCATE(pdata%EQPTR)
    IF (ALLOCATED(pdata%EQVAR)) DEALLOCATE(pdata%EQVAR)
    IF (ALLOCATED(pdata%VALUES)) DEALLOCATE(pdata%VALUES)
    IF (ALLOCATED(pdata%B)) DEALLOCATE(pdata%B)
    IF (ALLOCATED(pdata%X)) DEALLOCATE(pdata%X)
    ALLOCATE(pdata%NEQSB(1:pdata%NBLOCK),STAT=ST )
    ALLOCATE(pdata%EQPTR(1:pdata%NEQ+1),STAT=ST )
    ALLOCATE(pdata%EQVAR(1:pdata%NE),STAT=ST )
    ALLOCATE(pdata%VALUES(1:pdata%NE),STAT=ST )
    ALLOCATE(pdata%B(1:pdata%NEQ),STAT=ST )
    ALLOCATE(pdata%X(1:pdata%NEQ),STAT=ST )
    maxsbcols=0
    do i = 1, nblocks
      pdata%NEQSB(i)=rowptrin(i+1)-rowptrin(i)
      if (teems_verbosity()>=2) write(*,"('nomc block ',i4,' of dimension  ',i10,' X ',i10)") &
      i,pdata%NEQSB(i),colptrin(i+1)-colptrin(i)
      ncols=colptrin(i+1)-colptrin(i)
      if(maxsbcols.LT.ncols)maxsbcols=ncols
    end do
    maxsbcols=maxsbcols+m-colptrin(nblocks+1)
    if (teems_verbosity()>=2) print *, "row ", m,"nz ",nz,"maxcolsb ",maxsbcols
    h=1
    pdata%MAXSBCOLS=maxsbcols
    do j = 1, m
      pdata%EQPTR(j)=h
      do i = 1,neleperrow(j)
        pdata%EQVAR(h)=jcn1(h)
        pdata%VALUES(h)=values1(h)
        h=h+1
      end do
    end do
    pdata%EQPTR(m+1)=h
    do i = 1, m
      pdata%B(i)=b1(i)
    end do
  END IF
  CALL MPI_BARRIER(pdata%COMM,ERCODE)
  pdata%JOB = 23
  CALL MP48AD(pdata)
  IF (pdata%ERROR.LT.0) THEN
    WRITE (6,*) ' Unexpected MP48 analyse code on rank ',pdata%RANK
    redo(1)=pdata%ERROR
    RETURN
  END IF
  pdata%FACT_JOB=1
  pdata%JOB = 4
  CALL MP48AD(pdata)
  IF (pdata%ERROR.LT.0) THEN
    WRITE (6,*) ' Unexpected MP48 factorize code on rank ',pdata%RANK
    redo(1)=pdata%ERROR
    RETURN
  END IF
  IF (pdata%ERROR.EQ.2) THEN
    ! singular-matrix warning on the rebuild: fatal (see fast leg)
    IF (pdata%RANK.EQ.0) WRITE (6,*) ' MP48 reports the matrix is singular'
    redo(1)=-21
    RETURN
  END IF
  pinit=.true.
  pdata%JOB = 5
  CALL MP48AD(pdata)
  IF (pdata%ERROR.LT.0) THEN
    redo(1)=pdata%ERROR
    RETURN
  END IF
  IF (pdata%RANK.EQ.0) THEN
    do o=1,indata(2)
      x(o)=pdata%X(o)
    end do
  END IF
  redo(1)=0
END SUBROUTINE SPEC48_NOMC66_P

SUBROUTINE SPEC48_NOMC66_PFREE()
  use mp48_persist
  IMPLICIT NONE
  IF (pinit) THEN
    pdata%JOB = 6
    CALL MP48AD(pdata)
    pinit=.false.
  END IF
  ! user-supplied components survive JOB=6
  IF (ALLOCATED(pdata%NEQSB)) DEALLOCATE(pdata%NEQSB)
  IF (ALLOCATED(pdata%EQPTR)) DEALLOCATE(pdata%EQPTR)
  IF (ALLOCATED(pdata%EQVAR)) DEALLOCATE(pdata%EQVAR)
  IF (ALLOCATED(pdata%VALUES)) DEALLOCATE(pdata%VALUES)
  IF (ALLOCATED(pdata%B)) DEALLOCATE(pdata%B)
  IF (ALLOCATED(pdata%X)) DEALLOCATE(pdata%X)
END SUBROUTINE SPEC48_NOMC66_PFREE

SUBROUTINE SPEC51M_RANK(INSIZE,CNTL6,IRN,JCN,VA,IRNA,JCNA,KEEP,W,IW)
  use constants
  IMPLICIT NONE

  integer M,N,NE
  integer (4) INSIZE(*)
  integer (4) JCN(*),IRN(*),JCNA(*),IRNA(*),KEEP(*),IW(*)
  real (kind=DPC) VA(*),CNTL6(*),W(*)
  integer I,LA, MAXN,RANK1,SGNDET,T,NEFAC
  real (kind=DPC) LOGDET
  real (kind=DPC), pointer :: CNTL(:),RINFO(:)!,W(:),A(:),
  integer, pointer :: ICNTL(:),INFO(:)!,IW(:),KEEP(:),JCN1(:),IRN1(:),COLS(:),ROWS(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(5)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(6)) then
    LA=INSIZE(6)
  end if
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10))!A(LA),,W(5*MAXN)


  allocate(ICNTL(20),INFO(20))!JCN1(LA),IRN1(LA),,ROWS(M)COLS(N),,IW(6*M+3*N)

!     Factorize matrix
    IF (FSORD.EQ.1) THEN
      CALL MA48ID(CNTL,ICNTL)
    ELSE
      CALL MA48I(CNTL,ICNTL)
    ENDIF
    ! errors only below debug verbosity (silences duplicate-entry notes)
    if (teems_verbosity()<2) ICNTL(3)=1
    ! -3 workspace shortfalls are handled by caller-side growth, so
    ! MA48's own 'Error return ... because LA is' print would read as
    ! fatal in the logs; the wrappers print the diagnosis on genuinely
    ! fatal INFO(1) codes
    ICNTL(1)=0
    IF(CNTL6(1).EQ.0)THEN
      IF (FSORD.EQ.1) THEN
        CNTL(4)=1e-4
      else
        CNTL(4)=0.3
      ENDIF
    ELSE
      CNTL(4)=CNTL6(1)!1e-4!0.0000000001
    ENDIF
    ICNTL(7)=0
    ! INSIZE(7) (out): 0 = done, -3 = MA48 workspace too small
    ! (INSIZE(8) then holds the suggested LA; caller reallocates,
    ! re-stages IRN/JCN/VA and retries with INSIZE(6) = the new size)
    INSIZE(7)=0
    IF (FSORD.EQ.1) THEN
      CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      IF (INFO(1).EQ.-3) THEN
        INSIZE(7)=-3
        INSIZE(8)=max(INFO(3),INFO(4))
        deallocate(CNTL,RINFO,ICNTL,INFO)
        RETURN
      END IF
      CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
      IF (INFO(1).EQ.-3) THEN
        INSIZE(7)=-3
        INSIZE(8)=INFO(4)
        deallocate(CNTL,RINFO,ICNTL,INFO)
        RETURN
      END IF

!     Compute the determinant
      CALL MA51CD(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ELSE
      CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      IF (INFO(1).EQ.-3) THEN
        INSIZE(7)=-3
        INSIZE(8)=max(INFO(3),INFO(4))
        deallocate(CNTL,RINFO,ICNTL,INFO)
        RETURN
      END IF
      CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
      IF (INFO(1).EQ.-3) THEN
        INSIZE(7)=-3
        INSIZE(8)=INFO(4)
        deallocate(CNTL,RINFO,ICNTL,INFO)
        RETURN
      END IF
      CALL MA51C(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ENDIF

    IF(SGNDET.GT.0)THEN
      if (teems_verbosity()>=2) WRITE(6,'(/,A,F9.3)')&
          ' Determinant is positive; log(determinant) =',LOGDET
          INSIZE(4)=M
    ELSE IF(SGNDET.LT.0)THEN
        if (teems_verbosity()>=2) WRITE(6,'(/,A,F9.3)')&
          ' Determinant is negative; log(-determinant) =',LOGDET
          INSIZE(4)=M
      ELSE
        if (teems_verbosity()>=2) WRITE(6,'(/,A)')' Determinant is zero'
!     Determine the nonsingular submatrix of the factorization
        IF (FSORD.EQ.1) THEN
          CALL MA51AD(M,N,LA,IRN,KEEP,RANK1,IRNA,JCNA,W)!ROWS,COLS
        ELSE
          CALL MA51A(M,N,LA,IRN,KEEP,RANK1,IRNA,JCNA,W)!,ROWS,COLS
        ENDIF
        INSIZE(4)=RANK1
    END IF
  deallocate(CNTL,RINFO)!A,,W
  deallocate(ICNTL,INFO)!IRN1,JCN1,,ROWSCOLS,,KEEP,IW
END SUBROUTINE SPEC51M_RANK

!SUBROUTINE SPEC48_MSOL2(INSIZE,IRN,JCN,VA,B,X,IRNC,JCNC,VAC,IRNV,JCNV,VAV)
!  IMPLICIT NONE
!
!!  write(*,"(i5,i5)") IRN(1),JCN(1)
!!  write(*,"(i5,i5)") IRN(2),JCN(2)
!!  write(*,"(i5,i5)") IRN(3),JCN(3)
!!  write(*,"(i5,i5)") IRN(4),JCN(4)
!
!!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)
!!  write(*,"(i5,i5)") IRN(4),JCN(4)
!!  write(*,"(i5,i5)") IRN(4),JCN(4)
!    !write(*,"(A,i5,i5,F20.2)") 'Ain',IRN1(I),JCN1(i),VA(I)
!! Scale matrix
!    !write(*,"(A,i5,i5,F10.2,F10.2,F10.2)") 'Ainbs',I,J,R(I),C(J),A(II)
!    !write(*,"(A,i5,i5,F20.2,F20.2)") 'Ain',I,J,R(I)*C(J),A(II)
!!     Factorize matrix
!  !ICNTL(4)=0
!  !ICNTL(6)=1
!  !CNTL(4)=0.00000000000000001
!  !CNTL(2)=0
!
!  !do I=1,NE
!  !  write(*,"(A,i5,i5,i5,F20.2)") 'AAD',I,IRN1(I),JCN1(I),A(I)
!  !end do
!                RINFO)
!
!          ' Determinant is positive; log(determinant) =',LOGDET
!          ' Determinant is negative; log(-determinant) =',LOGDET
!
!    !write(*,"(i5,F10.2)") I,B(I)
!              RHS,SOL,ERROR,W,IW,INFO)
!    !write(*,"(A,i5,F20.2)") 'UI',I,X(I)
!    !I=J
!    !write(*,"(A,i10)") 'UI',I
!    !I=NC
!    !write(*,"(A,i10)") 'NC',I
!        !if(I.EQ.NEC) then
!          !write(*,"(A,i5,i5,F20.2)") 'RHS',IRNC(I),JCNC(I),VAC(I)
!        !end if
!        !if(RHS(J).NE.0) then
!        !end if
!              RHS,SOL,ERROR,W,IW,INFO)
!        !if(L.GE.9450) then
!          !if(J.EQ.7) then
!            !write(*,"(A,i5,i5,F20.10)") 'V',J,I,VAV(L)
!          !end if
!        !end if
!END SUBROUTINE SPEC48_MSOL2


SUBROUTINE SPEC48M_MSOL(INSIZE,IRN,JCN,VA,B,X,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0,KEEP)!,IRNV,JCNV,VAV
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZROW0(*),BIVINZCOL0(*),KEEP(*)!,JCNV(*),IRNV(*)
  integer(4)IRNB(*),JCNB(*)
  integer(4) JCNC(*),IRNC(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1!,NV,MV,NEV
  real(kind=DPC) VA(*),VAC(*),B(*),X(*),VALUESB(*),VECBIVI(*)
  !DOUBLE PRECISION LOGDET,SGNDET!,VAV(*)real (8)
  integer LA, MAXN,NBIVI,MBIVI
  integer(4) I,J,L,L1,L2,L3,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),SOL(:)!,VOUT(:)!,RHS(:)A(:),
  integer, pointer :: ICNTL(:),INFO(:),IW(:),JCNB1(:)!,IRNOUT(:),JCNOUT(:)!,JCN1(:),IRN1(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  MC=INSIZE(4)
  NC=INSIZE(5)
  NEC=INSIZE(6)
  MB=INSIZE(7)
  NB=INSIZE(8)
  NEB=INSIZE(9)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  MBIVI=INSIZE(14)
  NBIVI=INSIZE(15)
  JCN(1:NE)=JCN(1:NE)+1
  allocate(JCNB1(NEB))
  JCNB1(1:NEB)=JCNB(1:NEB)+1
  do J=1,M
    L=J+NE
    IRN((IRN(L)+1):IRN(L+1))=J
  end do
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(17)) then
    LA=INSIZE(17)
  end if
  ! INSIZE(18) (out): 0 = factorized, -3 = workspace too small
  ! (INSIZE(19) then holds MA48's suggested LA; caller reallocates and
  ! re-stages IRN/JCN/VA -- clobbered by the conversion above and by
  ! MA48 -- and retries with INSIZE(17) = the new size)
  INSIZE(18)=0
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3),SOL(M))!,RHS(M),W(5*MAXN)A(LA),

  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))!JCN1(LA),IRN1(LA),
!     Factorize matrix
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  ELSE
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  ! errors only below debug verbosity (silences duplicate-entry notes)
  if (teems_verbosity()<2) ICNTL(3)=1
  ! -3 workspace shortfalls are handled by caller-side growth, so
  ! MA48's own 'Error return ... because LA is' print would read as
  ! fatal in the logs; the wrappers print the diagnosis on genuinely
  ! fatal INFO(1) codes
  ICNTL(1)=0
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  ! KEEP is caller-allocated (bound M+9*N+7, i.e. ICNTL(6)=1); factors
  ! (VA/IRN/KEEP) return to C in place -- no scratch files
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ELSE
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).EQ.-3) THEN
    INSIZE(18)=-3
    INSIZE(19)=max(INFO(3),INFO(4))
    deallocate(CNTL,RINFO,ERROR1,SOL,JCNB1)
    deallocate(ICNTL,INFO,IW)
    RETURN
  END IF
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  CALL TEEMS_ONFAIL_DIAG(INFO(1))
  CALL TEEMS_ONFAIL_ABORT()
  END IF
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).EQ.-3) THEN
    INSIZE(18)=-3
    INSIZE(19)=INFO(4)
    deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)
    deallocate(ICNTL,INFO,IW)
    RETURN
  END IF
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    CALL TEEMS_ONFAIL_DIAG(INFO(1))
    CALL TEEMS_ONFAIL_ABORT()
  END IF
!     Compute the determinant

  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
  TRANS = .FALSE.
  if(sum(abs(B(1:N))).GT.0) then
    IF (FSORD.EQ.1) THEN
      CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
    ELSE
      CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
    ENDIF
  end if
  L=0
  TRANS = .FALSE.
  do J=1,MC-1!NC
      J3=IRNC(J+1)-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3!NEC
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ENDIF
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  end do

  J=MC
      J3=NEC-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3!NEC
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ELSE
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ENDIF
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)!,RHSA,
  deallocate(ICNTL,INFO,IW)!IRN1,,JCNOUT
END SUBROUTINE SPEC48M_MSOL


SUBROUTINE SPEC48M_MSOL_P(INSIZE,IRN,JCN,VA,B,X,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0,KEEP,REDO)
  ! Persistent-factor variant of SPEC48M_MSOL for the DBBD diagonal
  ! blocks (-fastrefac).  REDO(1) on entry: 0 = full analyse+factorize
  ! (IRN/JCN conversions + MA48A/AD + MA48B/BD JOB=1, as SPEC48M_MSOL);
  ! 1 = VA(1:NE) refilled by the caller in CSR order, IRN/JCN/KEEP are
  ! the persisted MA48 state, factorize with MA48B/BD JOB=2.  On exit:
  ! 0 = factorized and border solves done, <0 = fast factorize declined
  ! (MA48B/BD INFO(1)) BEFORE any border work — caller re-stages the
  ! block and retries with 0.
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZROW0(*),BIVINZCOL0(*),KEEP(*),REDO(*)
  integer(4)IRNB(*),JCNB(*)
  integer(4) JCNC(*),IRNC(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1
  real(kind=DPC) VA(*),VAC(*),B(*),X(*),VALUESB(*),VECBIVI(*)
  integer LA, MAXN,NBIVI,MBIVI
  integer(4) I,J,L,L1,L2,L3,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),SOL(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:),JCNB1(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  MC=INSIZE(4)
  NC=INSIZE(5)
  NEC=INSIZE(6)
  MB=INSIZE(7)
  NB=INSIZE(8)
  NEB=INSIZE(9)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  MBIVI=INSIZE(14)
  NBIVI=INSIZE(15)
  allocate(JCNB1(NEB))
  JCNB1(1:NEB)=JCNB(1:NEB)+1
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(17)) then
    LA=INSIZE(17)
  end if
  ! INSIZE(18) (out): 0 = factorized, -3 = workspace too small
  ! (INSIZE(19) then holds MA48's suggested LA; caller reallocates the
  ! persisted arrays, re-stages and retries with REDO(1)=0 and
  ! INSIZE(17) = the new size); REDO(1) is also set to -3 so the
  ! caller's decline channel sees the failure
  INSIZE(18)=0
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3),SOL(M))
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  ELSE
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  ! errors only below debug verbosity (silences duplicate-entry notes)
  if (teems_verbosity()<2) ICNTL(3)=1
  ! -3 workspace shortfalls are handled by caller-side growth, so
  ! MA48's own 'Error return ... because LA is' print would read as
  ! fatal in the logs; the wrappers print the diagnosis on genuinely
  ! fatal INFO(1) codes
  ICNTL(1)=0
  ! blocks whose entries turn unsuitable for the kept pivot sequence
  ! are refactorized as on a JOB=1 call
  ICNTL(11)=1
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  IF (REDO(1).EQ.0) THEN
    JCN(1:NE)=JCN(1:NE)+1
    do J=1,M
      L=J+NE
      IRN((IRN(L)+1):IRN(L+1))=J
    end do
    IF (FSORD.EQ.1) THEN
      CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
    ELSE
      CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
    ENDIF
    IF (INFO(1).EQ.-3) THEN
      INSIZE(18)=-3
      INSIZE(19)=max(INFO(3),INFO(4))
      REDO(1)=-3
      deallocate(CNTL,RINFO,ERROR1,SOL,JCNB1)
      deallocate(ICNTL,INFO,IW)
      RETURN
    END IF
    IF (INFO(1).LT.0) THEN
    WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
    CALL TEEMS_ONFAIL_DIAG(INFO(1))
    CALL TEEMS_ONFAIL_ABORT()
    END IF
    deallocate(IW)
    allocate(W(M),IW(2*M+2*N))
    IF (FSORD.EQ.1) THEN
      CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
    else
      CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
    ENDIF
    IF (INFO(1).EQ.-3) THEN
      INSIZE(18)=-3
      INSIZE(19)=INFO(4)
      REDO(1)=-3
      deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)
      deallocate(ICNTL,INFO,IW)
      RETURN
    END IF
    IF (INFO(1).NE.0) THEN
      WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
      INFO(1),'Solution not possible'
      write(*,"(A,i5)") 'RANK',INFO(5)
      CALL TEEMS_ONFAIL_DIAG(INFO(1))
      CALL TEEMS_ONFAIL_ABORT()
    END IF
  ELSE
    deallocate(IW)
    allocate(W(M),IW(2*M+2*N))
    IF (FSORD.EQ.1) THEN
      CALL MA48BD(M,N,NE,2,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
    else
      CALL MA48B(M,N,NE,2,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
    ENDIF
    IF (INFO(1).EQ.-3) THEN
      ! workspace too small for the kept factors: report -3 so the
      ! caller grows the arrays and redoes the analyse
      INSIZE(18)=-3
      INSIZE(19)=INFO(4)
      REDO(1)=-3
      deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)
      deallocate(ICNTL,INFO,IW)
      RETURN
    END IF
    IF (INFO(1).LT.0) THEN
      if (teems_verbosity()>=1) WRITE (6,'(A,I3)') &
        'Note: fast block refactorize declined, MA48B/BD INFO(1) =',INFO(1)
      REDO(1)=INFO(1)
      deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)
      deallocate(ICNTL,INFO,IW)
      RETURN
    END IF
  END IF

  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))
  TRANS = .FALSE.
  if(sum(abs(B(1:N))).GT.0) then
    IF (FSORD.EQ.1) THEN
      CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
    ELSE
      CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
    ENDIF
  end if
  L=0
  TRANS = .FALSE.
  do J=1,MC-1
      J3=IRNC(J+1)-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ENDIF
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)
              DO j2=IRNB(I)+1,IRNB(I+1)
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)
              do j2=IRNB(MB)+1,NEB
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  end do

  J=MC
      J3=NEC-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ELSE
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      ENDIF
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)
              DO j2=IRNB(I)+1,IRNB(I+1)
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)
              do j2=IRNB(MB)+1,NEB
                !$OMP ATOMIC
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  REDO(1)=0
  deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)
  deallocate(ICNTL,INFO,IW)
END SUBROUTINE SPEC48M_MSOL_P


SUBROUTINE SPEC48M_ESOL(INSIZE,IRN,VA,KEEP,B,SOL)
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) INSIZE(*),IRN(*),KEEP(*)
  integer M,N,NE
  real (kind=DPC) VA(*),B(*),SOL(*)
  integer LA, MAXN,I
  real (kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:)
  LOGICAL TRANS
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(17)) then
    LA=INSIZE(17)
  end if
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))
  allocate(ICNTL(20),INFO(20))

  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  ! errors only below debug verbosity (silences duplicate-entry notes)
  if (teems_verbosity()<2) ICNTL(3)=1
  JOB=1
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  else
    CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  ENDIF
  deallocate(CNTL,RINFO,W,ERROR1)
  deallocate(ICNTL,INFO,IW)
END SUBROUTINE SPEC48M_ESOL


SUBROUTINE SPEC48M_RPESOL(INSIZE,IRN,VA,KEEP,B,SOL,CNTL,RINFO,ERROR1,ICNTL,INFO,W,IW)
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) INSIZE(*),IRN(*),KEEP(*),ICNTL(*),INFO(*),IW(*)
  integer M,N,NE
  real(kind=DPC) VA(*),B(*),SOL(*)
  real(kind=DPC) CNTL(*),RINFO(*),W(*),ERROR1(*)
  integer LA!, MAXN,I
  LOGICAL TRANS
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
    LA=ceiling((NEFAC/100.0)*NE)
  if(LA.NE.INSIZE(17)) then
    LA=INSIZE(17)
  end if

  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  ! errors only below debug verbosity (silences duplicate-entry notes)
  if (teems_verbosity()<2) ICNTL(3)=1
  JOB=1
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  else
    CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  ENDIF
END SUBROUTINE SPEC48M_RPESOL



SUBROUTINE SPEC48_SSOL2LA(INSIZE,IRN,JCN,VA,B,X)
  ! One-shot LU solve.  INSIZE: 1=M 2=N 3=NE 4=LA percent (LA computed
  ! as ceiling(pct/100*NE) unless INSIZE(6)>0 gives an explicit LA),
  ! 5 (out): 0 = solved, -3 = workspace too small (INSIZE(6) then
  ! holds MA48's suggested LA; caller reallocates IRN/JCN/VA, refills
  ! them -- MA48 clobbers all three in place -- and retries with
  ! INSIZE(6) = the new size); other MA48 failures abort here.
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*)
  integer(4) INSIZE(*)
  integer M,N,NE,T
  real(kind=DPC) VA(*),B(*),X(*)
  integer I,J,L,LA, MAXN,SGNDET!,IFAIL,LP,II
  LOGICAL TRANS,checksol
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),LOGDET!A(:),,R(:),C(:)!,RHS(:),SOL(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:)!,JCN1(:),IRN1(:)
  ! -condest (MA60/MC71 solve-quality diagnostics; active only when the
  ! LU wrappers gate it on -- this kernel also serves the DBBD interface)
  LOGICAL CONDEST_ON
  integer(4) CDSTATUS,CDNOITER,CDKASE,CDIT
  integer(4) ICNTL60(5),KEEP60(10),JOB60(2)
  real(kind=DPC) OMEGA60(2),COND60(2),ERX60,RKEEP60(10)
  integer(4), allocatable :: CIRN(:),CJCN(:),C60IW(:)
  real(kind=DPC), allocatable :: CVA(:),CB(:),CX(:),CY(:),CRX(:),CD(:),C60W(:)
  integer(4), external :: TEEMS_CONDEST_ACTIVE
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(INSIZE(6).GT.0) then
    LA=INSIZE(6)
  endif
  INSIZE(5)=0
  ! pristine copy for the -condest diagnostics: MA48 clobbers the
  ! staged arrays in place, and MA60 needs the original system
  CONDEST_ON=.FALSE.
  IF (FSORD.EQ.1 .AND. TEEMS_CONDEST_ACTIVE().NE.0) THEN
    CDSTATUS=1
    DO I=1,M
      IF (B(I).NE.0.0) THEN
        CDSTATUS=0
        EXIT
      END IF
    END DO
    IF (CDSTATUS.EQ.1) THEN
      ! zero right-hand side (null-shock step): nothing to measure
      OMEGA60=0.0D0
      COND60=0.0D0
      ERX60=0.0D0
      CDNOITER=0
      CALL TEEMS_CONDEST_REPORT(CDSTATUS,OMEGA60(1),OMEGA60(2),ERX60,&
                                COND60(1),COND60(2),CDNOITER)
    ELSE
      CONDEST_ON=.TRUE.
      allocate(CIRN(NE),CJCN(NE),CVA(NE),CB(M))
      DO I=1,NE
        CIRN(I)=IRN(I)
        CJCN(I)=JCN(I)
        CVA(I)=VA(I)
      END DO
      DO I=1,M
        CB(I)=B(I)
      END DO
    END IF
  END IF
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))!,RHS(M),SOL(M),R(MAXN),C(MAXN),W(5*MAXN)A(LA),
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))!JCN1(LA),IRN1(LA),
!     Factorize matrix
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  ! errors only below debug verbosity (silences duplicate-entry notes)
  if (teems_verbosity()<2) ICNTL(3)=1
  ! -3 workspace shortfalls are handled by caller-side growth, so
  ! MA48's own 'Error return ... because LA is' print would read as
  ! fatal in the logs; the wrappers print the diagnosis on genuinely
  ! fatal INFO(1) codes
  ICNTL(1)=0
  T=M+5*N+4*N/ICNTL(6)+7
  allocate(KEEP(T))

  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).EQ.-3) THEN
    ! workspace too small: hand the suggested size back for a
    ! caller-side reallocate-and-retry
    INSIZE(5)=-3
    INSIZE(6)=max(INFO(3),INFO(4))
    deallocate(CNTL,RINFO,ERROR1,ICNTL,INFO,IW,KEEP)
    IF (CONDEST_ON) deallocate(CIRN,CJCN,CVA,CB)
    RETURN
  END IF
  ! the analyse itself fails on structural singularity (INFO(1)=-5);
  ! previously unchecked, which corrupted the follow-on factorize
  IF (INFO(1).LT.0) THEN
    WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
    CALL TEEMS_ONFAIL_DIAG(INFO(1))
    CALL TEEMS_ONFAIL_ABORT()
  END IF
  if (teems_verbosity()>=2) WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  if (teems_verbosity()>=2) WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
  IF (INFO(1).EQ.-3) THEN
    INSIZE(5)=-3
    INSIZE(6)=INFO(4)
    deallocate(CNTL,RINFO,W,ERROR1,ICNTL,INFO,IW,KEEP)
    IF (CONDEST_ON) deallocate(CIRN,CJCN,CVA,CB)
    RETURN
  END IF
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    WRITE (6,FMT='(A,i10)') 'INFO(5) =',INFO(5)
    CALL TEEMS_ONFAIL_DIAG(INFO(1))
    CALL TEEMS_ONFAIL_ABORT()
  END IF

  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(M,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
  else
    CALL MA48C(M,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
  endif
  if (teems_verbosity()>=2) WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  IF (CONDEST_ON) THEN
    ! MA60 iterative refinement + Arioli-Demmel-Duff backward/forward
    ! error and scaled condition numbers on the pristine copy, using
    ! the kept factors for the inv(A)/inv(A^T) products (MA48CD, trans
    ! per KASE).  Refinement runs on a COPY of X: the solution the
    ! caller receives is bit-identical with the flag on or off.
    allocate(CX(M),CY(M),CRX(M),CD(M),C60W(3*M),C60IW(2*M))
    DO I=1,M
      CX(I)=X(I)
      CD(I)=1.0D0
    END DO
    CALL MA60ID(ICNTL60,KEEP60,RKEEP60)
    ! no direct MA60 prints: the reporter owns the log line
    ICNTL60(1)=0
    JOB60(1)=1
    JOB60(2)=0
    CDKASE=0
    CDSTATUS=0
    DO CDIT=1,200
      CALL MA60AD(M,NE,CVA,CIRN,CJCN,CB,CX,CY,CD,C60W,C60IW,CDKASE,&
                  OMEGA60,ERX60,JOB60,COND60,CDNOITER,ICNTL60,KEEP60,RKEEP60)
      IF (CDKASE.EQ.0) EXIT
      IF (CDKASE.LT.0) THEN
        CDSTATUS=2
        IF (CDKASE.EQ.-3) CDSTATUS=3
        EXIT
      END IF
      TRANS=(CDKASE.EQ.1)
      CALL MA48CD(M,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  CY,CRX,ERROR1,W,IW,INFO)
      DO I=1,M
        CY(I)=CRX(I)
      END DO
    END DO
    CALL TEEMS_CONDEST_REPORT(CDSTATUS,OMEGA60(1),OMEGA60(2),ERX60,&
                              COND60(1),COND60(2),CDNOITER)
    deallocate(CX,CY,CRX,CD,C60W,C60IW)
    deallocate(CIRN,CJCN,CVA,CB)
  END IF
  deallocate(CNTL,RINFO,W,ERROR1)!A,,RHS,SOL
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,
END SUBROUTINE SPEC48_SSOL2LA

module ma48_persist
  ! State for the persistent-factor sequential LU path (-fastrefac).
  ! The caller keeps IRN/JCN/VA alive between calls; everything MA48
  ! needs beyond those (pivot sequence in KEEP, controls, workspaces)
  ! lives here so repeat steps can skip MA48A/AD entirely.
  use constants
  implicit none
  real(kind=DPC), allocatable, save :: pCNTL(:),pRINFO(:),pW(:),pERROR1(:)
  integer, allocatable, save :: pICNTL(:),pINFO(:),pIW(:),pKEEP(:)
  logical, save :: pready=.false.
  ! -condest pattern copy: on fast (JOB=2) steps IRN/JCN hold MA48
  ! state, not the original indices, so the pristine pattern is saved
  ! at rebuild time (VA(1:NE) is refilled in staging order each step,
  ! so values pair with this saved pattern)
  integer(4), allocatable, save :: pCIRN(:),pCJCN(:)
end module ma48_persist

SUBROUTINE SPEC48_SSOL2LA_P(INSIZE,IRN,JCN,VA,B,X)
  ! Sequential LU solve with a persistent pivot sequence.
  ! INSIZE: 1=M 2=N 3=NE 4=LA (absolute array length), 5 (inout): on
  ! entry 0 = full analyse+factorize, 1 = reuse the pivot sequence
  ! (MA48B/BD JOB=2); on exit 0 = solved, -3 = workspace too small
  ! (INSIZE(6) holds MA48's suggested LA; caller reallocates, refills
  ! and retries with 0), other <0 = fast factorize declined (MA48B/BD
  ! INFO(1)); caller refills IRN/JCN/VA and retries with 0.
  ! VA(1:NE) must be refilled by the caller each step in the same
  ! order it was filled for the analyse call; IRN/JCN and VA beyond NE
  ! are MA48 state and must not be touched between calls.
  use ma48_persist
  IMPLICIT NONE
  integer(4) INSIZE(*),IRN(*),JCN(*)
  real(kind=DPC) VA(*),B(*),X(*)
  integer M,N,NE,LA,MAXN,T
  LOGICAL TRANS
  ! -condest locals (see SPEC48_SSOL2LA; pattern copy lives in the module)
  LOGICAL CONDEST_ON
  integer(4) CDSTATUS,CDNOITER,CDKASE,CDIT,I
  integer(4) ICNTL60(5),KEEP60(10),JOB60(2)
  real(kind=DPC) OMEGA60(2),COND60(2),ERX60,RKEEP60(10)
  integer(4), allocatable :: C60IW(:)
  real(kind=DPC), allocatable :: CVA(:),CB(:),CX(:),CY(:),CRX(:),CD(:),C60W(:)
  integer(4), external :: TEEMS_CONDEST_ACTIVE
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  LA=INSIZE(4)
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  ! -condest: copy values and rhs before MA48 clobbers VA; the pattern
  ! is saved (module) on rebuild entries, where IRN/JCN are original
  CONDEST_ON=.FALSE.
  IF (FSORD.EQ.1 .AND. TEEMS_CONDEST_ACTIVE().NE.0) THEN
    CDSTATUS=1
    DO I=1,M
      IF (B(I).NE.0.0) THEN
        CDSTATUS=0
        EXIT
      END IF
    END DO
    IF (CDSTATUS.EQ.1) THEN
      OMEGA60=0.0D0
      COND60=0.0D0
      ERX60=0.0D0
      CDNOITER=0
      CALL TEEMS_CONDEST_REPORT(CDSTATUS,OMEGA60(1),OMEGA60(2),ERX60,&
                                COND60(1),COND60(2),CDNOITER)
    ELSE IF (INSIZE(5).NE.0 .AND. pready .AND. .NOT.allocated(pCIRN)) THEN
      ! fast step but no saved pattern (flag mid-run edge): no estimate
      CDSTATUS=2
      OMEGA60=0.0D0
      COND60=0.0D0
      ERX60=0.0D0
      CDNOITER=0
      CALL TEEMS_CONDEST_REPORT(CDSTATUS,OMEGA60(1),OMEGA60(2),ERX60,&
                                COND60(1),COND60(2),CDNOITER)
    ELSE
      CONDEST_ON=.TRUE.
      allocate(CVA(NE),CB(M))
      DO I=1,NE
        CVA(I)=VA(I)
      END DO
      DO I=1,M
        CB(I)=B(I)
      END DO
    END IF
  END IF
  IF (INSIZE(5).EQ.0 .OR. .NOT.pready) THEN
    IF (CONDEST_ON) THEN
      ! rebuild entry: IRN/JCN hold the original pattern -- save it
      if(allocated(pCIRN)) deallocate(pCIRN,pCJCN)
      allocate(pCIRN(NE),pCJCN(NE))
      DO I=1,NE
        pCIRN(I)=IRN(I)
        pCJCN(I)=JCN(I)
      END DO
    END IF
    if(allocated(pCNTL)) deallocate(pCNTL,pRINFO,pW,pERROR1,pICNTL,pINFO,pIW,pKEEP)
    allocate(pCNTL(10),pRINFO(10),pW(4*MAXN),pERROR1(3))
    allocate(pICNTL(20),pINFO(20),pIW(6*M+3*N))
    IF (FSORD.EQ.1) THEN
      CALL MA48ID(pCNTL,pICNTL)
    else
      CALL MA48I(pCNTL,pICNTL)
    endif
    ! errors only below debug verbosity (silences duplicate-entry notes)
    if (teems_verbosity()<2) pICNTL(3)=1
    ! -3 workspace shortfalls are handled by caller-side growth, so
    ! MA48's own 'Error return ... because LA is' print would read as
    ! fatal in the logs; the wrappers print the diagnosis on genuinely
    ! fatal INFO(1) codes
    pICNTL(1)=0
    ! any block whose entries turn unsuitable for the kept pivot
    ! sequence is refactorized as on a JOB=1 call instead of failing
    pICNTL(11)=1
    T=M+5*N+4*N/pICNTL(6)+7
    allocate(pKEEP(T))
    IF (FSORD.EQ.1) THEN
      CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,pKEEP,pCNTL,pICNTL,pIW,pINFO,pRINFO)
    else
      CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,pKEEP,pCNTL,pICNTL,pIW,pINFO,pRINFO)
    endif
    IF (pINFO(1).EQ.-3) THEN
      ! workspace too small: hand the suggested size back for a
      ! caller-side reallocate-and-retry
      INSIZE(5)=-3
      INSIZE(6)=max(pINFO(3),pINFO(4))
      pready=.false.
      IF (CONDEST_ON) deallocate(CVA,CB)
      RETURN
    END IF
    IF (pINFO(1).LT.0) THEN
      WRITE (6,'(A,I3)') 'Fatal STOP from MA48A/AD with INFO(1) =',pINFO(1)
      CALL TEEMS_ONFAIL_DIAG(pINFO(1))
      CALL TEEMS_ONFAIL_ABORT()
    END IF
    IF (FSORD.EQ.1) THEN
      CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,pKEEP,pCNTL,pICNTL,pW,pIW,pINFO,&
                  pRINFO)
    else
      CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,pKEEP,pCNTL,pICNTL,pW,pIW,pINFO,&
                  pRINFO)
    endif
    IF (pINFO(1).EQ.-3) THEN
      INSIZE(5)=-3
      INSIZE(6)=pINFO(4)
      pready=.false.
      IF (CONDEST_ON) deallocate(CVA,CB)
      RETURN
    END IF
    IF (pINFO(1).NE.0) THEN
      WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
      pINFO(1),'Solution not possible'
      WRITE (6,FMT='(A,i10)') 'INFO(5) =',pINFO(5)
      CALL TEEMS_ONFAIL_DIAG(pINFO(1))
      CALL TEEMS_ONFAIL_ABORT()
    END IF
    pready=.true.
  ELSE
    IF (FSORD.EQ.1) THEN
      CALL MA48BD(M,N,NE,2,LA,VA,IRN,JCN,pKEEP,pCNTL,pICNTL,pW,pIW,pINFO,&
                  pRINFO)
    else
      CALL MA48B(M,N,NE,2,LA,VA,IRN,JCN,pKEEP,pCNTL,pICNTL,pW,pIW,pINFO,&
                  pRINFO)
    endif
    IF (pINFO(1).EQ.-3) THEN
      INSIZE(5)=-3
      INSIZE(6)=pINFO(4)
      pready=.false.
      IF (CONDEST_ON) deallocate(CVA,CB)
      RETURN
    END IF
    IF (pINFO(1).LT.0) THEN
      ! fast factorize declined: caller refills the arrays and retries
      ! with a fresh analyse
      if (teems_verbosity()>=1) WRITE (6,'(A,I3)') &
        'Note: fast refactorize declined, MA48B/BD INFO(1) =',pINFO(1)
      INSIZE(5)=pINFO(1)
      IF (CONDEST_ON) deallocate(CVA,CB)
      RETURN
    END IF
  END IF
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(M,N,TRANS,1,LA,VA,IRN,pKEEP,pCNTL,pICNTL,&
              B,X,pERROR1,pW,pIW,pINFO)
  else
    CALL MA48C(M,N,TRANS,1,LA,VA,IRN,pKEEP,pCNTL,pICNTL,&
              B,X,pERROR1,pW,pIW,pINFO)
  endif
  IF (CONDEST_ON) THEN
    ! MA60 on the values/rhs copy + the module-saved pattern, solving
    ! with the persistent factors; X itself is never updated (see
    ! SPEC48_SSOL2LA)
    allocate(CX(M),CY(M),CRX(M),CD(M),C60W(3*M),C60IW(2*M))
    DO I=1,M
      CX(I)=X(I)
      CD(I)=1.0D0
    END DO
    CALL MA60ID(ICNTL60,KEEP60,RKEEP60)
    ICNTL60(1)=0
    JOB60(1)=1
    JOB60(2)=0
    CDKASE=0
    CDSTATUS=0
    DO CDIT=1,200
      CALL MA60AD(M,NE,CVA,pCIRN,pCJCN,CB,CX,CY,CD,C60W,C60IW,CDKASE,&
                  OMEGA60,ERX60,JOB60,COND60,CDNOITER,ICNTL60,KEEP60,RKEEP60)
      IF (CDKASE.EQ.0) EXIT
      IF (CDKASE.LT.0) THEN
        CDSTATUS=2
        IF (CDKASE.EQ.-3) CDSTATUS=3
        EXIT
      END IF
      TRANS=(CDKASE.EQ.1)
      CALL MA48CD(M,N,TRANS,1,LA,VA,IRN,pKEEP,pCNTL,pICNTL,&
                  CY,CRX,pERROR1,pW,pIW,pINFO)
      DO I=1,M
        CY(I)=CRX(I)
      END DO
    END DO
    CALL TEEMS_CONDEST_REPORT(CDSTATUS,OMEGA60(1),OMEGA60(2),ERX60,&
                              COND60(1),COND60(2),CDNOITER)
    deallocate(CX,CY,CRX,CD,C60W,C60IW)
    deallocate(CVA,CB)
  END IF
  INSIZE(5)=0
END SUBROUTINE SPEC48_SSOL2LA_P

SUBROUTINE SPEC48_PERSIST_FREE()
  use ma48_persist
  IMPLICIT NONE
  if(allocated(pCNTL)) deallocate(pCNTL,pRINFO,pW,pERROR1,pICNTL,pINFO,pIW,pKEEP)
  if(allocated(pCIRN)) deallocate(pCIRN,pCJCN)
  pready=.false.
END SUBROUTINE SPEC48_PERSIST_FREE

SUBROUTINE SPEC48M_SSOL2LA(INSIZE,IRN,JCN,VA,B,X)
  ! One-shot LU solve, explicit-LA variant (NDBBD final interface).
  ! INSIZE: 1=M 2=N 3=NE 4=LA percent, 6=explicit LA (always used),
  ! 5 (out): 0 = solved, -3 = workspace too small (INSIZE(6) then
  ! holds MA48's suggested LA; caller reallocates IRN/JCN/VA, refills
  ! them -- MA48 clobbers all three in place -- and retries);
  ! other MA48 failures abort here.
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*)
  integer(4) INSIZE(*)
  integer M,N,NE,T
  real(kind=DPC) VA(*),B(*),X(*)
  integer I,J,L,LA, MAXN,SGNDET!,IFAIL,LP,II
  LOGICAL TRANS,checksol
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),LOGDET!A(:),,R(:),C(:)!,RHS(:),SOL(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:)!,JCN1(:),IRN1(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(6)) then
    LA=INSIZE(6)
  end if
  INSIZE(5)=0
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))!,RHS(M),SOL(M),R(MAXN),C(MAXN),W(5*MAXN)A(LA),
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))!JCN1(LA),IRN1(LA),
!     Factorize matrix
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  ! errors only below debug verbosity (silences duplicate-entry notes)
  if (teems_verbosity()<2) ICNTL(3)=1
  ! -3 workspace shortfalls are handled by caller-side growth, so
  ! MA48's own 'Error return ... because LA is' print would read as
  ! fatal in the logs; the wrappers print the diagnosis on genuinely
  ! fatal INFO(1) codes
  ICNTL(1)=0
  T=M+5*N+4*N/ICNTL(6)+7
  allocate(KEEP(T))

  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).EQ.-3) THEN
    ! workspace too small: hand the suggested size back for a
    ! caller-side reallocate-and-retry
    INSIZE(5)=-3
    INSIZE(6)=max(INFO(3),INFO(4))
    deallocate(CNTL,RINFO,ERROR1,ICNTL,INFO,IW,KEEP)
    RETURN
  END IF
  ! the analyse itself fails on structural singularity (INFO(1)=-5);
  ! previously unchecked, which corrupted the follow-on factorize
  IF (INFO(1).LT.0) THEN
    WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
    CALL TEEMS_ONFAIL_DIAG(INFO(1))
    CALL TEEMS_ONFAIL_ABORT()
  END IF
  if (teems_verbosity()>=2) WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  if (teems_verbosity()>=2) WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
  IF (INFO(1).EQ.-3) THEN
    INSIZE(5)=-3
    INSIZE(6)=INFO(4)
    deallocate(CNTL,RINFO,W,ERROR1,ICNTL,INFO,IW,KEEP)
    RETURN
  END IF
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    WRITE (6,FMT='(A,i10)') 'INFO(5) =',INFO(5)
    CALL TEEMS_ONFAIL_DIAG(INFO(1))
    CALL TEEMS_ONFAIL_ABORT()
  END IF

  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(M,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
  else
    CALL MA48C(M,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,X,ERROR1,W,IW,INFO)
  endif
  if (teems_verbosity()>=2) WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  deallocate(CNTL,RINFO,W,ERROR1)!A,,RHS,SOL
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,
END SUBROUTINE SPEC48M_SSOL2LA

!SUBROUTINE SPEC48_SSOL0(INSIZE,IRN,JCN,VA,B,X)
!  IMPLICIT NONE
!
!  !NEFAC=int(INSIZE(4))
!
!    !write(*,"(A,i5,i5,F20.2)") 'Ain',IRN1(I),JCN1(i),VA(I)
!! Scale matrix
!    !write(*,"(A,i5,i5,F10.2,F10.2,F10.2)") 'Ainbs',I,J,R(I),C(J),A(II)
!    !write(*,"(A,i5,i5,F20.2,F20.2)") 'Ain',I,J,R(I)*C(J),A(II)
!!     Factorize matrix
!  !ICNTL(4)=0
!  !ICNTL(6)=1
!  !CNTL(4)=0.00000000000000001
!  !CNTL(2)=0
!
!  !do I=1,NE
!  !  write(*,"(A,i5,i5,i5,F20.2)") 'AAD',I,IRN1(I),JCN1(I),A(I)
!  !end do
!                RINFO)
!
!    !write(*,"(i5,F10.2)") I,B(I)
!              RHS,SOL,ERROR,W,IW,INFO)
!    !write(*,"(A,i5,F20.2)") 'XD',I,X(I)
!END SUBROUTINE SPEC48_SSOL0
!
!SUBROUTINE SPEC48_SSOL2(INSIZE,IRN,JCN,VA,B,X)
!  IMPLICIT NONE
!
!  !NEFAC=int(INSIZE(4))
!  !write(*,"(i5,i5)") IRN(1),JCN(1)
!  !write(*,"(i5,i5)") IRN(2),JCN(2)
!  !write(*,"(i5,i5)") IRN(3),JCN(3)
!  !write(*,"(i5,i5)") IRN(4),JCN(4)
!!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)
!    !write(*,"(i5,i5,i5,F20.2)") I,IRN1(I),JCN1(I),A(I)
!! Scale matrix
!    !write(*,"(A,i5,i5,F10.2,F10.2,F10.2)") 'Ainbs',I,J,R(I),C(J),A(II)
!
!!     Factorize matrix
!  !CNTL(2)=0.5
!  !CNTL(1)=0.0
!  !JOB=1
!              RHS,SOL,ERROR,W,IW,INFO)
!  !deallocate(A,CNTL,RINFO,W,ERROR)
!END SUBROUTINE SPEC48_SSOL2
!
!SUBROUTINE SPEC48_SSOL1(INSIZE,IRN,JCN,VA,B,X)
!  IMPLICIT NONE
!
!  !write(*,"(i5,i5)") IRN(1),JCN(1)
!  !write(*,"(i5,i5)") IRN(2),JCN(2)
!  !write(*,"(i5,i5)") IRN(3),JCN(3)
!  !write(*,"(i5,i5)") IRN(4),JCN(4)
!
!!  real (kind=myreal) A(LA),CNTL(10),LOGDET,RINFO(10),W(MAXN4)
!  !B(8)=10
!!     Factorize matrix
!  !CNTL(2)=0.5
!  !CNTL(1)=0.0
!  !JOB=1
!              RHS,SOL,ERROR,W,IW,INFO)
!  !deallocate(A,CNTL,RINFO,W,ERROR)
!END SUBROUTINE SPEC48_SSOL1


SUBROUTINE MY_SPAR_ADD3L(VECBIVI1,BIVIINDX1,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER(8) NZ1(*),NZ2(*)
  INTEGER(8) NZ0(*)
  INTEGER(8) BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ2(1))
  INTEGER(8) j,i,l
  j=NZ0(1)
  i=NZ1(1)
  do l=NZ2(1),1,-1
    if(BIVIINDX1(i).LT.BIVIINDX0(j)) then
      VECBIVI1(l)=VECBIVI0(j)
      BIVIINDX1(l)=BIVIINDX0(j)
      j=j-1
    else if(BIVIINDX1(i).EQ.BIVIINDX0(j)) then
      VECBIVI1(l)=VECBIVI1(i)+VECBIVI0(j)
      BIVIINDX1(l)=BIVIINDX0(j)
      j=j-1
      i=i-1
    else
      VECBIVI1(l)=VECBIVI1(i)
      BIVIINDX1(l)=BIVIINDX1(i)
      i=i-1
    end if
  end do
END SUBROUTINE MY_SPAR_ADD3L




SUBROUTINE MY_SPAR_ADD4L(VECBIVI1,BIVIINDX1,IRN,JCN,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2,NCOL1)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER IRN(*),JCN(*),NCOL1(*)
  INTEGER(8) NZ0(*),NZ1(*),NZ2(*)
  INTEGER(8) BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ2(1))
  INTEGER(8) j,i,l,ncol
  ncol=NCOL1(1)
  j=NZ0(1)
  i=NZ1(1)
  do l=NZ2(1),1,-1
    if(BIVIINDX1(i).LT.BIVIINDX0(j)) then
      VECBIVI1(l)=VECBIVI0(j)
          IRN(l)=BIVIINDX0(j)/ncol+1
          JCN(l)=MOD(BIVIINDX0(j),ncol)+1
      j=j-1
    else if(BIVIINDX1(i).EQ.BIVIINDX0(j)) then
      VECBIVI1(l)=VECBIVI1(i)+VECBIVI0(j)
          IRN(l)=BIVIINDX0(j)/ncol+1
          JCN(l)=MOD(BIVIINDX0(j),ncol)+1
      j=j-1
      i=i-1
    else
      VECBIVI1(l)=VECBIVI1(i)
          IRN(l)=BIVIINDX1(i)/ncol+1
          JCN(l)=MOD(BIVIINDX1(i),ncol)+1
      i=i-1
    end if
  end do
END SUBROUTINE MY_SPAR_ADD4L





SUBROUTINE MY_SPAR_COMPL(BIVIINDX1,NZ1,BIVIINDX0,NZ0,NZ2)
  INTEGER(8) BIVIINDX0(*),BIVIINDX1(*),NZ0(*),NZ1(*),NZ2(*)
  INTEGER(8) j,j1,i
  if(BIVIINDX1(NZ1(1)).LT.BIVIINDX0(NZ0(1))) then
      j=1
      j1=0
      do i=1,NZ1(1),1
        do while (BIVIINDX1(i).GT.BIVIINDX0(j))
          j=j+1
        end do
        if(BIVIINDX1(i).LT.BIVIINDX0(j)) then
          j1=j1+1
        end if
      end do
      NZ2(1)=NZ0(1)+j1
  else
      j=1
      j1=0
      do i=1,NZ0(1),1
        do while (BIVIINDX0(i).GT.BIVIINDX1(j))
          j=j+1
        end do
        if(BIVIINDX0(i).LT.BIVIINDX1(j)) then
          j1=j1+1
        end if
      end do
      NZ2(1)=NZ1(1)+j1
  end if
END SUBROUTINE MY_SPAR_COMPL



SUBROUTINE PREP48_ALU1(INSIZE,IRN,JCN,VA,W,IW,KEEP)
  use constants
  IMPLICIT NONE
  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*),INSIZE(*),IW(*),KEEP(*)
  integer M,N,NE,T,RANK,J1,J
  real(kind=DPC) VA(*),W(*)
  real(kind=DPC) LA1
  logical isopen
  !DOUBLE PRECISION LOGDET,SGNDET
  integer LA, MAXN
  real(kind=DPC), pointer :: CNTL(:),RINFO(:)
  integer, pointer :: ICNTL(:),INFO(:)
  character(len=1024) :: filename
  character(len=512) :: scrdir
  integer :: slen
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA1=(NEFAC/100.0)*NE
    LA=ceiling(LA1)
  endif
  if(LA.NE.INSIZE(17)) then
    LA=INSIZE(17)
  end if
  ! INSIZE(18) (out): 0 = factorized, -3 = workspace too small
  ! (INSIZE(19) then holds MA48's suggested LA; caller reallocates,
  ! re-stages IRN/JCN/VA and retries with INSIZE(17) = the new size)
  INSIZE(18)=0
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10))
  allocate(ICNTL(20),INFO(20))!,IW(6*M+3*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  ! errors only below debug verbosity (silences duplicate-entry notes)
  if (teems_verbosity()<2) ICNTL(3)=1
  ! -3 workspace shortfalls are handled by caller-side growth, so
  ! MA48's own 'Error return ... because LA is' print would read as
  ! fatal in the logs; the wrappers print the diagnosis on genuinely
  ! fatal INFO(1) codes
  ICNTL(1)=0
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).EQ.-3) THEN
    INSIZE(18)=-3
    INSIZE(19)=max(INFO(3),INFO(4))
    deallocate(CNTL,RINFO)
    deallocate(ICNTL,INFO)
    RETURN
  END IF
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  CALL TEEMS_ONFAIL_DIAG(INFO(1))
  CALL TEEMS_ONFAIL_ABORT()
  END IF
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).EQ.-3) THEN
    INSIZE(18)=-3
    INSIZE(19)=INFO(4)
    deallocate(CNTL,RINFO)
    deallocate(ICNTL,INFO)
    RETURN
  END IF
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    CALL TEEMS_ONFAIL_DIAG(INFO(1))
    CALL TEEMS_ONFAIL_ABORT()
  END IF
  ! factors (VA/IRN/KEEP, extents LA and T=INSIZE(13)) return to the
  ! caller in place -- C decides between resident handoff and the
  ! legacy _vav/_irnv/_keep scratch files
  deallocate(CNTL,RINFO)!,W)
  deallocate(ICNTL,INFO)!,IW,KEEP)
END SUBROUTINE PREP48_ALU1


SUBROUTINE PREP48M_MSOL(INSIZE,IRN,JCN,VA,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0,JCNB1,SOL,B,W,IW,KEEP)!,IRNV,JCNV,VAV
  use constants
  IMPLICIT NONE
  integer NEFAC,JOB
  integer(8) BIVINZROW0(*),L3
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZCOL0(*)!,JCNV(*),IRNV(*)
  integer(4)IRNB(*),JCNB(*),KEEP(*),IW(*)
  integer(4) JCNC(*),IRNC(*),JCNB1(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1!,NV,MV,NEV
  real(kind=DPC) VA(*),VAC(*),VALUESB(*),VECBIVI(*),SOL(*),B(*),W(*)
  !DOUBLE PRECISION LOGDET,SGNDET!,VAV(*)real (8)
  integer LA, MAXN,NBIVI,MBIVI
  integer(4) I,J,L,L1,L2,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol,isopen
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),ERROR1(:)!,VOUT(:)!,RHS(:)A(:),,SOL(:),B(:),W(:)
  integer, pointer :: ICNTL(:),INFO(:)!,IRNOUT(:),JCNOUT(:)!,JCN1(:),IRN1(:),JCNB1(:),IW(:),KEEP(:)
  character(len=1024) :: filename
  character(len=512) :: scrdir
  integer :: slen
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  MC=INSIZE(4)
  NC=INSIZE(5)
  NEC=INSIZE(6)
  MB=INSIZE(7)
  NB=INSIZE(8)
  NEB=INSIZE(9)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  MBIVI=INSIZE(14)
  NBIVI=INSIZE(15)
  JCN(1:NE)=JCN(1:NE)+1
  JCNB1(1:NEB)=JCNB(1:NEB)+1
  do J=1,M
    L=J+NE
    IRN((IRN(L)+1):IRN(L+1))=J
  end do
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(17)) then
    LA=INSIZE(17)
  end if
  ! INSIZE(18) (out): 0 = factorized, -3 = workspace too small
  ! (INSIZE(19) then holds MA48's suggested LA; caller reallocates,
  ! re-stages IRN/JCN/VA -- clobbered by the conversion above and by
  ! MA48 -- and retries with INSIZE(17) = the new size)
  INSIZE(18)=0
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))!,RHS(M),W(5*MAXN)A(LA),,SOL(M),B(M)
  allocate(ICNTL(20),INFO(20))!JCN1(LA),IRN1(LA),,IW(6*M+3*N)
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  ! errors only below debug verbosity (silences duplicate-entry notes)
  if (teems_verbosity()<2) ICNTL(3)=1
  ! -3 workspace shortfalls are handled by caller-side growth, so
  ! MA48's own 'Error return ... because LA is' print would read as
  ! fatal in the logs; the wrappers print the diagnosis on genuinely
  ! fatal INFO(1) codes
  ICNTL(1)=0
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).EQ.-3) THEN
    INSIZE(18)=-3
    INSIZE(19)=max(INFO(3),INFO(4))
    deallocate(CNTL,RINFO,ERROR1)
    deallocate(ICNTL,INFO)
    RETURN
  END IF
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  CALL TEEMS_ONFAIL_DIAG(INFO(1))
  CALL TEEMS_ONFAIL_ABORT()
  END IF
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).EQ.-3) THEN
    INSIZE(18)=-3
    INSIZE(19)=INFO(4)
    deallocate(CNTL,RINFO,ERROR1)
    deallocate(ICNTL,INFO)
    RETURN
  END IF
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    CALL TEEMS_ONFAIL_DIAG(INFO(1))
    CALL TEEMS_ONFAIL_ABORT()
  END IF
  JOB=1
  L=0
  TRANS = .FALSE.
  do J=1,MC-1!NC
      J3=IRNC(J+1)-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3!NEC
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      endif
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  end do
  J=MC
      J3=NEC-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3!NEC
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      endif
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              DO j2=IRNB(I)+1,IRNB(I+1)
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then !GE.1E-10.OR.SOL(I).LE.-1E-10) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)!(JCNB(j2)+1)!+NBIVI*(BIVINZROW0(JCNB(j2)+1)-1)
              do j2=IRNB(MB)+1,NEB
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  ! factors (VA/IRN/KEEP, extents LA and T=INSIZE(13)) return to the
  ! caller in place -- C decides between resident handoff and the
  ! legacy _vav/_irnv/_keep scratch files
  deallocate(CNTL,RINFO,ERROR1)!,RHSA,,SOL,JCNB1,B,W
  deallocate(ICNTL,INFO)!IRN1,,JCNOUT,IW,KEEP
END SUBROUTINE PREP48M_MSOL

SUBROUTINE PREP48M_MSOL_P(INSIZE,IRN,JCN,VA,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0,JCNB1,SOL,B,W,IW,&
  KEEP,REDO)
  ! Persistent-factor variant of PREP48M_MSOL for the NDBBD regional
  ! blocks (-fastrefac).  REDO(1) on entry: 0 = full analyse+factorize
  ! (as PREP48M_MSOL); 1 = VA(1:NE) refilled by the caller in CSR
  ! order, IRN/JCN/KEEP are the persisted MA48 state, factorize with
  ! MA48B/BD JOB=2.  On exit: 0 = factorized and border solves done,
  ! <0 = fast factorize declined (MA48B/BD INFO(1)) BEFORE any border
  ! work — caller re-stages the block and retries with 0.
  use constants
  IMPLICIT NONE
  integer NEFAC,JOB
  integer(8) BIVINZROW0(*),L3
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZCOL0(*),REDO(*)
  integer(4)IRNB(*),JCNB(*),KEEP(*),IW(*)
  integer(4) JCNC(*),IRNC(*),JCNB1(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1
  real(kind=DPC) VA(*),VAC(*),VALUESB(*),VECBIVI(*),SOL(*),B(*),W(*)
  integer LA, MAXN,NBIVI,MBIVI
  integer(4) I,J,L,L1,L2,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),ERROR1(:)
  integer, pointer :: ICNTL(:),INFO(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  MC=INSIZE(4)
  NC=INSIZE(5)
  NEC=INSIZE(6)
  MB=INSIZE(7)
  NB=INSIZE(8)
  NEB=INSIZE(9)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  MBIVI=INSIZE(14)
  NBIVI=INSIZE(15)
  JCNB1(1:NEB)=JCNB(1:NEB)+1
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  if(LA.NE.INSIZE(17)) then
    LA=INSIZE(17)
  end if
  ! INSIZE(18) (out): 0 = factorized, -3 = workspace too small
  ! (INSIZE(19) then holds MA48's suggested LA; caller reallocates the
  ! persisted arrays, re-stages and retries with REDO(1)=0 and
  ! INSIZE(17) = the new size); REDO(1) is also set to -3 so the
  ! caller's decline channel sees the failure
  INSIZE(18)=0
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3))
  allocate(ICNTL(20),INFO(20))
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  ! errors only below debug verbosity (silences duplicate-entry notes)
  if (teems_verbosity()<2) ICNTL(3)=1
  ! -3 workspace shortfalls are handled by caller-side growth, so
  ! MA48's own 'Error return ... because LA is' print would read as
  ! fatal in the logs; the wrappers print the diagnosis on genuinely
  ! fatal INFO(1) codes
  ICNTL(1)=0
  ! blocks whose entries turn unsuitable for the kept pivot sequence
  ! are refactorized as on a JOB=1 call
  ICNTL(11)=1
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  IF (REDO(1).EQ.0) THEN
    JCN(1:NE)=JCN(1:NE)+1
    do J=1,M
      L=J+NE
      IRN((IRN(L)+1):IRN(L+1))=J
    end do
    IF (FSORD.EQ.1) THEN
      CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
    else
      CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
    endif
    IF (INFO(1).EQ.-3) THEN
      INSIZE(18)=-3
      INSIZE(19)=max(INFO(3),INFO(4))
      REDO(1)=-3
      deallocate(CNTL,RINFO,ERROR1)
      deallocate(ICNTL,INFO)
      RETURN
    END IF
    IF (INFO(1).LT.0) THEN
    WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
    CALL TEEMS_ONFAIL_DIAG(INFO(1))
    CALL TEEMS_ONFAIL_ABORT()
    END IF
    IF (FSORD.EQ.1) THEN
      CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                  RINFO)
    else
      CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                  RINFO)
    endif
    IF (INFO(1).EQ.-3) THEN
      INSIZE(18)=-3
      INSIZE(19)=INFO(4)
      REDO(1)=-3
      deallocate(CNTL,RINFO,ERROR1)
      deallocate(ICNTL,INFO)
      RETURN
    END IF
    IF (INFO(1).NE.0) THEN
      WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
      INFO(1),'Solution not possible'
      write(*,"(A,i5)") 'RANK',INFO(5)
      CALL TEEMS_ONFAIL_DIAG(INFO(1))
      CALL TEEMS_ONFAIL_ABORT()
    END IF
  ELSE
    IF (FSORD.EQ.1) THEN
      CALL MA48BD(M,N,NE,2,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                  RINFO)
    else
      CALL MA48B(M,N,NE,2,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                  RINFO)
    endif
    IF (INFO(1).EQ.-3) THEN
      ! workspace too small for the kept factors: report -3 so the
      ! caller grows the arrays and redoes the analyse
      INSIZE(18)=-3
      INSIZE(19)=INFO(4)
      REDO(1)=-3
      deallocate(CNTL,RINFO,ERROR1)
      deallocate(ICNTL,INFO)
      RETURN
    END IF
    IF (INFO(1).LT.0) THEN
      if (teems_verbosity()>=1) WRITE (6,'(A,I3)') &
        'Note: fast block refactorize declined, MA48B/BD INFO(1) =',INFO(1)
      REDO(1)=INFO(1)
      deallocate(CNTL,RINFO,ERROR1)
      deallocate(ICNTL,INFO)
      RETURN
    END IF
  END IF
  JOB=1
  L=0
  TRANS = .FALSE.
  do J=1,MC-1
      J3=IRNC(J+1)-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      endif
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)
              DO j2=IRNB(I)+1,IRNB(I+1)
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)
              do j2=IRNB(MB)+1,NEB
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  end do
  J=MC
      J3=NEC-IRNC(J)
    if(J3.GT.0) then
      B(1:M)=0
      do I=1,J3
        L2=I+L
        B(JCNC(L2)+1)=VAC(L2)
      end do
      L=L+J3
      IF (FSORD.EQ.1) THEN
        CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      else
        CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
                  B,SOL,ERROR1,W,IW,INFO)
      endif
      do I=1,MB-1
        if(IRNB(I).NE.IRNB(I+1)) then
              L3=BIVINZCOL0(J)+BIVINZROW0(I)
              DO j2=IRNB(I)+1,IRNB(I+1)
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
      end do
        if(IRNB(MB).NE.NEB) then
              L3=BIVINZCOL0(J)+BIVINZROW0(MB)
              do j2=IRNB(MB)+1,NEB
                VECBIVI(L3)=VECBIVI(L3)-SOL(JCNB1(j2))*VALUESB(j2)
              end do
        end if
    end if
  REDO(1)=0
  ! factors (VA/IRN/JCN/KEEP) stay in the caller's persistent slot
  deallocate(CNTL,RINFO,ERROR1)
  deallocate(ICNTL,INFO)
END SUBROUTINE PREP48M_MSOL_P

SUBROUTINE SPAR_MULMIN(SOL,NROW,NZ,IRN,JCN,VA,RES)
  use constants
  IMPLICIT NONE
  integer I,J
  integer(4) JCN(*),IRN(*),NROW(*),NZ(*)
  real(kind=DPC) VA(*),RES(*),SOL(0:NROW(1)-1)
      do I=1,NROW(1)-1
              DO J=IRN(I)+1,IRN(I+1)
                !$OMP ATOMIC
                RES(I)=RES(I)-VA(j)*SOL(JCN(J))
              end do
      end do
              DO J=IRN(NROW(1))+1,NZ(1)
                !$OMP ATOMIC
                RES(I)=RES(I)-VA(j)*SOL(JCN(J))
              end do
END SUBROUTINE SPAR_MULMIN


SUBROUTINE SPAR_MULNOADD(SOL,NROW,NZ,IRN,JCN,VA,RES)
  use constants
  IMPLICIT NONE
  integer I,J
  integer(4) JCN(*),IRN(*),NROW(*),NZ(*)
  real(kind=DPC) VA(*),RES(*),SOL(0:NROW(1)-1)
  RES(1:NROW(1))=0
      do I=1,NROW(1)-1
              DO J=IRN(I)+1,IRN(I+1)
                RES(I)=RES(I)+VA(j)*SOL(JCN(J))
              end do
      end do
              DO J=IRN(NROW(1))+1,NZ(1)
                RES(NROW(1))=RES(NROW(1))+VA(j)*SOL(JCN(J))
              end do
END SUBROUTINE SPAR_MULNOADD

SUBROUTINE SPAR_VBIVIADD(SOL,BVCOL,BVROW,BVSIZE,NROW,NCOL,NZ,IRN,JCN,VA,RES)
  use constants
  IMPLICIT NONE
  integer I,J
  integer(4) BVCOL(*),JCN(*),IRN(*),NROW(*),NZ(*),NCOL(*)
  integer(8) BVROW(*),BVSIZE(*),L
  real(kind=DPC) VA(*),RES(0:(BVSIZE(1)-1)),SOL(0:NCOL(1))
      do I=1,NROW(1)-1
        L=BVCOL(1)+BVROW(I)
              DO J=IRN(I)+1,IRN(I+1)
                !$OMP ATOMIC
                RES(L)=RES(L)-VA(J)*SOL(JCN(J))
              end do
      end do
        L=BVCOL(1)+BVROW(NROW(1))
              DO J=IRN(NROW(1))+1,NZ(1)
                !$OMP ATOMIC
                RES(L)=RES(L)-VA(J)*SOL(JCN(J))
              end do
END SUBROUTINE SPAR_VBIVIADD

SUBROUTINE PATIO_MAT(insizeda,IRN,JCN,VBIVI,IRN1A,JCN1A)
  use constants
  IMPLICIT NONE
  integer I,J,colcut,rowcut
  integer(4) JCN(*),IRN(*),IRN1A(*),insizeda(*),JCN1A(*)
  real(kind=DPC) VBIVI(*)
  LOGICAL iterstop
        do I=1,insizeda(3)
          iterstop=.TRUE.
          colcut=0
          rowcut=0
          if(insizeda(1).NE.insizeda(4)) then
          do J=insizeda(1),insizeda(4)+1,-1
            if(irn1a(J).LT.irn(I)) then
              rowcut=rowcut+1
            else if(irn1a(j).EQ.irn(i)) then
              JCN(i)=1
              IRN(i)=1
              VBIVI(i)=0
              iterstop=.false.
              exit
            endif
          end do
          endif
          if(iterstop.EQV..true.) then
            if(insizeda(2).NE.insizeda(4)) then
            do J=insizeda(2),insizeda(4)+1,-1
            if(jcn1a(J).LT.jcn(I)) then
              colcut=colcut+1
            else if(JCN1A(j).EQ.JCN(i)) then
              JCN(i)=1
              IRN(i)=1
              VBIVI(i)=0
              iterstop=.false.
              exit
            endif
            end do
            endif
          endif
          if(iterstop.EQV..true.) then
            IRN(i)=IRN(i)-rowcut
            JCN(i)=JCN(i)-colcut
          endif
        end do
END SUBROUTINE PATIO_MAT
