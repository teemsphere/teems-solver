module constants
  implicit none
  INTEGER, PARAMETER :: DPC = KIND(0.0D0)
  INTEGER, PARAMETER :: FSORD = 1
  !switch to single by setting FSORD = 0, DPC and hsl_mp48d.f90 accordingly
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
    write(*,"(i10,i10,i10,i10)") nz,m,n,nblocks
    call mc66(m,n,nz,irn1,jcn1,nblocks,control,seed, &
      row_order,info,rowptr,column_order,colptr,&
      netcut,rowdiff,kblocks)
    write(*,"(a)") "generating the ordering =====1 "
    if (info /= 0) then
      call mc66_print_message(info)
      if (info < 0) stop "mc66 failed"
    end if
    write(*,"(a,I10)")       "netcut =                  ",netcut
    write(*,"(a,f10.2,'%')") "row dimension imbalance = ", rowdiff
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
      write(*,"('block ',i4,' of dimension  ',i10,' X ',i10)") &
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
  print *, "Process ", myid, " of ", numprocs, " is alive1"
  data%JOB = 25
  CALL MP48AD(data)
  IF (data%RANK.EQ.0) THEN
    IF (data%ERROR.LT.0) THEN
      WRITE (6,*) ' Unexpected error return'
    ELSE
      do o=1,m
        x(column_order(o))=data%X(o)
      end do
    END IF
    deallocate(row_order,rowptr,column_order, colptr)
  END IF
  call MPI_COMM_RANK( fcomm1, myid, ERCODE ) !MPI_COMM_WORLD
  call MPI_COMM_SIZE( fcomm1, numprocs, ERCODE ) !MPI_COMM_WORLD
  print *, "Process ", myid, " of ", numprocs, " is alive2"
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
      write(*,"('nomc block ',i4,' of dimension  ',i10,' X ',i10)") &
      i,data%NEQSB(i),colptrin(i+1)-colptrin(i)
      ncols=colptrin(i+1)-colptrin(i)
      if(maxsbcols.LT.ncols)maxsbcols=ncols
    end do
    maxsbcols=maxsbcols+m-colptrin(nblocks+1)
    print *, "row ", m,"nz ",nz,"maxcolsb ",maxsbcols
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
  print *, "Process ", myid, " of ", numprocs, " is alive1"
  data%JOB = 25
  CALL MP48AD(data)
  print *, "Processa ", myid, " of ", numprocs, " is alive1"
  IF (data%RANK.EQ.0) THEN
    IF (data%ERROR.LT.0) THEN
      WRITE (6,*) ' Unexpected error return'
    ELSE
      do o=1,m
        x(o)=data%X(o)
      end do
    END IF
  END IF
  print *, "Process ", myid, " of ", numprocs, " is alive2"
  data%JOB = 6
  CALL MP48AD(data)
END SUBROUTINE SPEC48_NOMC66

SUBROUTINE SPEC51_RANK(INSIZE,CNTL6,IRN,JCN,VA)
  use constants
  IMPLICIT NONE

  integer M,N,NE
  integer (4) INSIZE(*)
  integer (4) JCN(*),IRN(*)
  real (kind=DPC) VA(*),CNTL6(*)
  integer I,LA, MAXN,RANK1,SGNDET,T,NEFAC
  real (kind=DPC) LOGDET
  real (kind=DPC), pointer :: CNTL(:),RINFO(:),W(:)!A(:),
  integer, pointer :: COLS(:),ICNTL(:),INFO(:),IW(:),KEEP(:),ROWS(:)!JCN1(:),IRN1(:),
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(5)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),W(5*MAXN))!A(LA),


  allocate(COLS(N),ICNTL(20),INFO(20),IW(6*M+3*N),ROWS(M))!JCN1(LA),IRN1(LA),

!     Factorize matrix
    IF (FSORD.EQ.1) THEN
      CALL MA48ID(CNTL,ICNTL)
    ELSE
      CALL MA48I(CNTL,ICNTL)
    ENDIF
    IF(CNTL6(1).EQ.0)THEN
      IF (FSORD.EQ.1) THEN
        CNTL(4)=1e-4
      else
        CNTL(4)=0.3
      ENDIF
    ELSE
      CNTL(4)=CNTL6(1)!1e-4!0.0000000001
    ENDIF
    T=M+5*N+4*N/ICNTL(6)+7
    allocate(KEEP(T))
    ICNTL(7)=0
    IF (FSORD.EQ.1) THEN
      CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)

!     Compute the determinant
      CALL MA51CD(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ELSE
      CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
      CALL MA51C(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ENDIF

    IF(SGNDET.GT.0)THEN
      WRITE(6,'(/,A,F9.3)')&
          ' Determinant is positive; log(determinant) =',LOGDET
          INSIZE(4)=M
    ELSE IF(SGNDET.LT.0)THEN
        WRITE(6,'(/,A,F9.3)')&
          ' Determinant is negative; log(-determinant) =',LOGDET
          INSIZE(4)=M
      ELSE
        WRITE(6,'(/,A)')' Determinant is zero'
!     Determine the nonsingular submatrix of the factorization
        IF (FSORD.EQ.1) THEN
          CALL MA51AD(M,N,LA,IRN,KEEP,RANK1,ROWS,COLS,W)
        ELSE
          CALL MA51A(M,N,LA,IRN,KEEP,RANK1,ROWS,COLS,W)
        ENDIF
        INSIZE(4)=RANK1
        do I=1,M
          IRN(I)=ROWS(I)
        end do
        do I=1,N
          JCN(I)=COLS(I)
        end do
    END IF
  deallocate(CNTL,RINFO,W)!A,
  deallocate(COLS,ICNTL,INFO,IW,KEEP,ROWS)!IRN1,JCN1,
END SUBROUTINE SPEC51_RANK

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
    IF (FSORD.EQ.1) THEN
      CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)

!     Compute the determinant
      CALL MA51CD(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ELSE
      CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
      CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
      CALL MA51C(M,N,LA,VA,IRN,KEEP,SGNDET,LOGDET,W)
    ENDIF

    IF(SGNDET.GT.0)THEN
      WRITE(6,'(/,A,F9.3)')&
          ' Determinant is positive; log(determinant) =',LOGDET
          INSIZE(4)=M
    ELSE IF(SGNDET.LT.0)THEN
        WRITE(6,'(/,A,F9.3)')&
          ' Determinant is negative; log(-determinant) =',LOGDET
          INSIZE(4)=M
      ELSE
        WRITE(6,'(/,A)')' Determinant is zero'
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

SUBROUTINE SPEC48_MSOL(INSIZE,IRN,JCN,VA,B,X,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0)!,IRNV,JCNV,VAV
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZROW0(*),BIVINZCOL0(*)!,JCNV(*),IRNV(*)
  integer(4)IRNB(*),JCNB(*)
  integer(4) JCNC(*),IRNC(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1!,NV,MV,NEV
  real(kind=DPC) VA(*),VAC(*),B(*),X(*),VALUESB(*),VECBIVI(*)
  !DOUBLE PRECISION LOGDET,SGNDET!,VAV(*)real (8)
  integer LA, MAXN,NBIVI,MBIVI
  integer(4) I,J,L,L1,L2,L3,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol,isopen
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),SOL(:)!,VOUT(:)!,RHS(:)A(:),
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCNB1(:)!,IRNOUT(:),JCNOUT(:)!,JCN1(:),IRN1(:)
  character(len=1024) :: filename
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
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  allocate(KEEP(T))
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ELSE
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
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
  J=INSIZE(16)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)!,RHSA,
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,,JCNOUT
END SUBROUTINE SPEC48_MSOL

SUBROUTINE SPEC48M_MSOL(INSIZE,IRN,JCN,VA,B,X,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0)!,IRNV,JCNV,VAV
  use constants
  IMPLICIT NONE

  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZROW0(*),BIVINZCOL0(*)!,JCNV(*),IRNV(*)
  integer(4)IRNB(*),JCNB(*)
  integer(4) JCNC(*),IRNC(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1!,NV,MV,NEV
  real(kind=DPC) VA(*),VAC(*),B(*),X(*),VALUESB(*),VECBIVI(*)
  !DOUBLE PRECISION LOGDET,SGNDET!,VAV(*)real (8)
  integer LA, MAXN,NBIVI,MBIVI
  integer(4) I,J,L,L1,L2,L3,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol,isopen
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),SOL(:)!,VOUT(:)!,RHS(:)A(:),
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCNB1(:)!,IRNOUT(:),JCNOUT(:)!,JCN1(:),IRN1(:)
  character(len=1024) :: filename
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
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  allocate(KEEP(T))
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ELSE
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,RINFO)
  ENDIF
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
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
  J=INSIZE(16)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1)!,RHSA,
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,,JCNOUT
END SUBROUTINE SPEC48M_MSOL

SUBROUTINE SPEC48_ESOL(INSIZE,IRN,VA,KEEP,B,SOL)
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
END SUBROUTINE SPEC48_ESOL

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

SUBROUTINE SPEC48_RPESOL(INSIZE,IRN,VA,KEEP,B,SOL,CNTL,RINFO,ERROR1,ICNTL,INFO,W,IW)
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

  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  ENDIF
  JOB=1
  TRANS = .FALSE.
  IF (FSORD.EQ.1) THEN
    CALL MA48CD(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  else
    CALL MA48C(N,N,TRANS,JOB,LA,VA,IRN,KEEP,CNTL,ICNTL,&
              B,SOL,ERROR1,W,IW,INFO)
  ENDIF
END SUBROUTINE SPEC48_RPESOL

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


SUBROUTINE SPEC48_SSOL(INSIZE,IRN,JCN,VA,B,X)
  use constants
  IMPLICIT NONE

  INTEGER, PARAMETER :: myreal = SELECTED_REAL_KIND(8)!kind(real*8)!KIND( 0.0D0 )
  INTEGER, PARAMETER :: myint = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  INTEGER, PARAMETER :: myint32 = selected_int_kind(16)!kind(interger*4)!KIND( 1)
  integer NEFAC,JOB
  integer (kind=myint) JCN(*),IRN(*)
  integer (kind=myint32) INSIZE(*)
  integer M,N,NE,T
  real (kind=DPC) VA(*),B(*),X(*)
  integer I,J,L,LA, MAXN,SGNDET!,IFAIL,LP,II
  LOGICAL TRANS,checksol
  real (kind=DPC), pointer :: A(:),CNTL(:),RINFO(:),W(:),ERROR(:)!,LOGDET!,R(:),C(:)!,RHS(:),SOL(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCN1(:),IRN1(:)
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(4)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(A(LA),CNTL(10),RINFO(10),ERROR(3))!,RHS(M),SOL(M),R(MAXN),C(MAXN),W(5*MAXN)
  allocate(JCN1(LA),IRN1(LA),ICNTL(20),INFO(20),IW(6*M+3*N))

  do I=1,NE
    IRN1(I)=IRN(I)
    JCN1(I)=JCN(I)
    A(I)=VA(I)
  end do
  write(*,"(A,i10,i10,i15)") 'Matsize ',M,N,NE
!     Factorize matrix
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  T=M+5*N+4*N/ICNTL(6)+7
  allocate(KEEP(T))

  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
  STOP
  END IF
  JOB=1
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,A,IRN1,JCN1,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    STOP
  END IF

  deallocate(JCN1)
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
    CALL MA48CD(M,N,TRANS,JOB,LA,A,IRN1,KEEP,CNTL,ICNTL,&
              B,X,ERROR,W,IW,INFO)
  else
    CALL MA48C(M,N,TRANS,JOB,LA,A,IRN1,KEEP,CNTL,ICNTL,&
              B,X,ERROR,W,IW,INFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  deallocate(A,CNTL,RINFO,W,ERROR)!,RHS,SOL
  deallocate(IRN1,ICNTL,INFO,IW,KEEP)
END SUBROUTINE SPEC48_SSOL

SUBROUTINE SPEC48_SSOL2LA(INSIZE,IRN,JCN,VA,B,X)
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
  T=M+5*N+4*N/ICNTL(6)+7
  allocate(KEEP(T))

  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
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
  WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    WRITE (6,FMT='(A,i10)') 'INFO(5) =',INFO(5)
    STOP
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
  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
  deallocate(CNTL,RINFO,W,ERROR1)!A,,RHS,SOL
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,
END SUBROUTINE SPEC48_SSOL2LA

SUBROUTINE SPEC48M_SSOL2LA(INSIZE,IRN,JCN,VA,B,X)
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
  T=M+5*N+4*N/ICNTL(6)+7
  allocate(KEEP(T))

  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,JOB,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  WRITE (6,FMT='(A,I3/A)') 'INFO(3) =',INFO(3)/NE
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
  WRITE (6,FMT='(A,I3/A)') 'INFO(4) =',INFO(4)/NE
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    WRITE (6,FMT='(A,i10)') 'INFO(5) =',INFO(5)
    STOP
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
  WRITE (6,FMT='(A,I3/A)') 'INFO(1) =',INFO(1)
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

SUBROUTINE MY_SPAR_ADD(VECBIVI1,BIVIINDX1,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER NZ1(*),NZ0(*),NZ2(*)
  INTEGER BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  INTEGER j,j1,j1p,i,i1,j2,j3,j2p,ip,b1,e1,b2,e2,b3,e3,buffer,m,l
  buffer=100000!to avoid crash in fotran assignment
  j=NZ1(1)
  j1=NZ0(1)
  i=NZ2(1)
  do while(i.GT.0)!i=NZ2(1),1,-1
    i1=i
    j2=j
    j3=j1
    do while(BIVIINDX1(j).EQ.BIVIINDX0(j1).AND.BIVIINDX1(j).GT.-1)
          j1=j1-1
          j=j-1
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    j2p=j1+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)+VECBIVI0(j2p:j3)
        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          b3=j2p+buffer*m
          e3=b3+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)+VECBIVI0(b3:e3)
          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)+VECBIVI0(e3+1:j3)
          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)+VECBIVI0(j3-j2+ip:j3)
        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        e3=j3-j2+ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)+VECBIVI0(e3:j3)
        BIVIINDX1(e1:i1)=BIVIINDX1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)+VECBIVI0(j2p:j3-j2+i)
        BIVIINDX1(ip:i1-j2+i)=BIVIINDX1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        e3=j2p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)+VECBIVI0(e3:j3-j2+i)
        BIVIINDX1(e1:i1-j2+i)=BIVIINDX1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      end if
    end if
    i1=i
    j2=j
    do while(BIVIINDX1(j).GT.BIVIINDX0(j1))
          j=j-1
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)
        BIVIINDX1(e1:i1)=BIVIINDX1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
        BIVIINDX1(ip:ip+i-j1p)=BIVIINDX1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)
        BIVIINDX1(e1:ip+i-j1p)=BIVIINDX1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if

!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
    end if
    end if
    i1=i
    j2=j1
    do while (BIVIINDX1(j).LT.BIVIINDX0(j1))
    j1=j1-1
    end do
    i=i1-(j2-j1)
    if(i.NE.i1) then
    ip=i+1
    j1p=j1+1
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI0(j1p:j2)
        BIVIINDX1(ip:i1)=BIVIINDX0(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI0(b2:e2)
          BIVIINDX1(b1:e1)=BIVIINDX0(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI0(e2+1:j2)
          BIVIINDX1(e1+1:i1)=BIVIINDX0(e2+1:j2)
      end if
    end if
    !20 CONTINUE
  end do
END SUBROUTINE MY_SPAR_ADD

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

SUBROUTINE MY_SPAR_ADDL(VECBIVI1,BIVIINDX1,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER NZ1(*),NZ0(*),NZ2(*)
  INTEGER(8) BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  INTEGER j,j1,j1p,i,i1,j2,j3,j2p,ip,b1,e1,b2,e2,b3,e3,buffer,m,l
  buffer=100000!to avoid crash in fotran assignment
  j=NZ1(1)
  j1=NZ0(1)
  i=NZ2(1)
  do while(i.GT.0)!i=NZ2(1),1,-1
    i1=i
    j2=j
    j3=j1
    do while(BIVIINDX1(j).EQ.BIVIINDX0(j1).AND.BIVIINDX1(j).GT.-1)
          j1=j1-1
          j=j-1
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    j2p=j1+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)+VECBIVI0(j2p:j3)
        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          b3=j2p+buffer*m
          e3=b3+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)+VECBIVI0(b3:e3)
          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)+VECBIVI0(e3+1:j3)
          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)+VECBIVI0(j3-j2+ip:j3)
        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        e3=j3-j2+ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)+VECBIVI0(e3:j3)
        BIVIINDX1(e1:i1)=BIVIINDX1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)+VECBIVI0(j2p:j3-j2+i)
        BIVIINDX1(ip:i1-j2+i)=BIVIINDX1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        e3=j2p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)+VECBIVI0(e3:j3-j2+i)
        BIVIINDX1(e1:i1-j2+i)=BIVIINDX1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      end if
    end if
    i1=i
    j2=j
    do while(BIVIINDX1(j).GT.BIVIINDX0(j1))
          j=j-1
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
        BIVIINDX1(ip:i1)=BIVIINDX1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
          BIVIINDX1(b1:e1)=BIVIINDX1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
          BIVIINDX1(e1+1:i1)=BIVIINDX1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
        BIVIINDX1(i1-j2+ip:i1)=BIVIINDX1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)
        BIVIINDX1(e1:i1)=BIVIINDX1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
        BIVIINDX1(ip:ip+i-j1p)=BIVIINDX1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)
        BIVIINDX1(e1:ip+i-j1p)=BIVIINDX1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          BIVIINDX1(b1:e1-1)=BIVIINDX1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if

!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
!        !WRITE (*,"(i10,i10,i10,i10,i10)"),l,ip,i1,j1p,j2
!          !WRITE (*,"(i10,i10,i10,i10,i10)"),m,b1,e1,b2,e2
    end if
    end if
    i1=i
    j2=j1
    do while (BIVIINDX1(j).LT.BIVIINDX0(j1))
    j1=j1-1
    end do
    i=i1-(j2-j1)
    if(i.NE.i1) then
    ip=i+1
    j1p=j1+1
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI0(j1p:j2)
        BIVIINDX1(ip:i1)=BIVIINDX0(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI0(b2:e2)
          BIVIINDX1(b1:e1)=BIVIINDX0(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI0(e2+1:j2)
          BIVIINDX1(e1+1:i1)=BIVIINDX0(e2+1:j2)
      end if
    end if
    !20 CONTINUE
  end do
END SUBROUTINE MY_SPAR_ADDL

SUBROUTINE MY_SPAR_ADD1(VECBIVI1,BIVIINDX1,IRN,JCN,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2,NCOL1)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER IRN(*),JCN(*),NZ1(*),NZ0(*),NZ2(*),NCOL1(*)
  INTEGER BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  INTEGER j,j1,j1p,j2p,i,ncol,i1,j2,j3,ip,b1,e1,b2,b3,e3,e2,buffer,m,l
  buffer=100000!to avoid crash in fotran assignment
  j=NZ1(1)
  j1=NZ0(1)
  ncol=NCOL1(1)
  i=NZ2(1)
  do while(i.GT.0)!i=NZ2(1),1,-1
    i1=i
    j2=j
    j3=j1
    do while(BIVIINDX1(j).EQ.BIVIINDX0(j1).AND.BIVIINDX1(j).GT.-1)
          j1=j1-1
          j=j-1
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    j2p=j1+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)+VECBIVI0(j2p:j3)
        IRN(ip:i1)=BIVIINDX1(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX1(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          b3=j2p+buffer*m
          e3=b3+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)+VECBIVI0(b3:e3)
          IRN(b1:e1)=BIVIINDX1(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX1(b2:e2),ncol)+1
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)+VECBIVI0(e3+1:j3)
          IRN(e1+1:i1)=BIVIINDX1(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX1(e2+1:j2),ncol)+1
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)+VECBIVI0(j3-j2+ip:j3)
        IRN(i1-j2+ip:i1)=BIVIINDX1(ip:j2)/ncol+1
        JCN(i1-j2+ip:i1)=MOD(BIVIINDX1(ip:j2),ncol)+1
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        e3=j3-j2+ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)+VECBIVI0(e3:j3)
        IRN(e1:i1)=BIVIINDX1(e2:j2)/ncol+1
        JCN(e1:i1)=MOD(BIVIINDX1(e2:j2),ncol)+1
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          IRN(b1:e1-1)=BIVIINDX1(b2:e2-1)/ncol+1
          JCN(b1:e1-1)=MOD(BIVIINDX1(b2:e2-1),ncol)+1
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)+VECBIVI0(j2p:j3-j2+i)
        IRN(ip:ip+i-j1p)=BIVIINDX1(j1p:i)/ncol+1
        JCN(ip:ip+i-j1p)=MOD(BIVIINDX1(j1p:i),ncol)+1
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        e3=j2p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)+VECBIVI0(e3:j3-j2+i)
        IRN(e1:ip+i-j1p)=BIVIINDX1(e2:i)/ncol+1
        JCN(e1:ip+i-j1p)=MOD(BIVIINDX1(e2:i),ncol)+1
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          IRN(b1:e1-1)=BIVIINDX1(b2:e2-1)/ncol+1
          JCN(b1:e1-1)=MOD(BIVIINDX1(b2:e2-1),ncol)+1
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      end if
    end if
    i1=i
    j2=j
    do while(BIVIINDX1(j).GT.BIVIINDX0(j1))
          j=j-1
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1

    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
    end if

      if(i1-ip.LT.buffer) then
        IRN(ip:i1)=BIVIINDX1(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX1(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          IRN(b1:e1)=BIVIINDX1(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX1(b2:e2),ncol)+1
        end do
          IRN(e1+1:i1)=BIVIINDX1(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX1(e2+1:j2),ncol)+1
      end if
    end if
    i1=i
    j2=j1
    do while(BIVIINDX1(j).LT.BIVIINDX0(j1))
    j1=j1-1
    end do
    i=i1-(j2-j1)
    ip=i+1
    j1p=j1+1
    if(i.NE.i1) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI0(j1p:j2)
        IRN(ip:i1)=BIVIINDX0(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX0(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI0(b2:e2)
          IRN(b1:e1)=BIVIINDX0(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX0(b2:e2),ncol)+1
        end do
          VECBIVI1(e1+1:i1)=VECBIVI0(e2+1:j2)
          IRN(e1+1:i1)=BIVIINDX0(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX0(e2+1:j2),ncol)+1
      end if
    end if
    !20 CONTINUE
  end do
END SUBROUTINE MY_SPAR_ADD1

SUBROUTINE MY_SPAR_ADD1L(VECBIVI1,BIVIINDX1,IRN,JCN,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2,NCOL1)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*)
  INTEGER IRN(*),JCN(*),NZ1(*),NZ0(*),NZ2(*),NCOL1(*)
  INTEGER(8) BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  INTEGER j,j1,j1p,j2p,i,ncol,i1,j2,j3,ip,b1,e1,b2,b3,e3,e2,buffer,m,l
  buffer=100000!to avoid crash in fotran assignment
  j=NZ1(1)
  j1=NZ0(1)
  ncol=NCOL1(1)
  i=NZ2(1)
  do while(i.GT.0)!i=NZ2(1),1,-1
    i1=i
    j2=j
    j3=j1
    do while(BIVIINDX1(j).EQ.BIVIINDX0(j1).AND.BIVIINDX1(j).GT.-1)
          j1=j1-1
          j=j-1
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1
    j2p=j1+1
    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)+VECBIVI0(j2p:j3)
        IRN(ip:i1)=BIVIINDX1(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX1(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          b3=j2p+buffer*m
          e3=b3+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)+VECBIVI0(b3:e3)
          IRN(b1:e1)=BIVIINDX1(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX1(b2:e2),ncol)+1
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)+VECBIVI0(e3+1:j3)
          IRN(e1+1:i1)=BIVIINDX1(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX1(e2+1:j2),ncol)+1
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)+VECBIVI0(j3-j2+ip:j3)
        IRN(i1-j2+ip:i1)=BIVIINDX1(ip:j2)/ncol+1
        JCN(i1-j2+ip:i1)=MOD(BIVIINDX1(ip:j2),ncol)+1
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        e3=j3-j2+ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)+VECBIVI0(e3:j3)
        IRN(e1:i1)=BIVIINDX1(e2:j2)/ncol+1
        JCN(e1:i1)=MOD(BIVIINDX1(e2:j2),ncol)+1
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          IRN(b1:e1-1)=BIVIINDX1(b2:e2-1)/ncol+1
          JCN(b1:e1-1)=MOD(BIVIINDX1(b2:e2-1),ncol)+1
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)+VECBIVI0(j2p:j3-j2+i)
        IRN(ip:ip+i-j1p)=BIVIINDX1(j1p:i)/ncol+1
        JCN(ip:ip+i-j1p)=MOD(BIVIINDX1(j1p:i),ncol)+1
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        e3=j2p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)+VECBIVI0(e3:j3-j2+i)
        IRN(e1:ip+i-j1p)=BIVIINDX1(e2:i)/ncol+1
        JCN(e1:ip+i-j1p)=MOD(BIVIINDX1(e2:i),ncol)+1
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          b3=e3-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)+VECBIVI0(b3:e3-1)
          IRN(b1:e1-1)=BIVIINDX1(b2:e2-1)/ncol+1
          JCN(b1:e1-1)=MOD(BIVIINDX1(b2:e2-1),ncol)+1
          e1=b1
          e2=b2
          e3=b3
        end do
      end if
      end if
    end if
    i1=i
    j2=j
    do while(BIVIINDX1(j).GT.BIVIINDX0(j1))
          j=j-1
    end do
    i=i1-(j2-j)
    if(i.NE.i1) then
    ip=i+1
    j1p=j+1

    if(ip.GT.j2) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI1(j1p:j2)
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI1(b2:e2)
        end do
          VECBIVI1(e1+1:i1)=VECBIVI1(e2+1:j2)
      end if
    else
      if(j2-ip.LT.buffer) then
        VECBIVI1(i1-j2+ip:i1)=VECBIVI1(ip:j2)
      else
        l=(j2-ip)/buffer
        e1=i1-j2+ip+buffer*l
        e2=ip+buffer*l
        VECBIVI1(e1:i1)=VECBIVI1(e2:j2)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
      if(i-j1p.LT.buffer) then
        VECBIVI1(ip:ip+i-j1p)=VECBIVI1(j1p:i)
      else
        l=(i-j1p)/buffer
        e1=ip+buffer*l
        e2=j1p+buffer*l
        VECBIVI1(e1:ip+i-j1p)=VECBIVI1(e2:i)
        do m=l-1,0,-1
          b1=e1-buffer
          b2=e2-buffer
          VECBIVI1(b1:e1-1)=VECBIVI1(b2:e2-1)
          e1=b1
          e2=b2
        end do
      end if
    end if

      if(i1-ip.LT.buffer) then
        IRN(ip:i1)=BIVIINDX1(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX1(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          IRN(b1:e1)=BIVIINDX1(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX1(b2:e2),ncol)+1
        end do
          IRN(e1+1:i1)=BIVIINDX1(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX1(e2+1:j2),ncol)+1
      end if
    end if
    i1=i
    j2=j1
    do while(BIVIINDX1(j).LT.BIVIINDX0(j1))
    j1=j1-1
    end do
    i=i1-(j2-j1)
    ip=i+1
    j1p=j1+1
    if(i.NE.i1) then
      if(i1-ip.LT.buffer) then
        VECBIVI1(ip:i1)=VECBIVI0(j1p:j2)
        IRN(ip:i1)=BIVIINDX0(j1p:j2)/ncol+1
        JCN(ip:i1)=MOD(BIVIINDX0(j1p:j2),ncol)+1
      else
        l=(i1-ip)/buffer
        do m=0,l-1,1
          b1=ip+buffer*m
          e1=b1+buffer-1
          b2=j1p+buffer*m
          e2=b2+buffer-1
          VECBIVI1(b1:e1)=VECBIVI0(b2:e2)
          IRN(b1:e1)=BIVIINDX0(b2:e2)/ncol+1
          JCN(b1:e1)=MOD(BIVIINDX0(b2:e2),ncol)+1
        end do
          VECBIVI1(e1+1:i1)=VECBIVI0(e2+1:j2)
          IRN(e1+1:i1)=BIVIINDX0(e2+1:j2)/ncol+1
          JCN(e1+1:i1)=MOD(BIVIINDX0(e2+1:j2),ncol)+1
      end if
    end if
    !20 CONTINUE
  end do
END SUBROUTINE MY_SPAR_ADD1L

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



SUBROUTINE MY_SPAR_ADD2(VECBIVI1,BIVIINDX1,IRN,JCN,NZ1,VECBIVI0,BIVIINDX0,NZ0,NZ2,NCOL1,VECBIVI2,IRN2,JCN2,J2,CNTL3)
  use constants
  IMPLICIT NONE
  real(kind=DPC) VECBIVI1(*),VECBIVI0(*),VECBIVI2(*),CNTL3(*),cntl
  INTEGER IRN(*),JCN(*),IRN2(*),JCN2(*),NZ1(*),NZ0(*),NZ2(*),NCOL1(*),j,j1,i,ncol,l,J2(*)
  INTEGER BIVIINDX0(0:NZ0(1)),BIVIINDX1(0:NZ1(1))
  j=NZ1(1)
  j1=NZ0(1)
  ncol=NCOL1(1)
  l=NZ2(1)
  cntl=cntl3(1)
  do i=NZ2(1),1,-1
    IF(BIVIINDX1(j).EQ.BIVIINDX0(j1)) then
          VECBIVI1(i)=VECBIVI1(j)+VECBIVI0(j1)
          IRN(i)=BIVIINDX1(j)/ncol+1
          JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j1=j1-1
          j=j-1
          GO TO 20
    END IF
    IF(BIVIINDX1(j).GT.BIVIINDX0(j1)) then
          VECBIVI1(i)=VECBIVI1(j)
          IRN(i)=BIVIINDX1(j)/ncol+1
          JCN(i)=MOD(BIVIINDX1(j),ncol)+1
          j=j-1
          GO TO 20
    END IF
    VECBIVI1(i)=VECBIVI0(j1)
    IRN(i)=BIVIINDX0(j1)/ncol+1
    JCN(i)=MOD(BIVIINDX0(j1),ncol)+1
    j1=j1-1
    20 CONTINUE
    IF(VECBIVI1(i).GT.cntl.OR.VECBIVI1(i).LT.-cntl) then
          VECBIVI2(l)=VECBIVI1(i)
          IRN2(l)=IRN(i)
          JCN2(l)=JCN(i)
          l=l-1
    END IF
  end do
  J2(1)=l-1
END SUBROUTINE MY_SPAR_ADD2

SUBROUTINE MY_SPAR_COMP(BIVIINDX1,NZ1,BIVIINDX0,NZ0,NZ2)
  INTEGER BIVIINDX0(*),BIVIINDX1(*),NZ1(*),NZ0(*),NZ2(*),j,j1,i
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
END SUBROUTINE MY_SPAR_COMP

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

SUBROUTINE MY_VEC_COMZ(VECBIVI,BIVIINDX,col,row,colsize,NZ0,NZ1)
  use constants
  IMPLICIT NONE
  INTEGER, PARAMETER :: myint = selected_int_kind(4)
  real(kind=DPC) VECBIVI(*)
  INTEGER BIVIINDX(*),NZ1(*),col(*),row(*),colsize(*),i,j,l
  integer (kind=myint) NZ0(*)
  j=1
  do i=1,NZ0(1)
    if(VECBIVI(i).NE.0) then
      VECBIVI(j)=VECBIVI(i)
      l=i-1
      BIVIINDX(j)=col(MOD(l,colsize(1))+1)+row(l/colsize(1)+1)
      j=j+1
    end if
  end do
  NZ1(1)=j-1
END SUBROUTINE MY_VEC_COMZ

SUBROUTINE PREP48_ALU(INSIZE,IRN,JCN,VA)
  use constants
  IMPLICIT NONE
  integer NEFAC,JOB
  integer(4) JCN(*),IRN(*),INSIZE(*)
  integer M,N,NE,T,RANK,J1,J
  real(kind=DPC) VA(*)
  logical isopen
  !DOUBLE PRECISION LOGDET,SGNDET
  integer LA, MAXN
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:)
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:)
  character(len=1024) :: filename
  M=INSIZE(1)
  N=INSIZE(2)
  NE=INSIZE(3)
  NEFAC=INSIZE(10)
  RANK=INSIZE(11)
  J1=INSIZE(12)
  if(NEFAC.EQ.0) then
    LA=2*NE
  else
    LA=ceiling((NEFAC/100.0)*NE)
  endif
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10))
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  allocate(KEEP(T))
  KEEP=0
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
  END IF
  J=INSIZE(16)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO,W)
  deallocate(ICNTL,INFO,IW,KEEP)
END SUBROUTINE PREP48_ALU

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
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
  END IF
  J=INSIZE(16)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO)!,W)
  deallocate(ICNTL,INFO)!,IW,KEEP)
END SUBROUTINE PREP48_ALU1

SUBROUTINE PREP48_MSOL(INSIZE,IRN,JCN,VA,IRNC,JCNC,VAC,IRNB,JCNB,VALUESB,VECBIVI,BIVINZROW0,BIVINZCOL0)!,IRNV,JCNV,VAV
  use constants
  IMPLICIT NONE
  integer NEFAC,JOB
  integer(8) BIVINZROW0(*),L3
  integer(4) JCN(*),IRN(*),INSIZE(*),BIVINZCOL0(*)!,JCNV(*),IRNV(*)
  integer(4)IRNB(*),JCNB(*)
  integer(4) JCNC(*),IRNC(*)
  integer M,N,NE,T,NC,MC,NEC,RANK,J1!,NV,MV,NEV
  real(kind=DPC) VA(*),VAC(*),VALUESB(*),VECBIVI(*)
  !DOUBLE PRECISION LOGDET,SGNDET!,VAV(*)real (8)
  integer LA, MAXN,NBIVI,MBIVI
  integer(4) I,J,L,L1,L2,L4,L5,M0,M1,M2,M3,M4,M5,MB,NB,NEB,J2,J3
  LOGICAL TRANS,checksol,isopen
  real(kind=DPC), pointer :: CNTL(:),RINFO(:),W(:),ERROR1(:),SOL(:),B(:)!,VOUT(:)!,RHS(:)A(:),
  integer, pointer :: ICNTL(:),INFO(:),IW(:),KEEP(:),JCNB1(:)!,IRNOUT(:),JCNOUT(:)!,JCN1(:),IRN1(:)
  character(len=1024) :: filename
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
  MAXN=N
  IF (N.LT.M) THEN
    MAXN=M
  END IF
  allocate(CNTL(10),RINFO(10),ERROR1(3),SOL(M),B(M))!,RHS(M),W(5*MAXN)A(LA),
  allocate(ICNTL(20),INFO(20),IW(6*M+3*N))!JCN1(LA),IRN1(LA),
  IF (FSORD.EQ.1) THEN
    CALL MA48ID(CNTL,ICNTL)
  else
    CALL MA48I(CNTL,ICNTL)
  endif
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  allocate(KEEP(T))
  KEEP=0
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  deallocate(IW)
  allocate(W(M),IW(2*M+2*N))
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
  END IF
  JOB=1
  deallocate(W,IW)
  IF (JOB.EQ.1) THEN
    allocate(W(2*MAXN))
  else
    allocate(W(4*MAXN))
  END IF
  allocate(IW(MAXN))!A(LA),IRN1(LA),
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
  J=INSIZE(16)
    if(J.GT.98) then 
      J=7
    end if
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO,W,ERROR1,SOL,JCNB1,B)!,RHSA,
  deallocate(ICNTL,INFO,IW,KEEP)!IRN1,,JCNOUT
END SUBROUTINE PREP48_MSOL

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
  T=M+5*N+4*N/ICNTL(6)+7
  INSIZE(13)=T
  JOB=1
  IF (FSORD.EQ.1) THEN
    CALL MA48AD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  else
    CALL MA48A(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,RINFO)
  endif
  IF (INFO(1).LT.0) THEN
  WRITE (6,'(A,I3)') 'Error STOP from MA48A/AD with INFO(1) =',INFO(1)
  STOP
  END IF
  IF (FSORD.EQ.1) THEN
    CALL MA48BD(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  else
    CALL MA48B(M,N,NE,1,LA,VA,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,INFO,&
                RINFO)
  endif
  IF (INFO(1).NE.0) THEN
    WRITE (6,FMT='(A,I3/A)') 'STOP from MA48B/BD with INFO(1) =',&
    INFO(1),'Solution not possible'
    write(*,"(A,i5)") 'RANK',INFO(5)
    STOP
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
  J=INSIZE(16)
    if(J.GT.98) then 
      J=7
    end if
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A4,I4.4,I4.4,A4)") "_vav",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) VA(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_irnv",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) IRN(1:LA)
  CLOSE(J)
  inquire(unit=J, opened=isopen)
  DO WHILE (isopen)
    J=J+1
    if(J.EQ.99) then 
      J=7
    end if
    inquire(unit=J, opened=isopen)
  END DO
  write (filename, "(A5,I4.4,I4.4,A4)") "_keep",RANK,J1,".bin"
  OPEN(J, FILE=filename, STATUS='UNKNOWN', ACCESS='STREAM')
  WRITE (J) KEEP(1:T)
  CLOSE(J)
  deallocate(CNTL,RINFO,ERROR1)!,RHSA,,SOL,JCNB1,B,W
  deallocate(ICNTL,INFO)!IRN1,,JCNOUT,IW,KEEP
END SUBROUTINE PREP48M_MSOL

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

SUBROUTINE SPAR_MULADD(SOL,NROW,NZ,IRN,JCN,VA,RES)
  use constants
  IMPLICIT NONE
  integer I,J
  integer(4) JCN(*),IRN(*),NROW(*),NZ(*)
  real(kind=DPC) VA(*),RES(*),SOL(0:NROW(1)-1)
      do I=1,NROW(1)-1
              DO J=IRN(I)+1,IRN(I+1)
                RES(I)=RES(I)+VA(j)*SOL(JCN(J))
              end do
      end do
              DO J=IRN(NROW(1))+1,NZ(1)
                RES(NROW(1))=RES(NROW(1))+VA(j)*SOL(JCN(J))
              end do
END SUBROUTINE SPAR_MULADD

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
