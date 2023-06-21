C
C Ako je solver eksterni, sve rutine u ovom fajlu iskljuciti iz build-a
C
#ifndef MUMPS_EXTERNAL

! FAKTORIZACIJA SISTEMA
!
      SUBROUTINE MUMPS_FACTOR(mumps_par,IIZLAZ,INFF)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC) mumps_par
	INTEGER IERR, I, IIZLAZ,INFF
	
!  Iskljucivanje outputa
	mumps_par%ICNTL(1)=-1 
	mumps_par%ICNTL(2)=-1
	mumps_par%ICNTL(3)=-1
	mumps_par%ICNTL(4)=0

!	mumps_par%ICNTL(22)=1
!	mumps_par%DMUMPS_OOC_TMPDIR = 'C:\Users\miljan\Documents\3D remeshing
!     &\TOT REMESHING\000 MovingDomain\PAKSF\PakF\'

!  Call package for analysis and factorization
      mumps_par%ICNTL(14)=40
      IF(INFF.EQ.1) WRITE(*,*)'Starting analysis...'
	mumps_par%JOB = 1
      CALL DMUMPS(mumps_par)
    !  mumps_par%ICNTL(14)=1200
      mumps_par%JOB = 2
	IF(INFF.EQ.1)  WRITE(*,*)'Starting factorization...'
      CALL DMUMPS(mumps_par)

!  Greska u izvrsavanju solvera, dealokacija i kraj
	IF ( mumps_par%INFO(1).NE.0 ) THEN
        WRITE(*,*) 'INFO1:',mumps_par%INFO(1)
        WRITE(*,*) 'INFO2:',mumps_par%INFO(2)
        WRITE(*,*) 'INFO3:',mumps_par%INFO(3)
        WRITE(*,*) 'INFO7:',mumps_par%INFO(7)
        WRITE(*,*) 'INFO8:',mumps_par%INFO(8)
        WRITE(*,*) 'INFO9:',mumps_par%INFO(9)
        WRITE(*,*) 'INFO15:',mumps_par%INFO(15)
        WRITE(*,*) 'ICNTL14:',mumps_par%ICNTL(14)

        WRITE(IIZLAZ,*) 'INFO1  :',mumps_par%INFO(1)
        WRITE(IIZLAZ,*) 'INFO2  :',mumps_par%INFO(2)
        WRITE(IIZLAZ,*) 'INFO3  :',mumps_par%INFO(3)
        WRITE(IIZLAZ,*) 'INFO7  :',mumps_par%INFO(7)
        WRITE(IIZLAZ,*) 'INFO8  :',mumps_par%INFO(8)
        WRITE(IIZLAZ,*) 'INFO9  :',mumps_par%INFO(9)
        WRITE(IIZLAZ,*) 'INFO15 :',mumps_par%INFO(15)
        WRITE(IIZLAZ,*) 'ICNTL14:',mumps_par%ICNTL(14)

        CALL MUMPS_END(INFF)
	  WRITE(*,*) 'ERROR IN SOLVER EXECUTION - STOPPING CALCULATION'
	  STOP
      ENDIF
	END

!  ================================================================
!  ================================================================
!  ================================================================

      SUBROUTINE MUMPS_SOLUTION(RHS)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
	COMMON /MUMPS/ mumps_par
      TYPE (DMUMPS_STRUC) mumps_par
	DOUBLE PRECISION RHS
	DIMENSION RHS(*)
	INTEGER IERR, I
      
	DO I = 1, mumps_par%N
        mumps_par%RHS(I) = RHS(I)
      END DO

!  Call package for solution

      mumps_par%JOB = 3
      CALL DMUMPS(mumps_par)

	DO I = 1, mumps_par%N
        RHS(I) = mumps_par%RHS(I) 
      END DO
	END 


! =================================================================
! =================================================================
! =================================================================
      SUBROUTINE MUMPS_INIT(N,NZ,INFF)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
	COMMON /MUMPS/ mumps_par
      TYPE (DMUMPS_STRUC) mumps_par
	INTEGER IERR,N,NZ,INFF
	
	IF(INFF.EQ.1) WRITE(*,*) "Initializing solver..."
	
      CALL MPI_INIT(IERR)
      mumps_par%COMM = MPI_COMM_WORLD

!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
      mumps_par%JOB = -1
      mumps_par%SYM = 0
      mumps_par%PAR = 1

	CALL DMUMPS(mumps_par)

      mumps_par%N = N
	mumps_par%NZ = NZ

      ALLOCATE( mumps_par%IRN ( NZ ))
      ALLOCATE( mumps_par%JCN ( NZ ))
      ALLOCATE( mumps_par%A( NZ ))
      ALLOCATE( mumps_par%RHS ( N ))
	
	END

! =================================================================
! =================================================================
! =================================================================
      SUBROUTINE MUMPS_END(INFF)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
	COMMON /MUMPS/ mumps_par      
	TYPE (DMUMPS_STRUC) mumps_par
	INTEGER IERR,INFF
	
      IF(INFF.EQ.1) WRITE(*,*) "Killing solver..."

!  Destroy the instance (deallocate internal data structures)
      mumps_par%JOB = -2
      CALL DMUMPS(mumps_par)

      DEALLOCATE( mumps_par%IRN )
      DEALLOCATE( mumps_par%JCN )
      DEALLOCATE( mumps_par%A   )
      DEALLOCATE( mumps_par%RHS )

      CALL MPI_FINALIZE(IERR)

	END

! =========================================================================
! =========================================================================
! ROUTINE FOR READING THE SPARSE FILE FROM SPARSEASSEMBLER
      SUBROUTINE MUMPS_READ_SPARSE_FILE(NZ)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
	COMMON /MUMPS/ mumps_par      
	TYPE (DMUMPS_STRUC) mumps_par
      INTEGER NZ

!     Local variables
      INTEGER I,IFILE
      IFILE = 989

      OPEN (IFILE, FILE='Sparse.bin', FORM='BINARY')
      DO I = 1, NZ
        READ(IFILE) mumps_par%IRN(I),mumps_par%JCN(I),mumps_par%A(I)
      END DO

      CLOSE (IFILE)
      
      END

C=========================================================================
C=========================================================================
      SUBROUTINE MUMPSINIT(NEL,ID,JEDN,NET,NDIM,NETIP)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON /MUMPS_PAK/ IMUMPS,MUFILE,MUFILE2

      DIMENSION NEL(NDIM+1,*),ID(6,*)
	DIMENSION LM2(NETIP*NDIM)

      IF (IMUMPS.EQ.0) RETURN

       REWIND(MUFILE)

	NN=NETIP*NDIM

	WRITE(MUFILE) JEDN
	WRITE(MUFILE) NET

      n=1
	N1=0
      DO NBREL=1,NET
	N0=0
C	  WRITE(mufile,*) n
      DO II=1,NETIP*NDIM
        LM2(II)=0
      ENDDO
C=========================================================================
      DO KLM=1,NDIM
	  DO J=1,NETIP
          LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
 	  ENDDO
	ENDDO

      DO JJ=1,NN
	  IF (LM2(JJ).NE.0) THEN
	     N=N+1
	   N0=N0+1
        ENDIF
	ENDDO
	N1=N1+N0*N0  
	enddo
C	  WRITE(mufile,*) n

	WRITE(MUFILE) N-1
	WRITE(MUFILE) N1


      n=1
      DO NBREL=1,NET
	  WRITE(mufile) n
      DO II=1,NETIP*NDIM
        LM2(II)=0
      ENDDO
C=========================================================================
      DO KLM=1,NDIM
	  DO J=1,NETIP
          LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
 	  ENDDO
      ENDDO

      DO JJ=1,NN
	  IF (LM2(JJ).NE.0) N=N+1
	ENDDO
	  
	enddo
	  WRITE(mufile) n

C	WRITE(MUFILE,*) N-1

      DO 400 NBREL=1,NET

C TT21 is vector of unknowns values at element level
C TT210 is vector of unknowns values at element level at start of time step

      DO 125 I=1,NETIP*NDIM
        LM2(I)=0
 125  CONTINUE
C=========================================================================
      DO 130 KLM=1,NDIM
	 DO J=1,NETIP
        LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
	 ENDDO
 130  CONTINUE

      do i=1,nn
C	      WRITE(MUFILE,*) (LM2(I),I=1,NN)
	 IF (LM2(I).NE.0) WRITE(MUFILE) LM2(I)
      enddo
      
400   CONTINUE


	END
C=========================================================================
C
C=========================================================================
      SUBROUTINE MUMPSLEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,NETIP,NDIMID)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON /MUMPS_PAK/ IMUMPS,MUFILE,MUFILE2
      DIMENSION NEL(NDIM+1,*),ID(NDIMID,*)
c      DIMENSION NEL(9,*),ID(6,*)
	DIMENSION LM2(NETIP*NDIM)
      DOUBLE PRECISION elembuffer(67108864)
      INTEGER BUFFERSIZE

      DIMENSION SKEF(NDES,*)

      IF (IMUMPS.EQ.0) RETURN


      DO II=1,NETIP*NDIM
        LM2(II)=0
      ENDDO
C=========================================================================
      DO KLM=1,NDIM
	  DO J=1,NETIP
          LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
 	  ENDDO
	ENDDO

C      DO I=1,NETIP*ndim
C      DO J=1,NETIP*ndim
C	 IF (LM2(J).NE.0.AND.lm2(i).ne.0) then
C	   WRITE(MUFILE2) SKEF(J,I)
C       endif
C	ENDDO
C	ENDDO

C ========================================================================
C Milos dodao baferisanje od 256MB zbog performansi, Avgust 2009.
      BUFFERSIZE = 67108864
      DO KLM=1,NDIM
          DO J=1,NETIP
          LM2(KLM+(J-1)*NDIM)=ID(J,NEL(KLM,NBREL))
          ENDDO
        ENDDO

      buffercounter = 0
      DO I=1,NETIP*ndim
      DO J=1,NETIP*ndim
         if (LM2(J).NE.0.AND.lm2(i).ne.0) then
           buffercounter = buffercounter+1
           elembuffer(buffercounter) = SKEF(J,I)
         endif
         if (buffercounter.eq.BUFFERSIZE) then
           write(mufile2) (elembuffer(k), k=1,BUFFERSIZE)
           buffercounter = 0
         endif
      ENDDO
      ENDDO

      write(mufile2) (elembuffer(k), k=1, buffercounter)

	END         
C=========================================================================
C
C=========================================================================
      SUBROUTINE MUMPSRIGHT(A,N)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON /MUMPS_PAK/ IMUMPS,MUFILE,MUFILE2


      DIMENSION A(*)

      IF (IMUMPS.EQ.0) RETURN
	
C      do i=1,n
c	WRITE(MUFILE,*) (A(I),I=1,N)
 	WRITE(MUFILE2) (A(I), I=1,N)
C	enddo
C        CLOSE(MUFILE2)
	END
C=========================================================================
C
C=======================================================================
      SUBROUTINE checksolver(sile,jedn,right)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	DIMENSION SILE(*), right(*)


      iKc=0
      mufile2=107
      IIZLAZ=78
      OPEN (MUFILE2, FILE='FULLMATR.BIN',FORM='BINARY')
      OPEN (iizlaz, FILE='comparison.txt')

       
	read(mufile2) Jedn
	read(mufile2) iKc
	do i=1,jedn
	 right(i)=-right(i)
	enddo
      do k=1,iKc
       read(mufile2)i,j,a
       right(i)=right(i)+a*sile(j)
	enddo

 200  FORMAT (I5,1PD20.12)

      DO I=1,JEDN
        WRITE(IIZLAZ,200) i,right(i)
	ENDDO

      close (mufile2)
      close (iizlaz)

      END
C=======================================================================
C==========================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE FULLMA(ALEVO,DESNO,MaxKs,JEDN,SILE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION ALEVO(*),DESNO(*),sile(*)
	dimension MaxKs(*)


      iKc=0
      mufile2=107
	tol=1.0D-20
c      OPEN(iizlaz,FILE='PUNAMATR.TXT')
      OPEN (MUFILE2,FILE='/tmp/PAKF_MUMPS.OUT',FORM='BINARY')

      iKc=0

C	WRITE(mufile2) JEDN
C	WRITE(mufile2) iKc




	DO 20 I=1,JEDN
       DO 10 J=1,JEDN
        if(i.eq.1.and.j.eq.1) then
	     iKc=iKc+1
		 goto 10 
	  endif
            IF(J.GE.I) THEN
            iKs=MaxKs(j)-j+i
            if (iKs.gt.MaxKs(j-1)) then
              if (DABS(ALEVO(iKs)).GT.tol) then
                  iKc=iKc+1
C                  WRITE(mufile2,*)I,J,ALEVO(IKS)
C                  WRITE(IIZLAZ,*)I,J,ALEVO(IKS)
C                  Kc(iKc)=ALEVO(iKs)
               endif
            endif
            ELSE
             IKS=MAXKS(I)-I+J
             if (iKs.gt.MaxKs(I-1)) then
               if (DABS(DESNO(iKs)).GT.tol) then
                  iKc=iKc+1
C                  WRITE(mufile2,*)I,J,DESNO(IKS)
C                  WRITE(IIZLAZ,*)I,J,DESNO(IKS)
C                  Kc(iKc)=DESNO(iKs)
               endif
            endif
           ENDIF
C         write(iizlaz,200)i,j,vredn
  10  continue
  20  continue


C	WRITE(IIZLAZ,*) JEDN
C	WRITE(IIZLAZ,*) iKc
c      rewind (mufile2)        
	WRITE(mufile2) JEDN
	WRITE(mufile2) iKc
c	goto 50

      iKc=0

C
CE Subroutine FULLMATR IS USED FOR PRINTING FULL MATRIX
C
 
 
 
      
	DO 55 I=1,JEDN
       DO 45 J=1,JEDN
        if(i.eq.1.and.j.eq.1) then
	     iKc=iKc+1
          WRITE(mufile2)I,J,ALEVO(1)
		 goto 45 
	  endif

            IF(J.GE.I) THEN
            iKs=MaxKs(j)-j+i
            if (iKs.gt.MaxKs(j-1)) then
               if (DABS(ALEVO(iKs)).GT.tol) then
                  iKc=iKc+1
c                  WRITE(IIZLAZ,*)I,J,ALEVO(IKS)
                  WRITE(mufile2)I,J,ALEVO(IKS)
C                  Kc(iKc)=ALEVO(iKs)
               endif
            endif
            ELSE
             IKS=MAXKS(I)-I+J
             if (iKs.gt.MaxKs(I-1)) then
               if (DABS(DESNO(iKs)).GT.tol) then
                  iKc=iKc+1
c                  WRITE(IIZLAZ,*)I,J,DESNO(IKS)
                  WRITE(mufile2)I,J,DESNO(IKS)
C                  Kc(iKc)=DESNO(iKs)
               endif
            endif
           ENDIF
C         write(iizlaz,200)i,j,vredn
   45 continue
   55 continue

 200  FORMAT (2I5,1PD20.12)

    
 

 50     DO I=1,JEDN
c        WRITE(IIZLAZ,*) SILE(I)
        WRITE(mufile2) SILE(I)
	ENDDO

c      rewind (mufile2)        
c	WRITE(mufile2) JEDN
c	WRITE(mufile2) iKc

c      REWIND (UNIT=iizlaz, IOSTAT=IOS)
c     REWIND (IIZLAZ)
c	WRITE(IIZLAZ,*) JEDN
c	WRITE(IIZLAZ,*) iKc


      CLOSE (mufile2)
       
      END
C=======================================================================
#endif


C======================================================================
C SOLVER MUMPS poziv eksternog solvera
	SUBROUTINE SOLVER_EXTERNAL(SILE, JEDN)
	IMPLICIT NONE
        DOUBLE PRECISION SILE (*)
	INTEGER JEDN
	INTEGER I, IFILE, PFILE
	CHARACTER* 300 CMDLINE
	CHARACTER*55 JOBID
	IFILE = 57
C Citanje params.txt, $JOBID treci broj
	PFILE = 61
	OPEN (PFILE, FILE='params.txt')
	READ(PFILE,*) JOBID
	READ(PFILE,*) JOBID
	READ(PFILE,*) JOBID
	CLOSE (PFILE)
C Kraj $JOBID


!       CMDLINE = "mumps_convert"
!       CALL SYSTEM(CMDLINE)
 

!        CMDLINE = "mpiexec solver"
!        CMDLINE = "mpiexec elemental"
!        CMDLINE = "PakC_Elemental.exe"
	CMDLINE = "mpiexec -np 32 solver"


	CALL SYSTEM(CMDLINE)


	OPEN (IFILE, FILE='/tmp/RESENJE_'//JOBID//'.OUT',FORM='BINARY')
!	OPEN (IFILE, FILE='RESENJE.OUT',FORM='BINARY')
	DO I = 1, JEDN
	READ(IFILE) SILE(I)
	END DO
	CLOSE (IFILE)

	END
C=========================================================================
#ifndef MUMPS_EXTERNAL
C=========================================================================
C SOLVER MUMPS poziv internog solvera
	SUBROUTINE SOLVER_INTERNAL(SILE, JEDN, iSymmetric,IIZLAZ,INFF)
	IMPLICIT NONE

      INCLUDE 'dmumps_struc.h'

	COMMON /MUMPS/ mumps_par
      TYPE (DMUMPS_STRUC) mumps_par

      DOUBLE PRECISION SILE (*)
	INTEGER JEDN, ISPARSE_N, ISPARSE_NZ, iSymmetric,IIZLAZ,INFF

	ISPARSE_N = JEDN
        
      CALL sparseassembler_getnonzero(ISPARSE_NZ)
      CALL MUMPS_INIT(ISPARSE_N, ISPARSE_NZ,INFF) !, iSymmetric)
      CALL sparseassembler_getsparse(ISPARSE_NZ, mumps_par%IRN,
     1  mumps_par%JCN,mumps_par%A,INFF)
      CALL sparseassembler_kill(INFF)
      if(ISPARSE_N.gt.0.and.ISPARSE_NZ.gt.0)then
        CALL MUMPS_FACTOR(mumps_par,IIZLAZ,INFF)
        CALL MUMPS_SOLUTION(SILE)
        CALL MUMPS_END(INFF)
      else
        write(*,*) 'Error in solver calling!'
        stop
      endif

	END
C======================================================================



#endif
