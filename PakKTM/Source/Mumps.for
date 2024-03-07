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
