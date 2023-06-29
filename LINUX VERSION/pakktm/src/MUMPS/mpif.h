!
!  This file is part of MUMPS 5.0.2, released
!  on Fri Jul 15 09:12:54 UTC 2016
!
!
!  Copyright 1991-2016 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
!  University of Bordeaux.
!
!  This version of MUMPS is provided to you free of charge. It is
!  released under the CeCILL-C license:
!  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
!
!
!      Dummy mpif.h file including symbols used by MUMPS.
!
      INTEGER MPI_2DOUBLE_PRECISION
      INTEGER MPI_2INTEGER
      INTEGER MPI_2REAL
      INTEGER MPI_ANY_SOURCE
      INTEGER MPI_ANY_TAG
      INTEGER MPI_BYTE
      INTEGER MPI_CHARACTER
      INTEGER MPI_COMM_NULL
      INTEGER MPI_COMM_WORLD
      INTEGER MPI_COMPLEX
      INTEGER MPI_DOUBLE_COMPLEX
      INTEGER MPI_DOUBLE_PRECISION
      INTEGER MPI_INTEGER
      INTEGER MPI_LOGICAL
      INTEGER MPI_MAX
      INTEGER MPI_MAX_PROCESSOR_NAME
      INTEGER MPI_MAXLOC
      INTEGER MPI_MIN
      INTEGER MPI_MINLOC
      INTEGER MPI_PACKED
      INTEGER MPI_PROD
      INTEGER MPI_REAL
      INTEGER MPI_REPLACE
      INTEGER MPI_REQUEST_NULL
      INTEGER MPI_SOURCE
      INTEGER MPI_STATUS_SIZE
      INTEGER MPI_SUM
      INTEGER MPI_TAG
      INTEGER MPI_UNDEFINED
      INTEGER MPI_WTIME_IS_GLOBAL
      INTEGER MPI_LOR
      INTEGER MPI_LAND
      INTEGER MPI_INTEGER8
      INTEGER MPI_REAL8
      INTEGER MPI_BSEND_OVERHEAD
      PARAMETER (MPI_2DOUBLE_PRECISION=1)
      PARAMETER (MPI_2INTEGER=2)
      PARAMETER (MPI_2REAL=3)
      PARAMETER (MPI_ANY_SOURCE=4)
      PARAMETER (MPI_ANY_TAG=5)
      PARAMETER (MPI_BYTE=6)
      PARAMETER (MPI_CHARACTER=7)
      PARAMETER (MPI_COMM_NULL=8)
      PARAMETER (MPI_COMM_WORLD=9)
      PARAMETER (MPI_COMPLEX=10)
      PARAMETER (MPI_DOUBLE_COMPLEX=11)
      PARAMETER (MPI_DOUBLE_PRECISION=12)
      PARAMETER (MPI_INTEGER=13)
      PARAMETER (MPI_LOGICAL=14)
      PARAMETER (MPI_MAX=15)
      PARAMETER (MPI_MAX_PROCESSOR_NAME=31)
      PARAMETER (MPI_MAXLOC=16)
      PARAMETER (MPI_MIN=17)
      PARAMETER (MPI_MINLOC=18)
      PARAMETER (MPI_PACKED=19)
      PARAMETER (MPI_PROD=20)
      PARAMETER (MPI_REAL=21)
      PARAMETER (MPI_REPLACE=22)
      PARAMETER (MPI_REQUEST_NULL=23)
      PARAMETER (MPI_SOURCE=1)
      PARAMETER (MPI_STATUS_SIZE=2)
      PARAMETER (MPI_SUM=26)
      PARAMETER (MPI_TAG=2)
      PARAMETER (MPI_UNDEFINED=28)
      PARAMETER (MPI_WTIME_IS_GLOBAL=30)
      PARAMETER (MPI_LOR=31)
      PARAMETER (MPI_LAND=32)
      PARAMETER (MPI_INTEGER8=33)
      PARAMETER (MPI_REAL8=34)

      PARAMETER (MPI_BSEND_OVERHEAD=0)
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
