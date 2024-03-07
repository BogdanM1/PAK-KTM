C=======================================================================

C=======================================================================
      SUBROUTINE TABF(FN,NTAV,IBR,NFN,ARG,Y,NTMX,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
CE   PROGRAM TO DETERMINE VELUES OF FUNCTION GIVEN AS TABLE FOR
CE   ARGUMENT ARG
      DIMENSION FN(2,NFN,*),NTAV(*)
	tol=1.d-10
      IND=0
      NTMX=NTAV(IBR)
      IF(NTMX.GT.1) GO TO 10
      Y=FN(2,IBR,1)
      RETURN
CE   MULTILINEAR CURVE
   10 DO 50 I=2,NTMX
      I1=I-1
      IF(ARG.GE.FN(1,IBR,I1).AND.ARG.LT.FN(1,IBR,I)) GO TO 20
      GO TO 50
   20 Y=FN(2,IBR,I1)+(FN(2,IBR,I)-FN(2,IBR,I1))/
     *(FN(1,IBR,I)-FN(1,IBR,I1))*(ARG-FN(1,IBR,I1))
      GO TO 100
   50 CONTINUE
	if((dabs(arg-FN(1,IBR,1)).le.tol).or.
	1   (dabs(arg-FN(1,IBR,NTMX)).le.tol)) return
      IND=1
  100 RETURN
      END
C=======================================================================

C=======================================================================
      SUBROUTINE tabdiff(MATERP,AMATF,mat,arg,ddifus,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
CE   PROGRAM TO DETERMINE VELUES OF diffusion coeff given AS TABLE FOR
CE   ARGUMENT ARG
	common/diffusmain/indiffus,NDS,maxdt,NPRINTnds,lvreme,leprdt,
	1                  NUMMAT,maxPmat,LMATER,LMATFUN,LAREAD,LFLUXDX,
     1                  LFLUXDY,LFLUXDZ        
      DIMENSION MATERP(NUMMAT),AMATF(NUMMAT,2,*)
c      DIMENSION FN(2,NFN,*),NTAV(*)
      IND=0
	tol=1.d-10
      NTMX=MATERP(mat)
      IF(NTMX.GT.1) GO TO 10
      ddifus=AMATF(mat,2,1)
      RETURN
CE   MULTILINEAR CURVE
   10 DO 50 I=2,NTMX
      I1=I-1
      IF(ARG.GE.AMATF(mat,1,i1).AND.ARG.LT.AMATF(mat,1,i)) GO TO 20
      GO TO 50
   20 ddifus=AMATF(mat,2,i1)+(AMATF(mat,2,i)-AMATF(mat,2,i1))/
     *(AMATF(mat,1,i)-AMATF(mat,1,i1))*(ARG-AMATF(mat,1,i1))
      GO TO 100
   50 CONTINUE
      if(dabs(AMATF(mat,1,1)-arg).le.tol) then
	  ddifus=AMATF(mat,2,1)
	  go to 100
	endif
      if(dabs(AMATF(mat,1,NTMX)-arg).le.tol)then
	  ddifus=AMATF(mat,2,NTMX)
	  go to 100
	endif
      IND=1
  100 RETURN
      END
C===================================================================
C
C===================================================================
      SUBROUTINE MAXATE(MAXA,MHT,IDT,NEL,LM,NE,JEDN,INDMA,NTE,NWK,
     1 NDIMIDT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
CE    PROGRAM TO DETERMINE COLUMN HEIGHTS VECTOR AND MAXA
      DIMENSION MAXA(JEDN+1),MHT(*),IDT(NDIMIDT,*),NEL(NE,*),LM(*)
CE    ELEMENT LOOP
      DO 100 NLM=1,NE
      KK=0
      DO 2 I=1,NTE
      IF(NEL(NLM,I).EQ.0) GO TO 2
      N=NEL(NLM,I)
      DO J=1,2
        IF (IDT(J,N).GT.0) THEN
          KK=KK+1
          LM(KK)=IDT(J,N)
        ENDIF
      ENDDO
    2 CONTINUE
      LS=JEDN+1
      DO 10 I=1,KK
      IF (LM(I).LT.LS) LS=LM(I)
   10 CONTINUE
      DO 20 I=1,KK
      II=LM(I)
      ME=II-LS
      IF(ME.GT.MHT(II)) MHT(II)=ME
   20 CONTINUE
  100 CONTINUE
      IF (INDMA.EQ.0) RETURN
CE    VECTOR MAXA
      MAXA(1)=1
      MAXA(2)=2
      DO 200 I=2,JEDN
  200 MAXA(I+1)=MAXA(I)+MHT(I)+1
      NWK=MAXA(JEDN+1)-1
      LS = JEDN+1
      DO 210 I=1,LS
  210 MHT(I)=MAXA(I)
 3001 FORMAT(6(I5,2X,I5))
C
      RETURN    
      END
C=======================================================================

C=======================================================================
      SUBROUTINE RESEN(A,V,MAXA,NN,NWK,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /SRPSKI/ ISRPS
      DIMENSION A(*),V(*),MAXA(*)
CE     L*D*L(T) FACTORIZATION
      IF(KKK-2)40,150,150
   40 DO 140 N=1,NN
      KN=MAXA(N)
      KL=KN+1
      KU=MAXA(N+1)-1
      KH=KU-KL
      IF(KH)110,90,50
   50 K=N-KH
      IC=0
      KLT=KU
      DO 80 J=1,KH
      IC=IC+1
      KLT=KLT-1
      KI=MAXA(K)
      ND=MAXA(K+1)-KI-1
      IF(ND)80,80,60
   60 KK=MIN0(IC,ND)
      C=0.
      DO 70 L=1,KK
   70 C=C+A(KI+L)*A(KLT+L)
      A(KLT)=A(KLT)-C
   80 K=K+1
   90 K=N
      B=0.
      DO 100 KK=KL,KU
      K=K-1
      KI=MAXA(K)
      C=A(KK)/A(KI)
      B=B+C*A(KK)
  100 A(KK)=C
      A(KN)=A(KN)-B
  110 IF(A(KN).LE. 0.) THEN
        IF(ISRPS.EQ.0)
     *  WRITE(3,2000)N,A(KN)
        IF(ISRPS.EQ.1)
     *  WRITE(3,6000)N,A(KN)
        STOP
      ENDIF
  140 CONTINUE
      RETURN
CE     FORWARD REDUCTION
  150 DO 180 N=1,NN
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF(KU-KL)180,160,160
  160 K=N
      C=0.
      DO 170 KK=KL,KU
      K=K-1
  170 C=C+A(KK)*V(K)
      V(N)=V(N)-C
  180 CONTINUE
CE     BACK SUBSTITUTION
      DO 200 N=1,NN
      K=MAXA(N)
  200 V(N)=V(N)/A(K)
      IF(NN.EQ.1)RETURN
      N=NN
      DO 230 L=2,NN
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF (KU-KL)230,210,210
  210 K=N
      DO 220 KK=KL,KU
      K=K-1
  220 V(K)=V(K)-A(KK)*V(N)
  230 N=N-1
      RETURN
 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'
     1//' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I4,//' ','PIVOT=',
     2D20.12)
 6000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'
     1//' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I4,//' ','PIVOT=',
     2D20.12)
      END
C=======================================================================

C=======================================================================
      SUBROUTINE MINV(A,N,D,L,M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
CE    INVERT OF MATRIX
CE    INVERZION OF MATRIX
      DIMENSION L(*),M(*)
      DIMENSION A(*)
CE    SEARCHING MAX. MEMBER
      D=1.0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
      IF(DABS(BIGA)-DABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
CE    CHANGE ROW
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI)=HOLD
CE    CHANGE COLUMN
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI)=HOLD
   45 IF(BIGA) 48,46,48
   46 D=0.0
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
CE    REDUCTION OF MATRIX
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
      D=D*BIGA
      A(KK)=1.0/BIGA
   80 CONTINUE
CE    LAST
      K=N
  100 K=(K-1)
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI)=HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI)=HOLD
      GO TO 100
  150 RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UCZADP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
CE.    P R O G R A M
CE.        TO FORM POINTERS AND CALL ROUTINE FOR  READING PRESCRIBED
CE.        DISPLACEMENTS
CE.        NZADP - NUMBER OF PRESCRIBED DISPLACEMENTS
CE.    POINTERS
CE.        LNZADJ -  NUMBER OF EQUATION WITH PRESCRIBED DISPLACEMENTS
CE.        LNZADF -  NUMBER OF TIME FUNCTION
CE.        LZADFM -  MULTIPLICATION FACTOR
      COMMON A(1)
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /TEMPGL/ NPT,NGET,NMATT,NSTAC,MAXIT,NDT,NTCONS,
     1NET,NGE,JEDN,DT,TMAX,VREME,NPOC,NWK,NPRINT,LINTE,
     2KONF,INDSK,INDF,KOR,NTOTA,ITER,INDSC,INTEB,INDAX,IANIZ,NQE,IZIP,
     3ITOR,INDVEL,INDFL
      COMMON /REPEAT/ LCORD,LTABK,LTABC,LNTAKV,LNTACV,LIDTE,LIDT,LTCONS,
     1LMAX,LMAXA,LMAXM,LMHT,LNTAQE,LQEFN,LNTAQP,LQPFN,LNTAOK,LTOKFN,
     2LNTAHP,LHPFN,LNTAHR,LHRFN,LTT1,LTT2,LTT0,LSK,LCONM,LF,LFC,LNEL,
     3LIBFK,LFAKP,LISKC,LGUSM,LTOPM,LTMNM,LFOHR,LFVPOT
      COMMON/DARCYEQ/ ELDARCYX,ELDARCYY,ELDARCYZ,VDARCY,PINLET,
     1       DARCYX,JEDNDARCY,ISOLVER,IDARCY,LDARCYD,LDARCYMI,
     2       LIDARCYMA,NPTINDARCY 
CE    POINTERS FOR PRESCRIBED DISPLACEMENTS
      LNZADJ  = LMAX
      LZADFM = LNZADJ + 3*NZADP
      LMAX   = LZADFM + NZADP
      CALL ZADPOM(A(LNZADJ),A(LZADFM),A(LIDT),A(LIDARCYMA),NPT,NZADP)
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE ZADPOM(NZADJ,ZADFM,IDT,IDARCYMA,NPT,NZADP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*130 ACOZ
CE.    P R O G R A M
CE.        TO READ PRESCRIBED DISPLACEMENTS
CE.    I=1,NZADP  (NZADP - NUMBER OF PRESCRIBED DISPLACEMENTS)
CE.         NZADJ(I) -  NO. TIME FUNCTION FOR NODE - I
CE.         NZADF(I) -  NO. TIME FUNCTION FOR NODE - I
CE.         ZADFM(I) -  MULTIPLICATION FACTOR FOR NODE - I
      COMMON/CSMEAREDG/ISMEAREDG,NDIMIDG,ICONVG,INDIFFUSG,LAMIV,
     1     LDARCYCOEFV,MATNEWUSE,LHYDRCOEFV,LDPIPE,LDTISSUE,IPRINTVEL,
     2     IACTIVE(4),LZADATA1,LZADATA0,LZADAVXYZ,LWALLAREAV,
     3     LWALLPOROSV,LIPIPENODV,ILARGEVES
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
      DIMENSION NZADJ(3,NZADP),ZADFM(NZADP),IDT(NDIMIDG,NPT),
     1          IDARCYMA(NPT)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     1WRITE(3,2000)
      IF(ISRPS.EQ.1)
     1WRITE(3,6000)
      IF(ISRPS.EQ.0)
     1WRITE(3,2010)
      IF(ISRPS.EQ.1)
     1WRITE(3,6010)
      ENDIF
      I = 0
      NAUT=0
    5 I=I+1
      CALL ISPITA(ACOZ)
      IF(I.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(1,*) N,NC,FAK,KORC
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) N,IDOF,NC,FAK,KORC
      IF(ISMEAREDG.NE.1) GO TO 100
      IF(IDOF.LE.2) THEN
         IF(IDT(IDOF,N).GT.0) THEN
            WRITE(3,6050)N,IDOF,IDT(IDOF,N)
 6050       FORMAT(/' INPUT DATA ERROR'/
     1  ' FOR PRESCRIBED VALUES, DEGREE OF FREEDOM MUST BE RESTRAINED'/
     2  ' NODE NUMBER  =',I10/
     3  ' DEGREE OF FREEDOM, IDOF=',I5/
     4  ' GENERATED EQUATION NUMBER IN MATRIX IDT(IDOF,N)=',I10/
     5  ' CHECK DATA FEOR NODE CONSTRAINTM(MUST BE INPUT 1'/
     6  ' S T O P')
           STOP
         ENDIF
      ELSE
         IF(IDARCYMA(N).GE.0) THEN
            NN=0
            WRITE(3,6051)N,IDOF,NN
 6051       FORMAT(/' INPUT DATA ERROR'/
C     1  ' FOR PRESCRIBED NODAL VELOCITY, IDT(3,N) MUST BE EQ.1'/
     2  ' NODE NUMBER  =',I10/
     3  ' DEGREE OF FREEDOM, IDOF=',I5/
     4  ' INPUT VALUE FOR IDT(IDOF,N)=',I5/ 
     5  ' CHECK DATA FOR NODE CONSTRAINT (MUST BE INPUT 1)'/
     6  ' S T O P')
           STOP               
         ENDIF
      ENDIF
      GO TO 200
  100 CONTINUE
      IF(ISMEAREDG.EQ.0) THEN
         IND=1
         IF(IDOF.EQ.1.AND.IDT(1,N).EQ.0)IND=0
C         .AND.IDARCYMA(N).GE.0) IND=0
         IF(IDOF.EQ.2.AND.IDT(2,N).EQ.0)IND=0
C       .AND.IDARCYMA(N).GE.0) IND=0
         IF(IDOF.EQ.-1.AND.IDT(1,N).GT.0)IND=0
C        .AND.IDARCYMA(N).LT.0) IND=0
         IF(IDOF.EQ.-2.AND.IDT(2,N).GT.0)IND=0
C       .AND.IDARCYMA(N).LT.0) IND=0
         IF((IDOF.EQ.3.AND.IDT(3,N).EQ.0).OR.
     1      (IDOF.EQ.4.AND.IDT(4,N).EQ.0))IND=0
         IF(IND.EQ.1) THEN
            WRITE(3,6100)N,IDOF,(IDT(K,N),K=1,4)
 6100       FORMAT(/' INPUT ERROR FOR PRESCRIBED VALUE AT NODE:'/
     1     ' NODE NUMBER                                   N=',I10/
     2     ' DEGREE OF FREEDOM,                         IDOF=',I10/
     3     ' FOR IDOF.EQ.1,  IDT(1,N) MUST BE.EQ.1'/
     4     ' FOR IDOF.EQ.-1, IDT(1,N) MUST BE.EQ.0'/
     3     ' FOR IDOF.EQ.3, IDT(3,N) MUST BE.EQ.1'/
     3     ' FOR IDOF.EQ.4, IDT(4,N) MUST BE.EQ.1'/
     3     ' INPUT FOR IDT(1-4,N):'/   
     4     4(48X,'=',I10)/    
     3     ' S T O P')
           STOP
         ENDIF                
         GO TO 200  
      ENDIF
      IF(ISMEAREDG.EQ.2) THEN
         IND=1
         IF(IDOF.EQ.1.AND.IDT(1,N).EQ.0.AND.IDARCYMA(N).
     1       GE.0) IND=0
         IF(IDOF.EQ.-1.AND.IDT(1,N).GT.0.AND.IDARCYMA(N).LT.0) IND=0
         IF(IND.EQ.1) THEN
            WRITE(3,6102)N,IDOF,IDT(1,N)
 6102       FORMAT(/' INPUT ERROR FOR PRESCRIBED VALUE AT NODE:'/
     1     ' NODE NUMBER                                   N=',I10/
     2     ' DEGREE OF FREEDOM,                         IDOF=',I10/
     3     ' FOR IDOF.EQ.1,  IDT(1,N) MUST BE.EQ.1  IDT(1,N)=',I10/
     3     ' INPUT FOR IDT(1-4,N):'/   
     3     ' FOR IDOF.EQ.-1, IDT(1,N) MUST BE.EQ.0  IDT(1,N)=',I10/   
     3     ' S T O P')
           STOP
         ENDIF
         GO TO 200
      ELSE IF(ISMEAREDG.EQ.3) THEN
         IND=1
         IF(ICONVG.EQ.0) THEN
            IF(IDOF.EQ.3.AND.IDT(3,N).EQ.0) IND=0
         ELSE
             IF((IDOF.EQ.1.AND.IDT(1,N).EQ.0).OR.
     1       (IDOF.EQ.3.AND.IDT(3,N).EQ.0))IND=0
             IF(IDARCYMA(N).LT.0.AND.IDOF.EQ.-1.AND.IDT(1,N).GT.0)IND=0
         ENDIF         
         IF(IND.EQ.1) THEN
            IF(ICONVG.EQ.0)WRITE(3,6103)N,IDOF,IDT(3,N)
            IF(ICONVG.EQ.1)WRITE(3,6113)N,IDOF,IDT(1,N),IDT(3,N)
 6103       FORMAT(/' INPUT ERROR FOR PRESCRIBED VALUE AT NODE:'/
     1     ' NODE NUMBER                                 N=',I10/
     2     ' DEGREE OF FREEDOM,                       IDOF=',I10/
     3     ' MATRIX IDT......................     IDT(1,N)=',I10/
     3     ' INPUT FOR IDT(1-4,N):'/      
     3     ' FOR IDOF.EQ.3 MUST BE IDT(3,N).EQ.1'/
     4     ' FOR IDOF.EQ.-1 MUST BE IDT(1,N).EQ.0'/    
     3     ' S T O P')
 6113       FORMAT(/' INPUT ERROR FOR PRESCRIBED VALUE AT NODE:'/
     1     ' NODE NUMBER                                 N=',I10/
     2     ' DEGREE OF FREEDOM,                       IDOF=',I10/
     3     ' MATRIX IDT......................     IDT(1,N)=',I10/ 
     4     '           ......................     IDT(3,N)=',I10/
     3     ' INPUT FOR IDT(1-4,N):'/ 
     3     ' FOR IDOF.EQ. 1       MUST BE IDT(1,N).EQ.1'/
     4     ' FOR IDOF.EQ.-1       MUST BE IDT(1,N).EQ.0'/
     3     ' FOR IDOF.EQ. 3       MUST BE IDT(3,N).EQ.1'/  
     3     ' S T O P')
           STOP
         ENDIF 
         GO TO 200
      ELSE IF(ISMEAREDG.EQ.4) THEN
         IND=1
         IF((IDOF.EQ.1.AND.IDT(1,N).EQ.0).OR.
     1       (IDOF.EQ.2.AND.IDT(2,N).EQ.0))IND=0
         IF((IDOF.EQ.-1.AND.IDT(1,N).GT.0).OR.
     1       (IDOF.EQ.-2).AND.IDT(2,N).GT.0)IND=0
         IF(IND.EQ.1) THEN
            WRITE(3,6106)N,IDOF,IDT(1,N),IDT(3,N)
 6106       FORMAT(/' INPUT ERROR FOR PRESCRIBED VALUE AT NODE:'/
     1     ' NODE NUMBER                                 N=',I10/
     2     ' DEGREE OF FREEDOM,                       IDOF=',I10/
     3     ' INPUT FOR IDT(1-4,N):'/
     3     ' IDT(1,N) MUST BE.EQ.1 FOR IDOF.EQ.1, IDT(1,N)=',I10/
     3     ' IDT(2,N) MUST BE.EQ.1 FOR IDOF.EQ.2, IDT(2,N)=',I10/
     3     ' IDT(1,N) MUST BE.EQ.0 FOR IDOF.EQ.-1'/
     3     ' IDT(2,N) MUST BE.EQ.0 FOR IDOF.EQ.-2'/     
     3     ' S T O P')
           STOP
         ENDIF
         GO TO 200
      ELSE IF(ISMEAREDG.EQ.5) THEN
         IND=1
         IF(ICONVG.EQ.1) THEN
            IF(IDOF.GT.0) THEN
               IF((IDOF.EQ.1.AND.IDT(1,N).EQ.0).OR.
     2            (IDOF.EQ.2.AND.IDT(2,N).EQ.0).OR.
     3            (IDOF.EQ.3.AND.IDT(3,N).EQ.0).OR.
     1            (IDOF.EQ.4.AND.IDT(4,N).EQ.0))IND=0
            ELSE
               IF((IDOF.EQ.-1.AND.IDT(1,N).GT.0).OR. 
     1            (IDOF.EQ.-2.AND.IDT(2,N).GT.0))IND=0          
            ENDIF
         ELSE
            IF((IDOF.EQ.3.AND.IDT(3,N).EQ.0).OR.
     1         (IDOF.EQ.4.AND.IDT(4,N).EQ.0))IND=0 
         ENDIF   
         IF(IND.EQ.1) THEN
            
            IF(IDOF.GT.0) THEN
               WRITE(3,6107)N,IDOF,IDT(IDOF,N)
 6107          FORMAT(/' INPUT ERROR FOR PRESCRIBED VALUE AT NODE:'/
     1         ' NODE NUMBER                                 N=',I10/
     2         ' DEGREE OF FREEDOM,                       IDOF=',I10/
     3         ' INPUT FOR IDT(1-4,N):'/
     3         ' IDT(IDOF,N) MUST BE.EQ.0          IDT(IDOF,N)=',I10/ 
     3         ' S T O P')
            ELSE
               IDOFA=-IDOF
               WRITE(3,6117)N,IDOF,IDT(IDOFA,N)
 6117          FORMAT(/' INPUT ERROR FOR PRESCRIBED VALUE AT NODE:'/
     1        ' NODE NUMBER                                 N=',I10/
     2        ' DEGREE OF FREEDOM,                       IDOF=',I10/
     3        ' INPUT FOR IDT(1-4,N):'/
     3        ' IDT(-IDOF,N) MUST BE.GT.0            IDT(1,N)=',I10/ 
     3        ' S T O P')
            ENDIF
            STOP  
         ENDIF        
      ENDIF         
  200 CONTINUE    
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(3,5020) N,IDOF,NC,FAK
      NZADJ(1,I) = N
      NZADJ(2,I) =IDOF
      NZADJ(3,I) = NC
      ZADFM(I) = FAK
      IF(NAUT.GT.0) GO TO 30
      IF(KORC.GT.0) GO TO 20
      IF(I.EQ.NZADP) GO TO 50
      GO TO 5
   20 NAUT=1
      N1=N
      FAK1 = FAK
      KORA=KORC
      GO TO 5
CE    AUTOMATIC GENERATED DATA BETWEEN NODES N1 AND N2
   30 N2=N
      FAK2 = FAK
      RKORA = KORA
      RN1N2 = N2-N1
      DD = RKORA/RN1N2
      DF1F2 = DD*(FAK2 - FAK1)
      IAUT=(N2-N1)/KORA-1
      N = N1
      DO 34 J=1,IAUT
      I=I+1
      N=N+KORA
      NZADJ(1,I) = N
      NZADJ(2,I) = IDOF
      NZADJ(3,I) = NC
      ZADFM(I) = ZADFM(I-1)+DF1F2
   34 CONTINUE
      IF(I.EQ.NZADP) GO TO 50
      NAUT=0
      IF(KORC.EQ.0) GO TO 5
      KORA=KORC
      NAUT=1
      N1=N2
      FAK1=FAK2
      GO TO 5
CE    PRINT READ AND GENERATED NODAL DATA
   50 IPENALTY=0
      IF(NULAZ.NE.1.AND.NULAZ.NE.3.OR.IPENALTY.EQ.0) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(3,2040)
      IF(ISRPS.EQ.1)
     1WRITE(3,6040)
      DO 70 I=1,NZADP
      NJ=NZADJ(1,I)
      IDOF=NZADJ(2,I)
      IF(NJ.EQ.0) GO TO 70
      NC = NZADJ(3,I)
      FAK = ZADFM(I)
      WRITE(3,5050) NJ,IDOF,NC,FAK
   70 CONTINUE
      RETURN
 1000 FORMAT(I10,2I5,F10.2,I5)   
 5020 FORMAT(5X,I10,I5,2X,I5,7X,1PD12.5,I8) 
 5050 FORMAT(11X,I10,12X,I5,8X,I5,15X,1PD15.8)
 2000 FORMAT(///'1'/6X,'P O D A C I    O    Z A D A T I M    T E M P E R
     1 A T U R A M A'/6X,62('-'))
 2010 FORMAT(////6X,'U C I T A N I     P O D A C I'///11X,
     1'CVOR BR.    FUNK.BR.   FAKTOR MNOZ.   KORAK GEN.')
 2030 FORMAT(' ZA CVOR N =',I5,' POSTOJI OGRANICENJE,',
     1' UCITANA TEMPERATURA SE IGNORISE')
 2040 FORMAT(6X,
     1'G E N E R I S A N I     P O D A C I'///11X,'JEDNACINA BROJ',
     15X,'BROJ VREM. FUNKCIJE',5X,'FAKTOR MNOZENJA')
 6000 FORMAT(///'1'/6X,'D A T A    F O R   P R E S C R I B E D    T E M
     1P E R A T U R E S'/' C O N C E N T R A T I O N S',6X,65('-'))
 6010 FORMAT(/6X,'I N P U T    D A T A'/11X,
     1'NODE  IDOF FUNCTION NO  PROP FACTOR')
 6030 FORMAT(' ZA CVOR N =',I5,' POSTOJI OGRANICENJE,',
     1' UCITANA TEMPERATURA SE IGNORISE')
 6040 FORMAT(6X,
     1'G E N E R A T E D     D A T A'///11X,'EQUATIONS NO.',
     15X,'TIME FUNCTION NO.   ',5X,'MULTIPL. FACTOR')
      END
C =================================================================
C
C=====================================================================      
      SUBROUTINE ULZADVS(NZADNODVS,NZADFVS,ZADVALUEVS,IDVS,MULTIMOL,
     1           NDIMCELL,NDIMVS,NVSMAX,NPT,NTABFT,NZADVS) 
C
C     READ PRESCRIBED VALUES FOR CELL VESICLES
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*130 ACOZ
      DIMENSION NZADNODVS(4,NZADVS),NZADFVS(NZADVS),ZADVALUEVS(NZADVS),
     1          IDVS(MULTIMOL,NDIMVS,NVSMAX,NPT)
      WRITE(3,2000)NZADVS
 2000 FORMAT(/' PRESCRIBED VALUES FOR VESICLES WITHIN CELLS'/
     1        ' NUMBER OF PRESCRIBED VALUES, NZADVS=',I5)
      WRITE(3,2001)
 2001 FORMAT('      NODE  MOLEC  CELLL VESICLE  TIME FUNCT COEFFICIENT')
      DO NN=1,NZADVS
          CALL ISPITA(ACOZ) 
          READ(ACOZ,1000)NODE,MM,NCELL,NVS,NF,COEFF
 1000     FORMAT(I10,4I5,F10.2)
          INDERR=0
          IF(NODE.GT.NPT.OR.MM.GT.MULTIMOL.OR.NCELL.GT.NDIMCELL.
     1       OR.NVS.GT.5.OR.NF.GT.NTABFT)INDERR=1
          IF(IDVS(MM,NCELL,NVS,NODE).NE.0)INDERR=1 
          IF(INDERR.EQ.1) THEN
             WRITE(3,2010)NODE,MM,NCELL,NVS,NF
 2010        FORMAT(/' ERROR IN INPUT DATA FOR PRESCRIBED VALUES OF',
     1               ' VESICLES WITHIN CELLS'/' INPUT VALUES:'/
     2               ' NODE (LE. THAN NPT)            =',I7/
     3               ' MOLECULE NO(LE. THAN MULTIMOL) =',I7/
     4               ' CELL NUMBER(LE. THAN NDIMCELL) =',I7/
     5               ' VESICLE NO(LE. THAN 5)         =',I7/
     6               ' FUNCTION NO(LE. THAN NTABFT)   =',I7/ 
     7               ' DOF FOR VESICLE MUST BE RESTRAINED'/
     6               ' STOP')
             STOP
          ENDIF
          NZADNODVS(1,NN)=NODE
          NZADNODVS(2,NN)=MM
          NZADNODVS(3,NN)=NCELL
          NZADNODVS(4,NN)=NVS
          NZADFVS(NN)=NF
          ZADVALUEVS(NN)=COEFF
          WRITE(3,2020)NODE,MM,NCELL,NVS,NF,COEFF
 2020     FORMAT(I10,2I7,I8,I7,3X,F10.2)
      ENDDO
      
      RETURN
      END  
C =================================================================
C
C===================================================================== 
      SUBROUTINE UZADATLY(INPUTLY,ZADVALUELY,NZADPV,NZADLYMPH,
     1           NZADLYMPHMX,MULTIMOL,MM,KKK)    
C
C     READ PRESCRIBED VALUES FOR LYMPHATIC SYSTEM
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*130 ACOZ
      DIMENSION INPUTLY(MULTIMOL,3,NZADLYMPHMX),
     1          ZADVALUELY(MULTIMOL,NZADLYMPHMX),NZADPV(10)
      IF(MULTIMOL.EQ.1) THEN
          MMM=1
          WRITE(3,2000)NZADLYMPH
 2000     FORMAT(/' PRESCRIBED VALUES FOR LYMPHATIC SYSTEM'/
     1        ' NUMBER OF PRESCRIBED VALUES, NZADLYMPH=',I5/
     2   5X,'NODE',2X,'NDOF (1-PESS,2-CONC)',2X,'TIME F',3X,
     3   'MULT COEFF' )
      ELSE
          NZADLYMPH=NZADPV(MM)
          IF(NZADLYMPH.EQ.0) RETURN
          KKK=KKK+1
          IF(KKK.EQ.1)WRITE(3,2001)NZADLYMPHMX
 2001     FORMAT(/' PRESCRIBED VALUES FOR LYMPHATIC SYSTEM'/
     1        ' MAX NUMBER OF PRESCRIBED VALUES FOR A MODEL,',
     2        ' NZADLYMPHMX=',I5/
     3        5X,'NODE',2X,'NDOF (1-PESS,2-CONC)',2X,'TIME F',3X,
     3   'MULT COEFF' )      
          CALL ISPITA(ACOZ) 
          READ(ACOZ,1000)MMM
          NZADLYMPH=NZADPV(MMM)
          
 1000     FORMAT(I5)      
          WRITE(3,2002)MMM,NZADLYMPH
 2002     FORMAT(' MOLECULE NUMMBER MM=',I5,5X,' NZADLYMPH=',I5)        
      ENDIF
      DO NN=1,NZADLYMPH
          CALL ISPITA(ACOZ) 
          READ(ACOZ,1001)NODE,NDOF,NF,COEFF
 1001     FORMAT(I10,2I5,F10.2)
          INPUTLY(MMM,1,NN)=NODE
          INPUTLY(MMM,2,NN)=NDOF
          INPUTLY(MMM,3,NN)=NF
          ZADVALUELY(MMM,NN)=COEFF
          WRITE(3,2020)NODE,NDOF,NF,COEFF
 2020     FORMAT(I10,I7,12X,I7,3X,F10.2)
      ENDDO 
      RETURN
      END 
C =============================================================
C
C==============================================================