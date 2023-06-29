C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZT1DARCY(CORD,IDT,IDTV,IDTMIG,IDARCYMA,IPIPENODV,
     1 IACTIVE,VOLFRACTV,DIAMSMEV,THICKSMEV,VOLSMEV,WALLAREAV,
     2 DARCYCOEFV,AMIPV,AMICV,HYDRCOEFPV,HYDRCOEFCV,DARTENSOR,IDTVCELL,
     3 VFMULTIT,DIAMCELLV,IVMIGV,VMIGV,DPIPEEXV,DCELLSV,DPIPECAP,
     4 IPARTNODEV,PARTNODEV,IACTIVENODEV,ACTIVENODEV,ISOURCENODEV,
     5 SOURCENODEV,IADHESINODEV,ADHESINODEV,DWALLV,DWALLSMEV,DWALLCELLV,
     6 IDVS,IDLYMPHH,IDLYMPHC,RVLYMPH,HCOEFLYNOD,HCOEFLYNODW,CCOEFLYNOD,
     7 CCOEFLYNODW,DCOEFLYNOD,NPT,IPART,ISMEARED,IPARTSME,NZADVOL,JEDN,
     8 NZADP,ICONVG,ILARGEVES,NDIMIDG) 
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)   
      CHARACTER*130 ACOZ
      COMMON /REPEAT/ LCORD,LTABK,LTABC,LNTAKV,LNTACV,LIDTE,LIDT,LTCONS,
     1LMAX,LMAXA,LMAXM,LMHT,LNTAQE,LQEFN,LNTAQP,LQPFN,LNTAOK,LTOKFN,
     2LNTAHP,LHPFN,LNTAHR,LHRFN,LTT1,LTT2,LTT0,LSK,LCONM,LF,LFC,LNEL,
     3LIBFK,LFAKP,LISKC,LGUSM,LTOPM,LTMNM,LFOHR,LFVPOT
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
	common/diffusmain/indiffus,NDS,maxdt,NPRINTnds,lvreme,leprdt,
	1                  NUMMAT,maxPmat,LMATER,LMATFUN,LAREAD,LFLUXDX,
     1                  LFLUXDY 
      COMMON/DARCYEQ/ ELDARCYX,ELDARCYY,ELDARCYZ,VDARCY,PINLET,
     1       DARCYX,JEDNDARCY,ISOLVER,IDARCY,LDARCYD,LDARCYMI,
     2       LIDARCYMA,NPTINDARCY
 	COMMON/DARCYCONV/ICONV,IDINF,GUSTF,IDIFUSE,LMATCONV,LMATDIF,
     1       NMATTCONV,NMATTDIF,maxPmatC,maxPmatD,LMATFUNC,LMATFUND,
     2       LNELFICT,MAXCONNP,MAXCONNP2
      COMMON/DARCYHEAT/IHEAT,NMATTHEAT,maxPmatH,LMATHEAT,LMATFUNHRO,
     1      LMATFUNHSH,LMATFUNHK
      COMMON/CSMEARED2/IPASSME,JEDNPRES,JEDNCONC,LTT1P,LTT1C,LTT10P,
     1      LTT10C,LDIRCOEFF,LDPIPECA,LDARCYCA,I2D3D,IANIZS,NDIMTENSOR
      COMMON/CMULTITIS/MULTITIS,MULTIMOL,IACTIVECELL(4),NDIMCELL,LIDTV,
     1      LIDTVCELL,JEDNPV(10),JEDNCV(10),JEDNPMAX,JEDNCMAX,MM
      COMMON/MULTITISAD/LIPARTITV,LIACTIVEV,LIVSOURCEV,LIADHESIV,
     1 LVOLFRACT,LDARCYCOEFV,LHYDRCOEFPV,LHYDRCOEFCV,LAMIPV,LAMICV,
     2 LVOLSMEV,LDARTENSOR,LWALLAREAV,LVFMULTIT,LDIAMCELLV,
     3 LDPIPEEXV,LDCELLSV,LDPIPECAP,LIPARTNODEV,LPARTNODEV,
     4 LIACTIVENODEV,LACTIVENODEV,LISOURCENODEV,LSOURCENODEV,
     5 LADHESINODEV,LIADHESINODEV,LDWALLSMEV,LDWALLCELL,LDIAMSMEV,
     6 LTHICKSMEV,LDWALLV,LIPIPENODV,LICRITICALV,LCRITICALV,
     7 LTT1PM,LTT01PM,LTT1CM,LTT01CM,LFMULT  
      COMMON/CMATDATA/NDIMPIPEMAX,NDIMCAPMAX,NDIMCELLMAX,NDIMPMAX,
     4 NDIMAMAX,NDIMSMAX,NDIMADHMAX,NDIMWPIPEMAX,NDIMWCAPMAX,
     5 NDIMWCELLMAX,DARCYCOEF,AMICOEFP,AMICOEFC,HYDROCOEFP,HYDROCOEFC,
     6 DARCYTENSOR(6),DPIPE,DTISSUE,DIFFTENSOR(6),PARTPIPE1,PARTPIPE2,
     7 PARTCELL(2,4),ACTIVECOEF(2,4),SOURCECOEFF(4),ADHCOEFF(2,4),
     8 DWALLPIPE,DWALLCAP,DWALLCELL(4),INODEDTOT,INODEWTOT,LDPIPEEXA,
     9 LDPIPECAPA,LCELLA,LPARTA,LACTIVEA,LSOURCEA,LADHISEA,LWALLPIPEA,
     1 LWALLCAPA,LWALLCELLA,IREADDIFNODE
      COMMON/MATERCONST/ICMTABLE(10),IDMTABLE(10,7),IWMTABLE(10,12),
     1 NDIMPIPE(10),NDIMCAP(10),NDIMCELLV(10),NDIMPART(10),NDIMACT(10),
     2 NDIMSOURCE(10),NDIMADH(10),NDIMWPIPE(10),NDIMWCAP(10),
     3 NDIMWCELL(10),IDPIPEEXC(2,10),DPIPEEXC(2,10),IDPIPECAPC(10),
     1 DPIPECAPC(6,10),ICELLC(10,10),CELLC(10,10),IPARTC(12,10),
     2 PARTC(12,10),IACTIVEC(10,10),ACTIVEC(2,10,10),ISOURCEC(10,10),
     3 SOURCEC(10,10),IADHISEC(12,10),ADHISEC(2,12,10),IWALLPIPEC(10),
     4 WALLPIPEC(10),IWALLCAPC(10),WALLCAPC(10),IWALLCELLC(10,10),
     5 WALLCELLC(10,10),ICRITICALV(10,10),CRITICALV(10,10),
     6 WEIGHTMOLV(10)
      COMMON/CZADATMULT/NZADPV(10),NZADPMAX,LZADATM0,LZADATM1,LNZADNODF,
     1      LZADVALUE,NPTIMESYSV(10),NPTSYSMAX,LZADSYSV
      COMMON/CMIGRATION/IMIGTABLE(10),JEDNMIGV(10),VMIGC(10),
     1 INDCHTAXV(10),INDCHTAX,IMIGTOT,IMIGNODE,LIVMIGV,LVMIGV,LIDTMIG,
     2 INDMIG,LVFTT,JEDNMIG,INDCHTAXTOT,LCHTAXV
      COMMON/CCELLIN/INDCELLIN,INDCELLINV(10),VFVS(10,5),DIAMVS(10,5),
     1 RAVVS(10,5),DDIFFICY(10,10),DDIFFIVS(10,10,5),DDIFFICM(10,10,5)
     2 ,IPARTIVS(10,10,5),PARTIVS(10,10,5),NDIMVS,NVSMAX,LIDVS,NZADVS,
     3 LNZADNODVS,LNZADNFVS,LNZADFVS,LZADVALUEVS,LZADATAVS0,LZADATAVS
      COMMON/CLYMPH/INDLYMPH,INDLYMPHH,INDLYMPHC,LRVLYMPH,LPLYMHV,
     1 LCLYMHV,TABLYHW(10),TABLYCW(10,10),TABLYD(10),TENSORLYP(10,6),
     2 TENSORLYC(10,10,6),ITABLYP,ITABLYC,ITABLYD,LHCOEFLYNOD,
     3 LHCOEFLYNODW,LCCOEFLYNOD,LCCOEFLYNODW,LDCOEFLYNOD,LIDLYMPHH,
     4 LIDLYMPHC
      COMMON/CLYMHZADA/NZADLYMPH,LINPUTLY,LZADVALUELY,LZADALYMH1,
     1 LZADALYMH0,NZADPLYV(10),NZADLYMPHMX     
C    
      DIMENSION CORD(NPT,3),IDT(NDIMIDG,NPT),IDTV(MULTIMOL,NDIMIDG,NPT),
     1          IDTMIG(IMIGTOT,NPT),IDARCYMA(NPT),IPIPENODV(NPT),
     2          IACTIVE(4)  
      DIMENSION SME(5),IDTLOC(10),IDTLOCI(10),DXYZ(3),NODAUX(20),
     1DARCYCL(6),DIFFCAL(6),IDMULLOC(10),CELLOC1(12),CELLOC2(12),
     2   DIFLOC(12),IDTLOCIV(20)
     	DIMENSION VOLFRACTV(NPT),DIAMSMEV(NPT),THICKSMEV(NPT),
     2 VOLSMEV(NPT),WALLAREAV(NPT),DARCYCOEFV(NPT),AMIPV(NPT),
     3 AMICV(NPT),HYDRCOEFPV(NPT),HYDRCOEFCV(NPT),
     4 DARTENSOR(NDIMTENSOR,NPT),IDTVCELL(MULTIMOL,NDIMCELL,NPT),
     5 VFMULTIT(NDIMCELL,NPT),DIAMCELLV(NDIMCELL,NPT),
     6 IVMIGV(IMIGNODE,NPT),VMIGV(IMIGNODE,NPT), 
     6 DPIPEEXV(MULTIMOL,NPT,2),DCELLSV(MULTIMOL,NPT,NDIMCELLMAX),
     5 DPIPECAP(MULTIMOL,NDIMTENSOR,NPT),
     6 IPARTNODEV(MULTIMOL,NDIMPMAX,NPT),
     7 PARTNODEV(MULTIMOL,NDIMPMAX,NPT),
     8 IACTIVENODEV(MULTIMOL,NDIMAMAX,NPT),
     9 ACTIVENODEV(MULTIMOL,NDIMAMAX,NPT),
     7 ISOURCENODEV(MULTIMOL,NDIMSMAX,NPT),
     8 SOURCENODEV(MULTIMOL,NDIMSMAX,NPT),
     9 IADHESINODEV(MULTIMOL,NDIMADHMAX,NPT),   
     1 ADHESINODEV(MULTIMOL,NDIMADHMAX,NPT),DWALLV(MULTIMOL,NPT),
     9 DWALLSMEV(MULTIMOL,NPT),DWALLCELLV(MULTIMOL,NDIMWCELLMAX,NPT),
     1 IDVS(MULTIMOL,NDIMVS,NVSMAX,NPT),
     2 IDLYMPHH(NPT),IDLYMPHC(MULTIMOL,NPT),RVLYMPH(NPT),
     3 HCOEFLYNOD(NDIMTENSOR,NPT),HCOEFLYNODW(NPT),
     4 CCOEFLYNOD(MULTIMOL,NDIMTENSOR,NPT),CCOEFLYNODW(MULTIMOL,NPT),
     5 DCOEFLYNOD(NPT)  
C
C     READ NODES AT INTERFACE BETWEEN INLET MICROCHANNELS AND MACROCHANNEL
C     FOR NDS1-2 DEVICE OR 
C	READ NODES AT INTERFACE BETWEEN INLET RESERVOIR	AND NANOCHANNEL
C	FOR NDS-MD WITH RESERVOIRS 
C
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
         IF(ISMEARED.EQ.0) THEN
             WRITE(3,6010)
         ELSE
             WRITE(3,6110)
         ENDIF  
      ENDIF   
C WITH SMEARED MODEL     
 6010 FORMAT(/6X,'I N P U T   N O D A L   D A T A'/
     16X,31('-')/
     16X,' NODE',2X,'IND.PRE.CONC.',2X,'IND.PARTITION',14X,
     1'C O O R D I N A T E S',13X,'PARTITION COEFF',3X,
     2'IND.AUT.GEN.'/19X,'(0-NO, 1-YES)',13X,
     3'X',17X,'Y',17X,'Z',32X,'(0-NO,>0-YES)'/
     46X, ' NOTE ABOUT INPUT FOR MATRIX IDT(NDIMIDG,NPT):'/
     26X,'  1 - PRESCRIBED VALUES - PRESSURE,CONCENTRATION, OR HEAT'/
     3 6X,' -1 - PRESCRIBED VELOCITY (POSITIVE IN OUTWARD NORMAL)'/
     4 6X,' -2 - RESTRAINED FLUX'/
     56X,' FOR IDARCY=-2 (DIFF OR HEAT HOMOGENIZAT), ALL DOF ARE FREE')
     
 6110 FORMAT(//6X,'I N P U T   N O D A L   D A T A'/
     16X,31('-')/
     16X,'NOTE: LINES FOR MATERIAL DATA ARE NOT INPUT WHEN THE FLAG IN',
     1 ' MATERIAL TABLE IS NOT EQUAL TO 3'/' LINE NO'/
     13X,'1',' NODE',4X,'INDICATOR FOR PRESCR VALUES 1-4',
     1   4X,'C O O R D I N A T E S',5X,'PIPE FLAG',3X,
     1   'ADDIT CONC MULTIMODEL',1X,'INPUT',
     1   ' ALWAYS'/
     2    13X,'1,2-pressures; 3,4-conc',
     3    8X,'X',12X,'Y',12X,'Z'/
     2   11X,'  1 - prescribed; ',
     3   ' -1 - prescrib velocity (+ outward normal)')

       WRITE(3,6011)
 6011 FORMAT(3X,'2',6X,'CAPIL VOL FRACT',2X,'CAPIL DIAM',2X,'CAP WALL',
     1 ' THICK',
     1  2X,'NODAL VOLUME',4X,'NODAL AREA',27X,
     2   'INPUT ALWAYS')
       WRITE(3,6012)
 6012 FORMAT(3X,'3',6X,'DARCY C-TISS',2X 'VISCOS-PIPE',2X,
     1  2X 'VISCOS-CAP',2X,'H COEFF-PIPE',4X,'H COEFF-CPAILL',
     2  17X,'INPUT FOR ICONV.EQ.1')
       WRITE(3,6040)
 6040  FORMAT(3X,'4',6X,'DARCY TENSOR (MAX 6 COMPONENTS)',57X,
     1 'INCONV.EQ.1.AND.ICMTABLE(6).EQ.3')
       IF(MULTITIS.GT.0)
     1 WRITE(3,6013)
 6013  FORMAT(3X,'5',24X,'VOLUMETRIC FRACTIONS(1-NDIMCELL)',
     2 ' AND CELL DIAMETERS (1-NDIMCELL)',6X,
     2 'INPUT FOR COMPOSITE TISSUE')
       IF(IMIGNODE.GT.0) 
     1 WRITE(3,6014)
 6014  FORMAT (10X,'MIGRATION FLAGS(1-YES,0-NO)',5X,' VELOCITIES',
     1 45X,'INPUT WHEN IMIGTOT.GT.0')
       IF(IREADDIFNODE.EQ.0.AND.INDIFFUS.EQ.1) THEN
         WRITE(3,6029)
 6029    FORMAT(' NOTE: ALL DIFFUSION MATERIAL PARAMETERS ARE',
     1    ' INDEPENDENT OF NODES'/
     2           '       GIVEN AS CONSTANT OR CONCENTRATION DEPENDENT',
     3     ' AND ARE NOT INPUT AT NODAL POINTS')
       ENDIF  
       WRITE(3,6033)
 6033  FORMAT(4X,'---DIFFUSION DATA: INPUT WHEN',65X,
     1  'INDIFFUS.EQ.1.AND.IREADDIFNODE.EQ.1'/
     13X,'6',3X,' MOLECULE NO'/28X,'DIFF COEFF PIPE AND EXTRACELL',
     1  7X, 'DIFFUSION TENSOR',18X, 
     1  'INPUT WHEN INDIFFUS.EQ.1')
       WRITE(3,6034)
 6034  FORMAT(3X,'7',24X,'DIFF COEFFICENTS WITHIN CELLS',41X,
     1 'INPUT WHEN NDIMCELLMAX.GT.0')
       WRITE(3,6035)
 6035  FORMAT(3X,'8-11',9X,'INDICATORS',2X,'COEFFICIENTS:PARTIT/',
     1  'ACTIVE/','SOURCE/ADHESION')
       WRITE(3,6036)
 6036  FORMAT(3X,'12',11X,'DWALL PIPE',2X,'DWALLCAP',3X,
     1 'DIFF COEFFTS WALL CELL MEMRANES', 28X,'INPUT WHEN FLAG=3',
     2 ' FOR ANY WALL')
   
      IF(ICONV.EQ.1) CALL ICLEAR(IDARCYMA,NPT)
      CALL ICLEAR(IPIPENODV,NPT)
      IF(IMIGTOT.GT.0) THEN
          CALL ICLEAR(IDTMIG,IMIGTOT*NPT)
          IF(IMIGNODE.GT.0) THEN
             CALL ICLEAR(IVMIGV,IMIGNODE*NPT)
             CALL CLEAR(VMIGV,IMIGNODE*NPT)       
          ENDIF
      ENDIF    
      CALL CLEAR(DARCYCL,6)
      CALL CLEAR(DIFFCAL,6)
      NTCONS=0
      ILARGEVES=0
      I = 0
      NAUT=0
      MULTIMODA=2*(MULTIMOL-1)
C      
C     ------------------------- LOOP OVER NODES-----------------------      
C
      DO 400 NODE=1,NPT
C     LINE 1 ---------   
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) N,(IDTLOCI(J),J=1,4),(CORD(N,J),J=1,3),IPIPENOD,
     1      (IDTLOCIV(J),J=1,MULTIMODA)
 1000 FORMAT(I10,4I5,3e15.6,I5,20I2)
      IPIPENODV(N) = IPIPENOD
      DO J=1,4
         IF(IACTIVE(J).EQ.0) IDTLOCI(J)=1
      ENDDO
      IF(MULTIMOL.EQ.1) THEN
         DO J=1,4
            IDT(J,N)=IDTLOCI(J)
         ENDDO   
      ELSE
         DO MM=1,MULTIMOL
            IDTV(MM,1,N)=IDTLOCI(1)
            IDTV(MM,2,N)=IDTLOCI(2)
         ENDDO
         IDTV(1,3,N)=IDTLOCI(3)
         IDTV(1,4,N)=IDTLOCI(4)
         KK=1
         DO MM=2,MULTIMOL
            IF(IACTIVE(3).EQ.0) IDTLOCIV(KK)=1
            IF(IACTIVE(4).EQ.0) IDTLOCIV(KK+1)=1
            
            IDTV(MM,3,N)=IDTLOCIV(KK)
            IDTV(MM,4,N)=IDTLOCIV(KK+1)    
            KK=KK+2
         ENDDO
      ENDIF
      IF(MULTIMOL.EQ.1) THEN
         WRITE(3,5697)N,(IDTLOCI(J),J=1,4),(CORD(N,J),J=1,3),IPIPENOD
 5697    FORMAT(/I10,4X,4I5,6X,3E15.6,I5)
      ELSE
         WRITE(3,5698)N,(IDTLOCI(J),J=1,4),(CORD(N,J),J=1,3),IPIPENOD,
     1         (IDTLOCIV(J),J=1,MULTIMODA)
 5698    FORMAT(/I10,4I5,6X,3E15.6,I2,1X,10(2I2,2X))
      ENDIF           
C     LINE 2-----------------
      CALL ISPITA(ACOZ)
      READ(ACOZ,1010)VOLFRACTV(N),DIAMSMEV(N),THICKSMEV(N),VOLSME,
     1    WALLAREACC     
C     LINE 3 IN CASE OF ICONV.EQ.1 --------------
      DARCYC=0.D0
      AMIP=0.D0
      AMIC=0.D0
      HYDRCOEFP=0.D0
      HYDRCOEFC=0.D0
      IF(ICONV.EQ.1) THEN
         CALL ISPITA(ACOZ)
         READ(ACOZ,1010)DARCYC,AMIP,AMIC,HYDRCOEFP,HYDRCOEFC
 1010    FORMAT(7X,8E15.6)          
         IF(ICMTABLE(1).EQ.3)DARCYCOEFV(N)=DARCYC
         IF(ICMTABLE(2).EQ.3)AMIPV(N)=AMIP 
         IF(ICMTABLE(3).EQ.3)AMICV(N)=AMIC
         IF(ICMTABLE(4).EQ.3)HYDRCOEFPV(N)=HYDRCOEFP
         IF(ICMTABLE(5).EQ.3)HYDRCOEFCV(N)=HYDRCOEFC         
C     LINE 3A IN CASE OF ICMTABLE(6).EQ.3 NODAL VALUES FOR DARCY TENSOR --------------
         IF(ICMTABLE(6).EQ.3) THEN
            CALL ISPITA(ACOZ)
            READ(ACOZ,1010)(DARCYCL(K),K=1,6)
            DO K=1,NDIMTENSOR
               DARTENSOR(K,N)=DARCYCL(K)
            ENDDO 
         ENDIF  
      ENDIF
      IF(ISMEARED.EQ.1.OR.MULTITIS.GT.0)THEN
         VOLSMEV(N)=(1.-VOLFRACTV(N))*VOLSME
      ENDIF
      IF(WALLAREACC.LT.1E-10) THEN
         IF(DIAMSMEV(N).GT.1.E-10)  WALLAREACC=4./DIAMSMEV(N)
      ENDIF
      IF(IPIPENOD.EQ.0) THEN
         WALLAREAV(N)=VOLFRACTV(N)*WALLAREACC*VOLSME
      ELSE
        IF(WALLAREACC.LT.1E-10) WALLAREACC=1.D0
        WALLAREAV(N)= WALLAREACC*VOLSME
      ENDIF
      WRITE(3,5699)VOLFRACTV(N),DIAMSMEV(N),THICKSMEV(N),VOLSME,
     1    WALLAREAV(N),DARCYC,AMIP,AMIC,HYDRCOEFP,HYDRCOEFC
 5699 FORMAT(6X,5E15.6/6X,5E15.6)
      IF(ICONV.EQ.1.AND.ICMTABLE(6).EQ.3)
     1 WRITE(3,5699)(DARTENSOR(K,N),K=1,NDIMTENSOR)
C
c     LINE 4 - BASIC DATA FOR COMPOSITE-MULTITISSUE, VOLUME FRACTIONS
C
      IF(MULTITIS.GT.0) THEN
         CALL ISPITA(ACOZ)
         READ(ACOZ,1020)(CELLOC1(I),I=1,10)   
         CALL ISPITA(ACOZ)
         READ(ACOZ,1020)(CELLOC2(I),I=1,10)
         VFTOT=0.D0        
         DO I=1,NDIMCELL
            VFMULTIT(I,N)=CELLOC1(I)
            IF(CELLOC1(I).GT.0.999D0) THEN
                WRITE(3,2777)N,I,CELLOC1(I)
 2777           FORMAT(/' INPUT ERROR - VOLUME FRACTION OF CELL TYPE',
     1          ' MUST BE LT.0.999'/
     2          ' NODE NUMBER.....   ',I7/
     3          ' CELL TYPE NUMBER...',I7/
     4          ' VOLUME FRACTION... ',D12.3/
     5          ' S T O P')
                STOP
            ENDIF
            VFTOT=VFTOT+CELLOC1(I)
            DIAMCELLV(I,N)=CELLOC2(I)    
         ENDDO  
C
C        CHECK FOR VOLUME FRACTIONS - SUMM MUST BE EQUAL TO 1.
         IF(NDIMCELL.GT.1.AND.VFTOT.GT.0.999D0) THEN
            WRITE(3,2778)N
 2778           FORMAT(/' INPUT ERROR - VOLUME FRACTION OF CELLS',
     1          ' SUMM OF VOLUME FRACTIONS MUST BE LT.0.999'/
     2          ' NODE NUMBER.....   ',I7/
     5          ' S T O P')
            STOP                    
         ENDIF          
 1020    FORMAT(6X,12F10.4)
         WRITE(3,5702)(CELLOC1(I),I=1,10),(CELLOC2(I),I=1,10)
 5702    FORMAT(20X,10F10.4/20X,10F10.4)
 5706    FORMAT((12X,10E12.3))
C
         IF(IMIGNODE.GT.0) THEN
           CALL ISPITA(ACOZ)
           READ(ACOZ,1027)(IDTLOC(I),I=1,10)
 1027      FORMAT(I13,9I3)
           CALL ISPITA(ACOZ)
           READ(ACOZ,1020)(CELLOC1(I),I=1,10) 
           WRITE(3,2037)(IDTLOC(I),I=1,10),(CELLOC1(I),I=1,10)
 2037      FORMAT(I10,9I3/10X,10F10.4)
           KK=0
           DO I=1,NDIMCELL    
              IF(IMIGTABLE(I).EQ.3) THEN
                 KK=KK+1
                 IVMIGV(KK,N)=IDTLOC(I)
                 VMIGV(KK,N)=CELLOC1(I)
              ENDIF 
           ENDDO
        ENDIF  
      ENDIF                                                             
C
C     LINES FOR DIFFUSSION DATA
C     LOOP OVER MOLECULES
C
      IF(INDIFFUS.EQ.0) GO TO 440
      DO 420 MMI=1,MULTIMOL
C        MOLECULE NUMBER
         CALL ISPITA(ACOZ)
         READ(ACOZ,1023)MM, (IDTLOC(I),I=1,10),(IDTLOCIV(I),I=1,10)
 1023    FORMAT(I11,I4,9I3,3X,10I3)
C
         WRITE(3,5707)MM,(IDTLOC(I),I=1,10)
 5707    FORMAT(I12,I4,9I3/10X,10F10.3)
         IF(MULTITIS.GT.0) THEN       
             DO J=1,NDIMCELL
                IDTVCELL(MM,J,N)=IDTLOC(J)    
             ENDDO
             IF(INDCELLIN.GT.0) THEN
C               INTERNAL CELL DIFFUSION FOR EACH CELL WITH VESICLES ENTER AN IDVS LINE
                KIN=0
                DO J=1,NDIMCELL
                   NINI=INDCELLINV(J)
                   IF(NINI.GT.0) THEN
                      KIN=KIN+1
                      CALL ISPITA(ACOZ)
                      READ(ACOZ,5708)(IDTLOC(I),I=1,5)
 5708                 FORMAT(I15,4I3)
                      DO K=1,NINI
                        IDVS(MM,KIN,K,N)=IDTLOC(K)
                      ENDDO
                      WRITE(3,5709)(IDVS(MM,KIN,K,N),K=1,NINI)
 5709                 FORMAT(I16,4I3)
                   ENDIF
                ENDDO
             ENDIF                               
         ENDIF
        
        IF(IREADDIFNODE.EQ.0) GO TO 420
C
C        DIFFUSION COEFFICIENTS   
C
         CALL ISPITA(ACOZ)
         READ(ACOZ,1030)DDPIPE,DTISSUE,(DIFLOC(I),I=1,6)  
 1030    FORMAT(10X,12F10.2)  
         IF(NDIMPIPEMAX.GT.0) THEN
            DPIPEEXV(MM,N,1)=DDPIPE
            DPIPEEXV(MM,N,2)=DTISSUE
         ENDIF 
         IF(NDIMCAPMAX.GT.0) THEN
            DO K=1,NDIMTENSOR
               DPIPECAP(MM,K,N)=DIFLOC(K)
            ENDDO
         ENDIF  
         WRITE(3,5703)DPIPEEXV(MM,N,1),DPIPEEXV(MM,N,2),
     1               (DIFLOC(I),I=1,NDIMTENSOR)
5703     FORMAT(24X,2E12.3,12X,6E12.3)
         IF(MULTITIS.GT.0.AND.NDIMCELLMAX.GT.0) THEN
            CALL ISPITA(ACOZ)
            READ(ACOZ,1020)(DIFLOC(I),I=1,10)
            KK=0                  
            DO I=1,NDIMCELL    
                IF(ICELLC(I,MM).EQ.3) THEN
                   KK=KK+1              
                   DCELLSV(MM,N,KK)=DIFLOC(I)
                ENDIF
            ENDDO
            WRITE(3,5710)(DIFLOC(K),K=1,10)
5710        FORMAT(10X,10F10.4 )          
         ENDIF
C
C        PARTITIONING COEFFICIENTS  
C
         IF(NDIMPMAX.GT.0) THEN
            CALL ISPITA(ACOZ)
            READ(ACOZ,1021)NPARTN,(DIFLOC(I),I=1,12)
 1021       FORMAT(I7,3X,12F10.4)
            KK=0
            DO I=1,12 
              IF(IPARTC(I,MM).EQ.3)THEN
                 KK=KK+1
                 IPARTNODEV(MM,KK,N)=NPARTN
                 PARTNODEV(MM,KK,N)=DIFLOC(I)          
              ENDIF
            ENDDO 
            WRITE(3,5704)NPARTN,(DIFLOC(I),I=1,12)
 5704       FORMAT(I7,3X,12F10.4)
         ENDIF         
C
C        ACTIVE TRANSPORT COEFFICINETS  
C
         IF(NDIMAMAX.GT.0) THEN
            CALL ISPITA(ACOZ)
            READ(ACOZ,1021)NACTIV1,(DIFLOC(I),I=1,10)
            CALL ISPITA(ACOZ)
            READ(ACOZ,1022)(DARCYCL(I),I=1,10)
 1022       FORMAT(10X,12F10.2)
            KK=0 
            KK1=0                            
            DO I=1,10
               IF(IACTIVEC(I,MM).EQ.3)THEN
                   KK=KK+1
                   KK1=KK1+1
                   IACTIVENODEV(MM,KK1,N)=NACTIV1
                   ACTIVENODEV(MM,KK,N)=DIFLOC(I)
                   KK=KK+1
                   ACTIVENODEV(MM,KK,N)=DARCYCL(I)          
               ENDIF
            ENDDO 
            WRITE(3,5704)NACTIV1,(DIFLOC(I),I=1,10)
            WRITE(3,5714)(DARCYCL(I),I=1,10)
 5714       FORMAT(10X,12F10.4)   
         ENDIF
C
C        SOURCE TERMS                        
C
         IF(NDIMSMAX.GT.0) THEN
            CALL ISPITA(ACOZ)
            READ(ACOZ,1021)NSOURCE,(DIFLOC(I),I=1,10)
            KK=0
            DO I=1,10 
               IF(ISOURCEC(I,MM).EQ.3)THEN
                   KK=KK+1
                   ISOURCENODEV(MM,KK,N)=NSOURCE
                   SOURCENODEV(MM,KK,N)=DIFLOC(I)          
               ENDIF
            ENDDO
            WRITE(3,5704)NSOURCE,(DIFLOC(I),I=1,10) 
         ENDIF
C
C        ADHISION 
C
         IF(NDIMADHMAX.GT.0) THEN  
            CALL ISPITA(ACOZ)
            READ(ACOZ,1021)NADHIS1,(DIFLOC(I),I=1,10)
            CALL ISPITA(ACOZ)
            READ(ACOZ,1022)(DARCYCL(I),I=1,10)
            KK=0 
            KK1=0                             
            DO I=1,10       
               IF(IADHISEC(I,MM).EQ.3) THEN
                   KK=KK+1
                   KK1=KK1+1
                   IADHESINODEV(MM,KK1,N)=NADHIS1
                   ADHESINODEV(MM,KK,N)=DIFLOC(I)
                   KK=KK+1
                   ADHESINODEV(MM,KK,N)=DARCYCL(I)          
               ENDIF
            ENDDO 
            WRITE(3,5704)NADHIS1,(DIFLOC(I),I=1,10)
            WRITE(3,5714)(DARCYCL(I),I=1,10)
        ENDIF
C
C       WALL DIFFUSION DATA
C
        IF(NDIMWPIPEMAX.GT.0.OR.NDIMWCAPMAX.GT.0.OR.
     1     (MULTITIS.GT.0.AND.NDIMWCELLMAX.GT.0)) THEN
           CALL ISPITA(ACOZ)
           READ(ACOZ,1022)DWPIPE,DWALLSME,(DIFLOC(I),I=1,10) 
           IF(NDIMWPIPEMAX.GT.0) DWALLV(MM,N)=DWPIPE
           IF(NDIMWCAPMAX.GT.0)DWALLSMEV(MM,N)=DWALLSME   
           IF(MULTITIS.GT.0.AND.NDIMWCELLMAX.GT.0) THEN
              KK=0     
              DO I=1,NDIMCELL 
                 IF(IWALLCELLC(I,MM).EQ.3) THEN
                     KK=KK+1                      
                     DWALLCELLV(MM,KK,N)=DIFLOC(I)
                 ENDIF
              ENDDO
           ENDIF 
           WRITE(3,5714)DWPIPE,DWALLSME,(DIFLOC(I),I=1,10)
        ENDIF        
  420 CONTINUE

  440 CONTINUE  
C
      IF(INDLYMPH.EQ.0) GO TO 442 
C
C     READ DATA FOR LYMPHATIC SYSTEM
C
      IF(N.EQ.1) THEN
         WRITE(3,6100)
6100     FORMAT(' DATA FOR LYMPHATIC SYSTEM'/
     1   ' INDPRESC FOR CONV AND FOR ALL MOLECULES FOR DIFFUSION/',
     2   'VOLUME FRACTION, USED IF NOT IN THE TABLE:',
     2   ' DIAM,WALL HIDR COEF,DARCY TENSOR'/
     4   ' DIFFUSION: IF NOT IN THE TABLE: MOELEC NO/ WALL DIFF COEF,',
     5   'DIFF TENSOR') 
      ENDIF       
      CALL ISPITA(ACOZ)
      READ(ACOZ,1760)JH,RVLYMPH(N),DIAM,COEFWH,(DARCYCL(K),K=1,6),
     1 (IDTLOCIV(J),J=1,10)
      WRITE(3,6101)JH,(IDTLOCIV(J),J=1,MULTIMOL)
 6101 FORMAT(12X,I3,5X,10I2) 
      WRITE(3,6102)RVLYMPH(N),DIAM,COEFWH,(DARCYCL(K),K=1,NDIMTENSOR)
 6102 FORMAT(6X,9F10.4)
 1760 FORMAT(I5,9F10.4,10I2)

      IF(ISMEARED.EQ.1.OR.MULTITIS.GT.0)THEN
         VOLSMEV(N)=(1.-VOLFRACTV(N)-RVLYMPH(N))*VOLSME
      ENDIF 
      IF(INDLYMPHH.GT.0) THEN     
         IDLYMPHH(N)=JH
         IF(ITABLYP.EQ.0) THEN
            HCOEFLYNODW(N)=COEFWH
            DO K=1,NDIMTENSOR
               HCOEFLYNOD(K,N)=DARCYCL(K)
            ENDDO
         ENDIF
      ENDIF
      IF(ITABLYD.EQ.0) THEN
         DCOEFLYNOD(N)=DIAM
      ENDIF
      IF(INDLYMPHC.GT.0) THEN
         DO K=1,MULTIMOL
            IDLYMPHC(K,N)=IDTLOCIV(K)   
         ENDDO
         IF(ITABLYC.EQ.0) THEN
            DO J=1,MULTIMOL
               CALL ISPITA(ACOZ)
               READ(ACOZ,1761)MMM
 1761          FORMAT(I12)
               CALL ISPITA(ACOZ)
               READ(ACOZ,1762)COEFWC,(DARCYCL(K),K=1,6)
               WRITE(3,6103)MMM,COEFWC,(DARCYCL(K),K=1,NDIMTENSOR)
 6103          FORMAT(12X,I3,7F10.4)              
 1762          FORMAT(12X,7F10.4) 
               CCOEFLYNODW(MMM,N)=COEFWC 
               DO K=1,NDIMTENSOR
                  CCOEFLYNOD(MMM,K,N)=DARCYCL(K)
               ENDDO        
            ENDDO
         ENDIF
      ENDIF               
C     END FOR LYMPHATIC SYSTEM

  442 CONTINUE
C     END OF READING INPUT DATA  -----------------------       
 
 1011 FORMAT(7X,6E15.6,I3)   
C
      IF(ISMEARED.LE.1) THEN
        IDARCYMA(N) = IDTLOCI(1)
      ELSE
        IDARCYMA(N) = IDTLOCI(1)
      ENDIF
C
C    SET IDT MATRIX FOR ISEMAERD.EQ.1 AND MULTIMOL.EQ.1
C
      IF(ISMEARED.EQ.1.AND.MULTIMOL.EQ.1) THEN
         IF(ICONV.EQ.1) THEN
            IDT(1,N)=IDTLOCI(1)
            IF(INDIFFUS.EQ.1)IDT(2,N)=IDTLOCI(2)
            IDARCYMA(N)=IDTLOCI(1)
         ELSE
            IDT(1,N)=IDTLOCI(1)
            IDT(2,N)=IDTLOCI(2)
         ENDIF       
      ELSE
         DO K=1,4
            IDT(K,N)=IDTLOCI(K)
         ENDDO
         IF(ICONV.EQ.1.AND.IDTLOCI(1).LT.0)IDARCYMA(N)=-1
      ENDIF   
C      
  400 CONTINUE    
C
C     EQUATIONS FOR MIGRATION
C     
      JEDNMIG=0
      IF(MULTITIS.GT.0.AND.IMIGTOT.GT.0) THEN
          WRITE(3,2033)
 2033     FORMAT(/' EQUATIONS FOR  C E L L  M I G R A T I O N')         
          KID=0
          DO I=1,NDIMCELL 
              IF(IMIGTABLE(I).NE.0) THEN
                 KID=KID+1
                 DO N=1,NPT 
                    IF(IMIGTABLE(I).NE.3)THEN
                       JEDNMIG=JEDNMIG+1
                       IDTMIG(KID,N)=JEDNMIG
                    ELSE  
                       IF(IVMIGV(KID,N).GT.0) THEN
                          JEDNMIG=JEDNMIG+1
                          IDTMIG(KID,N)=JEDNMIG
                       ENDIF
                   ENDIF 
                ENDDO
                JEDNMIGV(I)=JEDNMIG
                WRITE(3,2034)I, JEDNMIG
 2034           FORMAT(' CELL NUMBER:',I5,5X,'NUMBER OF EQUATIONS:',I7/
     1         ' NODE  EQUATION, NODE   EQUATION.....')
                WRITE(3,2035)(N,IDTMIG(KID,N),N=1,NPT)
 2035           FORMAT(8(I3,I7,7X)) 
                IF(IMIGTABLE(I).EQ.3)THEN
                   WRITE(3,2036)
2036               FORMAT(' NODAL VELOCITIES'/
     1             ' NODE  VELOCITY,  NODE  VELOCITY....')
                   WRITE(3,2037)(N,VMIGV(KID,N),N=1,NPT)              
                ENDIF
            ENDIF
         ENDDO
      ENDIF  
C
      IF(MULTIMOL.GT.1) GO TO 500
C
C     ONE MOLECULE CASE
C
   50 JEDN  = 0
      NZADP = 0
      NZADLYMPH=0
      NZADVS=0
      JEDNPRES=0
      JEDNCONC=0     
      
      IF(NDIMIDG.GT.2) THEN
C
C       CASE WHEN IDIMIDG.GE.2-----------
C       
         DO 270 N=1,NPT
         DO 260 II=1,NDIMIDG
C
         IF(IPIPENODV(N).EQ.1) THEN
              IF(II.NE.1.AND.II.NE.3) THEN
                 IDT(II,N)=0
                 GO TO 260
              ENDIF
              IF(ICONV.EQ.1) THEN
                IF(INDIFFUS.EQ.0.AND.II.EQ.3) THEN
                   IDT(II,N)=0
                   GO TO 260  
                ENDIF
              ELSE
                IF(II.NE.3) THEN
                    IDT(II,N)=0
                    GO TO 260
                ENDIF
              ENDIF
         ELSE
            IF(ISMEARED.EQ.0) THEN
               IF(ICONV.EQ.1) THEN
                  IF(II.NE.2.AND.II.NE.4) THEN
                     IDT(II,N)=0
                     GO TO 260         
                  ENDIF
                  IF(INDIFFUS.EQ.0.AND.II.EQ.4) THEN
                     IDT(II,N)=0
                     GO TO 260
                  ENDIF
                  ! CELLS - NO PRESSURE:
                  IF(II.LE.2.AND.IDT(II,N).EQ.2) THEN
                     IDT(II,N)=0
                     GO TO 260
                  ENDIF
               ELSE
                  IF(II.NE.4) THEN
                      IDT(II,N)=0
                      GO TO 260
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF(IDT(II,N)) 252,253,254
  254    IDT(II,N)=0
         IF(IACTIVE(II).GT.0) THEN
              NZADP=NZADP+1
         ENDIF
         GO TO 260
  253    IF(II.LE.2) THEN 
            JEDNPRES=JEDNPRES+1
            IDT(II,N) = JEDNPRES 
         ELSE
            JEDNCONC=JEDNCONC+1
            IDT(II,N) = JEDNCONC            
         ENDIF  
         GO TO 260
  252    IF(II.LE.2) THEN
            JEDNPRES=JEDNPRES+1
            IDT(II,N) = JEDNPRES
            NZADP=NZADP+1 
         ELSE
            JEDNCONC=JEDNCONC+1
            IDT(II,N) = JEDNCONC            
         ENDIF 
C
  260    CONTINUE
C
C        ADD EQUATIONS FOR LYMPH
C
         IF(INDLYMPHH.GT.0) THEN
            IF(IPIPENODV(N).EQ.1) THEN
              IDLYMPHH(N)=0
            ELSE
              IF(IDLYMPHH(N).EQ.0) THEN
                JEDNPRES=JEDNPRES+1
                IDLYMPHH(N)=JEDNPRES
              ELSE
                NZADLYMPH=NZADLYMPH+1
                IDLYMPHH(N)=0
              ENDIF
            ENDIF
         ENDIF
         IF(INDLYMPHC.GT.0) THEN
            IF(IPIPENODV(N).EQ.1) THEN
              IDLYMPHC(1,N)=0
            ELSE
              IF(IDLYMPHC(1,N).EQ.0) THEN
                JEDNCONC=JEDNCONC+1
                IDLYMPHC(1,N)=JEDNCONC
              ELSE
                NZADLYMPH=NZADLYMPH+1
                IDLYMPHC(1,N)=0
              ENDIF
            ENDIF
         ENDIF
         IF(MULTITIS.EQ.0.OR.INDIFFUS.EQ.0.OR.MULTIMOL.GT.1) GO TO 270
C
C        COMPOSITE TISSUE AND DIFFUSION       
C
         KIN=0
         DO I=1,NDIMCELL
            IF(IPIPENODV(N).EQ.1) THEN 
               IDTVCELL(1,I,N)=0
            ELSE 
               IF(IDTVCELL(1,I,N).EQ.0) THEN
                  JEDNCONC=JEDNCONC+1
                  IDTVCELL(1,I,N)=JEDNCONC
               ELSE
                  NZADP=NZADP+1 
                  IDTVCELL(1,I,N)=0
               ENDIF
            ENDIF 
            IF(INDCELLIN.GT.0) THEN
               NINID=INDCELLINV(I)
               IF(NINID.GT.0) THEN
                  KIN=KIN+1
                  DO K=1,NINID
                     IF(IPIPENODV(N).EQ.1) THEN
                       IDVS(MM,KIN,K,N)=0
                     ELSE
                       IF(IDVS(MM,KIN,K,N).EQ.0) THEN
                         JEDNCONC=JEDNCONC+1
                         IDVS(MM,KIN,K,N)=JEDNCONC 
                       ELSE
                         IDVS(MM,KIN,K,N)=0
                         NZADVS=NZADVS+1
                       ENDIF
                     ENDIF
                  ENDDO 
               ENDIF
            ENDIF 
         ENDDO
C           
  270    CONTINUE
         GO TO 71
      ENDIF
C      
C     CASE WHEN IDIMIDG.LE.1 AND MULTIMOL.EQ.1
C
      DO 70 N=1,NPT
      DO 60 II=1,NDIMIDG
      IF(IDT(II,N)) 52,53,54
   54 IDT(II,N)=0
      IF(indiffus.EQ.0.AND.IHEAT.EQ.0) THEN
        IF(II.EQ.1) NZADP=NZADP+1
      ELSE
        IF(ICONV.NE.0) THEN
          IF(IDARCY.GE.0) NZADP=NZADP+1
        ELSEIF(II.EQ.2) THEN
          IF(IDARCY.GE.0) NZADP=NZADP+1
        ENDIF
      ENDIF
      GO TO 60
   53 JEDN=JEDN+1
      IDT(II,N) = JEDN
      GO TO 60
   52 JEDN      = JEDN + 1
      IDT(II,N) = JEDN
    
      IF(IDARCY.NE.0.OR.ISMEARED.NE.0) then
        IF(IDARCYMA(N).EQ.-1.AND.ICONV.NE.0) THEN
          NZADP=NZADP+1
        ENDIF
      ENDIF
   60 CONTINUE  
C        ADD EQUATIONS FOR LYMPH
      IF(INDLYMPH.EQ.0) GO TO 70  
      IF(INDLYMPHH.GT.0) THEN
          IF(IDLYMPHH(N).EQ.0) THEN
              JEDN=JEDN+1
              IDLYMPHH(N)=JEDN
          ELSE
              NZADLYMPH=NZADLYMPH+1
              IDLYMPHH(N)=0
          ENDIF
       ENDIF
       IF(INDLYMPHC.GT.0) THEN
          IF(IDLYMPHC(1,N).EQ.0) THEN
              JEDN=JEDN+1
              IDLYMPHC(1,N)=JEDN
          ELSE
              NZADLYMPH=NZADLYMPH+1
              IDLYMPHC(1,N)=0
          ENDIF
       ENDIF   
   70 CONTINUE
   
   71 CONTINUE
      IF(MULTIMOL.EQ.1) GO TO 171
      
  500 CONTINUE  
C
C SET MATRICES FOR MULTIMOL.GT.1 -----------------
C       CASE WHEN IDIMIDG.GE.2-----------
C  
      DO 550 MM=1,MULTIMOL
         JEDNPV(MM)=0 
         JEDNCV(MM)=0
         NZADPV(MM)=0  
         NZADPLYV(MM)=0                             
         DO 570 N=1,NPT
         JEDNP=0
         JEDNC=0
         DO 560 II=1,NDIMIDG   
         IF(IPIPENODV(N).EQ.1) THEN
              IF(II.NE.1.AND.II.NE.3) THEN
                 IDTV(MM,II,N)=0
                 GO TO 560
              ENDIF
              IF(ICONV.EQ.1) THEN
                IF(INDIFFUS.EQ.0.AND.II.EQ.3) THEN
                   IDTV(MM,II,N)=0
                   GO TO 560  
                ENDIF
              ELSE
                IF(II.NE.3) THEN
                    IDTV(MM,II,N)=0
                    GO TO 560
                ENDIF
              ENDIF
         ELSE
            IF(ISMEARED.EQ.0) THEN
               IF(ICONV.EQ.1) THEN
                  IF(II.NE.2.AND.II.NE.4) THEN
                      IDTV(MM,II,N)=0
                      GO TO 560                
                  ENDIF
                  IF(INDIFFUS.EQ.0.AND.II.EQ.4) THEN
                      IDTV(MM,II,N)=0
                      GO TO 560  
                  ENDIF
               ELSE
                  IF(II.NE.4) THEN
                      IDTV(MM,II,N)=0
                      GO TO 560                
                  ENDIF             
               ENDIF
            ENDIF
         ENDIF
         IF(IDTV(MM,II,N)) 552,553,554
  554    IDTV(MM,II,N)=0
         IF(IACTIVE(II).GT.0) THEN
              NZADPV(MM)=NZADPV(MM)+1 
         ENDIF
         GO TO 560
  553    IF(II.LE.2) THEN 
            JEDNPV(MM)=JEDNPV(MM)+1
            IDTV(MM,II,N) = JEDNPV(MM)
            JEDNP=JEDNPV(MM)
         ELSE
            JEDNCV(MM)=JEDNCV(MM)+1
            IDTV(MM,II,N) = JEDNCV(MM)
            JEDNC=JEDNCV(MM)            
         ENDIF  
         GO TO 560
  552    IF(II.LE.2) THEN
            JEDNPV(MM)=JEDNPV(MM)+1              
            IDTV(MM,II,N) = JEDNPV(MM)
            JEDNP=JEDNPV(MM)
            NZADPV(MM)=NZADPV(MM)+1 
         ELSE
            JEDNCV(MM)=JEDNCV(MM)+1                            
            IDTV(MM,II,N) = JEDNCV(MM)
            JEDNC=JEDNCV(MM)            
         ENDIF 
C
  560    CONTINUE
C---     ADD EQUATIONS FOR LYMPH  
         IF(INDLYMPHH.GT.0) THEN
            IF(IDLYMPHH(N).EQ.0) THEN
              JEDNP=JEDNP+1
              IDLYMPHH(N)=JEDNP
              JEDNPV(MM)=JEDNP
            ELSE
              NZADPLYV(MM)=NZADPLYV(MM)+1
              IDLYMPHH(N)=0
            ENDIF
         ENDIF
         IF(INDLYMPHC.GT.0) THEN
            IF(IDLYMPHC(MM,N).EQ.0) THEN
              JEDNC=JEDNC+1
              IDLYMPHC(MM,N)=JEDNC
              JEDNCV(MM)=JEDNC
            ELSE
              NZADPLYV(MM)=NZADPLYV(MM)+1
              IDLYMPHC(MM,N)=0
            ENDIF
         ENDIF 
C---     END OF LYMPH          
         IF(MULTITIS.EQ.0.OR.INDIFFUS.EQ.0) GO TO 570
C
C        COMPOSITE TISSUE AND DIFFUSION   
         KIN=0
         DO I=1,NDIMCELL
            IF(IPIPENODV(N).EQ.1) THEN
               IDTVCELL(1,I,N)=0
            ELSE          
               IF(IDTVCELL(MM,I,N).EQ.0) THEN   
                  JEDNCV(MM)=JEDNCV(MM)+1
                  IDTVCELL(MM,I,N)=JEDNCV(MM)
               ELSE
                  NZADPV(MM)=NZADPV(MM)+1 
                  IDTVCELL(MM,I,N)=0
               ENDIF
            ENDIF
            IF(INDCELLIN.GT.0) THEN
               NINID=INDCELLINV(I)
               IF(NINID.GT.0) THEN
                  KIN=KIN+1
                  DO K=1,NINID
                     IF(IDVS(MM,KIN,K,N).EQ.0) THEN
                        JEDNCV(MM)=JEDNCV(MM)+1
                        IDVS(MM,KIN,K,N)=JEDNCV(MM)
                     ELSE
                        IDVS(MM,KIN,K,N)=0
                        NZADVS=NZADVS+1                   
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
  570 CONTINUE
  550 CONTINUE
  
      WRITE(3,6060)
 6060 FORMAT(/6X,'G E N E R A T E D   E Q U A T I O N  N U M B E R S'/
     1        6X,'             MULTIMOLECULAR MODEL'/
     2 6X,'MOLECULE NUMBER'/
     3 6X,'NODE',10X,'CONVECTION',4X,'DIFFUSION:  PIPE (CAP)',
     4 ' EXTRACELL SPACE, CELLS  1,2,....NDIMCELL')
      DO MM=1,MULTIMOL
         WRITE(3,6061)MM
6061     FORMAT(13X,I2)
         DO N=1,NPT
           IF(MULTITIS.EQ.0) THEN
              WRITE(3,6062)N,(IDTV(MM,II,N),II=1,NDIMIDG)
6062          FORMAT(I10,4I5)   
           ELSE
              DO I=1,NDIMCELL 
                  IDTLOC(I)=IDTVCELL(MM,I,N)
              ENDDO
             WRITE(3,6063)N,(IDTV(MM,I,N),I=1,NDIMIDG),(IDTLOC(I),
     1                       I=1,NDIMCELL)
6063         FORMAT(I10,5X,4I7,3X,10I7)
           ENDIF
         ENDDO
      ENDDO
      IF(MULTITIS.GT.0.AND.IMIGTOT.GT.0) THEN
         WRITE(3,6070)
 6070    FORMAT(/6X,'EQUATIONS FOR CELL M I G R A T I O N'/
     1   6X,'NODE',12X,' C E L L S')
         DO N=1,NPT
            DO I=1,NDIMCELL
               IDTLOC(I)=0
            ENDDO
            KK=0
            DO I=1,NDIMCELL
              IF(IMIGTABLE(I).NE.0) THEN
                 IF(IVMIGV(KK,N).GT.0) THEN
                    KK=KK+1
                    IDTLOC(I)=IVMIGV(KK,N)
                 ENDIF
              ENDIF 
            ENDDO  
            WRITE(3,6071) N,(IDTLOC(I),I=1,NDIMCELL)
 6071       FORMAT(I10,5X,4I5)
         ENDDO
      ENDIF      
      RETURN
C END OF MULTIMOL.GT.1-------------------------------------------
C

C
CE    PRINT READ-IN AND GENERATED NODAL DATA ----------------
C
  171 IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
C      CALL WBROJK(KARTI,1)
C
      IF(NDIMIDG.LE.2) THEN
         WRITE(3,6027)
      ELSE
         WRITE(3,6028)
      ENDIF
C
      DO 80 I=1,NPT
          DO J=1,4
             IDTLOC(J)=0
          ENDDO
C          KK=0
          IF(ISMEARED.EQ.1) THEN
             IDTLOC(1)=IDT(1,I)
             IDTLOC(2)=IDT(2,I)
          ELSE
             DO J=1,4
                IDTLOC(J)=IDT(J,I)
             ENDDO
          ENDIF
          WRITE(3,5020) I,(IDTLOC(J),J=1,4),(CORD(I,J),J=1,3)
   80 CONTINUE
   
      IF(MULTITIS.EQ.0.OR.INDIFFUS.EQ.0) THEN
         IF(INDLYMPH.GT.0) GO TO 577      
         RETURN
      ENDIF
C
C     COMPOSITE TISSUE AND DIFFUSION
      WRITE(3,6400)
 6400 FORMAT(/6X,' EQUATIONS FOR DIFFUSION WITHIN COMPOSITE TISSUE'/
     1 6X,'MULECULE NUMBER'/
     1 16X,'NODE',20X,'C E L L S ....  1-10MAX') 
      DO MM=1,MULTIMOL
         WRITE(3,5033)MM
5033     FORMAT(/6X,'M O L E C U L E NUMBER.....',I5)
         DO N=1,NPT
            WRITE(3,5030) N,(IDTVCELL(MM,I,N),I=1,NDIMCELL)   
         ENDDO 
         IF(INDCELLIN.GT.0) THEN 
             WRITE(3,5037)
 5037       FORMAT(/' DIFFUSION WITHIN CELLS - VESICLES EQUATIONS')  
            KIN=0
            DO NC=1,NDIMCELL
               NINID=INDCELLINV(NC)          
               IF(NINID.GT.0) THEN
                  WRITE(3,5038)NC
 5038             FORMAT('    CELL NUMBER=',I5)
                  KIN=KIN+1
                  DO N=1,NPT
                     WRITE(3,5030) N,(IDVS(MM,KIN,K,N),K=1,NINID)  
                  ENDDO
               ENDIF
            ENDDO 
        ENDIF   
      ENDDO 
      IF(INDLYMPH.EQ.0) RETURN
C      
  577 CONTINUE
C     PRINT LYMPHATIC SYSTEM EQUATIONS
      WRITE(3,5130)
 5130 FORMAT(/' EQUATIONS FOR   L Y M P H A T I C  S Y S T E M'/
     1        ' NODE - EQUATION NUMBERS')
      IF(INDLYMPHH.GT.0) THEN
         WRITE(3,5131)
 5131    FORMAT(' CONVECTION:  NODE EQ NUMBER  NODE  EQ NUMBER......')
         KG=0
 5220    CONTINUE
         DO 5250 I = 1,8
            KG = KG +1
            IF(KG.GT.NPT) GO TO 5300     
            IDTLOC(I) = KG   
            IDTLOCI(I)=IDLYMPHH(KG)
 5250    CONTINUE
 5300    I = I - 1
         IF(I.NE.0) THEN
            WRITE(3,5132) (IDTLOC(J),IDTLOCI(J),J=1,I)
         ENDIF
         IF(KG.LT.NPT) GO TO 5220
      ENDIF
 5132 FORMAT(8(2I6,2X))    
      IF(INDLYMPHC.EQ.0) RETURN
C     DIFFUSION
      WRITE(3,5431)
 5431 FORMAT(' DIFFUSION:  NODE EQ NUMBER  NODE  EQ NUMBER......')  
      DO 5410 MM=1,MULTIMOL
         WRITE(3,5432)MM
 5432    FORMAT(/' MOLECULE NUMBER =',I5)
         KG=0
 5420    CONTINUE
         DO 5450 I = 1,8
            KG = KG +1
            IF(KG.GT.NPT) GO TO 5400     
            IDTLOC(I) = KG   
            IDTLOCI(I)=IDLYMPHC(MM,KG)
 5450    CONTINUE
 5400    I = I - 1
         IF(I.NE.0) THEN
            WRITE(3,5132) (IDTLOC(J),IDTLOCI(J),J=1,I)
         ENDIF
         IF(KG.LT.NPT) GO TO 5420
 5410 CONTINUE   
      RETURN      
 5030 FORMAT(11X,I10,5X,10I7)      
C
 2027 FORMAT(6X,'G E N E R I S A N I   P O D A C I   O   C V O R O V I M
     1 A'/6X,57('-')///11X,'CVOR',3X,'JEDN. BR.(>0)',
     19X,'K O O R D I N A T E'/11X,
     2'BROJ',3X,'KONS.TEMP.(=0)',4X,'X',10X,'Y',11X,'Z')
 6027 FORMAT(/6X,'G E N E R A T E D   N O D A L   D A T A'/
     16X,39('-')//
     116X,' NODE',29X,'EQUATION NUMBERS',28X,'C O O R D I N A T E S'
     2/86X,'X',16X,'Y',16X,'Z')
 6028 FORMAT(/6X,'G E N E R A T E D   N O D A L   D A T A'/
     16X,39('-')//
     116X,' NODE',29X,'EQUATION NUMBERS',28X,'C O O R D I N A T E S'
     2/43X,'PRESSURES        CONCENTRATIONS',12X,'X',16X,'Y',16X,'Z'/
     3 38X  'CAPILLAR  TISSUE CAPILLAR TISSUE(EX CELL SPACE)' )
C 5020 FORMAT(10X,2I10,6X,3E15.6)
 5020 FORMAT(11X,I10,12X,4I10,5X,E15.6,3X,E15.6,3X,E15.6)
C
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ULAZT2(TABK,TABC,NTAKV,NTACV,IBFK,
     1FAKP,GUSM,TOPM,TMNM,INTABK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*130 ACOZ
      COMMON /TEMPGL/ NPT,NGET,NMATT,NSTAC,MAXIT,NDT,NTCONS,
     1NET,NGE,JEDN,DT,TMAX,VREME,NPOC,NWK,NPRINT,LINTE,
     2KONF,INDSK,INDF,KOR,NTOTA,ITER,INDSC,INTEB,INDAX,IANIZ,NQE,IZIP,
     3ITOR,INDVEL,INDFL
      COMMON /REPEAT/ LCORD,LTABK,LTABC,LNTAKV,LNTACV,LIDTE,LIDT,LTCONS,
     1LMAX,LMAXA,LMAXM,LMHT,LNTAQE,LQEFN,LNTAQP,LQPFN,LNTAOK,LTOKFN,
     2LNTAHP,LHPFN,LNTAHR,LHRFN,LTT1,LTT2,LTT0,LSK,LCONM,LF,LFC,LNEL,
     3LIBFK,LFAKP,LISKC,LGUSM,LTOPM,LTMNM,LFOHR,LFVPOT
      COMMON /FUNKTE/ NANLK,NTABK,NANLC,NTABC,MAXTK,MAXTC,NQP,MAXTQP,
     2NHR,MAXTHR,NHP,MAXTHP,NTOK,MAXTOK,MAXTQE,KONK,KONC,NFCON,SBC,NFO
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
	common/diffusmain/indiffus,NDS,maxdt,NPRINTnds,lvreme,leprdt,
	1                  NUMMAT,maxPmat,LMATER,LMATFUN,LAREAD,LFLUXDX,
     1                  LFLUXDY,LFLUXDZ                    
      DIMENSION TABK(2,INTABK,*),TABC(2,NMATT,*),
     1NTAKV(*),NTACV(*),IBFK(NMATT,*),FAKP(NMATT,*),
     2GUSM(*),TOPM(*),TMNM(*)
CE     MATERIAL CONSTANTS TEMPERATURE DEPENDENCE
      DO 34 I=1,NMATT
      NTACV(I)=0
      DO 35 J=1,3
      IBFK(I,J)=0
   35 FAKP(I,J)=0.
   34 CONTINUE
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2036)
      IF(ISRPS.EQ.1)
     *WRITE(3,6036)
      ENDIF
 2036 FORMAT(/////' ','BROJEVI FUNKCIJA ZAVISNOSTI KOEFI. PROVODJENJA OD
     1 TEMPERATURE'/' ',61('-')///
     1' (ZA SVAKI OD "NMATT" MATERIJALA UCITAVA SE PO JEDNA KARTICA)'//
     1     /' MATER. BR.',3X,'BROJEVI FUNKCIJE',5X,'KOEFICIJENTI PROPORC
     2IONALNOSTI',5X,'GUST. MAT.',2X,'MAX.DOZ.TEMP.',2X,'MIN.DOZ.TEMP.')
c 6036 FORMAT(/////' ','DATA ABOUT CONDUCTION COEFFICIENT DEPENDENCE ON T
c     1EMPERATURES'/' ',60('-')///
c     1     /' MATER. No.',3X,'FUNCTION NUMBER ',7X,'MULTIPLIERS COEFICIE
c     2NT    ',5X,'DENS. MAT.',2X,'MAX.ALO.TEMP.',2X,'MIN.ALO.TEMP.')
 6036 FORMAT(/////' ','DATA ABOUT DIFFUSION COEFFICIENT DEPENDENCE ON C
     1ONCENTRATION'/' ',60('-')///
     1     /' MAT NO',4X,'FUN NOS',12X,'MULTIPLIER COEFICIENT
     2S',10X,'MAT DENS',5X,'MAX CONC ',4X,'MIN CONC')
      DO 37 I=1,NMATT
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) (IBFK(I,J),J=1,3),(FAKP(I,K),K=1,3),GUSM(I),TOPM(I),
     1TMNM(I)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (IBFK(I,J),J=1,3),(FAKP(I,K),K=1,3),GUSM(I),
     1TOPM(I),TMNM(I)
 1010 FORMAT(3I5,6F10.0)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(3,5035) I,(IBFK(I,J),J=1,3),(FAKP(I,K),K=1,3),GUSM(I),TOPM(I
     1),TMNM(I)
 5035 FORMAT(2X,I5,5X,3I5,3X,4(1PE12.4),2(1PE14.4))
   37 CONTINUE
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2002)
      IF(ISRPS.EQ.1)
     *WRITE(3,6002)
      ENDIF
 2002 FORMAT(//////' ZAVISNOST KOEFICIJENTA PROVODJENJA OD TEMPERATURE'/
     1' ',49('-')///' ZA SVAKU OD "NTABK" ZADATIH FUNKCIJA:'/
     1'  (FUNKCIJA BROJ I BROJ "NMAX" TACAKA ZA DEF. FUNKCIJE)'/
     2'  ("NMAX" PAROVA VREDN.: TEMPER. - KOEF. PROVODJENJA)')
C 6002 FORMAT(//////' CONDUCTION COEFFICIENT DEPENDENCE ON TEMPERATURE'/
C     1' ',48('-')///)
 6002 FORMAT(//////' DIFFUSION COEFFICIENT DEPENDENCE ON CONCENTRATION'/
     1' ',49('-')///' FOR EACH OF "NTABK" SPECIFIED FUNCTIONS:'/
     2'  (MAXIMUM "NMAX" PAIRS OF: CONCENTRATION - DIFFUSION COEFF)')
      IF (NTABK.EQ.0) GO TO 150
      DO 36 I=1,NTABK
   36 NTAKV(I)=0
CE    DIFFUSION COEFFICIENT CONCENTRATION DEPENDENCE - GIVEN AS TABLE
      DO 140 IK=1,NTABK
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) MM,NTAKV(MM)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1015) MM,NTAKV(MM)
 1015 FORMAT(2I5)
      IF(NTAKV(MM).EQ.0) NTAKV(MM)=1
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
         IF(ISRPS.EQ.0)
     *   WRITE(3,2003) MM,NTAKV(MM)
         IF(ISRPS.EQ.1.AND.IK.EQ.1)
     *   WRITE(3,6003) MM,NTAKV(MM)
      IF(ISRPS.EQ.1.AND.IK.GT.1)
     *WRITE(3,6004) MM,NTAKV(MM)
      ENDIF
 2003 FORMAT(//' FUNKCIJA BR.',I5,3X,'BROJ TACAKA NA KRIVOJ=',I5,3X,'(=1
     1,-CONST.KOEF.PROV.)'//'  TEMPERATURA',8X,'KOEFICIJENT PROVODJENJA'
     1)
 6003 FORMAT(//' FUNCTION NO=',I5,3X,'NUMBER OF PTS ON CURVE=',I5,5X,
     1 '(=1 - CONST DIFF COEF)'//'  CONCENTRATION',8X,
     1'DIFFUSION COEFICIENT')
 6004 FORMAT(//' FUNCTION NO=',I5,3X,'NUMBER OF PTS ON CURVE=',I5,5X,
     1//'  CONCENTRATION',8X,'DIFFUSION COEFICIENT')
      NMX=NTAKV(MM)
      IF(NMX.LE.MAXTK) GO TO 137
      IF(ISRPS.EQ.0)
     *WRITE(3,2005) MM,NMX,MAXTK
      IF(ISRPS.EQ.1)
     *WRITE(3,6005) MM,NMX,MAXTK
 2005 FORMAT(1X,'GRESKA U ULAZNIM PODACIMA :'/ ' MATERIJAL BROJ=',I3/
     1' BROJ TACAKA NA KRIVOJ ZAVISNOSTI KOEF.PROVODJ.OD TEMP.VECI OD MA
     2X.ZADATOG'/  ' BROJ TACAKA NA KRIVOJ=',I3/  ' ZADATI MAX.BROJ TACA
     3KA , MAXTK=',I3)
 6005 FORMAT(1X,'ERROR IN INPUT DATA FOR CONDUCTION COEFFICIENT:'/
     1' MATERIAL NUMBER =',I3/
     2' NUMBER POINTS ON CURVE =',I3/
     3' MAXIMAL NUMBER POINTS ON CURVE , MAXTK =',I3)
      STOP
  137 DO 138 J=1,NMX
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) (TABK(I,MM,J),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1016) (TABK(I,MM,J),I=1,2)
 1016 FORMAT(2E20.6)
 138  CONTINUE
 
 1001 FORMAT(2F10.0)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(3,5004) ((TABK(I,MM,J),I=1,2),J=1,NMX)
 5004 FORMAT(E20.6,5X,E20.6)
  140 CONTINUE
      GO TO 170
  150 IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2008)
      IF(ISRPS.EQ.1)
     *WRITE(3,6008)
      ENDIF
 2008 FORMAT(///' SVE ZAVISNOSTI KOEFICIJENATA PROVODJENJA TOPLOTE OD TE
     1MPERATURE - ANALITICKE')
 6008 FORMAT(///' ALL CONDUCTION COEFICIENT TEMPERATURE DEPENDENCE ARE',
     1' GIVEN IN ANALITICAL FORM ')
c  170 IF(NSTAC.EQ.0) RETURN
  170 if(NANLC.EQ.0 .AND. NTABC.EQ.0) RETURN
C
C     FOR UNSTEADY DIFFUSION, ONE TABLE FOR SPECIFIC HEAT, WITH UNIT C 
      IF(indiffus.EQ.1) THEN
	   DO I=1,NMATT
           NTACV(I)=1
	     TABC(1,I,1)=0.0
	     TABC(2,I,1)=1.D0
	   ENDDO
	   RETURN
	ENDIF
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2025)
      IF(ISRPS.EQ.1)
     *WRITE(3,6025)
      ENDIF
 6025 FORMAT(//////' SPECIFIC HEAT TEMPERATURE DEPEDENCE '/' ',35
     1('-')///)
 2025 FORMAT(//////' ZAVISNOST SPECIFICNE TOPLOTE OD TEMPERATURE'/' ',43
     1('-')///' ZA SVAKU OD "NTABC" ZADATIH FUNKCIJA:'/
     1'  (MATERIJAL BROJ I BROJ "NMX" TACAKA ZA DEF. FUNK.)'/
     2'  ("NMAX" PAROVA VREDN.: TEMPER. - SPECIF. TOPLOTA)')
      IF (NTABC.GT.0) GO TO 180
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2026)
      IF(ISRPS.EQ.1)
     *WRITE(3,6026)
      ENDIF
 6026 FORMAT(///' ALL SPECIFIC HEAT TEMPERATURE DEPEDENCE ARE GIVEN',
     1' ANALITICAL')
 2026 FORMAT(///' SVE ZAVISNOSTI SPECIFICNE TOPLOTE OD TEMPERATURE - ANA
     1LITICKE')
      RETURN
CE    SPECIFIC HEAT TEMPERATURE DEPEDENCE GIVEN AS TABLE
  180 DO 190 JC=1,NTABC
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) NN,NTACV(NN)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) NN,NTACV(NN)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2016) NN,NTACV(NN)
      IF(ISRPS.EQ.1)
     *WRITE(3,6016) NN,NTACV(NN)
      ENDIF
 2016 FORMAT(//' MATERIJAL BR.',I5,3X,'BROJ TACAKA NA KRIVOJ=',I3,3X,'(=
     11,-CONST.SPEC.TOPL.)'//'  TEMPERATURA',10X,'SPECIFICNA TOPLOTA')
 6016 FORMAT(//' MATERIAL No.=',I5,3X,'NUM. OF POINTS ON CURVE =',I3,3X,
     1'(=1,-CONST.SPEC.HEAT.)'//'  TEMPERATURE',12X,'SPECIFIC HEAT')
      NMX=NTACV(NN)
      IF (NMX.LE.MAXTC) GO TO 185
      IF(ISRPS.EQ.0)
     *WRITE(3,2017) NN,NMX,MAXTC 
      IF(ISRPS.EQ.1)
     *WRITE(3,6017) NN,NMX,MAXTC
 2017 FORMAT(/1X,'GRESKA U ULAZNIM PODACIMA :'/1X,'MATERIJAL BROJ=',I3
     1/1X,'BROJ TACAKA NA KRIVOJ ZAVISNOSTI SPEC. TOPLOTE OD TEMP. VECI
     2OD MAX.ZADATOG'/1X,'BROJ TACAKA NA KRIVOJ=',I3/1X,'ZADATI MAX.BROJ
     3 TACAKA , MAXTC=',I3)
 6017 FORMAT(1X,'ERROR IN INPUT DATA FOR SPECIFIC HEAT:'/
     1' MATERIAL NUMBER =',I3/
     2' NUMBER POINTS ON CURVE =',I3/
     3' MAXIMAL NUMBER POINTS ON CURVE , MAXTC =',I3)
      STOP
  185 DO 186 J=1,NMX
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) (TABC(I,NN,J),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1001) (TABC(I,NN,J),I=1,2)
  186 CONTINUE
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(3,5004) ((TABC(I,NN,J),I=1,2),J=1,NMX)
      DO 191 J=1,NMX
  191 TABC(2,NN,J)=TABC(2,NN,J)*GUSM(NN)
  190 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE KONVTEPD(F,DTT,TT1,IDT,IDTCD,IACTIVE,JEDN,ITER,TREFP,
     1  DTMODP,TREFD,DTMODD,EPSTA,EPSTR,KONVP,KONVD,IPRINT,NPRINT,KOR,
     2  NDT,VREME,NPT,NDIMIDT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON/MESSAGES/MESSAGE
      COMMON /SRPSKI/ ISRPS
      COMMON/FITNODES/INFF
C 
      DIMENSION F(JEDN),TT1(JEDN),DTT(JEDN),IDT(NDIMIDT,NPT),
     1       IDTCD(2,NPT),IACTIVE(4)
C
C      IDOF   = KK
	TOL    = 1.D-10
	TOLABS = 1.D-8
	TOLABSP=TOLABS
	TOLABSD=TOLABS
C 
      IDOF=0
      IPRTOT=0
      IDDTOT=0
      TT1ABSP=0.D0
      DTMODP=0.D0
      TT1ABSD=0.D0
      DTMODD=0.D0
      IF(ITER.EQ.0) THEN
         TREFP=0.D0
         TREFD=0.D0
      ENDIF
      DO 100 ID=1,4
      IF(IACTIVE(ID).EQ.0) GO TO 100
      IDOF=IDOF+1
      IPR=0
      IF(ID.EQ.1.OR.ID.EQ.3) THEN
        IPR=1
        IPRTOT=1
      ENDIF
      IDD=0
      IF(ID.EQ.2.OR.ID.EQ.4) THEN
         IDD=1
         IDDTOT=1 
      ENDIF
C
      DO 30 I=1,NPT
           II=IDT(IDOF,I)
           IF(II.GT.0) THEN
             IF(IPR.EQ.1) THEN
	           TT1ABSP = TT1ABSP + TT1(II)*TT1(II)   
	           TMAXP=TT1ABSP
                 DTMODP  = DTMODP  + F(II)*F(II)
             ELSEIF (IDD.EQ.1) THEN
                 TT1ABSD = TT1ABSD + TT1(II)*TT1(II)
	           TMAXD=TT1ABSD
                 DTMODD  = DTMODD  + F(II)*F(II)
             ENDIF
          ENDIF
   30 CONTINUE
  100 CONTINUE 
C     END OF LOOP OVER DOMAINS
C
      DTMODP  = DSQRT(DTMODP)
      DTMODD  = DSQRT(DTMODD)
	TT1ABSP = DSQRT(TT1ABSP)
	TT1ABSD = DSQRT(TT1ABSD)
	TMAXP=TT1ABSP
	TPRINTP = TT1ABSP
	TMAXD=TT1ABSD
	TPRINTD = TT1ABSD
C	IF(ITER.EQ.0) UNBREF=DTMOD
      IF(TMAXP.LE.TREFP) THEN
          TMAXP = TREFP
      ELSE
          TREFP = TMAXP
      ENDIF
      TOLTP = TMAXP*EPSTR
      IF(TMAXD.LE.TREFD) THEN
          TMAXD = TREFD
      ELSE
          TREFD = TMAXD
      ENDIF
      TOLTD = TMAXD*EPSTR      
   
C
C     CRITERIONS FOR PRESSUR
      IF(IPRTOT.EQ.1) THEN 
         IF(TMAXP.LT.TOL) THEN    
	      TPRINTP=TT1ABSP
            ERRORABSP=TOL
	      IF(TT1ABSP.GT.TOL)ERRORABSP=TMAXP/TT1ABSP
	      ERRORRELP=TOL 
	  ELSE
	     ERRORABSP=DTMODP/TPRINTP
		   ERRORRELP=DTMODP/TMAXP
           IF(ERRORRELP.GT.EPSTR) KONVP = 0  
           if((EPSTA.gt.1.d-10).and.(dtmodP.le.EPSTA))konvP=1
	     IF(DTMODP.LT.TOL)KONVP=1 
c         dodat kriterijum za absolutnu toleranciju:TOLABS=0.05D0           
          INDABSP=0
	    if(DTMODP.lt.TOLABSP.AND.KONVP.EQ.0) THEN
	       INDABSP=1  
	       konvP=1
	       IF(INFF.EQ.0) RETURN 
             WRITE(MESSAGE,9040) TIME,KOR,ITER,dtmodP,TOLABS
 9040        FORMAT(/
     1  ' CONVERGENCE REACHED FOR INCREMENT OF PRESSURE SMALLER',
     2  ' THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1  ' TIME=',1PD12.5/
     1  ' STEP NO. =',I5/
	2  ' ITERATION NO. =',I5/
     3  ' INCERMENT =',1PD12.5/
     4  ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/) 
           ENDIF  
         ENDIF        
         IF(TPRINTP.LT.TOL) THEN
            ERRORABSP=0.D0
	   ELSE
	      ERRORABSP=DTMODP/TPRINTP
	   ENDIF 
      ENDIF	
C
C     CRITERIONS FOR DIFFUSION
      IF(IDDTOT.EQ.1) THEN 
         IF(TMAXD.LT.TOL) THEN    
	      TPRINTD=TT1ABSD
            ERRORABSD=TOL
	      IF(TT1ABSD.GT.TOL)ERRORABSD=TMAXD/TT1ABSD
	      ERRORRELD=TOL 
	  ELSE
	     ERRORABSD=DTMODD/TPRINTD
		   ERRORRELD=DTMODD/TMAXD
           IF(ERRORRELD.GT.EPSTR) KONVD = 0  
           if((EPSTA.gt.1.d-10).and.(dtmodD.le.EPSTA))konvD=1
	     IF(DTMODD.LT.TOL)KONVD=1 
c         dodat kriterijum za absolutnu toleranciju:TOLABS=0.05D0           
          INDABSD=0
	    if(DTMODD.lt.TOLABSD.AND.KONVD.EQ.0) THEN
	       INDABSD=1  
	       konvD=1
	       IF(INFF.EQ.0) RETURN 
             WRITE(MESSAGE,9041) TIME,KOR,ITER,dtmodD,TOLABS
 9041        FORMAT(/
     1  ' CONVERGENCE REACHED FOR INCREMENT OF CONCENTRATION SMALLER',
     2  ' THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1  ' TIME=',1PD12.5/
     1  ' STEP NO. =',I5/
	2  ' ITERATION NO. =',I5/
     3  ' INCERMENT =',1PD12.5/
     4  ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/) 
           ENDIF  
         ENDIF        
         IF(TPRINTD.LT.TOL) THEN
            ERRORABSD=0.D0
	   ELSE
	      ERRORABSD=DTMODD/TPRINTD
	   ENDIF 
      ENDIF	
	 
C------------------- LOOP FOR PRINTOUT OF PRESSURES AND CONCENTRATIONS ----------
C
      KK=0
      DO 400 IDOF=1,2
      IF(IDOF.EQ.1) THEN
         IF(IPRTOT.EQ.0) GO TO 400
         INDABS=INDABSP
      ELSEIF (IDOF.EQ.2) THEN
         IF(IDDTOT.EQ.0) GO TO 400
         INDABS=INDABSD
      ENDIF
      KK=KK+1 
      IF(ITER.EQ.0.AND.KK.EQ.1) THEN
          WRITE(*,2000)KOR,NDT
      ENDIF
      IF(IDOF.EQ.1) THEN
        DTMOD=DTMODP
        TPRINT=TPRINTP
        ERRORABS=ERRORABSP  
      ELSE
        DTMOD=DTMODD
        TPRINT=TPRINTD
        ERRORABS=ERRORABSD 
      ENDIF
      IF(KK.EQ.1) THEN
          WRITE(*,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS 
      ELSE
          WRITE(*,2021) IDOF,DTMOD,TPRINT,ERRORABS 
      ENDIF             
      IF(INDABS.EQ.1) THEN
        WRITE(*,2002) TOLABS
      ENDIF
	IF(NPRINT.GT.1) THEN
	  IF(KOR.GT.1.AND.KOR.LT.NDT) THEN
           IF(IPRINT.NE.(NPRINT-1)) GO TO 120
	  ENDIF    
      ENDIF
      IF(ISRPS.EQ.0) THEN
        WRITE(3,1000) ITER,DTMOD,ITER,ENE
        WRITE(*,1000) ITER,DTMOD,ITER,ENE
      ELSE
	  IF(ITER.EQ.0.AND.KK.EQ.1) THEN
	      IF(ILAMDA.EQ.0) THEN
	         WRITE(3,2000)KOR,NDT
	         IF(KOR.EQ.1) WRITE(3,2030)  
	      ELSE
	         WRITE(3,2020)KOR,NDT
	      ENDIF
        ENDIF
        IF(KK.EQ.1) THEN
            WRITE(3,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS
        ELSE
           WRITE(3,2021) IDOF,DTMOD,TPRINT,ERRORABS
        ENDIF
        IF(INDABS.EQ.1) THEN
          WRITE(3,2002) TOLABS
        ENDIF
	ENDIF
  120 CONTINUE
	
  400 CONTINUE
      RETURN
 1000 FORMAT(/'      ITERACIJA = ',I5,5X,'NORMA  = ',1PD12.4/
     1        '        NORMA (',I3,') / / NORMA (',I3,' )  = ',
     2               1PD12.4)
 2010 FORMAT(/' ITERATION =',I3,1X,': INCREMENT      = ',2X,
     1               1PD12.4,3X,'NORM     = ',1PD12.4/
     115X,'ICR/ICREMSTEP',5X,'ICR/NORM=',21PD12.4)
 2002 FORMAT(/' UNBALANCED NORM.LT.ABS-TOLERANCE = ',1PD12.4)
 2000 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM')
 2020 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM',
     2  5X,'UNBALANCED REF RHS',5X,'UBALANCED/REF RHS LAMDA') 
 2001 FORMAT(I5,8X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4)
 2021 FORMAT(13X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4)
 2030 FORMAT(' IDOF=1 - PRESSURES; IDOF=2 - CONCENTRATIONS')   
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE KONVTEPD2(F,DTT,TT1,IDT,IDTCD,IACTIVE,JEDN,ITER,TREF,
     1  DTMOD,EPSTA,EPSTR,KONV,IPRINT,NPRINT,KOR,NDT,VREME,NPT,NDIMIDT,
     2  IPASSME,ICONV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON/MESSAGES/MESSAGE
      COMMON /SRPSKI/ ISRPS
      COMMON/FITNODES/INFF
C 
      DIMENSION F(JEDN),TT1(JEDN),DTT(JEDN),IDT(NDIMIDT,NPT),
     1       IDTCD(2,NPT),IACTIVE(4)
C
C      IDOF   = KK
	TOL    = 1.D-10
	TOLABS = 1.D-8
C 
      IDOF=0
      IPRTOT=0
      IDDTOT=0
      TT1ABS=0.D0
      DTMOD=0.D0
      IF(ITER.EQ.0) TREF=0.D0
      IDIN=0
      IDOF=1
      IF(IPASSME.EQ.2) THEN
         IDIN=2
         IDOF=2
      ENDIF
C
      DO 30 I=1,JEDN
	     TT1ABS = TT1ABS + TT1(I)*TT1(I)   
	     TMAX=TT1ABS
           DTMOD  = DTMOD  + F(I)*F(I)
   30 CONTINUE
C
      DTMOD  = DSQRT(DTMOD)
	TT1ABS = DSQRT(TT1ABS)
	TMAX=TT1ABS
	TPRINT = TT1ABS
C	IF(ITER.EQ.0) UNBREF=DTMOD
      IF(TMAX.LE.TREF) THEN
          TMAX = TREF
      ELSE
          TREF = TMAX
      ENDIF
      TOLT = TMAX*EPSTR
C
      IF(TMAX.LT.TOL) THEN    
	      TPRINT=TT1ABS
            ERRORABS=TOL
	      IF(TT1ABS.GT.TOL)ERRORABS=TMAX/TT1ABS
	      ERRORREL=TOL 
      ELSE
	     ERRORABS=DTMOD/TPRINT
		   ERRORREL=DTMOD/TMAX
           IF(ERRORREL.GT.EPSTR) KONV = 0  
           if((EPSTA.gt.1.d-10).and.(dtmod.le.EPSTA))konv=1
	     IF(DTMOD.LT.TOL)KONV=1 
c          dodat kriterijum za absolutnu toleranciju:TOLABS=0.05D0           
           INDABS=0
	     if(DTMOD.lt.TOLABS) THEN
	         INDABS=1  
	         konv=1
C	         IF(INFF.EQ.0) RETURN 
	         IF(IPASSME.EQ.1) THEN
                  WRITE(MESSAGE,9040) VREME,KOR,ITER,dtmod,TOLABS
 9040             FORMAT(/
     1         ' CONVERGENCE REACHED FOR INCREMENT OF PRESSURE SMALLER',
     2         ' THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1         ' TIME=',1PD12.5/
     1         ' STEP NO. =',I5/
	2         ' ITERATION NO. =',I5/
     3         ' INCERMENT =',1PD12.5/
     4         ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/) 
            ELSE
                 WRITE(MESSAGE,9041) VREME,KOR,ITER,dtmod,TOLABS
 9041            FORMAT(/
     1         ' CONVERGENCE REACHED FOR INCREMENT OF CONCENTRATION',
     2         ' SMALLER THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1         ' TIME=',1PD12.5/
     1         ' STEP NO. =',I5/
	2         ' ITERATION NO. =',I5/
     3         ' INCERMENT =',1PD12.5/
     4         ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/) 
            ENDIF
          ENDIF               
      ENDIF        
      IF(TPRINT.LT.TOL) THEN
         ERRORABS=0.D0
	ELSE
	   ERRORABS=DTMOD/TPRINT
	ENDIF	 
C------------------- LOOP FOR PRINTOUT OF PRESSURES AND CONCENTRATIONS ----------
C
C      KK=0 
      IF(ITER.EQ.0) THEN
          IF(IPASSME.EQ.1.OR.(IPASSME.EQ.2.AND.
     1  ICONV.EQ.0)) THEN
             WRITE(*,2000)KOR,NDT
          ENDIF
      ENDIF
      IF(ITER.EQ.0) THEN
          WRITE(*,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS 
      ELSE
          WRITE(*,2021) IDOF,DTMOD,TPRINT,ERRORABS 
      ENDIF             
      IF(INDABS.EQ.1) THEN
        WRITE(*,2002) TOLABS
      ENDIF
	IF(NPRINT.GT.1) THEN
	  IF(KOR.GT.1.AND.KOR.LT.NDT) THEN
           IF(IPRINT.NE.(NPRINT-1)) GO TO 120
	  ENDIF    
      ENDIF
      IF(ISRPS.EQ.0) THEN
        WRITE(3,1000) ITER,DTMOD,ITER,ENE
        WRITE(*,1000) ITER,DTMOD,ITER,ENE
      ELSE
	  IF(ITER.EQ.0.AND.(IPASSME.EQ.1.OR.(IPASSME.EQ.2
     1  .AND.ICONV.EQ.0)))THEN
	      IF(ILAMDA.EQ.0) THEN
	         WRITE(3,2000)KOR,NDT
	         IF(KOR.EQ.1) WRITE(3,2030)  
	      ELSE
	         WRITE(3,2020)KOR,NDT
	      ENDIF
        ENDIF
C        IF(ITER.EQ.0) THEN
            WRITE(3,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS
C        ELSE
C           WRITE(3,2021) IDOF,DTMOD,TPRINT,ERRORABS
C        ENDIF
        IF(INDABS.EQ.1) THEN
          WRITE(3,2002) TOLABS
        ENDIF
	ENDIF
  120 CONTINUE
	
  400 CONTINUE
      RETURN
 1000 FORMAT(/'      ITERACIJA = ',I5,5X,'NORMA  = ',1PD12.4/
     1        '        NORMA (',I3,') / / NORMA (',I3,' )  = ',
     2               1PD12.4)
 2010 FORMAT(/' ITERATION =',I3,1X,': INCREMENT      = ',2X,
     1               1PD12.4,3X,'NORM     = ',1PD12.4/
     115X,'ICR/ICREMSTEP',5X,'ICR/NORM=',21PD12.4)
 2002 FORMAT(/' UNBALANCED NORM.LT.ABS-TOLERANCE = ',1PD12.4)
 2000 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM')
 2020 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM',
     2  5X,'UNBALANCED REF RHS',5X,'UBALANCED/REF RHS LAMDA') 
 2001 FORMAT(I5,8X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4)
 2021 FORMAT(13X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4)
 2030 FORMAT(' IDOF=1 - PRESSURES; IDOF=2 - CONCENTRATIONS')   
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE KONVTEPD2M(F,DTT,TT1,IDT,IDTCD,IACTIVE,JEDN,ITER,TREF,
     1  DTMOD,EPSTA,EPSTR,KONV,IPRINT,NPRINT,KOR,NDT,VREME,NPT,NDIMIDT,
     2  IPASSME,ICONV,MULTIMOL,MM) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON/MESSAGES/MESSAGE
      COMMON /SRPSKI/ ISRPS
      COMMON/FITNODES/INFF
C 
      DIMENSION F(JEDN),TT1(MULTIMOL,JEDN),DTT(JEDN),IDT(NDIMIDT,NPT),
     1       IDTCD(2,NPT),IACTIVE(4)
C
C      IDOF   = KK
	TOL    = 1.D-10
	TOLABS = 1.D-8
C 
      IDOF=0
      IPRTOT=0
      IDDTOT=0
      TT1ABS=0.D0
      DTMOD=0.D0
      IF(ITER.EQ.0) TREF=0.D0
      IDIN=0
      IDOF=1
      IF(IPASSME.EQ.2) THEN
         IDIN=2
         IDOF=2
      ENDIF
C
        DO 30 I=1,JEDN
	     TT1ABS = TT1ABS + TT1(MM,I)*TT1(MM,I)   
	     TMAX=TT1ABS
           DTMOD  = DTMOD  + F(I)*F(I)
   30   CONTINUE
C
      DTMOD  = DSQRT(DTMOD)
	TT1ABS = DSQRT(TT1ABS)
	TMAX=TT1ABS
	TPRINT = TT1ABS
C	IF(ITER.EQ.0) UNBREF=DTMOD
      IF(TMAX.LE.TREF) THEN
          TMAX = TREF
      ELSE
          TREF = TMAX
      ENDIF
      TOLT = TMAX*EPSTR
C
      IF(TMAX.LT.TOL) THEN    
	      TPRINT=TT1ABS
            ERRORABS=TOL
	      IF(TT1ABS.GT.TOL)ERRORABS=TMAX/TT1ABS
	      ERRORREL=TOL 
      ELSE
	     ERRORABS=DTMOD/TPRINT
		   ERRORREL=DTMOD/TMAX
           IF(ERRORREL.GT.EPSTR) KONV = 0  
           if((EPSTA.gt.1.d-10).and.(dtmod.le.EPSTA))konv=1
	     IF(DTMOD.LT.TOL)KONV=1 
c          dodat kriterijum za absolutnu toleranciju:TOLABS=0.05D0           
           INDABS=0
	     if(DTMOD.lt.TOLABS) THEN
	         INDABS=1  
	         konv=1
C	         IF(INFF.EQ.0) RETURN 
	         IF(IPASSME.EQ.1) THEN
                  WRITE(MESSAGE,9040) VREME,KOR,ITER,dtmod,TOLABS
 9040             FORMAT(/
     1         ' CONVERGENCE REACHED FOR INCREMENT OF PRESSURE SMALLER',
     2         ' THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1         ' TIME=',1PD12.5/
     1         ' STEP NO. =',I5/
	2         ' ITERATION NO. =',I5/
     3         ' INCERMENT =',1PD12.5/
     4         ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/) 
            ELSE
                 WRITE(MESSAGE,9041) VREME,KOR,ITER,dtmod,TOLABS
 9041            FORMAT(/
     1         ' CONVERGENCE REACHED FOR INCREMENT OF CONCENTRATION',
     2         ' SMALLER THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1         ' TIME=',1PD12.5/
     1         ' STEP NO. =',I5/
	2         ' ITERATION NO. =',I5/
     3         ' INCERMENT =',1PD12.5/
     4         ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/) 
            ENDIF
          ENDIF               
      ENDIF        
      IF(TPRINT.LT.TOL) THEN
         ERRORABS=0.D0
	ELSE
	   ERRORABS=DTMOD/TPRINT
	ENDIF	 
C------------------- LOOP FOR PRINTOUT OF PRESSURES AND CONCENTRATIONS ----------
C
C      KK=0 
      IF(ITER.EQ.0) THEN
          IF(IPASSME.EQ.1.OR.(IPASSME.EQ.2.AND.ICONV.EQ.0)) THEN
             WRITE(*,2000)KOR,NDT
          ENDIF
      ENDIF
      IF(ITER.EQ.0) THEN
          WRITE(*,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS 
      ELSE
          WRITE(*,2021) IDOF,DTMOD,TPRINT,ERRORABS 
      ENDIF             
      IF(INDABS.EQ.1) THEN
        WRITE(*,2002) TOLABS
      ENDIF
	IF(NPRINT.GT.1) THEN
	  IF(KOR.GT.1.AND.KOR.LT.NDT) THEN
           IF(IPRINT.NE.(NPRINT-1)) GO TO 120
	  ENDIF    
      ENDIF
      IF(ISRPS.EQ.0) THEN
        WRITE(3,1000) ITER,DTMOD,ITER,ENE
        WRITE(*,1000) ITER,DTMOD,ITER,ENE
      ELSE
	  IF(ITER.EQ.0.AND.(IPASSME.EQ.1.OR.(IPASSME.EQ.2
     1  .AND.ICONV.EQ.0)))THEN
	      IF(ILAMDA.EQ.0) THEN
	         WRITE(3,2000)KOR,NDT
	         IF(KOR.EQ.1) WRITE(3,2030)  
	      ELSE
	         WRITE(3,2020)KOR,NDT
	      ENDIF
        ENDIF
C        IF(ITER.EQ.0) THEN
            WRITE(3,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS
C        ELSE
C           WRITE(3,2021) IDOF,DTMOD,TPRINT,ERRORABS
C        ENDIF
        IF(INDABS.EQ.1) THEN
          WRITE(3,2002) TOLABS
        ENDIF
	ENDIF
  120 CONTINUE
	
  400 CONTINUE
      RETURN
 1000 FORMAT(/'      ITERACIJA = ',I5,5X,'NORMA  = ',1PD12.4/
     1        '        NORMA (',I3,') / / NORMA (',I3,' )  = ',
     2               1PD12.4)
 2010 FORMAT(/' ITERATION =',I3,1X,': INCREMENT      = ',2X,
     1               1PD12.4,3X,'NORM     = ',1PD12.4/
     115X,'ICR/ICREMSTEP',5X,'ICR/NORM=',21PD12.4)
 2002 FORMAT(/' UNBALANCED NORM.LT.ABS-TOLERANCE = ',1PD12.4)
 2000 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM')
 2020 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM',
     2  5X,'UNBALANCED REF RHS',5X,'UBALANCED/REF RHS LAMDA') 
 2001 FORMAT(I5,8X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4)
 2021 FORMAT(13X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4)
 2030 FORMAT(' IDOF=1 - PRESSURES; IDOF=2 - CONCENTRATIONS')   
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE KONVTE(F,DTT,TT1,IDT,IDTCD,JEDN,ITER,TREF,DTMOD,EPSTA,
     1  EPSTR,KONV,IPRINT,NPRINT,KOR,NDT,TIME,NPT,ILAMDA,KK,NDIMIDT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON/MESSAGES/MESSAGE
      COMMON /SRPSKI/ ISRPS
      COMMON/FITNODES/INFF
C 
      DIMENSION F(JEDN),TT1(JEDN),DTT(JEDN),IDT(NDIMIDT,NPT)
      DIMENSION IDTCD(2,NPT)
C
      KONV   = 1
      IDOF   = KK
	TOL    = 1.D-10
	TOLABS = 1.D-8
	milos  = 1
	if(milos.eq.1) go to 10
      IF(DABS(EPSTA).LT.TOL) GO TO 10
      IF(ITER.EQ.0) GO TO 50
      DO 20 I = 1,JEDN
      IF(DABS(F(I)).LT.EPSTA) GO TO 20
      KONV = 0
      GO TO 50
   20 CONTINUE
   50 IF(KONV.EQ.0) RETURN
      IF(DABS(EPSTR).LT.TOL) RETURN
C
   10 DTMOD  = 0.
      TMAX   = 0.
	TT1ABS = 0.D0
	if(iter.eq.0)tref=0.d0
C
      DO 40 N=1,NPT
          I   = IDT(KK,N)
          IF(I.LE.0) GO TO 40
          TMAX   = TMAX   + DTT(I)*DTT(I)
	    TT1ABS = TT1ABS + TT1(I)*TT1(I)
          DTMOD  = DTMOD  + F(I)*F(I)
   40   CONTINUE
C
      DTMOD  = DSQRT(DTMOD)
C      TMAX   = DSQRT(TMAX)
	TT1ABS = DSQRT(TT1ABS)
	TMAX=TT1ABS
	TPRINT = TT1ABS
C	IF(ITER.EQ.0) UNBREF=DTMOD
      IF(TMAX.GT.TREF) GO TO 70
      TMAX = TREF
      GO TO 80
   70 TREF = TMAX
   80 TOLT = TMAX*EPSTR
C
C     CORRECTION FOR CONVERGENCE (OCT. 16,09)
C      TORELAT=UNBREF*EPSTR
      IF(TMAX.LT.TOL) THEN
	   TPRINT=TT1ABS
         ERRORABS=TOL
	   IF(TT1ABS.GT.TOL)ERRORABS=TMAX/TT1ABS
	   ERRORREL=TOL
         GO TO 300
	ENDIF
	ERRORREL=DTMOD/TMAX
      IF(ERRORREL.GT.EPSTR) KONV = 0
C      IF(DTMOD.GT.TOLT) KONV = 0
	if((EPSTA.gt.1.d-10).and.(dtmod.le.EPSTA))konv=1
	IF(DTMOD.LT.TOL)KONV=1
C
c     dodat kriterijum za absolutnu toleranciju:TOLABS=0.05D0      
C      
      INDABS=0
	if(DTMOD.lt.TOLABS.AND.KONV.EQ.0) THEN
	  INDABS=1  
	  konv=1
	  IF(INFF.EQ.0) RETURN 
        WRITE(MESSAGE,9040) TIME,KOR,ITER,dtmod,TOLABS
 9040   FORMAT(/
     1  ' CONVERGENCE REACHED FOR INCREMENT SMALLER',
     2  ' THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1  ' TIME=',1PD12.5/
     1  ' STEP NO. =',I5/
	2  ' ITERATION NO. =',I5/
     3  ' INCERMENT =',1PD12.5/
     4  ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/)   
      ENDIF
C      ENE=DTMOD/TMAX
 300  IF(TPRINT.LT.TOL) THEN
        ERRORABS=0.D0
	ELSE
	  ERRORABS=DTMOD/TPRINT
	ENDIF
C
      IF(INFF.EQ.0) RETURN
C
      IF(ITER.EQ.0.AND.KK.EQ.1) THEN
        IF(ILAMDA.EQ.0) THEN
            WRITE(*,2000)KOR,NDT
        ELSE
            WRITE(*,2020)KOR,NDT
        ENDIF
      ENDIF
      WRITE(*,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS              
      IF(INDABS.EQ.1) THEN
        WRITE(*,2002) TOLABS
      ENDIF
	IF(NPRINT.GT.1) THEN
	  IF(KOR.GT.1.AND.KOR.LT.NDT) THEN
           IF(IPRINT.NE.(NPRINT-1)) GO TO 120
	  ENDIF    
      ENDIF
      IF(ISRPS.EQ.0) THEN
        WRITE(3,1000) ITER,DTMOD,ITER,ENE
        WRITE(*,1000) ITER,DTMOD,ITER,ENE
      ELSE
	  IF(ITER.EQ.0.AND.KK.EQ.1) THEN
	      IF(ILAMDA.EQ.0) THEN
	         WRITE(3,2000)KOR,NDT
	      ELSE
	         WRITE(3,2020)KOR,NDT
	      ENDIF
        ENDIF
        WRITE(3,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS
        IF(INDABS.EQ.1) THEN
          WRITE(3,2002) TOLABS
        ENDIF
	ENDIF
  120 CONTINUE
      RETURN
 1000 FORMAT(/'      ITERACIJA = ',I5,5X,'NORMA  = ',1PD12.4/
     1        '        NORMA (',I3,') / / NORMA (',I3,' )  = ',
     2               1PD12.4)
 2010 FORMAT(/' ITERATION =',I3,1X,': INCREMENT      = ',2X,
     1               1PD12.4,3X,'NORM     = ',1PD12.4/
     115X,'ICR/ICREMSTEP',5X,'ICR/NORM=',21PD12.4)
 2002 FORMAT(/' UNBALANCED NORM.LT.ABS-TOLERANCE = ',1PD12.4)
 2000 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM')
 2020 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM',
     2  5X,'UNBALANCED REF RHS',5X,'UBALANCED/REF RHS LAMDA') 
 2001 FORMAT(I5,8X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4)  
      END   
C=========================================================================
C
C=======================================================================
      SUBROUTINE KONVTEMIG(F,TT1,JEDN,ITER,TREF,DTMOD,EPSTA,
     1  EPSTR,KONV,IPRINT,NPRINT,KOR,NDT,TIME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON/MESSAGES/MESSAGE
      COMMON /SRPSKI/ ISRPS
C 
      DIMENSION F(JEDN),TT1(JEDN)
C
      KONV   = 1
      IDOF=1
      INFF=0
	TOL    = 1.D-10
	TOLABS = 1.D-8
C
   10 DTMOD  = 0.
      TMAX   = 0.
	TT1ABS = 0.D0
	if(iter.eq.0)tref=0.d0
C
      DO 40 I=1,JEDN
	    TT1ABS = TT1ABS + TT1(I)*TT1(I)
          DTMOD  = DTMOD  + F(I)*F(I)
   40 CONTINUE 
C
      DTMOD  = DSQRT(DTMOD)
	TT1ABS = DSQRT(TT1ABS)
	TMAX=TT1ABS
	TPRINT = TT1ABS
C	IF(ITER.EQ.0) UNBREF=DTMOD
      IF(TMAX.GT.TREF) GO TO 70
      TMAX = TREF
      GO TO 80
   70 TREF = TMAX
   80 TOLT = TMAX*EPSTR
C
C     CORRECTION FOR CONVERGENCE (OCT. 16,09)
C      TORELAT=UNBREF*EPSTR
      IF(TMAX.LT.TOL) THEN
	   TPRINT=TT1ABS
         ERRORABS=TOL
	   IF(TT1ABS.GT.TOL)ERRORABS=TMAX/TT1ABS
	   ERRORREL=TOL
         GO TO 300
	ENDIF
	ERRORREL=DTMOD/TMAX
      IF(ERRORREL.GT.EPSTR) KONV = 0
C      IF(DTMOD.GT.TOLT) KONV = 0
	if((EPSTA.gt.1.d-10).and.(dtmod.le.EPSTA))konv=1
	IF(DTMOD.LT.TOL)KONV=1
C
c     dodat kriterijum za absolutnu toleranciju:TOLABS=0.05D0      
C      
      INDABS=0
	if(DTMOD.lt.TOLABS.AND.KONV.EQ.0) THEN
	  INDABS=1  
	  konv=1
	  IF(INFF.EQ.0) RETURN 
        WRITE(MESSAGE,9040) TIME,KOR,ITER,dtmod,TOLABS
 9040   FORMAT(/
     1  ' CONVERGENCE REACHED FOR INCREMENT SMALLER',
     2  ' THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1  ' TIME=',1PD12.5/
     1  ' STEP NO. =',I5/
	2  ' ITERATION NO. =',I5/
     3  ' INCERMENT =',1PD12.5/
     4  ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/)   
      ENDIF
C      ENE=DTMOD/TMAX
 300  IF(TPRINT.LT.TOL) THEN
        ERRORABS=0.D0
	ELSE
	  ERRORABS=DTMOD/TPRINT
	ENDIF
C
      IF(ITER.EQ.0) THEN
          WRITE(*,2000)KOR,NDT
      ENDIF
      WRITE(*,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS              
      IF(INDABS.EQ.1) THEN
        WRITE(*,2002) TOLABS
      ENDIF
	IF(NPRINT.GT.1) THEN
	  IF(KOR.GT.1.AND.KOR.LT.NDT) THEN
           IF(IPRINT.NE.(NPRINT-1)) GO TO 120
	  ENDIF    
      ENDIF
	IF(ITER.EQ.0) THEN
	     WRITE(3,2000)KOR,NDT
      ENDIF
      WRITE(3,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS
      IF(INDABS.EQ.1) THEN
          WRITE(3,2002) TOLABS
      ENDIF
  120 CONTINUE
      RETURN
 2010 FORMAT(/' ITERATION =',I3,1X,': INCREMENT      = ',2X,
     1               1PD12.4,3X,'NORM     = ',1PD12.4/
     115X,'ICR/ICREMSTEP',5X,'ICR/NORM=',21PD12.4)
 2002 FORMAT(/' UNBALANCED NORM.LT.ABS-TOLERANCE = ',1PD12.4)
 2000 FORMAT (/' ITERATIONS FOR STEP ON CELL MIGRATION:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM')
 2020 FORMAT (/' ITERATIONS FOR STEP ON CELL MIGRATION:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM',
     2  5X,'UNBALANCED REF RHS',5X,'UBALANCED/REF RHS LAMDA') 
 2001 FORMAT(I5,8X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4)  
      END
C =========================================================================
C
C=========================================================================
      SUBROUTINE KONVTEMULT(F,TT1,JEDN,ITER,TREF,DTMOD,EPSTA,EPSTR,
     1  KONV,IPRINT,NPRINT,KOR,NDT,VREME,NPT,KK,MULTIMOL,MM)
C
C     CONVERGENCE CHECK FOR MULTIMOLECULAR MODEL
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON/MESSAGES/MESSAGE
      COMMON /SRPSKI/ ISRPS
      COMMON/FITNODES/INFF
C 
      DIMENSION F(JEDN),TT1(MULTIMOL,JEDN)
C
      KONV   = 1
      IDOF   = KK
	TOL    = 1.D-10
	TOLABS = 1.D-8
	milos  = 1
	if(milos.eq.1) go to 10
      IF(DABS(EPSTA).LT.TOL) GO TO 10
      IF(ITER.EQ.0) GO TO 50
      DO 20 I = 1,JEDN
      IF(DABS(F(I)).LT.EPSTA) GO TO 20
      KONV = 0
      GO TO 50
   20 CONTINUE
   50 IF(KONV.EQ.0) RETURN
      IF(DABS(EPSTR).LT.TOL) RETURN
C
   10 DTMOD  = 0.
      TMAX   = 0.
	TT1ABS = 0.D0
	if(iter.eq.0)tref=0.d0
C
      DO 40 I=1,JEDN
	    TT1ABS = TT1ABS + TT1(MM,I)*TT1(MM,I)
          DTMOD  = DTMOD  + F(I)*F(I)
   40   CONTINUE 
      DTMOD  = DSQRT(DTMOD)
C      TMAX   = DSQRT(TMAX)
	TT1ABS = DSQRT(TT1ABS)
	TMAX=TT1ABS
	TPRINT = TT1ABS
C	IF(ITER.EQ.0) UNBREF=DTMOD
      IF(TMAX.GT.TREF) GO TO 70
      TMAX = TREF
      GO TO 80
   70 TREF = TMAX
   80 TOLT = TMAX*EPSTR
C
C     CORRECTION FOR CONVERGENCE (OCT. 16,09)
C      TORELAT=UNBREF*EPSTR
      IF(TMAX.LT.TOL) THEN
	   TPRINT=TT1ABS
         ERRORABS=TOL
	   IF(TT1ABS.GT.TOL)ERRORABS=TMAX/TT1ABS
	   ERRORREL=TOL
         GO TO 300
	ENDIF
	ERRORREL=DTMOD/TMAX
      IF(ERRORREL.GT.EPSTR) KONV = 0
C      IF(DTMOD.GT.TOLT) KONV = 0
	if((EPSTA.gt.1.d-10).and.(dtmod.le.EPSTA))konv=1
	IF(DTMOD.LT.TOL)KONV=1
C
c     dodat kriterijum za absolutnu toleranciju:TOLABS=0.05D0      
C      
      INDABS=0
	if(DTMOD.lt.TOLABS.AND.KONV.EQ.0) THEN
	  INDABS=1  
	  konv=1
	  IF(INFF.EQ.0) RETURN 
        WRITE(MESSAGE,9040) TIME,KOR,ITER,dtmod,TOLABS
 9040   FORMAT(/
     1  ' CONVERGENCE REACHED FOR INCREMENT SMALLER',
     2  ' THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1  ' TIME=',1PD12.5/
     1  ' STEP NO. =',I5/
	2  ' ITERATION NO. =',I5/
     3  ' INCERMENT =',1PD12.5/
     4  ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/)   
      ENDIF
C      ENE=DTMOD/TMAX
 300  IF(TPRINT.LT.TOL) THEN
        ERRORABS=0.D0
	ELSE
	  ERRORABS=DTMOD/TPRINT
	ENDIF
C
      IF(INFF.EQ.0) RETURN
C
      IF(ITER.EQ.0.AND.KK.EQ.1) THEN
         WRITE(*,2000)KOR,NDT
      ENDIF
      WRITE(*,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS              
      IF(INDABS.EQ.1) THEN
        WRITE(*,2002) TOLABS
      ENDIF
	IF(NPRINT.GT.1) THEN
	  IF(KOR.GT.1.AND.KOR.LT.NDT) THEN
           IF(IPRINT.NE.(NPRINT-1)) GO TO 120
	  ENDIF    
      ENDIF
	  IF(ITER.EQ.0.AND.KK.EQ.1) THEN
	     WRITE(3,2000)KOR,NDT
        ENDIF
        WRITE(3,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS
        IF(INDABS.EQ.1) THEN
          WRITE(3,2002) TOLABS
        ENDIF
  120 CONTINUE
      RETURN
 2010 FORMAT(/' ITERATION =',I3,1X,': INCREMENT      = ',2X,
     1               1PD12.4,3X,'NORM     = ',1PD12.4/
     115X,'ICR/ICREMSTEP',5X,'ICR/NORM=',21PD12.4)
 2002 FORMAT(/' UNBALANCED NORM.LT.ABS-TOLERANCE = ',1PD12.4)
 2000 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM')
 2020 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM',
     2  5X,'UNBALANCED REF RHS',5X,'UBALANCED/REF RHS LAMDA') 
 2001 FORMAT(I5,8X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4)  
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE KONVTL(F,TT1,LAMDAV,JEDN,ITER,TREF,DTMOD,EPSTA,EPSTR,
	1                  KONV,IPRINT,NPRINT,KOR,NDT,TIME,FLAMREF,NPT,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	COMMON/MESSAGES/MESSAGE
      COMMON /SRPSKI/ ISRPS
      COMMON/FITNODES/INFF
      COMMON/AULCOMM/ ALUNBALANCE
      
      DIMENSION F(JEDN),TT1(JEDN),LAMDAV(NPT)
      KONV = 1
      IDOF=3
	TOL=1.D-10
	TOLABS=1.D-8
	TOLLAM=1.D-3
	KK=2
   10 DTMOD = 0.
      TMAX = 0.
	TT1ABS=0.D0
	if(iter.eq.0)tref=0.d0
      DO  NODE=1,NPT
         I=LAMDAV(NODE)
         IF(I.GT.0) THEN
            IF(IND.EQ.1) THEN
               IF(TT1(I).LT.0.D0) THEN
                  TT1(I)=TOLLAM
               ENDIF
	         TT1ABS=TT1ABS+TT1(I)*TT1(I)
               DTMOD = DTMOD + F(I)*F(I)
            ELSE
              IF(ITER.EQ.0) THEN
                 FLAMREF=FLAMREF + F(I)*F(I)
              ELSE
                DTMOD = DTMOD + F(I)*F(I)
              ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF(IND.EQ.0) THEN
         IF(ITER.EQ.0) THEN
             FLAMREF=DSQRT(FLAMREF)
             IF(FLAMREF.LT.TOL) FLAMREF=TOL
             ALUNBALANCE=FLAMREF
         ELSE 
           ALUNBALANCE=DSQRT(DTMOD)
         ENDIF
         RETURN
      ENDIF
      DTMOD = DSQRT(DTMOD)
	TT1ABS=DSQRT(TT1ABS)
	TMAX=TT1ABS
	TPRINT=TT1ABS
C	IF(ITER.EQ.0) UNBREF=DTMOD
      IF(TMAX.GT.TREF) GO TO 70
      TMAX = TREF
      GO TO 80
   70 TREF = TMAX
   80 TOLT = TMAX*EPSTR
      TOLRHS=FLAMREF*EPSTR
C
C     CORRECTION FOR CONVERGENCE (OCT. 16,09)
C      TORELAT=UNBREF*EPSTR
      IF(TMAX.LT.TOL) THEN
	   TPRINT=TT1ABS
         ERRORABS=TOL
	   IF(TT1ABS.GT.TOL)ERRORABS=TMAX/TT1ABS
	   ERRORREL=TOL
         GO TO 300
	ENDIF
	ERRORREL=DTMOD/TMAX
	IF(ITER.EQ.0) THEN
	   ERRORRELRHS=1.D0 
	ELSE
	   ERRORRELRHS=DABS(ALUNBALANCE/FLAMREF)
	ENDIF
	
      !IF(ERRORREL.GT.EPSTR.OR.ERRORRELRHS.GT.EPSTRL) KONV = 0   !miljan
      IF(ERRORRELRHS.GT.EPSTR) KONV = 0
      
C      IF(DTMOD.GT.TOLT) KONV = 0
	if((EPSTA.gt.1.d-10).and.(dtmod.le.EPSTA))konv=1
	IF(DTMOD.LT.TOL)KONV=1
c     dodat kriterijum za absolutnu toleranciju:TOLABS=0.05D0      
C      
      INDABS=0
	if(DTMOD.lt.TOLABS.AND.KONV.EQ.0) THEN
	  INDABS=1  
	  konv=1
	  IF(INFF.EQ.0) RETURN 
        WRITE(MESSAGE,9040) TIME,KOR,ITER,dtmod,TOLABS
 9040   FORMAT(/
     1  ' CONVERGENCE REACHED FOR INCREMENT SMALLER',
     2  ' THAN NUMERICAL TOLERANCE IN PAK-T  '/
	1  ' TIME=',1PD12.5/
     1  ' STEP NO. =',I5/
	2  ' ITERATION NO. =',I5/
     3  ' INCERMENT =',1PD12.5/
     4  ' NUMERICAL TOLERANCE IN PAK-T, TOLABS=',1PD12.5/)   
      ENDIF
C      ENE=DTMOD/TMAX
 300  IF(TPRINT.LT.TOL) THEN
        ERRORABS=0.D0
	ELSE
	  ERRORABS=DTMOD/TPRINT
	ENDIF
C
      IF(INFF.EQ.0) RETURN
C
      IF(ITER.EQ.0.AND.KK.EQ.1) THEN
        WRITE(*,2000)KOR,NDT
      ENDIF
        WRITE(*,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS,FLAMREF,
     1                ERRORRELRHS                
      IF(INDABS.EQ.1) THEN
        WRITE(*,2002) TOLABS
      ENDIF
	IF(NPRINT.GT.1) THEN
	  IF(KOR.GT.1.AND.KOR.LT.NDT) THEN
           IF(IPRINT.NE.(NPRINT-1)) GO TO 120
	  ENDIF    
      ENDIF
      IF(ISRPS.EQ.0) THEN
        WRITE(3,1000) ITER,DTMOD,ITER,ENE
        WRITE(*,1000) ITER,DTMOD,ITER,ENE
      ELSE
	  IF(ITER.EQ.0.AND.KK.EQ.1) THEN
	    WRITE(3,2000)KOR,NDT
        ENDIF
        WRITE(3,2001) ITER,IDOF,DTMOD,TPRINT,ERRORABS,FLAMREF,
     1                ERRORRELRHS
        IF(INDABS.EQ.1) THEN
          WRITE(3,2002) TOLABS
        ENDIF
	ENDIF
  120 CONTINUE
      RETURN
 1000 FORMAT(/'      ITERACIJA = ',I5,5X,'NORMA  = ',1PD12.4/
     1        '        NORMA (',I3,') / / NORMA (',I3,' )  = ',
     2               1PD12.4)
 2010 FORMAT(/' ITERATION =',I3,1X,': INCREMENT      = ',2X,
     1               1PD12.4,3X,'NORM     = ',1PD12.4/
     115X,'ICR/ICREMSTEP',5X,'ICR/NORM=',21PD12.4)
 2002 FORMAT(/' UNBALANCED NORM.LT.ABS-TOLERANCE = ',1PD12.4)
 2000 FORMAT (/' ITERATIONS FOR STEP:', I10,' /',I7/
     1 ' ITERATION',5X,'IDOF',5X,'INCREMENT',11X,'NORM',10X,'ICR/NORM')
 2001 FORMAT(I5,8X,I5,3X,1PD12.4,6X,1PD12.4,4X,1PD12.4,5X,1PD12.4,10X,
     1       1PD12.4)  
      END
C=======================================================================

C=======================================================================
      SUBROUTINE ELEMT(NETIP,IND)
      GO TO (100,200,300),IND
  100 GO TO (110,120,130,110),NETIP
 
  130 CALL UL3DGL
      RETURN
  120 CALL UL2DGL
      RETURN
  110 CALL UL1DGL
      RETURN
  200 GO TO (210,220,230,210),NETIP
  230 CALL K21MGL
      RETURN
  220 CALL K09MGL
      RETURN
  210 CALL K03MGL
      RETURN
  300 GO TO (310,320,330,310),NETIP
  330 CONTINUE
      RETURN
  320 CONTINUE
      RETURN
  310 CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE ULTABF(NTAV,FN,NBR,MAXT,IBR,IMAX,ITMAX,TFMX,TMAX,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*130 ACOZ
CE READ INPUT DATA FOR FUNCTIONS GIVEN AS TABLE
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
      DIMENSION NTAV(*),FN(2,NBR,*)
      DO 100 N=1,NBR
      ITMAX=1
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) IBR,IMAX
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) IBR,IMAX
      IF(IMAX.LE.MAXT) GO TO 20
      RETURN
   20 NTAV(IBR)=IMAX
      ITMAX=0
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2000) IBR,NTAV(IBR)
      IF(ISRPS.EQ.1)
     *WRITE(3,6000) IBR,NTAV(IBR)
      ENDIF
 2000 FORMAT(//11X,'FUNKCIJA BR.=',I3,5X,'BROJ TACAKA NA KRIVOJ=',I3/)
 6000 FORMAT(//11X,'FUNCTION NO.=',I3,5X,'NUMBER OF POINTS ON CURVE=',
     1 I3/)
      DO 25 J=1,IMAX
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*)(FN(I,IBR,J),I=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1001)(FN(I,IBR,J),I=1,2)
  25  CONTINUE
      IF(IND.EQ.0) GO TO 30
      IF(MAXT.EQ.1) GO TO 30
      IF(FN(1,IBR,IMAX).GE.TMAX) GO TO 30
      ITMAX=1
      TFMX=FN(1,IBR,IMAX)
      RETURN
   30 IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2001)
      IF(ISRPS.EQ.1)
     *WRITE(3,6001)
 2001 FORMAT(11X,'ARGUMENT',6X,'FUNKCIJA')
 6001 FORMAT(11X,'ARGUMENT',6X,'FUNCTION')
      WRITE(3,5002)((FN(I,IBR,J),I=1,2),J=1,IMAX)
 5002 FORMAT(10X,E10.4,5X,E10.4)
      ENDIF
  100 CONTINUE
      IBR=0
      RETURN
 1000 FORMAT(10I5)
 1001 FORMAT(2E10.4)
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE POCTE(IDT,TT1,TCONS,CINITIAL,NDIMIDT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*130 ACOZ
CE  READ INITIAL TEMPERATURES
      COMMON /TEMPGL/ NPT,NGET,NMATT,NSTAC,MAXIT,NDT,NTCONS,
     1NET,NGE,JEDN,DT,TMAX,VREME,NPOC,NWK,NPRINT,LINTE,
     2KONF,INDSK,INDF,KOR,NTOTA,ITER,INDSC,INTEB,INDAX,IANIZ,NQE,IZIP,
     3ITOR,INDVEL,INDFL
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
	common/diffusmain/indiffus,NDS,maxdt,NPRINTnds,lvreme,leprdt,
	1                  NUMMAT,maxPmat,LMATER,LMATFUN,LAREAD,LFLUXDX,
     1                  LFLUXDY,LFLUXDZ
C
      DIMENSION IDT(NDIMIDT,NPT),TT1(*),TCONS(*),T14(7),N14(14)
      IF(NPOC.GT.1) GO TO 50  
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2000)
      IF(ISRPS.EQ.1)
     *WRITE(3,6000)
      ENDIF
 6000 FORMAT(/////6X,'I N I T I A L   C O N C E N T R A T I O N S'/
     16X,39('-')//)
 2000 FORMAT(/////6X,'P O C E T N E   T E M P E R A T U R E'/
     16X,37('-')//)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) TPOC
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) TPOC
 1000 FORMAT(E10.4,2I5)
      CINITIAL=TPOC
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2001) TPOC
       IF(ISRPS.EQ.1) THEN
           WRITE(3,6002) TPOC
	  ENDIF
       ENDIF
 2001 FORMAT(///
     111X,'SVE TACKE IMAJU ISTU POCETNU TEMPERATURU ... TPOC =',1PD10.3/
     216X,'(ZA NSTAC.EQ.1)')
 6001 FORMAT(///
     111X,'ALL POINTS HAVE SAME INITIAL CONCENTRATION.. CINITAL =',
     1 E10.4/
     216X,'(FOR NSTAC.EQ.1)')
 6002 FORMAT(//' ALL POINTS HAVE INTIAL CONCENTRATION,',
     1     '  ...CINITAL =',E10.4/)
      J=0
      DO 20 I=1,NPT
C      IF(NTCONS.EQ.0) GO TO 10
       II=IDT(2,I)
       IF(II.EQ.0) GO TO 20
   10 TT1(II)=TPOC
   20 CONTINUE
 
CE  POINTS ARE WITH DIFERENT TEMPERATURES  - NOT TESTED HERE
50    IF(ISRPS.EQ.0)
     *WRITE(3,2002)
      IF(ISRPS.EQ.1)
     *WRITE(3,6012)
 2002 FORMAT(//1X,'TACKE IMAJU RAZLICITE TEMPERATURE:')
 6012 FORMAT(//1X,'POINTS ARE WITH DIFERENT TEMPERATURES:')
 
      DO 300 IN=1,NPOC
      READ(1,1000)TPOC,NN1,NNA
      IF(ISRPS.EQ.0)
     *WRITE(3,2000)TPOC,NN1,NNA
      IF(ISRPS.EQ.1)
     *WRITE(3,6000)TPOC,NN1,NNA
      IF(NN1.EQ.0) GO TO 90
      IF(ISRPS.EQ.0)
     *WRITE(3,2005)
      IF(ISRPS.EQ.1)
     *WRITE(3,6005)
 2005 FORMAT(1X,'POJEDINACNO ZADATI CVOROVI:')
 6005 FORMAT(1X,'POJEDINACNO ZADATI CVOROVI:')
   55 READ(1,1001)(N14(I),I=1,14)
      DO 70 I=1,14
      IF(N14(I).EQ.0) GO TO 90
      J=N14(I)
      IF(IDT(2,J).EQ.0) GO TO 60
      K=IDT(2,J)
      TT1(K)=TPOC
      GO TO 70
   60 IF(DABS(TCONS(J)-TPOC).LT.1.D-10) GO TO 70
      IF(ISRPS.EQ.0)
     *WRITE(3,2004) J,TCONS(J),TPOC
      IF(ISRPS.EQ.1)
     *WRITE(3,6004) J,TCONS(J),TPOC
 2004 FORMAT(/1X,'GRESKA U ULAZNIM PODACIMA:'/1X,'CVOT BROJ=',I5/
     11X,'UCITANA KONSTANTNA TEMPERATURA=',D12.4/
     21X,'UCITANA POCETNA    TEMPERATURA=',D12.4/1X,'NISU JEDNAKE'/
     31X,'ZA CVOR CIJA JE TEMPERATURA KONSTANTNA NE TREBA UCITAVATI POCE
     4TNU TEMPERATURU')
 6004 FORMAT(/1X,'GRESKA U ULAZNIM PODACIMA:'/1X,'CVOT BROJ=',I5/
     11X,'UCITANA KONSTANTNA TEMPERATURA=',D12.4/
     21X,'UCITANA POCETNA    TEMPERATURA=',D12.4/1X,'NISU JEDNAKE'/
     31X,'ZA CVOR CIJA JE TEMPERATURA KONSTANTNA NE TREBA UCITAVATI POCE
     4TNU TEMPERATURU')
      STOP
   70 CONTINUE
      GO TO 55
   90 IF(NNA.EQ.0) GO TO 300
   92 READ(1,1001)N1,KOR,N2
 1001 FORMAT(14I5)
      IF(N1.EQ.0) GO TO 300
      IF(ISRPS.EQ.0)
     *WRITE(3,2008) N1,KOR,N2
      IF(ISRPS.EQ.1)
     *WRITE(3,6008) N1,KOR,N2
 2008 FORMAT(1X,'PRVI CVOR ,N1=',I5,3X,'KORAK=',I3,3X,'POSLEDNJI CVOR',
     1'N2=',I5)
 6008 FORMAT(1X,'FIRST NODE ,N1=',I5,3X,'STEP=',I3,3X,'LAST NODE',
     1'N2=',I5)
      I14=0
      N=N1-KOR
   95 N=N+KOR
      IF(N.LE.N2) GO TO 100
      IF(ISRPS.EQ.0)
     *WRITE(3,2009)
      IF(ISRPS.EQ.1)
     *WRITE(3,6009)
 2009 FORMAT(1X,'GRESKA U ULAZNIM PODACIMA: POSLEDNJI CVOR MORA SE ',
     1'DOBITI SABIRANJEM KORAKA KOR POCEVSI OD PRVOG CVORA N1')
 6009 FORMAT(1X,'GRESKA U ULAZNIM PODACIMA: POSLEDNJI CVOR MORA SE ',
     1'DOBITI SABIRANJEM KORAKA KOR POCEVSI OD PRVOG CVORA N1')
      STOP
  100 I14=I14+1
      N14(I14)=N
      IF(N.EQ.N2) GO TO 110
      IF(I14.EQ.14) GO TO 105
      GO TO 95
  105 J14=14
      I14=0
      GO TO 115
  110 J14=I14
  115 DO 120 I=1,J14
      J=N14(I)
      IF(IDT(2,J).EQ.0) GO TO 130
      K=IDT(2,J)
      TT1(K)=TPOC
      GO TO 120
  130 IF(DABS(TCONS(J)-TPOC).LT.1.D-10) GO TO 120
      WRITE(3,2004)J,TCONS(J),TPOC
      STOP
  120 CONTINUE
      IF(J.EQ.N2) GO TO 92
      GO TO 95
  300 CONTINUE
      IF(ISRPS.EQ.0)
     *WRITE(3,2010)
      IF(ISRPS.EQ.1)
     *WRITE(3,6010)
 2010 FORMAT(/1X,'POCETNE TEMPER. I KONSTANTNE TEMPER. CVOROVA -',
     1' UCITANE I GENERISANE :'/1X,'CVOROVI I TEMPERATURE')
 6010 FORMAT(/1X,'POCETNE TEMPER. I KONSTANTNE TEMPER. CVOROVA -',
     1' UCITANE I GENERISANE :'/1X,'CVOROVI I TEMPERATURE')
      N=0
      NN=0
  205 DO 250 I=1,7
      N=N+1
      IF(N.GT.NPT) GO TO 230
      N14(I)=N
      IF(IDT(2,N).EQ.0) GO TO 210
      NN=NN+1
      T14(I)=TT1(NN)
      GO TO 250
  210 T14(I)=TCONS(N)
  250 CONTINUE
      WRITE(3,5011)(N14(I),T14(I),I=1,7)
      GO TO 205
  230 I=I-1
      IF(I.EQ.0) RETURN
      WRITE(3,5011)(N14(J),T14(J),J=1,I)
 5011 FORMAT(7(I5,1X,D11.4))
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE POCTE2(IDT,TT1,TCONS,CINITIAL,NDIMIDT,ISMEARED)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*130 ACOZ
CE  READ INITIAL TEMPERATURES
      COMMON /TEMPGL/ NPT,NGET,NMATT,NSTAC,MAXIT,NDT,NTCONS,
     1NET,NGE,JEDN,DT,TMAX,VREME,NPOC,NWK,NPRINT,LINTE,
     2KONF,INDSK,INDF,KOR,NTOTA,ITER,INDSC,INTEB,INDAX,IANIZ,NQE,IZIP,
     3ITOR,INDVEL,INDFL
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
	common/diffusmain/indiffus,NDS,maxdt,NPRINTnds,lvreme,leprdt,
	1                  NUMMAT,maxPmat,LMATER,LMATFUN,LAREAD,LFLUXDX,
     1                  LFLUXDY,LFLUXDZ
C
      DIMENSION IDT(NDIMIDT,NPT),TT1(*),TCONS(*),T14(7),N14(14)
      IF(NPOC.GT.1) GO TO 50  
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2000)
      IF(ISRPS.EQ.1)
     *WRITE(3,6000)
      ENDIF
 6000 FORMAT(/////6X,'I N I T I A L   C O N C E N T R A T I O N S'/
     16X,39('-')//)
 2000 FORMAT(/////6X,'P O C E T N E   T E M P E R A T U R E'/
     16X,37('-')//)
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) TPOC
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) TPOC
 1000 FORMAT(E10.4,2I5)
      CINITIAL=TPOC
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2001) TPOC
       IF(ISRPS.EQ.1) THEN
           WRITE(3,6002) TPOC
	  ENDIF
       ENDIF
 2001 FORMAT(///
     111X,'SVE TACKE IMAJU ISTU POCETNU TEMPERATURU ... TPOC =',1PD10.3/
     216X,'(ZA NSTAC.EQ.1)')
 6001 FORMAT(///
     111X,'ALL POINTS HAVE SAME INITIAL CONCENTRATION.. CINITAL =',
     1 E10.4/
     216X,'(FOR NSTAC.EQ.1)')
 6002 FORMAT(//' ALL POINTS HAVE INTIAL CONCENTRATION,',
     1     '  ...CINITAL =',E10.4/)
      J=0
      DO 20 I=1,NPT
C        IF(NTCONS.EQ.0) GO TO 10
         II=IDT(3,I)
         IF(II.GT.0) THEN
            TT1(II)=TPOC
         ENDIF
         II=IDT(4,I)
         IF(II.GT.0) THEN
            TT1(II)=TPOC
         ENDIF
   20 CONTINUE
 
CE  POINTS ARE WITH DIFERENT TEMPERATURES  - NOT TESTED HERE
50    IF(ISRPS.EQ.0)
     *WRITE(3,2002)
      IF(ISRPS.EQ.1)
     *WRITE(3,6012)
 2002 FORMAT(//1X,'TACKE IMAJU RAZLICITE TEMPERATURE:')
 6012 FORMAT(//1X,'POINTS ARE WITH DIFERENT TEMPERATURES:')
C
      RETURN  ! PRIVREMENO - MILJAN
C      
      IF(ISMEARED.NE.1) THEN
         WRITE(3,2500) ISMEARED
 2500    FORMAT(' NOT ALOWED DIFFERENT CONCENTRATONS FOR',
     1      ' ISMEARED.NE.1'/
     2      ' ISMEARED=',I5/
     3      ' STOP IN POCTE2 SUBROUTINE')
         STOP     
      ENDIF
      DO 300 IN=1,NPOC
      READ(1,1000)TPOC,NN1,NNA
      IF(ISRPS.EQ.0)
     *WRITE(3,2000)TPOC,NN1,NNA
      IF(ISRPS.EQ.1)
     *WRITE(3,6000)TPOC,NN1,NNA
      IF(NN1.EQ.0) GO TO 90
      IF(ISRPS.EQ.0)
     *WRITE(3,2005)
      IF(ISRPS.EQ.1)
     *WRITE(3,6005)
 2005 FORMAT(1X,'POJEDINACNO ZADATI CVOROVI:')
 6005 FORMAT(1X,'POJEDINACNO ZADATI CVOROVI:')
   55 READ(1,1001)(N14(I),I=1,14)
      DO 70 I=1,14
      IF(N14(I).EQ.0) GO TO 90
      J=N14(I)
      IF(IDT(2,J).EQ.0) GO TO 60
      K=IDT(2,J)
      TT1(K)=TPOC
      GO TO 70
   60 IF(DABS(TCONS(J)-TPOC).LT.1.D-10) GO TO 70
      IF(ISRPS.EQ.0)
     *WRITE(3,2004) J,TCONS(J),TPOC
      IF(ISRPS.EQ.1)
     *WRITE(3,6004) J,TCONS(J),TPOC
 2004 FORMAT(/1X,'GRESKA U ULAZNIM PODACIMA:'/1X,'CVOT BROJ=',I5/
     11X,'UCITANA KONSTANTNA TEMPERATURA=',D12.4/
     21X,'UCITANA POCETNA    TEMPERATURA=',D12.4/1X,'NISU JEDNAKE'/
     31X,'ZA CVOR CIJA JE TEMPERATURA KONSTANTNA NE TREBA UCITAVATI POCE
     4TNU TEMPERATURU')
 6004 FORMAT(/1X,'GRESKA U ULAZNIM PODACIMA:'/1X,'CVOT BROJ=',I5/
     11X,'UCITANA KONSTANTNA TEMPERATURA=',D12.4/
     21X,'UCITANA POCETNA    TEMPERATURA=',D12.4/1X,'NISU JEDNAKE'/
     31X,'ZA CVOR CIJA JE TEMPERATURA KONSTANTNA NE TREBA UCITAVATI POCE
     4TNU TEMPERATURU')
      STOP
   70 CONTINUE
      GO TO 55
   90 IF(NNA.EQ.0) GO TO 300
   92 READ(1,1001)N1,KOR,N2
 1001 FORMAT(14I5)
      IF(N1.EQ.0) GO TO 300
      IF(ISRPS.EQ.0)
     *WRITE(3,2008) N1,KOR,N2
      IF(ISRPS.EQ.1)
     *WRITE(3,6008) N1,KOR,N2
 2008 FORMAT(1X,'PRVI CVOR ,N1=',I5,3X,'KORAK=',I3,3X,'POSLEDNJI CVOR',
     1'N2=',I5)
 6008 FORMAT(1X,'FIRST NODE ,N1=',I5,3X,'STEP=',I3,3X,'LAST NODE',
     1'N2=',I5)
      I14=0
      N=N1-KOR
   95 N=N+KOR
      IF(N.LE.N2) GO TO 100
      IF(ISRPS.EQ.0)
     *WRITE(3,2009)
      IF(ISRPS.EQ.1)
     *WRITE(3,6009)
 2009 FORMAT(1X,'GRESKA U ULAZNIM PODACIMA: POSLEDNJI CVOR MORA SE ',
     1'DOBITI SABIRANJEM KORAKA KOR POCEVSI OD PRVOG CVORA N1')
 6009 FORMAT(1X,'GRESKA U ULAZNIM PODACIMA: POSLEDNJI CVOR MORA SE ',
     1'DOBITI SABIRANJEM KORAKA KOR POCEVSI OD PRVOG CVORA N1')
      STOP
  100 I14=I14+1
      N14(I14)=N
      IF(N.EQ.N2) GO TO 110
      IF(I14.EQ.14) GO TO 105
      GO TO 95
  105 J14=14
      I14=0
      GO TO 115
  110 J14=I14
  115 DO 120 I=1,J14
      J=N14(I)
      IF(IDT(2,J).EQ.0) GO TO 130
      K=IDT(2,J)
      TT1(K)=TPOC
      GO TO 120
  130 IF(DABS(TCONS(J)-TPOC).LT.1.D-10) GO TO 120
      WRITE(3,2004)J,TCONS(J),TPOC
      STOP
  120 CONTINUE
      IF(J.EQ.N2) GO TO 92
      GO TO 95
  300 CONTINUE
      IF(ISRPS.EQ.0)
     *WRITE(3,2010)
      IF(ISRPS.EQ.1)
     *WRITE(3,6010)
 2010 FORMAT(/1X,'POCETNE TEMPER. I KONSTANTNE TEMPER. CVOROVA -',
     1' UCITANE I GENERISANE :'/1X,'CVOROVI I TEMPERATURE')
 6010 FORMAT(/1X,'POCETNE TEMPER. I KONSTANTNE TEMPER. CVOROVA -',
     1' UCITANE I GENERISANE :'/1X,'CVOROVI I TEMPERATURE')
      N=0
      NN=0
  205 DO 250 I=1,7
      N=N+1
      IF(N.GT.NPT) GO TO 230
      N14(I)=N
      IF(IDT(2,N).EQ.0) GO TO 210
      NN=NN+1
      T14(I)=TT1(NN)
      GO TO 250
  210 T14(I)=TCONS(N)
  250 CONTINUE
      WRITE(3,5011)(N14(I),T14(I),I=1,7)
      GO TO 205
  230 I=I-1
      IF(I.EQ.0) RETURN
      WRITE(3,5011)(N14(J),T14(J),J=1,I)
 5011 FORMAT(7(I5,1X,D11.4))
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE POCTEMULTMOL(IDT,IDTV,IDTVCELL,TT1,TT1C,TT1CM,
     1     IACTIVECELL,ICINITV,CINITV,NDIMIDT,ISMEARED,NDIMCELL,
     2     NDIMCELL2,MULTIMOL,JEDNMAX,JEDN,NPT,MM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE    ASSIGN INITIAL CONCENTRATIONS TO NODES FOR MULTIMOLAR CASE,
C     INCLUDING CELLS  
C
      DIMENSION IDT(NDIMIDT,NPT),IDTV(MULTIMOL,NDIMIDT,NPT),
     5 IDTVCELL(MULTIMOL,NDIMCELL,NPT),TT1(JEDN),TT1C(JEDN),
     7 TT1CM(MULTIMOL,JEDNMAX),IACTIVECELL(4)
      DIMENSION ICINITV(10,NDIMCELL2),CINITV(10,NDIMCELL2)
      DIMENSION CC(12),JEDNLOC(12)
      IF(MM.EQ.1) THEN
         WRITE(3,2000)
 2000    FORMAT(/6X,'I N I T I A L  CONCENTRATIONS'/
     1  ' MOLECULE MM',5X,'PIPE  EX.CELL SPACE AND CELLS'/
     2  8X,'NODE',5X,'EQUATIONS',5X,'CONCENTRATIONS')
      ENDIF   
      DO 20 N=1,NPT
         DO I=1,2
            IF(ICINITV(MM,I).GT.0) THEN
               IF(MULTIMOL.EQ.1) THEN
                  NN=IDT(2+I,N)
                  IF(ISMEARED.EQ.1) THEN
                     IF(NN.GT.0) TT1(NN)=CINITV(MM,I)
                  ELSE
                     IF(NN.GT.0) TT1C(NN)=CINITV(MM,I)
                  ENDIF
               ELSE
                  NN=IDTV(MM,2+I,N)
                  IF(NN.GT.0)TT1CM(MM,NN)= CINITV(MM,I) 
               ENDIF
            ENDIF
         ENDDO
         IF(NDIMCELL.EQ.0) GO TO 20 
         DO J=1,NDIMCELL
               IF(ICINITV(MM,2+J).GT.0) THEN
                  NN=IDTVCELL(MM,J,N)
                  IF(MULTIMOL.EQ.1) THEN
                     IF(ISMEARED.EQ.1) THEN
                        IF(NN.GT.0) TT1(NN)=CINITV(MM,2+J)
                     ELSE
                        IF(NN.GT.0) TT1C(NN)=CINITV(MM,2+J)
                     ENDIF            
                  ELSE
                     IF(NN.GT.0)TT1CM(MM,NN)= CINITV(MM,2+J)
                  ENDIF
               ENDIF   
         ENDDO
   20 CONTINUE
      WRITE(3,2010)MM
 2010 FORMAT(I5)
      DO N=1,NPT
         CALL CLEAR(CC,12)
         CALL ICLEAR(JEDNLOC,12)
         DO I=1,2
            IF(ICINITV(MM,I).GT.0) THEN
                IF(MULTIMOL.EQ.1) THEN
                  NN=IDT(2+I,N)
                  IF(ISMEARED.EQ.1) THEN
                     IF(NN.GT.0) CC(I)=TT1(NN)
                  ELSE
                     IF(NN.GT.0) CC(I)=TT1C(NN)
                  ENDIF
               ELSE
                  NN=IDTV(MM,2+I,N)
                  IF(NN.GT.0) CC(I)=TT1CM(MM,NN)
               ENDIF
               JEDNLOC(I)=NN 
            ENDIF
         ENDDO        
         IF(NDIMCELL.GT.0) THEN
            DO J=1,NDIMCELL
               IF(ICINITV(MM,2+J).GT.0) THEN
                  NN=IDTVCELL(MM,J,N)
                  IF(NN.GT.0) THEN
                      JEDNLOC(2+J)=NN
                      IF(MULTIMOL.EQ.1) THEN
                         IF(ISMEARED.EQ.1) THEN
                            IF(NN.GT.0) CC(2+J)=TT1(NN)
                         ELSE
                            IF(NN.GT.0) CC(2+J)=TT1C(NN)
                         ENDIF
                      ELSE
                         CC(2+J)=TT1CM(MM,NN)
                      ENDIF
                   ENDIF
               ENDIF   
            ENDDO
         ENDIF
         WRITE(3,2020)N,(JEDNLOC(J),J=1,2+NDIMCELL)
 2020    FORMAT(I10,12I5)
         WRITE(3,2021)(CC(J),J=1,2+NDIMCELL)
 2021    FORMAT(10X,12F10.4)          
      ENDDO
C 
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INDST(ISKC,KSK,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*130 ACOZ
CE  READ AND GENERATE INDICATORS FOR PRINTOUT NODES
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
      DIMENSION ISKC(*),KSK(*)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2007)
      IF(ISRPS.EQ.1)
     *WRITE(3,6007)
      ENDIF
 2007 FORMAT(/////'1',
     1' UCITANI ZELJENI CVOROVI ZA KOJE SE STAMPAJU REZULTATI'/' ',53('-
     1')///
     1'  (UCITAVANJE SE ZAVRSAVA PRAZNOM KARTICOM)'//
     1' CVOR BROJ',10X,'KORAK')
 6007 FORMAT(/////'1',
     1' NODES FOR TEMERATURES PRINT'/' ',27('-')///
     1' NODE NUMB.',10X,'STEP')
      NAUT=0
    5 CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(1,*) N,KORC
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) N,KORC
 1000 FORMAT(14I5)
      ISKC(N)=1
      IF(N.EQ.0) GO TO 22
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(3,1008) N,KORC
 1008 FORMAT(I5,12X,I5)
      IF(NAUT.GT.0) GO TO 30
      IF(KORC.GT.0) GO TO 20
      GO TO 5
   20 NAUT=1
      N1=N
      KORA=KORC
      GO TO 5
   30 N2=N
      IF(N1.GT.N2) KORA=-KORA
      IAUT=(N2-N1)/KORA-1
      N=N1
      DO 34 J=1,IAUT
      N=N+KORA
   34 ISKC(N)=ISKC(N1)
      NAUT=0
      IF(KORC.EQ.0) GO TO 5
      KORA=KORC
      NAUT=1
      N1=N2
      GO TO 5
   22 IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      IF(ISRPS.EQ.0)
     *WRITE(3,2020)
      IF(ISRPS.EQ.1)
     *WRITE(3,6020)
      ENDIF
 2020 FORMAT(//////
     1' GENERISANI ZELJENI CVOROVI ZA KOJE SE STAMPAJU REZULTATI'/' ',56
     1('-')///)
 6020 FORMAT(//////
     1' GENERATED NODES FOR TEMPERATURES PRINTING               '/' ',41
     1('-')///)
      J=0
      DO 23 I=1,NPT
      IF(ISKC(I).EQ.0) GO TO 23
      J=J+1
      KSK(J)=I
   23 CONTINUE
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(3,1000) (KSK(I),I=1,J)
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE MATSTE(SK,MAXA,F,SKE,FE,LM,NCV,INDSK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
CE  INSERTING ELEMENT MATRIXES AND VECTORS
CE  INTO SYSTEM MATRIXES AND VECTORS
      DIMENSION SK(*),MAXA(*),F(*),SKE(*),FE(*),LM(*)
      K=0
      DO 200 I=1,NCV
      IVR=LM(I)
      IF(INDSK.EQ.0) GO TO 110
      DO 100 J=I,NCV
      K=K+1
      KOL=LM(J)
      IF(IVR.LE.0.OR.KOL.LE.0) GO TO 100
      IF(IVR-KOL) 10,10,20
   10 KS=MAXA(KOL) + KOL - IVR
      GO TO 50
   20 KS=MAXA(IVR)+IVR-KOL
   50 SK(KS)=SK(KS)+SKE(K)
  100 CONTINUE
  110 IF(IVR.LE.0) GO TO 200
      F(IVR)=F(IVR)+FE(I)
  200 CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE ANALKF(MAT,TEMP,AKX,AKY,AKZ,TOPM,TMNM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
CE  FUNCTION CONDUCTION COEFICIENT-TEMPERATURE GIVEN ANALITICAL
      DIMENSION TOPM(*),TMNM(*)
      GO TO (10,20,30,40,50,60,70,80,90,100),MAT
   10 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=-471.87+3.14*TEMP-4.83D-3*TEMP**2+2.38D-7*TEMP**3
      AKY=AKX
      AKZ=AKX
      RETURN
   20 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=71.75-0.08*TEMP+7.08D-5*TEMP**2-3.86D-8*TEMP**3
      AKY=AKX
      AKZ=AKX
      RETURN
   30 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=404.67-0.08*TEMP+5.03D-5*TEMP**2-3.81D-8*TEMP**3
      AKY=AKX
      AKZ=AKX
      RETURN
   40 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=204.79-0.27*TEMP+1.42D-4*TEMP**2-2.50D-8*TEMP**3
      AKY=AKX
      AKZ=AKX
      RETURN
   50 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=459.67-0.21*TEMP+2.60D-5*TEMP**2+1.40D-9*TEMP**3
      AKY=3.23-2.83D-3*TEMP+1.13D-6*TEMP**2-1.44D-10*TEMP**3
      AKZ=AKX
      RETURN
   60 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=3.23-2.83D-3*TEMP+1.13D-6*TEMP**2-1.44D-10*TEMP**3
      AKY=459.67-0.21*TEMP+2.60D-5*TEMP**2+1.40D-9*TEMP**3
      AKZ=AKY
      RETURN
   70 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=142.65-0.01*TEMP-1.04D-5*TEMP**2+1.81D-9*TEMP**3
      AKY=AKX
      AKZ=AKX
      RETURN
   80 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=169-0.14*TEMP+6.1D-5*TEMP*TEMP
      AKY=169-0.14*TEMP+6.1D-5*TEMP*TEMP
      AKZ=169-0.14*TEMP+6.1D-5*TEMP*TEMP
      RETURN
   90 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=442.0-0.31*TEMP+1.39D-4*TEMP*TEMP-3.43D-8*TEMP*TEMP*TEMP
      AKY=2.81-3.29D-3*TEMP+2.08D-6*TEMP*TEMP-4.64D-10*TEMP*TEMP*TEMP
      AKZ=AKX
      RETURN
  100 IF(TEMP.GT.TOPM(MAT)) GO TO 1000
      IF(TEMP.LT.TMNM(MAT)) GO TO 2000
      AKX=2.81-3.29D-3*TEMP+2.08D-6*TEMP*TEMP-4.64D-10*TEMP*TEMP*TEMP
      AKY=442.0-0.31*TEMP+1.39D-4*TEMP*TEMP-3.43D-8*TEMP*TEMP*TEMP
      AKZ=AKX
      RETURN
 1000 WRITE(3,1001) MAT,TEMP,TOPM(MAT)
 1001 FORMAT(/' PREKORACENA MAXIMALNA VREDNOST DOZVOLJENE TEMPERATURE'/
     1' U MATERIJALU BROJ, MAT=',I5/' TEKUCA TEMPERATURA, TEMP=',D14.5/
     2' MAXIMALNA DOZVOLJENA TEMPERATURA U MATERIJALU, TOPM=',D14.5)
      STOP
 2000 WRITE(3,2001) MAT,TEMP,TMNM(MAT)
 2001 FORMAT(/' PREKORACENA MINIMALNA VREDNOST DOZVOLJENE TEMPERATURE'/
     1' U MATERIJALU BROJ, MAT=',I5/' TEKUCA TEMPERATURA, TEMP=',D14.5/
     2' MANIMALNA DOZVOLJENA TEMPERATURA U MATERIJALU, TMNM=',D14.5)
      STOP
      END
C=======================================================================
C=======================================================================
      SUBROUTINE ANALCF(MAT,TEMP,CC,GUSM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
CE  FUNCTION SPECIFIC HEAT-TEMPERATURE GIVEN ANALITICAL
      DIMENSION GUSM(*)
      GO TO (10,20,30,40,50,60,70,80,90,100),MAT
   10 CP=1016.48-1.17*TEMP+2.53D-3*TEMP**2-7.49D-7*TEMP**3
      CC=CP*GUSM(MAT)
      RETURN
   20 CP=-85.28+2.92*TEMP-4.85D-3*TEMP**2+3.02D-6*TEMP**3
      CC=CP*GUSM(MAT)
      RETURN
   30 CP=329.33+0.15*TEMP-3.68D-5*TEMP**2+1.60D-8*TEMP**3
      CC=CP*GUSM(MAT)
      RETURN
   40 CP=416.11+2.14*TEMP-1.02D-3*TEMP**2+1.70D-7*TEMP**3
      CC=CP*GUSM(MAT)
      RETURN
   50 CP=416.11+2.14*TEMP-1.02D-3*TEMP**2+1.70D-7*TEMP**3
      CC=CP*GUSM(MAT)
      RETURN
   60 CP=416.11+2.14*TEMP-1.02D-3*TEMP**2+1.70D-7*TEMP**3
      CC=CP*GUSM(MAT)
      RETURN
   70 CP=189.63+0.24*TEMP-1.05D-4*TEMP**2+1.59D-8*TEMP**3
      CC=CP*GUSM(MAT)
      RETURN
   80 CP=133.6+0.02*TEMP-7.1D-6*TEMP*TEMP
      CC=CP*GUSM(MAT)
      RETURN
   90 CP=755.6+2.29*TEMP-1.59D-3*TEMP*TEMP+3.94D-7*TEMP*TEMP*TEMP
      CC=CP*GUSM(MAT)
      RETURN
  100 CP=755.6+2.29*TEMP-1.59D-3*TEMP*TEMP+3.94D-7*TEMP*TEMP*TEMP
      CC=CP*GUSM(MAT)
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE TESTLNEL(NEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NEL(NE2D,4)
      
      COMMON /EL2K09/NMAT2D,MAT2D,NE2D,NGAU2X,NGAU2Y,NP2DMX,
     1IPR2DC,IPTG2,N2DGU,MEL2D,IQE2,NBR2,NCV2,LM2C(9),IDEL2(5),
     2NOD5(5),NND5,IELX2,NLM2,NTHIC,THIC
      
      DO 200 NLM2= 1,NE2D
      DO 50 I=1,NBR2
      J=NEL(NLM2,I)
   50 CONTINUE 
  200 CONTINUE
  
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ULTABFSUB(NTAV,FN,NBR,MAXT,TMAX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*100 ACOZ
CE READ INPUT DATA FOR FUNCTIONS GIVEN AS TABLE
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
      DIMENSION NTAV(NBR),FN(2,NBR,MAXT)
      DO 100 N=1,NBR
      CALL ISPITA(ACOZ)
      READ(ACOZ,1000) IBR,IMAX
      IF(IMAX.GT.MAXT) THEN
         WRITE(3,2000)IBR,IMAX,MAXT
2000     FORMAT(/' ERROR IN TIME FUNCTION: NUMBER OF POINTS GT.MAXIMUM'/
     1    ' FUNCTION NUMBER.........  IBR =',I5/
     1    ' NUMBER OF POINTS....      IMAX=',I5/
     2    ' MAX NUMBER OF POINTS....  MAXT=',I5/
     3    ' S T O P')
         STOP
      ENDIF 
   20 NTAV(IBR)=IMAX
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
         WRITE(3,6000) IBR,NTAV(IBR)
      ENDIF
 6000 FORMAT(//11X,'FUNCTION NO.=',I3,5X,'NUMBER OF POINTS ON CURVE=',
     1 I3/)
      DO 25 J=1,IMAX
         CALL ISPITA(ACOZ)
         READ(ACOZ,1001)(FN(I,IBR,J),I=1,2)
  25  CONTINUE
      IF(FN(1,IBR,IMAX).LE.TMAX.AND.IMAX.GT.1)THEN
         WRITE(3,2001)IBR,FN(1,IBR,IMAX),TMAX
2001     FORMAT(/' ERROR IN TIME FUNCTION: FUNCTION NOT DEFINED IN THE'/
     1    ' ENTIRE TIME INTERVAL'/ 
     1    ' FUNCTION NUMBER.........               IBR =',I5/
     1    ' LAST VALUE OF TIME-ARGUMENT  FN(1,IBR,IMAX)=',D12.4/
     2    ' TIME INTERVAL.......................   TMAX=',D12.4/
     3    ' S T O P')
         STOP
      ENDIF
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
         WRITE(3,6001)
 6001    FORMAT(11X,'ARGUMENT',6X,'FUNCTION')
         WRITE(3,5002)((FN(I,IBR,J),I=1,2),J=1,IMAX)
 5002    FORMAT(10X,E10.4,5X,E10.4)
      ENDIF
  100 CONTINUE
      RETURN
 1000 FORMAT(10I5)
 1001 FORMAT(2E10.4)
      END
C===========================================================================
C
C===========================================================================