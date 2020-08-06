C   14/06/90 102061623  MEMBER NAME  TMP      (#F11JET)     FVS
      SUBROUTINE VGDAT
C
C  INITIAL SETTING OF THE IO CHANNELS
C  ALL INPUT CHANNELS ARE SET TO 5 ( STANDARD INPUT ) FOR VEGAS AND INPL
C  ALL OUTPUT CHANNELS ARE SET TO 6 ( STANDARD OUTPUT ) FOR VEGAS AND IN
C  INPUT AND OUTPUT CHANNELS FOR SAVE AND RESTR ARE SET TO 7
C  INPUT AND OUTPUT CHANNELS FOR SAVE2 AND RESTR2 ARE SET TO 9
C  INPUT AND OUTPUT CHANNELS FOR VGSAVE AND VGRSTR ARE SET TO 8
C  INPUT AND OUTPUT CHANNELS FOR THE RANDOM NUMBER GENERATOR STATE IS SE
C  AN ADDITIONAL IO CHANNEL IS SET TO 18 ( SO THIS CHANNEL IS RESERVED F
C                                         ALSO )
C
C  AUTHOR : S. DE JONG
C
      COMMON/VGASIO/NINP,NOUTP
      COMMON/VGPLIO/NIN,NOUT
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
      LOGICAL FIRST
C
      DATA FIRST /.TRUE./
C
      IF(FIRST)THEN
         FIRST=.FALSE.
         NINP =5
         NIN  =5
CCC - I02GAU
CCC      NOUTP=6
CCC      NOUT =6
         NOUT = NOUTP
CCC - I02GAU
         LUN1 =7
         LUN2 =9
         LUN3 =8
         LUN4 =17
         LUN5 =18
      ENDIF
C
      RETURN
C
      END
      SUBROUTINE VEGAS(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH,IENTRY)
C
C  SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG*N
C    - BY G.P.  LEPAGE  SEPT 1976/ (REV) APR 1978
C
C    -FTN5 VERSION 21-8-1984
C    -HBOOK/HPLOT INTERFACE 6-1-1985
C
C  AUTHOR                                       : G. P. LEPAGE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FXN
      DOUBLE PRECISION ZZF1,ZZW
      DIMENSION XIN(50),R(50),DX(10),IA(10),KG(10),DT(10)
      DIMENSION XL(10),XU(10),QRAN(10),X(10)
      COMMON/VGASIO/NINP,NOUTP
      COMMON/VGB2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     + ,D(50,10),DI(50,10)
      COMMON/VGRES/S1,S2,S3,S4
      DOUBLE PRECISION S1,S2,S3,S4
CCC I02GAU THIS NUMBER IS NEEDED.
      COMMON /CCALL/ CALLS,ITMXCOM
C --- CHANGES ACCORDING TO MANUAL BY G.P.LEPAGE
      COMMON /BVEG3/ ALPH,NDMX,MDS
CCC I02GAU
C
C
      DATA XL,XU/10*0.,10*1./
CCC I02GAU @@@
      DATA NDMX/50/,ALPH/1.5/,ONE/1./,MDS/1/
CCC I02GAU @@@
C
C
CCC I02GAU @@@
C
C --- SWITCH FOR VEGAS-ENTRIES (IF ENTRY POINTS DON'T WORK)
      IF (IENTRY .EQ. 0) GOTO 7890
      GOTO (2001,2002,2003) IENTRY
7890  CONTINUE
C      NDMX=50
C      ALPH=1.5        !!! C A U T I O N !!! 3-JET-MC
C      ONE=1.
C      MDS=1
CCC I02GAU @@@
C
C      WRITE(6,*) ACC
C      WRITE(6,*) NDIM
C      WRITE(6,*) NCALL
C      WRITE(6,*) ITMX
C      WRITE(6,*) NPRN
C      WRITE(6,*) IGRAPH
C      WRITE(6,*) IENTRY

      CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'0VEGAS CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
      IF(NPRN.GT.0)THEN
         IPR=0
      ELSE
         IPR=1
      ENDIF
      NDO=1
C      WRITE(6,*) 'VEGAS HI',NDIM
      DO 1 J=1,NDIM
         XI(1,J)=ONE
1     CONTINUE
C
 2001 CONTINUE
C
C      ENTRY VEGAS1(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C...INITIALIZES CUMMULATIVE VARIABLES,BUT NOT GRID
      CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'0VEGAS1 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
      IF(MOD(IGRAPH,100) .GT. 0 )THEN
CCC I02GAU
C        IF(MOD(IGRAPH,10) .GT. 0)THEN
C           NOW=MOD(IGRAPH,100)/10
C        ELSE
C           NOW=MOD(IGRAPH,10)
C        ENDIF
         IF (IGRAPH .GE. 1) THEN
            NOW=1
         ELSE
            NOW=0
         ENDIF
         F1=0.D0
         W=0.D0
CCC I02GAU
         ZZF1  = SNGL(F1)
         ZZW   = SNGL(W)
         CALL INPLOT(NOW,ZZF1,ZZW,0)
      ENDIF
      IF(IGRAPH .GE. 100) THEN
         NOW=MOD(IGRAPH,1000)/100
C        CALL HVBOOK(NOW,F1,W)
         CONTINUE
      ENDIF
C
      IT=0
      SI=0.
      SI2=SI
      SWGT=SI
      SCHI=SI
      SCALLS=SI
C
 2002 CONTINUE
C      ENTRY VEGAS2(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C...NO INITIALIZATION
      CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'0VEGAS2 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
      ND=NDMX
      NG=1
      IF(MDS.NE.0)THEN
         NG=(NCALL/2.)**(1./NDIM)
         MDS=1
         IF((2*NG-NDMX).GE.0)THEN
            MDS=-1
            NPG=NG/NDMX+1
            ND=NG/NPG
            NG=NPG*ND
         ENDIF
      ENDIF
CDGR
      WRITE(6,*) 'MDS=',MDS
CDGR
C
      K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2) NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=DXG**(2*NDIM)/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE
      DO 3 J=1,NDIM
         DX(J)=XU(J)-XL(J)
         XJAC=XJAC*DX(J)
3     CONTINUE
C
C  REBIN PRESERVING BIN DENSITY
C
      IF(ND.NE.NDO)THEN
         RC=NDO/XND
         DO 7 J=1,NDIM
            K=0
            XN=0.
            DR=XN
            I=K
4           K=K+1
            DR=DR+ONE
            XO=XN
            XN=XI(K,J)
5           IF(RC.GT.DR) GO TO 4
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR
            IF(I.LT.NDM) GO TO 5
            DO 6  I=1,NDM
               XI(I,J)=XIN(I)
6           CONTINUE
            XI(ND,J)=ONE
7        CONTINUE
         NDO=ND
      ENDIF
C
       IF(NPRN.NE.0.AND.NPRN.NE.10)WRITE(NOUTP,200)NDIM,CALLS,IT,ITMX
     + ,ACC,MDS,ND
       IF(NPRN.EQ.10)WRITE(NOUTP,290)NDIM,CALLS,ITMX,ACC,MDS,ND
C
 2003 CONTINUE
C      ENTRY VEGAS3(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C     - MAIN INTEGRATION LOOP
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'0VEGAS3 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
9     CONTINUE
      IT=IT+1
      TI=0.
      TSI=TI
      IF(MOD(IGRAPH,100) .GT. 0 )THEN
         ZZF1  = SNGL(F1)
         ZZW   = SNGL(W)
         CALL INPLOT(NOW,ZZF1,ZZW,1)
      ENDIF
      IF(IGRAPH .GE. 100) THEN
C        CALL HVRSET(NOW,F1,W)
         CONTINUE
      ENDIF
C
      DO 10 J=1,NDIM
         KG(J)=1
         DO 10 I=1,ND
            D(I,J)=TI
            DI(I,J)=TI
10    CONTINUE
C
11    FB=0.
      F2B=FB
      K=0
C
12    CONTINUE
      K=K+1
      DO 121 J=1,NDIM
         QRAN(J)=VGRAN(0.0)
121   CONTINUE
      WGT=XJAC
      DO 15 J=1,NDIM
         XN=(KG(J)-QRAN(J))*DXG+ONE
         IA(J)=XN
         IAJ=IA(J)
         IAJ1=IAJ-1
         IF(IAJ.LE.1)THEN
            XO=XI(IAJ,J)
            RC=(XN-IAJ)*XO
         ELSE
            XO=XI(IAJ,J)-XI(IAJ1,J)
            RC=XI(IAJ1,J)+(XN-IAJ)*XO
         ENDIF
         X(J)=XL(J)+RC*DX(J)
         WGT=WGT*XO*XND
15    CONTINUE
C
CCC I02GAU --> NOW THE WEIGHT IS AVAILABLE IN THE INTEGRATED FUNCTION...
CCC   F=FXN(X)*WGT
CC      WRITE(6,*) 'NPG,K=',NPG,K
CCC I02GAU
      ITMXCOM=ITMX
      F=FXN(X,WGT)*WGT
CCC I02GAU
      F1=F/CALLS
      W=WGT/CALLS
      IF(MOD(IGRAPH,100) .GT. 0 )THEN
         ZZF1  = SNGL(F1)
         ZZW   = SNGL(W)
         CALL INPLOT(NOW,ZZF1,ZZW,2)
      ENDIF
      IF(IGRAPH .GE. 100) THEN
C        CALL HVFILL(NOW,F1,W)
         CONTINUE
      ENDIF
C
      F2=F*F
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
         IAJ=IA(J)
         DI(IAJ,J)=DI(IAJ,J)+F/CALLS
         IF(MDS.GE.0)  D(IAJ,J)=D(IAJ,J)+F2
16    CONTINUE
      IF(K.LT.NPG) GO TO 12
C
      F2B=F2B*NPG
      F2B=DSQRT(F2B)
      F2B=DABS((F2B-FB)*(F2B+FB))
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.LT.0)THEN
         DO 17 J=1,NDIM
            IAJ=IA(J)
            D(IAJ,J)=D(IAJ,J)+F2B
17       CONTINUE
      ENDIF
      K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
C
C  FINAL RESULTS FOR THIS ITERATION
C
      TI=TI/CALLS
      TSI=TSI*DV2G
      TI2=TI*TI
      IF(TSI .EQ. 0.)THEN
         WGT = 0.
      ELSE
         WGT=TI2/TSI
      ENDIF
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      IF(SWGT .EQ. 0.)THEN
         AVGI=TI
      ELSE
         AVGI=SI/SWGT
      ENDIF
      IF(SI2 .EQ. 0.)THEN
         SD=TSI
      ELSE
         SD=SWGT*IT/SI2
      ENDIF
      SCALLS=SCALLS+CALLS
      CHI2A=0.
      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
      IF(SD .NE. 0.)THEN
         SD=ONE/SD
         SD=DSQRT(SD)
      ELSE
         SD=TSI
      ENDIF
      IF(NPRN.NE.0)THEN
         TSI=DSQRT(TSI)
CCCI02GAU
C        IF(NPRN.NE.10)WRITE(NOUTP,201)IPR,IT,TI,TSI,AVGI,SD,CHI2A

         IF (NPRN.NE.10) THEN
C           IF (IT .EQ. 1) WRITE(NOUTP,3218)
            WRITE(NOUTP,201)0  ,IT,TI,TSI,AVGI,SD,CHI2A
3218     FORMAT('1')
         ENDIF
CCCI02GAU
         IF(NPRN.EQ.10)WRITE(NOUTP,203)IT,TI,TSI,AVGI,SD,CHI2A
         IF(NPRN.LT.0)THEN
            DO 20 J=1,NDIM
               WRITE(NOUTP,202)J
               WRITE(NOUTP,204)(XI(I,J),DI(I,J),D(I,J),I=1,ND)
20         CONTINUE
         ENDIF
      ENDIF
C
C   REFINE GRID
C
CC      WRITE(6,*) 'SD, AVGI=', SD, AVGI
21    IF(SD .NE. 0.)THEN
         REL = DABS(SD/AVGI)
      ELSE
         REL = 0.
      ENDIF
CC      WRITE(6,*) REL, ACC, IT, ITMX
CC      WRITE(6,*) AVGI, SD, TI, TSI, IGRAPH
      IF(REL.LE.DABS(ACC).OR.IT.GE.ITMX)NOW=2
      S1=AVGI
      S2=SD
      S3=TI
      S4=TSI
      IGRPH=MOD(IGRAPH,100)
      IF(IGRPH .GT. 0 .AND. IGRPH .LT. 10)THEN
         ZZF1  = SNGL(F1)
         ZZW   = SNGL(W)
         CALL INPLOT(NOW,ZZF1,ZZW,3)
      ELSE IF(IGRPH .GE. 10 )THEN
         ZZF1  = SNGL(F1)
         ZZW   = SNGL(W)
         CALL INPLOT(NOW,ZZF1,ZZW,4)
      ENDIF
      IF(IGRAPH .GE. 100) THEN
C        CALL HVEDIT(NOW,F1,W)
         CONTINUE
      ENDIF
C
      DO 23 J=1,NDIM
         XO=D(1,J)
         XN=D(2,J)
         D(1,J)=(XO+XN)/2.
         DT(J)=D(1,J)
         DO 22 I=2,NDM
            D(I,J)=XO+XN
            XO=XN
            XN=D(I+1,J)
            D(I,J)=(D(I,J)+XN)/3.
            DT(J)=DT(J)+D(I,J)
22       CONTINUE
         D(ND,J)=(XN+XO)/2.
         DT(J)=DT(J)+D(ND,J)
23    CONTINUE
C
      DO 28 J=1,NDIM
         RC=0.
         DO 24 I=1,ND
            R(I)=0.
            IF(D(I,J).GT.0.)THEN
               XO=DT(J)/D(I,J)
               R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
            ENDIF
            RC=RC+R(I)
24       CONTINUE
         RC=RC/XND
         K=0
         XN=0.
         DR=XN
         I=K
25       K=K+1
         DR=DR+R(K)
         XO=XN
         XN=XI(K,J)
26       IF(RC.GT.DR) GO TO 25
         I=I+1
         DR=DR-RC
         IF(DR .EQ. 0.)THEN
            XIN(I)=XN
         ELSE
            XIN(I)=XN-(XN-XO)*DR/R(K)
         ENDIF
         IF(I.LT.NDM) GO TO 26
         DO 27 I=1,NDM
            XI(I,J)=XIN(I)
27       CONTINUE
         XI(ND,J)=ONE
28    CONTINUE
C
      IF(IT.LT.ITMX.AND.DABS(ACC).LT.REL)GO TO 9
C
      S1=AVGI
      S2=SD
      S3=CHI2A
      RETURN
C
199   FORMAT(A)
200   FORMAT('0INPUT PARAMETERS FOR VEGAS   NDIM=',I3
     +,'  NCALL=',F9.0/28X,'  IT=',I5,'  ITMX =',I5/28X
     +,'  ACC=',D9.3/28X,'  MDS=',I3,'   ND=',I4//)
290    FORMAT('0VEGAS  NDIM=',I3,'  NCALL=',F9.0,'  ITMX =',I5
     + ,'  ACC=',D9.3,'  MDS=',I3,'   ND=',I4)
201   FORMAT(/I1,'INTEGRATION BY VEGAS'/'0ITERATION NO',I3,
     +'.   INTEGRAL =',D14.8/20X,'STD DEV  =',D10.4/
     +' ACCUMULATED RESULTS.   INTEGRAL =',D14.8/
     +24X,'STD DEV  =',D10.4 / 24X,'CHI**2 PER ITN   =',D10.4)
202   FORMAT('0DATA FOR AXIS',I2 / 7X,'X',7X,'  DELT I  ',
     +2X,' CONVCE    ',11X,'X',7X,'  DELT I  ',2X,' CONVCE     '
     +,11X,'X',7X,'  DELT I  ',2X,' CONVCE     '/)
204   FORMAT(1X,3D12.4,5X,3D12.4,5X,3D12.4)
203   FORMAT(1X,I3,D20.8,D12.4,D20.8,D12.4,D12.4)
C
      END
C
C --- I02GAU
C --- SET PARAMETERS FOR HISTOGRAMS
      SUBROUTINE VHIS(XMIN,XMAX,NLPX,LTOPX,LLX,
     &                NLSX,NDDX,NAVEX,IHIS,TEXTX)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*32 TEXTX
      CHARACTER*32 TEXT(10),VTEXT(2)
      CHARACTER*1  CHR(40)
      CHARACTER*1  HMIN,HPLUS,HBLANK,HSTAR
      COMMON/VGPLIO/NIN,NOUT
      COMMON/VGLPLT/XL(10),V1(2),V2(2),AV(10)
      DOUBLE PRECISION XL,V1,V2,AV
      COMMON/VGPLTB/XKS(40,10),YKS(40,10)
      COMMON/VGRES/S1,S2,S3,S4
      DOUBLE PRECISION S1,S2,S3,S4
      DIMENSION ZAV(10),YAV(10),ZSV(10),YSV(10),ZTV(10)
      DIMENSION XLMAX(10),XLMIN(10),NLP(10),LTOP(10),LL(10)
      DIMENSION NUMB(12)
      DIMENSION XLS(42,10),YLS(42,10),NLSN(42,10),MLSN(42,10),DLS(10)
     +,XLAV(10),XLSQ(10),XLAVA(10),SXA(10),TLIM(6),TOP(10),XLTQ(10)
      DIMENSION NBIN(41),NLOG(41),SLOG(41),TLOG(41),HV(12)
      DIMENSION V1MAX(2),V1MIN(2),V2MAX(2),V2MIN(2),NV1(2)
     +,NV2(2)
      DIMENSION VM(12,12,2),NVM(12,12,2),BIN1(2),BIN2(2),VOL(2)
     +,WM(12,12,2),MVM(12,12,2)
C
C --- I02GAU
C
      COMMON /MIS1/ TEXT,VTEXT,
     &              CHR,
     &              HMIN,HPLUS,HBLANK,HSTAR
      COMMON /MIS2/ ZAV,YAV,ZSV,YSV,ZTV,
     &              XLMAX,XLMIN,NLP,LTOP,LL,
     &              NUMB
      COMMON /MIS3/ XLS,YLS,NLSN,MLSN
      COMMON /MIS4/ DLS,
     &              XLAV,XLSQ,XLAVA,SXA,TLIM,TOP,XLTQ,
     &              NBIN,NLOG
      COMMON /MIS5/ SLOG,TLOG,HV,
     &              V1MAX,V1MIN,V2MAX,V2MIN,NV1,
     &              NV2
      COMMON /MIS6/ VM,BIN1,BIN2,VOL,
     &              WM,MVM,NVM
      COMMON /MIS7/ FSQA,WTOW,
     &              KK,ITT,NLS,NDD,NAVE,NLPS,KT,I1,I2
C
      NLS=NLSX
      NDD=NDDX
      NAVE=NAVEX
      IF (IHIS .GT. 0) THEN
         XLMIN(IHIS)=XMIN
         XLMAX(IHIS)=XMAX
         NLP(IHIS)=NLPX
         LTOP(IHIS)=LTOPX
         LL(IHIS)=LLX
         DO 2, I=1,32
            TEXT(IHIS)(I:I)=' '
 2       CONTINUE
         I=1
 1       CONTINUE
            IF (TEXTX(I:I) .EQ. '$') GOTO 3
            TEXT(IHIS)(I:I)=TEXTX(I:I)
            I=I+1
         GOTO 1
 3       CONTINUE
      ENDIF
C
      RETURN
      END
C
C --- I02GAU
C
C --- I02GAU
      SUBROUTINE INPLOT(NOW,FF,PDX,IENTRY)
C --- I02GAU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C --- I02GAU
      DOUBLE PRECISION FF,PDX
      CALL INPL1OT(NOW,FF,PDX,IENTRY)
      RETURN
      END
C
      SUBROUTINE INPL1OT(NOW,FF,PDX,IENTRY)
C
C  AUTHOR                       : J. VERMASEREN
C  TERMINAL WIDTH ADDITION      : J. VAN DER HORST
C
C --- I02GAU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C --- I02GAU
      DOUBLE PRECISION FF,PDX
      CHARACTER*32 TEXT(10),VTEXT(2)
      CHARACTER*1  CHR(40)
      CHARACTER*1  HMIN,HPLUS,HBLANK,HSTAR
      COMMON/VGPLIO/NIN,NOUT
      COMMON/VGLPLT/XL(10),V1(2),V2(2),AV(10)
      DOUBLE PRECISION XL,V1,V2,AV
      COMMON/VGPLTB/XKS(40,10),YKS(40,10)
      COMMON/VGRES/S1,S2,S3,S4
      DOUBLE PRECISION S1,S2,S3,S4
      DIMENSION ZAV(10),YAV(10),ZSV(10),YSV(10),ZTV(10)
      DIMENSION XLMAX(10),XLMIN(10),NLP(10),LTOP(10),LL(10)
      DIMENSION NUMB(12)
      DIMENSION XLS(42,10),YLS(42,10),NLSN(42,10),MLSN(42,10),DLS(10)
     +,XLAV(10),XLSQ(10),XLAVA(10),SXA(10),TLIM(6),TOP(10),XLTQ(10)
      DIMENSION NBIN(41),NLOG(41),SLOG(41),TLOG(41),HV(12)
      DIMENSION V1MAX(2),V1MIN(2),V2MAX(2),V2MIN(2),NV1(2)
     +,NV2(2)
      DIMENSION VM(12,12,2),NVM(12,12,2),BIN1(2),BIN2(2),VOL(2)
     +,WM(12,12,2),MVM(12,12,2)
C
C --- I02GAU
C
      COMMON /MIS1/ TEXT,VTEXT,
     &              CHR,
     &              HMIN,HPLUS,HBLANK,HSTAR
      COMMON /MIS2/ ZAV,YAV,ZSV,YSV,ZTV,
     &              XLMAX,XLMIN,NLP,LTOP,LL,
     &              NUMB
      COMMON /MIS3/ XLS,YLS,NLSN,MLSN
      COMMON /MIS4/ DLS,
     &              XLAV,XLSQ,XLAVA,SXA,TLIM,TOP,XLTQ,
     &              NBIN,NLOG
      COMMON /MIS5/ SLOG,TLOG,HV,
     &              V1MAX,V1MIN,V2MAX,V2MIN,NV1,
     &              NV2
      COMMON /MIS6/ VM,BIN1,BIN2,VOL,
     &              WM,MVM,NVM
      COMMON /MIS7/ FSQA,WTOW,
     &              KK,ITT,NLS,NDD,NAVE,NLPS,KT,I1,I2
C
C --- I02GAU
C
      DATA TLIM/1.6,2.5,4.0,6.666666667E+00,10.,16./
      DATA MLS,MAV,NDMAX/10,10,2/
      DATA NGRAPH/0/
      DATA HMIN /'-'/ , HPLUS /'+'/ , HBLANK /' '/ , HSTAR /'*'/
C
C --- I02GAU
      GOTO (2101,2102,2103,2104) IENTRY
C --- I02GAU
C
      CALL VGDAT
C
      IGRAPH=NOW
      NOW=0
      KK=0
      ITT=0
C
      IF(IGRAPH.GT.0)THEN
         IF(IGRAPH.NE.NGRAPH)THEN
C --- I02GAU
C            READ(NIN,*)NLS
C --- I02GAU
            WRITE(NOUT,814)NLS
         ENDIF
         IF(NLS.LT.0) NLS=0
         IF(NLS.NE.0)THEN
            IF(NLS.GT.MLS)THEN
               WRITE(NOUT,813)MAV,MLS,NDMAX
               STOP
            ENDIF
            IF(IGRAPH.NE.NGRAPH) WRITE(NOUT,815)
            DO 801 I=1,NLS
               IF(IGRAPH.NE.NGRAPH)THEN
C --- I02GAU
C                  READ(NIN,*)XLMIN(I),XLMAX(I),NLP(I),LTOP(I),
C     +                         LL(I),TEXT(I)
C --- I02GAU
                  WRITE(NOUT,816)I,XLMIN(I),XLMAX(I),NLP(I),LTOP(I),
     +                           LL(I),TEXT(I)
               ENDIF
               IF(NLP(I).LT.1) NLP(I)=1
               IF(NLP(I).GT.40) NLP(I)=40
               DLS(I)=(XLMAX(I)-XLMIN(I))/NLP(I)
               NLPS=NLP(I)+2
               DO 300 J=1,NLPS
                  YLS(J,I)=0
                  MLSN(J,I)=0
300            CONTINUE
801         CONTINUE
         ENDIF
         IF(IGRAPH.NE.NGRAPH)THEN
C --- I02GAU
C            READ(NIN,*)NDD
C --- I02GAU
            WRITE(NOUT,817)NDD
         ENDIF
         IF(NDD.LT.0) NDD=0
         IF(NDD.NE.0)THEN
            IF(NDD.GT.NDMAX)THEN
               WRITE(NOUT,813)MAV,MLS,NDMAX
               STOP
            ENDIF
            IF(IGRAPH.NE.NGRAPH) WRITE(NOUT,818)
            DO 803 I=1,NDD
               IF(IGRAPH.NE.NGRAPH)THEN
                  READ(NIN,*)V1MIN(I),V1MAX(I),NV1(I),V2MIN(I),
     +                         V2MAX(I),NV2(I),VTEXT(I)
                  WRITE(NOUT,819)I,V1MIN(I),V1MAX(I),NV1(I),V2MIN(I),
     +                         V2MAX(I),NV2(I),VTEXT(I)
               ENDIF
               IF(NV1(I).LT.1)NV1(I)=1
               IF(NV2(I).LT.1) NV2(I)=1
               IF(NV1(I).GT.10) NV1(I)=10
               IF(NV2(I).GT.10) NV2(I)=10
               BIN1(I)=(V1MAX(I)-V1MIN(I))/NV1(I)
               BIN2(I)=(V2MAX(I)-V2MIN(I))/NV2(I)
               VOL(I)=BIN1(I)*BIN2(I)
803         CONTINUE
            WTOW=0.
            DO 805 I=1,NDD
            DO 805 J=1,12
            DO 805 K=1,12
               WM(K,J,I)=0.
               MVM(K,J,I)=0
805         CONTINUE
         ENDIF
         IF(IGRAPH.NE.NGRAPH)THEN
C --- I02GAU
C            READ(NIN,*)NAVE
C --- I02GAU
            WRITE(NOUT,820)NAVE
         ENDIF
         IF(NAVE.LT.0) NAVE=0
         IF(NAVE.GT.MAV)THEN
            WRITE(NOUT,813)MAV,MLS,NDMAX
            STOP
         ENDIF
         DO 11 I=1,MAV
            YAV(I)=0.
            YSV(I)=0.
11       CONTINUE
         KT=0
      ELSE
         NAVE=0
         NLS=0
         NDD=0
      ENDIF
      NGRAPH=IGRAPH
      RETURN
C
C --- I02GAU
 2101 CONTINUE
C --- I02GAU
C
C      ENTRY REPLOT(NOW,FF,PDX)
C
      IF(NAVE.NE.0)THEN
         DO 62 I=1,NAVE
            ZAV(I)=0.
            ZTV(I)=0.
            ZSV(I)=0.
62       CONTINUE
      ENDIF
      FSQA=0.
      KT=KT+1
      IF(NLS.NE.0)THEN
         DO 302 I=1,NLS
            NLPS=NLP(I)+2
            XLAV(I)=0.
            XLTQ(I)=0.
            XLSQ(I)=0.
            DO 302 J=1,NLPS
               XLS(J,I)=0.
               NLSN(J,I)=0
302      CONTINUE
      ENDIF
      IF(NDD.NE.0)THEN
         DO 402 I=1,NDD
            N1=NV1(I)+2
            N2=NV2(I)+2
            DO 402 I1=1,N1
            DO 402 I2=1,N2
               VM(I1,I2,I)=0.
               NVM(I1,I2,I)=0
402      CONTINUE
      ENDIF
      RETURN
C
C --- I02GAU
 2102 CONTINUE
C --- I02GAU
C
C      ENTRY XPLOT(NOW,FF,PDX)
C
      FSQA=FSQA+FF*FF/PDX
      ITT=ITT+1
      IF(NLS.NE.0)THEN
         DO 304  I=1,NLS
            NLPS=(XL(I)-XLMIN(I))/DLS(I)+1.
            IF(NLPS.LT.0) NLPS=0
            IF(NLPS.GT.NLP(I))NLPS=NLP(I)+1
            NLPS=NLPS+1
            XLS(NLPS,I)=XLS(NLPS,I)+FF/DLS(I)
            NLSN(NLPS,I)=NLSN(NLPS,I)+1
            XLAV(I)=XLAV(I)+FF*XL(I)
            XLTQ(I)=XLTQ(I)+FF*FF*XL(I)/PDX
            XLSQ(I)=XLSQ(I)+(FF*XL(I))**2/PDX
304      CONTINUE
      ENDIF
      IF(NDD.NE.0)THEN
         DO 404 I=1,NDD
            I1=(V1(I)-V1MIN(I))/BIN1(I)+2
            IF(I1.LT.1) I1=1
            IF(I1.GT.NV1(I)+2) I1=NV1(I)+2
            I2=(V2(I)-V2MIN(I))/BIN2(I)+2
            IF(I2.LT.1) I2=1
            IF(I2.GT.NV2(I)+2) I2=NV2(I)+2
            VM(I1,I2,I)=VM(I1,I2,I)+FF/VOL(I)
            NVM(I1,I2,I)=NVM(I1,I2,I)+1
404      CONTINUE
      ENDIF
      IF(NAVE.EQ.0)RETURN
      DO 22 I=1,NAVE
         ZAV(I)=ZAV(I)+AV(I)*FF
         ZTV(I)=ZTV(I)+FF*FF*AV(I)/PDX
         ZSV(I)=ZSV(I)+(AV(I)*FF)**2/PDX
22    CONTINUE
99    RETURN
C
C --- I02GAU
 2103 CONTINUE
C --- I02GAU
C
C      ENTRY PLOTIT(NOW,FF,PDX)
C
      IF(NLS.EQ.0) GO TO 315
      IF(KK.LE.0)THEN
         DO 306 I=1,NLS
            NLPS=NLP(I)+2
            DO 306 J=1,NLPS
               MLSN(J,I)=NLSN(J,I)
               YLS(J,I)=XLS(J,I)
306      CONTINUE
      ELSE
         VBEF=VTOT
         VU=(S4/S3)**2
         DO 309 I=1,NLS
            NLPS=NLP(I)+2
            DO 309 J=1,NLPS
               IF(NLSN(J,I).EQ.0) GO TO 309
               IF(MLSN(J,I).NE.0)THEN
                  AL1=VU/NLSN(J,I)
                  AL2=VBEF/MLSN(J,I)
                  MLSN(J,I)=MLSN(J,I)+NLSN(J,I)
                  YLS(J,I)=(AL2*XLS(J,I)+AL1*YLS(J,I))/(AL1+AL2)
               ELSE
                  MLSN(J,I)=NLSN(J,I)
                  YLS(J,I)=XLS(J,I)
               ENDIF
309      CONTINUE
      ENDIF
      DO 311 I=1,NLS
         SXF=XLSQ(I)-XLAV(I)*XLAV(I)
         SXT=XLTQ(I)-XLAV(I)*S3
         SX2=XLSQ(I)/XLAV(I)**2+FSQA/S3**2-2.*XLTQ(I)/(XLAV(I)*S3)
         SX2=SX2*(XLAV(I)/S3)**2
         IF(KT.EQ.1)THEN
            XLAVA(I)=XLAV(I)/S3
            SXA(I)=SX2
         ELSE
            XHELP=SX2+SXA(I)
            IF(XHELP.NE.0)THEN
               XLAVA(I)=(XLAV(I)*SXA(I)/S3+XLAVA(I)*SX2)/XHELP
               SXA(I)=SXA(I)*SX2/XHELP
            ENDIF
         ENDIF
311   CONTINUE
      VTOT=(S2/S1)**2
      IF(NOW.EQ.2)THEN
         DO 341 I=1,NLS
            TOP(I)=0.
            NLPS=NLP(I)+1
            DO 341 J=2,NLPS
               XLS(J,I)=YLS(J,I)/S1
               IF(XLS(J,I).GT.TOP(I))TOP(I)=XLS(J,I)
341         CONTINUE
         DO 342 I=1,NLS
            IF(LTOP(I).LE.0) LTOP(I)=I
            LTO=LTOP(I)
            IF(TOP(I).GT.TOP(LTO))TOP(LTO)=TOP(I)
342      CONTINUE
         YLOG=DLOG10(S1*S1)*.5
         DO 314 I=1,NLS
            WRITE(NOUT,321)I
            NLPS=NLP(I)+1
            LTO=LTOP(I)
            TOP(I)=TOP(LTO)
            IF(TOP(I).EQ.0) TOP(I)=1.
            AN1=DLOG10(TOP(I))
            N1=AN1
            IF(N1.GT.AN1) N1=N1-1
            Z1=TOP(I)*10.**(-N1)
            DO 343 L=1,4
               IF(Z1.LT.TLIM(L)) GO TO 344
343         CONTINUE
            L=5
344         IF(TOP(I).LT.1.6/(XLMAX(I)-XLMIN(I)))L=L+1
            TOPM=TLIM(L)*10.**N1
            DO 345 J=2,NLPS
               NBIN(J)=( XLS(J,I)*20./TOPM+1.5 )
               IF(LL(I).LT.0)NBIN(J)=0
               IF(XLS(J,I).GT.0)THEN
                  TLOG(J)=DLOG10(XLS(J,I))
                  SLOG(J)=TLOG(J)+YLOG
                  NLOG(J)=( (TLOG(J)-N1)*8.+33.5 )  * .5
                  IF(LL(I).GT.0) NLOG(J)=0
               ELSE
                  SLOG(J)=0
                  TLOG(J)=0
                  NLOG(J)=0
               ENDIF
345         CONTINUE
            WRITE(NOUT,322)TEXT(I)
            N1P1=N1+1
            N1M4=N1 -4
            WRITE(NOUT,323)
            XMIN=XLMIN(I)
            XMAX=XMIN+DLS(I)
            WRITE(NOUT,324)XMIN,XMAX,YLS(2,I),SLOG(2),XLS(2,I),
     +                     TLOG(2),MLSN(2,I)
            XKS(1,I)  = (XMIN+XMAX)/2.
            YKS(1,I)  = YLS(2,I)
            DO 348 J=3,NLPS
               XMIN=XMAX
               XMAX=XMIN+DLS(I)
               XKS(J-1,I) = (XMIN+XMAX)/2.
               YKS(J-1,I) = YLS(J,I)
               WRITE(NOUT,324)XMIN,XMAX,YLS(J,I),SLOG(J),XLS(J,I),
     +                        TLOG(J),MLSN(J,I)
348         CONTINUE
            WRITE(NOUT,325)
            EL1=YLS(1,I)*DLS(I)
            EL2=EL1/S1
            WRITE(NOUT,327)EL1,EL2,MLSN(1,I)
            EL1=YLS(42,I)*DLS(I)
            EL2=EL1/S1
C --- I02GAU
CC          WRITE(NOUT,328)EL1,EL2,MLSN(NLPS1,I)
            WRITE(NOUT,328)EL1,EL2,MLSN(NLPS+1,I)
C --- I02GAU
            SXSQ=SQRT(SXA(I)/ITT)
            WRITE(NOUT,329)XLAVA(I),SXSQ
            WRITE(NOUT,723)TLIM(L),N1,N1P1,N1M4
            DO 991 L=1,40
               CHR(L)=HMIN
               IF(NLOG(L+1).EQ.21) CHR(L)=HPLUS
               IF(NBIN(L+1).EQ.21) CHR(L)=HSTAR
991         CONTINUE
            WRITE(NOUT,724)CHR
            NLPSV  =  NLPS
            IF(NLPSV .GT. 21)NLPSV=21
            DO 993 J=3,NLPSV
               DO 992 L=1,40
                  CHR(L)=HBLANK
                  IF(NLOG(L+1).EQ.23-J) CHR(L)=HPLUS
                  IF(NBIN(L+1).EQ.23-J) CHR(L)=HSTAR
992            CONTINUE
               WRITE(NOUT,724)CHR
993         CONTINUE
            NLPS1=NLPSV+1
            IF(NLPSV.EQ.21) GO TO 996
            DO 995 J=NLPS1,21
               DO 994 L=1,40
                  CHR(L)=HBLANK
                  IF(NLOG(L+1).EQ.23-J) CHR(L)=HPLUS
                  IF(NBIN(L+1).EQ.23-J) CHR(L)=HSTAR
994            CONTINUE
               WRITE(NOUT,724)CHR
995         CONTINUE
996         DO 997 L=1,40
               CHR(L)=HMIN
               IF(NLOG(L+1).EQ.1) CHR(L)=HPLUS
               IF(NBIN(L+1).EQ.1) CHR(L)=HSTAR
997         CONTINUE
            WRITE(NOUT,724)CHR
314      CONTINUE
      ENDIF
315   IF(NDD.NE.0)THEN
         WBEF=WTOT
         DO 500 I=1,NDD
            NX=NV1(I)+2
            NY=NV2(I)+2
            IF(KK.LE.0)THEN
               DO 501 J=1,NX
               DO 501 K=1,NY
                  WM(J,K,I)=VM(J,K,I)
                  MVM(J,K,I)=NVM(J,K,I)
501            CONTINUE
            ELSE
               VU=(S4/S3)**2
               DO 503 J=1,NX
               DO 503 K=1,NY
                  IF(NVM(J,K,I).NE.0)THEN
                     IF(MVM(J,K,I).NE.0)THEN
                      AL1=VU/NVM(J,K,I)
                      AL2=VBEF/MVM(J,K,I)
                      MVM(J,K,I)=MVM(J,K,I)+NVM(J,K,I)
                      WM(J,K,I)=(AL2*VM(J,K,I)+AL1*WM(J,K,I))/(AL1+AL2)
                     ELSE
                      MVM(J,K,I)=NVM(J,K,I)
                      WM(J,K,I)=VM(J,K,I)
                     ENDIF
                  ENDIF
503            CONTINUE
            ENDIF
500      CONTINUE
         WTOT=(S2/S1)**2
         IF(NOW.EQ.2)THEN
            DO 408 I=1,NDD
               WRITE(NOUT,481)I,VTEXT(I)
               VVV=V2MAX(I)
               MVV=NV1(I)+2
               NVV=NV2(I)+1
               SIZE=VOL(I)/S1
               DO 406 I2=1,NVV
                  J2=NVV+2-I2
                  DO 410 I1=1,MVV
                     NUMB(I1)=1000.*WM(I1,J2,I)*SIZE+.5
410               CONTINUE
                  WRITE(NOUT,486)(NUMB(I1),I1=1,MVV)
                  WRITE(NOUT,483)(WM(I1,J2,I),I1=1,MVV)
                  WRITE(NOUT,486)(MVM(I1,J2,I),I1=1,MVV)
                  WRITE(NOUT,484)VVV
                  VVV=VVV-BIN2(I)
                  IF(ABS(VVV/BIN2(I)).LT.1.E-10)VVV=0.
406            CONTINUE
               DO 411 I1=1,MVV
                  NUMB(I1)=1000.*WM(I1,1,I)*SIZE+.5
411            CONTINUE
               WRITE(NOUT,486)(NUMB(I1),I1=1,MVV)
               WRITE(NOUT,483)(WM(I1,1,I),I1=1,MVV)
               WRITE(NOUT,486)(MVM(I1,1,I),I1=1,MVV)
               WRITE(NOUT,482)
               MVV=MVV-1
               DO 407 I1=1,MVV
                  HV(I1)=V1MIN(I)+(I1-1)*BIN1(I)
                  IF(ABS(HV(I1)/BIN1(I)).LT.1.E-10)HV(I1)=0.
407            CONTINUE
               WRITE(NOUT,485)(HV(I1),I1=1,MVV)
408         CONTINUE
         ENDIF
      ENDIF
      IF(NAVE.NE.0)THEN
         IF(NOW.EQ.2) WRITE(NOUT,26)
         DO 24 I=1,NAVE
            SXF=ZSV(I)-ZAV(I)*ZAV(I)
            SXT=ZSV(I)/ZAV(I)**2+FSQA/S3**2-2.*ZTV(I)/(ZAV(I)*S3)
            SX2=SXT*(ZAV(I)/S3)**2
            IF(KT.EQ.1)THEN
               YAV(I)=ZAV(I)/S3
               YSV(I)=SX2
            ELSE
               XHELP=SX2+YSV(I)
               IF(XHELP.NE.0)THEN
                  YAV(I)=(YSV(I)*ZAV(I)/S3+YAV(I)*SX2)/XHELP
                  YSV(I)=YSV(I)*SX2/XHELP
               ENDIF
            ENDIF
            YSSQ=SQRT(YSV(I)/ITT)
            IF(NOW.EQ.2) WRITE(NOUT,27)I,YAV(I),YSSQ
24       CONTINUE
      ENDIF
      NOW=1
      KK=KK+1
      RETURN
C
27    FORMAT(12X,I2,9X,E15.5,5X,E15.3)
26    FORMAT('1',10X,'THE FOLLOWING ARE AVERAGES WITH ERROR ESTIMATE'/)
321   FORMAT('1',10X,'SINGLE DIFFERENTIAL CROSS-SECTION NUMBER',I3///)
322   FORMAT(' SINGLE DIFFERENTIAL CROSS SECTION OF ',A32/)
323   FORMAT(11X,'LIMITS',9X,1HI,16X,'ACCUMULATED RESULTS '
     +      /26X,'I'/5X,'LOWER',7X,'UPPER',4X,1HI,5X,'DS/DX',
     +       4X,'DLOG10   (DS/DX)/S  DLOG10  POINTS'
     +      /2X,24('-'),'I',50('-'))
324   FORMAT(E12.4,E12.4,3H  I,2(E12.4,F8.2),I8)
325   FORMAT(2X,24('-'),'I',50('-'))
723   FORMAT(//41X,'UPPER BIN',6X,'LOWER BIN'/17X,'* LINEAR      PLOT',
     +       F8.2,'*10**',I3,8X,'0'/17X,'+ LOGARITHMIC PLOT',6X,
     +       '10**',I3,8X,'10**',I3)
724   FORMAT(19X,'I',40A1,'I')
327   FORMAT(7X,'TOTAL UNDERFLOW',4X,1HI,E12.4,E20.4,I16)
328   FORMAT(7X,'TOTAL  OVERFLOW',4X,1HI,E12.4,E20.4,I16)
329   FORMAT(//19X,'ACCUMULATED AVERAGE =',E12.5
     +/19X,'ESTIMATED ERROR     =',E12.5)
481   FORMAT('1',45X,'DOUBLE DIFFERENTIAL CROSS-SECTION NUMBER',I3/
     +/60X,'X-AXIS ',3A4/60X,'Y-AXIS ',A32/)
482   FORMAT(20X,11('I',9X))
483   FORMAT(11X,E9.3,11('I',E9.3))
484   FORMAT(1X,E10.3,  '---------',11('I---------'))
485   FORMAT('0',14X,11E10.3)
486   FORMAT(11X,I8,1X,11('I',I8,1X))
810   FORMAT(I2)
811   FORMAT(2E12.4,3I2,8A4)
812   FORMAT(2(2E10.3,I4),6A4)
813   FORMAT('1***ERROR***',10X,'TOO MANY PLOTS REQUESTED'//
     +22X,'THE UPPER LIMITS ARE '//19X,I2,' AVERAGES'//19X,I2,
     +' ONE DIMENSIONAL PLOTS'//19X,I2,' TWO DIMENSIONAL PLOTS'////
     +22X,'***EXECUTION IS HALTED***')
814   FORMAT('1NUMBER OF SINGLE DIFFERENTIAL CROSS SECTIONS ',
     +'REQUESTED =',I3/)
C --- I02GAU
815   FORMAT(' INFORMATION ON THE DATA CARDS'//
     +'  I',10X,'XLMIN',12X,'XLMAX',7X,
     +'BINS  CORRELATION  TYPE',5X,'TEXT'/)
C815   FORMAT(' INFORMATION ON THE DATA CARDS'//
C     +'  I',10X,'XLMIN',12X,'XLMAX',7X,
C     +'BINS  CORRELLATION  TYPE',5X,'TEXT'/)
C --- I02GAU
816   FORMAT(I3,2E17.4,I8,I9,I10,5X,A32)
817   FORMAT('0NUMBER OF DOUBLE DIFFERENTIAL CROSS SECTIONS',
     +' REQUESTED =',I3/)
818   FORMAT(' INFORMATION ON THE DATA CARDS'//
     +'  I',7X,'V1MIN',7X,'V1MAX',2X,'BINS',7X,'V2MIN'
     +,7X,'V2MAX',2X,'BINS',5X,'TEXT 1',8X,'TEXT 2'/)
819   FORMAT(I3,2E12.3,I5,1X,2E12.3,I5,6X,3A4,2X,A32)
820   FORMAT('0NUMBER OF AVERAGES REQUESTED =',I3)
C
C
C --- I02GAU
 2104 CONTINUE
C --- I02GAU
C
C      ENTRY PLOTTA(NOW,FF,PDX)
C
      IF(NLS.EQ.0) GO TO 1315
      IF(KK.LE.0)THEN
         DO 1306 I=1,NLS
            NLPS=NLP(I)+2
            DO 1306 J=1,NLPS
               MLSN(J,I)=NLSN(J,I)
               YLS(J,I)=XLS(J,I)
1306     CONTINUE
      ELSE
         VBEF=VTOT
         VU=(S4/S3)**2
         DO 1309 I=1,NLS
            NLPS=NLP(I)+2
            DO 1309 J=1,NLPS
               IF(NLSN(J,I).EQ.0) GO TO 1309
C              IF(MLSN(J,I).EQ.0)THEN
CHH
               IF(MLSN(J,I).NE.0)THEN
CHH
                  AL1=VU/NLSN(J,I)
                  AL2=VBEF/MLSN(J,I)
                  MLSN(J,I)=MLSN(J,I)+NLSN(J,I)
                  YLS(J,I)=(AL2*XLS(J,I)+AL1*YLS(J,I))/(AL1+AL2)
               ELSE
                  MLSN(J,I)=NLSN(J,I)
                  YLS(J,I)=XLS(J,I)
               ENDIF
1309     CONTINUE
      ENDIF
      DO 1311 I=1,NLS
         SXF=XLSQ(I)-XLAV(I)*XLAV(I)
         SXT=XLTQ(I)-XLAV(I)*S3
         SX2=XLSQ(I)/XLAV(I)**2+FSQA/S3**2-2.*XLTQ(I)/(XLAV(I)*S3)
         SX2=SX2*(XLAV(I)/S3)**2
         IF(KT.EQ.1)THEN
            XLAVA(I)=XLAV(I)/S3
            SXA(I)=SX2
         ELSE
            XHELP=SX2+SXA(I)
            IF(XHELP.NE.0)THEN
               XLAVA(I)=(XLAV(I)*SXA(I)/S3+XLAVA(I)*SX2)/XHELP
               SXA(I)=SXA(I)*SX2/XHELP
            ENDIF
         ENDIF
1311  CONTINUE
      VTOT=(S2/S1)**2
      IF(NOW.EQ.2)THEN
         DO 1341 I=1,NLS
            TOP(I)=0.
            NLPS=NLP(I)+1
            DO 1341 J=2,NLPS
               XLS(J,I)=YLS(J,I)/S1
               IF(XLS(J,I).GT.TOP(I))TOP(I)=XLS(J,I)
1341     CONTINUE
         DO 1342 I=1,NLS
            IF(LTOP(I).LE.0) LTOP(I)=I
            LTO=LTOP(I)
            IF(TOP(I).GT.TOP(LTO))TOP(LTO)=TOP(I)
1342     CONTINUE
         YLOG=DLOG10(S1*S1)*.5
         DO 1314 I=1,NLS
            WRITE(NOUT,321)I
            NLPS=NLP(I)+1
            LTO=LTOP(I)
            TOP(I)=TOP(LTO)
            IF(TOP(I).EQ.0) TOP(I)=1.
            AN1=DLOG10(TOP(I))
            N1=AN1
            IF(N1.GT.AN1) N1=N1-1
            Z1=TOP(I)*10.**(-N1)
            DO 1343 L=1,4
               IF(Z1.LT.TLIM(L)) GO TO 1344
1343        CONTINUE
            L=5
1344        IF(TOP(I).LT.1.6/(XLMAX(I)-XLMIN(I)))L=L+1
            TOPM=TLIM(L)*10.**N1
            DO 1345 J=2,NLPS
               NBIN(J)=XLS(J,I)*40./TOPM+1.5
               IF(LL(I).LT.0)NBIN(J)=0
               IF(XLS(J,I).GT.0)THEN
                  TLOG(J)=DLOG10(XLS(J,I))
                  SLOG(J)=TLOG(J)+YLOG
                  NLOG(J)=(TLOG(J)-N1)*8.+33.5
                  IF(LL(I).GT.0) NLOG(J)=0
               ELSE
                  SLOG(J)=0
                  TLOG(J)=0
                  NLOG(J)=0
               ENDIF
1345        CONTINUE
            WRITE(NOUT,322)TEXT(I)
            N1P1=N1+1
            N1M4=N1 -4
            WRITE(NOUT,1323)TLIM(L),N1,N1P1,N1M4
            DO 1347 L=1,40
               CHR(L)=HMIN
               IF(NLOG(L+1).EQ.41) CHR(L)=HPLUS
               IF(NBIN(L+1).EQ.41) CHR(L)=HSTAR
1347        CONTINUE
            XMIN=XLMIN(I)
            XMAX=XMIN+DLS(I)
            WRITE(NOUT,1324)XMIN,XMAX,YLS(2,I),SLOG(2),XLS(2,I),
     +                      TLOG(2),MLSN(2,I),CHR
            XKS(1,I)  = (XMIN+XMAX)/2.
            YKS(1,I)  = YLS(2,I)
            DO 1348 J=3,NLPS
               XMIN=XMAX
               XMAX=XMIN+DLS(I)
               XKS(J-1,I) = (XMIN+XMAX)/2.
               YKS(J-1,I) = YLS(J,I)
               DO 1349 L=1,40
                  CHR(L)=HBLANK
                  IF(NLOG(L+1).EQ.43-J) CHR(L)=HPLUS
                  IF(NBIN(L+1).EQ.43-J) CHR(L)=HSTAR
1349           CONTINUE
               WRITE(NOUT,1324)XMIN,XMAX,YLS(J,I),SLOG(J),XLS(J,I),
     +                         TLOG(J),MLSN(J,I),CHR
1348        CONTINUE
            NLPS1=NLPS+1
            IF(NLPS.EQ.41) GO TO 1352
            DO 1351 J=NLPS1,41
               DO 1350 L=1,40
                  CHR(L)=HBLANK
                  IF(NLOG(L+1).EQ.43-J) CHR(L)=HPLUS
                  IF(NBIN(L+1).EQ.43-J) CHR(L)=HSTAR
1350           CONTINUE
1351        CONTINUE
            WRITE(NOUT,1325)CHR
1352        DO 1353 L=1,40
               CHR(L)=HMIN
               IF(NLOG(L+1).EQ.1) CHR(L)=HPLUS
               IF(NBIN(L+1).EQ.1) CHR(L)=HSTAR
1353        CONTINUE
            WRITE(NOUT,1326)CHR
            EL1=YLS(1,I)*DLS(I)
            EL2=EL1/S1
            WRITE(NOUT,1327)EL1,EL2,MLSN(1,I)
            EL1=YLS(42,I)*DLS(I)
            EL2=EL1/S1
            WRITE(NOUT,1328)EL1,EL2,MLSN(NLPS1,I)
            SXSQ=SQRT(SXA(I)/ITT)
            WRITE(NOUT,1329)XLAVA(I),SXSQ
1314     CONTINUE
      ENDIF
1315  CONTINUE
      IF(NDD.NE.0)THEN
         WBEF=WTOT
         DO 1500 I=1,NDD
            NX=NV1(I)+2
            NY=NV2(I)+2
            IF(KK.LE.0)THEN
               DO 1501 J=1,NX
                DO 1501 K=1,NY
                   WM(J,K,I)=VM(J,K,I)
                   MVM(J,K,I)=NVM(J,K,I)
1501           CONTINUE
            ELSE
               VU=(S4/S3)**2
               DO 1503 J=1,NX
                DO 1503 K=1,NY
                   IF(NVM(J,K,I).NE.0)THEN
                      IF(MVM(J,K,I).NE.0)THEN
                         AL1=VU/NVM(J,K,I)
                         AL2=VBEF/MVM(J,K,I)
                         MVM(J,K,I)=MVM(J,K,I)+NVM(J,K,I)
                         WM(J,K,I)=(AL2*VM(J,K,I)+AL1*WM(J,K,I))/
     +                             (AL1+AL2)
                      ELSE
                         MVM(J,K,I)=NVM(J,K,I)
                         WM(J,K,I)=VM(J,K,I)
                      ENDIF
                   ENDIF
1503           CONTINUE
            ENDIF
1500     CONTINUE
         WTOT=(S2/S1)**2
         IF(NOW.EQ.2)THEN
            DO 1408 I=1,NDD
               WRITE(NOUT,481)I,VTEXT(I)
               VVV=V2MAX(I)
               MVV=NV1(I)+2
               NVV=NV2(I)+1
               SIZE=VOL(I)/S1
               DO 1406 I2=1,NVV
                  J2=NVV+2-I2
                  DO 1410 I1=1,MVV
                     NUMB(I1)=1000.*WM(I1,J2,I)*SIZE+.5
1410              CONTINUE
                  WRITE(NOUT,486)(NUMB(I1),I1=1,MVV)
                  WRITE(NOUT,483)(WM(I1,J2,I),I1=1,MVV)
                  WRITE(NOUT,486)(MVM(I1,J2,I),I1=1,MVV)
                  WRITE(NOUT,484)VVV
                  VVV=VVV-BIN2(I)
                  IF(ABS(VVV/BIN2(I)).LT.1.E-10)VVV=0.
1406           CONTINUE
               DO 1411 I1=1,MVV
                  NUMB(I1)=1000.*WM(I1,1,I)*SIZE+.5
1411           CONTINUE
               WRITE(NOUT,486)(NUMB(I1),I1=1,MVV)
               WRITE(NOUT,483)(WM(I1,1,I),I1=1,MVV)
               WRITE(NOUT,486)(MVM(I1,1,I),I1=1,MVV)
               WRITE(NOUT,482)
               MVV=MVV-1
               DO 1407 I1=1,MVV
                  HV(I1)=V1MIN(I)+(I1-1)*BIN1(I)
                  IF(ABS(HV(I1)/BIN1(I)).LT.1.E-10)HV(I1)=0.
1407           CONTINUE
               WRITE(NOUT,485)(HV(I1),I1=1,MVV)
1408        CONTINUE
         ENDIF
      ENDIF
      IF(NAVE.NE.0)THEN
         IF(NOW.EQ.2) WRITE(NOUT,26)
         DO 124 I=1,NAVE
            SXF=ZSV(I)-ZAV(I)*ZAV(I)
            SXT=ZSV(I)/ZAV(I)**2+FSQA/S3**2-2.*ZTV(I)/(ZAV(I)*S3)
            SX2=SXT*(ZAV(I)/S3)**2
            IF(KT.EQ.1)THEN
               YAV(I)=ZAV(I)/S3
               YSV(I)=SX2
            ELSE
               XHELP=SX2+YSV(I)
               IF(XHELP.NE.0)THEN
                  YAV(I)=(YSV(I)*ZAV(I)/S3+YAV(I)*SX2)/XHELP
                  YSV(I)=YSV(I)*SX2/XHELP
               ENDIF
            ENDIF
            YSSQ=SQRT(YSV(I)/ITT)
            IF(NOW.EQ.2) WRITE(NOUT,27)I,YAV(I),YSSQ
124      CONTINUE
      ENDIF
      NOW=1
      KK=KK+1
      RETURN
C
1323  FORMAT(11X,'LIMITS',9X,'I',16X,'ACCUMULATED RESULTS',15X,'I',24X
     +,'UPPER BIN',6X,'LOWER BIN'/26X,'I',50X,'I * LINEAR      PLOT',
     +F8.2,'*10**',I3,8X,'0'/5X,'LOWER',7X,'UPPER',4X,'I',5X,'DS/DX',4X
     +,'DLOG10   (DS/DX)/S  DLOG10  POINTS  I + LOGARITHMIC PLOT',6X
     +,'10**',I3,8X,'10**',I3/2X,24('-'),'I',50('-'),'I')
1324  FORMAT(E12.4,E12.4,'  I',2(E12.4,F8.2),I8,'  I',4X,'I',40A1,'I')
1325  FORMAT(26X,'I',50X,'I',4X,'I',40A1,'I')
1326  FORMAT(2X,24('-'),'I',50('-'),'I',4X,'I',40A1,'I')
1327  FORMAT(7X,'TOTAL UNDERFLOW',4X,'I',E12.4,E20.4,I16,2X,'I')
1328  FORMAT(7X,'TOTAL  OVERFLOW',4X,'I',E12.4,E20.4,I16,2X,'I')
1329  FORMAT(//19X,'ACCUMULATED AVERAGE =',E12.5
     +/19X,'ESTIMATED ERROR     =',E12.5)
C
      END
       SUBROUTINE PICKIN(S,V1,V2,V3,V4,V5,DJ,NOPT,Y)
C
C  AUTHOR          : J. VERMASEREN
C  DEBUGGED        : G.-J. VAN OLDENBORGH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( PI = 3.14159265 )
       DIMENSION Y(4)
      COMMON/VGPZZ/W1,W2,W3,W4,W5,D1,D2,D5,D7,SL1
      COMMON/VGEXTR/TS1,TS2,T1,T2
      COMMON/VGLEVI/GRAM,DD1,DD2,DD3,DD4,DD5,DELTA,G4,SA1,SA2
      COMMON/VGDOTP/P12,P13,P14,P15,P23,P24,P25,P34,P35,P45,P1K2,P2K1
C
C
       X1=Y(1)
       X2=Y(2)
       X3=Y(3)
       W1=V1*V1
       W2=V2*V2
       W3=V3*V3
       W4=V4*V4
       W5=V5*V5
       SIG=V4+V5
       SIG1=SIG*SIG
       SIG2=SIG1
       D1=W3-W1
       D2=W5-W2
       D5=W1-W2
       D6=W4-W5
       SS=S+D5
       RL1=SS*SS-4.*W1*S
       IF(RL1.LE.0)THEN
          DJ=0.0
          RETURN
       ENDIF
       SL1=DSQRT(RL1)
       IF(NOPT.EQ.0)THEN
          SMAX=S+W3-2.*V3*DSQRT(S)
          CALL MAPS2(TS2,X3,SIG1,SMAX,DS2)
          SIG1=TS2
          SP=S+W3-SIG1
       ENDIF
       D3=SIG1-W2
       RL2=SP*SP-4.*S*W3
       IF(RL2.LE.0)THEN
          DJ=0.0
          RETURN
       ENDIF
       SL2=DSQRT(RL2)
       T1MAX=W1+W3-(SS*SP+SL1*SL2)/(2.*S)
       T1MIN=(D1*D3+(D3-D1)*(D3*W1-D1*W2)/S)/T1MAX
C********************************************
C  DEBUG ADDED BY GEERT-JAN VAN OLDENBORGH
C
       IF(NOPT.NE.0)THEN
          S2TMIN=(S*(V1-V3)-V3*(D5-V1*V3))/V1
          IF(S2TMIN.GT.S2MIN)T1MIN=(V1-V3)**2
       ENDIF
C********************************************
       CALL MAPT1(T1,X1,T1MIN,T1MAX,DT1)
       D4=W4-T1
       D8=T1-W2
       T13=T1-W1-W3
       SA1=-(T1-D1)*(T1-D1)/4.+W1*T1
       IF(SA1.GE.0)THEN
          DJ=0.0
          RETURN
       ENDIF
       SL3=DSQRT(-SA1)
       IF(W1.NE.0)THEN
          SB=(S*(T1-D1)+D5*T13)/(2.*W1)+W3
          SD=SL1*SL3/W1
          SE=(S*(T1*(S+T13-W2)-W2*D1)+W3*(D5*D8+W2*W3))/W1
          IF(DABS((SB-SD)/SD).LT.1.0)THEN
             S2MAX=SB+SD
             SPLUS=SE/S2MAX
          ELSE
             SPLUS=SB-SD
             S2MAX=SE/SPLUS
          ENDIF
       ELSE
          S2MAX=(S*(T1*(S+D8-W3)-W2*W3)+W2*W3*(W2+W3-T1))/SS/T13
          SPLUS=SIG2
       ENDIF
       S2X=S2MAX
       IF(NOPT)5,6,7
5      IF(SPLUS.GT.SIG2)SIG2=SPLUS
       IF(NOPT.LT.-1)CALL MAPS2(TS2,X3,SIG2,S2MAX,DS2)
       IF(NOPT.EQ.-1)CALL MAPLA(TS2,T1,W2,X3,SIG2,S2MAX,DS2)
6      S2X=TS2
7      R1=S2X-D8
       R2=S2X-D6
       RL4=(R1*R1-4.*W2*S2X)*(R2*R2-4.*W5*S2X)
       IF(RL4.LE.0)THEN
          DJ=0.0
          RETURN
       ENDIF
       SL4=DSQRT(RL4)
       T2MAX=W2+W5-(R1*R2+SL4)/(2.*S2X)
       T2MIN=(D2*D4+(D4-D2)*(D4*W2-D2*T1)/S2X)/T2MAX
       CALL MAPT2(T2,X2,T2MIN,T2MAX,DT2)
       D7=T1-T2
       R3=D4-T2
       R4=D2-T2
       B=R3*R4-2.*(T1+W2)*T2
       C=T2*D6*D8+(D6-D8)*(D6*W2-D8*W5)
       T25=T2-W2-W5
       SA2=-R4*R4/4.+W2*T2
       IF(SA2.GE.0)THEN
          DJ=0.0
          RETURN
       ENDIF
       SL6=2.*DSQRT(-SA2)
       G4=-R3*R3/4.+T1*T2
       IF(G4.GE.0)THEN
          DJ=0.0
          RETURN
       ENDIF
       SL7=DSQRT(-G4)*2.
       SL5=SL6*SL7
       IF(DABS((SL5-B)/SL5).LT.1.0)THEN
          S2MIN=(-SL5-B)/(2.*T2)
          S2P=C/(T2*S2MIN)
       ELSE
          S2P=(SL5-B)/(2.*T2)
          S2MIN=C/(T2*S2P)
       ENDIF
       IF(NOPT.GT.1)CALL MAPS2(TS2,X3,S2MIN,S2MAX,DS2)
       IF(NOPT.EQ.1)CALL MAPLA(TS2,T1,W2,X3,S2MIN,S2MAX,DS2)
       AP=-(TS2+D8)*(TS2+D8)/4.+TS2*T1
       IF(W1.NE.0)THEN
          DD1=-W1*(TS2-S2MAX)*(TS2-SPLUS)/4.
       ELSE
          DD1=SS*T13*(TS2-S2MAX)/4.
       ENDIF
       DD2=-T2*(TS2-S2P)*(TS2-S2MIN)/4.
       YY4=COS(PI*Y(4))
       DD=DD1*DD2
       P12=0.5*(S-W1-W2)
       ST=TS2-T1-W2
       DELB=(2.*W2*R3+R4*ST)*(4.*P12*T1-(T1-D1)*ST)/(16.*AP)
       IF(DD.LE.0)THEN
          DJ=0.0
          RETURN
       ENDIF
       DELTA=DELB-YY4*ST*DSQRT(DD)/(2.*AP)
       TS1=(2.*P12*R3-4.*DELTA)/ST+T2+W1
       IF(AP.GE.0)THEN
          DJ=0.0
          RETURN
       ENDIF
       DJ=DS2*DT1*DT2*PI*PI/(8.*SL1*DSQRT(-AP))
       GRAM=(1.-YY4)*(1.+YY4)*DD/AP
       P13=-T13/2.
       P14=(D7+TS1-W3)/2.
       P15=(S+T2-TS1-W2)/2.
       P23=(S+T1-TS2-W1)/2.
       P24=(TS2-D7-W5)/2.
       P25=-T25/2.
       P34=(TS1-W3-W4)/2.
       P35=(S+W4-TS1-TS2)/2.
       P45=(TS2-W4-W5)/2.
       P1K2=(TS1-T2-W1)/2.
       P2K1=ST/2.
       IF(W2.NE.0)THEN
          SBB=(S*(T2-D2)-D5*T25)/(2.*W2)+W5
          SDD=SL1*SL6/(2.*W2)
          SEE=(S*(T2*(S+T25-W1)-W1*D2)+W5*(W1*W5-D5*(T2-W1)))/W2
          IF(SBB/SDD.GE.0)THEN
             S1P=SBB+SDD
             S1M=SEE/S1P
          ELSE
             S1M=SBB-SDD
             S1P=SEE/S1M
          ENDIF
          DD3=-W2*(S1P-TS1)*(S1M-TS1)/4.
       ELSE
          S1P=(S*(T2*(S-W5+T2-W1)-W1*W5)+W1*W5*(W1+W5-T2))/T25/(S-D5)
          DD3=-T25*(S-D5)*(S1P-TS1)/4.
       ENDIF
       ACC3=(S1P-TS1)/(S1P+TS1)
       SSB=(-R3*(D1-T1)+2.*T1*(T2+W1))/(2.*T1)
       SSD=SL3*SL7/T1
       SSE=(T2-W1)*(W4-W3)+(T2-W4+D1)*((T2-W1)*W3-(W4-W3)*W1)/T1
       IF(SSB/SSD.GE.0)THEN
          S1PP=SSB+SSD
          S1PM=SSE/S1PP
       ELSE
          S1PM=SSB-SSD
          S1PP=SSE/S1PM
       ENDIF
       DD4=-T1*(TS1-S1PP)*(TS1-S1PM)/4.
       ACC4=(TS1-S1PM)/(TS1+S1PM)
       DD5=DD1+DD3+((P12*(T1-D1)/2.-W1*P2K1)*(P2K1*(T2-D2)-W2*R3)
     +     -DELTA*(2.*P12*P2K1-W2*(T1-D1)))/P2K1
C
       RETURN
       END
       SUBROUTINE ORIENT(S,V1,V2,V3,V4,V5,DJ,NOPT,Y)
C
C  AUTHOR          : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VGVARB/E,E1,E2,E3,E4,E5,P,P3,P4,P5,CT3,ST3,CT4,ST4,CT5,
     +ST5,CP3,SP3,CP5,SP5
      COMMON/VGVARC/AL3,AL4,BE4,BE5,DE3,DE5,PP3,PP4,PP5
      COMMON/VGPZZ/W1,W2,W3,W4,W5,D1,D2,D5,D7,SL1
      COMMON/VGEXTR/TS1,TS2,T1,T2
      COMMON/VGLEVI/GRAM,DD1,DD2,DD3,DD4,DD5,DELTA,G4,SA1,SA2
      COMMON/VGDOTP/P12,P13,P14,P15,P23,P24,P25,P34,P35,P45,P1K2,P2K1
       DIMENSION Y(4)
C
C
       CALL PICKIN(S,V1,V2,V3,V4,V5,DJ,NOPT,Y)
       IF(DJ.EQ.0)RETURN
       E=DSQRT(S)
       RE=0.5/E
       E1=RE*(S+D5)
       E2=RE*(S-D5)
       P=RE*SL1
       DE3=RE*(TS2-W3+D5)
       DE5=RE*(TS1-W5-D5)
       E3=E1-DE3
       E4=DE3+DE5
       E5=E2-DE5
       IF(E4.LT.V4)THEN
          DJ=0.0
          RETURN
       ENDIF
       P3=DSQRT(E3*E3-W3)
       P4=DSQRT((E4-V4)*(E4+V4))
       P5=DSQRT(E5*E5-W5)
       PP3=DSQRT(DD1/S)/P
       PP5=DSQRT(DD3/S)/P
       ST3=PP3/P3
       ST5=PP5/P5
       IF(ST3.GT.1..OR.ST5.GT.1.)THEN
          DJ=0.0
          RETURN
       ENDIF
       CT3=DSQRT(1.-ST3*ST3)
       CT5=DSQRT(1.-ST5*ST5)
       IF(E1*E3.LT.P13)CT3=-CT3
       IF(E2*E5.GT.P25)CT5=-CT5
       AL3=ST3*ST3/(1.+CT3)
       BE5=ST5*ST5/(1.-CT5)
       IF(DD5.LT.0)THEN
          DJ=0.0
          RETURN
       ENDIF
       PP4=DSQRT(DD5/S)/P
       ST4=PP4/P4
       IF(ST4.GT.1.)THEN
          DJ=0.0
          RETURN
       ENDIF
       CT4=DSQRT(1.-ST4*ST4)
       IF(E1*E4.LT.P14)CT4=-CT4
       AL4=1.-CT4
       BE4=1.+CT4
       IF(CT4.LT.0)BE4=ST4*ST4/AL4
       IF(CT4.GE.0)AL4=ST4*ST4/BE4
       RR=DSQRT(-GRAM/S)/(P*PP4)
       SP3=RR/PP3
       SP5=-RR/PP5
       IF(DABS(SP3).GT.1..OR.DABS(SP5).GT.1.)THEN
          DJ=0.0
          RETURN
       ENDIF
       CP3=-DSQRT(1.-SP3*SP3)
       CP5=-DSQRT(1.-SP5*SP5)
       A1=PP3*CP3-PP5*CP5
       IF(DABS(PP4+PP3*CP3+CP5*PP5).LT.DABS(DABS(A1)-PP4))RETURN
       IF(A1.LT.0)CP5=-CP5
       IF(A1.GE.0)CP3=-CP3
       RETURN
       END
      DOUBLE PRECISION FUNCTION TREAT(F,X,NDIM)
C
C  AUTHOR      : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F
      DIMENSION X(10),Z(10)
      COMMON/VGB2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     + ,D(50,10),DI(50,10)
C
C
      DATA NCALL/0/
C
      IF(NCALL.EQ.0)THEN
         NCALL=1
         R=NDO
         R=R**NDIM
      ENDIF
      W=R
      DO 4 I=1,NDIM
         XX=X(I)*NDO
         J=XX
         JJ=J+1
         Y=XX-J
         IF(J.LE.0)THEN
            DD=XI(1,I)
         ELSE
            DD=XI(JJ,I)-XI(J,I)
         ENDIF
         Z(I)=XI(JJ,I)-DD*(1.-Y)
         W=W*DD
4     CONTINUE
      TREAT=W*F(Z)
C
      RETURN
      END
      SUBROUTINE SETGEN(F,NDIM,NPOIN,NPRIN,NTREAT)
C
C  AUTHOR      : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F
      DIMENSION X(10),N(10)
      COMMON/VGMAXI/MDUM,MBIN,FFMAX,FMAX(2500),NM(2500)
      COMMON/VGASIO/NINP,NOUTP
C
C
      CALL VGDAT
C
      MBIN=3
      FFMAX=0.
      SUM=0.
      SUM2=0.
      SUM2P=0.
      MAX=MBIN**NDIM
      IF(NPRIN.GE.2)WRITE(NOUTP,200)MBIN,MAX,NPOIN
      DO 5 J=1,MAX
         NM(J)=0
         FMAX(J)=0.
5     CONTINUE
      DO 1 J=1,MAX
         JJ=J-1
         DO 2 K=1,NDIM
            JJJ=JJ/MBIN
            N(K)=JJ-JJJ*MBIN
            JJ=JJJ
2        CONTINUE
         FSUM=0.
         FSUM2=0.
         DO 3 M=1,NPOIN
            DO 4 K=1,NDIM
               X(K)=(VGRAN(0.0)+N(K))/MBIN
4           CONTINUE
            IF(NTREAT.GT.0)Z=TREAT(F,X,NDIM)
            IF(NTREAT.LE.0)Z=F(X)
            IF(Z.GT.FMAX(J))FMAX(J)=Z
            FSUM=FSUM+Z
            FSUM2=FSUM2+Z*Z
3        CONTINUE
         AV=FSUM/NPOIN
         AV2=FSUM2/NPOIN
         SIG2=AV2-AV*AV
         SIG=SQRT(SIG2)
         SUM=SUM+AV
         SUM2=SUM2+AV2
         SUM2P=SUM2P+SIG2
         IF(FMAX(J).GT.FFMAX)FFMAX=FMAX(J)
         EFF=10000.
         IF(FMAX(J).NE.0)EFF=FMAX(J)/AV
         IF(NPRIN.GE.3)WRITE(NOUTP,100)J,AV,SIG,FMAX(J),EFF,
     +                                 (N(KJ),KJ=1,NDIM)
1     CONTINUE
      SUM=SUM/MAX
      SUM2=SUM2/MAX
      SUM2P=SUM2P/MAX
      SIG=SQRT(SUM2-SUM*SUM)
      SIGP=SQRT(SUM2P)
      EFF1=0.
      DO 6 J=1,MAX
         EFF1=EFF1+FMAX(J)
6     CONTINUE
      EFF1=EFF1/(MAX*SUM)
      EFF2=FFMAX/SUM
      IF(NPRIN.GE.1)WRITE(NOUTP,101)SUM,SIG,SIGP,FFMAX,EFF1,EFF2
C
100   FORMAT(I6,3X,G13.6,G12.4,G13.6,F8.2,3X,10I1)
101   FORMAT('0THE AVERAGE FUNCTION VALUE =',G14.6/
     +       ' THE OVERALL STD DEV        =',G14.4/
     +       ' THE AVERAGE STD DEV        =',G14.4/
     +       ' THE MAXIMUM FUNCTION VALUE =',G14.6/
     +       ' THE AVERAGE INEFFICIENCY   =',G14.3/
     +       ' THE OVERALL INEFFICIENCY   =',G14.3/)
200   FORMAT('1SUBROUTINE SETGEN USES A',I3,'**NDIM DIVISION'/
     + ' THIS RESULTS IN ',I7,' CUBES'/
     + ' THE PROGRAM PUT ',I5,' POINTS IN EACH CUBE TO FIND',
     + ' STARTING VALUES FOR THE MAXIMA'//)
C
      RETURN
      END
      SUBROUTINE GENERA(F,NDIM,NEVENT,NSTRAT,NTREAT)
C
C  AUTHOR      : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       EXTERNAL F
      COMMON/VGMAXI/MDUM,MBIN,FFMAX,FMAX(2500),NM(2500)
      DIMENSION X(10),N(10)
C
C
       AM=MBIN
       AMI=1./AM
       MAX=MBIN**NDIM
       NCALL=0
       MCALL=0
       NEV=0
       MEV=0
       M=0
       IF(NSTRAT.GT.0)GO TO 10
          NN=0
1         CONTINUE
          NN=NN+1
          J=VGRAN(0.0)*MAX+1.
          Y=VGRAN(0.0)*FFMAX
          NM(J)=NM(J)+1
          IF(Y.GT.FMAX(J))GO TO 1
          JJ=J-1
          DO 2 K=1,NDIM
             JJJ=JJ/MBIN
             N(K)=JJ-JJJ*MBIN
             X(K)=(VGRAN(0.0)+N(K))*AMI
             JJ=JJJ
2         CONTINUE
          IF(NTREAT.GT.0)G=TREAT(F,X,NDIM)
          IF(NTREAT.LE.0)G=F(X)
          NCALL=NCALL+1
          IF(Y.GT.G)GO TO 1
          CALL ACCEPT(1)
          NEV=NEV+1
          IF(NEV.GE.NEVENT)THEN
             RETURN
          ENDIF
          IF(G.LE.FMAX(J))GO TO 1
3         CONTINUE
          IF(NM(J).EQ.1)GO TO 7
          A=NM(J)*(G-FMAX(J))/FFMAX
          M=M+1
          T=A
4         CONTINUE
          IF(T.LT.1.)THEN
             IF(VGRAN(0.0).LE.T)GO TO 7
             T=1.
          ENDIF
          T=T-1.
          DO 6 K=1,NDIM
             X(K)=(VGRAN(0.0)+N(K))*AMI
6         CONTINUE
          IF(NTREAT.GT.0)Z=TREAT(F,X,NDIM)
          IF(NTREAT.LE.0)Z=F(X)
          MCALL=MCALL+1
          NCALL=NCALL+1
          IF(Z.LT.FMAX(J))GO TO 4
          CALL ACCEPT(1)
          MEV=MEV+1
          NEV=NEV+1
          IF(NEV.LT.NEVENT)THEN
             IF(Z.LE.G)GO TO 4
             B=NM(J)*(Z-FMAX(J))/FFMAX
             T=T+(B-A)
             A=B
             G=Z
             GO TO 4
          ENDIF
7         CONTINUE
          FMAX(J)=G
          IF(NSTRAT.GT.0)GO TO 18
          IF(G.GT.FFMAX)FFMAX=G
          IF(NEV.LT.NEVENT)GO TO 1
8         CONTINUE
          RETURN
10     CONTINUE
       NCYCLE=0
11     CONTINUE
       NCYCLE=NCYCLE+1
       H=FFMAX
       TF=NSTRAT/FFMAX
       J=0
12     CONTINUE
       JJ=J
       J=J+1
       DO 13 K=1,NDIM
          JJJ=JJ/MBIN
          N(K)=JJ-JJJ*MBIN
          JJ=JJJ
13     CONTINUE
       TOT=FMAX(J)*TF
       NM(J)=NM(J)+NSTRAT
       G=FMAX(J)
14     CONTINUE
       IF(TOT.GE.1.)GO TO 15
       IF(TOT.LE.0)GO TO 17
       IF(VGRAN(0.0).GT.TOT)GO TO 17
15     CONTINUE
       TOT=TOT-1.
       DO 16 K=1,NDIM
          X(K)=(VGRAN(0.0)+N(K))*AMI
16     CONTINUE
       IF(NTREAT.GT.0)Z=TREAT(F,X,NDIM)
       IF(NTREAT.LE.0)Z=F(X)
       NCALL=NCALL+1
       IF(Z.LT.VGRAN(0.0)*FMAX(J))GO TO 14
       IF(Z.GT.G)G=Z
       NEV=NEV+1
       CALL ACCEPT(1)
       GO TO 14
17     CONTINUE
       IF(G.GT.FMAX(J))GO TO 3
18     CONTINUE
       IF(G.GT.H)H=G
       IF(J.LT.MAX)GO TO 12
       CALL ACCEPT(0)
       FFMAX=H
       IF(NEV.LT.NEVENT)GO TO 11
C
       CALL ACCEPT(-1)
       RETURN
C
       END
      SUBROUTINE SAVE1(NDIM)
C
C   STORES VEGAS DATA (UNIT 7) FOR LATER RE-INITIALIZATION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VGB2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     + ,D(50,10),DI(50,10)
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C  AUTHOR : J. VERMASEREN
C
C
      CALL VGDAT
      REWIND LUN1
C
      WRITE(LUN1,200) NDO,IT,SI,SI2,SWGT,SCHI,
     1      ((XI(I,J),I=1,NDO),J=1,NDIM)
     2     ,((DI(I,J),I=1,NDO),J=1,NDIM)
C
200   FORMAT(2I8,4E16.10E2/(5E16.10E2))
C
      RETURN
      END
      SUBROUTINE RESTR1(NDIM)
C
C   ENTERS INITIALIZATION DATA FOR VEGAS
C
C  AUTHOR : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VGB2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     + ,D(50,10),DI(50,10)
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C
      CALL VGDAT
      REWIND LUN1
C
      READ(LUN1,200) NDO,IT,SI,SI2,SWGT,SCHI,
     1      ((XI(I,J),I=1,NDO),J=1,NDIM)
     2     ,((DI(I,J),I=1,NDO),J=1,NDIM)
C
200   FORMAT(2I8,4E16.10E2/(5E16.10E2))
C
      RETURN
      END
      SUBROUTINE SAVE2(NDIM)
C
C  AUTHOR      :J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VGMAXI/MDUM,MBIN,FFMAX,FMAX(2500),NM(2500)
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C
      CALL VGDAT
      REWIND LUN2
C
      MAX=MBIN**NDIM
      WRITE(LUN2,100)MBIN,FFMAX
      WRITE(LUN2,101)(FMAX(I),I=1,MAX)
      WRITE(LUN2,102)(NM(I),I=1,MAX)
C
100   FORMAT(I10,E16.10E2)
101   FORMAT(5E16.10E2)
102   FORMAT(8I10)
C
      RETURN
      END
      SUBROUTINE RESTR2(NDIM)
C
C  AUTHOR      : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VGMAXI/MDUM,MBIN,FFMAX,FMAX(2500),NM(2500)
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C
      CALL VGDAT
      REWIND LUN2
C
      READ(LUN2,100)MBIN,FFMAX
      MAX=MBIN**NDIM
      READ(LUN2,101)(FMAX(I),I=1,MAX)
      READ(LUN2,102)(NM(I),I=1,MAX)
C
100   FORMAT(I10,E16.10E2)
101   FORMAT(5E16.10E2)
102   FORMAT(8I10)
C
      RETURN
      END
C**********HS
C     SUBROUTINE VGSAVE(NDIM)
C
C   STORES VEGAS DATA (UNIT LUN3, DEFAULT UNIT 8) FOR LATER RE-INITIALIZ
C
C   AUTHOR : S. DE JONG
C
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     COMMON/VGB2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
C    + ,D(50,10),DI(50,10)
C     COMMON/VGMAXI/MDUM,MBIN,FFMAX,FMAX(2500),NM(2500)
C
C  SELECTION PARAMETERS FOR PHOTO PRODUCTION
C
C  FOR A DESCRIPTION OF THE VECTOR COEFFICIENTS SEE A SEPARATE SHEET AVA
C  FROM THE AUTHORS (M. POLETIEK OR S. DE JONG).
C
C     COMMON /VGSLCT/ IVGPAR(40),RVGPAR(20)
C     INTEGER IVGPAR
C     COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C
C     CALL VGDAT
C     CALL VGPDAT
C     REWIND LUN3
C
C     WRITE(LUN3)NDIM
C     WRITE(LUN3)(IVGPAR(I),I=1,20)
C     WRITE(LUN3)(RVGPAR(I),I=1,20)
C
C     WRITE(LUN3) NDO,IT,SI,SI2,SWGT,SCHI,
C    1      ((XI(I,J),I=1,NDO),J=1,NDIM)
C    2     ,((DI(I,J),I=1,NDO),J=1,NDIM)
C
C     MAX=MBIN**NDIM
C     WRITE(LUN3)MBIN,FFMAX
C     WRITE(LUN3)(FMAX(I),I=1,MAX)
C     WRITE(LUN3)(NM(I),I=1,MAX)
C
C  FOR TEST UNFORMATTED
C  IF TEST SUCCEEDS REMOVE FORMAT STATEMENTS
C
C100   FORMAT(I16)
C101   FORMAT((5I16))
C102   FORMAT((5E16.9E3))
C
C200   FORMAT(2I8,4E16.10E2/(5E16.10E2))
C
C300   FORMAT(I10,E16.10E2)
C301   FORMAT(5E16.10E2)
C302   FORMAT(8I10)
C
C     RETURN
C     END
C***********HS
C
C     SUBROUTINE VGRSTR(NDM)
C
C  RESTORES VEGAS DATA (UNIT LUN3, DEFAULT UNIT 8) FOR LATER RE-INITIALI
C
C   AUTHOR : S. DE JONG
C
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     COMMON/VGB2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
C    + ,D(50,10),DI(50,10)
C     COMMON/VGMAXI/MDUM,MBIN,FFMAX,FMAX(2500),NM(2500)
C
C  SELECTION PARAMETERS FOR PHOTO PRODUCTION
C
C  FOR A DESCRIPTION OF THE VECTOR COEFFICIENTS SEE A SEPARATE SHEET AVA
C  FROM THE AUTHORS (M. POLETIEK OR S. DE JONG).
C
C     COMMON /VGSLCT/ IVGPAR(40),RVGPAR(20)
C     INTEGER IVGPAR
C     COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C
C     CALL VGDAT
C     CALL VGPDAT
C     REWIND LUN3
C
C     READ(LUN3)NDIM
C     IF(NDM.NE.NDIM)THEN
C        PRINT*,'YOU GAVE THE WRONG DIMENSION'
C        PRINT*,'THE DIMENSION YOU GAVE WAS : ',NDM
C        PRINT*,'THE DIMENSION FOUND ON THE FILE WAS : ',NDIM
C        PRINT*,'IS THIS REALLY WHAT YOU WANT ???'
C        PRINT*,'YOUR DIMENSION WILL BE RETURNED WITH THE SAME',
C    1          'VALUE AS FOUND ON THE FILE'
C        PRINT*,'THE PROGRAM CONTINUES'
C     ENDIF
C     NDM=NDIM
C     READ(LUN3)(IVGPAR(I),I=1,20)
C     READ(LUN3)(RVGPAR(I),I=1,20)
C
C     READ(LUN3) NDO,IT,SI,SI2,SWGT,SCHI,
C    1      ((XI(I,J),I=1,NDO),J=1,NDIM)
C    2     ,((DI(I,J),I=1,NDO),J=1,NDIM)
C
C     READ(LUN3)MBIN,FFMAX
C     MAX=MBIN**NDIM
C     READ(LUN3)(FMAX(I),I=1,MAX)
C     READ(LUN3)(NM(I),I=1,MAX)
C
C  FOR TEST UNFORMATTED
C  IF TEST SUCCEEDS REMOVE FORMAT STATEMENTS
C
C100   FORMAT(I16)
C101   FORMAT((5I16))
C102   FORMAT((5E16.9E3))
C
C200   FORMAT(2I8,4E16.10E2/(5E16.10E2))
C
C300   FORMAT(I10,E16.10E2)
C301   FORMAT(5E16.10E2)
C302   FORMAT(8I10)
C
C     RETURN
C     END
      SUBROUTINE MAPLIN(T,X,TMIN,TMAX,DT)
C
C  SIMPLE LINEAR MAPPING
C
C  AUTHOR : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DT=TMIN-TMAX
      T=TMIN-DT*X
      RETURN
      END
      SUBROUTINE MAPLOG(T,X,TMIN,TMAX,DT)
C
C  SIMPLE LOGARITHMIC MAPPING
C
C  AUTHOR : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Y=TMAX/TMIN
      T=TMIN*Y**X
      DT=-T*DLOG(Y)
      RETURN
      END
      SUBROUTINE MAPINV(W2,X,W2MIN,W2MAX,DW)
C
C  SIMPLE INVERSE MAPPING
C
C  AUTHOR : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      WMIN=1./W2MIN
      WMAX=1./W2MAX
      DW=WMIN-WMAX
      W2=1./(WMAX+DW*X)
      DW=DW*W2*W2
      RETURN
      END
       SUBROUTINE ACCEPT(N)
C
C  THIS ROUTINE SHOULD IN PRINCIPLE BE PROVIDED BY THE USER.
C
C  ONCE A GOOD EVENT HAS BEEN PRODUCED BY "GENERA" ACCEPT IS CALLED AND
C  THE USER CAN PROCESS THIS EVENT.
C
C  IN THE COMMON BLOCKS "VARIAB", "VARIAC" AND "VARIAD" ARE ALL THE
C  VARIABLES CONTAINING THE PHYSICS PARAMETERS OF THE EVENT.
C
C  IF THE PATCH "LUND" IS SELECTED THE EVENTS WILL BE PUT IN THE LUND
C  COMMON BLOCK "LUJETS", AND THE PARTICLES ARE FRAGMENTED
C  USING THE LUND FRAGMENTATION ALGORITHM.
C
C         ---====>  THE LUND VERSION USED IS JETSET 6.2  <====---
C                                           ============
C  OTHERWISE THE ROUTINE IS JUST A DUMMY.
C
C  AUTHOR : S. DE JONG
C
C --- I02GAU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C --- I02GAU
      COMMON/VGVARB/E,E1,E2,E3,E4,E5,P,P3,P4,P5,CT3,ST3,CT4,ST4,CT5,
     +ST5,CP3,SP3,CP5,SP5
      COMMON/VGVARC/AL3,AL4,BE4,BE5,DE3,DE5,PP3,PP4,PP5
      COMMON/VGVARD/E6,E7,P6,P7,CT6,ST6,CT7,ST7,CP6,SP6,CP7,SP7,W
C
C
       DATA NCALL/0/
C
       NCALL=NCALL+1


       RETURN
       END

C----------------------------------------------------------------------
      FUNCTION VGRAN(S)
C
C   PICK YOUR CHOICE FOR A RANDOM NUMBER GENERATOR !
C
C  AUTHOR : J. VERMASEREN
C
      DOUBLE PRECISION VGRAN,DRANF
C
C  Modified GPS
C      VGRAN  = DRANF(S)
      DOUBLE PRECISION RVEC(1)
      call RM48(RVEC,1)
      VGRAN = RVEC(1)
      RETURN
      END

C----------------------------------------------------------------------
C GPS: intialise your favourite generator at start of a given independent 
C      sequence 
      SUBROUTINE VGINIT_RAN(iseq)
      implicit none
      integer :: iseq
      if (iseq .lt. 1) then
         write(0,*) 'VGINIT_RAN: iseq should be > 0, but was', iseq
         stop
      else
         write(0,*) 'Initialising RM48 with iseq=',iseq
         call RM48IN(iseq,0,0)
      end if
      END SUBROUTINE VGINIT_RAN
      SUBROUTINE VGINIT_RAN3(iseed)
      implicit none
      integer :: iseed(3)
      if (iseed(1) .lt. 1 .or. iseed(2) .lt. 0 .or. iseed(3) .lt.0) then
         write(0,*) 'VGINIT_RAN: iseed should be > 0,1,1, but was',iseed
         stop
      else
         call RM48IN(iseed(1),iseed(2),iseed(3))
      end if
      END SUBROUTINE VGINIT_RAN3




      SUBROUTINE RANDAT
C
C  INITIALISES THE NUMBER NCALL TO 0 TO FLAG THE FIRST CALL
C  OF THE RANDOM NUMBER GENERATOR
C
C  AUTHOR : S. DE JONG
C
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
C
      LOGICAL FIRST
C
C
      DATA FIRST /.TRUE./
C
      IF(FIRST)THEN
         FIRST=.FALSE.
         NCALL=0
      ENDIF
C
      RETURN
C
      END
      FUNCTION RANF(DUMMY)
C
C   RANDOM NUMBER FUNCTION TAKEN FROM KNUTH
C   (SEMINUMERICAL ALGORITHMS).
C   METHOD IS X(N)=MOD(X(N-55)-X(N-24),1/FMODUL)
C   NO PROVISION YET FOR CONTROL OVER THE SEED NUMBER.
C
C   RANF GIVES ONE RANDOM NUMBER BETWEEN 0 AND 1.
C   IRN55 GENERATES 55 RANDOM NUMBERS BETWEEN 0 AND 1/FMODUL.
C   IN55  INITIALIZES THE 55 NUMBERS AND WARMS UP THE SEQUENCE.
C
C  AUTHOR                     : J. VERMASEREN
C  EXTENSION TO START THROUGH : S. DE JONG
C
      PARAMETER (FMODUL=1.E-09)
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
C
C
      CALL RANDAT
C
      IF( NCALL.EQ.0 ) THEN
          CALL IN55 ( IA,234612947 )
          MCALL = 55
          NCALL = 1
      ENDIF
      IF ( MCALL.EQ.0 ) THEN
          CALL IRN55(IA)
          MCALL=55
      ENDIF
      RANF=IA(MCALL)*FMODUL
      MCALL=MCALL-1
      RETURN
      END
      DOUBLE PRECISION FUNCTION DRANF(DUMMY)
C
C  AUTHOR : J. VERMASEREN
C
      EXTERNAL RANF
C
C
      DRANF=DBLE(RANF(DUMMY))
      RETURN
      END
      SUBROUTINE IN55(IA,IX)
C
C  AUTHOR : J. VERMASEREN
C
      PARAMETER (MODULO=1000000000)
      INTEGER IA(55)
C
      IA(55)=IX
      J=IX
      K=1
      DO 10 I=1,54
         II=MOD(21*I,55)
         IA(II)=K
         K=J-K
         IF(K.LT.0)K=K+MODULO
         J=IA(II)
10    CONTINUE
      DO 20 I=1,10
         CALL IRN55(IA)
20    CONTINUE
      RETURN
      END
      SUBROUTINE IRN55(IA)
C
C  AUTHOR : J. VERMASEREN
C
      PARAMETER (MODULO=1000000000)
      INTEGER IA(55)
C
      DO 10 I=1,24
         J=IA(I)-IA(I+31)
         IF(J.LT.0)J=J+MODULO
         IA(I)=J
10    CONTINUE
      DO 20 I=25,55
         J=IA(I)-IA(I-24)
         IF(J.LT.0)J=J+MODULO
         IA(I)=J
20    CONTINUE
      RETURN
      END
      SUBROUTINE SIRN55
C
C  THIS ROUTINE SAVES THE STATE OF THE RANDOM NUMBER GENERATOR RANF
C  FROM LOGICAL UNIT LUN4
C
C  AUTHOR : S. DE JONG
C
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C
      CALL RANDAT
      CALL VGDAT
C
      WRITE(LUN4,1)NCALL,MCALL,(IA(I),I=1,55)
C
1     FORMAT(2I20,(5I16))
C
      RETURN
      END
      SUBROUTINE RIRN55
C
C  THIS ROUTINE READS THE STATE OF THE RANDOM NUMBER GENERATOR RANF
C  FROM LOGICAL UNIT LUN4
C
C  AUTHOR : S. DE JONG
C
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C
      CALL RANDAT
      CALL VGDAT
C
      READ(LUN4,1)NCALL,MCALL,(IA(I),I=1,55)
C
1     FORMAT(2I20,(5I16))
C
      RETURN
      END
      SUBROUTINE MAPT1(T,X,TMIN,TMAX,DT)
C
C  T1 MAPPING FOR PHOTO PRODUCTION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VGASIO/NINP,NOUTP
C
      DT=TMIN-TMAX
      T=TMIN-DT*X
C
      RETURN
      END
      SUBROUTINE MAPT2(T,X,TMIN,TMAX,DT)
C
C  T2 MAPPING FOR PHOTO PRODUCTION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VGASIO/NINP,NOUTP
C
      Y=TMAX/TMIN
      T=TMIN*Y**X
      DT=-T*DLOG(Y)
C
      RETURN
      END
      SUBROUTINE MAPS2(S2,X,SMIN,SMAX,DS)
C
C  S2 MAPPING FOR PHOTO PRODUCTION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VGASIO/NINP,NOUTP
C
      Y=SMAX/SMIN
      S2=SMIN*Y**X
      DS=S2*DLOG(Y)
C
      RETURN
      END
      SUBROUTINE MAPLA(X,Y,Z,U,XM,XP,D)
C
C  SELECT THE RIGHT MAPPING ROUTINE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VGASIO/NINP,NOUTP
C
      XMB=XM-Y-Z
      XPB=XP-Y-Z
      C=-4.*Y*Z
      ALP=DSQRT(XPB*XPB+C)
      ALM=DSQRT(XMB*XMB+C)
      AM=XMB+ALM
      AP=XPB+ALP
      YY=AP/AM
      ZZ=YY**U
      X=Y+Z+(AM*ZZ-C/(AM*ZZ))/2.
      AX=DSQRT((X-Y-Z)**2+C)
      D=AX*DLOG(YY)
C
      RETURN
      END
