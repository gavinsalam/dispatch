C ---------------------------------------------------------------------
C
C --- MJADE CLUSTERING
C
C ---------------------------------------------------------------------

      SUBROUTINE CLUSTERCUT(DMCUTLOC2, IFPRINT)

C --- COMMON BLOCKS FOR MC INTEGRATION
C
C --- COMMON BLOCKS FOR JET MATRIX ELEMENTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C

C---------------------------------------------------------------------

      INTEGER NMAXENT
      PARAMETER (NMAXENT=20)

      DOUBLE PRECISION POUT,PPOL
      INTEGER NENTRY,IPTYPE,IIDENT,INOUT
      COMMON /EVTRECORD/
     &       POUT(NMAXENT,4),
     &       PPOL(NMAXENT,4),
     &       NENTRY,IPTYPE(NMAXENT),IIDENT(NMAXENT),INOUT(NMAXENT)

      COMMON /PVPRECORD/
     &       PVP(NMAXENT,4),
     &       PVPPOL(NMAXENT,4)

C---------------------------------------------------------------------

CC      INTEGER NMAXENT
CC      PARAMETER (NMAXENT=20)

CC      DOUBLE PRECISION POUT,PPOL
CC      INTEGER NENTRY,IPTYPE,IIDENT,INOUT
CC      COMMON /EVTRECORD/
CC     &       POUT(NMAXENT,4),
CC     &       PPOL(NMAXENT,4),
CC     &       NENTRY,IPTYPE(NMAXENT),IIDENT(NMAXENT),INOUT(NMAXENT)

CC      COMMON /PVPRECORD/
CC     &       PVP(NMAXENT,4),
CC     &       PVPPOL(NMAXENT,4)

C --- BACKUP CLUSTERS

      DOUBLE PRECISION POUTBU,PPOLBU
      INTEGER NENTRYBU,IPTYPEBU,IIDENTBU,INOUTBU
      COMMON /EVTRECORDBU/
     &       POUTBU(NMAXENT,4),
     &       PPOLBU(NMAXENT,4),
     &       NENTRYBU,IPTYPEBU(NMAXENT),IIDENTBU(NMAXENT),
     &       INOUTBU(NMAXENT)

      COMMON /PVPRECORDBU/
     &       PVPBU(NMAXENT,4),
     &       PVPPOLBU(NMAXENT,4)

      DOUBLE PRECISION EVWGT,EVXSECT,EVERROR
      INTEGER IEVACCPT
      COMMON /USERWGT/
     &       EVWGT,
     &       EVXSECT,EVERROR,
     &       IEVACCPT

      DOUBLE PRECISION DMFACT2,DMREN2
      COMMON /USCALES/
     &       DMFACT2,DMREN2

      INTEGER IUSERINT,IUSERFLOAT,IUSERSTRING
      PARAMETER (IUSERINT=100,IUSERFLOAT=100,IUSERSTRING=100)
      DOUBLE PRECISION DUSERPAR
      INTEGER IUSERPAR
      CHARACTER*100 CUSERPAR
      COMMON /USERPAR/
     &       DUSERPAR(IUSERFLOAT),
     &       IUSERPAR(IUSERINT),
     &       CUSERPAR(IUSERSTRING)

C
      COMMON /INOUT/
     &       NOUTD
C
      COMMON /CUTS/
     &       CUT,DMM02,QQ02,WW02,XIIMIN,
     &       DTHEKPMIN,DTHEKPMAX,
     &       DTHEPIMIN,DTHEPIMAX

      COMMON /SMALLCUTS/
     &       CUTLARGE

      COMMON /ZCUTS/
     &       ZZMIN,ZZMAX

      COMMON /ZCCUTIF/
     &       IZCUT

      COMMON /RSCHEME/
     &       IREC

      COMMON /RSCHEME1/
     &       IRECFRAME

      COMMON /MORECUTS1/
     &       EELMIN

      COMMON /MORECUTS2/
     &       PTPVP

      COMMON /MORECUTS/
     &       WW2MAXC,QQ2MAXC,XHMINC,XHMAXC
C
      COMMON /MOREMORE/
     &       DMCUT2,
     &       ICUTTYPE

      COMMON /XTRAFLAGS/
     &       INOTOPIN

      COMMON /PROTON/
     &       DMPROTON

      COMMON /LABFRAME/
     &       EPEE,DTHETAKP,DTHETAP1,DTHETAP2,DTHETAP3
C
      COMMON /LABFRAME1/
     &       EPEE1
C
      COMMON /EVENT/
     &       IC,ICHOICE,NPARTONS,IPROC,IFLAOUT(4)
C
      COMMON /MCONST/
     &       PI,ZETA2,TWOPI
C
      COMMON /QUARKS/
     &       PDENSITY(13),P1DENSITY(13),CHARGE(13),
     &       NOUTFLAV
C
      COMMON /COLOUR/
     &       XNC,XCF,XNQ,XNG
C
      COMMON /COUPLINGS/
     &       ALPHAS,ALPHAEM
C
      COMMON /EMRUNC/
     &       ALPHAEMF,
     &       IAEMRUN
C
      COMMON /PSECT/
     &       AA(12,8,8),BB(6,8,8),BB1(3,2,2),CC(8,8),
     &       FF(18),FFU(18),FFD(18),
     &       FF3(18),FF3U(18),FF3D(18),
     &       FF5(18),FF5U(18),FF5D(18),
     &       FFQ(3,12,12),PXSECT(90),
     &       XS(17),XSIS(18),PROB(90),PROBSUM(90),TXSECT,DXSECT(4),
     &       IXSECT(4)
C
      COMMON /MATRIXEL/
     &       I0LOOP,I1LOOP,IREALF,IREALI,IREAL4
C
      COMMON /MATRIXEL1/
     &       INLO
C
      COMMON /MATRIXEL1/
     &       I0LOOPLT
C
      COMMON /PROCESS/
     &       KPROC(3,90),KFLAV(3,90,4)
C
      COMMON /KINEMA1/
     &       W0,RHO,
     &       PHI,PHIL,PHI1,PHI2,Z,B,E,U12,R12,
     &       EEK,EEKP,EEP0,EEP1,EEP2,EEP3,
     &       ZK,ZKP,ZP1,ZP2,ZP3,
     &       E13,E23,
     &       COSPHI1,SINPHI1,
     &       COSPHI2,SINPHI2,
     &       COSPHIL,SINPHIL,
     &       COSPHI,SINPHI,
     &       U,U1,T,DNY
C
      COMMON /KINEMA/
     &       SH,XH,Y,
     &       VPS1,VPS2,VPS3,VPS4,VPS5,VPS6,RANDN,
     &       QQ2,WW2,XII,SP,XP,
     &       SR1,SR2,SR3,
     &       DPSWGT
C
      COMMON /ISS/
     &       DV(7,3),DVS(7,3),
     &       ZPRIME,ZPRIME1,SIGMA,D,ZETA,TSTAR,T1STAR,T2STAR
C
      COMMON /INVR/
     &       S01,S02,S03,S12,S13,S23,
     &       SL00,SL01,SL02,SL03,
     &       SL10,SL11,SL12,SL13,
     &       SL
C
      COMMON /INVRR/
     &       SL0R,SL1R
C
      COMMON /YINVR/
     &       YIA,YIB,YAB
C
      COMMON /XQINT/
     &       XHMI,XHMA,XHMIP,XHMAP,
     &       QQMI,QQMA,QQMIP,QQMAP,
     &       YMI,YMA,
     &       WWMI,WWMA,
     &       ZX,ZQ,
     &       DNORMXY
C
      COMMON /DEBUG/
     &       TR(20),IDEBUG
C
      COMMON /ANALYSE/
     &       IANFLAG(10)
C
      COMMON /FSHSKIN/
     &       E1,E2,E3,ER,
     &       V12,V13,V23,VR1,VR2,VR3,
     &       EPERM1,EPERM2,EPERM3,EPERM4,
     &       XPERM12,XPERM14,
     &       QUOT1,QUOT2,QUOT4,QUOT12,QUOT14
C
      COMMON /JSINTERFACE/
     &       QMASS(7),PMASS(13),
     &       ECM,PA1,PA2,PA4,PM1,PM2,PM3,PM4,X1,X2,X4,X12,X14,
     &       KFCODE(13),IFLAV1,IFLAV2,IFLAV3,IFLAV4,
     &       IKVALID
C
      COMMON /FACTSCALE/
     &       RHOSCALE,RHOR,DNF,DNFR,
     &       ISCALE,ISCALER

      COMMON /E665CLUSTER/
     &       ICLUSTYPE

      PARAMETER (NMPART=4)
      COMMON /ETRANSV/
     &       ET(NMPART),EL(NMPART),EP(NMPART),
     &       EET(NMPART),EEP(NMPART),EXF(NMPART)
C
      COMMON /NEWSCALES/
     &       DMFMIN,DMFMAX,DMRMIN,DMRMAX
C
      COMMON /PAR01/
     &       PRENA,PRENB,PFACTA,PFACTB,PCUTA,PCUTB,
     &       DMCMIN,DMCMAX
C
      PARAMETER (NINITFLAG=1000)
      COMMON /INITBLOCK/
     &       INITCOLD(NINITFLAG)
C
      COMMON /DGR/
     &       IDGR

      COMMON /DMI/
     &       IDMI

      COMMON /HAMPEL/
     &       IHAMPEL

      COMMON /DGR1/
     &       IDOINT

      COMMON /DGR2/
     &       IOLDVERSION

      PARAMETER (MAXCLUSTER=2000)
      COMMON /CLUSTER/
     &       PP(MAXCLUSTER,4),
     &       IFREMNANT(MAXCLUSTER),
     &       IFTOUCH(MAXCLUSTER),
     &       NP(MAXCLUSTER),NCLUSTER

      COMMON /CLUSTER1/
     &       NJETS

      PARAMETER (NERRPOSS=50)

      COMMON /DEBUGCAT/
     &       I1ERRCALL(NERRPOSS),
     &       I1ERRCOUNT(NERRPOSS),
     &       I1ERROLD(NERRPOSS)

      COMMON /HELICITY/
     &       DHELFAC(4),
     &       IHELFLAG(4),IOFFHEL(4)

      COMMON /HELCONST/
     &       H30,H3R,H3X,H3G,
     &       H40,H4R,H4X,H4G,
     &       HPHI3,H2PHI3,
     &       HPHI4,H2PHI4

      COMMON /PROGRESS/
     &       IPROGRESS

      COMMON /INTPROGRESS/
     &       IINTPROGRESS

      COMMON /CHECKVIRT/
     &       VCANCEL

      COMMON /SUBTRACT4/
     &       ZM,CEFFECTIVE,Z0P0,ZPHAT,SIGMAM,ZPLOWER,ZPUPPER,
     &       ISUB4

      COMMON /XCUTS/
     &       NJVIOLATE

      COMMON /HATVAR/
     &       SHAT, XIHAT, XIHATCUT

      COMMON /M2JADE/
     &       CUTL

C --- GENERAL

      COMMON /GEN1/
     &       IFGENERAL,ILEVEL

      COMMON /XSECTS/
     &       TQ(4),TG(4),TQSUM,TGSUM,
     &       TQM1(4),TGM1(4),TQM1SUM,TGM1SUM,
     &       VTQ1(4),VTQ2(4),VTG1(4),VTG2(4),
     &       RFTQSUM,RFTGSUM,
     &       RFTQSUMRPLUS,RFTGSUMRPLUS,
     &       RITQSUMDELTA,RITQSUMREG,RITQSUMPLUS,
     &       RITPSUMDELTA,RITPSUMREG,RITPSUMPLUS,
     &       RITG1SUMDELTA,RITG1SUMREG,RITG1SUMPLUS,
     &       RITG2SUMDELTA,RITG2SUMREG,RITG2SUMPLUS,
     &       RITQSUMZPLUS,
     &       RITPSUMZPLUS,
     &       RITG1SUMZPLUS,
     &       RITG2SUMZPLUS

      COMMON /PDS/
     &       PDQ,PDG,PDP,
     &       PDQU,PDGU,PDPU

      COMMON /RSCALEDGR/
     &       SCALELOGF,SCALELOGG

      COMMON /NEWKIN/
     &       BETAMAX,ZVAR,XIPRIME,V,ALPHAMAX,R,RVAR,BVAR,
     &       XIPART1,XIPART2,XPPRIME,
     &       ZNEW,PHINEW

      PARAMETER (NLEVEL=10)
      COMMON /GWEIGHT/
     &       GDXSECT(NLEVEL,3),
     &       AGDXSECT(NLEVEL,3),
     &       UWEIGHT(NLEVEL),
     &       ICUTF(NLEVEL)

      COMMON /NEWTEST/
     &       XT01,XT02,XT03,XT12,XT13,XT23,XT1R,XT2R,XT3R

      COMMON /GLUINOMASS/
     &       GLMASS

      COMMON /MECOMMON/
     &       CROSSNUM(8,8,4,4),
     &       CROSSDEN(8,8,4,4)

      COMMON /METYPE/
     &       IFNEWME

      COMMON /GLUINOFLAV/
     &       FNEW1,FNEW2,FNEW3,FNEW4,FNEW5,FNEW6,FNEW7

      COMMON /GLUINOSUM/
     &       GLSUM
C
C --- COMMON BLOCKS USED BY VEGAS AND BY THE MAIN PROGRAM
C
C --- COMMON BLOCKS USED BY VEGAS
C
      REAL AVGT,ERRT,AVGI,ERRI
      COMMON /VGRES/
     &       AVGT,ERRT,AVGI,ERRI
C
      COMMON /VGB2/
     &       NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS,
     &       DI02GAU(50,10),DI(50,10)
C
      COMMON /VGASIO/
     &       NINP,NOUTP
C
      COMMON /BVEG3/
     &       ALPH,NDMX,MDS
C
      COMMON /CCALL/
     &       CALLS,ITMX
C
      COMMON /WEIGHT/
     &       WGTCOM
C
      COMMON /DVRES/
     &       VGINT,VGERR

      COMMON /VARMINV/
     &       BPARPR
C
C --- COMMON BLOCKS USED IN THE MAIN PROGRAM
C
      CHARACTER*32 PLOTNAME
      CHARACTER*32 JOBNAME
      COMMON /IO/
     &       PLOTNAME,JOBNAME
C
      COMMON /PROJETSTAT/
     &       INTGFLAG,NUMTOT,NUMUSE
C
      COMMON /VINT/
     &       NPOINTS,NIITMX,IALPH,NPOINTS1,NIITMX1

      COMMON /VINTSCALE/
     &       DMULT

      COMMON /WEIGHTI/
     &       TOTWGT
C
      COMMON /PDFCOMM/
     &       IFPDFLIB

      COMMON /PDFLOC/
     &       DLAMBDALOC,XMINLOC,XMAXLOC,Q2MINLOC,Q2MAXLOC,
     &       LOCNFL,LOCORD

      COMMON /PARCOUP/
     &       DLAMBDA,
     &       IARUN,IINF,
     &       ILAMBDA
C
      COMMON /PAR001/
     &       JOBNR,KPAR,IIXH,IIY,IPD,IPDSET,NINTDIM,MPLOT

      COMMON /INTFN/
     &       DSUM(100)
C
C --- ARRAY CONTAINING ADDITIONAL CUTS FOR (XH,Y)-GRID
C
      COMMON /GCUTS/
     &       ACUT(10),SUMCUT(10),DMAXCUT(10),
     &       NACUTS,NSURV(10)
C
C --- GRID IN XH, QQ2
C
      COMMON /GRID/
     &       XHMIN,XHMAX,QQ2MIN,QQ2MAX
C
C --- COMMON BLOCKS USED FOR BINNING THE DATA
C
      PARAMETER (NHISTMAX=20,NBINMAX=1000)
      CHARACTER*32 CBINNAME
      COMMON /BINCOM/
     &       BINMIN(NHISTMAX),BINMAX(NHISTMAX),
     &       BINSIZE(NHISTMAX),
     &       BINL(NHISTMAX,NBINMAX),BINU(NHISTMAX,NBINMAX),
     &       BINCENTER(NHISTMAX,NBINMAX),
     &       XBINDATA(NHISTMAX,NBINMAX),
     &       IBIN(NHISTMAX),IHUSED(NHISTMAX),
     &       CBINNAME(NHISTMAX)

      COMMON /W50510/ IFLPRTW
      COMMON /W50511/ NPTYPEW,NGROUPW,NSETW,MODEW,NFLW,LOW,TMASW
      COMMON /W50512/ QCDL4W,QCDL5W
      COMMON /W50513/ XMINW,XMAXW,Q2MINW,Q2MAXW


      DOUBLE PRECISION PPSTAR(4),P1K(4),P2K(4)

C --- COPY MOMENTA TO NEW ARRAY ((REMNANT AND) PARTON MOMENTA)
      NCLUSTER=0
      DO 10, I=1,NENTRY
         IF (IIDENT(I) .GE. 0) THEN
            NCLUSTER=NCLUSTER+1
            DO 20, K=1,4
               IF (IRECFRAME .EQ. 0) THEN
                  PP(NCLUSTER,K)=POUT(I,K)
               ELSE
                  PP(NCLUSTER,K)=PVP(I,K)
               ENDIF
20          CONTINUE
C --------- SET FLAGS
            IF (IIDENT(I) .EQ. 0) THEN
               IFREMNANT(NCLUSTER)=1
               IF (ICLUSTYPE .EQ. 1) THEN
                  IFTOUCH(NCLUSTER)=1
               ELSEIF (ICLUSTYPE .EQ. 2) THEN
                  IFTOUCH(NCLUSTER)=0
               ELSE
                  STOP
               ENDIF
            ELSE
               IFTOUCH(NCLUSTER)=1
               IFREMNANT(NCLUSTER)=0
            ENDIF
         ENDIF
10    CONTINUE

C --- THE CLUSTERING OF MOMENTA
30    CONTINUE
C ------ FIND THE SMALLEST INVARIANT MASS
         DSMALL=1.D37
         I1SMALL=0
         I2SMALL=0
         DO 40, I1=1,NCLUSTER-1
            DO 50, I2=I1+1,NCLUSTER
               IF (IFTOUCH(I1) .EQ. 1 .AND. IFTOUCH(I2) .EQ. 1) THEN
                  DO 55, K=1,4
                     P1K(K)=PP(I1,K)
                     P2K(K)=PP(I2,K)
                     IF (     IREC .EQ. -1
     &                   .OR. IREC .EQ. 0
     &                   .OR. IREC .EQ. 1) THEN
C --------------------- SUM OF FOUR-VECTORS
                        PPSTAR(K)=P1K(K)+P2K(K)
                     ELSEIF (IREC .EQ. 2) THEN
C --------------------- E0
                        PPSTAR(K)=P1K(K)+P2K(K)
                     ELSE
C --------------------- P
                        PPSTAR(K)=P1K(K)+P2K(K)
                     ENDIF
55                CONTINUE
                  IF (     IREC .EQ. -1
     &                .OR. IREC .EQ. 0) THEN
C ------------------ E
                     DPPSTAR2= PPSTAR(4)**2
     &                       -(PPSTAR(1)**2+PPSTAR(2)**2+PPSTAR(3)**2)
                  ELSEIF (IREC .EQ. 1) THEN
C ------------------ JADE
                     COSTHETA=(P1K(1)*P2K(1)
     &                        +P1K(2)*P2K(2)
     &                        +P1K(3)*P2K(3))
     &                      /DSQRT(
     &                         (P1K(1)**2+P1K(2)**2+P1K(3)**2)
     &                        *(P2K(1)**2+P2K(2)**2+P2K(3)**2)
     &                            )
                     DPPSTAR2=2.*P1K(4)*P2K(4)*(1.-COSTHETA)
                  ELSEIF (IREC .EQ. 2) THEN
C ------------------ E0
                     DPPSTAR2= PPSTAR(4)**2
     &                       -(PPSTAR(1)**2+PPSTAR(2)**2+PPSTAR(3)**2)
                  ELSE
C ------------------ P
                     DPPSTAR2= PPSTAR(4)**2
     &                       -(PPSTAR(1)**2+PPSTAR(2)**2+PPSTAR(3)**2)
                  ENDIF
                  IF (     IFREMNANT(I1) .EQ. 1
     &                .OR. IFREMNANT(I2) .EQ. 1) THEN
                     DPPSTAR2=DPPSTAR2*CUTL
                  ENDIF
                  IF (DPPSTAR2 .LT. DSMALL) THEN
                     DSMALL=DPPSTAR2
                     I1SMALL=I1
                     I2SMALL=I2
                  ENDIF
               ENDIF
50          CONTINUE
40       CONTINUE

C         WRITE(6,*) 'DSMALL=',DSMALL,DMCUTLOC2,I1SMALL, I2SMALL

         IF (IFPRINT .NE. 0) THEN
            WRITE(6,*) 'RATIO=', DSMALL/DMCUTLOC2
         ENDIF

C ------ CHECK IF COMBINATION OF TWO CLUSTERS IS REQUIRED
         IF (DSMALL .GE. DMCUTLOC2) THEN
            GOTO 100
         ENDIF
C ------ REPLACE THE TWO MOMENTA BY THE CLUSTER MOMENTUM
         DO 72, K=1,4
            P1K(K)=PP(I1SMALL,K)
            P2K(K)=PP(I2SMALL,K)
            IF (     IREC .EQ. -1
     &          .OR. IREC .EQ. 0
     &          .OR. IREC .EQ. 1) THEN
C ------------ SUM OF FOUR-VECTORS: E AND JADE
               PPSTAR(K)=P1K(K)+P2K(K)
            ELSEIF (IREC .EQ. 2) THEN
C ------------ E0
               PPSTAR(K)=P1K(K)+P2K(K)
            ELSE
C ------------ P
               PPSTAR(K)=P1K(K)+P2K(K)
            ENDIF
72       CONTINUE
C ------ RESCALE IF NECESSARY
         IF (     IREC .EQ. -1
     &       .OR. IREC .EQ. 0
     &       .OR. IREC .EQ. 1) THEN
C ---------- SUM OF FOUR-VECTORS
         ELSEIF (IREC .EQ. 2) THEN
C ---------- E0
             ALPHA=DSQRT(PPSTAR(4)**2
     &                  /(PPSTAR(1)**2
     &                   +PPSTAR(2)**2
     &                   +PPSTAR(3)**2))
             PPSTAR(1)=ALPHA*PPSTAR(1)
             PPSTAR(2)=ALPHA*PPSTAR(2)
             PPSTAR(3)=ALPHA*PPSTAR(3)
         ELSE
C ---------- P
             ALPHA=DSQRT(PPSTAR(4)**2
     &                  /(PPSTAR(1)**2
     &                   +PPSTAR(2)**2
     &                   +PPSTAR(3)**2))
             PPSTAR(4)=PPSTAR(4)/ALPHA
         ENDIF
         DO 70, K=1,4
            PP(I1SMALL,K)=PPSTAR(K)
70       CONTINUE
         IFTOUCH(I1SMALL)=1
         IF (     IFREMNANT(I1SMALL) .EQ. 1
     &       .OR. IFREMNANT(I2SMALL) .EQ. 1) THEN
            IFREMNANT(I1SMALL)=1
         ELSE
            IFREMNANT(I1SMALL)=0
         ENDIF
         DO 71, K=1,4
            PP(I2SMALL,K)=PP(NCLUSTER,K)
71       CONTINUE
         IFTOUCH(I2SMALL)=IFTOUCH(NCLUSTER)
         IFREMNANT(I2SMALL)=IFREMNANT(NCLUSTER)
         NCLUSTER=NCLUSTER-1
         IF (NCLUSTER .GE. 2) THEN
            GOTO 30
         ENDIF
100   CONTINUE

      NJETS=NCLUSTER

      GOTO 999


 999  CONTINUE
      RETURN
      END

