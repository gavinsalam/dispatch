C ---------------------------------------------------------------------
C
C --- NLO MATRIX ELEMENTS FOR (2+1) PARTON FINAL STATES
C
C ---------------------------------------------------------------------

      SUBROUTINE EVALVIRT(
     &              PAR,
     &              I0LOOPIN,I1LOOPIN,I0LOOPLTIN,
     &              XSECTQ,XSECTG
     &           )

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

      DOUBLE PRECISION PAR(11)

C ---------------------------------------------------------------------

      S01IN=PAR(1)
      S02IN=PAR(2)
      S12IN=PAR(3)
      SL01IN=PAR(4)
      XBIN=PAR(5)
      YIN=PAR(6)
      Q2IN=PAR(7)
      XIIN=PAR(8)
CC      ALPHASIN=PAR(9)
      ALPHASIN=1.D0
      DMREN2IN=PAR(10)
      DNFRIN=PAR(11)

      XH=XBIN
      Y=YIN
      QQ2=Q2IN
      XII=XIIN

      ALPHAS=ALPHASIN
      DMREN2=DMREN2IN
      DNFR=DNFRIN

      S01=S01IN
      S02=S02IN
      S12=S12IN
      SL01=SL01IN

      SH=QQ2/XH/Y
      SP=SH*XII

      Z=S01/SP/Y
      XP=1.-S12/SP/Y

      I0LOOP=I0LOOPIN
      I1LOOP=I1LOOPIN

      I0LOOPLT=I0LOOPLTIN

C ---------------------------------------------------------------------

      XNC=3.
      XCF=4./3.

      PI=3.1415927
      ZETA2=PI**2/6.

C ---------------------------------------------------------------------

      DO 10, I=1,1
         DO 20, J=1,2
            GDXSECT(I,J)=0.0
20       CONTINUE
10    CONTINUE

C .....................................................................

C --- (2+1) BORN, VIRTUAL
      ILEVEL=1

C --- BORN MATRIX ELEMENTS
      CALL DMNLO2BORN()

C --- PARTIAL CROSS SECTIONS, SUM IN TXSECT
      CALL DMNLO2A()

C .....................................................................

      DO 30, J=1,2
         DXSECT(J)=GDXSECT(1,J)
30    CONTINUE
      TXSECT=DXSECT(1)+DXSECT(2)

      XSECTQ=DXSECT(1)
      XSECTG=DXSECT(2)

      RETURN
      END

C ---------------------------------------------------------------------
C
C --- PARTIAL CROSS SECTIONS
C
C ---------------------------------------------------------------------
C
C --- (2+1) BORN
C
C ---------------------------------------------------------------------

      SUBROUTINE DMNLO2BORN()

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

      CALL HELICITYF(XP,Z,SL01)
C --- RELATIVE TO OUTGOING QUARK
      CALL BORNXPZQ(XP,Z,TQ,TG)
C --- RELATIVE TO OUTGOING GLUON
      CALL BORNXPZQ(XP,1.-Z,TQM1,TGM1)
      TQSUM=0.
      TGSUM=0.
      TQM1SUM=0.
      TGM1SUM=0.
      DO 5, I=1,4
         TQSUM=TQSUM+DHELFAC(I)*TQ(I)
         TGSUM=TGSUM+DHELFAC(I)*TG(I)
         TQM1SUM=TQM1SUM+DHELFAC(I)*TQM1(I)
         TGM1SUM=TGM1SUM+DHELFAC(I)*TGM1(I)
5     CONTINUE

      RETURN
      END

C ---------------------------------------------------------------------
C
C --- (2+1) BORN, VIRTUAL
C
C ---------------------------------------------------------------------

      SUBROUTINE DMNLO2A()

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

C --- FACTOR IN FRONT OF THE CROSS SECTION
      PF=1.
     &  *32*QQ2
     &  *(1.+(1.-Y)**2)/2./Y**2

C .....................................................................

C --- BORN TERMS

      IF (I0LOOP .EQ. 1) THEN

         GDXSECT(ILEVEL,1)=GDXSECT(ILEVEL,1)
     &            +PF*XCF*TQSUM
         GDXSECT(ILEVEL,2)=GDXSECT(ILEVEL,2)
     &            +PF*0.5*TGSUM

      ENDIF

C .....................................................................

C --- VIRTUAL CORRECTIONS

      IF (I1LOOP .EQ. 1) THEN

         CALL VIRTCORR()

         VTQ1SUM=0.
         VTQ2SUM=0.
         VTG1SUM=0.
         VTG2SUM=0.
         DO 15, I=1,4
            VTQ1SUM=VTQ1SUM+DHELFAC(I)*VTQ1(I)
            VTQ2SUM=VTQ2SUM+DHELFAC(I)*VTQ2(I)
            VTG1SUM=VTG1SUM+DHELFAC(I)*VTG1(I)
            VTG2SUM=VTG2SUM+DHELFAC(I)*VTG2(I)
15       CONTINUE

         SCALELOGF=DLOGDG(5001,DMREN2/QQ2)
     &            *(11./6.*XNC-1./3.*DNFR)

         GDXSECT(ILEVEL,1)=GDXSECT(ILEVEL,1)
     &            +PF*ALPHAS/2./PI
     &            *(
     &              XCF**2*VTQ1SUM-0.5*XNC*XCF*VTQ2SUM
     &             +SCALELOGF*XCF*TQSUM
     &             )
         GDXSECT(ILEVEL,2)=GDXSECT(ILEVEL,2)
     &            +PF*ALPHAS/2./PI
     &            *(
     &              0.5*XCF*VTG1SUM-0.25*XNC*VTG2SUM
     &             +SCALELOGF*0.5*TGSUM
     &             )

      ENDIF

C .....................................................................

      RETURN
      END

C ---------------------------------------------------------------------
C
C --- HELICITY FACTORS (2+1)
C
C ---------------------------------------------------------------------

      SUBROUTINE HELICITYF(XP1,Z1,SL01VAR)

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

      IDUMMY=I0LOOPLT
      DO 17, ICNT=1,4
         IHELFLAG(ICNT)=MOD(IDUMMY,2)
         IDUMMY=IDUMMY/2
         DHELFAC(ICNT)=0.
17    CONTINUE

      IF (     XP1    .LT. 1.D-10
     &    .OR. 1.-XP1 .LT. 1.D-10
     &    .OR. 1.-Y   .LT. 1.D-10
     &    .OR. Z1     .LT. 1.D-10
     &    .OR. 1.-Z1  .LT. 1.D-10
     &   ) THEN
         DCOSPHIL=1.
      ELSE
         DCOSPHIL= ((1.-XP1*Y)*Z1+(1.-XP1)*(1.-2.*Z1)-XP1*Y*SL01VAR/QQ2)
     &            /(2.*DSQRTDG(1302,XP1*(1.-XP1)*(1.-Y)*Z1*(1.-Z1)))
      ENDIF

      IF (IHELFLAG(1) .EQ. 1) THEN
         DHELFAC(1)=1.D0
      ENDIF

      IF (IHELFLAG(2) .EQ. 1) THEN
         DHELFAC(2)=(4*(1.-Y)+1.+(1.-Y)**2)/(1.+(1.-Y)**2)
      ENDIF

      IF (IHELFLAG(3) .EQ. 1) THEN
         DHELFAC(3)= 2.*(2.-Y)*DSQRTDG(1301,1.-Y)
     &              *DCOSPHIL/(1.+(1.-Y)**2)
      ENDIF

      IF (IHELFLAG(4) .EQ. 1) THEN
         DHELFAC(4)= 2.*2.*(1.-Y)
     &              *(2.*DCOSPHIL**2-1.)/(1.+(1.-Y)**2)
      ENDIF

      RETURN
      END

C ---------------------------------------------------------------------
C
C --- BORN TERMS AS A FUNCTION OF XP AND ZQ
C
C ---------------------------------------------------------------------

      SUBROUTINE BORNXPZQ(XP1,Z1,TQ1,TG1)

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

      DIMENSION TQ1(4),TG1(4)

      TQ1(1)= (1.-XP1)/(1.-Z1)+(1.-Z1)/(1.-XP1)
     &       +2.*XP1*Z1/(1.-XP1)/(1.-Z1)

      TG1(1)= (1.-Z1)/Z1+Z1/(1.-Z1)
     &       -2.*XP1*(1.-XP1)/Z1/(1.-Z1)

      TQ1(2)=2.*XP1*Z1

      TG1(2)=4.*XP1*(1.-XP1)

      TQ1(3)= -2.*DSQRTDG(1303,Z1*XP1/(1.-Z1)/(1.-XP1))
     &       *((1.-XP1)*(1.-Z1)+XP1*Z1)

      TG1(3)= -2.*DSQRTDG(1303,XP1*(1.-XP1)/Z1/(1.-Z1))
     &       *(1.-2.*XP1)*(1.-2.*Z1)

      TQ1(4)=XP1*Z1

      TG1(4)=2.*XP1*(1.-XP1)

      RETURN
      END

C ---------------------------------------------------------------------
C
C --- VIRTUAL CORRECTIONS AS A FUNCTION OF S01,S02,S12
C
C ---------------------------------------------------------------------

      SUBROUTINE VIRTCORR()

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

      YIA=S01/QQ2
      YIB=S02/QQ2
      YAB=S12/QQ2

C .....................................................................

C --- QUARK INITIATED

      Y13=YAB
      Y03=YIB

      Y01=1.D0+Y13-Y03

      DLOG01=DLOGDG(4001,Y01)
      DLOG03=DLOGDG(4004,Y03)
      DLOG13=DLOGDG(4003,Y13)
      RFUNCT0103=RFUNCT(Y01,Y03)
      RFUNCT0113=RFUNCT(Y01,-Y13)
      RFUNCT0313=RFUNCT(Y03,-Y13)

C --- BORN-TERM
      T3JETQ = Y13/Y03 + Y03/Y13 + 2.D0*Y01/(Y13*Y03)

C --- CORRECTIONS

      ADDCON=T3JETQ*(-2.*ZETA2-DLOG01**2-8.)
     & -( 4.*DLOG01*(2.*Y01/(-Y13+Y03)+Y01**2/(-Y13+Y03)**2)
     &   +DLOG13
     &    *((4.*Y01-2.*Y13)/(Y01+Y03)+Y13*Y03/(Y01+Y03)**2)
     &   +DLOG03
     &    *((4.*Y01+2.*Y03)/(Y01-Y13)+Y13*Y03/(Y01-Y13)**2)
     &   +2.*(Y01**2+(Y01+Y03)**2)/Y13/Y03*RFUNCT0113
     &   +2.*(Y01**2+(Y01-Y13)**2)/Y13/Y03*RFUNCT0103
     &   +Y01*(4./(-Y13+Y03)+1./(Y01-Y13)+1./(Y01+Y03))
     &   +Y01/Y13-Y01/Y03+Y13/Y03+Y03/Y13
     &  )
      EDVIRT1Q=ADDCON
C      EDVIRT1Q=T3JETQ
C      EDVIRT1Q=DLOG13+DLOG03+DLOG01
C     &        +( Y01/Y03+Y03/Y13+Y13/Y01
C     &          +Y03/Y01+Y13/Y03+Y01/Y13)
C     &        *(DLOG13**2+DLOG03**2+DLOG01**2)

      ADDCON=T3JETQ*(( 2.*ZETA2-DLOG01**2
     &               +(DLOG13**2-PI**2)
     &               +DLOG03**2)+2.*RFUNCT0313)
     & -( DLOG13*2.*Y13/(Y01+Y03)
     &   +DLOG03*(-2.)*Y03/(Y01-Y13)
     &   +4.*DLOG01*(Y01**2/(-Y13+Y03)**2+2.*Y01/(-Y13+Y03))
     &   -2.*(-Y13/Y03-Y03/Y13-Y01/Y13+Y01/Y03-2.*Y01/(-Y13+Y03))
     &   +2.*RFUNCT0113*(Y01**2+(Y01+Y03)**2)/Y13/Y03
     &   +2.*RFUNCT0103*(Y01**2+(Y01-Y13)**2)/Y13/Y03
     &  )
      EDVIRT2Q=ADDCON

C --- BORN-TERM
      T3JETQ = 4.*XP**2*(0.5*Y01)

C --- CORRECTIONS

      ADDCON=T3JETQ*(
     &   (-2.*ZETA2
     &   -DLOG01**2-8)
     &   +8
     &   +2.*DLOG01*(-Y03/(Y13-Y03)+Y01*Y03/(Y13-Y03)**2)
     &   +2.*DLOG13
     &   +2.*RFUNCT0113*Y01/Y03
     &   -2.*RFUNCT0103
     &   -7.-2.*Y13/(Y13-Y03)
     &              )
      EDVIRT1LQ=ADDCON

      ADDCON=T3JETQ*(
     &   (2.*ZETA2
     &   -DLOG01**2
     &   +(DLOG13**2-PI**2)
     &   +DLOG03**2)
     &   +2.*DLOG01*(-Y03/(Y13-Y03)+Y01*Y03/(Y13-Y03)**2)
     &   +2.*RFUNCT0113*Y01/Y03
     &   -2.*RFUNCT0103
     &   +2.*RFUNCT0313
     &   -2.*Y03/(Y13-Y03)-2.
     &              )
      EDVIRT2LQ=ADDCON

C --- BORN-TERM
      T3JETQ = 4.*XP**2*(0.5*Y01)

C --- CORRECTIONS

      ADDCON=T3JETQ*(
     &   -2.*ZETA2
     &   -DLOG01**2
     &   +2.*DLOG01*(Y13/(Y13-Y03)-Y01*Y13/(Y13-Y03)**2)
     &   +2.*DLOG03
     &   +2.*RFUNCT0103*(-Y01/Y13)
     &   -2.*RFUNCT0113
     &   -7.+2.*Y03/(Y13-Y03)
     &              )
      EDVIRT1RQ=ADDCON

      ADDCON=T3JETQ*(
     &   2.*ZETA2
     &   -DLOG01**2
     &   +(DLOG13**2-PI**2)
     &   +DLOG03**2
     &   +2.*DLOG01*(Y13/(Y13-Y03)-Y01*Y13/(Y13-Y03)**2)
     &   +2.*RFUNCT0103*(-Y01/Y13)
     &   -2.*RFUNCT0113
     &   +2.*RFUNCT0313
     &   +2.*Y13/(Y13-Y03)-2.
     &              )
      EDVIRT2RQ=ADDCON

C --- BORN-TERM
      T3JETQ = 4.*XP**2*Y01

C --- CORRECTIONS

      ADDCON=T3JETQ*(
     &   -2.*ZETA2
     &   +2.*DLOG01
     &   +0.5*DLOG13*(2.*Y03/(Y01+Y03)-Y03*Y13/(Y01+Y03)**2)
     &   +0.5*DLOG03*(2.*Y13/(-Y01+Y13)-Y03*Y13/(-Y01+Y13)**2)
     &   -DLOG01**2
     &   -RFUNCT0113
     &   -RFUNCT0103
     &   -0.5*(14.+Y01/(Y01-Y13)+Y01/(Y01+Y03))
     &              )
      EDVIRT1XQ=ADDCON

      ADDCON=T3JETQ*(
     &   2.*ZETA2
     &   -DLOG01**2
     &   +(DLOG13**2-PI**2)
     &   +DLOG03**2
     &   +2.*DLOG01
     &   +DLOG13*(-Y13/(Y01+Y03))
     &   +DLOG03*(Y03/(Y01-Y13))
     &   -RFUNCT0113
     &   -RFUNCT0103
     &   +2.*RFUNCT0313
     &              )
      EDVIRT2XQ=ADDCON

C .....................................................................

C --- GLUON INITIATED

      Y01=YIA
      Y02=YIB

      Y12 = -1.D0+Y01+Y02

      DLOG01=DLOGDG(4011,Y01)
      DLOG02=DLOGDG(4012,Y02)
      DLOG12=DLOGDG(4013,Y12)
      RFUNCT0102=RFUNCT(Y01,Y02)
      RFUNCT0112=RFUNCT(Y01,-Y12)
      RFUNCT0212=RFUNCT(Y02,-Y12)

C --- BORN-TERM
      T3JETG = Y01/Y02 + Y02/Y01 - 2.D0*Y12/(Y01*Y02)

C --- CORRECTIONS

      ADDCON=T3JETG*(-2.*ZETA2
     &               -(DLOG12**2-PI**2)-8.)
     & +4.*DLOG12*(-2.*Y12/(Y01+Y02)+Y12**2/(Y01+Y02)**2)
     & +DLOG01
     &  *((-4.*Y12+2.*Y01)/(-Y12+Y02)-Y01*Y02/(-Y12+Y02)**2)
     & +DLOG02
     &  *((-4.*Y12+2.*Y02)/(-Y12+Y01)-Y01*Y02/(-Y12+Y01)**2)
     & -2.*(Y12**2+(-Y12+Y02)**2)/Y01/Y02*RFUNCT0112
     & -2.*(Y12**2+(-Y12+Y01)**2)/Y01/Y02*RFUNCT0212
     & -Y12*(4./(Y01+Y02)+1./(-Y12+Y01)+1./(-Y12+Y02))
     & +Y12/Y01+Y12/Y02-Y01/Y02-Y02/Y01
      EDVIRT1G=ADDCON

      ADDCON=T3JETG*((2.*ZETA2
     &               -(DLOG12**2-PI**2)
     &               +DLOG01**2
     &               +DLOG02**2)+2.*RFUNCT0102)
     & +DLOG01*(-2.)*Y01/(-Y12+Y02)
     & +DLOG02*(-2.)*Y02/(-Y12+Y01)
     & +4.*DLOG12*(Y12**2/(Y01+Y02)**2-2.*Y12/(Y01+Y02))
     & -2.*(Y01/Y02+Y02/Y01-Y12/Y01-Y12/Y02+2.*Y12/(Y01+Y02))
     & -2.*RFUNCT0112*(Y12**2+(-Y12+Y02)**2)/Y01/Y02
     & -2.*RFUNCT0212*(Y12**2+(-Y12+Y01)**2)/Y01/Y02
      EDVIRT2G=ADDCON

C --- BORN-TERM
      T3JETG = 4.*XP**2*(Y12)

C --- CORRECTIONS

      ADDCON=T3JETG*(
     & (-2.*ZETA2
     & -(DLOG12**2-PI**2)-8)
     & +8
     & +2.*DLOG12
     & +0.5*DLOG01*(-2.*Y02/(Y12-Y02)+Y01*Y02/(Y12-Y02)**2)
     & +0.5*DLOG02*(-2.*Y01/(Y12-Y01)+Y01*Y02/(Y12-Y01)**2)
     & -RFUNCT0112
     & -RFUNCT0212
     & -0.5*(14.+Y12/(Y12-Y01)+Y12/(Y12-Y02))
     &              )
      EDVIRT1LG=ADDCON

      ADDCON=T3JETG*(
     & (2.*ZETA2
     & -(DLOG12**2-PI**2)
     & +DLOG01**2
     & +DLOG02**2)
     & +2.*DLOG12
     & +DLOG01*(-Y01/(Y12-Y02))
     & +DLOG02*(-Y02/(Y12-Y01))
     & -RFUNCT0112
     & -RFUNCT0212
     & +2.*RFUNCT0102
     &              )
      EDVIRT2LG=ADDCON

C --- BORN-TERM
      T3JETG = 4.*XP**2*(0.5*Y12)

C --- CORRECTIONS

      ADDCON=T3JETG*(
     & -2.*ZETA2
     & -(DLOG12**2-PI**2)
     & +2.*DLOG12*(Y01/(Y01+Y02)-Y01*Y12/(Y01+Y02)**2)
     & +2.*DLOG02
     & +2.*RFUNCT0212*(-Y12/Y01)
     & -2.*RFUNCT0112
     & -7.-2.*Y02/(Y01+Y02)
     &              )
      EDVIRT1RG=ADDCON

      ADDCON=T3JETG*(
     & 2.*ZETA2
     & -(DLOG12**2-PI**2)
     & +DLOG01**2
     & +DLOG02**2
     & +2.*DLOG12*(Y01/(Y01+Y02)-Y01*Y12/(Y01+Y02)**2)
     & +2.*RFUNCT0212*(-Y12/Y01)
     & -2.*RFUNCT0112
     & +2.*RFUNCT0102
     & +2.*Y01/(Y01+Y02)-2.
     &              )
      EDVIRT2RG=ADDCON

C --- BORN-TERM
      T3JETG = 4.*XP**2*(0.5*Y12)

C --- CORRECTIONS

      ADDCON=T3JETG*(
     & -2.*ZETA2
     & -(DLOG12**2-PI**2)
     & +2.*DLOG12*(Y02/(Y01+Y02)-Y02*Y12/(Y01+Y02)**2)
     & +2.*DLOG01
     & +2.*RFUNCT0112*(-Y12/Y02)
     & -2.*RFUNCT0212
     & -7.-2.*Y01/(Y01+Y02)
     &              )
      EDVIRT1XG=ADDCON

      ADDCON=T3JETG*(
     & 2.*ZETA2
     & -(DLOG12**2-PI**2)
     & +DLOG01**2
     & +DLOG02**2
     & +2.*DLOG12*(Y02/(Y01+Y02)-Y02*Y12/(Y01+Y02)**2)
     & +2.*RFUNCT0112*(-Y12/Y02)
     & -2.*RFUNCT0212
     & +2.*RFUNCT0102
     & +2.*Y02/(Y01+Y02)-2.
     &              )
      EDVIRT2XG=ADDCON

C .....................................................................

C --- LINEAR COMBINATIONS (HELICITY CROSS SECTIONS)

C --- PROJECTION OPERATORS
      COMDEN=4.*XP*(1-XP)**2*Z**2*(1-Z)**2
      H30=(-Z*(1.-Z)-XP*(1.-XP)+2.*XP*(1.-XP)*Z*(1.-Z))/COMDEN
      H3R=(XP+Z-2.*XP*Z)/COMDEN
      H3X=((1.-XP)*(1.-Z)+XP*Z)/COMDEN
      H3G=(-1)*(2.*XP*(1.-XP)*Z*(1.-Z))/COMDEN
      H40=( 1.5*XP**2*(1.-XP)**2
     &     +1.5*Z**2*(1.-Z)**2
     &     +0.5*(1.-XP-Z)**2
     &         *(2.*XP*(1.-XP)+2.*Z*(1.-Z)-(Z-XP)**2)
     &    )/COMDEN
      H4R=(-XP*(1.-XP)-Z*(1.-Z)+3.*XP*(1.-XP)*Z*(1.-Z))/COMDEN
      H4X=-((1.-XP)**2*(1.-Z)**2+XP**2*Z**2+XP*(1.-XP)*Z*(1.-Z))/COMDEN
      H4G=(-1)*(-2.*XP*(1.-XP)*Z*(1.-Z)
     &          *((1.-XP)*(1.-Z)+XP*Z))/COMDEN

      HPHI3  = -2.*DSQRTDG(8010,XP*(1.-XP)*Z*(1.-Z))/XP
     &            *((1.-XP)*(1.-Z)+XP*Z)
      HPHI4  = -2.*DSQRTDG(8010,XP*(1.-XP)*Z*(1.-Z))/XP
      H2PHI3 = Z*(1.-Z)*(1.-XP)
      H2PHI4 = 0.

      EDVIRT1H3Q  = 0.
     &             +H30*EDVIRT1LQ
     &             +H3R*EDVIRT1RQ
     &             +H3X*EDVIRT1XQ
     &             +H3G*EDVIRT1Q

      EDVIRT1H3G  = 0.
     &             +H30*EDVIRT1LG
     &             +H3R*EDVIRT1RG
     &             +H3X*EDVIRT1XG
     &             +H3G*EDVIRT1G

      EDVIRT1H4Q  = 0.
     &             +H40*EDVIRT1LQ
     &             +H4R*EDVIRT1RQ
     &             +H4X*EDVIRT1XQ
     &             +H4G*EDVIRT1Q

      EDVIRT1H4G  = 0.
     &             +H40*EDVIRT1LG
     &             +H4R*EDVIRT1RG
     &             +H4X*EDVIRT1XG
     &             +H4G*EDVIRT1G

      EDVIRT2H3Q  = 0.
     &             +H30*EDVIRT2LQ
     &             +H3R*EDVIRT2RQ
     &             +H3X*EDVIRT2XQ
     &             +H3G*EDVIRT2Q

      EDVIRT2H3G  = 0.
     &             +H30*EDVIRT2LG
     &             +H3R*EDVIRT2RG
     &             +H3X*EDVIRT2XG
     &             +H3G*EDVIRT2G

      EDVIRT2H4Q  = 0.
     &             +H40*EDVIRT2LQ
     &             +H4R*EDVIRT2RQ
     &             +H4X*EDVIRT2XQ
     &             +H4G*EDVIRT2Q

      EDVIRT2H4G  = 0.
     &             +H40*EDVIRT2LG
     &             +H4R*EDVIRT2RG
     &             +H4X*EDVIRT2XG
     &             +H4G*EDVIRT2G

      EDVIRT1PHIQ  = 0.
     &              +HPHI3*EDVIRT1H3Q
     &              +HPHI4*EDVIRT1H4Q

      EDVIRT12PHIQ = 0.
     &              +H2PHI3*EDVIRT1H3Q
     &              +H2PHI4*EDVIRT1H4Q

      EDVIRT1PHIG  = 0.
     &              +HPHI3*EDVIRT1H3G
     &              +HPHI4*EDVIRT1H4G

      EDVIRT12PHIG = 0.
     &              +H2PHI3*EDVIRT1H3G
     &              +H2PHI4*EDVIRT1H4G

      EDVIRT2PHIQ  = 0.
     &              +HPHI3*EDVIRT2H3Q
     &              +HPHI4*EDVIRT2H4Q

      EDVIRT22PHIQ = 0.
     &              +H2PHI3*EDVIRT2H3Q
     &              +H2PHI4*EDVIRT2H4Q

      EDVIRT2PHIG  = 0.
     &              +HPHI3*EDVIRT2H3G
     &              +HPHI4*EDVIRT2H4G

      EDVIRT22PHIG = 0.
     &              +H2PHI3*EDVIRT2H3G
     &              +H2PHI4*EDVIRT2H4G

      VTQ1(1)=EDVIRT1Q
      VTQ1(2)=EDVIRT1LQ
      VTQ1(3)=EDVIRT1PHIQ
      VTQ1(4)=EDVIRT12PHIQ

      VTQ2(1)=EDVIRT2Q
      VTQ2(2)=EDVIRT2LQ
      VTQ2(3)=EDVIRT2PHIQ
      VTQ2(4)=EDVIRT22PHIQ

      VTG1(1)=EDVIRT1G
      VTG1(2)=EDVIRT1LG
      VTG1(3)=EDVIRT1PHIG
      VTG1(4)=EDVIRT12PHIG

      VTG2(1)=EDVIRT2G
      VTG2(2)=EDVIRT2LG
      VTG2(3)=EDVIRT2PHIG
      VTG2(4)=EDVIRT22PHIG

      RETURN
      END

C ---------------------------------------------------------------------
C
      FUNCTION DSQRTDG(ILOC,XIN1)
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
      IF (XIN1 .LT. 0.D0) THEN
         WRITE(6,*) 'ERROR SQRT'
         WRITE(6,*) '**ARGUMENT:',XIN1
         WRITE(6,*) '**LOCATION:',ILOC
         DSQRTDG=0.
      ELSE
         DSQRTDG=DSQRT(XIN1)
      ENDIF
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
      FUNCTION DLOGDG(ILOC,XIN1)
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
      IF (XIN1 .LE. 0.D0) THEN
         WRITE(6,*) 'ERROR LN'
         WRITE(6,*) '**ARGUMENT:',XIN1
         WRITE(6,*) '**LOCATION:',ILOC
         DLOGDG=-1.D30
      ELSE
         DLOGDG=DLOG(XIN1)
      ENDIF
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C === THE "R"-FUNCTION
C
C     THIS FUNCTION IS ONLY CORRECT, IF AT LEAST ONE ARGUMENT IS NOT
C     SMALLER THAN 0.
C
C ---------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION RFUNCT(XPAR,YPAR)

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

      COMPLEX*16 CS2

      RFUNCT = DLOGDG(4021,DABS(XPAR))*DLOGDG(4022,DABS(YPAR))
     &        -DLOGDG(4023,DABS(XPAR))*DLOGDG(4024,DABS(1.D0-XPAR))
     &        -DLOGDG(4025,DABS(YPAR))*DLOGDG(4026,DABS(1.D0-YPAR))
     &        -DREAL(CS2(DCMPLX(XPAR,1.D-10)))
     &        -DREAL(CS2(DCMPLX(YPAR,1.D-10)))
     &        +ZETA2

      RETURN

      END

C ---------------------------------------------------------------------
