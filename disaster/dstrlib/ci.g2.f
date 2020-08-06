C ---------------------------------------------------------------------
C
C INTERFACE TO THE OLD CLUSTER ROUTINE
C
C ---------------------------------------------------------------------

      SUBROUTINE CI(PIN, NPA, DCSCALE2, NJETSOUT, ETAIN, ICUT, IFLAGS)

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


      DOUBLE PRECISION PIN(4,7)
      INTEGER ICUT

      PTMIN=1.0
      PTMAX=300.0
      RAPMIN=-3.5
      RAPMAX=3.5
      RATIO = 29.8845
      ECMLOC = 300

      CALL CIFULL(
     &        PIN,
     &        NPA,
     &        DCSCALE2,
     &        NJETSOUT,
     &        ETAIN,
     &        ICUT,
     &        IFLAGS,
     &        PTMIN,
     &        PTMAX,
     &        RAPMIN,
     &        RAPMAX,
     &        RATIO,
     &        ECMLOC
     &     )

      RETURN
      END

C ---------------------------------------------------------------------
C
C INTERFACE TO THE OLD CLUSTER ROUTINE
C
C ---------------------------------------------------------------------

      SUBROUTINE CIFULL(
     &              PIN,
     &              NPA,
     &              DCSCALE2,
     &              NJETSOUT,
     &              ETAIN,
     &              ICUT,
     &              IFLAGS,
     &              PTMIN,
     &              PTMAX,
     &              RAPMIN,
     &              RAPMAX,
     &              RATIO,
     &              ECMLOC
     &           )

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


      DOUBLE PRECISION PIN(4,7)
      INTEGER ICUT

      CUTL=1.0

      IRECFRAME=1
      IREC=-1
      ICLUSTYPE=1

      NENTRY=NPA

      IIDENT(1)=0
      DO K=1,4
         PVP(1,K)=(1-ETAIN)/ETAIN*PIN(K,1)
      ENDDO

      DO I=2,NPA
      IIDENT(I)=I-1
         DO K=1,4
            PVP(I,K)=PIN(K,I)
         ENDDO
      ENDDO

C      WRITE(6,*) DCSCALE2, ETAIN
C      DO I=1,NENTRY
C         WRITE(6,11) PVP(I,1),PVP(I,2),PVP(I,3),PVP(I,4),IIDENT(I)
C      ENDDO
C11    FORMAT (4E16.6,1I3)

      CALL CLUSTERCUT(DCSCALE2, IFLAGS)

      IF (IFLAGS .NE. 0) THEN
         DO I=1,NJETS
            WRITE(6,111) PP(I,1),PP(I,2),PP(I,3),PP(I,4),IIDENT(I),
     &                   IFREMNANT(I)
         ENDDO
111      FORMAT (4E16.6,1I3,1I3)
      ENDIF
C      WRITE(6,112) PIN(1,1), PIN(2,1), PIN(3,1), PIN(4,1)
C      WRITE(6,112) PIN(1,6), PIN(2,6), PIN(3,6), PIN(4,6)
C112   FORMAT (4E16.6)

      NJETSOUT=NJETS

C --- NOW CHECK RAPIDITY AND PT CUTS
CC      PTMIN=1.0
CC      PTMAX=300.0
CC      RAPMIN=-3.5
CC      RAPMAX=3.5
CC      RATIO = 29.8845
CC      ECM = 300
      ICUT=1
      DO I=1, NJETSOUT
C         WRITE(6,*) I, IFREMNANT(I)
         IF (IFREMNANT(I) .EQ. 0) THEN
C --------- CALCULATE INNER PRODUCTS
C --------- LEPTON
            RL = PP(I,4) * PIN(4,6)
     &          -PP(I,1) * PIN(1,6)
     &          -PP(I,2) * PIN(2,6)
     &          -PP(I,3) * PIN(3,6)
C --------- INCIDENT
            RP = PP(I,4) * PIN(4,1)
     &          -PP(I,1) * PIN(1,1)
     &          -PP(I,2) * PIN(2,1)
     &          -PP(I,3) * PIN(3,1)
            RP = RP / ETAIN
            RR2 = PP(I,4) ** 2
     &           -PP(I,1) ** 2
     &           -PP(I,2) ** 2
     &           -PP(I,3) ** 2
            RE =   2 * DSQRT(RATIO) / ECMLOC / (RATIO+RATIO)
     &               * (RP + RATIO * RL)
            RZ = - 2 * DSQRT(RATIO) / ECMLOC / (RATIO+RATIO)
     &               * (RP - RATIO * RL)
            RT2 = RE**2 - RZ**2 - RR2
            RT = DSQRT(RT2)
            RABS2 = RE**2 - RR2
            RABS = DSQRT(RABS2)
            COSTH = - RZ / RABS
            V = 0.5*(1-COSTH)
C            WRITE(6,*) V, 1-V
            IF (V .GT. 1.D-40 .AND. 1.D0-V .GT. 1.D-40) THEN
               PRAP = 0.5 * DLOG((1-V)/V)
            ELSEIF (V .LE. 0.5 .AND. V .GT. -1.D-2) THEN
               PRAP = 20
            ELSEIF (V .GE. 0.5 .AND. 1.D0-V .GT. -1.D-2) THEN
               PRAP = -20
            ELSE
               WRITE(6,*) 'RAPIDITY OUT OF RANGE'
               STOP
            ENDIF
            IF (RT .LT. PTMIN .OR. RT .GT. PTMAX) THEN
               ICUT=0
            ENDIF
            IF (PRAP .LT. RAPMIN .OR. PRAP .GT. RAPMAX) THEN
               ICUT=0
            ENDIF
            IF (IFLAGS .EQ. 1) THEN
               WRITE(6,223) I, RT, PRAP
223            FORMAT(1I3,2E16.6)
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END

C ---------------------------------------------------------------------

C --- INTERFACE VIA AN ARRAY
C     ORDER: (P1..P3, ENERGY) ... () ... ()
C     REQUIRES THE MOMENTUM FRACTION, TOO!

      SUBROUTINE CIARRAY(PARRAY, NPA, DCSCALE2, NJETSOUT, ETAIN, ICUT,
     &                   IFLAGS)

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


      DOUBLE PRECISION PARRAY(*)
      DOUBLE PRECISION PIN(4,7)

      INTEGER ICUT

C      WRITE(6,*) DCSCALE2, ETAIN

C --- COPY INTO PIN

      INDE=1
      DO I=1, NPA
         DO K=1,4
            PIN(K,I)=PARRAY(INDE)
            INDE=INDE+1
         ENDDO
C         WRITE(6,*) NPA, I, PIN(1,I), PIN(2,I), PIN(3,I), PIN(4,I)
      ENDDO
      DO K=1,4
         PIN(K,6)=PARRAY(INDE)
         INDE=INDE+1
      ENDDO
C      WRITE(6,*) PIN(1,6), PIN(2,6), PIN(3,6), PIN(4,6)

      CALL CI(PIN, NPA, DCSCALE2, NJETSOUT, ETAIN, ICUT, IFLAGS)

      RETURN
      END

C ---------------------------------------------------------------------
