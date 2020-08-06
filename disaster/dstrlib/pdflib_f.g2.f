C --------------------------------------------------------------------
C
C --> INTERFACE FOR PDFLIB TO C++
C
C     FILE:              PDFLIB_F.F
C     CREATED:           16.04.1997
C     LAST MODIFICATION: 16.10.1997
C
C --------------------------------------------------------------------

C --- SET PDFLIB PARAMETERS, AND GET BACK BOUNDS

      SUBROUTINE PDFSETPAR(IPD,IPDSET,
     &                     NFL,NORDER,DLAMBDA4,XMIN,XMAX,Q2MIN,Q2MAX)

      IMPLICIT NONE

      INTEGER IPD, IPDSET
      INTEGER NFL,NORDER
      DOUBLE PRECISION DLAMBDA4,XMIN,XMAX,Q2MIN,Q2MAX

      DOUBLE PRECISION TMASW, QCDL4W, QCDL5W,
     &                 XMINW, XMAXW, Q2MINW, Q2MAXW

      INTEGER IFLPRTW,
     &        NPTYPEW, NGROUPW, NSETW, MODEW, NFLW, LOW

      COMMON /W50510/ IFLPRTW
      COMMON /W50511/ NPTYPEW,NGROUPW,NSETW,MODEW,NFLW,LOW,TMASW
      COMMON /W50512/ QCDL4W,QCDL5W
      COMMON /W50513/ XMINW,XMAXW,Q2MINW,Q2MAXW

      CHARACTER*20 PARM(20)
      DOUBLE PRECISION VAL(20)

      INTEGER I

CC      IFLPRTW = 2

      PARM(1) = 'NPTYPE'
      VAL(1)  = 1

      PARM(2) = 'NGROUP'
      VAL(2)  = IPD

      PARM(3) = 'NSET'
      VAL(3)  = IPDSET

      DO 10, I = 4, 20
         PARM(I) = ' '
         VAL(I)  = 0.0
10    CONTINUE

      CALL PDFSET(PARM,VAL)

      NFL = IABS(NFLW)

      IF (LOW .EQ. 1) THEN
         NORDER=1
      ELSE
         NORDER=2
      ENDIF

      DLAMBDA4 = QCDL4W
      XMIN     = XMINW
      XMAX     = XMAXW
      Q2MIN    = Q2MINW
      Q2MAX    = Q2MAXW

C      WRITE(6,*) NPTYPEW,NGROUPW,NSETW,MODEW
C      WRITE(6,*) 'INITIALIZATION PDFLIB DONE!'

      RETURN
      END

C ---------------------------------------------------------------------

C --- GET PARTON DISTRIBUTION FOR GIVEN X, SCALE

      SUBROUTINE PDFGET(X,SCALE,PDEN)

      IMPLICIT NONE

      DOUBLE PRECISION X, SCALE, PDEN(13)

      DOUBLE PRECISION PDF(-6:6)
      INTEGER I, J

      CALL PFTOPDG(X, SCALE, PDF)

      J = 1
      DO 10, I = -6, 6
         PDEN(J) = PDF(I)
C         WRITE(6,*) I, PDF(I)/X, J, PDEN(J)/X
         J = J+1
10    CONTINUE

      RETURN
      END

C --------------------------------------------------------------------

C --- SET PDFLIB PARAMETERS FOR ALPHA_S CALCULATION

      SUBROUTINE PDFASSETPAR(ILO, DLAMBDA4)

      IMPLICIT NONE

      INTEGER ILO
      DOUBLE PRECISION DLAMBDA4

      DOUBLE PRECISION QCDL4W, QCDL5W

      COMMON /W50512/ QCDL4W,QCDL5W

      CHARACTER*20 PARM(20)
      DOUBLE PRECISION VAL(20)

      INTEGER I

      PARM(1) = 'LO'
      VAL(1)  = ILO

      DO 10, I = 2, 20
         PARM(I) = ' '
         VAL(I)  = 0.0
10    CONTINUE

      CALL PDFSET(PARM, VAL)

      DLAMBDA4 = QCDL4W

CC      WRITE(6,*) 'L4 = ', QCDL4W
CC      WRITE(6,*) 'L5 = ', QCDL5W

      RETURN
      END

C --------------------------------------------------------------------

C --- CALCULATE ALPHA_S VIA PDFLIB ROUTINE

      FUNCTION PDFAS(SCALE)

      IMPLICIT NONE

      DOUBLE PRECISION PDFAS, SCALE, ALPHAS2

      PDFAS = ALPHAS2(SCALE)

      RETURN
      END

C ---------------------------------------------------------------------
