C --------------------------------------------------------------------
C
C --> FORTRAN INTERFACES (USER ROUTINE) FOR DISASTER++
C
C     FILE:              DISASTER_F.F
C     CREATED:           01.10.1997
C     LAST MODIFICATION: 16.12.1997
C
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C
C --> FORTRAN MAIN PROGRAM
C
C --------------------------------------------------------------------

      PROGRAM DISASTER

      IMPLICIT NONE

      INTEGER IFLAG
      INTEGER PPCMAIN

      WRITE(6,*) 'FORTRAN PROGRAM START...'

C --- PPCMAIN STARTS THE C++ INTERFACE
      IFLAG = 1
C      IFLAG = PPCMAIN()
      IF (IFLAG .EQ. 0) THEN
         GOTO 110
      ENDIF

C --- OPEN THE LOG FILE
      OPEN(FILE='LOGFILE', UNIT=42, STATUS='UNKNOWN')

C --- START DISASTER++
      CALL DISASTER_CA()

C --- SET SOME PARAMETERS

      CALL DISASTER_CD('ECM', DBLE(300.0))

      CALL DISASTER_CI('LEPTON_INTEGRATION', 1)

      CALL DISASTER_CD('XBMIN', DBLE(0.005))
      CALL DISASTER_CD('XBMAX', DBLE(0.01))

      CALL DISASTER_CD('YMIN', DBLE(0.01))
      CALL DISASTER_CD('YMAX', DBLE(0.03))

      CALL DISASTER_CI('PROCESS_INDEX', 2)
      CALL DISASTER_CI('NUMBER_OF_FINAL_STATE_PARTONS_IN_BORN_TERM',
     &                 2)

      CALL DISASTER_CI('POINTS', 100)
      CALL DISASTER_CD('FACT_PREP',  DBLE(0.5))
      CALL DISASTER_CD('FACT_FINAL', DBLE(1.5))
      CALL DISASTER_CI('ITERATIONS', 10)

      CALL DISASTER_CC('RUN_MC', 0)

C --- END DISASTER++
      CALL DISASTER_CB()

C --- CLOSE THE LOGFILE
      CLOSE(UNIT=42)

110   CONTINUE

      WRITE(6,*) '... FORTRAN PROGRAM END'

      STOP
      END

C --------------------------------------------------------------------
C
C USER ROUTINES
C
C --------------------------------------------------------------------

      SUBROUTINE USER1(IACTION)

      IMPLICIT NONE

      INTEGER IFINAL, NOBS, NEVT
      DOUBLE PRECISION SUM, SUM2, VAR(100), VAR2(100)
      COMMON /DISSTAT/ SUM(100), SUM2(100), IFINAL, NOBS, NEVT

      INTEGER IACTION

      INTEGER I

C      WRITE(6,*) 'IACTION=', IACTION

      IF (IACTION .EQ. 2) THEN

         DO I = 1, NOBS

            VAR2(I)=(NEVT * SUM2(I) - SUM(I) ** 2) / (NEVT - 1)

            IF (VAR2(I) .GE. 0.) THEN
               VAR(I) = DSQRT(VAR2(I))
            ELSE
               WRITE(6,*) 'NEGATIVE VARIANCE = ',I,' ', VAR2(I)
               VAR(I) = 0.
            ENDIF

            WRITE(6,100) I, SUM(I), VAR(I)

100         FORMAT('SUM (', 1I3, ') =', 1E16.6, ' +- ', 1E16.6)
         ENDDO

      ELSEIF (IACTION .EQ. 3) THEN

         IFINAL = 0

      ELSEIF (IACTION .EQ. 4) THEN

         IFINAL = 1

         NOBS = 3
         NEVT = 0
         DO I = 1, NOBS
            SUM(I)  = 0.
            SUM2(I) = 0.
         ENDDO

      ENDIF

      RETURN
      END

C --------------------------------------------------------------------

      FUNCTION USER2(
     &            NR,
     &            NL,
     &            FVECT,
     &            NPARTONS,
     &            XB,
     &            Q2,
     &            XI,
     &            WEIGHT,
     &            IRPS,
     &            IALPHAS,
     &            IALPHAEM,
     &            LOGNF
     &         )

      IMPLICIT NONE

      DOUBLE PRECISION USER2

      INTEGER NR
      INTEGER NL
      DOUBLE PRECISION FVECT(0:3, -10:10, 1:*)
      INTEGER NPARTONS(1:*)
      DOUBLE PRECISION XB(1:*)
      DOUBLE PRECISION Q2(1:*)
      DOUBLE PRECISION XI(1:*)
      DOUBLE PRECISION WEIGHT(1:11, 1:*)
      INTEGER IRPS(1:*)
      INTEGER IALPHAS(1:*)
      INTEGER IALPHAEM(1:*)
      INTEGER LOGNF(1:*)

      INTEGER IFINAL, NOBS, NEVT
      DOUBLE PRECISION SUM, SUM2, VAR(100), VAR2(100)
      COMMON /DISSTAT/ SUM(100), SUM2(100), IFINAL, NOBS, NEVT

      INTEGER IFILLED(20)
      INTEGER NFILLED

      INTEGER I, J, K

      INTEGER ICLUST
      INTEGER NJETS
      DOUBLE PRECISION CUTSCALE2

      INTEGER IPDF_COLLECTION
      INTEGER IPDF_PARAMETRIZATION
      INTEGER IPDF_SET
      INTEGER IALPHAS_VARIANT
      INTEGER IALPHAS_ORDER
      DOUBLE PRECISION DALPHAS_LAMBDAQCD4
      INTEGER IALPHAEM_VARIANT

      DOUBLE PRECISION AOBS
      DOUBLE PRECISION DISASTER_CAO

      DOUBLE PRECISION OBS(30, 10)
      INTEGER IOBS

      DOUBLE PRECISION FFACTIN(4)
      DOUBLE PRECISION FFACTOUT(13)

      DOUBLE PRECISION AFACT

      DOUBLE PRECISION ALPHAS, ALPHAEM

      DOUBLE PRECISION TSUM

      INTEGER NFLOC

C      CALL PFPX_F()

C --- PARAMETERS FOR MRSD-', ALPHA_S = 0.1, ALPHA_EM = 1/137
      IPDF_COLLECTION      = 1
      IPDF_PARAMETRIZATION = 3
      IPDF_SET             = 31
      IALPHAS_VARIANT      = -1
      IALPHAS_ORDER        = 0
      DALPHAS_LAMBDAQCD4   = 0.2
      IALPHAEM_VARIANT     = 2

C --- PARAMETERS FOR MRSD-', ALPHA_S = 2-LOOP, ALPHA_EM = F.S.C.
C      IPDF_COLLECTION      = 1
C      IPDF_PARAMETRIZATION = 3
C      IPDF_SET             = 31
C      IALPHAS_VARIANT      = 1
C      IALPHAS_ORDER        = 2
C      DALPHAS_LAMBDAQCD4   = 0.2
C      IALPHAEM_VARIANT     = 1

      NFLOC = 5

C --- PRINT EVENT RECORDS AND CONTRIBUTIONS

      IF (0 .EQ. 0) THEN
         GOTO 200
      ENDIF

      IFILLED(1) = -8
      IFILLED(2) = -7
      IFILLED(3) = -5
      IFILLED(4) = -4
      IFILLED(5) = -1
      IFILLED(6) =  0
      IFILLED(7) =  1
      IFILLED(8) =  2
      NFILLED = 5

      WRITE(6,*) 'USER2 CALLED'

      DO I = 1, NR
         WRITE(6,110) I, NPARTONS(I), Q2(I), XI(I)
         DO J = 1, NFILLED + NPARTONS(I)
            WRITE(6,100) IFILLED(J),
     &                   FVECT(0, IFILLED(J), I),
     &                   FVECT(1, IFILLED(J), I),
     &                   FVECT(2, IFILLED(J), I),
     &                   FVECT(3, IFILLED(J), I)
         ENDDO
      ENDDO

      DO I = 1, NL
         WRITE(6,120) I,
     &                IRPS(I), IALPHAS(I), IALPHAEM(I), LOGNF(I)
         DO J = 1, 11
            WRITE(6,130) J, WEIGHT(J, I)
         ENDDO
      ENDDO

100   FORMAT(1I5, 4E16.6)
110   FORMAT(2I5, '-->', 2E16.6)
120   FORMAT(1I5, '++>', 4I5)
130   FORMAT(1I5, 1E16.6)

200   CONTINUE

C --------------------------------------------------------------------

C --- CALCULATE THE ADAPTATION OBSERVABLE
      AOBS = DISASTER_CAO(
     &          IPDF_COLLECTION,
     &          IPDF_PARAMETRIZATION,
     &          IPDF_SET,
     &          IALPHAS_VARIANT,
     &          IALPHAS_ORDER,
     &          DALPHAS_LAMBDAQCD4,
     &          IALPHAEM_VARIANT
     &       )

C      WRITE(6,*) 'AOBS=', AOBS

C --------------------------------------------------------------------

C --- EVALUATE OBSERVABLES IF THIS IS THE FINAL RUN

      IF (IFINAL .EQ. 1) THEN

C --------------------------------------------------------------------

C ------ SUM...

         NEVT = NEVT + 1

C ------ CLUSTER THE EVENTS
         DO I = 1, NR
            CUTSCALE2 = 0.02 * Q2(I) * (1.-XB(I)) / XB(I)
C --------- CAREFUL! ICLUST MODIFIES THE EVENT RECORD!
            NJETS = ICLUST(CUTSCALE2, FVECT, NPARTONS(I), I)
C            WRITE(6,*) 'NJETS=', NJETS
C --------- BIN FOR (1+1) JETS
            IF (NJETS .EQ. 1) THEN
               OBS(I, 1) = 1.
            ELSE
               OBS(I, 1) = 0.
            ENDIF
C --------- BIN FOR (2+1) JETS
            IF (NJETS .EQ. 2) THEN
               OBS(I, 2) = 1.
            ELSE
               OBS(I, 2) = 0.
            ENDIF
C --------- BIN FOR (3+1) JETS
            IF (NJETS .EQ. 3) THEN
               OBS(I, 3) = 1.
            ELSE
               OBS(I, 3) = 0.
            ENDIF
         ENDDO

C         IOBS = 2
C
C         DO I = 1, NR
C            OBS(I, IOBS) = 1.
C         ENDDO

C --------------------------------------------------------------------

C ------ CONVOLUTE OBSERVABLES, WEIGHTS AND FLAVOUR FACTORS

         DO IOBS = 1, NOBS

            TSUM = 0.

            DO I = 1, NL

C ------------ CALCULATE COUPLING CONSTANTS AND PRODUCTS
C              OF PARTON DENSITIES AND QUARK CHARGES

               FFACTIN(1) = XI(IRPS(I))
               FFACTIN(2) = DSQRT(Q2(IRPS(I)))
               FFACTIN(3) = DSQRT(Q2(IRPS(I)))
               FFACTIN(4) = DSQRT(Q2(IRPS(I)))

               CALL DISASTER_CFF(
     &                 IPDF_COLLECTION,
     &                 IPDF_PARAMETRIZATION,
     &                 IPDF_SET,
     &                 IALPHAS_VARIANT,
     &                 IALPHAS_ORDER,
     &                 DALPHAS_LAMBDAQCD4,
     &                 IALPHAEM_VARIANT,
     &                 NFLOC,
     &                 FFACTIN,
     &                 FFACTOUT
     &              )

               ALPHAS  = FFACTOUT(12)
               ALPHAEM = FFACTOUT(13)

C         WRITE(6,210) I, FFACTIN(1), FFACTIN(2)
C210      FORMAT(1I5, 2E16.6)

C ------------ ALL SCALES ASSUMED TO BE Q**2 --> LOGARITHMS VANISH

               IF (LOGNF(I) .EQ. 0) THEN
                  AFACT = 1.
               ELSE
                  AFACT = 0.
               ENDIF

               AFACT =   AFACT
     &                 * ALPHAS  ** IALPHAS(I)
     &                 * ALPHAEM ** IALPHAEM(I)
     &                 * OBS(IRPS(I), IOBS)

               DO J = 1, 11

                  TSUM = TSUM + AFACT * WEIGHT(J, I) * FFACTOUT(J)

               ENDDO

            ENDDO

C            WRITE(6,*) 'TSUM=', TSUM

            SUM(IOBS)  = SUM(IOBS)  + TSUM
            SUM2(IOBS) = SUM2(IOBS) + TSUM * TSUM

         ENDDO

C --------------------------------------------------------------------

      ENDIF

      USER2 = AOBS

      RETURN
      END

C --------------------------------------------------------------------

      SUBROUTINE USER3(AVERAGE, ERROR)

      IMPLICIT NONE

      DOUBLE PRECISION AVERAGE, ERROR

100   FORMAT('RESULT OF INTEGRATION =', 1E16.6, ' +- ', 1E16.6)

      WRITE(6, 100) AVERAGE, ERROR

      RETURN
      END

C --------------------------------------------------------------------
C
C --> PRINTING ROUTINE (CHAR* FROM C++)
C
C --------------------------------------------------------------------

      SUBROUTINE PRINT_F(STRING, LENGTH, UNIT)

      IMPLICIT NONE

      CHARACTER*1000 STRING
      INTEGER LENGTH, UNIT

      CHARACTER*7 FMTS

      INTEGER I1, I2, I3

      CHARACTER*1 NUMBER(0:9)
      DATA
     &   NUMBER(0) /'0'/,
     &   NUMBER(1) /'1'/,
     &   NUMBER(2) /'2'/,
     &   NUMBER(3) /'3'/,
     &   NUMBER(4) /'4'/,
     &   NUMBER(5) /'5'/,
     &   NUMBER(6) /'6'/,
     &   NUMBER(7) /'7'/,
     &   NUMBER(8) /'8'/,
     &   NUMBER(9) /'9'/

      I3 = LENGTH / 100
      I2 = (LENGTH - I3 * 100) / 10
      I1 = (LENGTH - I3 * 100 - I2 * 10)

      IF (LENGTH .EQ. 0) THEN

         WRITE(UNIT=UNIT, FMT=*)

      ELSE

         FMTS(1:1) = '('
         FMTS(2:2) = '1'
         FMTS(3:3) = 'A'
         FMTS(4:4) = NUMBER(I3)
         FMTS(5:5) = NUMBER(I2)
         FMTS(6:6) = NUMBER(I1)
         FMTS(7:7) = ')'

         WRITE(UNIT=UNIT, FMT=FMTS) STRING(1:LENGTH)

      ENDIF

      RETURN
      END

C --------------------------------------------------------------------
C
C --> PROVOKE A FLOATING POINT EXCEPTION
C
C --------------------------------------------------------------------

      SUBROUTINE PFPX_F()

      IMPLICIT NONE

      DOUBLE PRECISION D

      WRITE(6,*) 'PROVOKED FPE (FORTRAN)'
      D = DSQRT(-1.D0)
      WRITE(6,*) 'FPE (FORTRAN) VALUE = ', D

      RETURN
      END

C --------------------------------------------------------------------
