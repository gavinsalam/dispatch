C---------------------------------------------------------------------
C
C ---> POLYLOGARITHMS
C
C ---> DIRK GRAUDENZ
C
C      DERIVED FROM DPOLY.F
C
C ---> 02\10\1991
C ---> 17\03\1995
C
C ---------------------------------------------------------------------
C
C ===> MAP COMPLEX NUMBER INTO |Z| <= 1, RE(Z) <= 1/2
C      WITH THE EXCEPTION OF Z=1.D0+0.D0*I
C
C ---------------------------------------------------------------------
C
      SUBROUTINE CMAP(Z,ZHAT,IMAP)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 Z,ZHAT
C
      AZ=CDABS(Z)
      A1MZ=CDABS(1.D0-Z)
      DREZ=DREAL(Z)
      DIMZ=DIMAG(Z)
C
      IF (DREZ .EQ. 1.D0 .AND. DIMZ .EQ. 0.D0) THEN
         IMAP=0
         ZHAT=Z
      ELSE
         IF (DREZ .LE. .5D0) THEN
            IF (AZ .LE. 1.D0) THEN
               IMAP=1
               ZHAT=Z
            ELSE
               IMAP=2
               ZHAT=1.D0/Z
            ENDIF
         ELSE
            IF (A1MZ .LE. 1.D0) THEN
               IMAP=3
               ZHAT=1.D0-Z
            ELSE
               IMAP=4
               ZHAT=1.D0/(1.D0-Z)
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> SIGN FUNCTION
C      DNEWSIGN(X,Y)=|X|*SGN(Y)
C
C ---------------------------------------------------------------------
C
      FUNCTION DNEWSIGN(X,Y)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      IF (Y .LT. 0.D0) THEN
         DNEWSIGN=-DABS(X)
      ELSEIF (Y .EQ. 0.D0) THEN
         DNEWSIGN=0.D0
      ELSE
         DNEWSIGN=DABS(X)
      ENDIF
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> INITIALIZE CONSTANTS
C
C ---------------------------------------------------------------------
C
      SUBROUTINE DMATHINIT()
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      DATA ICOLD/0/
C
      ICOLD=1
C
      DBERN(1)= 1.D0/6.D0
      DBERN(2)=-1.D0/30.D0
      DBERN(3)= 1.D0/42.D0
      DBERN(4)=-1.D0/30.D0
      DBERN(5)= 5.D0/66.D0
      DBERN(6)=-691.D0/2730.D0
      DBERN(7)= 7.D0/6.D0
      DBERN(8)=-3617.D0/510.D0
      DBERN(9)= 43867.D0/798.D0

      DO K=1,9
         DKCONST(K)=1.D0/DFLOAT(2*K*(2*K+1))
      ENDDO

      DBERNOULLI( 0)=1.0
      DBERNOULLI( 1)=-.5
      DBERNOULLI( 2)=DBERN(1)
      DBERNOULLI( 3)=0.
      DBERNOULLI( 4)=DBERN(2)
      DBERNOULLI( 5)=0.
      DBERNOULLI( 6)=DBERN(3)
      DBERNOULLI( 7)=0.
      DBERNOULLI( 8)=DBERN(4)
      DBERNOULLI( 9)=0.
      DBERNOULLI(10)=DBERN(5)
      DBERNOULLI(11)=0.
      DBERNOULLI(12)=DBERN(6)
      DBERNOULLI(13)=0.
      DBERNOULLI(14)=DBERN(7)
      DBERNOULLI(15)=0.
      DBERNOULLI(16)=DBERN(8)
      DBERNOULLI(17)=0.
      DBERNOULLI(18)=DBERN(9)
      DBERNOULLI(19)=0.
      DBERNOULLI(20)=-174611./330.
      DBERNOULLI(21)=0.
C
      IPMAX=21
C
      DFACT(0)=1.
      DO 10, I=1,IPMAX+2
         DFACT(I)=I*DFACT(I-1)
 10   CONTINUE
C
      DO 20, I=0,IPMAX
         DS12C(I)=0.5*(DFLOAT(I)+1.)*DBERNOULLI(I)/DFACT(I+2)
         DSUM=0.
         DO 30, J=0,I
            DSUM=DSUM
     &          +DBERNOULLI(J)*DBERNOULLI(I-J)/DFACT(J+1)
     &           *DFACT(I)/DFACT(I-J)/DFACT(I+1)
 30      CONTINUE
         DLI3C(I)=DSUM
 20   CONTINUE
C
      DPI=3.1415926535897932384D0
      DZETA2=DPI*DPI/6.D0
C
C --- ???
C
      DZETA3=1.20205690315959D0
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> ERROR MESSAGE
C
C ---------------------------------------------------------------------
C
      SUBROUTINE DMATHERR(IERR)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMMON /DMATHDIS/ NERR
      DATA NERR /0/
C
      NERRMAX=10
      IF (NERR .LT. NERRMAX) THEN
         NERR=NERR+1
         WRITE(6,*) '-----------------------------------'
         WRITE(6,*) '! DMATH ERROR #', IERR
         WRITE(6,*) '-----------------------------------'
         IF (NERR .EQ. NERRMAX) THEN
            WRITE(6,*) '*** --> ... LAST WARNING! <-- ***'
         ENDIF
      ENDIF
C
      IF (IERR .LT. 0) THEN
         STOP
      ENDIF
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> COMPLEX LOGARITHM WITH A CUT ALONG (-INF,0]
C
C ---------------------------------------------------------------------
C
      FUNCTION CNEWLOG(Z)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 Z,CNEWLOG
C
      IF (ICOLD .EQ. 0) THEN
         CALL DMATHINIT()
      ENDIF
C
      X=DREAL(Z)
      Y=DIMAG(Z)
C
C --- CHECK THE CUT OF THE LOG
C
      IF (Y .EQ. 0.D0 .AND. X .LE. 0.D0) THEN
         CALL DMATHERR( 1)
      ENDIF
C
      R=CDABS(Z)
C
      IF (X .LT. 0.D0) THEN
         IF (Y .NE. 0.D0) THEN
            PHI=DATAN(Y/X)+DNEWSIGN(DPI,Y)
         ELSE
            PHI=DPI
         ENDIF
      ELSEIF (X .EQ. 0.D0) THEN
         PHI=DNEWSIGN(DPI/2.D0,Y)
      ELSE
         PHI=DATAN(Y/X)
      ENDIF
C
      CNEWLOG=DCMPLX(DLOG(R),PHI)
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> ELEMENTARY SPENCE FUNCTION FOR |Z|<=1, RE(Z)<=0.5
C
C :_MOD_: CORRECT ONLY FOR FABS(Z) <~ 1.E-16; == 0 BELOW! (LOGARITHM)
C --> EXPLICIT EXPANSION!
C ... ALSO CHECK THE BOUND KMAX...
C
C ---------------------------------------------------------------------
C
      FUNCTION CELSPENCE(Z)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 Z,CELSPENCE,ZW,ZU,ZSUM,CNEWLOG
C
      IF (ICOLD .EQ. 0) THEN
         CALL DMATHINIT()
      ENDIF
C
C --- CHECK THE RANGE OF Z
C
      IF (CDABS(Z) .GT. 1.D0 .OR. DREAL(Z) .GT. .5D0) THEN
         CALL DMATHERR( 2)
      ENDIF
C
C      WRITE(6,*) Z
      ZW=-CNEWLOG(1.D0-Z)
      ZSUM=ZW-.25D0*ZW*ZW
      ZU=ZW
C
C --- HERE: SUM FROM 9 TO 1 TO OBTAIN HIGHER ACCURACY ???
      IF (DREAL(ZW) .LT. 1.E-10 .AND. DIMAG(ZW) .LT. 1.E-10) THEN
         KMAX = 5
CC         KMAX = 9
      ELSE
         KMAX = 9
      ENDIF

      DO 1, K=1,KMAX
CC         WRITE(6,*) ZU,ZW,K
CCCC         ZU=ZU*ZW*ZW/DFLOAT(2*K*(2*K+1))
         ZU=ZU*ZW*ZW*DKCONST(K)
         ZSUM=ZSUM+ZU*DBERN(K)
1     CONTINUE
C
      CELSPENCE=ZSUM
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> COMPLEX SPENCE FUNCTION
C
C ---------------------------------------------------------------------
C
      FUNCTION CS2(Z)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 Z,CS2,ZHAT,ZRESULT,CNEWLOG,CELSPENCE
C
      IF (ICOLD .EQ. 0) THEN
         CALL DMATHINIT()
      ENDIF
C
C --- CHECK THE CUT OF THE SPENCE FUNCTION
C
      IF (DIMAG(Z) .EQ. 0.D0 .AND. DREAL(Z) .GT. 1.D0) THEN
         CALL DMATHERR( 3)
      ENDIF
C
      CALL CMAP(Z,ZHAT,IMAP)
C
      IF (IMAP .EQ. 0) THEN
         ZRESULT=DZETA2
      ELSEIF (IMAP .EQ. 1) THEN
         ZRESULT=CELSPENCE(ZHAT)
      ELSEIF (IMAP .EQ. 2) THEN
         ZRESULT=-CELSPENCE(ZHAT)-.5D0*CNEWLOG(-Z)**2-DZETA2
      ELSEIF (IMAP .EQ. 3) THEN
         ZRESULT=-CELSPENCE(ZHAT)-CNEWLOG(Z)*CNEWLOG(1.D0-Z)+DZETA2
      ELSEIF (IMAP .EQ. 4) THEN
         ZRESULT= CELSPENCE(ZHAT)-CNEWLOG(1.D0-Z)*CNEWLOG(-Z)
     &           +.5D0*CNEWLOG(1.D0-Z)**2-DZETA2
      ENDIF
C
      CS2=ZRESULT
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> REAL SPENCE FUNCTION FOR ARGUMENTS NOT ON THE CUT
C
C ---------------------------------------------------------------------
C
      FUNCTION SPENCE(X)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 CS2
C
      SPENCE=DREAL(CS2(DCMPLX(X,0.D0)))
C
C      WRITE(6,*) 'FORTRAN: ARG, VALUE=', X, SPENCE
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> ELEMENTARY S1,2(Z) FOR |Z|<=1, RE(Z)<=0.5
C
C ---------------------------------------------------------------------
C
      FUNCTION CELS12(Z)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 Z,CELS12,ZW,ZSUM,CNEWLOG
C
      IF (ICOLD .EQ. 0) THEN
         CALL DMATHINIT()
      ENDIF
C
C --- CHECK THE RANGE OF Z
C
      IF (CDABS(Z) .GT. 1.D0 .OR. DREAL(Z) .GT. .5D0) THEN
         CALL DMATHERR( 2)
      ENDIF
C
      ZW=-CNEWLOG(1.D0-Z)
      ZSUM=0.
      DO 10,I=IPMAX,0,-1
         ZSUM=ZSUM+DS12C(I)*ZW**(I+2)
 10   CONTINUE
C
      CELS12=ZSUM
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> COMPLEX S1,2
C
C ---------------------------------------------------------------------
C
      FUNCTION CS12(Z)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 Z,CS12,ZHAT,ZRESULT,CNEWLOG,CS2,
     &           CELS12,CELS3
C
      IF (ICOLD .EQ. 0) THEN
         CALL DMATHINIT()
      ENDIF
C
C --- CHECK THE CUT OF THE SPENCE FUNCTION
C
      IF (DIMAG(Z) .EQ. 0.D0 .AND. DREAL(Z) .GT. 1.D0) THEN
         CALL DMATHERR( 3)
      ENDIF
C
      CALL CMAP(Z,ZHAT,IMAP)
CC      WRITE(6,*) 'IMAP=',IMAP
CC      WRITE(6,*) 'Z   =',Z
CC      WRITE(6,*) 'ZHAT=',ZHAT
C
      IF (IMAP .EQ. 0) THEN
         ZRESULT=DZETA3
      ELSEIF (IMAP .EQ. 1) THEN
         ZRESULT=CELS12(ZHAT)
      ELSEIF (IMAP .EQ. 2) THEN
         ZRESULT=-CELS12(ZHAT)+CELS3(ZHAT)
     &          -1./3.*CNEWLOG(-Z)**3-DZETA2*CNEWLOG(-Z)
     &          -CNEWLOG(-Z)*CS2(Z)+DZETA3
      ELSEIF (IMAP .EQ. 3) THEN
         ZRESULT=-CELS3(ZHAT)-CNEWLOG(ZHAT)*CS2(Z)
     &          -0.5*CNEWLOG(Z)*CNEWLOG(ZHAT)**2
     &          +DZETA2*CNEWLOG(ZHAT)+DZETA3
      ELSEIF (IMAP .EQ. 4) THEN
         ZRESULT=-CELS3(ZHAT)-0.5*CNEWLOG(1.-Z)**2*CNEWLOG(-Z)
     &          +1./6.*CNEWLOG(1.-Z)**3-CNEWLOG(1.-Z)*CS2(Z)
     &          +DZETA3-DZETA2*CNEWLOG(1.-Z)
      ENDIF
C
      CS12=ZRESULT
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> REAL S1,2 FOR ARGUMENTS NOT ON THE CUT
C
C ---------------------------------------------------------------------
C
      FUNCTION SP12(X)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 CS12
C
      SP12=DREAL(CS12(DCMPLX(X,0.D0)))
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> ELEMENTARY LI3(Z) FOR |Z|<=1, RE(Z)<=0.5
C
C ---------------------------------------------------------------------
C
      FUNCTION CELS3(Z)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 Z,CELS3,ZW,ZSUM,CNEWLOG
C
      IF (ICOLD .EQ. 0) THEN
         CALL DMATHINIT()
      ENDIF
C
C --- CHECK THE RANGE OF Z
C
      IF (CDABS(Z) .GT. 1.D0 .OR. DREAL(Z) .GT. .5D0) THEN
         CALL DMATHERR( 2)
      ENDIF
C
      ZW=-CNEWLOG(1.D0-Z)
      ZSUM=0.
      DO 10,I=IPMAX,0,-1
         ZSUM=ZSUM+DLI3C(I)*ZW**(I+1)
 10   CONTINUE
C
      CELS3=ZSUM
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> COMPLEX LI3
C
C ---------------------------------------------------------------------
C
      FUNCTION CS3(Z)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 Z,CS3,ZHAT,ZRESULT,CNEWLOG,CS2,
     &           CELS12,CELS3
C
      IF (ICOLD .EQ. 0) THEN
         CALL DMATHINIT()
      ENDIF
C
C --- CHECK THE CUT OF THE SPENCE FUNCTION
C
      IF (DIMAG(Z) .EQ. 0.D0 .AND. DREAL(Z) .GT. 1.D0) THEN
         CALL DMATHERR( 3)
      ENDIF
C
      CALL CMAP(Z,ZHAT,IMAP)
CC      WRITE(6,*) 'IMAP=',IMAP
CC      WRITE(6,*) 'Z   =',Z
CC      WRITE(6,*) 'ZHAT=',ZHAT
C
      IF (IMAP .EQ. 0) THEN
         ZRESULT=DZETA3
      ELSEIF (IMAP .EQ. 1) THEN
         ZRESULT=CELS3(ZHAT)
      ELSEIF (IMAP .EQ. 2) THEN
         ZRESULT=-(-CELS3(1./Z)+1./6.*CNEWLOG(-Z)**3+DZETA2*CNEWLOG(-Z))
      ELSEIF (IMAP .EQ. 3) THEN
         ZRESULT=-CELS12(1.-Z)+CNEWLOG(Z)*CS2(Z)
     &          +0.5*CNEWLOG(1.-Z)*CNEWLOG(Z)**2+DZETA3
      ELSEIF (IMAP .EQ. 4) THEN
         ZRESULT=CELS12(1./(1.-Z))-CELS3(1./(1.-Z))
     &          +CNEWLOG(1.-1./(1.-Z))*CS2(1./(1.-Z))
     &          +1./6.*CNEWLOG(1./(1.-Z))**3+DZETA2*CNEWLOG(-1./Z)
     &          +0.5*CNEWLOG(1./(1.-Z))
     &           *CNEWLOG(1.-1./(1.-Z))*CNEWLOG(-Z)
      ENDIF
C
      CS3=ZRESULT
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
C ===> REAL LI3 FOR ARGUMENTS NOT ON THE CUT
C
C ---------------------------------------------------------------------
C
      FUNCTION SP3(X)
C
C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C
C
      COMPLEX*16 CS3
C
      SP3=DREAL(CS3(DCMPLX(X,0.D0)))
C
      RETURN
      END
C

C ---------------------------------------------------------------------
C
C ===> DEI(X)=E_1(X)+GAMMA+LN(X)
C      CF. ABRAMOWITZ & STEGUN
C
C ---------------------------------------------------------------------

      FUNCTION DEI(X)

C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C

      IF (X .LT. 0.D0) THEN
         WRITE(6,*) 'TRIED TO EVALUATE DEI(X) AT',X
         STOP
      ELSEIF (X .LE. 1) THEN
         DEI= .57721566490153286061D0
     &      +(-.57721566 + X * (
     &           .99999193 + X * (
     &            -.24991055 + X * (
     &               .05519968 + X * (
     &                -.00976004 + X * (
     &                   .00107857
     &                                 )
     &                               )
     &                             )
     &                           )
     &                         )
     &       )
      ELSE
         DEI= .57721566490153286061D0 + DLOG(X)
     &       +DEXP(-X)/X
     &      *( .2677737343 + X * (
     &           8.6347608925 + X * (
     &             18.0590169730 + X * (
     &               8.5733287401 + X
     &                                 )
     &                              )
     &                           )
     &       )
     &      /( 3.9584969228 + X * (
     &           21.0996530827 + X * (
     &             25.6329561486 + X * (
     &               9.5733223454 + X
     &                                 )
     &                               )
     &                            )
     &       )

      ENDIF

      RETURN
      END

C ---------------------------------------------------------------------
C
C ===> DEI(X)=E_1(X)+GAMMA+LN(X)
C      BY DIRECT SUMMATION OF N TERMS
C
C ---------------------------------------------------------------------

      FUNCTION DEISUM(X,NMAX)

C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C

      SUM=0.D0

      DO 20, N=NMAX,1,-1

      FACT=1.
      DO 10, I=1,N
         FACT=FACT*(X/DFLOAT(I))
10    CONTINUE

      SUM = SUM - (-1)**N/DFLOAT(N)*FACT

20    CONTINUE

      DEISUM=SUM

      RETURN
      END

C ---------------------------------------------------------------------
C
C ===> GAMMA FUNCTION
C      CF. ABRAMOWITZ & STEGUN
C
C ---------------------------------------------------------------------

      FUNCTION DGGAMMA(X)

C
C --- COMMON BLOCKS FOR DPOLY.F
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /CDMATH/
     &       DBERN(9),DPI,DZETA2,DZETA3,
     &       DBERNOULLI(0:21),
     &       DFACT(0:23),
     &       DLI3C(0:21),
     &       DS12C(0:21),
     &       DKCONST(9),
     &       ICOLD,IPMAX
C
C

      IF (X .GT. 10000. .OR. X .LT. -10000.) THEN
         WRITE(6,*) 'ARE YOU KIDDING???'
         STOP
      ENDIF

      FACT=1.D0

10    CONTINUE
         IF (X .LT. 1.9999D0) GOTO 20
         FACT=FACT*(X-1.D0)
         X=X-1.D0
         GOTO 10
20    CONTINUE

30    CONTINUE
         IF (X .GT. 0.9998D0) GOTO 40
         FACT=FACT/X
         X=X+1.D0
         GOTO 30
40    CONTINUE

C --- NOW X IS ROUGHLY IN BETWEEN OF 1 AND 2...

      Y=X-1.D0

      VALUE=1.+ Y *
     &      (-.577191652 + Y *
     &      ( .988205891 + Y *
     &      (-.897056937 + Y *
     &      ( .918206857 + Y *
     &      (-.756704078 + Y *
     &      ( .482199394 + Y *
     &      (-.193527818 + Y *
     &      ( .035868343
     &      ))))))))

      DGGAMMA=VALUE*FACT

CC      WRITE(6,*) X,Y,FACT,VALUE

      RETURN
      END

C ---------------------------------------------------------------------
