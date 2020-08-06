C --------------------------------------------------------------------
C
C --> INTERFACE TO VEGAS
C
C --------------------------------------------------------------------

      SUBROUTINE VEGAS1_F(DPAR, IPAR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION AVGT, ERRT, AVGI, ERRI

      COMMON /VGRES/
     &       AVGT,ERRT,AVGI,ERRI

      COMMON /VGB2/
     &       NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS,
     &       DI02GAU(50,10),DI(50,10)

      COMMON /VGASIO/
     &       NINP,NOUTP

      COMMON /BVEG3/
     &       ALPH,NDMX,MDS

      COMMON /CCALL/
     &       CALLS,ITMX

      EXTERNAL FFVEGAS1

      DIMENSION DPAR(0 : *)
      DIMENSION IPAR(0 : *)

      ACC      = DPAR(0)
      AVGTOUT  = DPAR(1)
      ERRTOUT  = DPAR(2)
      AVGIOUT  = DPAR(3)
      ERRIOUT  = DPAR(4)
      CALLSOUT = DPAR(5)
      ALPHIN   = DPAR(6)

      NDIM     = IPAR(0)
      NCALL    = IPAR(1)
      ITMXIN   = IPAR(2)
      NPRN     = IPAR(3)
      IGRAPH   = IPAR(4)
      IENTRY   = IPAR(5)
      NINPIN   = IPAR(6)
      NOUTPIN  = IPAR(7)
      ITMXOUT  = IPAR(8)
      MDSIN    = IPAR(9)

      NINP  = NINPIN
      NOUTP = NOUTPIN
      ALPH  = ALPHIN
      MDS   = MDSIN

      CALL VEGAS(FFVEGAS1,ACC,NDIM,NCALL,ITMXIN,NPRN,IGRAPH,IENTRY)

      AVGTOUT  = AVGT
      ERRTOUT  = ERRT
      AVGIOUT  = AVGI
      ERRIOUT  = ERRI
      CALLSOUT = CALLS
      ITMXOUT  = ITMX

      DPAR(0) = ACC
      DPAR(1) = AVGTOUT
      DPAR(2) = ERRTOUT
      DPAR(3) = AVGIOUT
      DPAR(4) = ERRIOUT
      DPAR(5) = CALLSOUT
      DPAR(6) = ALPHIN

      IPAR(0) = NDIM
      IPAR(1) = NCALL
      IPAR(2) = ITMXIN
      IPAR(3) = NPRN
      IPAR(4) = IGRAPH
      IPAR(5) = IENTRY
      IPAR(6) = NINPIN
      IPAR(7) = NOUTPIN
      IPAR(8) = ITMXOUT
      IPAR(9) = MDSIN

      RETURN
      END

C --------------------------------------------------------------------
C
C --> FUNCTION INTEGRATED BY VEGAS
C
C --------------------------------------------------------------------

      FUNCTION FFVEGAS1(X, WGT)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION X(10)

      COMMON /CCALL/
     &       CALLS,ITMX

C      WRITE(6,*) 'X, WGT = ', X(1), X(2), X(3), WGT

C --- CALL THE C ROUTINE
      FFVEGAS1 = FVEGAS1(X, WGT, CALLS, ITMX)

      RETURN
      END

C --------------------------------------------------------------------
