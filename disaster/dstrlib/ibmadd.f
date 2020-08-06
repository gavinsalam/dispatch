c --------------------------------------------------------------------
c
c source: ibmadd.f (by Th. Hadig)
c
c --------------------------------------------------------------------

      SUBROUTINE PDFSET(PARM, VAL)

      CHARACTER*20 PARM(20)
      DOUBLE PRECISION VAL(20)

      CALL PDFSET_(PARM, VAL)
      RETURN
      END


      SUBROUTINE PFTOPDG(X, SCALE, PDF)

      DOUBLE PRECISION X, SCALE
      DOUBLE PRECISION PDF(-6:6)

      CALL PFTOPDG_(X, SCALE, PDF)
      RETURN
      END

      FUNCTION ALPHAS2 (SCALE)

      DOUBLE PRECISION SCALE, ALPHAS2

      ALPHAS2 = ALPHAS2_(SCALE)

      RETURN
      END

c --------------------------------------------------------------------
