      REAL*8 FUNCTION TSTEP(F,FDOT,F2DOT,F3DOT,ETA,BODY,RS)
*
*
*       General time-step criterion.
*       ----------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  F(3),FDOT(3),F2DOT(3),F3DOT(3)
*
*       Obtain new integration step using composite expression.
*       STEP = (ETA*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
*
      F2 = F(1)**2 + F(2)**2 + F(3)**2
      FDOT2 = FDOT(1)**2 + FDOT(2)**2 + FDOT(3)**2
      F2DOT2 = F2DOT(1)**2 + F2DOT(2)**2 + F2DOT(3)**2
      F3DOT2 = F3DOT(1)**2 + F3DOT(2)**2 + F3DOT(3)**2
*
      TSTEP = SQRT(F2*F3DOT2)+F2DOT2
      IF (TSTEP.NE.0.0D0) THEN
         TSTEP =  (SQRT(F2*F2DOT2)+FDOT2)/TSTEP
         TSTEP = SQRT(ETA*TSTEP)
      END IF
*
      RETURN
*
      END
