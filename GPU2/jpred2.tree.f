      SUBROUTINE JPRED2(I,IPRED,CMX,CMV)
*
*
*       Neighbour prediction of single particle.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  CMX(3),CMV(3)
*
*
*       Choose between copy and standard prediction.
      IF (IPRED.EQ.0) THEN
          DO 10 K = 1,3
              CMX(K) = X(K,I)
              CMV(K) = XDOT(K,I)
   10     CONTINUE
      ELSE
          S = MAX(TIME - T0(I), 0.0D0)
          S2 = 0.5*S
          S3 = S / 3.0D0
          CMX(1) = X0(1,I)+ S*(X0DOT(1,I)+ S2*(FI(1,I)+ S3*FIDOT(1,I)))
          CMX(2) = X0(2,I)+ S*(X0DOT(2,I)+ S2*(FI(2,I)+ S3*FIDOT(2,I)))
          CMX(3) = X0(3,I)+ S*(X0DOT(3,I)+ S2*(FI(3,I)+ S3*FIDOT(3,I)))
          CMV(1) = X0DOT(1,I) + S*(FI(1,I) + S2*FIDOT(1,I))
          CMV(2) = X0DOT(2,I) + S*(FI(2,I) + S2*FIDOT(2,I))
          CMV(3) = X0DOT(3,I) + S*(FI(3,I) + S2*FIDOT(3,I))
      END IF
*
      RETURN
*
      END
