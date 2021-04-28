      SUBROUTINE JPRED(I)
*
*
*       Neighbour prediction of single particle.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/PREDICT/ TPRED(NMAX)
*
*
*       Skip prediction if done at current time.
      IF (TPRED(I).EQ.TIME) GO TO 10
*
*       Adopt low order prediction for standard single particle.
      S = TIME - T0(I)
      S2 = 0.5*S
      S3 = S / 3.0D0
      X(1,I) = X0(1,I) + S*(X0DOT(1,I) + S2*(FI(1,I) + S3*FIDOT(1,I)))
      X(2,I) = X0(2,I) + S*(X0DOT(2,I) + S2*(FI(2,I) + S3*FIDOT(2,I)))
      X(3,I) = X0(3,I) + S*(X0DOT(3,I) + S2*(FI(3,I) + S3*FIDOT(3,I)))
      XDOT(1,I) = X0DOT(1,I) + S*(FI(1,I) + S2*FIDOT(1,I))
      XDOT(2,I) = X0DOT(2,I) + S*(FI(2,I) + S2*FIDOT(2,I))
      XDOT(3,I) = X0DOT(3,I) + S*(FI(3,I) + S2*FIDOT(3,I))
*
*       Resolve the components of any perturbed pair.
      IF (I.GT.N) THEN
          JPAIR = I - N
          IF (LIST(1,2*JPAIR-1).GT.0) THEN
!$omp critical
              ZZ = 1.0
*       Distinguish between low and high-order prediction of U & UDOT.
              IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
              CALL KSRES2(JPAIR,J1,J2,ZZ)
!$omp end critical
          END IF
      END IF
*
      TPRED(I) = TIME
*
   10 RETURN
*
      END
