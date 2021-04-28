      SUBROUTINE CMFIRR2(I,NNB,KLIST,XI,XID,FIRR,FD)
*
*
*       Irregular force correction from J > N (singles or unpert KS).
*       -------------------------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/PREDICT/ TPRED(NMAX)
      REAL*8  XI(3),XID(3),FIRR(3),FD(3),DX(3),DV(3),
     &        FP(6),FPD(6),FCM(3),FCMD(3),CMX(3),CMV(3),XK(6),VK(6)
      INTEGER KLIST(LMAX)
*
*
*       Set neighbour number & list index of the last single particle.
      NNB1 = NNB + 1
      NNB2 = NNB1
   25 IF (KLIST(NNB2).LE.N) GO TO 30
      NNB2 = NNB2 - 1
      IF (NNB2.GT.1) GO TO 25
*       Include special case of only c.m. neighbours.
      GO TO 40
*
*       Exit if no binaries on neighbour list.
   30 IF (NNB2.EQ.NNB1) GO TO 60
*       Perform differential correction for regularized c.m. neighbours.
 40   DO 50 LL = NNB2+1,NNB1
          J = KLIST(LL)
*       See whether to sum over binary components.
          JPAIR = J - N
          J1 = 2*JPAIR - 1
*       Skip unperturbed binary (treated as single particle).
          IF (LIST(1,J1).EQ.0) THEN
              GO TO 50
          END IF
*
*       Obtain coordinates & velocity by prediction or copying.
          IPRED = 1
          IF (TPRED(J).EQ.TIME) IPRED = 0
          CALL KSRES3(JPAIR,J1,J2,IPRED,CMX,CMV,XK,VK)
*
*     Obtain individual c.m. force with single particle approximation.
          DO L = 1,3
             FCM(L) = 0
             FCMD(L) = 0
          END DO
          CALL CUTOFF(X(1,I),XDOT(1,I),CMX,CMV,FCM,FCMD,BODY(J),RS(I))
*
*     Evaluate perturbation on first component of body #J.
          DO L = 1,3
             FP(L) = 0
             FP(L+3) = 0
             FPD(L) = 0
             FPD(L+3) = 0
          END DO
          CALL CUTOFF(X(1,I),XDOT(1,I),XK,VK,
     &                FP, FPD, BODY(J1), RS(I))
          CALL CUTOFF(X(1,I),XDOT(1,I),XK(4),VK(4),
     &                FP(4), FPD(4), BODY(J2), RS(I))
*
*     Accumulate individual correction terms together.
          DO 49 L = 1,3
              FIRR(L) = FIRR(L) + (FP(L) + FP(L+3) - FCM(L))
              FD(L) = FD(L) + (FPD(L) + FPD(L+3) - FCMD(L))
   49     CONTINUE
   50 CONTINUE
*
   60 RETURN
*
      END
