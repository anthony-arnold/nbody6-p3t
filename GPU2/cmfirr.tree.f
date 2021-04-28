      SUBROUTINE CMFIRR(I,I1,FIRR,FD)
*
*
*       Irregular force corrections for active KS.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/PREDICT/ TPRED(NMAX)
      REAL*8  FIRR(3),FD(3),DX(3),DV(3),FP(6),FPD(6),FCM(3),FCMD(3)
      REAL*8  CMX(3),CMV(3),XK(6),VK(6)
*
*
*       Set perturber number and KS indices.
      NP = LIST(1,I1)
      I2 = I1 + 1
      J1 = -1
      BODYIN = 1.0/BODY(I)
*
*       Form irregular force components for perturbed KS pair.
      DO 50 LL = 2,NP+1
          J = LIST(LL,I1)
*       Obtain coordinates & velocity by prediction or copying.
          IPRED = 1
          IF (TPRED(J).EQ.TIME) IPRED = 0
          IF (J.LE.N) THEN
              CALL JPRED2(J,IPRED,CMX,CMV)
              K = J
          ELSE
              JPAIR = J - N
              J1 = 2*JPAIR - 1
              IF (LIST(1,J1).EQ.0) THEN
                  CALL JPRED2(J,IPRED,CMX,CMV)
                  K = J
              ELSE
                  CALL KSRES3(JPAIR,J1,J2,IPRED,CMX,CMV,XK,VK)
                  K = J1
              END IF
          END IF
*
*       Copy c.m. values for case of single or unperturbed KS.
          IF (K.EQ.J) THEN
              DO 10 L = 1,3
                  XK(L) = CMX(L)
                  VK(L) = CMV(L)
   10         CONTINUE
          END IF
*
*       Obtain individual c.m. force with single particle approximation.
          DO L = 1,3
             FCM(L) = 0
             FCMD(L) = 0
          END DO
          CALL CUTOFF(X(1,I),XDOT(1,I),CMX,CMV,FCM,FCMD,BODY(J),RS(I))
*
*     Evaluate perturbation on first component due to body #K.
 20       DO L = 1,3
             FP(L) = 0
             FP(L+3) = 0
             FPD(L) = 0
             FPD(L+3) = 0
          END DO
          CALL CUTOFF(X(1,I1),XDOT(1,I1),XK,VK,
     &                FP, FPD, BODY(K), RS(I))
          CALL CUTOFF(X(1,I2),XDOT(1,I2),XK,VK,
     &                FP(4), FPD(4), BODY(K), RS(I))
*
*     Accumulate individual correction terms.
c          WRITE(*,*) 'CHECK', (FIRR(K),K=1,3), (FD(K),K=1,3),
c     &         (FCM(K),K=1,3), (FCMD(K),K=1,3),
c     &         (FP(K),K=1,6), (FPD(K),K=1,6)
          DO 40 L = 1,3
              FIRR(L) = FIRR(L) + ((BODY(I1)*FP(L) +
     &                              BODY(I2)*FP(L+3))*BODYIN - FCM(L))
              FD(L) = FD(L) + ((BODY(I1)*FPD(L) +
     &                          BODY(I2)*FPD(L+3))*BODYIN - FCMD(L))
   40     CONTINUE
*
*       Include possible second KS component without c.m. contribution.
          IF (K.EQ.J1) THEN
              DO 45 L = 1,3
                  FCM(L) = 0.0
                  FCMD(L) = 0.0
                  XK(L) = XK(L+3)
                  VK(L) = VK(L+3)
   45         CONTINUE
              K = K + 1
              GO TO 20
          END IF
   50 CONTINUE
*
*       Check force correction due to regularized chain.
      IF (NCH.GT.0) THEN
          CALL KCPERT(I,I1,FIRR,FD)
      END IF
*
      RETURN
*
      END
