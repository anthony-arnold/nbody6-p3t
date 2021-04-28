      SUBROUTINE KCPERT(II,I1,FIRR,FD)
*
*
*       Differential force correction on active KS from chain.
*       ------------------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  FIRR(3),FD(3),DX(3),DV(3),FP(3),FPD(3),FCM(3),FCMD(3)
*
*
*       See whether perturber list contains chain c.m. body #ICH.
      J = 0
      NNB1 = LIST(1,I1) + 1
      DO 1 L = 2,NNB1
          JJ = LIST(L,I1)
          IF (JJ.EQ.ICH) J = JJ
    1 CONTINUE
*       Exit if chain c.m. not identified.
      IF (J.EQ.0) GO TO 30
*
*       Obtain chain c.m. force on each KS component for diff correction.
      BODYIN = 1.0/BODY(II)
      I = I1
 5    DO K=1,3
         FCM(K) = 0.0D0
         FCMD(K) = 0.0D0
      END DO
      CALL CUTOFF(X(1,I),XDOT(1,I),X(1,J),XDOT(1,J),
     &            FCM, FCMD, BODY(J), RS(I))
*
      DO 10 K = 1,3
          FP(K) = 0.0
          FPD(K) = 0.0
   10 CONTINUE
*
*       Form contributions from each chain member for each component.
      DO 20 JJ = 1,NCH
         CALL CUTOFF(X(1,I),XDOT(1,I),XC(1,JJ),UC(1,JJ),
     &               FP,FPD,BODYC(JJ),RS(I))
   20 CONTINUE
*
*       Add differential corrections (chain contribution minus c.m. term).
      DO 25 K = 1,3
          FIRR(K) = FIRR(K) + BODY(I)*(FP(K) - FCM(K))*BODYIN
          FD(K) = FD(K) + BODY(I)*(FPD(K) - FCMD(K))*BODYIN
   25 CONTINUE
*
*       Treat the second component in a similar way.
      IF (I.EQ.I1) THEN
          I = I1 + 1
          GO TO 5
      END IF
*
   30 RETURN
*
      END
