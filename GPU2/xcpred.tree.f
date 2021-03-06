      SUBROUTINE XCPRED(KCASE)
*
*
*       Prediction of global chain variables.
*       -------------------------------------
*
      INCLUDE 'common6.h'
        REAL*8  M,MASS,MC
        PARAMETER (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      SAVE
*
*
*       Check prediction indicator for perturbers & c.m. (new skip 05/16).
      IF (KCASE.EQ.0.OR.LISTC(1).EQ.0) GO TO 4
*
*       Check adding chain c.m. #ICH to perturber list for prediction.
      IF (KCASE.EQ.1) THEN
          NNB2 = LISTC(1) + 2
          LISTC(NNB2) = ICH
      ELSE
*      Note KCASE = 2 for chain c.m. prediction at new block-time.
          NNB2 = LISTC(1) + 1
      END IF
*
*       Predict X & XDOT of perturbers and (if needed) c.m. to order FDOT.
      DO 1 L = 2,NNB2
          J = LISTC(L)
          S = TIME - T0(J)
*       Accept positive interval outside range (no bad effects in DIFSY1).
*         S = MIN(S,STEP(J))
          S = MAX(STEP(J),0.0D0)   ! Note TIME may not be updated (05/16).
          S2 = 0.5*S
          S3 = S / 3.0D0
          X(1,J)= X0(1,J) +S*(X0DOT(1,J) +S2*(FI(1,J) +S3*FIDOT(1,J)))
          X(2,J)= X0(2,J) +S*(X0DOT(2,J) +S2*(FI(2,J) +S3*FIDOT(2,J)))
          X(3,J)= X0(3,J) +S*(X0DOT(3,J) +S2*(FI(3,J) +S3*FIDOT(3,J)))
          XDOT(1,J) = X0DOT(1,J) +S*(FI(1,J) +S2*FIDOT(1,J))
          XDOT(2,J) = X0DOT(2,J) +S*(FI(2,J) +S2*FIDOT(2,J))
          XDOT(3,J) = X0DOT(3,J) +S*(FI(3,J) +S2*FIDOT(3,J))
    1 CONTINUE
*
*       Obtain global coordinates & velocities from current chain & c.m.
   4  LK = 0
      DO 10 I = 1,NN
          DO 5 K = 1,3
              LK = LK + 1
              XC(K,I) = XCH(LK) + X(K,ICH)
              UC(K,I) = VCH(LK) + XDOT(K,ICH)
    5     CONTINUE
   10 CONTINUE
*
      RETURN
*
      END
