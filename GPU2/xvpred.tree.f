      SUBROUTINE XVPRED(I1,NNB)
*
*
*       Prediction of coordinates & velocities.
*       ---------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Perform full N prediction to highest order (defined by NNB >= N).
      IF (NNB.GE.N) GO TO 20
*       Use high order for special singles (NNB < 0) or c.m. bodies (I1 > N).
      IF (NNB.LT.0.OR.(NNB.EQ.0.AND.I1.GT.N)) GO TO 20
*
*       Adopt low order for standard single particles (I1 <= N & NNB = 0).
      IF (NNB.EQ.0) THEN
          I = I1
          S = TIME - T0(I)
          IF (S.EQ.0.0D0) GO TO 40
          S2 = 0.5*S
          S3 = S / 3.0D0
          X(1,I) = X0(1,I)+ S*(X0DOT(1,I)+ S2*(FI(1,I)+ S3*FIDOT(1,I)))
          X(2,I) = X0(2,I)+ S*(X0DOT(2,I)+ S2*(FI(2,I)+ S3*FIDOT(2,I)))
          X(3,I) = X0(3,I)+ S*(X0DOT(3,I)+ S2*(FI(3,I)+ S3*FIDOT(3,I)))
          XDOT(1,I) = X0DOT(1,I)+ S*(FI(1,I)+ S2*FIDOT(1,I))
          XDOT(2,I) = X0DOT(2,I)+ S*(FI(2,I)+ S2*FIDOT(2,I))
          XDOT(3,I) = X0DOT(3,I)+ S*(FI(3,I)+ S2*FIDOT(3,I))
          GO TO 40
      END IF
*
*       Predict all neighbours to order FDOT.
      NNB1 = NNB+ 1
      DO 5 L = 2,NNB1
          I = LIST(L,I1)
          S = MAX(TIME - T0(I), 0.0D0)
c          IF (S.EQ.0.0D0) GO TO 5
          S2 = 0.5*S
          S3 = S / 3.0D0
          X(1,I) = X0(1,I)+ S*(X0DOT(1,I)+ S2*(FI(1,I)+ S3*FIDOT(1,I)))
          X(2,I) = X0(2,I)+ S*(X0DOT(2,I)+ S2*(FI(2,I)+ S3*FIDOT(2,I)))
          X(3,I) = X0(3,I)+ S*(X0DOT(3,I)+ S2*(FI(3,I)+ S3*FIDOT(3,I)))
          XDOT(1,I) = X0DOT(1,I)+ S*(FI(1,I)+ S2*FIDOT(1,I))
          XDOT(2,I) = X0DOT(2,I)+ S*(FI(2,I)+ S2*FIDOT(2,I))
          XDOT(3,I) = X0DOT(3,I)+ S*(FI(3,I)+ S2*FIDOT(3,I))
    5 CONTINUE
*
*       Resolve the components of any perturbed pairs (last neighbours).
   10 I = LIST(NNB1,I1)
      IF (I.GT.N) THEN
          JPAIR = I - N
          IF (LIST(1,2*JPAIR-1).GT.0) THEN
              CALL RESOLV(JPAIR,1)
          END IF
          NNB1 = NNB1 - 1
          IF (NNB1.GT.1) GO TO 10
      END IF
      GO TO 40
*
*       Set index of first body (#I1 if NNB = 0 or NNB >= N).
   20 I = I1
*       Predict coordinates & velocities of body #I to order F3DOT.
   25 DT = TIME - T0(I)
      IF ((DT.EQ.0.0D0.AND.IPHASE.LT.4).OR.BODY(I).EQ.0.0D0) GO TO 35
      CALL JPRED(I, 1, X(1,I), XDOT(1,I))
      GO TO 35
      A1 = 0.05*DT
      A2 = 0.25*DT
c      A3 = T0(I) - T0R(I) - SMAX*2
c      A4 = 0.5*A3
*
      DO 30 K = 1,3
          FK = D0(K,I)
          F1DOTK = D1(K,I)
          F2DOTK = 0.5*D2(K,I)
          F3DOTK = ONE6*D3(K,I)
          X(K,I) = ((((F3DOTK*A1 + ONE12*F2DOTK)*DT +
     &                   ONE6*F1DOTK)*DT + 0.5*FK)*DT + X0DOT(K,I))*DT +
     &                   X0(K,I)
          XDOT(K,I)  = (((F3DOTK*A2 + ONE3*F2DOTK)*DT +
     &                            0.5*F1DOTK)*DT + FK)*DT + X0DOT(K,I)
   30 CONTINUE
*
*       Resolve the components of perturbed pairs.
 35   IF (I.GT.N) THEN
          JPAIR = I - N
          IF (LIST(1,2*JPAIR-1).GT.0) THEN
              CALL RESOLV(JPAIR,1)
          END IF
      END IF
*
*       Check whether prediction of single particle or full N.
      IF (NNB.GT.0) THEN
          I = I + 1
          IF (I.LE.NNB) GO TO 25
      END IF
*
   40 RETURN
*
      END
