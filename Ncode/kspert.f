      SUBROUTINE KSPERT(I1,NNB0,XI,VI,FP,FD)
*
*
*       Perturbation on KS pair.
*       ------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,OMEGA,
     &       MBULGE,RBULGE,MDISK,SCALEA,SCALEB,MHALO,RHALO,VCIRC
      REAL*8  XI(6),VI(6),FP(6),FD(6),TF(3),TD(3)
      REAL*8  XG(3),XGDOT(3),FM(3),FMD(3),FS(3),FSD(3)
*
*
*       Initialize the perturbing force & first derivative.
      DO 10 K = 1,6
          FP(K) = 0.0D0
          FD(K) = 0.0D0
   10 CONTINUE
*
*       Set index of the last single perturber.
      NNB2 = NNB0 + 1
   15 IF (LIST(NNB2,I1).LE.N) GO TO 20
      NNB2 = NNB2 - 1
      IF (NNB2.GT.1) GO TO 15
*       Include special case of only c.m. perturbers.
      GO TO 30
*
*       Obtain the perturbation from single particles.
   20 DO 25 L = 2,NNB2
          K = LIST(L,I1)
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(1) = FP(1) + A1*A6
          FP(2) = FP(2) + A2*A6
          FP(3) = FP(3) + A3*A6
          V1 = XDOT(1,K) - VI(1)
          V2 = XDOT(2,K) - VI(2)
          V3 = XDOT(3,K) - VI(3)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(1) = FD(1) + (V1 - A1*A9)*A6
          FD(2) = FD(2) + (V2 - A2*A9)*A6
          FD(3) = FD(3) + (V3 - A3*A9)*A6
*       Perturbation on first component.
*
          A1 = X(1,K) - XI(4)
          A2 = X(2,K) - XI(5)
          A3 = X(3,K) - XI(6)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(4) = FP(4) + A1*A6
          FP(5) = FP(5) + A2*A6
          FP(6) = FP(6) + A3*A6
          V1 = XDOT(1,K) - VI(4)
          V2 = XDOT(2,K) - VI(5)
          V3 = XDOT(3,K) - VI(6)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(4) = FD(4) + (V1 - A1*A9)*A6
          FD(5) = FD(5) + (V2 - A2*A9)*A6
          FD(6) = FD(6) + (V3 - A3*A9)*A6
*       Perturbation on second component.
   25 CONTINUE
*
*       See whether to include any remaining c.m. perturbers.
      IF (NNB2.GT.NNB0) GO TO 40
*
   30 KDUM = 0
*       Dummy index to enable summation of c.m. or resolved components.
      NNB3 = NNB2 + 1
      DO 35 L = NNB3,NNB0+1
          K = LIST(L,I1)
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*       See whether c.m. approximation applies (ignore unperturbed case).
          J = K - N
          IF (RIJ2.GT.CMSEP2*R(J)**2) THEN
              V1 = XDOT(1,K) - VI(1)
              V2 = XDOT(2,K) - VI(2)
              V3 = XDOT(3,K) - VI(3)
              A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
              GO TO 33
          END IF
*
*       Resolve pair #J and sum over individual components.
          CALL KSRES2(J,J1,J2,RIJ2)
          KDUM = J1
          K = KDUM
   32     A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          V1 = XDOT(1,K) - VI(1)
          V2 = XDOT(2,K) - VI(2)
          V3 = XDOT(3,K) - VI(3)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
   33     A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(1) = FP(1) + A1*A6
          FP(2) = FP(2) + A2*A6
          FP(3) = FP(3) + A3*A6
          FD(1) = FD(1) + (V1 - A1*A9)*A6
          FD(2) = FD(2) + (V2 - A2*A9)*A6
          FD(3) = FD(3) + (V3 - A3*A9)*A6
*       Perturbation on first component.
*
          A1 = X(1,K) - XI(4)
          A2 = X(2,K) - XI(5)
          A3 = X(3,K) - XI(6)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(4) = FP(4) + A1*A6
          FP(5) = FP(5) + A2*A6
          FP(6) = FP(6) + A3*A6
          V1 = XDOT(1,K) - VI(4)
          V2 = XDOT(2,K) - VI(5)
          V3 = XDOT(3,K) - VI(6)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(4) = FD(4) + (V1 - A1*A9)*A6
          FD(5) = FD(5) + (V2 - A2*A9)*A6
          FD(6) = FD(6) + (V3 - A3*A9)*A6
*       Perturbation on second component.
*
          IF (K.EQ.KDUM) THEN
              K = K + 1
              GO TO 32
          END IF
   35 CONTINUE
*
*       Check perturbation correction due to regularized chain.
   40 IF (NCH.GT.0) THEN
          DO 45 L = 2,NNB2
              J = LIST(L,I1)
              IF (J.GT.ICH) GO TO 50
              IF (J.EQ.ICH) THEN
                  J1 = I1
                  CALL FCHAIN(J1,1,XI(1),VI(1),FP(1),FD(1))
                  J1 = J1 + 1
                  CALL FCHAIN(J1,1,XI(4),VI(4),FP(4),FD(4))
                  GO TO 50
              END IF
   45     CONTINUE
      END IF 
*
*       Set the relative perturbing force and first derivative.
   50 DO 55 K = 1,3
          FP(K) = FP(K) - FP(K+3)
          FD(K) = FD(K) - FD(K+3)
          TF(K) = 0.0D0
          TD(K) = 0.0D0
   55 CONTINUE
*
*       See whether the linearized perturbation should be included.
      IF (KZ(14).GT.0.AND.KZ(14).LT.3) THEN
          Q1 = XI(1) - XI(4)
          Q3 = XI(3) - XI(6)
          CALL XTRNLP(Q1,Q3,TF)
*
*       Use same formalism for the first derivative (omit Coriolis force).
          VX = VI(1) - VI(4)
          VZ = VI(3) - VI(6)
          CALL XTRNLP(VX,VZ,TD)
          DO 60 K = 1,3
             FP(K) = FP(K) + TF(K)
             FD(K) = FD(K) + TD(K)
   60     CONTINUE
      END IF
*
*       Check optional Plummer potential.
      IF (MP.GT.0.0) THEN
          RI2 = AP2
          RRDOT = 0.0
*       Form one central distance and scalar product of relative motion.
          DO 65 K = 1,3
              RI2 = RI2 + XI(K)**2
              RRDOT = RRDOT + (XI(K) - XI(K+3))*(VI(K) - VI(K+3))
   65     CONTINUE
          ZF = 1.0/RI2
*       Write current mass inside RI as MP*R3*ZF^{3/2} (Heggie & Hut p.73).
          FMP = MP*ZF*SQRT(ZF)
          DO 70 K = 1,3
              XREL = XI(K) - XI(K+3)
              VREL = VI(K) - VI(K+3)
              FP(K) = FP(K) - XREL*FMP
              FD(K) = FD(K) - (VREL - 3.0*RRDOT*ZF*XREL)*FMP
   70     CONTINUE
      END IF
*
*       Include possible contribution from central point-mass.
      IF (GMG.GT.0.0D0) THEN
          DO 80 K = 1,3
              XG(K) = RG(K) + XI(K)   ! Note XG used for each KS component.
              XGDOT(K) = VG(K) + VI(K)
   80     CONTINUE
          CALL FNUC(XG,XGDOT,FS,FSD)
          DO 85 K = 1,3
              XG(K) = RG(K) + XI(K+3)
              XGDOT(K) = VG(K) + VI(K+3)
   85     CONTINUE
          CALL FNUC(XG,XGDOT,FM,FMD)
          DO 90 K = 1,3
              FP(K) = FP(K) + (FS(K) - FM(K))
              FD(K) = FD(K) + (FSD(K) - FMD(K))
   90     CONTINUE
      END IF
*
      RETURN
*
      END
