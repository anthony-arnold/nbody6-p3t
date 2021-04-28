      SUBROUTINE GPUINT(I,FIRR,FD,FREG,FDR)
*
*
*     Regular GPU leapfrog force integrator.
*     --------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tree.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &     LISTC(LMAX)
      REAL*8 FRX(3),FDX(3)
      REAL*8 FIRR(3),FREG(3),FD(3),FDR(3),XI(3),XIDOT(3)
*
      DO K = 1,3
         XI(K) = X0(K,I)
         XIDOT(K) = X0DOT(K,I)
      END DO
*
*     Update the neighbour list.
      NNB = LIST(1,I)
c     DTR = TIME - T0R(I)
      DTR = SMAX*4
*
*     Choose appropriate force loop for single particle or c.m. body.
      IF (I.GT.N) THEN
*     Treat unperturbed KS in the single particle approximation.
         I1 = 2*(I - N) - 1
         IF (LIST(1,I1).GT.0) THEN
*     Correct irregular force on active c.m. body (same form as NBINT).
            CALL CMFIRR(I,I1,FIRR,FD)
            GO TO 20
         END IF
      END IF
*
*     Add corrections from c.m. particles (singles or unperturbed).
      IF (NPAIRS.GT.0) THEN
         CALL CMFIRR2(I,NNB,LIST(1,I),XI,XIDOT,FIRR,FD)
      END IF
*
*     Include treatment for regularized subsystem.
      IF (NCH.GT.0) THEN
*     Distinguish between chain c.m. and any other particle.
         IF (NAME(I).EQ.0) THEN
            CALL CHFIRR(I,1,XI,XIDOT,FIRR,FD)
         ELSE
*     Search the chain perturber list for #I.
            NP1 = LISTC(1) + 1
            DO 15 L = 2,NP1
               J = LISTC(L)
               IF (J.GT.I) GO TO 20
               IF (J.EQ.I) THEN
                  CALL FCHAIN(I,1,XI,XIDOT,FIRR,FD)
                  GO TO 20
               END IF
 15         CONTINUE
         END IF
      END IF
*
*     Check optional interstellar clouds.
 20   IF (KZ(13).GT.0) THEN
         CALL FCLOUD(I,FREG,FDR,1)
      END IF
*
*     See whether an external force should be added.
      IF (KZ(14).GT.0) THEN
*     Save current values for deriving work done by tides (#14 = 3).
         IF (KZ(14).EQ.3) THEN
            DO 22 K = 1,3
               FRX(K) = FREG(K)
               FDX(K) = FDR(K)
 22         CONTINUE
         END IF
*
*     !$omp critical
*     Obtain the tidal perturbation (force and first derivative).
         CALL XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDR,2)
*     !$omp end critical
*
*     Form rate of tidal energy change during last regular step.
         IF (KZ(14).EQ.3) THEN
            WDOT = 0.0
            W2DOT = 0.0
*     W3DOT = 0.0
            DO 24 K = 1,3
               PX = FREG(K) - FRX(K)
               DPX = FDR(K) - FDX(K)
               WDOT = WDOT + XIDOT(K)*PX
               W2DOT = W2DOT + (FREG(K) + FIRR(K))*PX + XIDOT(K)*DPX
*     Suppress incomplete third-order term which gives worse result.
*     W3DOT = W3DOT + 2.0*(FREG(K) + FIRR(K))*DPX +
*     &                            (FDR(K) + FD(K))*PX
 24         CONTINUE
*     Note: second-order term derived by Douglas Heggie (Aug/03).
         END IF
*
*     Include the force from optional gaseous Plummer potential.
         IF (KZ(14).GE.3) THEN
            CALL XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDR,-1)
         END IF
      END IF
*
*
*     Assume small mass at centre for special case of no neighbours.
      ZFI = FIRR(1)**2 + FIRR(2)**2 + FIRR(3)**2
      IF (NNB.EQ.0.OR.ZFI.EQ.0.0D0) THEN
         NBVOID = NBVOID + 1
      END IF
*
*     Accumulate tidal energy change for general galactic potential.
 70   IF (KZ(14).EQ.3) THEN
*     Note: Taylor series at end of interval with negative argument.
!     $omp critical
         ETIDE = ETIDE + BODY(I)*(0.5*W2DOT*DTR - WDOT)*DTR
!     $omp end critical
*     Note: integral of Taylor series for V*P using final values.
      END IF
*
*     Set forces.
      DTSQ = DTR**2
      DT6 = 6.0D0/(DTR*DTSQ)
      DT2 = 2.0D0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DTR13 = ONE3*DTR
      DO K = 1,3
         FDR(K) = (FRDOT(K,I) - (FREG(K) - FIRR(K)))*DTR
*       Subtract change of neighbour force to form actual first derivative.
         DFR = FR(K,I) - (FIRR(K) - FI(K,I)) - (FREG(K) - FIRR(K))
         FDR0 = FDR(K) - (FIDOT(K,I) - FD(K))
*
         FRD = FRDOT(K,I)
         SUM = FRD + FDR0
         AT3 = 2.0D0*DFR + DTR*SUM
         BT2 = -3.0D0*DFR - DTR*(SUM + FRD)
*
         FRDOT(K,I) = FDR(K)
         FR(K,I) = FREG(K) - FIRR(K)
         D0R(K,I) = FREG(K) - FIRR(K)
         D1R(K,I) = FRDOT(K,I)
         D2R(K,I) = (3.0D0*AT3 + BT2)*DT2
         D3R(K,I) = AT3*DT6
*
         FI(K,I) = FIRR(K)
         FIDOT(K,I) = FD(K)
         D0(K,I) = FIRR(K)
         D1(K,I) = FD(K)
*
         F(K,I) = 0.5D0*(FR(K,I) + FIRR(K))
         FDOT(K,I) = ONE6*(FRDOT(K,I) + FD(K))
      END DO
*
      RETURN
*
      END
