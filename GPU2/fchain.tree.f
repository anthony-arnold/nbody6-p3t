      SUBROUTINE FCHAIN(I,IR,XI,XIDOT,FIRR,FD)
*
*
*       Irregular force & derivative due to chain.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  XI(3),XIDOT(3),DX(3),DV(3),FIRR(3),FD(3),XIS(3),VIS(3)
      REAL*8  FCI(3),FDC(3)
*
*
*       Use c.m. values for correction of perturbed KS (call from KCPERT).
      I2 = 0
      IF (I.LT.IFIRST) THEN
          IPAIR = KVEC(I)
          IF (I.EQ.2*IPAIR) I2 = 1
          ICM = N + IPAIR
*       Save local variables for individual chain contributions.
          DO 1 K = 1,3
              XIS(K) = XI(K)
              VIS(K) = XIDOT(K)
              XI(K) = X(K,ICM)
              XIDOT(K) = XDOT(K,ICM)
    1     CONTINUE
      END IF
*
*
*     Form terms for the original chain c.m. interaction.
      DO K = 1,3
         FCI(K) = 0
         FDC(K) = 0
      END DO
      CALL CUTOFF(XI,XIDOT,X(1,ICH),XDOT(1,ICH),FCI,FDC,BODY(ICH),RSMIN)
*
*       Subtract force and first derivative from current values.
      DO 10 K = 1,3
          FIRR(K) = FIRR(K) - FCI(K)
          FD(K) = FD(K) - FDC(K)
   10 CONTINUE
*
*       Restore XI & XIDOT for KS components (perturbed case).
      IF (I.LT.IFIRST) THEN
          DO 12 K = 1,3
              XI(K) = XIS(K)
              XIDOT(K) = VIS(K)
   12     CONTINUE
      END IF
*
*       Resolve chain coordinates & velocities using the predicted c.m.
*     IF (IR.EQ.0) THEN
*         CALL XCPRED(0)
*     END IF
*
*       Obtain contributions from all members of the chain.
      DO 20 J = 1,NCH
         CALL CUTOFF(XI,XIDOT,XC(1,J),UC(1,J),FIRR,FD,BODYC(J),RSMIN)
   20 CONTINUE
*
      RETURN
*
      END
