      SUBROUTINE SEARCH(I,IKS)
*
*
*       Close encounter search.
*       -----------------------

      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*

*
*       Increase counter for regularization attempts and set critical step.
      NKSTRY = NKSTRY + 1
      RJMIN2 = 1.0
      DTS = MAX(SMIN,4.0*STEP(I))
*
      FMAX = 0.0
      NCLOSE = 0
*       Find dominant neighbours by selecting all STEP(J) <= 4*STEP(I).
      L = LIST(1,I) + 2
    2 L = L - 1
      IF (L.LT.2) GO TO 10
      IF (LIST(L,I).LE.N) GO TO 4
*
*       Check first whether any c.m. with small step is within range.
      J = LIST(L,I)
*       Include mass condition (STEP may be large).
      IF (STEP(J).GT.DTS.AND.BODY(J).LT.10.0*BODY(I)) GO TO 2
      A1 = X(1,J) - X(1,I)
      A2 = X(2,J) - X(2,I)
      A3 = X(3,J) - X(3,I)
      RIJ2 = A1*A1 + A2*A2 + A3*A3
      IF (RIJ2.GT.RMIN22) GO TO 2
*
      FIJ = BODY(J)/RIJ2
      IF (FMAX.LT.FIJ) FMAX = FIJ
*       Abandon further search if c.m. force exceeds half total force.
      IF (FMAX**2.LT.F(1,I)**2 + F(2,I)**2 + F(3,I)**2) THEN
          NCLOSE = NCLOSE + 1
          JLIST(NCLOSE) = J
          GO TO 2
      ELSE
          GO TO 10
      END IF
*
*       Continue searching single particles with current value of FMAX.
    4 JCOMP = 0
*
      DO 6 K = L,2,-1
          J = LIST(K,I)
          IF (STEP(J).GT.DTS) GO TO 6
          A1 = X(1,J) - X(1,I)
          A2 = X(2,J) - X(2,I)
          A3 = X(3,J) - X(3,I)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          IF (RIJ2.LT.8.0*RMIN22) THEN
              NCLOSE = NCLOSE + 1
              JLIST(NCLOSE) = J
*       Record index of every single body with small step inside 2*RMIN.
              FIJ = (BODY(I) + BODY(J))/RIJ2
              IF (FIJ.GT.FMAX) THEN
                  FMAX = FIJ
*       Save square distance and global index of dominant body.
                  RJMIN2 = RIJ2
                  JCOMP = J
              END IF
          END IF
    6 CONTINUE
*
*       See whether dominant component is a single particle inside RMIN.
      IF (JCOMP.LT.IFIRST.OR.JCOMP.GT.N) GO TO 10
*       Accept one single candidate inside RMIN (which makes PERT = 0).
      IF (RJMIN2.GT.RMIN2) GO TO 10
*
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
*       Only select approaching particles (include nearly circular case).
      RIJMIN = SQRT(RJMIN2)
      IF (RDOT.GT.0.02*SQRT((BODY(I) + BODY(JCOMP))*RIJMIN)) GO TO 10
*
*       Ensure a massive neighbour is included in perturbation estimate.
      BCM = BODY(I) + BODY(JCOMP)
      IF (BODY1.GT.10.0*BCM) THEN
          JBIG = 0
          BIG = BCM
          NNB1 = LIST(1,I) + 1
          DO 20 L = 2,NNB1
              J = LIST(L,I)
              IF (BODY(J).GT.BIG) THEN
                  JBIG = J
                  BIG = BODY(J)
              END IF
   20     CONTINUE
*       Check whether already present, otherwise add to JLIST.
          DO 25 L = 1,NCLOSE
              IF (JLIST(L).EQ.JBIG) THEN
                  JBIG = 0
              END IF
   25     CONTINUE
          IF (JBIG.GT.0) THEN
              NCLOSE = NCLOSE + 1
              JLIST(NCLOSE) = JBIG
          END IF
      END IF
*
*       Evaluate vectorial perturbation due to the close bodies.
      CALL FPERT(I,JCOMP,NCLOSE,PERT)
*
*       Accept #I & JCOMP if the relative motion is dominant (GI < 0.01).
      GI = PERT*RJMIN2/BCM
      IF (GI.GT.0.01) THEN
*         IF (KZ(4).GT.0.AND.TIME-TLASTT.GT.4.44*TCR/FLOAT(N))
*    &                                             CALL EVOLVE(JCOMP,0)
          GO TO 10
      END IF
*
*       Exclude any c.m. body of compact subsystem (TRIPLE, QUAD or CHAIN).
      DO 8 ISUB = 1,NSUB
          NAMEI = NAMES(1,ISUB)
          IF (NAMEI.EQ.NAME(I).OR.NAMEI.EQ.NAME(JCOMP)) GO TO 10
    8 CONTINUE
*
*       Also check possible c.m. body of chain regularization (NAME = 0).
      IF (NCH.GT.0) THEN
          IF (NAME(I).EQ.0.OR.NAME(JCOMP).EQ.0) GO TO 10
      END IF
*
*     Check for equal neighbour lists.
*     Neighbour lists should be sorted so checking element-wise is OK.
c      IF (LIST(1,I).NE.LIST(1,JCOMP)) GO TO 10
c      LJ = 2
c      DO LI=2,LIST(1,I)
c         IF (LIST(LI,I).NE.JCOMP) THEN
c            IF (LIST(LJ,JCOMP).EQ.I) LJ = LJ + 1
c            IF (LIST(LJ,JCOMP).NE.LIST(LI,I)) GO TO 10
c            LJ = LJ + 1
c         END IF
c      END DO
*
*       Save index and increase indicator to denote new regularization.
      ICOMP = I
      IKS = IKS + 1
*
   10 RETURN
*
      END
