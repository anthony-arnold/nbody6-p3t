      SUBROUTINE KSTERM
*
*
*       Termination of KS regularization.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8  SAVE(15),FREG(3),FDREG(3),F1(3),F1DOT(3),A(9)
*
*
*       Copy pair index from COMMON save and define KS components & c.m.
      IPAIR = KSPAIR
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ICM = N + IPAIR
      JMIN = 0
*
*       Save regular c.m. step and form irregular time-step (cf. ADJUST).
      CMSTEP = STEPR(ICM)
      DTCL = 0.04*SQRT(R(IPAIR)**3/BODY(ICM))
*       Switch to arbitrary small value for (CHAIN) collision.
      IF (R(IPAIR).LT.1.0D-06) DTCL = DTMIN
*       Include first-order prediction in regular force & derivative.
      DTR = TIME - T0R(ICM)
      DO 200 K = 1,3
          FREG(K) = FR(K,ICM) + D1R(K,ICM)*DTR
          FDREG(K) = FRDOT(K,ICM) + D2R(K,ICM)*DTR
  200 CONTINUE
*
*       Prepare termination at block time (KS, triple, quad, chain or merge).
      IF (TIME.LE.TBLOCK.AND.IPHASE.LT.9) THEN
          TIME0 = TBLOCK
*       Skip KS integration for unperturbed orbit or T0(I1) at end of block.
          IF (LIST(1,I1).EQ.0.OR.T0(I1).EQ.TIME0) GO TO 3
*
*       See whether the time interval should be modified by KSLOW procedure.
          IF (KSLOW(IPAIR).GT.1) THEN
              IMOD = KSLOW(IPAIR)
              ZMOD = FLOAT(ISLOW(IMOD))
          ELSE
              ZMOD = 1.0
          END IF
*
    1     DT0 = TIME0 - T0(I1)
*       Integrate up to current block time in case interval is too large.
          IF (DT0.GT.STEP(I1)) THEN
              TIME = T0(I1) + STEP(I1)
              H0(IPAIR) = H(IPAIR)
*             Z = -0.5D0*H(IPAIR)*DTAU(IPAIR)**2
*             CALL STUMPF(IPAIR,Z)
              CALL KSINT(I1,1)   ! Dummy argument; hopefully nothing happens.
              DTU = DTAU(IPAIR)
              STEP(I1) = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU
     &                                                   + R(IPAIR))*DTU
              STEP(I1) = ZMOD*STEP(I1)
              DT = STEP(I1)
              CALL STEPK(DT,DTN)
              STEP(I1) = MIN(DTN,SMAX)
*       Restrict increase of R for superfast particles in one block-step.
              IF (H(IPAIR).GT.100.0.AND.R(IPAIR).GT.RMIN) GO TO 3
              GO TO 1
          END IF
*
          IF (NPAIRS.GT.0) GO TO 3
*       Determine the last regularized step by Newton-Raphson iteration.
          DTU = DT0/(R(IPAIR)*ZMOD)
          DTU = MIN(DTU,DTAU(IPAIR))
*       Include rare case of zero interval due to subtracting large values.
          DTU = MAX(DTU,1.0D-10)
          ITER = 0
    2     Y0 = DT0 - ZMOD*((ONE6*TDOT3(IPAIR)*DTU +
     &                             0.5*TDOT2(IPAIR))*DTU + R(IPAIR))*DTU
          YPR = -((0.5*TDOT3(IPAIR)*DTU + TDOT2(IPAIR))*DTU + R(IPAIR))
          YPR = ZMOD*YPR
          DTU = DTU - Y0/YPR
          DT1 = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU +
     &                                                     R(IPAIR))*DTU
          DT1 = ZMOD*DT1
          ITER = ITER + 1
          IF (ABS(DT0 - DT1).GT.1.0E-10*STEP(I1).AND.ITER.LT.10) GO TO 2
*
*       Advance the KS solution to next block time and terminate at TIME0.
          DTAU(IPAIR) = DTU
          STEP(I1) = DT1
          TIME = T0(I1) + DT1
          H0(IPAIR) = H(IPAIR)
*         Z = -0.5D0*H(IPAIR)*DTU**2
*         CALL STUMPF(IPAIR,Z)
          CALL KSINT(I1,1)
          DT = STEP(I1)
          CALL STEPK(DT,DTN)
          STEP(I1) = MIN(DTN,SMAX)
    3     TIME = TIME0
*
*       Predict X & XDOT for body #JCOMP (note TIME = TBLOCK if second call).
          IF (N.LT.5000) THEN
              IF (JCOMP.GE.IFIRST) THEN   ! procedure suppressed for fast code.
                  CALL XVPRED(JCOMP,-1)
                  IF (GAMMA(IPAIR).GT.0.2.AND.JCOMP.LE.N) THEN
                      JMIN = JCOMP
*       Initialize T0, X0 & X0DOT for XVPRED & FPOLY on large perturbation.
                      T0(JCOMP) = TIME
                      DO 4 K = 1,3
                          X0(K,JCOMP) = X(K,JCOMP)
                          X0DOT(K,JCOMP) = XDOT(K,JCOMP)
    4                 CONTINUE
                  END IF
              END IF
          END IF
      END IF
*
*       Predict coordinates and evaluate potential energy w.r.t. perturbers.
      CALL KSRES(IPAIR,J1,J2,0.0D0)
      NP = LIST(1,I1)
      DO 5 L = 1,NP
          JPERT(L) = LIST(L+1,I1)
    5 CONTINUE
      JLIST(1) = I1
      JLIST(2) = I2
      CALL NBPOT(2,NP,POT1)
*
*       Rectify the orbit to yield U & UDOT consistent with binding energy.
      CALL KSRECT(IPAIR)    ! Dangerous to omit unperturbed case (9/16).
*
      CALL RESOLV(IPAIR,2)
      RIJ2 = 0
      VIJ2 = 0
      DO 300 K = 1,3
      RIJ2 = RIJ2 + (X(K,I1) - X(K,I2))**2
      VIJ2 = VIJ2 + (XDOT(K,I1) - XDOT(K,I2))**2
  300 CONTINUE
      RIJ = SQRT(RIJ2)
      HH = 0.5*VIJ2 - BODY(ICM)/RIJ
      WRITE (6,303)  HH, (HH - H(IPAIR))/HH
  303 FORMAT (' ERR!   HH ERR ',1P,3E10.2)
*       Retain final KS variables for explicit restart at merge termination.
      IF (TIME.LE.TBLOCK.AND.IPHASE.EQ.6) THEN
          HM(NMERGE) = H(IPAIR)
          DO 6 K = 1,4
              UM(K,NMERGE) = U(K,IPAIR)
              UMDOT(K,NMERGE) = UDOT(K,IPAIR)
    6     CONTINUE
      END IF
*
*       Check optional diagnostic output for disrupted new hard binary.
      IF (KZ(8).EQ.0) GO TO 10
      IF (LIST(2,I1+1).NE.0.OR.H(IPAIR).GT.0.0) GO TO 10
      IF (GAMMA(IPAIR).GT.0.5.AND.JCOMP.GT.0.OR.IPHASE.EQ.7) THEN
          IF (JCOMP.EQ.0.OR.IPHASE.EQ.7) JCOMP = I1
          K = 0
          IF (JCOMP.GT.N) THEN
              J2 = 2*(JCOMP - N)
              K = LIST(2,J2)
          END IF
          SEMI = -0.5*BODY(ICM)/H(IPAIR)
          EB = -0.5*BODY(I1)*BODY(I2)/SEMI
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          WRITE (8,8)  TIME+TOFF, NAME(I1), NAME(I2), K, NAME(JCOMP),
     &                 BODY(JCOMP), EB, SEMI, R(IPAIR), GAMMA(IPAIR), RI
    8     FORMAT (' END BINARY   T =',F8.1,'  NAME = ',2I6,I3,I6,
     &                      '  M(J) =',F8.4,'  EB =',F10.5,'  A =',F8.5,
     &                          '  R =',F8.5,'  G =',F5.2,'  RI =',F5.2)
      END IF
*
   10 IF (KZ(10).GT.1) THEN
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
*       Include error of the bilinear relation (cf. Book Eqn.(4.30)).
          ERR = U(4,IPAIR)*UDOT(1,IPAIR) - U(3,IPAIR)*UDOT(2,IPAIR) +
     &          U(2,IPAIR)*UDOT(3,IPAIR) - U(1,IPAIR)*UDOT(4,IPAIR)
          WRITE (6,15)  TIME+TOFF, BODY(I1), BODY(I1+1), DTAU(IPAIR),
     &                  R(IPAIR), RI, H(IPAIR), IPAIR, GAMMA(IPAIR),
     &                  STEP(I1), LIST(1,I1), LIST(1,ICM), ERR
   15     FORMAT (/,' END KSREG    TIME =',F8.2,1P,2E9.1,0P,F6.3,1P,
     &                E10.1,0P,F7.2,F9.2,I5,F8.3,1P,E10.1,2I5,1P,E10.2)
      END IF
*
*       Obtain global coordinates & velocities.
      CALL RESOLV(IPAIR,2)
*
*       Correct for differential potential energy due to rectification.
      CALL NBPOT(2,NP,POT2)
*       Add correction term with opposite sign for conservation.
*     ECOLL = ECOLL + (POT2 - POT1)
*     IF (ABS(POT1-POT2).GT.0.0001) WRITE (6,16)  POT1,BE(3),POT1-POT2
*  16 FORMAT (' CORRECT:    POT1 BE3 POT1-POT2  ',2F10.6,F10.6)
*
*       Modify c.m. neighbour radius by density contrast and set new values.
      NNB = LIST(1,ICM) + 1
*       Check predicted neighbour number and form volume ratio.
      NBP = MIN(ALPHA*SQRT(FLOAT(NNB)*RS(ICM))/(RS(ICM)**2),ZNBMAX)
      NBP = MAX(NBP,INT(ZNBMIN))
      A0 = FLOAT(NBP)/FLOAT(NNB)
*       Re-determine neighbour list on zero membership for distant binary.
      IF (LIST(1,ICM).EQ.0) THEN
          RS0 = 0.1*(ABS(X(1,ICM)) + ABS(X(2,ICM)))
          CALL NBLIST(ICM,RS0)
          NNB = LIST(1,ICM) + 1
      END IF
*       Copy all neighbours in the case of merger.
      IF (IPHASE.EQ.6) THEN
          A0 = 1.0
          NBP = NNB - 1
      END IF
      IF (RS(ICM).GT.-100.0*BODY(ICM)/H(IPAIR)) A0 = 1.0
*       Accept old c.m. values for small length scale ratio or H > 0.
      RS(I1) = RS(ICM)*A0**0.3333
      RS(I1+1) = RS(I1)
      RS2 = RS(I1)**2
*
*       Select neighbours for components inside the modified c.m. sphere.
   20 NNB1 = 1
      DO 25 L = 2,NNB
          J = LIST(L,ICM)
          RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
     &                                    (X(3,ICM) - X(3,J))**2
*       Ensure that at least the predicted neighbour number is reached.
          IF (RIJ2.LT.RS2.OR.(NNB1 + NNB - L.LE.NBP.AND.
     &        NNB1.LT.NNBMAX-1)) THEN
              NNB1 = NNB1 + 1
              ILIST(NNB1) = J
          END IF
   25 CONTINUE
*
*       Check that there is space for adding dominant component later.
      IF (NNB1.GE.NNBMAX.AND.IPHASE.NE.6) THEN
          RS2 = 0.9*RS2
          GO TO 20
      END IF
*
*       Reduce pair index, total number & single particle index.
      NPAIRS = NPAIRS - 1
      NTOT = N + NPAIRS
      IFIRST = 2*NPAIRS + 1
      IF (NPAIRS.EQ.0) NNTB = 0
*
*       Save name of components & flag for modifying LISTD in UPDATE.
      JLIST(1) = NAME(I1)
      JLIST(2) = NAME(I1+1)
      JLIST(3) = LIST(2,I1+1)
*
*       Skip adjustment of tables if last or only pair being treated.
      IF (IPAIR.EQ.NPAIRS + 1) GO TO 60
*
*       Move the second component before the first.
      DO 50 KCOMP = 2,1,-1
          I = 2*IPAIR - 2 + KCOMP
*
          DO 30 K = 1,3
              SAVE(K) = X(K,I)
              SAVE(K+3) = X0DOT(K,I)
   30     CONTINUE
*       Current velocity has been set in routine RESOLV.
          SAVE(7) = BODY(I)
          SAVE(8) = RS(I)
          SAVE(9) = RADIUS(I)
          SAVE(10) = TEV(I)
          SAVE(11) = BODY0(I)
          SAVE(12) = TEV0(I)
          SAVE(13) = EPOCH(I)
          SAVE(14) = SPIN(I)
          SAVE(15) = ZLMSTY(I)
          NAMEI = NAME(I)
          KSI = KSTAR(I)
          LAST = 2*NPAIRS - 1 + KCOMP
*
*       Move up global variables of other components.
          DO 40 J = I,LAST
              J1 = J + 1
              DO 35 K = 1,3
                  X(K,J) = X(K,J1)
*       Copy latest X & X0DOT (= 0) of single components for predictor.
                  X0(K,J) = X(K,J)
                  X0DOT(K,J) = X0DOT(K,J1)
                  XDOT(K,J) = XDOT(K,J1)
   35         CONTINUE
              BODY(J) = BODY(J1)
              RS(J) = RS(J1)
              RADIUS(J) = RADIUS(J1)
              TEV(J) = TEV(J1)
              TEV0(J) = TEV0(J1)
              BODY0(J) = BODY0(J1)
              EPOCH(J) = EPOCH(J1)
              SPIN(J) = SPIN(J1)
              ZLMSTY(J) = ZLMSTY(J1)
              NAME(J) = NAME(J1)
              KSTAR(J) = KSTAR(J1)
              STEP(J) = STEP(J1)
              T0(J) = T0(J1)
              K = LIST(1,J1) + 1
              IF (K.EQ.1) K = 2
*       Transfer unmodified neighbour lists (include flag in 2nd comp).
              DO 38 L = 1,K
                  LIST(L,J) = LIST(L,J1)
   38         CONTINUE
   40     CONTINUE
*
*       Set new component index and copy basic variables.
          I = LAST + 1
          DO 45 K = 1,3
              X(K,I) = SAVE(K)
              X0DOT(K,I) = SAVE(K+3)
              XDOT(K,I) = SAVE(K+3)
   45     CONTINUE
          BODY(I) = SAVE(7)
          RS(I) = SAVE(8)
          RADIUS(I) = SAVE(9)
          TEV(I) = SAVE(10)
          BODY0(I) = SAVE(11)
          TEV0(I) = SAVE(12)
          EPOCH(I) = SAVE(13)
          SPIN(I) = SAVE(14)
          ZLMSTY(I) = SAVE(15)
          NAME(I) = NAMEI
          KSTAR(I) = KSI
   50 CONTINUE
*
*       Include removal of the circularization name NAMEC from chaos table.
      IF (KSTAR(ICM).GE.10.AND.NCHAOS.GT.0.AND.IPHASE.NE.6) THEN
*       Note that NAMEC may remain dormant for hierarchical systems.
          II = -ICM
          CALL SPIRAL(II)
      END IF
*
*       Update all regularized variables.
      CALL REMOVE(IPAIR,2)
*
*       Remove old c.m. from all COMMON tables (no F & FDOT correction).
      CALL REMOVE(ICM,3)
*
*       Set new global index of first & second component.
   60 ICOMP = 2*NPAIRS + 1
      JCOMP = ICOMP + 1
*
*       Save c.m. neighbour list for routine FPOLY1/2 (may be renamed below).
      ILIST(1) = NNB1 - 1
      DO 65 L = 1,NNB1
          LIST(L,ICOMP) = ILIST(L)
   65 CONTINUE
*
*       Modify all relevant COMMON list arrays.
      CALL UPDATE(IPAIR)
*
*       Check replacing of single KS component by corresponding c.m.
   70 IF (NNB1.GT.1.AND.LIST(2,ICOMP).LT.ICOMP) THEN
          J2 = LIST(2,ICOMP)
          J = KVEC(J2) + N
          IF (NNB1.EQ.2) THEN
              LIST(2,ICOMP) = J
          ELSE
              L = 2
   80         JNEXT = LIST(L+1,ICOMP)
              IF (JNEXT.LT.J) THEN
                  LIST(L,ICOMP) = JNEXT
                  L = L + 1
                  IF (L.LT.NNB1) GO TO 80
                  LIST(NNB1,ICOMP) = J
              ELSE IF (JNEXT.EQ.J) THEN
                  NNB1 = NNB1 - 1
                  DO 85 LL = L,NNB1
                      LIST(LL,ICOMP) = LIST(LL+1,ICOMP)
   85             CONTINUE
              ELSE
                  LIST(L,ICOMP) = J
              END IF
*       Check again until first neighbour > ICOMP.
              GO TO 70
          END IF
      END IF
*
*       Make space for dominant component and copy members to JCOMP list.
      DO 90 L = NNB1,2,-1
          LIST(L+1,ICOMP) = LIST(L,ICOMP)
          LIST(L+1,JCOMP) = LIST(L,ICOMP)
   90 CONTINUE
*
*       Set dominant component in first location and specify membership.
      LIST(2,ICOMP) = JCOMP
      LIST(2,JCOMP) = ICOMP
      LIST(1,ICOMP) = NNB1
      LIST(1,JCOMP) = NNB1
*
      TIME = TBLOCK
*       Initialize T0, T0R, X0 & X0DOT for both components.
      T0(ICOMP) = TIME
      T0(JCOMP) = TIME
      T0R(ICOMP) = TIME
      T0R(JCOMP) = TIME
      DO 95 K = 1,3
          X0(K,ICOMP) = X(K,ICOMP)
          X0(K,JCOMP) = X(K,JCOMP)
          X0DOT(K,ICOMP) = XDOT(K,ICOMP)
          X0DOT(K,JCOMP) = XDOT(K,JCOMP)
   95 CONTINUE
*
*       Predict current coordinates & velocities for the neighbours.
      CALL XVPRED(ICOMP,NNB1)
*
*       Form new force polynomials (skip triple, quad, merge & chain).
      IF (N.LT.5000) THEN
*       Use old method for non-standard termination (not from nbody6.f).
*
          CALL FPOLY1(ICOMP,JCOMP,2)
          CALL FPOLY2(ICOMP,JCOMP,2)
*
*       Improve force polynomials of strong perturber after rectification.
          IF (JMIN.GE.IFIRST) THEN
              CALL FPOLY1(JMIN,JMIN,0)
              CALL FPOLY2(JMIN,JMIN,0)
          END IF
*
      ELSE
*
*       Treat each component in turn.
          I = ICOMP
  100     DO 105 K = 1,3
              FI(K,I) = 0.0
              D1(K,I) = 0.0
              FR(K,I) = FREG(K)
              FRDOT(K,I) = FDREG(K)
              D0R(K,I) = FREG(K)
              D1R(K,I) = FDREG(K)  ! ignore old tidal force.
  105     CONTINUE
*
*       Obtain irregular force & first derivative for body #I.
          KDUM = 0
          NNB = LIST(1,I)
*       Loop over neighbours only (cf. FPOLY1).
          DO 130 L = 2,NNB+1
              J = LIST(L,I)
              IF (J.GT.N) THEN
                  JPAIR = J - N
*       Use c.m. approximation for unperturbed binary.
                  IF (LIST(1,2*JPAIR-1).GT.0) THEN
                      KDUM = 2*JPAIR - 1
                      J = KDUM
                  END IF
              END IF
*
  110         DO 115 K = 1,3
                  A(K) = X(K,J) - X(K,I)
                  A(K+3) = XDOT(K,J) - XDOT(K,I)
  115         CONTINUE
*
              A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
              A(8) = BODY(J)*A(7)*SQRT(A(7))
              A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
*       Accumulate irregular force and first derivative.
              DO 120 K = 1,3
                  F1(K) = A(K)*A(8)
                  F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
                  FI(K,I) = FI(K,I) + F1(K)
                  D1(K,I) = D1(K,I) + F1DOT(K)
  120         CONTINUE
*
*       Check for KS component.
              IF (J.EQ.KDUM) THEN
                  J = J + 1
                  GO TO 110
              END IF
  130     CONTINUE
*
*       Check option for external force.
          IF (KZ(14).GT.0) THEN
              CALL XTRNLD(I,I,1)
          END IF
*
*       Truncate the orbital time-step (formed from Kepler period).
          CALL STEPK(DTCL,DTN)
          STEP(I) = MIN(DTN,SMAX)
          ITER = 0
*       Perform commensurability check.
  140     IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
              STEP(I) = 0.5D0*STEP(I)
              ITER = ITER + 1
              IF (ITER.LT.16.OR.STEP(I).GT.DTK(40)) GO TO 140
              STEP(I) = DTK(40)
          END IF
*
*       Initialize next irregular time and copy c.m. regular time-step.
          TNEW(I) = TIME + STEP(I)
          VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
          DT0 = 0.5*RS(I)/SQRT(VI2)
          DT = MIN(DT0,SMAX,CMSTEP)
          CALL STEPK(DT,DTN)
          STEPR(I) = DTN
          ITER = 0
*       Perform commensurability check.
  145     IF (DMOD(TIME,STEPR(I)).NE.0.0D0) THEN
              STEPR(I) = 0.5D0*STEPR(I)
              ITER = ITER + 1
              IF (ITER.LT.16.OR.STEPR(I).GT.DTK(40)) GO TO 145
              STEPR(I) = DTK(40)
          END IF
*       Set force and all derivatives.
          DO 150 K = 1,3
              F(K,I) = FI(K,I) + FR(K,I)
              FDOT(K,I) = D1(K,I) + D1R(K,I)
              F(K,I) = 0.5*F(K,I)
              FDOT(K,I) = ONE6*FDOT(K,I)
              D0(K,I) = FI(K,I)
              FIDOT(K,I) = D1(K,I)
              D2(K,I) = 0.0
              D3(K,I) = 0.0
              D0R(K,I) = FREG(K)
              FR(K,I) = FREG(K)
              FRDOT(K,I) = FDREG(K)
              D1R(K,I) = FDREG(K)
              D2R(K,I) = 0.0
              D3R(K,I) = 0.0
  150     CONTINUE
*
*       Treat second component in the same way.
          IF (I.EQ.ICOMP) THEN
              I = JCOMP
              GO TO 100
          END IF
*
*       Improve force polynomials of strong perturber after rectification.
*         IF (JMIN.GE.IFIRST) THEN
*             CALL FPOLY1(JMIN,JMIN,0)
*             CALL FPOLY2(JMIN,JMIN,0)
*         END IF
      END IF
*
*       Check updating of chain perturber list.
      IF (NCH.GT.0) THEN
          CALL CHFIND
      END IF
*
      RETURN
*
      END
