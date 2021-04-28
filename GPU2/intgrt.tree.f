      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tree.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/PREDICT/ TPRED(NMAX)
      INTEGER ISTAT(NMAX)
      PARAMETER (NPMAX=16)
      REAL*4 TINTG(2)
      REAL*8   XNEW(3,NMAX),XDNEW(3,NMAX),RCUT(NMAX)
      INTEGER  NXTLST(NMAX),LISTQ(NMAX),NL(20)
      INTEGER  IRR(NMAX)
      LOGICAL LOOP,LSTEPM,STEPXCOMM
      SAVE IQ,ICALL,NQ,LQ,LOOP,LSTEPM,STEPM,ISAVE,JSAVE,ISTART,TBH
      DATA IQ,ICALL,LQ,LOOP,LSTEPM,STEPM /0,2,11,.TRUE.,.FALSE.,0.03125/
      DATA ISAVE,JSAVE,ISTART /0,0,0/
      SAVE TLISTQ,TMIN,ICOMP0,JCOMP0,NPACT,XNEW,XDNEW,RCUT
      SAVE NXTLST,LISTQ,IRR,ISTAT
*     SAVE CPRED,CPRED2,CPRED3
*     DATA CPRED,CPRED2,CPRED3 /0.0D0,0.0D0,0.0D0/
*
*
*       Update quantized value of STEPM for large N (first time only).
      IF (.NOT.LSTEPM.AND.NZERO.GT.1024) THEN
          K = INT((FLOAT(NZERO)/1024.0)**0.333333)
          STEPM = 0.03125D0/2**(K-1)
          LSTEPM = .TRUE.
       END IF
*
*       Open the GPU libraries on each new run (note nnbmax = NN is printed).
      IF (ISTART.EQ.0) THEN
          NN = NTOT + 100
          CALL GPUNB_OPEN(NN, NNBMAX)
*       Set larger value for GPUIRR (note further possible increase of NTOT).
          NNN = NTOT + 100
*       Increase NNN further for cold start to deal with many KS.
          IF (QVIR.LT.0.01.AND.TIME.EQ.0.0D0) NNN = NTOT + 500
*       Allow for missing primordial binaries.
          NNN = NNN + NBIN0 - NPAIRS
          CALL GPUIRR_OPEN(NNN,LMAX)
          ISTART = 1
*       Define parameter for predicting active particles in library.
          IF (N.LE.5000) THEN
              NPACT = 75
          ELSE IF (N.LE.50000) THEN
              NPACT = 100
          ELSE IF (N.LE.100000) THEN
              NPACT = 150
          ELSE
              NPACT = 200
          END IF
          TBH = TTOT
          CALL STOPWATCH(TBEG)
          WRITE (6,2) TBEG
 2        FORMAT(//,9X,'START INTGRT ', E11.3)
      END IF
*
*       Search for high velocities after escape or KS/chain termination.
  999 IF (KZ(37).GT.0.AND.(IPHASE.EQ.-1.OR.IPHASE.GE.2)) THEN
          CALL HIVEL(0)
      END IF
*
*       Initialize new force times and predict to previous TIME.
      IF (IPHASE.EQ.1) THEN
*       Treat case of new KS with components from ICOMP0 & JCOMP0.
          I = ICOMP0
*       Avoid updating unchanged locations.
 1000     IF (I.GE.IFIRST) THEN
             TNEW(I) = T0(I) + STEP(I)
             TPRED(I) = -1.0
             CALL JPRED(I)
             CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),
     &            FI(1,I),FIDOT(1,I),
     &            BODY(I),T0(I))
          END IF
*       Select second component or new c.m. (not connected with ICOMP0).
          IF (I.EQ.ICOMP0) THEN
              I = JCOMP0
              GO TO 1000
          ELSE IF (I.EQ.JCOMP0) THEN
              I = NTOT
              GO TO 1000
          END IF
!$omp parallel do private(J)
          DO 1010 J = IFIRST,NTOT
              CALL GPUIRR_SET_LIST(J,LIST(1,J))
 1010     CONTINUE
!$omp end parallel do
*
*       Remove ICOMP0 & JCOMP0 from LISTQ and add NTOT.
          CALL REPAIR(ICOMP0,JCOMP0,NTOT,0,NQ,LISTQ,TMIN)
*       Check KS termination (locations IFIRST & IFIRST + 1).
      ELSE IF (IPHASE.EQ.2) THEN
*       Update relevant variables for the KS components and later c.m.
          I1 = IFIRST
          I2 = I1 + 1
 1015     DO 1030 I = I1,I2
             TNEW(I) = T0(I) + STEP(I)
             TPRED(I) = -1.0
             CALL JPRED(I)
             CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),
     &            FI(1,I),FIDOT(1,I),
     &            BODY(I),T0(I))
 1030     CONTINUE
*
*       Include all the more recent c.m. bodies.
          IF (I1.EQ.IFIRST) THEN
              I1 = N + KSPAIR
              I2 = NTOT
              GO TO 1015
*       Note loop is ignored for case of last c.m. (KSPAIR = NPAIRS + 1).
          END IF
*
*       Re-send all neighbour lists using parallel directive (cheap loop).
!$omp parallel do private(I)
          DO 1035 I = IFIRST,NTOT
              CALL GPUIRR_SET_LIST(I,LIST(1,I))
 1035     CONTINUE
!$omp end parallel do
*
*       Add two first single particles to LISTQ and remove terminated c.m.
          CALL REPAIR(N+KSPAIR,0,IFIRST,IFIRST+1,NQ,LISTQ,TMIN)
      END IF
*
*       Update all SET_JP & SET_LIST except for frequent IPHASE = 1 or 2.
      IF (IPHASE.LE.0.OR.IPHASE.GE.3) THEN
*       Include at TIME = 0, CHAIN (7, 8), COLL & COAL, escape & MERGE/RESET.
!$omp parallel do private(I)
*       Note IPHASE = 3 also after ENFORCED KS in ADJUST (if no escape).
          DO 1040 I = IFIRST,NTOT
             CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),
     &            FI(1,I),FIDOT(1,I),
     &            BODY(I),T0(I))
             CALL GPUIRR_SET_LIST(I,LIST(1,I))
*     Initialize end-point of new times and predict all at previous TIME.
              TNEW(I) = T0(I) + STEP(I)
              TPRED(I) = -1.0
              CALL JPRED(I)
*       Note that this prediction is strictly not necessary.
 1040     CONTINUE
!$omp end parallel do
*       Prescribe level search on return, except for new and terminated KS.
          LOOP = .TRUE.
          TLISTQ = TIME
*       Refresh chain perturber predictions and LISTC after each RETURN.
*         IF (NCH.GT.0) THEN   ! Switched to routine CHAIN (02/16).
*             CALL XCPRED(2)
*             CALL CHLIST(ICH)
*         END IF
       END IF
*
*       Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
      TPREV = TIME
*       Take close encounter step as typical lower value for time-step list.
      CALL STEPK(DTMIN,DTM)
*
*       Determine level for the smallest step (ignore extreme values).
      LQS = 20
      DO 1050 L = 1,20
          IF (DTM.EQ.DTK(L)) THEN
              LQS = L
          END IF
 1050 CONTINUE
*
*       Specify upper level for optimized membership and reset IQ.
      LQB = MAX(1, LQS - 8)
      LQS = LQB + 8
      IF (IQ.LT.0) ICALL = 0
      IQ = 0
*
*       Check updating new list of block steps with T0 + STEP =< TLISTQ.
    1 ICALL = ICALL + 1
*       Reset TMIN second & third time after change to catch new chain step.
      IF (TIME.GE.TLISTQ.OR.ICALL.LE.3) THEN
*       Update interval by optimization at major times (sqrt of N-NPAIRS).
          IF (DMOD(TLISTQ,2.0D0).EQ.0.0D0.OR.LOOP) THEN
              LOOP = .FALSE.
              DO 10 L = 1,20
                  NL(L) = 0
   10         CONTINUE
              DO 14 I = IFIRST,NTOT
*       Count steps at five different levels for the smallest values.
                 DO 12 L = LQB,LQS
                    STEMP = STEP(I)
                    IF (STEMP.LT.DTK(L)) NL(L) = NL(L) + 1
   12             CONTINUE
   14         CONTINUE
              NLSUM = 0
*       Determine interval by summing smallest steps until near sqrt(N-N_b).
              NSQ = INT(SQRT(FLOAT(N - NPAIRS)))
              LQ = LQS
              DO 15 L = LQS,LQB,-1
                  NLSUM = NLSUM + NL(L)
                  IF (NLSUM.LE.NSQ) LQ = L
   15         CONTINUE
*             WRITE (6,16)  TIME+TOFF,NQ,NLSUM,LQ,(NL(K),K=LQB,LQS)
*  16         FORMAT (' LEVEL CHECK:    T NQ NLSUM LQ NL  ',
*    &                                  F9.3,3I5,2X,7I4)
          END IF
*
*       Increase interval by optimized value.
          NQ = 0
          TMIN = 1.0D+10
   18     TLISTQ = TLISTQ + DTK(LQ)
          DO 20 I = IFIRST,NTOT
              IF (TNEW(I).LE.TLISTQ.AND.BODY(I).GT.0.0D0) THEN
                  NQ = NQ + 1
                  LISTQ(NQ) = I
                  TMIN = MIN(TNEW(I),TMIN)
              ELSE IF (TNEW(I).LT.TADJ.AND.BODY(I).EQ.0.0D0) THEN
              WRITE (76,19)  IPHASE, NAME(I), TNEW(I), BODY(I)*SMU
   19         FORMAT (' ZERO-ZERO!    IPH NM TNEW M ',I4,I7,1P,2E10.2)
              TNEW(I) = 1.0D+06
              STEP(I) = 1.0D+06
              CALL FLUSH(76)
              END IF
   20     CONTINUE
*       Increase interval in rare case of zero membership.
          IF (NQ.EQ.0) GO TO 18
*       Make a slight adjustment for high levels and small membership.
          IF (LQ.LE.15.AND.NQ.LE.2) LQ = MAX(1, LQ - 1)
      END IF
*
*     TTT = TIME     ! Used for profiling (needs wtime.o).
*     TLAST = TIME   ! Activate on negative time test.
*
      DO 1220 L = 1,NQ
          J = LISTQ(L)
          IF (TNEW(J).LT.TADJ.AND.BODY(J).EQ.0.0D0) THEN
              WRITE (77,1215)  L, NQ, NAME(J), BODY(J)*SMU
 1215         FORMAT (' ZERO!!   L NQ NM M  ',2I5,I7,1P,E10.2)
              CALL FLUSH(77)
              TNEW(J) = 1.0D06
          END IF
 1220  CONTINUE
*     Find all particles in next block (TNEW = TMIN) and set TIME.
      CALL INEXT(NQ,LISTQ,TMIN,NXTLEN,NXTLST)
*
*       Set new time and save block time (for regularization terminations).
      I = NXTLST(1)
      TMIN = 0
      TIME = T0(I) + STEP(I)
      TBLOCK = TIME
      TTOT = TIME + TOFF
*
*       Include diagnostics for negative TIME and block-step diagnostics.
*     IF (TIME.LT.TLAST.AND.TIME.LT.TLISTQ) THEN
*     WRITE (6,1300) NXTLEN, NQ, TIME, TIME-TLAST,STEP(I)
*1300 FORMAT (' NEGATIVE!    LEN NQ T T-TL S ',2I6,F10.5,1P,2E10.2)
*     CALL FLUSH(6)
*     END IF
c      IN = NXTLST(NXTLEN)
c      WRITE (6,22)  IN,NAME(IN),NXTLEN,NSTEPU,NSTEPI,TIME,
c     &     STEP(IN),STEPR(IN)
c 22   FORMAT (' INTGRT   I NM LEN #U #I T S SR  ',
c     &  3I6,2I11,F10.4,1P,2E10.2)
c      CALL FLUSH(6)
*
*       Terminate on small irregular time-step (means problems).
      IF (STEP(I).LT.1.0D-11) THEN
          WRITE (6,24)  I, NAME(I), NXTLEN, NSTEPR, STEP(I), STEPR(I)
   24     FORMAT (' SMALL STEP!!    I NAME LEN #R SI SR ',
     &                              3I6,I11,1P,2E10.2)
          WRITE (6,1800)  TIME, T0R(I)+STEPR(I)
 1800     FORMAT (' CHECK!!   T T0R+SR  ',2F14.6)
          STOP
      END IF
*
*     TT2 = DBLTIM()
*     CPRED3 = CPRED3 + (TT2 - TT1)
*       Re-determine list if current time exceeds boundary.
      IF (TIME.GT.TLISTQ) GO TO 1
*
*       Check option for advancing interstellar clouds.
      IF (KZ(13).GT.0) THEN
          CALL CLINT
      END IF
*
*       Check optional integration of cluster guiding centre.
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
          IF (KZ(14).EQ.3.AND.STEPXCOMM(TIME,STEPX)) THEN
              CALL GCINT
          END IF
*       Include mass loss by gas expulsion (Kroupa et al. MN 321, 699).
          IF (MPDOT.GT.0.0D0.AND.TIME + TOFF.GT.TDELAY) THEN
              CALL PLPOT1(PHI1)
              MP = MP0/(1.0 + MPDOT*(TIME + TOFF - TDELAY))
*       Adjust the tidal radius for strict correctness.
              RTIDE = RTIDE0*((ZMASS + MP)/(1.0 + MP0))**0.3333
*       Replace by exponential mass loss for faster decrease.
*             DT = TIME + TOFF - TDELAY
*             MP = MP0*EXP(-MPDOT*DT)
              CALL PLPOT1(PHI2)
*       Add differential correction for energy conservation.
              EMDOT = EMDOT + (PHI1 - PHI2)
          END IF
      END IF
*
*       Include commensurability test (may be suppressed if no problems).
c      IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
c          WRITE (6,25)  I, NAME(I), NSTEPI, TIME, STEP(I), TIME/STEP(I)
c   25     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
c     &                        2I6,I11,F12.5,1P,E9.1,0P,F16.4)
c          STOP
c      END IF
*
*       Check for new regularization at end of previous block.
      IF (IKS.GT.0) THEN
          IPHASE = 1
*       Copy the saved KS component indices and time.
          ICOMP = ISAVE
          JCOMP = JSAVE
          TIME = TSAVE
          ICOMP0 = ICOMP
          JCOMP0 = JCOMP
          GO TO 100
      END IF
*
*     TT1 = DBLTIM()
*       Check next adjust time before beginning a new block.
      IF (TIME.GT.TADJ) THEN
          TIME = TADJ
          IPHASE = 3
          GO TO 100
      END IF
*
*       Check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
         TIME = TNEXT
          CALL OUTPUT
          GO TO 1
      END IF
*
*     See whether to perform a leapfrog step.
      IF (TBLCKR.EQ.TPREV) THEN
         CALL STOPWATCH(TBEG)
*
*     Give a velocity kick to regular force candidates.
!     $omp parallel do schedule(runtime)
         DO J = IFIRST,NTOT
            CALL RKICK(J)
            T0R(J) = TIME + SMAX*2
         END DO
!     $omp end parallel do
*
*     Predict particles using new velocities and positions.
!     $omp parallel do schedule(runtime)
         DO J = IFIRST,NTOT
            CALL GPUIRR_SET_JP(J,X0(1,J),X0DOT(1,J),
     &           FI(1,J),FIDOT(1,J),
     &           BODY(J),T0(J))
            TPRED(J) = -1.0
         END DO
!     $omp end parallel do
         CALL STOPWATCH(TEND)
         TINTRE = TINTRE + TEND - TBEG
      END IF
*
*       See whether to advance ARchain or KS at first new time.
      IF (TIME.GT.TPREV) THEN
         CALL SUBINT(IQ,I10)
*       Check collision/coalescence indicator.
          IF (IQ.LT.0) GO TO 999
       END IF
*
*     Form lists of candidates for new irregular and regular force.
      IF (TIME.GE.TBLCKR + SMAX*4) THEN
         DO L = 1,NXTLEN
            IRR(L) = NXTLST(L)
         END DO
         NFR = NXTLEN
      ELSE
         DO L = 1,NXTLEN
            IRR(L) = 0
         END DO
         NFR = 0
      END IF
c      NFR = 0
c      DO 28 L = 1,NXTLEN
c          J = NXTLST(L)
c          IF (TNEW(J).GE.T0R(J) + SMAX*4) THEN
c              NFR = NFR + 1
c              IRR(L) = J
c          ELSE
c              IRR(L) = 0
c          END IF
c 28    CONTINUE
*
      CALL GPUIRR_PRED_ALL(IFIRST,NTOT,TIME)
      NNPRED = NNPRED + 1
*
*       Save new time (output time at TIME > TADJ) and increase # blocks.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
      TMIN = 1.0D+10
      IKS0 = IKS
*
*       Predict chain c.m. only at new block-time.
      IF (NCH.GT.0) THEN
          CALL JPRED(ICH)
*         CALL XCPRED(2)
      END IF
*
*       Evaluate all irregular forces & derivatives in the GPUIRR library.
*     DO 46 II = 1,NXTLEN
*         I = NXTLST(II)
*         CALL GPUIRR_FIRR(I,GF(1,II),GFD(1,II))
*     46 CONTINUE
c      DO II = 1,NXTLEN
c         I = NXTLST(II)
c         IF(RS(I).NE.RS(IFIRST)) THEN
c            WRITE(*,*) 'CHECK RS', I, RS(I)
c         END IF
c      END DO
      CALL GPUIRR_FIRR_VEC(NXTLEN,NXTLST,GF,GFD,RS)
*
*
*       Choose between standard and parallel irregular integration.
      CALL STOPWATCH(TBEG)
      IF (NXTLEN.LE.NPMAX) THEN
*
*       Correct the irregular steps sequentially.
          DO 48 II = 1,NXTLEN
              I = NXTLST(II)
*       Advance the irregular step (no copy of X0 to GPUIRR needed here).
              CALL NBINT(I,IKS,IRR(II),GF(1,II),GFD(1,II))
*             CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),FI(1,I),FIDOT(1,I),
*    &                                                BODY(I),T0(I))
*
*       Save indices and TIME of first KS candidates in the block.
              IF (IKS0.EQ.0.AND.IKS.GT.0) THEN
                  ISAVE = ICOMP
                  JSAVE = JCOMP
                  TSAVE = TIME
                  IKS0 = IKS
              END IF
   48     CONTINUE
*
      ELSE
*       Perform irregular correction in parallel (flag for perturbed binary).
!$omp parallel do private(I)
          DO 50 II = 1,NXTLEN
*       Initialize array for repeated calls to ensure thread-safe procedure.
              ISTAT(II) = -1
              I = NXTLST(II)
              CALL NBINTP(I,IRR(II),GF(1,II),GFD(1,II),ISTAT(II),
     &                    XNEW(1,II),XDNEW(1,II))
*       Note that rejected parallel part contains few operations.
   50     CONTINUE
!$omp end parallel do
          IPHASE = 0
          IKS = 0
*       Examine all block members for irregular steps or X0 & X0DOT copy.
          DO 500 II = 1,NXTLEN
              I = NXTLST(II)
              IF (ISTAT(II).GT.0) THEN
*       Correct exceptional irregular steps in serial.
                  CALL NBINT(I,IKS,IRR(II),GF(1,II),GFD(1,II))
              END IF
 500       CONTINUE
!     $omp parallel do private(I,K)
           DO II = 1,NXTLEN
              I = NXTLST(II)
              IF(ISTAT(II).LE.0) THEN
*       Copy new position and velocity.
                 T0(I) = TIME
                 DO K=1,3
                    X0(K,I) = XNEW(K,II)
                    X0DOT(K,I) = XDNEW(K,II)
                 END DO
              END IF
           END DO
!     $omp end parallel do
        END IF
      CALL STOPWATCH(TEND)
      TINTIR = TINTIR + TEND - TBEG
      RNSTEPI = RNSTEPI + NXTLEN
*
*
*       Check regular force updates (NFR members on block-step #NBLCKR).
      IF (NFR.GT.0) THEN
*
      NBLCKR = NBLCKR + 1
      NN = NTOT - IFIRST + 1
      TBLCKR = TIME
*       Obtain irregular & regular force and determine current neighbours.
c     RSA = RSNEXT
      RSA = RSMIN
      NOFL(1) = 0
 512  NOFL2 = 0
*
      CALL STOPWATCH(TBEG)
!     $omp parallel do
      DO I = IFIRST,NTOT
         RCUT(I) = RSA+RBUFF
      END DO
!     $omp end parallel do
      CALL STOPWATCH(TEND)
      TINTRE = TINTRE + TEND - TBEG
*
      CALL GPUNB_REGF(NN,BODY(IFIRST),
     &     X0(1,IFIRST),X0DOT(1,IFIRST),
     &     GPUACC(1,IFIRST), GPUJRK(1,IFIRST),
     &     LMAX, LIST(1,IFIRST), RCUT(IFIRST))
*
      CALL STOPWATCH(TBEG)
      DO I = IFIRST,NTOT
         NNB = LIST(1, I)
         IF (NNB.LT.0) THEN
            RI2 = (X(1,I)-RDENS(1))**2 + (X(2,I)-RDENS(2))**2 +
     &           (X(3,I)-RDENS(3))**2
            WRITE (41,556)  NSTEPR, NAME(I), LIST(1,I), NNB,
     &           RS(I), SQRT(RI2)
 556        FORMAT (' OVERFLOW!   #R NAME NB0 NB RS ri ',
     &           I11,I6,2I5,2F8.3)
            CALL FLUSH(41)
            NOFL2 = NOFL2 + 1
         END IF
      END DO
*
      IF (NOFL2.GT.0) THEN
         RSA = 0.9 * RSA
         NOFL(1) = NOFL(1) + NOFL2
         GO TO 512
      END IF
*
!     $omp parallel do
      DO I = IFIRST,NTOT
         CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),
     &        FI(1,I),FIDOT(1,I),
     &        BODY(I),T0(I))
      END DO
!     $omp end parallel do
*
!     $omp parallel do private(NNB,LP,L,ITEMP) schedule(guided)
      DO I = IFIRST,NTOT
         NNB = LIST(1,I)
         LP = 1
         DO L = 2,NNB+1
            ITEMP = LIST(L,I) + IFIRST
            IF (ITEMP.NE.I) THEN
               LP = LP + 1
               LIST(LP,I) = ITEMP
            END IF
         END DO
         LIST(1,I) = LP - 1
         CALL GPUIRR_SET_LIST(I, LIST(1,I))
         NXTLST(I - IFIRST + 1) = I
      END DO
      NXTLEN = NTOT - IFIRST + 1
!     $omp end parallel do
      CALL STOPWATCH(TEND)
      TINTRE = TINTRE + TEND - TBEG
*
c!     $omp parallel do
c      DO I = IFIRST,NTOT
c         CALL GPUIRR_FIRR(I,GF(1,I),GFD(1,I),RSA)
c      END DO
c     !     $omp end parallel do
      CALL GPUIRR_FIRR_VEC(NXTLEN,NXTLST,GF,GFD,RS)
      CALL STOPWATCH(TBEG)
!     $omp parallel do
      DO I = IFIRST,NTOT
         TPRED(I) = -1
         X(1,I) = X0(1,I)
         X(2,I) = X0(2,I)
         X(3,I) = X0(3,I)
         XDOT(1,I) = X0DOT(1,I)
         XDOT(2,I) = X0DOT(2,I)
         XDOT(3,I) = X0DOT(3,I)
         RS(I) = RSA
      END DO
!     $omp end parallel do
c!     $omp parallel do
      DO I = IFIRST,NTOT
         CALL GPUINT(I,GF(1,I),GFD(1,I),GPUACC(1,I),GPUJRK(1,I))
      END DO
c!     $omp end parallel do
!     $omp parallel do
      DO I = IFIRST,NTOT
         CALL GPUINT2(I,GF(1,I),GFD(1,I))
      END DO
!     $omp end parallel do
      CALL STOPWATCH(TEND)
      TINTRE = TINTRE + TEND - TBEG
      RSMIN = RSA
*
*       Accumulate the sum of overflows (NOFL(1) holds current number).
      NOFL(2) = NOFL(2) + NOFL(1)
      RNSTEPR = RNSTEPR + NFR
*
      END IF
*
      CALL STOPWATCH(TBEG)
*     Determine next block time (note STEP may shrink in GPUCOR).
!     $omp parallel do reduction(min:TMIN)
      DO 350 L = 1,NXTLEN
         I = NXTLST(L)
         TI = TNEW(I)
         TMIN = MIN(TMIN, TI)
 350  CONTINUE
!     $omp end parallel do
*
*       Enlarge search to LISTQ members in extremely rare case of NXTLEN = 1.
      IF (NXTLEN.EQ.1.OR.((TMIN-TIME).GT.(TIME-TMIN0))) THEN
!     $omp parallel do reduction(min:TMIN)
          DO 355 L = 1,NQ      ! Note this actually happend in 2 runs (08/16).
             J = LISTQ(L)
             TJ = TNEW(J)
             TMIN = MIN(TJ, TMIN)
 355      CONTINUE
!     $omp end parallel do
      END IF
*
*       Copy current coordinates & velocities from corrected values, IQ.NE.0.
      IF (IKS.GT.0.OR.IQ.NE.0.OR.TIME.GE.TADJ.OR.TIME.GE.TNEXT.OR.
     &   (KZ(19).GE.3.AND.STEPXCOMM(TIME,STEPX))) THEN
*       Note: need copy X & XDOT at output (skip in XVPRED; also cf. #31 >0).
!$omp parallel do private(I, L, K)
      DO 360 L = 1,NXTLEN
          I = NXTLST(L)
          DO 58 K = 1,3
              X(K,I) = X0(K,I)
              XDOT(K,I) = X0DOT(K,I)
   58     CONTINUE
  360 CONTINUE
!$omp end parallel do
      END IF
*
*       Send corrected active particles to GPUIRR library.
!$omp parallel do private(I, L)
      DO 60 L = 1,NXTLEN
          I = NXTLST(L)
          CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),
     &         FI(1,I),FIDOT(1,I),
     &         BODY(I),T0(I))
   60 CONTINUE
!$omp end parallel do
*
      CALL STOPWATCH(TEND)
      TINTIR = TINTIR + TEND - TBEG
*     TT2 = DBLTIM()
*     CPRED2 = CPRED2 + (TT2 - TT1)
*     IF (TIME.EQ.TCRIT) WRITE (6,380)  CPRED2, TT2, CPRED, CPRED3
* 380 FORMAT (' TIMING     CPRED2 TT2 CPRED CPRED3 ',1P,4E12.2)
*
*       Check integration of tidal tail members.
      IF (NTAIL.GT.0) THEN
*       Allow large quantized interval with internal iteration.
          IF (DMOD(TIME,0.25D0).EQ.0.0D0) THEN
              DO 65 J = ITAIL0,NTTOT
                  IF (TNEW(J).LE.TIME) THEN
                      CALL NTINT(J)
                  END IF
   65         CONTINUE
          END IF
      END IF
*
*       Exit on KS termination, new multiple regularization or merger.
      IF (IQ.NE.0) THEN
          NBPREV = 0
          IF (IQ.GE.4.AND.IQ.NE.7) THEN
              CALL DELAY(IQ,-1)
          ELSE
*       Ensure correct KS index (KSPAIR may denote second termination).
              KSPAIR = KVEC(I10)
              IPHASE = IQ
          END IF
          GO TO 100
      END IF
*
*       Perform optional check on high-velocity particles at major times.
c      IF (KZ(37).GT.0.AND.LISTV(1).GT.0) THEN
c          IF (DMOD(TIME,STEPM).EQ.0.0D0) THEN
c              CALL SHRINK(TMIN)
c              IF (LISTV(1).GT.0) THEN
c                  CALL HIVEL(-1)
c              END IF
c          END IF
c      END IF
*
*       Check optional mass loss time at end of block-step.
      IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 100-year step (new polynomials).
          IF (TIME.GT.TMDOT.AND.STEPXCOMM(TIME,STEPX)) THEN
              IF (KZ(19).GE.3) THEN
                  CALL MDOT
               ELSE
                  CALL MLOSS
              END IF
              IF (IPHASE.LT.0) GO TO 999
          END IF
      END IF
*
*       Include optional plotting of single BH orbits (1 or 2).
      IF (KZ(45).GT.0.AND.TBH.LT.TIME+TOFF) THEN
          CALL BHPLOT
*       Update time interval (try 5 points per time unit).
          TBH = TIME + TOFF + 2.0D-01
      END IF
*
*       Advance counters and check timer & optional COMMON save (NSUB = 0).
      NTIMER = NTIMER + NXTLEN
      IF (NTIMER.LT.NMAX) GO TO 1

      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.10*N.AND.N.GT.1000.AND.NSUB.EQ.0) THEN
          NSTEPS = 0
          IF (KZ(1).EQ.3) CALL MYDUMP(1,1)
      END IF
*
*       Check option for general binary search.
*     IF (KZ(4).GT.0.AND.TIME - TLASTS.GT.DELTAS) THEN
*         CALL EVOLVE(0,0)
*     END IF
*
*       Include facility for termination of run (create dummy file STOP).
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          IF (NSUB.EQ.0)  WRITE (6,70)
   70     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 1
*
*       Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*       Specify zero step to enforce termination.
          DO 75 L = 1,NSUB
              STEPS(L) = 0.0D0
   75     CONTINUE
          NTIMER = NMAX
          GO TO 1
      END IF
*
*       Obtain elapsed wall-clock time (hours, minutes & seconds).
      CALL WTIME(IHOUR,IMIN,ISEC)
      SECS = 3600.0*IHOUR + 60.0*IMIN + ISEC
      WTOT = WTOT + SECS - WTOT0
      WTOT0 = SECS
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          WT = WTOT/3600.0
          CALL MYDUMP(1,1)
          WRITE (6,80)  TIME+TOFF, TCOMP, CPUTOT/60.0, ERRTOT, DETOT, WT
   80     FORMAT (/,9X,'COMMON SAVED AT TIME =',F8.2,'  TCOMP =',F7.1,
     &                 '  CPUTOT =',F6.1,'  ERRTOT =',F10.6,
     &                 '  DETOT =',F10.6,'  WTOT =',F7.1)
      END IF
*
*       Close GPU & GPUIRR library and stop.
      CALL GPUNB_CLOSE
      CALL GPUIRR_CLOSE
      STOP
*
  100 RETURN
*
      END
*
      LOGICAL FUNCTION STEPXCOMM(T,X)
      REAL*8 T,X
      IF (X.EQ.0.0D0) THEN
         STEPXCOMM = .FALSE.
      ELSE
         STEPXCOMM = DMOD(T,X).EQ.0.0D0
      END IF
      END
