      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*     -------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tree.h'
      REAL*4 TINTG(2)
      REAL*8 RCUT(NMAX)
      REAL*8 FRX(3),FDX(3)
      SAVE ISTART, RCUT
*
*       Open the GPU libraries on each new run (note nnbmax = NN is printed).
      IF (ISTART.EQ.0) THEN
          NN = NTOT + 100
          CALL GPUNB_OPEN(NN, NNBMAX)
          ISTART = 1
          TBH = TTOT
*
*     Set empty neighbour radius.
!     $omp parallel do
          DO I = IFIRST,NTOT
             RCUT(I) = 0.0D0
          END DO
!     $omp end parallel do
*
          CALL STOPWATCH(TBEG)
          WRITE (6,2) TBEG
 2        FORMAT(//,9X,'START INTGRT ', E11.3)
       END IF
      NN = NTOT - IFIRST + 1
      STEPLF = SMAX
      STEPLF2 = SMAX * 0.5D0
      IPHASE = 0
      TPREV = TIME
!     $omp parallel do schedule(runtime)
      DO J = IFIRST,NTOT
         X0(1,J) = X(1,J)
         X0(2,J) = X(2,J)
         X0(3,J) = X(3,J)
      END DO
!     $omp end parallel do
*
*     Set new time and save block time
 1    TIME = TIME + SMAX
      TTOT = TIME + TOFF
*
*     Half Drift
***   T_pos = T_vel + 0.5
*
      CALL STOPWATCH(TBEG)
!     $omp parallel do schedule(runtime)
      DO J = IFIRST,NTOT
         DO K=1,3
            X0(K,J) = X0(K,J) + X0DOT(K,J)*STEPLF2
         END DO
      END DO
!     $omp end parallel do
!     $omp parallel do schedule(runtime)
      DO J = IFIRST,NTOT
         T0(J) = TIME
      END DO
!     $omp end parallel do
      CALL STOPWATCH(TEND)
      TINTIR = TINTRE + TEND - TBEG
      RNSTEPI = RNSTEPI + NN
*
* Update galaxy guiding centre
      IF (KZ(14).GT.1) THEN
         CALL GCINT
      END IF

*     New forces
      CALL GPUNB_LF(NN,BODY(IFIRST),
     &     X0(1,IFIRST),X0DOT(1,IFIRST),
     &     GPUACC(1,IFIRST), GPUJRK(1,IFIRST),
     &     GPUPHI(IFIRST), RCUT(IFIRST))
*
!     $omp parallel do schedule(runtime)
      DO J = IFIRST,NTOT
         FR(1,J) = GPUACC(1,J)
         FR(2,J) = GPUACC(2,J)
         FR(3,J) = GPUACC(3,J)
      END DO
!     $omp end parallel do
*
      IF (KZ(13).GT.0) THEN
!     $omp parallel do schedule(runtime)
         DO J = IFIRST,NTOT
*     Check optional interstellar clouds.
            CALL FCLOUD(I,FR(1,J),FRDOT(1,J),1)
         END DO
!     $omp end parallel do
      END IF
      IF (KZ(14).GT.0) THEN
      DO J = IFIRST,NTOT
*     Save current values for deriving work done by tides (#14 = 3).
         IF (KZ(14).EQ.3) THEN
            DO 22 K = 1,3
               FRX(K) = FR(K,J)
               FDX(K) = FRDOT(K,J)
 22         CONTINUE
         END IF

         CALL XTRNLF(
     &        X0(1,J), X0DOT(1,J),
     &        FI(1,J), FR(1,J),
     &        FIDOT(1,J), FRDOT(1,J),  2)
*
*     Form rate of tidal energy change during last regular step.
         IF (KZ(14).EQ.3) THEN
            WDOT = 0.0
            W2DOT = 0.0
            DO 24 K = 1,3
               PX = FR(K,J) - FRX(K)
               DPX = FRDOT(K,J) - FDX(K)
               WDOT = WDOT + X0DOT(K,J)*PX
               W2DOT = W2DOT + FR(K,J)*PX + X0DOT(K,J)*DPX
 24         CONTINUE
*     Note: second-order term derived by Douglas Heggie (Aug/03).
         END IF
*
*     Include the force from optional gaseous Plummer potential.
         IF (KZ(14).GE.3) THEN
            CALL XTRNLF(
     &           X0(1,J), X0DOT(1,J),
     &           FI(1,J), FR(1,J),
     &           FIDOT(1,J), FRDOT(1,J), -1)
*
*     Accumulate tidal energy change for general galactic potential.
*     Note: Taylor series at end of interval with negative argument.
            ETIDE = ETIDE + BODY(J)*(0.5*W2DOT*SMAX - WDOT)*SMAX
*     Note: integral of Taylor series for V*P using final values.
         END IF
      END DO
      END IF
*
*     Full-kick
***   T_vel = T_pos + 0.5
*
      CALL STOPWATCH(TBEG)
!     $omp parallel do schedule(runtime)
      DO J = IFIRST,NTOT
         DO K=1,3
            X0DOT(K,J) = X0DOT(K,J) + FR(K,J)*STEPLF
         END DO
         T0R(J) = TIME + SMAX*0.5D0
      END DO
!     $omp end parallel do
      CALL STOPWATCH(TEND)
      TINTRE = TINTRE + TEND - TBEG
*
*
*     Half Drift
      CALL STOPWATCH(TBEG)
!     $omp parallel do schedule(runtime)
      DO J = IFIRST,NTOT
         DO K=1,3
            X0(K,J) = X0(K,J) + X0DOT(K,J)*STEPLF2
         END DO
      END DO
*
*
*** HERE: T_pos = T_vel
*
*       Check next adjust time before beginning a new block.
      IF (TIME.GE.TADJ) THEN
          TIME = TADJ
          IPHASE = 3
          GO TO 100
      END IF
*
*       Advance counters and check timer & optional COMMON save (NSUB = 0).
      NTIMER = NTIMER + NN
      IF (NTIMER.LT.NMAX) GO TO 1
      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.10*N.AND.N.GT.1000.AND.NSUB.EQ.0) THEN
          NSTEPS = 0
          IF (KZ(1).EQ.3) CALL MYDUMP(1,1)
      END IF
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
*       Obtain elapsed wall-clock time (hours, minutes & seconds).
      CALL WTIME(IHOUR,IMIN,ISEC)
      SECS = 3600.0*IHOUR + 60.0*IMIN + ISEC
      WTOT = WTOT + SECS - WTOT0
      WTOT0 = SECS
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
!     $omp parallel do
         DO I = IFIRST,NTOT
            X(1,I) = X0(1,I)
            X(2,I) = X0(2,I)
            X(3,I) = X0(3,I)
            XDOT(1,I) = X0DOT(1,I)
            XDOT(2,I) = X0DOT(2,I)
            XDOT(3,I) = X0DOT(3,I)
         END DO
!     $omp end parallel do
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
      STOP
*       Check optional mass loss time at end of block-step.
 100  IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 100-year step (new polynomials).
          IF (TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              IF (KZ(19).GE.3) THEN
                  CALL MDOT
              ELSE
                  CALL MLOSS
              END IF
          END IF
      END IF
*     COPY for adjust
!     $omp parallel do
 105  DO I = IFIRST,NTOT
         X(1,I) = X0(1,I)
         X(2,I) = X0(2,I)
         X(3,I) = X0(3,I)
         XDOT(1,I) = X0DOT(1,I)
         XDOT(2,I) = X0DOT(2,I)
         XDOT(3,I) = X0DOT(3,I)
      END DO
!     $omp end parallel do
*
 110  RETURN
*
      END
