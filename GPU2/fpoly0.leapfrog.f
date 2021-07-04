      SUBROUTINE FPOLY0(RS0)
*
*
*       GPU initialization of forces and neighbour lists.
*       -------------------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tree.h'
      REAL*8  RCUT(NMAX)
      SAVE RCUT
*       Open the GPU libraries on new run (note nnbmax = NN is printed).
      NN = NTOT + 100
      CALL GPUNB_OPEN(NN, NNBMAX)
*
*       Set provisional neighbour radius and regular time-step for GPUNB.
      NFI = NTOT - IFIRST + 1
!     $omp parallel do
      DO 10 I = IFIRST,NTOT
         T0(I) = 0.0
         RS(I) = 0.0
*       Initialize F & FDOT.
          DO 5 K = 1,3
              F(K,I) = 0.0
              FDOT(K,I) = 0.0
    5     CONTINUE
   10 CONTINUE
!     $omp end parallel do
*       Define maximum GPU neighbour number and initialize counters.
      RSMIN = 0.0D0
 59   NBMAX = MIN(NNBMAX + 150,LMAX-5)
      NOFL2 = 0
!     $omp parallel do
      DO I = IFIRST,NTOT
         RCUT(I) = 0.0D0
      END DO
!     $omp end parallel do
*
      CALL GPUNB_REGF(NFI,BODY(IFIRST),
     &     X(1,IFIRST),XDOT(1,IFIRST),
     &     FR(1,IFIRST), D1R(1,IFIRST),
     &     LMAX, LIST(1,IFIRST), RCUT(IFIRST))
*
      DO I = IFIRST,NTOT
         NNB = LIST(1,I)
         IF (NNB.LT.0) THEN
            RI2 = (X(1,I)-RDENS(1))**2 + (X(2,I)-RDENS(2))**2 +
     &           (X(3,I)-RDENS(3))**2
            WRITE (41,40)  NSTEPR, NAME(I), NNB,
     &           RS(I), SQRT(RI2)
 40         FORMAT (' OVERFLOW!   #R NAME NB RS ri ',
     &           I11,I7,I5,2F8.3)
            CALL FLUSH(41)
            STOP
         END IF
      END DO
*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNLD(IFIRST,NTOT,1)
      END IF
*
*       Form total force & force derivative and extra variables for XVPRED.
!$omp parallel do private(I)
      DO 110 I = IFIRST,NTOT
         RS(I) = RSMIN
         DO 105 K = 1,3
*
            D0(K,I) = 0.0D0
            D0R(K,I) = FR(K,I)
            FIDOT(K,I) = 0.0D0
            FRDOT(K,I) = D1R(K,I)
*
            F(K,I) = FR(K,I) + FI(K,I)
            FDOT(K,I) = FRDOT(K,I) + FIDOT(K,I)
  105     CONTINUE
  110 CONTINUE
!     $omp end parallel do
*
*       Close the GPU libraries (limits change in INTGRT).
      CALL GPUNB_CLOSE
*
      RETURN
*
      END
