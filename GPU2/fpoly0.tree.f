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
*
      TBLCKR = TIME
*     Read optimal average neighbour number and buffer shell radius
      READ (5,*) NNBOPT, RBUFF
      WRITE (*,*) 'NNBOPT, RBUFF, TBLCKR', NNBOPT, RBUFF, TBLCKR
c     RSNEXT = RS0
*
*
*       Open the GPU libraries on new run (note nnbmax = NN is printed).
      NN = NTOT + 100
      CALL GPUNB_OPEN(NN, NNBMAX)
*
*       Set provisional neighbour radius and regular time-step for GPUNB.
      NFI = NTOT - IFIRST + 1
!     $omp parallel do
      DO 10 I = IFIRST,NTOT
         T0(I) = 0.0
         RS(I) = RS0
*       Initialize F & FDOT.
          DO 5 K = 1,3
              F(K,I) = 0.0
              FDOT(K,I) = 0.0
    5     CONTINUE
   10 CONTINUE
!     $omp end parallel do
*       Set larger value for GPUIRR (note further possible increase of NTOT).
      NNN = NN + 10
      CALL GPUIRR_OPEN(NNN,LMAX)
*
*       Send all particles (X & XDOT) and zero F & FDOT to GPU.
!$omp parallel do private(I)
      DO 15 I = IFIRST,NTOT
          CALL GPUIRR_SET_JP(I,X(1,I),XDOT(1,I),F(1,I),FDOT(1,I),
     &                                          BODY(I),T0(I))
   15 CONTINUE
!$omp end parallel do
*
      CALL GPUIRR_PRED_ALL(IFIRST,NTOT,TIME)  ! may not be needed.
*
*       Send all single particles and c.m. bodies (NFI) to GPU.
c      CALL GPUNB_SEND(NFI,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
*
*       Define maximum GPU neighbour number and initialize counters.
      RSMIN = RS0
 59   NBMAX = MIN(NNBMAX + 150,LMAX-5)
      NOFL2 = 0
!     $omp parallel do
      DO I = IFIRST,NTOT
         RCUT(I) = RSMIN + RBUFF
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
            RSMIN = 0.9 * RSMIN
            NOFL2 = NFOL2 + 1
         END IF
      END DO
*
      IF (NOFL2.GT.0) THEN
         NOFL(1) = NOFL(1) + NOFL2
         GO TO 59
      END IF
*
!     $omp parallel do private(NNB,LP,L,ITEMP)
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
      END DO
!     $omp end parallel do
*
!     $omp parallel do schedule(guided)
      DO I = IFIRST,NTOT
         CALL GPUIRR_FIRR(I,FI(1,I),D1(1,I),RSMIN)
      END DO
!     $omp end parallel do
*
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
            FR(K,I) = FR(K,I) - FI(K,I)
*
            D0(K,I) = FI(K,I)
            D0R(K,I) = FR(K,I)
            FIDOT(K,I) = D1(K,I)
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
      CALL GPUIRR_CLOSE
*
      RETURN
*
      END
