      SUBROUTINE START
*
*
*       Initialization of data & polynomials.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      EXTERNAL SCALE,MERGE
      PARAMETER  (NS=12)
*
*
*       Initialize global scalars, counters & useful constants.
      CALL ZERO
*
*       Read input parameters.
      CALL INPUT
*
*       Set initial conditions: BODY(I), X(K,I), XDOT(K,I); I=1,N & K=1,3.
      CALL DATA
*
*       Scale initial conditions to new units.
      CALL SCALE
*
*       Set total mass in case routines DATA & SCALE are not used.
      ZMASS = 0.0D0
      DO 10 I = 1,N
          ZMASS = ZMASS + BODY(I)
   10 CONTINUE
*
*       Define mean mass in scaled units and solar mass conversion factor.
      BODYM = ZMASS/FLOAT(N)
      IF (KZ(5).NE.3) THEN
          ZMBAR = ZMBAR/BODYM
      END IF
*
*       Introduce scaling factors DAYS, YRS, SU, RAU, SMU, TSTAR & VSTAR.
      CALL UNITS
*
*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNL0
      END IF
*
*       Check optional scaling to hot system.
      IF (KZ(29).GT.0) THEN
          CALL HOTSYS
      END IF
*
*       Check option for initial binaries.
      IF (KZ(8).EQ.1.OR.KZ(8).GE.3) THEN
          CALL BINPOP
      END IF
*
*       Include stable primordial triples.
      IF (KZ(18).GT.1.AND.KZ(8).GT.0) THEN
          CALL HIPOP
      END IF
*
*       Check optional initialization for tidal two-body capture (suppress).
*     IF (KZ(27).GT.0) THEN
*         CALL INTIDE
*     END IF
*
*       Set sequential name, maximum mass & primary velocity.
      BODY1 = 0.0
      DO 20 I = 1,N
          NAME(I) = I
          BODY1 = MAX(BODY1,BODY(I))
          DO 15 K = 1,3
              X0DOT(K,I) = XDOT(K,I)
   15     CONTINUE
   20 CONTINUE
*
*       Check renaming BH to avoid confusion with primordial binaries.
      IF (KZ(24).GT.0.AND.NBIN0.GT.0) THEN
          DO 22 I = 1,N
              IF (BODY(I).GT.0.5*BODY1) II = I  ! catches the last BH.
   22     CONTINUE
          NAM2 = NAME(II)
          NAM1 = NAME(II-1)
          NAME(II) = 2
          NAME(II-1) = 1
          WRITE (6,25)  II, NAME(II), NAME(II-1), BODY(II)*SMU
   25     FORMAT (' RENAME BH    I NM BODY M ',3I5,1P,4E10.2)
      END IF
*
*       Randomize particle indices (dummy routine for standard code).
      CALL SWAP
*
*       Initialize fixed block steps (40 levels).
      CALL IBLOCK
*
*       Create table of inverse Stumpff coefficients.
      DO 30 I = 1,NS
          SCOEFF(I) = 1.0D0/((I + 1)*(I + 2))
   30 CONTINUE
*
*       Set optional stellar evolution parameters or define STEPX.
      IF (KZ(19).GT.2) THEN
          CALL INSTAR
      ELSE IF (KZ(14).GT.1) THEN
          DT = 1.0E-03/TSTAR
          CALL STEPK(DT,DTN)
          STEPX = DTN
      END IF
*
*       Initialize optional cloud parameters.
      IF (KZ(13).GT.0) THEN
          CALL CLOUD0
      END IF
*
*       Regularize any hard primordial binaries (assume sequential ordering).
      IF (NBIN0.GT.0) THEN
          DO 50 IPAIR = 1,NBIN0
              ICOMP = 2*IPAIR - 1
              JCOMP = 2*IPAIR
              RIJ2 = 0.0
*       Apply standard distance criterion.
              DO 45 K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICOMP) - X(K,JCOMP))**2
   45         CONTINUE
              IF (RIJ2.LT.RMIN**2) THEN
                  CALL KSPREG
              END IF
   50     CONTINUE
      END IF
*
*       Perform fast initialization on GPU.
      RS0 = RC
      CALL FPOLY0(RS0)
*
*       Initialize force polynomials and time-steps (standard way).
!$omp parallel do private(I)
      DO 55 I=IFIRST,NTOT
          CALL FPOLY2(I,I,0)
   55 ENDDO
!$omp end parallel do
*
*     Generate perturber lists for primordial binaries.
      CALL KSPINIT
*
*       Search other binaries
      RMIN2 = RMIN**2
      RMIN22 = 4.0*RMIN2
      I = N
   99 IKS = 0
      IF (STEP(I).LT.DTMIN) THEN
          CALL SEARCH(I,IKS)
          IF (IKS.GT.0) THEN
              IF (ICOMP.GT.JCOMP) THEN
                  ITMP = ICOMP
                  ICOMP = JCOMP
                  JCOMP = ITMP
              END IF
              CALL KSREG
              GO TO 99
          END IF
      END IF
      I = I - 1
      IF (I.GT.IFIRST) GO TO 99
*
*       Include optional regularization of primordial triples.
      IF (KZ(18).GT.1.AND.NHI0.GT.0) THEN
          KSPAIR = 1
*       Note that each KS pair will move to the top of the queue.
   60     ICOMP = 2*KSPAIR - 1
          ICM = KSPAIR + N
          RX2 = 1.0
*       Find index of closest outer component without any assumption.
          DO 70 J = IFIRST,N
              RIJ2 = 0.0
              DO 65 K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICM) - X(K,J))**2
   65         CONTINUE
              IF (RIJ2.LT.RX2) THEN
                  RX2 = RIJ2
                  JCOMP = J
              END IF
   70     CONTINUE
*       Increase KS pointer for rejected separation (include possible exit).
          IF (RX2.GT.RMIN**2) THEN
              KSPAIR = KSPAIR + 1
              IF (KSPAIR.GT.NPAIRS) GO TO 80
              GO TO 60
          END IF
*       Evaluate PCRIT for R0(NPAIRS) in MERGE since IMPACT is bypassed.
          CALL HISTAB(KSPAIR,JCOMP,PMIN,RSTAB)
*       Initialize the triple (satisfies R < RMIN without stability test).
          IPHASE = 6
          CALL MERGE
*       Examine the same ICM (moved up after successful MERGE).
          IF (NMERGE.LT.NHI0) THEN
              GO TO 60
          END IF
      END IF
*
*       Check the average neighbour number.
 80   NNB0 = 0
      NNBHI = 0
      DO 85 I = IFIRST,NTOT
         NNB0 = NNB0 + LIST(1,I)
         NNBHI = MAX(NNBHI, LIST(1,I))
   85 CONTINUE
      ZNB = FLOAT(NNB0)/FLOAT(NTOT-NPAIRS)
      IF (ZNB.LT.0.25*ZNBMAX.OR.ZNB.LT.0.25*SQRT(FLOAT(N))) THEN
          WRITE (6,90)  ZNB, ZNBMAX
 90       FORMAT (/,12X,'WARNING!   SMALL NEIGHBOUR NUMBERS',
     &         ' <NNB> =', F8.3, ' <MAX> =', F8.3)
      END IF
*
      WRITE (6,91) NNBHI
 91   FORMAT(/,12X,'NNBHI = ' I6)
      RETURN
*
      END
