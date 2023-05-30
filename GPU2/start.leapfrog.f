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
*       Special evolution code
      if (kz(19).eq.5) then
         zmbar = zmbar*n
         call instar
         call evolve_mf
      end if

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
 30    CONTINUE
*
*       Set optional stellar evolution parameters or define STEPX.
      IF (KZ(19).EQ.5) THEN
        kz(19)=4
        toff = toff/tstar
        DT = 1.0d-04/TSCALE
        CALL STEPK(DT,DTN)
        IF(DTN*TSCALE.LT.100.0) DTN = 2.d0*DTN
        STEPX = DTN
        TMDOT = STEPX
c convert TEV values from Myr to NBODY units and subtract TOFF
         do i=1,n
           TEV(I)  = TEV(I)/TSTAR-TOFF
           TEV0(I) = TEV0(I)/TSTAR-TOFF
           EPOCH(I) = EPOCH(I)-TOFF*TSTAR
         end do
c set gal. cent. time step to toff
         if (KZ(14).gt.0) TG=TOFF
         SPNFAC = ZMBAR*SU**2/(1.0D+06*TSTAR)
      ELSE

*
*     Set optional stellar evolution parameters or define STEPX.
         IF (KZ(19).GT.2) THEN
            CALL INSTAR
         ELSE IF (KZ(14).GT.1) THEN
            DT = 1.0E-03/TSTAR
            CALL STEPK(DT,DTN)
            STEPX = DTN
         END IF
ENDIF
*
*       Initialize optional cloud parameters.
      IF (KZ(13).GT.0) THEN
          CALL CLOUD0
      END IF
*
*       Perform fast initialization on GPU.
      RS0 = RC
      CALL FPOLY0(RS0)
*
*       Check the average neighbour number.
 80   NNB0 = 0
      NNBHI = 0
      DO 85 I = IFIRST,NTOT
         NNB0 = NNB0 + LIST(1,I)
         NNBHI = MAX(NNBHI, LIST(1,I))
   85 CONTINUE
      ZNB = FLOAT(NNB0)/FLOAT(NTOT-NPAIRS)
      IF (ZNB.NE.0.0) THEN
          WRITE (6,90)  ZNB, ZNBMAX
 90       FORMAT (/,12X,'WARNING!   NON-ZERO NEIGHBOUR NUMBERS',
     &         ' <NNB> =', F8.3, ' <MAX> =', F8.3)
      END IF
*
      WRITE (6,91) NNBHI
 91   FORMAT(/,12X,'NNBHI = ' I6)
      RETURN
*
      END
