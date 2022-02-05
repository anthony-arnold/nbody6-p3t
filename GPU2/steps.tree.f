      SUBROUTINE STEPS(I1,I2)
*
*
*     Initialization of time-steps & prediction variables.
*     ----------------------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tree.h'
      REAL*8  FDUM(3)
*
*
*     Set new steps and initialize prediction variables.
      DO 40 I = I1,I2

         ZFI = FI(1,I)**2 + FI(2,I)**2 + FI(3,I)**2
         ZFID = D1(1,I)**2 + D1(2,I)**2 + D1(3,I)**2
         IF (ZFI.EQ.0.0D0.OR.ZFID.EQ.0.0D0) THEN
            DT = SMAX
         ELSE
*
*     Include regular force in the irregular step (cf. Makino & SJA 1992).
            DO 5 K = 1,3
               FDUM(K) = FI(K,I) + FR(K,I)
 5          CONTINUE
*
*     Determine irregular and regular steps by the general criterion.
            DT = TSTEP(FDUM,D1(1,I),D2(1,I),D3(1,I),ETAI,
     &           BODY(I), RS(I)+RBUFF)
         END IF
         DTR = SMAX * 4
*
*     Reduce irregular step for case of triple, quad, chain or merger.
         IF (IPHASE.GE.4) DT = 0.5*DT
*
*     Initialize the times and obtain discrete steps (block-step version).
         T0(I) = TIME
         T0R(I) = TIME
*
*     Convert predicted step to nearest block time-step (truncated down).
         CALL STEPK(DT,DTN)
         DTRN = DTR
         IF (TIME.LE.0.0D0) THEN
            STEP(I) = DTN
            STEPR(I) = DTRN
         ELSE
*     Reduce steps by factor 2 until commensurate with current time.
            STEP(I) = DTN
            STEPR(I) = DTRN
            ITER = 0
 10         IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
               STEP(I) = 0.5D0*STEP(I)
               ITER = ITER + 1
               IF (ITER.LT.16.OR.STEP(I).GT.DTK(40)) GO TO 10
               STEP(I) = DTK(40)
               WRITE (6,15)  I, ITER, TIME/STEP(I), DT, STEP(I)
 15            FORMAT (' WARNING!   I ITER T/STEP DT STEP ',
     &              I5,I4,F16.4,1P,2E9.1)
            END IF
         END IF
*
*     Reduce irregular step if STEPR < STEP.
 25      IF (STEPR(I).LT.STEP(I)) THEN
            STEP(I) = 0.5D0*STEP(I)
            GO TO 25
         END IF
*
*     Set maximum irregular step if no neighbours.
         IF (LIST(1,I).EQ.0) THEN
            DT = SMAX*4
         END IF
*
*     Initialize or update array for new block times.
         TNEW(I) = STEP(I) + T0(I)
*     IF (TIME.GT.0.0) THEN
*     WRITE (7,28)  I, NAME(I), TIME, DT, STEP(I), STEPR(I)
*     28          FORMAT (' STEPS:   I NM TIME DT STEP DTR ',
*     &                           2I5,F12.6,1P,3E9.1)
*     CALL FLUSH(7)
*     END IF
*
*     Set prediction variables (X0DOT set by START, KSREG or KSTERM).
         DO 30 K = 1,3
            X0(K,I) = X(K,I)
            F(K,I) = 0.5D0*F(K,I)
            FDOT(K,I) = ONE6*FDOT(K,I)
 30      CONTINUE
 40   CONTINUE
*
      RETURN
*
      END
