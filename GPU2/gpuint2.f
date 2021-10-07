      SUBROUTINE GPUINT2(I,FIRR,FD)
*
*
*     Regular GPU leapfrog force integrator.
*     --------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tree.h'
      REAL*8 FIRR(3),FD(3)
*
      NNB = LIST(1,I)
*
*     Kick velocity again.
      CALL RKICK(I)
*
      T0(I) = TIME
      T0R(I) = TIME
*
      ZFI = FIRR(1)**2 + FIRR(2)**2 + FIRR(3)**2
      ZFID = FD(1)**2 + FD(2)**2 + FD(3)**2
      IF (ZFI.EQ.0.0D0.OR.ZFID.EQ.0.0D0) THEN
         TTMP = SMAX
      ELSE
*
*     Specify new time-step by standard criterion (STEPI version not good).
         TTMP = SQRT(ETAI * ZFI/ZFID)
      END IF
c      DT0 = TTMP
c      IF (I.GT.N.AND.NNB.LT.20) THEN
c         VI2 = X0DOT(1,I)**2+X0DOT(2,I)**2+X0DOT(3,I)**2
c         DT0 = 0.1*RS(I)/SQRT(VI2)
c         TTMP = MIN(TTMP,DT0)
cc         IF (RI2.GT.0.1) TTMP = 0.25*TTMP
c      END IF
*
*     Select discrete value (increased by 2, decreased by 2 or unchanged).
c      IF (TTMP.GT.2.0*STEP(I)) THEN
c         IF (DMOD(TIME,2.0*STEP(I)).EQ.0.0D0) THEN
c              TTMP = MIN(2.0*STEP(I),SMAX)
c          ELSE
c              TTMP = STEP(I)
c          END IF
c      ELSE IF (TTMP.LT.STEP(I)) THEN
c          TTMP = 0.5*STEP(I)
c          IF (TTMP.GT.DT0) THEN
c              TTMP = 0.5*TTMP
c          END IF
c      ELSE
c          TTMP = STEP(I)
c     END IF
c      IPOW = LOG(TTMP) / LOG(2.0D0)
c     TTMP = 2.0D0**(IPOW - 1)
      DT = TTMP
      CALL STEPK(DT, TTMP)
*
*     Set new block step and update next time.
 87   STEP(I) = TTMP
      ITER = 0
 89   IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
         STEP(I) = 0.5D0*STEP(I)
         ITER = ITER + 1
         IF (ITER.LT.16.OR.STEP(I).GT.DTK(40)) GO TO 89
         STEP(I) = DTK(40)
      END IF
      TNEW(I) = STEP(I) + T0(I)
*
      RETURN
*
      END
