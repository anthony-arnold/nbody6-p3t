      SUBROUTINE GCINT
*
*
*       Integration of 3D cluster orbit.
*       --------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,OMEGA,
     &       MBULGE,RBULGE,MDISK,SCALEA,SCALEB,MHALO,RHALO,VCIRC
      REAL*8  FM(3),FD(3),FS(3),FSD(3)
      REAL*8 DTGC
      parameter (DTGC=0.000244140625)  ! 1/4096
*
*
*       Predict coordinates and velocities to order FDOT.
      If (KZ(14).EQ.5 .and. tg.eq.0.d0) then
         call forceir13 (rg, vg, fg, fgd)
      end if
      IF (KZ(14).EQ.6 .and. tg.eq.0.d0) THEN
         call forcebo15 (rg, vg, fg, fgd)
      END IF

      DT = DTGC

*       Note: integration step may exceed STEPX (depends on block-step).
 100  DT2 = 0.5*DT
      DT3 = ONE3*DT
      RG2 = 0.0
      RGVG = 0.0

      DO 10 K = 1,3
          RG(K) = ((FGD(K)*DT3 + FG(K))*DT2 + VG(K))*DT + RG(K)
          VG(K) = (FGD(K)*DT2 + FG(K))*DT + VG(K)
          RG2 = RG2 + RG(K)**2
          RGVG = RGVG + RG(K)*VG(K)
   10 CONTINUE
*
*       Obtain force and first derivative of point-mass galaxy.
      IF (KZ(14).EQ.3.AND.GMG.GT.0.0D0) THEN
          CALL FNUC(RG,VG,FM,FD)
      ELSE
          DO 15 K = 1,3
              FM(K) = 0.0
              FD(K) = 0.0
   15     CONTINUE
      END IF
*
      if (kz(14).eq.999) then
*       Check bulge force.
      IF (GMB.GT.0.0D0) THEN
          CALL FBULGE(RG,VG,FS,FSD)
          DO 20 K = 1,3
              FM(K) = FM(K) + FS(K)
              FD(K) = FD(K) + FSD(K)
   20     CONTINUE
      END IF
*
*       Include optional Miyamoto disk component (Book eq. 8.52).
      IF (DISK.GT.0.0D0) THEN
          CALL FDISK(RG,VG,FS,FSD)
          DO 25 K = 1,3
              FM(K) = FM(K) + FS(K)
              FD(K) = FD(K) + FSD(K)
   25     CONTINUE
      END IF
*
*       Check addition of logarithmic galaxy potential (not linearized).
      IF (V02.GT.0.0D0) THEN
          CALL FHALO(RG,VG,FS,FSD)
          DO 30 K = 1,3
              FM(K) = FM(K) + FS(K)
              FD(K) = FD(K) + FSD(K)
   30     CONTINUE
      END IF
      end if

      IF (KZ(14).EQ.5) THEN
          call forceir13 (rg, vg, fm, fd)
       END IF
       IF (KZ(14).EQ.6) THEN
          call forcebo15 (rg, vg, fm, fd)
       END IF

*
*       Set time factors for corrector.
      DT13 = ONE3*DT
      DTSQ12 = ONE12*DT**2
      TG = TG + DT
*
*       Include the Hermite corrector and update F & FDOT.
      DO 40 K = 1,3
          DF = FG(K) - FM(K)
          SUM = FGD(K) + FD(K)
          AT3 = 2.0*DF + SUM*DT
          BT2 = -3.0*DF - (SUM + FGD(K))*DT
          RG(K) = RG(K) + (0.6*AT3 + BT2)*DTSQ12
          VG(K) = VG(K) + (0.75*AT3 + BT2)*DT13
          FG(K) = FM(K)
          FGD(K) = FD(K)
   40 CONTINUE
*
*       Update angular velocity in case of non-circular orbit.
*     OM1 = (RG(2)*VG(3) - RG(3)*VG(2))/RG2
*     OM2 = (RG(3)*VG(1) - RG(1)*VG(3))/RG2
      OM3 = (RG(1)*VG(2) - RG(2)*VG(1))/RG2
*     OMEGA2 = OM1**2 + OM2**2 + OM3**2
      OMEGA2 = OM3**2
      OMEGA = SQRT(OMEGA2)
      TIDAL(4) = 2.0*OMEGA

      if (tg.lt.time+toff) goto 100

c Check for end condition
      IF (KZ(14).GE.5.and.(time+toff)*tstar-tcrit.gt.-30.d0) THEN
          call checkend (rg, vg)
      END IF
*
      RETURN
*
      END
