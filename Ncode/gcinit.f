      SUBROUTINE GCINIT
*
*
*       Initialization of 3D cluster orbit.
*       -----------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,OMEGA,
     &       MBULGE,RBULGE,MDISK,SCALEA,SCALEB,MHALO,RHALO,VCIRC
      REAL*8  FD(3),FDD(3)
*
*
      IF (KZ(14).EQ.5) then
         WRITE(*,*) "Init force Irrgang"
         CALL FORCEIR13(RG,VG,FG,FGD)
         GOTO 999
      END IF

      IF (KZ(14).EQ.6) then
         WRITE(*,*) "Init force Bovy"
         CALL FORCEBO15(RG,VG,FG,FGD)
         GOTO 999
      END IF
*       Obtain initial force and first derivative from point-mass galaxy.
      IF (GMG.GT.0.0D0) THEN
          RIN2 = 1.0/(RG(1)**2 + RG(2)**2 + RG(3)**2)
          RIN3 = RIN2*SQRT(RIN2)
          RGVG = 3.0*(RG(1)*VG(1) + RG(2)*VG(2) + RG(3)*VG(3))*RIN2
*
          DO 1 K = 1,3
              FG(K) = -GMG*RG(K)*RIN3
              FGD(K) = -GMG*(VG(K) - RGVG*RG(K))*RIN3
    1     CONTINUE
      ELSE

          DO 5 K = 1,3
              FG(K) = 0.0
              FGD(K) = 0.0
    5     CONTINUE
      END IF
*
*      Add possible effect of the bulge (for determining tidal radius).
      IF (GMB.GT.0.0D0) THEN
          CALL FBULGE(RG,VG,FD,FDD)
          DO 10 K = 1,3
              FG(K) = FG(K) + FD(K)
              FGD(K) = FGD(K) + FDD(K)
   10     CONTINUE
      END IF
*
*       Include optional disk contribution (may dominate tidal field).
      IF (DISK.GT.0.0D0) THEN
          CALL FDISK(RG,VG,FD,FDD)
          DO 15 K = 1,3
              FG(K) = FG(K) + FD(K)
              FGD(K) = FGD(K) + FDD(K)
   15     CONTINUE
      END IF
*
*       Check addition of logarithmic galaxy potential (not linearized).
      IF (V02.GT.0.0D0) THEN
          CALL FHALO(RG,VG,FD,FDD)
          DO 20 K = 1,3
               FG(K) = FG(K) + FD(K)
               FGD(K) = FGD(K) + FDD(K)
   20     CONTINUE
      END IF
*
*       Initialize the time of GC integration.
 999  TG = TIME + TOFF
*
      RETURN
*
      END
