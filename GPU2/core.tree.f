      SUBROUTINE CORE
*
*
*       DENSITY CENTRE & CORE RADIUS.
*       -----------------------------
*
      INCLUDE 'common6.h'
      PARAMETER(NOPT=6)
      REAL*8 RAD(NOPT,NMAX), DENS(NMAX)

c      NMX = MIN(5000+IFIRST, NTOT) ! NUMBER OF STARS WHICH ARE BEING USED
      NMX = MIN(NTOT, INT(5000**0.5*N**0.5))
*
*     Find 6 nearest neighbours.
      CALL KNN(NMX, NTOT-IFIRST+1, X(1,IFIRST), NOPT, RAD(1,IFIRST))
!     $omp parallel do
      DO I=IFIRST,NMX
         DENS(I) = BODY(I)/SQRT(RAD(NOPT,I))**3.D0
      END DO
!     $omp end parallel do
*
      DO K=1,3
         RDENS(K) = 0.D0
      END DO
*
      RHO1 = 0.0D0
      DO J=IFIRST,NMX
         DO K=1,3
            RDENS(K) = RDENS(K) + X(K,J) * DENS(J)
         END DO
         RHO1 = DENS(J) + RHO1
      END DO
*
      RHO1 = MAX(RHO1,ZMASS/RSCALE**3)
      DO K=1,3
         RDENS(K) = RDENS(K) / RHO1
      END DO
*
      RHOS = 0.0D0
      RHO2 = 0.0D0
      RC = 0.0D0
      DO J=IFIRST,NMX
         RID2 = 0.0D0
         DO K=1,3
            RID = X(K,J) - RDENS(K)
            RID2 = RID2 + RID**2
         END DO
*
         RHO2 = DENS(J)**2 + RHO2
         RC = RC + DENS(J)**2*RID2
         RHOS = MAX(RHOS, DENS(J))
      END DO
*
      IF (RHO2.GT.0.0D0) RC = SQRT(RC/RHO2)
      RHOD = RHO2/RHO1
      RHOM = RHOS
*       Set square core radius and its inverse for routine REGINT.
      RC22 = RC**2
      RC2IN = 1.0/RC22
*
*       Sum particles & mass inside the core radius and set rms velocity.
      NC1 = 0
      VC = 0.D00
      ZMC = 0.D00
      DO J=IFIRST,NTOT
         XID = X(1,J) - RDENS(1)
         YID = X(2,J) - RDENS(2)
         ZID = X(3,J) - RDENS(3)
         RID2 = XID**2 + YID**2 + ZID**2
         IF (RID2.LT.RC22) THEN
            NC1 = NC1 + 1
            VC = VC + XDOT(1,J)**2 + XDOT(2,J)**2 + XDOT(3,J)**2
            ZMC = ZMC + BODY(J)
         END IF
      END DO
*
*       Set core membership & rms velocity.
      NC = MAX(NC1,2)
      VC = SQRT(VC/FLOAT(NC))
*
*
      RETURN
*
      END
