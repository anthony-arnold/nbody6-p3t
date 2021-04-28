      SUBROUTINE CUTOFF(XI,XIDOT,XJ,XJDOT,D0,D1,BODYJ,RCUT)
*
*
*     Cutoff function for neighbour forces.
*     -------------------------------------
*
      REAL*8 XI(3), XIDOT(3), XJ(3), XJDOT(3), D0(3), D1(3),
     &     BODYJ, RCUT
      REAL*8 DX(3), DVX(3)
      REAL*8 RIN, DRICUT, R2, RV, R, RINV, ALPHA, MRINV3
      REAL*8 X, XV, PK, DPK, PM, DPM, FAC, FJK
*
*     Compute the inner and cutoff radii
      RIN = RCUT * 0.1
      DRICUT = 1.0 / (RCUT - RIN)
*
*     Get distance and velocity difference
      R2 = 0
      RV = 0
      DO K = 1,3
         DX(K) = XJ(K) - XI(K)
         DVX(K) = XJDOT(K) - XIDOT(K)
         R2 = R2 + DX(K)**2
         RV = RV + DX(K)*DVX(K)
      END DO
      R = SQRT(R2)
      RINV = 1.0/R
*     Get constants
      ALPHA = -3.0D0 * RINV**2 * RV
      MRINV3 = BODYJ * RINV**3
*     Determine the hard cutoffs
      X = (R - RIN) * DRICUT
      XV = RV * RINV * DRICUT
*     Clamp X to [0,1]
      IF (X.LE.0.0D0) THEN
         X = 0.0D0
         XV = 0.0D0
      ELSE IF (X.GT.1.0D0) THEN
         X = 1.0D0
         XV = 0.0D0
      END IF
*
*     Set up polynomials
      PK = 1.0 - (((-20.0*X + 70.0)*X - 84.0)*X + 35.0)*(X**4)
      DPK = ((((-140.0*X + 420.0)*X - 420.0)*X + 140.0)*(X**3))*XV
      PM = PK * MRINV3
      DPM = DPK * MRINV3
*
*     Accumulate forces
      DO K = 1,3
         FAC = PM * DX(K)
         FJK = PM * DVX(K) + ALPHA * FAC - DPM * DX(K)
         D0(K) = D0(K) + FAC
         D1(K) = D1(K) + FJK
      END DO
*
      RETURN
*
      END
*
