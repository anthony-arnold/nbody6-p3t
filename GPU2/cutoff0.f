      SUBROUTINE CUTOFF0(XI,XJ,D0,BODYJ,RCUT)
*
*
*     Cutoff function for neighbour force without derivative.
*     -------------------------------------
*
      REAL*8 XI(3), XJ(3), D0(3), BODYJ, RCUT
      REAL*8 DX(3)
      REAL*8 RIN, DRICUT, R2, R, RINV, MRINV3
      REAL*8 X, PK, PM, FAC
*
*     Compute the inner and cutoff radii
      RIN = RCUT * 0.1
      DRICUT = 1.0 / (RCUT - RIN)
*
*     Get distance and velocity difference
      R2 = 0
      DO K = 1,3
         DX(K) = XJ(K) - XI(K)
         R2 = R2 + DX(K)**2
      END DO
      R = SQRT(R2)
      RINV = 1.0/R
*     Get constants
      MRINV3 = BODYJ * RINV**3
*     Determine the hard cutoffs
      X = (R - RIN) * DRICUT
*     Clamp X to [0,1]
      X = MAX(X,0.0D0)
      X = MIN(X,1.0D0)
*
*     Set up polynomials
      PK = 1.0 - (((-20.0*X + 70.0)*X - 84.0)*X + 35.0)*(X**4)
      PM = PK * MRINV3
*
*     Accumulate forces
      DO K = 1,3
         FAC = PM * DX(K)
         D0(K) = D0(K) + FAC
      END DO
*
      RETURN
*
      END
*
