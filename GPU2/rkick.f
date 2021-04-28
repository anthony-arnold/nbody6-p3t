      SUBROUTINE RKICK(I)
*
*
*       Leapfrog velocity kick.
*       -----------------------
*
      INCLUDE 'common6.h'
*
*     Kick velocity.
      DTR2 = SMAX * 2
      DO K = 1,3
         X0DOT(K,I) = X0DOT(K,I) + DTR2*FR(K,I)
      END DO
*
      RETURN
*
      END
