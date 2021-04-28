      SUBROUTINE NBLIST(I,RS0)
*
*
*       Neighbour list.
*       ------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'tree.h'
*
*
*       Form square neighbour radius.
      RS2 = (RS0 + RBUFF)**2
      ILIST(1) = 0
      ITER = 0
*
      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
      IF (VI2.GT.0.0) THEN
          DTR = 0.1*RS0/SQRT(VI2)
          VRFAC = -0.1*RS2/DTR
*       Estimated radial velocity factor for outer sphere.
      ELSE
          VRFAC = 0.0
      END IF
*
*       Search all single particles & c.m. bodies.
    1 RCRIT2 = 1.59*RS2
      NNB = 1
      DO 10 J = IFIRST,NTOT
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          IF (RIJ2.GT.RCRIT2) GO TO 10
          IF (RIJ2.GT.RS2) THEN
*       Skip velocity test for primordial binaries.
              IF (TIME.LE.0.0D0.AND.J.LE.2*NBIN0) GO TO 6
              RIJDOT = 0.0
              DO 4 K = 1,3
                  RIJDOT = RIJDOT +
     &                         (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
    4         CONTINUE
*       Accept member if maximum penetration factor is significant.
              IF (RIJDOT.GE.VRFAC) GO TO 10
          ELSE
              IF (J.EQ.I) GO TO 10
          END IF
*       Include both components in case of primordial binaries.
    6     IF (J.LE.2*NBIN0.AND.TIME.EQ.0.0D0) THEN
              JPAIR = KVEC(J)
*       Consider odd and even locations and skip body $I.
              IF (J.EQ.2*JPAIR-1) THEN
                  NNB = NNB + 1
                  ILIST(NNB) = J
                  IF (I.NE.J + 1) THEN
                      NNB = NNB + 1
                      ILIST(NNB) = J + 1
                  END IF
*       Exclude previously saved component.
              ELSE IF (ILIST(NNB).LT.J-1) THEN
                  IF (I.NE.J - 1) THEN
                      NNB = NNB + 1
                      ILIST(NNB) = J - 1
                  END IF
                  NNB = NNB + 1
                  ILIST(NNB) = J
              END IF
          ELSE
              NNB = NNB + 1
              ILIST(NNB) = J
          END IF
   10 CONTINUE
*
*       Set neighbour radius and copy list membership.
      RS(I) = RS0
      LIST(1,I) = NNB - 1
      DO 20 L = 2,NNB
          LIST(L,I) = ILIST(L)
   20 CONTINUE
*
   25 RETURN
*
      END
