      SUBROUTINE MDOT
*
*
*       Mass loss from evolving stars.
*       Updated 6/1/98 by J. Hurley
*       ------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  BM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX),MLIST(NMAX)
      INTEGER JX(2)
      REAL*8 MASS(2),MASSC(2),RAD(2),RADC(2),LUMIN(2)
      REAL*8 AGE0(2),TM0(2),TBGB(2),MENV(2),RENV(2),K2STR(2)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN
      REAL*8 M0,M1,RM,AGE,LUM,MC,RCC,ME,RE,RNEW,LNEW,MCH
      PARAMETER (MCH=1.44d0)
      REAL*8 M10,RXL1,RXL2,Q,ROL(2),RLPERI,RM0
      REAL*8 EPS,ALPHA2,TOL
      PARAMETER(EPS=1.0d-06,ALPHA2=0.09d0,TOL=1.0d-10)
      REAL*8 DTX(2),DMS(2),DMA(2),DMX(2),EPCH0(2),TEVK,DIFF
      REAL*8 DT,DSEP,DJMB,DJGR,DJORB,DJT,DJTI,DTGR,DSPIN,JSPBRU,OSPBRU
      REAL*8 RJ,HJ,TJ2,VORB2,VWIND2,IVSQM,OMV2,DTXMIN,DMXMAX,DELET
      REAL*8 JSPIN(2),OSPIN(2),DJSPIN(2),JORB,OORB,AURSUN
      PARAMETER(AURSUN=214.95d0)
      REAL*8 K2,K3,ACC1,ACC2,BETA,XI
      PARAMETER(K3=0.21d0,ACC1=3.920659d+08,ACC2=1.5d0)
      PARAMETER(BETA=0.125d0,XI=1.d0,MAXBLK=40)
      REAL*8 MLWIND,RL
      EXTERNAL MLWIND,RL
      CHARACTER*8  WHICH1
      LOGICAL ICORR,IKICK
      SAVE IWARN
      DATA IWARN /0/
*
*       Define current time (unit of million years) and update counter.
      TTOT = TIME + TOFF
      TPHYS = TTOT*TSTAR
      TIME0 = TIME
      IBLUE = 0
      TMDOT = 1.0D+10
*
*       Form list of look-up times and determine next TMDOT.
    1 ML = 0
      ML0 = 0
      DO 2 J = 1,NTOT
          IF (TEV(J).LE.TIME) THEN
              ML0 = ML0 + 1
              MLIST(ML0) = J
          ELSE
              TMDOT = MIN(TMDOT,TEV(J))
          END IF
    2 CONTINUE
      IF (ML0.EQ.0) GO TO 100
      IF (ML0.GE.NMAX) THEN
          WRITE (6,3)  ML0, TIME
    3     FORMAT (' DANGER!    MDOT LIMIT   ML T',I5,F8.2)
          STOP
      END IF
*
*       Treat list members sequentially until GO TO 5 at end.
    4 ML = ML + 1
      I = MLIST(ML)
*
*       Avoid the same binary or merger
      IF (TEV(I).GT.TIME) THEN
          IF (ML.LT.ML0) THEN
              GO TO 4
           ELSE
              GO TO 100
          END IF
      END IF
*
*       Restore synchronized TIME in case of RESET, CHRECT & KSPERI.
    5 TIME = TIME0
*
*       Copy relevant mass (standard case or merger ghost member).
*       Also determine if the star is in a standard binary or is the inner
*       component of a hierachy and if so set the mass accretion flag.
      JX(1) = I

      MASS(1) = BODY(I)*ZMBAR
*
      NMDOT = NMDOT + 1
*
      DO 200 K = 1,1
*
         I = JX(K)
*       Set interval since last mass update.
         DTX(K) = 1.0D+06*(TEV(I) - TEV0(I))*TSTAR
*       Set the initial mass and current type.
         M0 = BODY0(I)*ZMBAR
         M1 = MASS(K)
         MC = 0.D0
         KW = KSTAR(I)
*
*       Obtain stellar parameters at previous epoch.
         AGE = TEV0(I)*TSTAR - EPOCH(I)
         AGE = MAX(AGE,0.D0)
         AGE0(K) = AGE
         EPCH0(K) = EPOCH(I)
         CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
         CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM0,LUM,KW,MC,RCC,ME,RE,K2)

         RAD(K) = RM0
         LUMIN(K) = LUM
         TM0(K) = TM
         TBGB(K) = TSCLS(1)
         MASSC(K) = MC
         RADC(K) = RCC
         MENV(K) = ME
         RENV(K) = RE
         K2STR(K) = K2
*
*       Ensure that type change occurs at time TEV.
         IF(KW.NE.KSTAR(I))THEN
            IF (KW.GE.13) THEN
                WRITE (6,190)  I, NAME(I), KSTAR(I), KW, MASS(K)
  190           FORMAT (' NS/BH FORMATION    ',2I7,2I4,F7.2)
            END IF
            KW = KSTAR(I)
            M1 = MASS(K)
         ENDIF
*
*       Evaluate mass loss due to stellar wind.
         RLPERI = 0.D0
         DMX(K) = MLWIND(KW,LUM,RM0,M1,MC,RLPERI,ZMET)
         DMA(K) = 0.D0
*
*       Convert the spin angular momentum of the star to physical units.
         JSPIN(K) = MAX(SPIN(I)*SPNFAC,1.0D-10)
*
*       Evaluate the spin of the star.
         OSPIN(K) = JSPIN(K)/(K2*RM0*RM0*(M1-MC)+K3*RCC*RCC*MC)
*
 200  CONTINUE
*
      DT = DTX(1)
*
      DO 210 K = 1,1
*
         I = JX(K)
         KW = KSTAR(I)
*
*       Calculate the change in spin angular momentum due to mass loss and
*       also due to mass accretion and/or tides if the star is in a binary.
         DJSPIN(K) = (2.D0/3.D0)*DMX(K)*RAD(K)*RAD(K)*OSPIN(K)
*       Restrict the time-step for extreme conditions (JH 05/09).
         IF(DJSPIN(K).GT.0.0D0)THEN
            DT = MIN(DT,0.3D0*JSPIN(K)/DJSPIN(K))
         ENDIF
*
*       Include magnetic braking for stars with convective envelopes.
         CALL MAGBRK(KW,MASS(K),MENV(K),RAD(K),OSPIN(K),DJMB)
*       Limit to a 3% angular momentum change for the star owing to MB.
         IF(DJMB.GT.0.0D0)THEN
            DT = MIN(DT,0.03D0*JSPIN(K)/DJMB)
            DJSPIN(K) = DJSPIN(K) + DJMB
         ENDIF
*
*       Evaluate the mass loss from the star in the interval DT.
         DMS(K) = (DMX(K) - DMA(K))*DT
         DMR = ABS(DMS(K)/(MASS(K) + 1.0d-10))
*
*       Restrict accumulated mass loss to maximum of 2%.
         IF(DMR.GT.0.02)THEN
            DT = DT*0.02/DMR
            DMS(K) = 0.02*MASS(K)
         ENDIF
*
*       Check that mass loss does not exceed the envelope mass.
         IF(KSTAR(I).LT.10)THEN
            DML = MAX(MASS(K) - MASSC(K),1.0D-07)
            IF(DML.LT.DMS(K))THEN
               DT = (DML/DMS(K))*DT
               DMS(K) = DML
            ENDIF
         ENDIF
         DTX(K) = DT
*
 210  CONTINUE
*
      DO 220 K = 1,1
*       Set the initial mass and current type.
         I = JX(K)
         DT = DTX(K)
         M0 = BODY0(I)*ZMBAR
         M10 = M0
         M1 = MASS(K)
         MC = MASSC(K)
         KW = KSTAR(I)
*       Set indicator for mass loss correction.
         DMS(K) = (DMX(K) - DMA(K))*DT
         DMR = ABS(DMS(K)/(MASS(K) + 1.0d-10))
         IF(DMR.GT.TINY)THEN
            M1 = M1 - DMS(K)
*       Check rejuvenation of MS, HG or HE star.
            IF(KW.LE.2.OR.KW.EQ.7)THEN
               M0 = M1
               CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
               IF(KW.EQ.2)THEN
                  IF(GB(9).LT.MASSC(K).OR.M10.GT.ZPARS(3))THEN
                     M0 = M10
                  ELSE
                     EPCH0(K) = TM + (TSCLS(1)-TM)*(AGE0(K)-TM0(K))/
     &                               (TBGB(K)-TM0(K))
                     EPCH0(K) = TEV0(I)*TSTAR - EPCH0(K)
                  ENDIF
               ELSE
                  EPCH0(K) = TEV0(I)*TSTAR - AGE0(K)*TM/TM0(K)
               ENDIF
            ENDIF
         ENDIF
         DJSPIN(K) = DJSPIN(K)*DT
*
*       Set current age to actual look-up value (allows looping).
         TEVK = TEV0(I) + DT/(1.0D+06*TSTAR)
*       Set indicator for mass loss correction.
         AGE = TEVK*TSTAR - EPCH0(K)
         AGE0(K) = AGE
*
*       Determine stellar evolution time scales and luminosity.
         CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Include temporary error check.
         IF((AGE-TN).GT.1.0d-02)THEN
*           IF(I.GE.IFIRST.OR.(I.LT.IFIRST.AND.NAME(ICM).GT.0))THEN
               WRITE(6,994)NAME(I), KW, DMS(K), AGE, TN
*           ENDIF
 994        FORMAT(' MDOT WARNING! AGE > TN   NM KW DMS AGE TN ',
     &                                   I6,I4,F7.3,1P,2E9.1)
            IF(KW.LE.6)THEN
               AGE = MIN(AGE,0.9999D0*TSCLS(11))
            ELSE
               AGE = MIN(AGE,0.9999D0*TSCLS(5))
               AGE = MIN(AGE,1.0001D0*TN)
            ENDIF
         ENDIF
*##
*       Obtain stellar parameters at current epoch.
         CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW,MC,RCC,ME,RE,K2)
*
*       Check for change of CO WD to ONe on rapid accretion (21/11/08).
         IF(DMR.GT.TINY.AND.KW.EQ.11)THEN
            DME = 2.08D-03*(1.D0/(1.D0 + ZPARS(11)))*RM
            IF(ABS(DMA(K)).GT.0.4d0*DME) KW = 12
         ENDIF
*
*       Check various aspects related to change of type.
         IF(KW.NE.KSTAR(I))THEN
            IF(KW.EQ.15)THEN
               WRITE(6,915)I,NAME(I),M0,KSTAR(I),TEVK*TSTAR-EPCH0(K)
 915           FORMAT (' MDOT WARNING! SN KW=15   I NM M0 KSTAR AGE ',
     &                                            2I6,F7.3,I4,E9.1)
            ELSEIF(KW.GE.13)THEN
*
*       Force new NS or BH to have a period of one second.
               OSPIN(K) = 2.0D+08
               JSPIN(K) = K3*RCC*RCC*MC*OSPIN(K)
            ENDIF
*       Check inclusion of mass loss on type change.
*       Note DM should be zero if ICORR not already true and no SN.
            DMS(K) = MASS(K) - M1
            DMR = ABS(DMS(K)/(MASS(K) + 1.0d-10))
         ENDIF
*
*       Set mass loss and new radius in N-body units.
         DMSUN = DMS(K)
         DM = DMS(K)/ZMBAR
         RNEW = RM/SU
         LNEW = LUM
         KW0 = KSTAR(I)

         BODY(I) = M1/ZMBAR
         BODY0(I) = M0/ZMBAR

*       Accumulate total mass loss (solar units) and reduce cluster mass.
         ZMDOT = ZMDOT + DMSUN
         ZMASS = ZMASS - DM

         IF(DMSUN.LT.TOL) GOTO 250

         VI2 = XDOT(1,I)**2+XDOT(2,I)**2+XDOT(3,I)**2

         POTJ = 0.d0
         NMX=MIN(N,10000)

!     $omp end parallel do private(RIJ) reduction(+:POTJ)
         DO J = 1,NMX
           IF (J.NE.I) THEN
             RIJ = (X(1,I)-X(1,J))**2+(X(2,I)-X(2,J))**2+
     *          (X(3,I)-X(3,J))**2
*
             POTJ = POTJ + BODY(J)/SQRT(RIJ+EPS)
           END IF
         END DO 
!     $omp end parallel do

         POTJ=(POTJ*N)/(1.d0*NMX)

         EMDOT = EMDOT - DM*POTJ + 0.5*DM*VI2
*
*       Update the maximum single body mass but skip compact subsystems.
         IF(MASS(K)/ZMBAR.GE.0.99*BODY1.AND.NSUB.EQ.0)THEN
            BODY1 = 0.d0
            DO 35 J = 1,N
               BODY1 = MAX(BODY1,BODY(J))
  35        CONTINUE
         ENDIF

*       Update the mass loss counters for types > 2.
         IF(KW.EQ.3)THEN
            ZMRG = ZMRG + DMSUN
         ELSEIF(KW.EQ.4)THEN
            ZMHE = ZMHE + DMSUN
         ELSEIF(KW.EQ.5.OR.KW.EQ.6)THEN
            ZMRS = ZMRS + DMSUN
         ELSEIF(KW.GE.7.AND.KW.LE.9)THEN
            ZMNH = ZMNH + DMSUN
         ELSEIF(KW.GE.10.AND.KW.LE.12)THEN
            ZMWD = ZMWD + DMSUN
         ELSEIF(KW.EQ.13.OR.KW.EQ.15)THEN
            ZMSN = ZMSN + DMSUN
         ELSEIF(KW.EQ.14)THEN
            ZMBH = ZMBH + DMSUN
         ENDIF

         IKICK = .FALSE.
         IF (KZ(25).GT.0.AND.KW.GE.10.AND.KW.LE.12) IKICK = .TRUE.
*       Ensure all NS/BH are assigned a kick (might depend on DM).
         IF (KW.EQ.13.OR.KW.EQ.14) IKICK = .TRUE.
         if (IKICK) THEN
              CALL KICK(I,1,KW,DM)
         END IF
*
*       Do not Perform neighbour force corrections on significant mass loss.
*       Ensure that massless supernova remnant will escape next output.
         IF (KW.EQ.15) THEN
            T0(I) = TADJ + DTADJ
            STEP(I) = 1.0D+06
            STEPR(I) = 1.0D+06
            RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
            VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
            DO 55 L = 1,3
*       Ensure that ghost will escape next output (far from fast escapers).
               X0(L,I) = MIN(1.0d+04+X(L,I),1000.0*RSCALE*X(L,I)/RI)
               X(L,I) = X0(L,I)
               X0DOT(L,I) = SQRT(0.004*ZMASS/RSCALE)*XDOT(L,I)/VI
               XDOT(L,I) = X0DOT(L,I)
               F(L,I) = 0.0D0
               FDOT(L,I) = 0.0D0
               D0(L,I) = 0.0
               D1(L,I) = 0.0
               D2(L,I) = 0.0D0
               D3(L,I) = 0.0D0
               D1R(L,I) = 0.0
               D2R(L,I) = 0.0D0
               D3R(L,I) = 0.0D0
   55       CONTINUE
         END IF
*
 250     CONTINUE
*
*       Update event counters for types > 2.
         IF(KW.NE.KW0)THEN
            IF(KW.EQ.3)THEN
               NRG = NRG + 1
            ELSEIF(KW.EQ.4)THEN
               NHE = NHE + 1
            ELSEIF(KW.EQ.5)THEN
               NRS = NRS + 1
            ELSEIF(KW.GE.7.AND.KW.LE.9.AND.KW0.LE.6)THEN
               NNH = NNH + 1
            ELSEIF(KW.GE.10.AND.KW.LE.12.AND.KW0.LE.9)THEN
               NWD = NWD + 1
            ELSEIF((KW.EQ.13.OR.KW.EQ.15).AND.KW0.LE.12)THEN
               NSN = NSN + 1
            ELSEIF(KW.EQ.14.AND.KW0.LE.13)THEN
               NBH = NBH + 1
            ENDIF
         ENDIF
*
*       Include consistency warnings.
         IF (RNEW - RADIUS(I).GT.0.5*RADIUS(I)) THEN
            WRITE(6,924) I, NAME(I), TPHYS, DT/1.0d+06,
     &                    KSTAR(I), KW, M0, M1, RADIUS(I)*SU, RNEW*SU
 924        FORMAT(' EXPAND!    I NM TP DTP K* KW M0 M R RN ',
     &                          2I6,F7.1,F7.3,2I4,2F7.1,2F7.1)
            CALL FLUSH(6)
         END IF
*
*       Update R, L, classification type & spin angular momentum.
         RADIUS(I) = RNEW
         ZLMSTY(I) = LNEW
         KSTAR(I) = KW
         JSPIN(K) = MAX(JSPIN(K) - DJSPIN(K),1.0D-10)
*
*       Ensure that the star does not spin up beyond break-up.
         OSPBRU = TWOPI*SQRT(MASS(K)*AURSUN**3/RAD(K)**3)
         JSPBRU = (K2STR(K)*(MASS(K)-MASSC(K))*RAD(K)*RAD(K) +
     &             K3*MASSC(K)*RADC(K)*RADC(K))*OSPBRU
         IF(JSPIN(K).GT.JSPBRU) JSPIN(K) = JSPBRU
         SPIN(I) = JSPIN(K)/SPNFAC
*
*       Update epoch and check binary diagnostics for transition to new type.
         IF (KW.NE.KW0) THEN
            IF(KW.LT.0)THEN
               WRITE(6,925)I,NAME(I),KW0,KW,M10,M1,DMS(K),RM,AGE
 925           FORMAT(' MDOT CHANGE: I NM K* M0 M1 DM R AGE',
     &                               2I6,2I4,2F6.1,F7.3,F7.2,F8.1)
            ENDIF
         END IF
*
*       Include optional diagnostics (new type or significant mass loss).
         IF (KZ(19).GT.3.AND.KW0.NE.KW) THEN
            WHICH1 = ' TYPE   '
            WRITE(6,926)WHICH1, TPHYS, I, NAME(I), DMR, KW0, KW, M0, M1,
     &                  RADIUS(I)*SU, EMDOT
 926        FORMAT(' NEW',A8,' TPHYS I NAM DM/M KW0 KW M0 M R EMD ',
     &                         F7.1,2I7,F6.2,2I3,2F7.2,F7.1,F10.5)
         END IF
*
*       Base new time scale for changes in radius & mass on stellar type.
         EPOCH(I) = AGE0(K) + EPCH0(K) - AGE
*        TEV(I) = (AGE0(K) + EPCH0(K))/TSTAR
         TEV(I) = TEVK
         TEV0(I) = TEV(I)
         CALL TRDOT(I,DTM,M1)
         TEV(I) = TEV(I) + DTM
*
*       Note any formation of black holes or TZ object.
         IF(KW.EQ.14.AND.KSTAR(I).LT.14.AND.I.LE.N)THEN
            WRITE(6,930) I, NAME(I), KW0, KW, KSTAR(I), M0, M1, DMR
 930        FORMAT(' NEW BH/TZ    I NM K0 KW K* M0 M1 DM/M ',
     &                            2I6,3I4,3F7.2)
         END IF
*
*       Perform consistency check on massive WD (skip HMDOT treatment).
         IF ((KSTAR(I).GE.10.AND.KSTAR(I).LE.12))THEN
             IF (BODY(I)*ZMBAR.GT.MCH)THEN
               WRITE(6,932)I,KW,BODY0(I)*ZMBAR,BODY(I)*ZMBAR,RADIUS(I)
 932           FORMAT(' DANGER!  MDOT   I K* M0 M R ',2I5,2F7.2,1P,E9.1)
               WRITE(6,934) NAME(I), TTOT, TEV0(I), TEV(I)
 934           FORMAT(' NAM T TEV0 TEV ',I6,3F10.2)
               STOP
            ENDIF
         ENDIF

*       Include stellar radius check in case of NaN.
         IF (RADIUS(I).GE.0.0.AND.RADIUS(I).LT.1.0) GO TO 105
         WRITE(6,936) I, KSTAR(I), M1, RADIUS(I)
 936     FORMAT(' DANGER!    MDOT    I K* M1 R* ',I5,I4,F7.2,1P,E10.1)
         STOP
 105     CONTINUE
 220  CONTINUE
*
*       See whether any other stars should be considered.
      IF (TEV(I).GT.TIME)  THEN
          TMDOT = MIN(TMDOT,TEV(I))
      ELSE
          GO TO 5
      END IF

      IF (ML.LT.ML0) GO TO 4
*
 100  RETURN
*
      END
