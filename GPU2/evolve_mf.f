      SUBROUTINE EVOLVE_MF
      INCLUDE 'common6.h'
      integer wdflag,nsflag
      common /bhtreat/ wdflag,nsflag
      PARAMETER  (NS=12)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN,LUM
      REAL*8 MENV,RENV,ME,RE,K2,K3,TEND
      PARAMETER(K3=0.21D0)
      REAL*8 MLWIND,RL
      EXTERNAL MLWIND,RL
      real*8 age,tphys,mc,m0,m1,rm,mass,dtm,zmet
      integer*4 i,kw,kw0
      real*8 age0,epch0,rad,TM0,MASSC,mcltotnb,mcltotph
      real*8 DMX,DMA,BODYI(NMAX)

      read (5,*) tend

      write (*,'(/a,f7.1,a/)') "Evolving all stars to an end time of ",
     *    tend," Myr"

      tstar = 1.d0

      do i=1,n
        mass     = body(i)*zmbar
        body0(i) = body(i)

        bodyi(i) = body(i)
        spin(i)  = 0.d0

c       tphys    = 0.d0
c       epoch(i) = 0.d0
c       tev0(i)  = 0.d0
        tev(i)   = 0.d0

c       kw  = 1
c       m0  = mass
c       m1  = mass
c       age = 0.d0

c       CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
c       CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
c    &       RM,LUM,KW,MC,RCC,MENV,RENV,K2)

c       kstar(i) = kw

c       CALL TRDOT(I,DTM,M1)
c       tev(i) = dtm

  502   mass = body(i)*zmbar
c       if (tev(i).lt.tend) then
            dt = 1E6*(TEV(i) - TEV0(i))*TSTAR
c       else
c           dt = 1E6*(TEND - TEV0(i))*TSTAR
c       end if

        m0 = body0(i)*zmbar
        m1 = mass
        mc = 0.D0
        kw = kstar(i)

        age = tev0(i)*tstar - epoch(i)
        age0 = age
        epch0 = epoch(i)

        CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
        CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &       RM0,LUM,KW,MC,RCC,MENV,RENV,K2)


        RAD = RM0
        RADC = RCC
        LUMIN = LUM
        TM0 = TM
        MASSC = MC
        MENV = ME
        RENV = RE
        TBGB = TSCLS(1)
        K2STR = K2

        IF (KW.NE.KSTAR(I))THEN
           IF (KW.GE.13) THEN
              WRITE(6,'(a,I7,2I4,F7.2)')"NS/BH Form ",I,KSTAR(I),KW,MASS
          END IF
          KW = KSTAR(I)
          M1 = MASS
        ENDIF

        DMX = MLWIND(KW,LUM,RM0,M1,MC,0.d0,ZMET)
        DMA = 0.D0

        DMS = (DMX - DMA)*DT
        DMR = ABS(DMS/(MASS + 1.0d-10))

        IF (DMR.GT.0.02) THEN
           DT = DT*0.02/DMR
           DMS = 0.02*MASS
        ENDIF

        IF(KSTAR(I).LT.10)THEN
           DML = MAX(MASS - MASSC,1.0D-07)
           IF (DML.LT.DMS)THEN
              DT = (DML/DMS)*DT
              DMS = DML
           ENDIF
        ENDIF

        M0 = BODY0(I)*ZMBAR
        M10 = M0
        M1 = MASS
        MC = MASSC
        KW = KSTAR(I)

        DMS = (DMX - DMA)*DT
        DMR = ABS(DMS/(MASS + 1.0d-10))

        IF (DMR.GT.TINY) THEN
           M1 = M1 - DMS

*       Check rejuvenation of MS, HG or HE star.
c          IF (KW.LE.2.OR.KW.EQ.7) THEN
c              M0 = M1
c              CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
c              IF (KW.EQ.2) THEN
c                 IF(GB(9).LT.MASSC.OR.M10.GT.ZPARS(3))THEN
c                    M0 = M10
c                 ELSE
c                    EPCH0 = TM + (TSCLS(1)-TM)*(AGE0-TM0)/(TBGB-TM0)
c                    EPCH0 = TEV0(I)*TSTAR - EPCH0
c                 ENDIF
c              ELSE
c                 EPCH0 = TEV0(I)*TSTAR - AGE0*TM/TM0
c              ENDIF
c           ENDIF
        END IF

        TEVK = TEV0(I) + DT/(1.0D+06*TSTAR)
        AGE = TEVK*TSTAR - EPCH0
        AGE0 = AGE

        CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
        CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW,MC,RCC,ME,RE,K2)

        IF (KW.NE.KSTAR(I)) THEN
           DMS = MASS - M1
           DMR = ABS(DMS/(MASS + 1.0d-10))
        ENDIF

        DMSUN = DMS
        DM = DMS/ZMBAR
        KW0 = KSTAR(I)

        BODY(I)   = M1/ZMBAR
        BODY0(I)  = M0/ZMBAR
        KSTAR(I)  = KW
        RADIUS(I) = RM

        IF (KZ(19).GT.3.AND.KW0.NE.KW) THEN
           WRITE(6,926)' TYPE   ', TEV(I), I, NAME(I), DMR, KW0, KW,
     &                  M0,M1,RADIUS(I)*SU, EMDOT
 926       FORMAT(' NEW',A8,' TPHYS I NAM DM/M KW0 KW M0 M R EMD ',
     &                         F7.1,2I7,F6.2,2I3,2F8.3,F7.1,F10.5)
        END IF

        EPOCH(I) = AGE0 + EPCH0 - AGE
        TEV(I) = TEVK
        TEV0(I) = TEV(I)

        CALL TRDOT(I,DTM,M1)
c       write(*,*)"in teff: ",(time+toff)*tstar,
c    *       name(i),(tev(i)+toff)*tstar,(tev0(i)+toff)*tstar,epoch(i)
c       if (tev(i).gt.400)    stop

        TEV(I) = TEV(I) + DTM
        if (tev(i).lt.tend) goto 502

        write (*,'(f8.2,i7,i3,3f11.3,2e13.4)') tend,i,kstar(i),
     *      bodyi(i)*zmbar,body(i)*zmbar,tev0(i),LUM,RM
      end do

      toff = tend

      mcltotnb = 0.d0
      mcltotph = 0.d0

      do i=1,n
        mcltotnb=mcltotnb+body(i)
        mcltotph=mcltotph+body(i)*zmbar
      end do

      write (*,*) "Total cluster mass (NBODY): ",mcltotnb
      write (*,*) "Total cluster mass (PHYS):  ",mcltotph

      do i=1,n
         body(i)  = body(i)/mcltotnb    ! *zmbar/(mcltotnb*mcltotph)
         body0(i) = body0(i)/mcltotnb   ! *zmbar/(mcltotnb*mcltotph)
      end do

      zmbar=mcltotph  ! mcltotnb*mcltotph

      return
      end
