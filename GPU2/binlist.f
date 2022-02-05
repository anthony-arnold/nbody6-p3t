      SUBROUTINE BINLIST (NBINA)
      include 'common6.h'
      real*8 mmin
      integer nbina

      ttot = time + toff

      nbina = 0      ! Reset counter for total binary number

      DO IPAIR = 1,NPAIRS
*       Skip pairs with zero mass of c.m. particle (merged binary ghost).
          IF (BODY(N+IPAIR).GT.0.0D0) THEN
*       Predict coordinates, velocities & binding energy.
            CALL RESOLV(IPAIR,1)

            i1  = 2*ipair-1
            i2  = 2*ipair

            call bparest (i1, i2, nbina, 1)

          END IF
      END DO

c search for further non-regularised binaries

      do i1=ifirst,n
        nnb = list(1,i1)
        dmin = 1.E10
        merk = 0
        do j=1,nnb
          i2 = list(j+1,i1)
          dis=sqrt((x(1,i1)-x(1,i2))**2+(x(2,i1)-x(2,i2))**2+
     *         (x(3,i1)-x(3,i2))**2)
          if (dis.lt.dmin) then
            merk=i2
            dmin=dis
          end if
        end do

        if (merk.eq.0) then
          dmin = 1.E10
          do j=ifirst,ntot
            dis=sqrt((x(1,i1)-x(1,j))**2+(x(2,i1)-x(2,j))**2+
     *         (x(3,i1)-x(3,j))**2)
            if (dis.lt.dmin.and.i1.ne.j) then
              merk=j
              dmin=dis
            end if
          end do
        end if

        i2=merk
        if (i1.lt.i2.and.i2.ne.0) then

        rdis = sqrt((x(1,i1)-x(1,i2))**2+(x(2,i1)-x(2,i2))**2+
     *      (x(3,i1)-x(3,i2))**2)
        vcmx = (body(i1)*xdot(1,i1)+body(i2)*xdot(1,i2))/
     *             (body(i1)+body(i2))
        vcmy = (body(i1)*xdot(2,i1)+body(i2)*xdot(2,i2))/
     *             (body(i1)+body(i2))
        vcmz = (body(i1)*xdot(3,i1)+body(i2)*xdot(3,i2))/
     *             (body(i1)+body(i2))
        eb1  = 0.5d0*body(i1)*((xdot(1,i1)-vcmx)**2+
     *           (xdot(2,i1)-vcmy)**2+(xdot(3,i1)-vcmz)**2)
        eb2  = 0.5d0*body(i2)*((xdot(1,i2)-vcmx)**2+
     *           (xdot(2,i2)-vcmy)**2+(xdot(3,i2)-vcmz)**2)
        ebini= -body(i1)*body(i2)/rdis+eb1+eb2

        rav = 1.d0/(0.5d0*n/(4.d0*3.1415/3.d0*0.8**3))**(1.d0/3.d0)

        if (ebini.lt.0.d0.and.rdis.lt.rav/3.d0) then

          xcmx = (body(i1)*x(1,i1)+body(i2)*x(1,i2))/
     *             (body(i1)+body(i2))
          xcmy = (body(i1)*x(2,i1)+body(i2)*x(2,i2))/
     *             (body(i1)+body(i2))
          xcmz = (body(i1)*x(3,i1)+body(i2)*x(3,i2))/
     *             (body(i1)+body(i2))

          mmin = 1.E10
          do j=ifirst, ntot
            if (j.ne.i1.and.j.ne.i2) then
              dis = sqrt((x(1,j)-xcmx)**2+(x(2,j)-xcmy)**2+
     *               (x(3,j)-xcmz)**2)
              if (dis.lt.mmin) then
                mmin = dis
                merk = j
              end if
            end if
          end do

          fr1 = body(i2)/rdis**2
          fr2 = body(i1)/rdis**2

          fd2 = body(merk)/mmin**3*rdis

          if (fr1.gt.10.0*fd2.and.fr2.gt.10.0*fd2) then
c         write (*,*) 'Bb: ',name(i1),name(i2),name(merk),fr1,
c    *           fr2,fd,fd2
            call bparest (i1, i2, nbina, 2)
          end if
        end if
        end if
      end do

      call flush(38)
      isend = -1     ! Force new send in intgrt.f

      RETURN
      END


      SUBROUTINE BPAREST (i1, i2, nbina, mode)
      include 'common6.h'
      COMMON/BINARY/  ZM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      integer i1,i2,nbina,mode
      character*3 comms
      real*8 xx,lbin,lx,ly,lz,mred,rapo,rper,amax,mges,x1,vper,vapo

      rdis = sqrt((x(1,i1)-x(1,i2))**2+(x(2,i1)-x(2,i2))**2
     *          +(x(3,i1)-x(3,i2))**2)

      vcmx = (body(i1)*xdot(1,i1)+body(i2)*xdot(1,i2))/
     *             (body(i1)+body(i2))
      vcmy = (body(i1)*xdot(2,i1)+body(i2)*xdot(2,i2))/
     *             (body(i1)+body(i2))
      vcmz = (body(i1)*xdot(3,i1)+body(i2)*xdot(3,i2))/
     *             (body(i1)+body(i2))

      eb1  = 0.5d0*body(i1)*((xdot(1,i1)-vcmx)**2+
     *           (xdot(2,i1)-vcmy)**2+(xdot(3,i1)-vcmz)**2)
      eb2  = 0.5d0*body(i2)*((xdot(1,i2)-vcmx)**2+
     *           (xdot(2,i2)-vcmy)**2+(xdot(3,i2)-vcmz)**2)
      ebini= -body(i1)*body(i2)/rdis+eb1+eb2

      if (ebini.lt.0.d0) then

      x1cm=(body(i1)*x(1,i1)+body(i2)*x(1,i2))/(body(i1)+body(i2))  ! CM
      x2cm=(body(i1)*x(2,i1)+body(i2)*x(2,i2))/(body(i1)+body(i2))
      x3cm=(body(i1)*x(3,i1)+body(i2)*x(3,i2))/(body(i1)+body(i2))

      rdc = sqrt((rdens(1)-x1cm)**2+(rdens(2)-x2cm)**2+
     *              (rdens(3)-x3cm)**2)

      mred = body(i1)*body(i2)/(body(i1)+body(i2))        ! Mass of reduced particle
      lx = (x(2,i1)-x(2,i2))*(xdot(3,i1)-xdot(3,i2))-
     *         (x(3,i1)-x(3,i2))*(xdot(2,i1)-xdot(2,i2))
      ly = (x(3,i1)-x(3,i2))*(xdot(1,i1)-xdot(1,i2))-
     *         (x(1,i1)-x(1,i2))*(xdot(3,i1)-xdot(3,i2))
      lz = (x(1,i1)-x(1,i2))*(xdot(2,i1)-xdot(2,i2))-
     *         (x(2,i1)-x(2,i2))*(xdot(1,i1)-xdot(1,i2))
      lbin = sqrt(lx**2+ly**2+lz**2)

      ebin = 0.5*mred*((xdot(1,i1)-xdot(1,i2))**2+
     *    (xdot(2,i1)-xdot(2,i2))**2+(xdot(3,i1)-xdot(3,i2))**2)-
     *     body(i1)*body(i2)/rdis

      xx = body(i1)*body(i2)/(2.d0*ebin)
      xx2 = xx**2+0.5*mred*lbin**2/ebin
      if (xx2.lt.0.0) xx2=0.0
      rapo = -xx+sqrt(xx2)
      rper = -xx-sqrt(xx2)

      mges = (body(i1)+body(i2))*zmbar
      fac1 = rbar*PARSEC/ASTUNIT

      if (rapo.ne.rper) then
        x1 = (rper/rapo)**2-1.d0
        vper = sqrt(2.d0*(GAU*mges/rapo/fac1-GAU*mges/rper/fac1)/x1)
        vapo = rper/rapo*vper
      else
        vper = sqrt(GAU*mges/rapo/fac1)
        vapo = vper
      end if

      amax = body(i1)*body(i2)/(-2.d0*ebin)
      ecc  = (rapo-rper)/(rapo+rper)

      if (mode.eq.1) comms = ' '
      if (mode.eq.2) comms = 'E'

      if (i1.gt.n) comms = comms(1:1)//'T1'
      if (i2.gt.n) comms = comms(1:1)//'T2'

      if (nmerge.ne.0) then
        do ii=1,nmerge
          if (-namem(ii).eq.name(i1)) comms = comms(1:1)//'T1'
          if (-namem(ii).eq.name(i2)) comms = comms(1:1)//'T2'
        end do
      end if

      write (38,'(f10.3,2i7,2i3,a3,2f8.4,2f11.5,2e15.6,1x,a3)')
     *  ttot,name(i1),name(i2),kstar(i1),kstar(i2),' N ',body(i1)
     *    *zmbar,body(i2)*zmbar,rdc,ecc,amax,ebini,comms

      v1 = sqrt((xdot(1,i1)-vcmx)**2+(xdot(2,i1)-vcmy)**2+
     *    (xdot(3,i1)-vcmz)**2)

      v2 = sqrt((xdot(1,i2)-vcmx)**2+(xdot(2,i2)-vcmy)**2+
     *    (xdot(3,i2)-vcmz)**2)

      write (38,'(f10.3,2i7,2i3,a3,2f8.4,2e13.5,2f13.7,1x,a3)')
     *  ttot*tstar,name(i1),name(i2),kstar(i1),kstar(i2),' P ',body(i1)
     *    *zmbar,body(i2)*zmbar,rper*fac1,rapo*fac1,vper,
     *     vapo,comms

      if (kz(14).eq.0.or.rdc.lt.rtide) nbina = nbina + 1

      end if

      return
      end
