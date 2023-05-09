      subroutine auswert (xa, va, bodya, phi)
*
* Berechnet Lag. Radii, Bound Mass + Cluster Centre Position...
*
      include 'common6.h'

      integer i,j,k,bbn,nstar,nst(19),nbina
      integer bin,nbin(50),nn,ntyg(15)
      real*8 rad(nmax),mass(nmax),mbin(50),bodya(nmax)
      real*8 mtyg(15),mex,mg,find_rt
      real*8 xa(3,nmax),va(3,nmax),phi(nmax)
      real*8 ro,rl(0:19),rh2,dis2,mlow,mup,
     *   inm,inma,mbound,meanm(19),rc2
c      real*8 MBULGE,MDISK,MHALO
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,OMEGA,
     &       MBULGE,RBULGE,MDISK,SCALEA,SCALEB,MHALO,RHALO,VCIRC

      mlow = min(body1,bodyn)
      mup  = max(body1,bodyn)

      if (kz(14).gt.1) then
        call gcint
      end if

      nstar = ntot-(ifirst-1)/2

      inm = 0.d0
      do k=1,n
        inm = inm + bodya(k)
      end do

 501  rtide=find_rt(kz(14),inm,rdens)

      rt2 = rtide**2
      inma=inm
      inm = 0.d0

      do k=1,n
        dis2=((rdens(1)-xa(1,k))**2+(rdens(2)-xa(2,k))**2+
     *     (rdens(3)-xa(3,k))**2)
        if (dis2.lt.rt2) then
          inm=inm+bodya(k)
        end if
      end do

      if (abs(inm-inma).gt.1.E-10) goto 501

      write(*,'(/a,2f13.6/)')'Tidal radius: NBODY/phys ',rtide,
     * rtide*rbar

c     dvr(1) = 0.d0
c     dvr(2) = 0.d0
c     dvr(3) = 0.d0

      nin = 0

      do j=ifirst,ntot
        dis2=((rdens(1)-xa(1,j))**2+(rdens(2)-xa(2,j))**2+
     *     (rdens(3)-xa(3,j))**2)
        if (dis2.lt.rt2) then
c         do k=1,3
c           dvr(k) = dvr(k) + va(k,j)
c         end do
          nin = nin + 1
        end if
      end do

c     if (nin.gt.0) then
c       do k=1,3
c         dvr(k) = dvr(k) / (1.d0*nin)
c       end do
c     end if

c get binary statistics
      if (nbin0.gt.0) then
        call binlist(nbina)
      else
        nbina = 0
      end if

      nimbh = 0

      if (kz(50).eq.1) then
        do ii=1,n
          if (name(ii).eq.1) then
             nimbh=ii
          end if
        end do
        if (nimbh.eq.0) then
          write (*,*) "BH Particle not identified in auswert !"
          stop
        end if
      end if

      nbound = 0
      mbound = 0.d0
      nn = 0
      bbn = 0

      do k=1,n
        dis2=((rdens(1)-xa(1,k))**2+(rdens(2)-xa(2,k))**2+
     *     (rdens(3)-xa(3,k))**2)
        if (dis2.lt.rt2) then
            nbound=nbound+1
            mbound=mbound+bodya(k)

            if (k.ne.nimbh.and.bodya(k)*zmbar.lt.150.0.
     *               and.bodya(k).gt.0.0) then
              nn=nn+1
              rad(nn)  = sqrt(dis2)
              mass(nn) = bodya(k)
          end if
        end if
      end do

      call sort2(nn, rad, mass)

C Kumulieren der Masse
      do k=2,nn
        mass(k)=mass(k)+mass(k-1)
      end do

      rl(0) = 0.d0
      rl(19) = 0.d0

      do k=1,19
        if (k.lt.11) then
          ro = 1.d0*k*mass(nn)/100.d0
        else
          ro = 1.d0*(k-9)*mass(nn)/10.d0
        end if

        if (mass(1).gt.ro) then
          rl(k)=rad(1)*ro/mass(1)
        else
          do j=2,nn
            if (mass(j-1).lt.ro.and.mass(j).ge.ro) then
               rl(k)=rad(j-1)+(rad(j)-rad(j-1))*(ro-mass(j-1))/
     *               (mass(j)-mass(j-1))
            end if
          end do
        end if
      end do
      if (rl(19).eq.0.d0) then
       if (nn.gt.0) then
        rl(19)=rad(nn)
       else
        rl(19) = 0.0
       end if
      end if

      if (body1.ne.bodyn) then

      do k=1,19
        meanm(k) = 0.d0
        nst(k) = 0
      end do

      n = ntot-(ifirst-1)/2

      do i=1,n
        dis2=sqrt((rdens(1)-xa(1,i))**2+(rdens(2)-xa(2,i))**2+
     *     (rdens(3)-xa(3,i))**2)
        if (i.ne.nimbh.and.bodya(i)*zmbar.lt.150.0.and.
     *              bodya(i)*zmbar.gt.0.0) then
        do k=1,19
          if (dis2.le.rl(k).and.dis2.gt.rl(k-1)) then
            if (i.le.n) then
              meanm(k) = meanm(k)+bodya(i)
              nst(k) = nst(k) + 1
            end if
          end if
        end do
        end if
      end do

      do k=1,19
        if (nst(k).gt.0) then
          meanm(k)=meanm(k)/(1.d0*nst(k))*zmbar
        else
          meanm(k) = 0.d0
        end if
      end do

      write(34,'(2f13.4,20(f13.6,2x))') ttot,ttot*tstar,
     *     (meanm(k),k=1,19),mass(nn)/nn*zmbar
      call flush(34)


      do k=1,19
       if (nst(k).eq.0) then
          write (*,*) 'Radial bin without star: ',ttot,k,rl(k),rl(k-1)
         if (k.gt.1) meanm(k) = meanm(k-1)
       end if
      end do

      rc2 = rl(3)**2
      rh2 = rl(14)**2
      rt2 = rtide**2

      do bin=1,10
        mbin(bin)=log10(mlow)+(log10(mup)-log10(mlow))*(bin-0.5)/10.d0
        nbin(bin) = 0
      end do

      do i=ifirst,ntot
        dis2=((rdens(1)-xa(1,i))**2+(rdens(2)-xa(2,i))**2+
     *     (rdens(3)-xa(3,i))**2)
       if (dis2.lt.rc2.and.bodya(i).gt.0.0) then
        if (i.le.n) then
         xm = bodya(i)*zmbar
        bin=1+int((log10(xm)-log10(mlow))/(log10(mup)-log10(mlow))*10.0)
         if (bin.le.10.and.bin.ge.1) nbin(bin) = nbin(bin)+1
        else
         i1 = 2*(i-n)-1
         i2 = 2*(i-n)
          xm = bodya(i1)*zmbar
          bin = 1 + int((log10(xm)-log10(mlow))/
     *             (log10(mup)-log10(mlow))*10.0)
          if (bin.le.10.and.bin.ge.1) nbin(bin) = nbin(bin)+1
          xm = bodya(i2)*zmbar
          bin = 1 + int((log10(xm)-log10(mlow))/
     *             (log10(mup)-log10(mlow))*10.0)
          if (bin.le.10.and.bin.ge.1) nbin(bin) = nbin(bin)+1
        end if
       end if
      end do

      write(35,'(2f13.4,i4,10(f12.6,i6))') ttot,ttot*tstar,1,
     *     (10**mbin(k),nbin(k),k=1,10)
      call flush(35)

      do bin=1,10
        nbin(bin) = 0
      end do

      do i=ifirst,ntot
        dis2=((rdens(1)-xa(1,i))**2+(rdens(2)-xa(2,i))**2+
     *     (rdens(3)-xa(3,i))**2)
       if (dis2.lt.rh2.and.bodya(i).gt.0.0) then
        if (i.le.n) then
         xm = bodya(i)*zmbar
        bin=1+int((log10(xm)-log10(mlow))/(log10(mup)-log10(mlow))*10.0)
         if (bin.le.10.and.bin.ge.1) nbin(bin) = nbin(bin)+1
        else
         i1 = 2*(i-n)-1
         i2 = 2*(i-n)
          xm = bodya(i1)*zmbar
          bin = 1 + int((log10(xm)-log10(mlow))/
     *             (log10(mup)-log10(mlow))*10.0)
          if (bin.le.10.and.bin.ge.1) nbin(bin) = nbin(bin)+1
          xm = bodya(i2)*zmbar
          bin = 1 + int((log10(xm)-log10(mlow))/
     *             (log10(mup)-log10(mlow))*10.0)
          if (bin.le.10.and.bin.ge.1) nbin(bin) = nbin(bin)+1
        end if
       end if
      end do

      write(35,'(2f13.4,i4,10(f12.6,i6))') ttot,ttot*tstar,2,
     *     (10**mbin(k),nbin(k),k=1,10)
      call flush(35)

      do bin=1,10
        nbin(bin) = 0
      end do

      do i=ifirst,ntot
        dis2=((rdens(1)-xa(1,i))**2+(rdens(2)-xa(2,i))**2+
     *     (rdens(3)-xa(3,i))**2)
       if (dis2.lt.rt2.and.bodya(i).gt.0.0) then
        if (i.le.n) then
         xm = bodya(i)*zmbar
        bin=1+int((log10(xm)-log10(mlow))/(log10(mup)-log10(mlow))*10.0)
        if (bin.le.10.and.bin.ge.1) nbin(bin) = nbin(bin)+1
        else
         i1 = 2*(i-n)-1
         i2 = 2*(i-n)
          xm = bodya(i1)*zmbar
          bin = 1 + int((log10(xm)-log10(mlow))/
     *             (log10(mup)-log10(mlow))*10.0)
          if (bin.le.10.and.bin.ge.1) nbin(bin) = nbin(bin)+1
          xm = bodya(i2)*zmbar
          bin = 1 + int((log10(xm)-log10(mlow))/
     *             (log10(mup)-log10(mlow))*10.0)
          if (bin.le.10.and.bin.ge.1) nbin(bin) = nbin(bin)+1
        end if
       end if
      end do

      write(35,'(2f13.4,i4,10(f12.6,i6))') ttot,ttot*tstar,3,
     *     (10**mbin(k),nbin(k),k=1,10)
      call flush(35)

      end if

      write(*,'(a,2f11.4,i9,a,i6,a,f8.5)')  'Auswertung: T = ',ttot,
     * ttot*tstar,nbound,'  Sterne gebunden ',nbina,'    Mtot= ',mbound

      write(*,'(/a,10(e11.4)/)') 'Lagr.-Radien (geb. Teilchen):',
     *   (rl(2*k),k=1,5),(rl(2*k-1),k=6,10)

      dmtot=ZMDOT/ZMBAR

      write(31,'(2f13.4,3i7,2f13.4)')
     *       ttot,ttot*tstar,nbound,nbina,bbn,dmtot,mbound
      call flush(31)

      write(32,'(2f13.4,19(e14.6,2x))') ttot,ttot*tstar,
     *     (rl(k),k=1,19)
      call flush(32)

      if (KZ(19).ne.0) then

      mg  = 0.d0
      mex = 0.d0

      do j=1,15
       ntyg(j) = 0
       mtyg(j) = 0.d0
      end do

      do i=1,n
        dis2=((rdens(1)-xa(1,i))**2+(rdens(2)-xa(2,i))**2+
     *     (rdens(3)-xa(3,i))**2)
       if (dis2.lt.rt2.and.bodya(i).gt.0.0) then
        do j=0,14
          if (kstar(i).eq.j) then
            ntyg(j+1) = ntyg(j+1) + 1
            mtyg(j+1) = mtyg(j+1) + bodya(i)
          end if
        end do

        mg = mg + bodya(i)
        if (kstar(i).gt.1) then
          mex = mex + bodya(i)
        end if
       end if
      end do

      write (*,*) 'Mex/Mg: ',mex/mg

      do j=1,15
       mtyg(j) = mtyg(j) / mg
      end do

      write (*,'(/a,15(i3,i7,f9.6))') 'Mver: ',
     *                  (j-1,ntyg(j),mtyg(j),j=1,15)

      write (36,'(2f13.6,15(i3,i7,f9.6))')
     *   ttot,ttot*tstar,(j-1,ntyg(j),mtyg(j),j=1,15)
      call flush(36)
      end if

      write(37,'(2f13.4,a3,6e14.06)') ttot,ttot*tstar,' 1 ',
     *    rdens(1),rdens(2),rdens(3),dvr(1),dvr(2),dvr(3)
      call flush(37)

      if (kz(14).gt.1) then
        write(37,'(2f13.4,a3,6e14.06)') ttot,ttot*tstar,' 2 ',
     *   rg(1),rg(2),rg(3),vg(1),vg(2),vg(3)
        call flush(37)
      end if

      return
      end

c -------------------------------------------------------------------------

      real*8 function find_rt(gmpar,mcl,rdens)
      IMPLICIT REAL*8  (A-H,O-Z)
      integer gmpar,i
      real*8 mcl,rgal,xg(3),ff,ffa,dcl,fs(3),fm(3),fmd(3),fsd(3)
      real*8 vrad,vper,omg,rdens(3)
      real*8 MBULGE,MDISK,MHALO
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,OMEGA,
     &       MBULGE,RBULGE,MDISK,SCALEA,SCALEB,MHALO,RHALO,VCIRC

      if (gmpar.gt.1) then
         rgal = sqrt((rg(1)-rdens(1))**2+(rg(2)-rdens(2))**2+
     *         (rg(3)-rdens(3))**2)
      end if

      if (gmpar.eq.1) then
        find_rt = (mcl/(3.d0*omega**2))**(1.d0/3.d0)
      else if (gmpar.eq.2) then
        find_rt = (mcl/(3.0*gmg))**(1.0/3.0)*rgal
      else if (gmpar.eq.3) then
        find_rt = (mcl/(2.d0*vcirc**2.0))**(1.0/3.0)*rgal**(2.0/3.0)
      else if (gmpar.eq.5 .or. gmpar.eq.6) then
         find_rt = 0.d0
         ffa = 0.d0
         if (gmpar.eq.5) then
            CALL FORCEIR13(RG,VG,FS,FSD)
         else
            CALL FORCEBO15(RG,VG,FS,FSD)
         endif
         xg(2) = rg(2)
         xg(3) = rg(3)
         vrad=(vg(1)*rg(1)+vg(2)*rg(2)+vg(3)*rg(3))
     *           /sqrt(rg(1)**2+rg(2)**2+rg(3)**2)
         vper = vg(1)**2+vg(2)**2+vg(3)**2-vrad**2
         if (vper.gt.0.d0) then
            omg=sqrt(vper)/sqrt(rg(1)**2+rg(2)**2+rg(3)**2)
            do i=1,10000
              dcl=0.05*i
              xg(1)=rg(1)+dcl
              if (gmpar.eq.5) then
                 CALL FORCEIR13(RG,VG,FM,FMD)
              else
                 CALL FORCEBO15(RG,VG,FM,FMD)
              endif
              ff=fm(1)-fs(1)+omg**2.0*dcl-mcl/dcl**2.0
              if (ffa.le.0.d0.and.ff.gt.0.d0) then
                find_rt = dcl
              end if
              ffa = ff
            end do
         end if
         if (find_rt.eq.0.d0) then
           find_rt=1.E10
           write (*,*) "Couldn't determine tidal radius !"
         end if
      end if

      return
      end
