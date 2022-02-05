      subroutine auswert_iso (xa, va, bodya, phi)
*
* Berechnet Lag. Radii, Bound Mass + Cluster Centre Position...
*
      include 'common6.h'

      integer i,j,k,bbn,nst(19),bin,nbin(50),nn,ntyg(15),nimbh,nbina
      real*8 rad(nmax),mass(nmax),mbin(50),mtyg(15),mex,mg
      real*8 bodya(nmax)
      real*8 xa(3,nmax),va(3,nmax),phi(nmax),ro,rl(0:19),rh2,dis2,
     *   mbound,meanm(19),rc2,estar,mtin,mlow,mup
*

      mlow = min(body1,bodyn)
      mup  = max(body1,bodyn)

      mtin = 0.d0

      do j=ifirst,ntot
        estar = phi(j)+0.5*((va(1,j)-dvr(1))**2+(va(2,j)-dvr(2))**2
     *      +(va(3,j)-dvr(3))**2)

c Assume all stars bound if cluster has reflective boundary
        if (estar.lt.0.d0.or.kz(29).ne.0) then
          mtin = mtin + bodya(j)
        end if
      end do

c get binary statistics
      if (nbin0.gt.0) then
        call binlist(nbina)
      else
        nbina = 0
      end if

      nbound = 0
      mbound = 0.d0
      nn = 0
      bbn = 0

C Check which star is IMBH
      nimbh = 0

      if (kz(50).eq.1) then
        do k=ifirst,ntot
          if (k.le.n) then
            if (name(k).eq.1) nimbh=k
          else
            i1 = 2*(k-n)-1
            i2 = 2*(k-n)
            if (name(i1).eq.1.or.name(i2).eq.1) nimbh=k
          end if
        end do
      end if

      do k=ifirst,ntot
        dis2=((rdens(1)-xa(1,k))**2+(rdens(2)-xa(2,k))**2+
     *     (rdens(3)-xa(3,k))**2)
        estar = phi(k)+0.5*((va(1,k)-dvr(1))**2+(va(2,k)-dvr(2))**2
     *      +(va(3,k)-dvr(3))**2)

c Assume all stars bound if cluster has reflective boundary
        if (estar.lt.0.d0.or.kz(29).ne.0) then
          if (k.le.n) then
            nbound=nbound+1
            mbound=mbound+bodya(k)
          else
            nbound=nbound+2
            i1 = 2*(k-n)-1
            i2 = 2*(k-n)
            mbound=mbound+bodya(i1)
            mbound=mbound+bodya(i2)
            bbn = bbn+1
          end if
          if (k.ne.nimbh) then
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
      if (rl(19).eq.0.d0) rl(19)=rad(nn)

      if (body1.ne.bodyn) then

      do k=1,19
        meanm(k) = 0.d0
        nst(k) = 0
      end do

      n = ntot-(ifirst-1)/2

      do i=ifirst,ntot
        dis2=sqrt((rdens(1)-xa(1,i))**2+(rdens(2)-xa(2,i))**2+
     *     (rdens(3)-xa(3,i))**2)
        if (i.ne.nimbh) then
          do k=1,19
            if (dis2.le.rl(k).and.dis2.gt.rl(k-1)) then
              if (i.le.n) then
                meanm(k) = meanm(k)+bodya(i)
                nst(k) = nst(k) + 1
              else
                i1 = 2*(i-n)-1
                i2 = 2*(i-n)
                meanm(k) = meanm(k)+bodya(i1)
                meanm(k) = meanm(k)+bodya(i2)
                nst(k) = nst(k) + 2
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

      do bin=1,10
        mbin(bin) = log10(mlow)+(log10(mup)-log10(mlow))*(bin-0.5)/10.d0
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
        estar = phi(i)+0.5*((va(1,i)-dvr(1))**2+(va(2,i)-dvr(2))**2
     *     +(va(3,i)-dvr(3))**2)
       if (estar.lt.0.d0.and.bodya(i).gt.0.0) then
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

      write(*,'(/a,2f11.4,i9,a,i6,a,f8.5)')  'Time: T = ',ttot,
     * ttot*tstar,nbound,'  stars bound ',nbina,'    Mtot= ',mbound

      write(*,'(/a,10(e11.4)/)') 'Lagr.-Radii (bound stars):',
     *   (rl(2*k),k=1,5),(rl(2*k-1),k=6,10)

      dmtot=ZMDOT/ZMBAR

c     dmtot=ZMDOT/(ZMASS*ZMBAR)

      write(31,'(2f13.4,3i7,2f13.4)')
     *       ttot,ttot*tstar,nbound,nbina,bbn,dmtot,mbound
      call flush(31)

      write(32,'(2f13.4,19(e14.6,2x))') ttot,ttot*tstar,
     *     (rl(k),k=1,19)
      call flush(32)

      if (KZ(20).gt.0) then

      mg  = 0.d0
      mex = 0.d0

      do j=1,15
       ntyg(j) = 0
       mtyg(j) = 0.d0
      end do

      do i=1,n
       estar = phi(i)+0.5*((va(1,i)-dvr(1))**2+(va(2,i)-dvr(2))**2
     *      +(va(3,i)-dvr(3))**2)
       if (estar.lt.0.d0.and.bodya(i).gt.0.0) then
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

      return
      end
