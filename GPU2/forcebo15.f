      subroutine forcebo15 (rp, vp, fp, fdp)
c
c calculates force and 1st derivative for a particle moving
c in the galactic model of Bovy (2015), galpy
c units: physical units
c
c     implicit none
      include "common6.h"
      real*8 rp(3), vp(3), ru(3), vu(3), fp(3), fdp(3)
      real*8 mb,rhobm,bulgecut,m2,a2,b2,rhodm,rsdm, rv, xd3, xd5
      parameter(mb=5.0E9, rhobm=-1.8d0, bulgecut=1900.d0, m2=6.8E10,
     *  a2=3000.d0, b2=280.d0, rhodm=7.1111E-3, rsdm=16000.d0)
      real*8 f1(3),f2(3),f3(3),fd1(3),fd2(3),fd3(3)
      real*8 r,mbulge,mhalo,mhalo2,dis3,dis5,xf

c     a+b parameters in pc, m parameters in Msun
      ru(1) = rp(1)*rbar
      ru(2) = rp(2)*rbar
      ru(3) = rp(3)*rbar

      vu(1) = vp(1)*vstar
      vu(2) = vp(2)*vstar
      vu(3) = vp(3)*vstar

      dis  = sqrt(ru(1)**2+ru(2)**2+ru(3)**2)
      rv = ru(1)*vu(1) + ru(2)*vu(2) + ru(3)*vu(3)

      dis3 = dis*dis*dis;
      dis5 = dis3*dis*dis;

      do k=1,3
        f3(k)=0.d0
        fd3(k)=0.d0
        f2(k)=0.d0
        fd2(k)=0.d0
      end do

c bulge
      mbulge = mb*gammp(1.5+rhobm/2.0,dis*dis/bulgecut/bulgecut)
      xf = rv/dis3/dis*(dis/bulgecut)**(3.0-rhobm)*
     *   exp(-dis*dis/bulgecut/bulgecut)

      f1(1) = -mbulge*ru(1)/dis3
      f1(2) = -mbulge*ru(2)/dis3
      f1(3) = -mbulge*ru(3)/dis3

      fd1(1) = -mbulge*vu(1)/dis3+(3.d0*mbulge*rv/dis5-mbulge*xf)*ru(1)
      fd1(2) = -mbulge*vu(2)/dis3+(3.d0*mbulge*rv/dis5-mbulge*xf)*ru(2)
      fd1(3) = -mbulge*vu(3)/dis3+(3.d0*mbulge*rv/dis5-mbulge*xf)*ru(3)

c disc
      xc1=sqrt(ru(1)**2+ru(2)**2+(a2+sqrt(ru(3)**2+b2**2))**2)
      xc2 = a2+sqrt(ru(3)**2+b2**2)
      xc3 = sqrt(ru(3)**2+b2**2)

      f2(1) = -m2*ru(1)/xc1**3
      f2(2) = -m2*ru(2)/xc1**3
      f2(3) = -m2*ru(3)/xc1**3/
     *   sqrt(ru(3)**2+b2**2)*(a2+sqrt(ru(3)**2+b2**2))

      fd2(1) = -m2*vu(1)/xc1**3+3*m2*ru(1)/xc1**5*(rv+
     *    ru(3)*vu(3)*a2/xc3)
      fd2(2) = -m2*vu(2)/xc1**3+3*m2*ru(2)/xc1**5*(rv+
     *    ru(3)*vu(3)*a2/xc3)
      fd2(3) = -m2*vu(3)/xc1**3+
     *   m2*xc2/xc3/xc1**5*(3*ru(3)*rv+ru(3)**2*vu(3)*a2/xc3)+
     *    m2*vu(3)/xc1**3*(ru(3)**2*a2/xc3**3-xc2/xc3)

c halo
      mhalo  = 4.0*PI*rhodm*rsdm**3.0*(log(1.0+dis/rsdm)-dis/(dis+rsdm))
      mhalo2 = 4.0*PI*rhodm*rsdm**3.0/(dis+rsdm)**2

      f3(1) =-mhalo*ru(1)/dis3
      f3(2) =-mhalo*ru(2)/dis3
      f3(3) =-mhalo*ru(3)/dis3

      fd3(1) = -mhalo*vu(1)/dis3+3.0*mhalo*rv*ru(1)/dis5
     *   -mhalo2*rv*ru(1)/dis3
      fd3(2) = -mhalo*vu(2)/dis3+3.0*mhalo*rv*ru(2)/dis5
     *   -mhalo2*rv*ru(2)/dis3
      fd3(3) = -mhalo*vu(3)/dis3+3.0*mhalo*rv*ru(3)/dis5
     *   -mhalo2*rv*ru(3)/dis3

      fp(1) = (f1(1) + f2(1) + f3(1))  ! Transform force to Msun/pc^2
      fp(2) = (f1(2) + f2(2) + f3(2))
      fp(3) = (f1(3) + f2(3) + f3(3))

      fdp(1) = (fd1(1)+fd2(1)+fd3(1))  ! Fdot in Msun/pc^2/Myr
      fdp(2) = (fd1(2)+fd2(2)+fd3(2))
      fdp(3) = (fd1(3)+fd2(3)+fd3(3))

      fp(1) = fp(1)/zmbar*rbar**2
      fp(2) = fp(2)/zmbar*rbar**2
      fp(3) = fp(3)/zmbar*rbar**2

      fdp(1) = fdp(1)/zmbar*rbar**2/vstar
      fdp(2) = fdp(2)/zmbar*rbar**2/vstar
      fdp(3) = fdp(3)/zmbar*rbar**2/vstar

      return
      end


      subroutine orb_energchk_bo15(rp, vp)
c
c calculates energy and angular momentum of a particle moving
c in the galactic model of Bovy (2015), galpot
c units: pc, pc/myr
c
c
      include 'common6.h'
      real*8 rp(3), vp(3), phi(3)
      real*8 mb,rhobm,bulgecut,m2,a2,b2,rhodm,rsdm, rv, xd3, xd5
      parameter(mb=5.0E9, rhobm=-1.8d0, bulgecut=1900.d0, m2=6.8E10,
     *  a2=3000.d0, b2=280.d0, rhodm=7.1111E-3, rsdm=16000.d0)
      real*8 lx,ly,lz,ltot,lde, drat2, dgamm
      real*8 ru(3),vu(3)
      real*8 gamrat
      parameter (gamrat=6.388367725022)  ! for rhobm = -1.8

      ru(1) = rp(1)*rbar
      ru(2) = rp(2)*rbar
      ru(3) = rp(3)*rbar

      vu(1) = vp(1)*vstar
      vu(2) = vp(2)*vstar
      vu(3) = vp(3)*vstar

      vg = sqrt(vu(1)**2+vu(2)**2+vu(3)**2)
      rgal = sqrt(ru(1)**2+ru(2)**2+ru(3)**2)

      drat2 = rgal*rgal/bulgecut/bulgecut;

c bulge
      phi(1) = mb*gamrat*gammp(1.0+rhobm/2.0,drat2)/bulgecut-
     *    mb*gammp(1.5+rhobm/2.0,drat2)/rgal;

      phi(2) = -m2/sqrt(ru(1)**2+ru(2)**2+(a2+
     *      sqrt(ru(3)**2+b2**2))**2)

c halo
      phi(3) = -4.0*PI*RHODM*RSDM*RSDM*RSDM*log(1.0+rgal/RSDM)/rgal;
      etot = G*(phi(1)+phi(2)+phi(3))+0.5d0*vg**2

      lx = ru(2)*vu(3)-ru(3)*vu(2)
      ly = ru(3)*vu(1)-ru(1)*vu(3)
      lz = ru(1)*vu(2)-ru(2)*vu(1)
      ltot = sqrt(lx**2+ly**2+lz**2)

      if (time+toff.eq.0.d0) then
         e0_orb=etot
         v0_orb=vg
         l0_orb=lz
      end if

      if (l0_orb.ne.0.d0) then
         lde = (lz-l0_orb)/l0_orb
      else
         lde = lz
      end if

      write (*,'(/a,f13.5,a,f12.5,2(a,f14.9))')
     *   'Cluster orbit: Dis= ',rgal,' Vel= ',vg,' Delta E = ',
     *    (etot-e0_orb)/v0_orb**2, ' Delta L = ',lde

      return
      end

      FUNCTION gammp(a,x)
      REAL*8 a,gammp,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)write(*,*)'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END


      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      write(*,*) 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END


      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)write(*,*) 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      write (*,*) 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
