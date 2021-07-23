      subroutine forceir13 (rp, vp, fp, fdp) 
c 
c calculates force and 1st derivative for a particle moving
c in the galactic model of Irrgang et al. (2013), Table 1
c units: physical units 
c
c     implicit none
      include "common6.h"
      real*8 rp(3), vp(3), ru(3), vu(3), fp(3), fdp(3)
      real*8 m1,b1,m2,a2,b2,m3,a3, rv, xd3, xd5
      parameter (m1=9.5E9 , b1=230.0, m2 = 6.6E10, a2 = 4220.0, 
     *   b2 = 292.0, m3=2.335E10, a3=2562.0)
      real*8 f1(3),f2(3),f3(3),fd1(3),fd2(3),fd3(3)
      real*8 r,xc1,xc2,xc3

c     a+b parameters in pc, m parameters in Msun

      ru(1) = rp(1)*rbar
      ru(2) = rp(2)*rbar
      ru(3) = rp(3)*rbar

      vu(1) = vp(1)*vstar
      vu(2) = vp(2)*vstar
      vu(3) = vp(3)*vstar

      dis  = sqrt(ru(1)**2+ru(2)**2+ru(3)**2)
      rv = ru(1)*vu(1) + ru(2)*vu(2) + ru(3)*vu(3)

c bulge
      xc1 = sqrt(b1**2+dis**2)
      xd3 = xc1**3
      xd5 = xc1**5

      f1(1) = -m1*ru(1)/xd3
      f1(2) = -m1*ru(2)/xd3
      f1(3) = -m1*ru(3)/xd3

      fd1(1) = -m1*(vu(1)/xd3 - 3.d0 * rv * ru(1) / xd5)
      fd1(2) = -m1*(vu(2)/xd3 - 3.d0 * rv * ru(2) / xd5)
      fd1(3) = -m1*(vu(3)/xd3 - 3.d0 * rv * ru(3) / xd5)

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
      f3(1) =-m3*ru(1)/(dis*a3*a3+dis*dis*a3)
      f3(2) =-m3*ru(2)/(dis*a3*a3+dis*dis*a3)
      f3(3) =-m3*ru(3)/(dis*a3*a3+dis*dis*a3)

      fd3(1) = -m3*vu(1)/(dis*a3*a3+dis*dis*a3)
     *   +m3*ru(1)*rv/(dis*a3*a3+dis*dis*a3)**2*(a3*a3/dis+2.0*a3)
      fd3(2) = -m3*vu(2)/(dis*a3*a3+dis*dis*a3)
     *   +m3*ru(2)*rv/(dis*a3*a3+dis*dis*a3)**2*(a3*a3/dis+2.0*a3)
      fd3(3) = -m3*vu(3)/(dis*a3*a3+dis*dis*a3)
     *   +m3*ru(3)*rv/(dis*a3*a3+dis*dis*a3)**2*(a3*a3/dis+2.0*a3)

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


      subroutine orb_energchk_irr13(rp, vp)
c
c calculates energy and angular momentum of a particle moving 
c in the galactic model of Irrgang et al. (2013), Table 1
c units: pc, pc/myr
c
c
      include 'common6.h'
      real*8 rp(3), vp(3), phi(3)
      real*8 m1,b1,m2,a2,b2,m3,a3
      parameter (m1=9.5E9 , b1=230.0, m2 = 6.6E10, a2 = 4220.0,
     *   b2 = 292.0, m3=2.335E10, a3=2562.0)
      real*8 lx,ly,lz,ltot,lde
      real*8 ru(3),vu(3)

      ru(1) = rp(1)*rbar
      ru(2) = rp(2)*rbar
      ru(3) = rp(3)*rbar

      vu(1) = vp(1)*vstar
      vu(2) = vp(2)*vstar
      vu(3) = vp(3)*vstar

      vg = sqrt(vu(1)**2+vu(2)**2+vu(3)**2)
      rgal = sqrt(ru(1)**2+ru(2)**2+ru(3)**2)

      phi(1) = -m1/sqrt(ru(1)**2+ru(2)**2+ru(3)**2+b1**2)

      phi(2) = -m2/sqrt(ru(1)**2+ru(2)**2+(a2+
     *      sqrt(ru(3)**2+b2**2))**2)

      phi(3) = m3/a3*log(1.0+rgal/a3) 
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


     subroutine checkend (rp, vp)
      include "common6.h"
      real*8 rp(3), vp(3), ru(3), vu(3), dis, vr

      ru(1) = rp(1)*rbar
      ru(2) = rp(2)*rbar
      ru(3) = rp(3)*rbar

      vu(1) = vp(1)*vstar
      vu(2) = vp(2)*vstar
      vu(3) = vp(3)*vstar

      dis=sqrt((rfin(1)-ru(1))**2+(rfin(2)-ru(2))**2+(rfin(3)-ru(3))**2)
      vr = (ru(1)-rfin(1))*vu(1)+(ru(2)-rfin(2))*vu(2)+
     *       (ru(3)-rfin(3))*vu(3)

      if (vrfinold.le.0.d0.and.vr.gt.0.d0) then
             vrfinold=vr
             call output
             stop
      end if

      vrfinold = vr

      return
      end
