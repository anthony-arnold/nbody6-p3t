      subroutine dcen (xa, va, bodya)
      include 'common6.h'
      integer i,j,k,l,nstar,merk
      real*8 dens(nmax),xa(3,nmax),va(3,nmax),bodya(nmax)
      real*4 rad(nmax)
      real*8 sd,mind
      INTEGER casopt
      real*8 eps2
      parameter (eps2 = 0.d0)

      casopt = max(min(7,ntot/2),3)

      nmx = min(10000+ifirst, ntot)
      rhom = 0.d0

!$omp parallel do private(J,K,L,RAD,MIND,MERK)
      do j=ifirst,nmx
        do k=ifirst,ntot
          rad(k-ifirst+1) = sqrt((xa(1,j)-xa(1,k))**2+
     *         (xa(2,j)-xa(2,k))**2+(xa(3,j)-xa(3,k))**2)
          if (rad(k-ifirst+1).eq.0.d0) rad(k-ifirst+1)=1.E10
        end do
        rad(j-ifirst+1) = 2.e30   ! Abstand von Teilchen zu sich selbst

        do l=1,casopt
          mind = 1.E30
          do i=1,ntot-ifirst+1
            if (rad(i).lt.mind) then
              mind = rad(i)
              merk = i
            end if
          end do
          rad(merk) = 1.E30
        end do
        dens(j) = bodya(j)/mind**3.d0  ! rad(casopt)
      end do
!$omp end parallel do

      nstar = ntot-(ifirst-1)/2

      do k=1,3
        rdens(k) = 0.d0
        dvr(k) = 0.d0
      end do

      sd = 0.d0
      do j=ifirst,nmx
        do k=1,3
          rdens(k) = rdens(k) + xa(k,j) * dens(j)
c         dvr(k) = dvr(k) + va(k,j) * dens(j)
        end do
        sd = sd + dens(j)
      end do

      do k=1,3
        rdens(k) = rdens(k) / sd
c       dvr(k) = dvr(k) / sd
      end do

      write (*,*) "in dcen: ",rdens(1),rdens(2),rdens(3)

c use all stars to calculate cluster velocity
      do j=ifirst,ntot
        do k=1,3
          dvr(k) = dvr(k) + va(k,j)
        end do
      end do

      do k=1,3
        dvr(k) = dvr(k) / (1.0*(ntot-ifirst))
      end do

      return
      end
