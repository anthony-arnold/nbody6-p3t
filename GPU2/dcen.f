      subroutine dcen (xa, va, bodya)
      USE OMP_LIB
      include 'common6.h'
      integer i,j,k,l,nstar,merk(0:100)
      real*8 dens(nmax),xa(3,nmax),va(3,nmax),bodya(nmax)
      real*8 rad(nmax)
      real*8 sd,mind
      INTEGER casopt
      real*8 eps2
      parameter (eps2 = 0.d0)
      integer*4 jj,istart,iend,merkk,ithreadm

      casopt = max(min(7,ntot/2),3)

      ithreadm = OMP_get_max_threads()

      nmx = min(10000+ifirst, ntot)
      rhom = 0.d0

      do j=ifirst,nmx
!$omp parallel do
        do k=ifirst,ntot
          rad(k-ifirst+1) = (xa(1,j)-xa(1,k))**2+
     *         (xa(2,j)-xa(2,k))**2+(xa(3,j)-xa(3,k))**2
          if (rad(k-ifirst+1).eq.0.d0) rad(k-ifirst+1)=1.E10
        end do
!$omp end parallel do
        rad(j-ifirst+1) = 2.e30   ! Abstand von Teilchen zu sich selbst

        do l=1,casopt
!$omp parallel private (jj,merkk,i,istart,iend)
          jj = OMP_GET_THREAD_NUM()

          istart = jj*(ntot-ifirst+2)/ithreadm
          iend   = (jj+1)*(ntot-ifirst+2)/ithreadm-1
          merkk = istart

          do i=istart+1,iend
            if (rad(i).lt.rad(merkk)) merkk=i
          end do

          merk(jj)=merkk
!$omp end parallel

          merkk = merk(0)
          do i=1,ithreadm-1
            if (rad(merk(i)).lt.rad(merkk)) merkk=merk(i)
          end do
          mind = rad(merkk)
          rad(merkk) = 1.E30
        end do
        dens(j) = bodya(j)/mind**1.5d0  ! rad(casopt)
      end do

      nstar = ntot-(ifirst-1)/2

      do k=1,3
        rdens(k) = 0.d0
        dvr(k) = 0.d0
      end do

      sd = 0.d0
!$omp parallel do private(k) reduction(+:sd,rdens)
      do j=ifirst,nmx
        do k=1,3
          rdens(k) = rdens(k) + xa(k,j) * dens(j)
        end do
        sd = sd + dens(j)
      end do
!$omp end parallel do

      do k=1,3
        rdens(k) = rdens(k) / sd
      end do

c calculate half-mass radius
!$omp parallel do
      do j=ifirst,ntot
       rad(j-ifirst+1)=(xa(1,j)-rdens(1))**2+(xa(2,j)-rdens(2))**2+
     *      (xa(3,j)-rdens(3))**2
      end do
!$omp end parallel do

      call sort(ntot-ifirst+1,rad)

      rh = sqrt(rad(int(0.5*(ntot-ifirst+1))))
      rh2 = rh*rh

      sd = 0.0
!$omp parallel do reduction(+:sd,dvr)
      do j=ifirst,ntot
        rad(j)=(xa(1,j)-rdens(1))**2+(xa(2,j)-rdens(2))**2
     *      +(xa(3,j)-rdens(3))**2
        if (rad(j).lt.rh2) then
          dvr(1) = dvr(1) + va(1,j) * bodya(j)
          dvr(2) = dvr(2) + va(2,j) * bodya(j)
          dvr(3) = dvr(3) + va(3,j) * bodya(j)
          sd = sd + bodya(j)
        end if
      end do
!$omp end parallel do

      do k=1,3
        dvr(k) = dvr(k) / sd
      end do

      if (kz(50).eq.1) then
        nbh = 0
        do ii=1,ntot
          if (name(ii).eq.1) then
            if (bodya(ii).eq.0.d0) then
              ifd = 0
              do iii=1,ifirst-1
                if(name(n+kvec(iii)).lt.0.and.(iii.EQ.2*KVEC(iii)-1.or.
     *            NAME(2*kvec(iii)).GT.NZERO)) THEN
                  call FINDJ(III,IGHOST,IM)
                  if (ighost.eq.ii) then
                      ifd=1
                       write (*,*) "Identified BH ghost: ",n+kvec(iii)
                      nbh=n+kvec(iii)
                  end if
                end if
              end do
              if (ifd.eq.0) then
                write (*,*) "Unidentified BH ghost !",ifirst,ighost
                goto 100
              end if
            else if (ii.ge.ifirst) then
               nbh=ii
            else
               nbh=n+kvec(ii)
            end if
          end if
        end do
        if (nbh.eq.0) then
          write (*,*) "BH Particle not identified in dcen !"
          stop
        end if

        dis = sqrt((rdens(1)-xa(1,nbh))**2+(rdens(2)-xa(2,nbh))**2+
     *       (rdens(3)-xa(3,nbh))**2)

        write (*,*) "Distance of BH particle to centre: ",dis

        if (dis.gt.1E2) then
           write (*,*) "BH Particle outside core !"
           write (*,*) "Core position: ",rdens(1),rdens(2),rdens(3)
           write (*,*) "BH position: ",xa(1,nbh),xa(2,nbh),xa(3,nbh)
           stop
        end if

c use IMBH position as density centre
        rdens(1)=xa(1,nbh)
        rdens(2)=xa(2,nbh)
        rdens(3)=xa(3,nbh)
      end if

 100  return
      end
