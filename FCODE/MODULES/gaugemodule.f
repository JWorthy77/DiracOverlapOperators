!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module dynamicmodule
      use pacc
      use arraysizes
      implicit none
      
      complex(prc) :: phi(Nv,4)
      real(prc) :: phir(Nv)

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setDynamicNoise()
      use rvmodule
      use options
      use gaugefield
      use basicdiracopsmod
      implicit none
      integer id
      real(prc) rv(Nv,2)
      complex(prc) :: eta(Nv,4)

      do id=1,4
        call gauss0(rv)
        eta(:,id)=cmplx(rv(:,1),rv(:,2),prc)
      end do
      call DWilson(eta,phi,u,.true.,baremass,czero)

      return
      end subroutine setDynamicNoise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module dynamicmodule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module gaugemodule
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makeGaugeField(GZERO)
      use gaugefield
      implicit none

      logical GZERO
      integer isw,Nsw
      real(prc) dH
      real(prc) thetat(Nv,3)

      if (GZERO) then
        theta=0
        call coef(u,theta)
        return
      end if

      Nsw=50
!     loop over Nsweep Hybrid MC steps
      do isw=1,Nsw
        print *,"sweep:",isw," of ",Nsw
        thetat=theta
        call march(dH,thetat)
        call accept(dH,thetat)
      end do
      call coef(u,theta)

      return
      end subroutine makeGaugeField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makeQuenchedGaussianThirringField()
      use gaugefield
      implicit none

      print *,"Make rvs for quenched Gaussian Thirring"
!     grv ~ 1/s.rt(2.pi).exp(-1/2.(x/s)^2)
!     s=1/beta^(1/2) =>
!     theta_i ~ (beta/(2.pi))^(1/2)*exp(-beta/2*theta^2)"
!     z=s.x+mu ~ N(mu,s)
      call setGRVs(Nv*3,theta)
      theta=theta/sqrt(gbeta)
      call coef(u,theta)

      return
      end subroutine makeQuenchedGaussianThirringField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makeQuenchedCosineThirringField()
      use gaugefield
      implicit none
      logical COSINE
      integer,parameter :: Ncentre=100
      integer,parameter :: Nedge=Ncentre+1
      real(prc),dimension(Ncentre) :: ctheta,dist
      real(prc),dimension(Nedge) :: etheta,cdist
      real(prc) dtheta,rv,mult,crv
      integer i,idx,mu

      print *,"Make rvs for quenched  cosine Thirring"
      dtheta=2*pi/Ncentre
      etheta(1)=-pi
      do i=2,Nedge
        etheta(i)=etheta(i-1)+dtheta
        ctheta(i-1)=(etheta(i-1)+etheta(i))/two
        dist(i-1)=cospdf(gbeta,ctheta(i-1))
      end do
      dist=dist/sum(dist)
      cdist(1)=zero
      do i=2,Nedge
        cdist(i)=cdist(i-1)+dist(i-1)
      end do

      do i=1,Nv
        do mu=1,3
          rv=urv()
          idx=distindex(rv,Nedge,cdist)
          mult=(rv-cdist(idx))/(cdist(idx+1)-cdist(idx))
          crv = etheta(idx)+mult*dtheta
          theta(i,mu)=crv
        end do
      end do
          
      call coef(u,theta)

      return
      end subroutine makeQuenchedCosineThirringField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      real(prc) function cospdf(beta,theta)
!      implicit none
!      real(prc) beta,theta
!
!      cospdf = exp(beta*(cos(theta)-one))
!      return
!      end function cospdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      integer function distindex(rv,Nd,dist)
c      implicit none
c      real(prc) rv
c      integer Nd
c      real(prc),dimension(Nd) :: dist
c      integer i
c
c      do i=1,Nd-1
c        if ((rv.ge.dist(i)).and.rv.le.dist(i+1)) then
c          distindex=i
c          return
c        endif
c      enddo
c      print *,"index not found"
c      stop
c      return
c      end function distindex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine force(theta,dSdpi)
c      use options
      implicit none
      real(prc) theta(Nv,3)
      real(prc) dSdpi(Nv,3)

      if (GAUGETYPE.eq.1) then 
        call gaugeforce(theta,dSdpi)
      elseif (GAUGETYPE.eq.2) then
        call gaugeforceThirring3(theta,dSdpi)
      elseif (GAUGETYPE.eq.3) then
        print *,"should be in forceSigma"
        stop
c        call gaugeforceGN(theta,dSdpi)
      endif

      if (DYNAMIC.gt.0) then
        print *,'DYNAMIC - force'
        call fermionforce(dSdpi) 
      endif

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gaugeforce(thetat,dSdpi)
c      use arraysizes
c      use indices
c      use options
      implicit none
      real(prc) thetat(Nv,3),dSdpi(Nv,3)
      integer mu,nu,i,ish(Nv)
      real(prc) Sigma(Nv),a(Nv)

      do mu=1,3

        Sigma=zero
        dSdpi(:,mu)=zero
        do nu=1,3
          if (nu.ne.mu) then
c           first the +nu staple...
            do i=1,Nv
              a(i)=thetat(iu(i,mu),nu)-thetat(iu(i,nu),mu)-thetat(i,nu)
            end do
            Sigma=Sigma+a

c  ... and then the -nu staple
            call gather(Nv,ish,id(1,nu),iu(1,mu))
            do i=1,Nv
         a(i)=-thetat(ish(i),nu)-thetat(id(i,nu),mu)+thetat(id(i,nu),nu)
            end do
            Sigma=Sigma+a
          end if
        end do
c (plaquette=staple +theta_link) and since we have added a forward
c and a backward staple in each of the 2 planes (rem 3 dimensions)
c we must add 4 theta_links to each staplesum.
        dSdpi(:,mu)=gbeta*(Sigma(:)+4.0*thetat(:,mu))
      end do

      return
      end subroutine gaugeforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gather(n,a,b,c)
      integer n,a(n),b(n),c(n)
      integer i
c
      do i=1,n
        a(i)=b(c(i))
      enddo
c
      return
      end subroutine gather
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gaugeforceThirring3(theta,dSdpi)
      implicit none
      real(prc) theta(Nv,3),dSdpi(Nv,3)

      dSdpi=Nferms*gbeta*theta

      return
      end subroutine gaugeforceThirring3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gaugeforceGN(sigma,dSdpi)
      implicit none
      real(prc) sigma(Nv),dSdpi(Nv)

      dSdpi=Nferms*gbeta*sigma

      return
      end subroutine gaugeforceGN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforce(dSdpi)
      use gammas
      use axbmodule1
      use gaugefield
      use dynamicmodule
      use basicdiracopsmod
      implicit none
      real(prc) dSdpi(Nv,3)
      complex(prc) eta(Nv,4)
      complex tzi
      integer mu,idirac,i
      integer igork1
      procedure(),pointer :: Mptr

      tzi=cmplx(zero,two)
      Mptr=>DWilson
      call IMdagM(phi,eta,u,.false.,baremass,czero,Mptr)

      do mu=1,3
        do idirac=1,4
          igork1=gamin(mu,idirac)
          do i=1,Nv
!           Wilson terms
            dSdpi(i,mu)=dSdpi(i,mu)+real(tzi*
     &  (conjg(eta(i,idirac))*u(i,mu)*eta(iu(i,mu),idirac)
     &  -conjg(eta(iu(i,mu),idirac))*conjg(u(i,mu))*eta(i,idirac)))
!           Dirac terms
            dSdpi(i,mu)=dSdpi(i,mu)+real(tzi*gamval(mu,idirac)*
     &  (conjg(eta(i,idirac))*u(i,mu)*eta(iu(i,mu),igork1)
     &  +conjg(eta(iu(i,mu),idirac))*conjg(u(i,mu))*eta(i,igork1)))
          end do
        end do
      end do

      return
      end subroutine fermionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine hamilton(theta,pp,h,hg,hp,hf)
c      use options
c      use arraysizes
      implicit none
      real(prc) theta(Nv,3),pp(Nv,3)
      real(prc) h,hg,hp,hf

      hp=0.5*sum(pp*pp)
      if (GAUGETYPE.eq.1) then
        call SGqed3(theta,hg)
      elseif (GAUGETYPE.eq.2) then
        call SGthirring3(theta,hg)
      elseif (GAUGETYPE.eq.3) then
        print *,"GAUGETYPE is 3: GN model - should be in hamiltonSigma"
        stop
!        call SG_GN(theta,hg)
      endif
      h=hg+hp
      if (DYNAMIC.gt.0) then
        print *,'DYNAMIC - hamilton'
        call hamiltonFermion(hf)
        h=h+hf
      endif

      return
      end              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SGqed3(theta,hg)
c      use arraysizes
c      use numbers
c      use options
      implicit none
      real(prc) theta(Nv,3)
      real(prc) hg
      integer mu,nu,i
      real(prc) Sigma(Nv),a(Nv)

      hg=zero
      do mu=1,3
        do nu=1,mu-1
 
          do i=1,Nv
           Sigma(i)=theta(iu(i,mu),nu)-theta(iu(i,nu),mu)-theta(i,nu)
          end do
          a=theta(:,mu)+Sigma
          hg=hg+sum(a*a)

        end do
      end do
      hg=gbeta*hg/2
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SGthirring3(theta,hg)
      implicit none
      real(prc) theta(Nv,3)
      real(prc) hg

      hg=Nferms*gbeta*sum(theta*theta)

      return
      end subroutine SGthirring3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine hamiltonFermion(hf)
      use basicdiracopsmod
      use gaugefield
      use dynamicmodule
      use axbmodule1
      implicit none
      real(prc) hf
      procedure(),pointer :: Mptr
      complex(prc) TMP(Nv,4)

      print *,'hamiltonFermion (dynamic)'
      Mptr=>DWilson
      call IMdagM(phi,TMP,u,.false.,baremass,czero,Mptr)
      hf=sum(conjg(phi)*TMP)

      return
      end subroutine hamiltonFermion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine march(dH,thetat)
c      use arraysizes
c      use numbers
c      use rvmodule
c      use options
      use dynamicmodule
      implicit none
      real(prc) dH,H0,H1
      real(prc) thetat(Nv,3)
      real(prc) dSdpi(Nv,3)
      real(prc) dt,hg,hp,hf
      real(prc) pp(Nv,3),ps(Nv,4)
      real(prc) etime,proby,ytest,avsteps
      integer mu,ts,tsmax

      etime=one
      tsmax=10000
      dt=0.001
      avsteps=etime/dt
      proby=one/avsteps

!     set eta for dynamic component
      if(DYNAMIC.gt.0) then
        call setDynamicNoise()
      end if

!     set heatbath (initialise) for p 
      do mu=1,3
        call gaussp(ps)
        pp(:,mu)=ps(:,1)
      enddo
c      print *,pp

!     half time-step for pp
      call hamilton(thetat,pp,H0,hg,hp,hf)
c      print *,H0,hg,hp
      call force(thetat,dSdpi)
      pp=pp-dt/2*dSdpi

c      print *,dSdpi

!     hybrid loop
      proby=1.0/avsteps
      do ts=1,tsmax
c        print *,dt*pp
        thetat=thetat+dt*pp
        call force(thetat,dSdpi)
c      print *,dSdpi
        ytest=urv()
c        print *,"ytest:",ytest,"proby:",proby
        if (ytest.lt.proby) then
          pp=pp-dt/2*dSdpi
          goto 501
        else
          pp=pp-dt*dSdpi
        endif
        call hamilton(thetat,pp,H1,hg,hp,hf)
c        print *,H0,H1,hg/Nv,hp/Nv
      end do

501   continue
      call hamilton(thetat,pp,H1,hg,hp,hf)
      print *,hg/Nv,hp/Nv
      dH=H0-H1
      return
      end subroutine march
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine accept(dH,thetat)
c      use arraysizes
c      use rvmodule
      use gaugefield
      implicit none
      real(prc) dH
      real(prc) thetat(Nv,3)
      real(prc) y
      real x
c**********************************************************************
c  Monte Carlo step: accept new fields with probability=
c              min(1,exp(H0-H1))
c**********************************************************************
      y=exp(dH)
      x=urv()
      if (x.lt.y) then
        theta=thetat
      end if
      return
      end subroutine accept
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module gaugemodule
