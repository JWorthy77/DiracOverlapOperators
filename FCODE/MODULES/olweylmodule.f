      module olweylmod
      use pacc
      use arraysizes
      use numbers
      use options
      use ratfuncs
      use weylmod
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLWeyl(R,DR,u,DAGGER,mass,SRF)
      implicit none
!     calculates DR = DOl*R for 2 component fermion
!     DOl = (1+m)/2+(1-m)/2.V(Dw) although mass should be 0
      complex(prc),intent(in) :: R(Nv,2)
      complex(prc),intent(out) :: DR(Nv,2)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) :: TMP(Nv,2),VR(Nv,2),DL(Nv,2)

      call VOLWeyl(R,DR,u,DAGGER,-MDW,SRF)
      DR=(one+mass)/two*R + (one-mass)/two*DR

      return
      end subroutine DOLWeyl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLWeyl(RR,S,u,DAGGER,dwmass,SRF)
!     approximate Voverlap using partial fraction rational functions
      implicit none
      complex(prc),intent(in) :: RR(Nv,2)
      complex(prc),intent(out) :: S(Nv,2)
      complex(prc),intent(in) ::  u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,2),TMP2(Nv,2)
      integer j,p,nn,nd,nf
c      procedure(),pointer :: Mptr => NULL()
      real(prc) mult

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      nf=nn-nd+1
      S=zero
      if (nf.ge.1) then ! it should only ever be 0 or 1
        do j=0,nf-1
          mult=front(j+1)
!         T = m.DWilson^j.RR
          TMP1=mult*RR
          do p=1,j
            call WdagW(TMP1,TMP2,u,DAGGER,dwmass,czero)
            TMP1=TMP2
          end do
          S=S+TMP1
        end do
      end if

c      Mptr => WdagWpC
      do j=1,nd
         add = -denom(j)
         call IMWeyl(RR,TMP1,u,DAGGER,dwmass,add)
         TMP2 = pf(j)*TMP1
         S=S+TMP2
      end do

      call DWeyl(S,TMP1,u,DAGGER,dwmass,czero)
      S=SRF%mult*TMP1

      end associate
      return
      end subroutine VOLWeyl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMWeyl(RR,DR,u,DAGGER,mass,add)
!     solve Mptr.DR = RR
!     requires M to be symmetric, positive definite
      implicit none
      integer,parameter :: kferm = 2*Nv
      integer,parameter :: niterc=10*Nv
      complex(prc),intent(in) :: RR(kferm)
      complex(prc),intent(out) :: DR(kferm)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call WdagWpC(RR,x2,u,DAGGER,mass,add)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call WdagWpC(p,x2,u,DAGGER,mass,add)
        alphan=sum(conjg(r)*r)
        alphad=sum(conjg(p)*x2)
        alpha=alphan/alphad

        DR=DR+alpha*p
        r=r-alpha*x2
        betan=sum(conjg(r)*r)
        if (betan.lt.resid) goto 8
        beta=betan/alphan
        p=r+beta*p
      end do

8     continue
c      print *,itercg,niterc,betan
      return
      end subroutine IMWeyl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLWeyl(RR,DR,u,DAGGER,mass,SRF) 
      use ratfuncs
      implicit none
      complex(prc),intent(in) :: RR(Nv,2)
      complex(prc),intent(out) :: DR(Nv,2)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) TMP(Nv,2)

      if (.not. DAGGER) then
        call DOLWeyl(RR,TMP,u,.not.DAGGER,mass,SRF)
        call IDOLWdagDOLW(TMP,DR,u,DAGGER,mass,SRF)
      elseif (DAGGER) then
        call IDOLWdagDOLW(RR,TMP,u,DAGGER,mass,SRF)
        call DOLWeyl(TMP,DR,u,.not.DAGGER,mass,SRF)
      end if

      return
      end subroutine IDOLWeyl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLWdagDOLW(RR,DR,u,DAGGER,mass,SRF)
      use ratfuncs
      implicit none
      integer,parameter :: kferm = 2*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc),intent(in) :: RR(kferm)
      complex(prc),intent(out) :: DR(kferm)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc),intent(in) :: SRF

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call DOLWeyl(RR,x1,u,DAGGER,mass,SRF)
      call DOLWeyl(x1,x2,u,.not.DAGGER,mass,SRF)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call DOLWeyl(p,x1,u,DAGGER,mass,SRF)
        call DOLWeyl(x1,x2,u,.not.DAGGER,mass,SRF)
        alphan=sum(conjg(r)*r)
        alphad=sum(conjg(p)*x2)
        alpha=alphan/alphad

        DR=DR+alpha*p
        r=r-alpha*x2
        betan=sum(conjg(r)*r)
        if (betan.lt.resid) goto 8
        beta=betan/alphan
        p=r+beta*p
      end do

8     continue
      print *,itercg,niterc,betan
      return
      end subroutine IDOLWdagDOLW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testolweyl()
      use rvmodule
      use gaugefield
      implicit none
      complex(prc) RR(Nv,2),DR(Nv,2),TMP(Nv,2),TMP2(Nv,2)
      type(sgnratfunc) :: SRF
      integer l
     
      call setRVs(2*Nv,RR)

      call setHTcoeffs(10,SRF)

      call DOLWeyl(RR,DR,u,.true.,zero,SRF)
      do l=2,6
        call setHTcoeffs(10*l,SRF)
        call DOLWeyl(RR,TMP2,u,.true.,zero,SRF)
        print *,"l:",maxval(abs(TMP-TMP2))
        TMP=TMP2
      end do

      call DOLWeyl(RR,TMP,u,.false.,zero,SRF)
      call DOLWeyl(TMP,TMP2,u,.true.,zero,SRF)
      call IDOLWdagDOLW(TMP2,DR,u,.false.,zero,SRF)
      print *,"DOLWdag.DOLW.IDOLWdagDOLW:",maxval(abs(DR-RR))

      call setHTcoeffs(40,SRF)
      call DOLWeyl(RR,TMP,u,.false.,zero,SRF)
      call IDOLWeyl(TMP,DR,u,.false.,zero,SRF)
      print *,"OLW.IOLW:",maxval(abs(DR-RR))

      return
      end subroutine testolweyl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine estimateWeylOLExtrema(Nmax,MAXMIN,eig)
!     the basic power method
      use rvmodule
      use gaugefield
      use utilsmod
      implicit none
      integer,intent(in) :: Nmax,MAXMIN
      real(prc) eig
      complex(prc),dimension(Nv,2) :: R,TMP,DR
      type(sgnratfunc) :: SRF
      integer i

      call setRVs(Nv*2,R)
      call normalise2(R)
      call setHTcoeffs(30,SRF)
      open(unit=81,file='EigDetails.dat',access='append',
     &           status='unknown',form='formatted')
      do i=1,Nmax
        if (MAXMIN.eq.1) then
          call DOLWeyl(R,TMP,u,.false.,zero,SRF)
          call DOLWeyl(TMP,DR,u,.true.,zero,SRF)
        elseif (MAXMIN.eq.-1) then
          call IDOLWdagDOLW(R,DR,u,.false.,zero,SRF)
        endif
        R=DR
        eig=mag2(R)
        R=R/eig
        print *,i,eig
        write(81,*) i,eig
      end do
      close(81)
      print *,"eig:",eig
      eig=sqrt(eig)
      print *,"eig:",eig
      if (MAXMIN.eq.-1) then
        eig=one/eig
      endif
      print *,"eig:",eig
      print *,"Lam:",eig/sqrt(1-eig*eig)

      return
      end subroutine estimateWeylOLExtrema
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcWOLExtrema()
      use options
      use gaugemodule
      use statsmod
      implicit none
      integer,parameter :: Nav=10
      real(prc) eig,eigs(Nav),eigav,eigsd
      real(prc) Lam,Lams(Nav),Lamav,Lamsd
      integer l,Nl,j
     
      Nl=15
      do l=-Nl,Nl
        open(unit=11,file='EigsAP3.dat',access='append',
     &           status='unknown',form='formatted')
        gbeta=sqrt(sqrt(two))**(l)
        do j=1,Nav
c          call makeQuenchedGaussianThirringField()
          call makeQuenchedCosineThirringField()
          call estimateWeylOLExtrema(20,-1,eig)
          eigs(j)=eig
          Lams(j)=eig/sqrt(one-eig*eig)
        end do
        call calcVarReal(Nav,eigs,eigav,eigsd) 
        call calcVarReal(Nav,Lams,Lamav,Lamsd) 

        write(11,*) gbeta,eigav,Lamav
        print *,gbeta,eigav,Lamav,eigsd,Lamsd,Nav
        close(11)
      end do

      return
      end subroutine calcWOLExtrema
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module olweylmod
