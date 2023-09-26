!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module gaugefieldSU3
      use pacc
      use arraysizes
      use numbers
      use determinants
      implicit none
      integer,parameter :: Nxset=50
      real(prc) thetaSU3(3,3,Nv,3) ! turn these into a derived type
      complex(prc) uSU3(3,3,Nv,3)
      complex(prc) GellManns(3,3,8)
      complex(prc) Xset(3,3,2*Nxset) ! set of matrices to choose link update

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initSU3
      use paulimodule
      use determinants
      implicit none
      integer i,n,mu
   
      uSU3=cone
      do mu=1,3
        do n=1,Nv
          do i=1,3
            uSU3(i,i,n,mu)=cone
          end do
        end do
      end do
      call setPauliMatrices
c      print *,"sig1"
c      print *,sigPauli(:,:,1)
c      print *,"sig2"
c      print *,sigPauli(:,:,2)
c      print *,"sig3"
c      print *,sigPauli(:,:,3)
      call createXset()
      do i=1,Nxset
        print *,"i",i
        print *,Xset(:,:,i)
        print *,det3(Xset(:,:,i))
        print *,det3(Xset(:,:,i+Nxset))
      end do
      return
      end subroutine initSU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine metropolisSU3(u)
      use rvmodule
      implicit none
      complex(prc) :: u(3,3,Nv,3)
      real(prc) SG0,Sl0,SlT,DS
      integer m,Nm,idx
      integer n,mu
      integer nkeep
      real(prc) rv,pkeep
      complex(prc) :: U0link(3,3),UTlink(3,3)

      call SG_SU3(u,SG0)
      Nm=100*Nv*3
c      Nm=100
      nkeep=0
      do m=1,Nm
        ! choose link for potential update
        rv=urv()
        n=ceiling(Nv*rv)
        mu=ceiling(3*rv)
        ! store current link
        U0link=u(:,:,n,mu)
        ! calc local link action
        Sl0=linkAction(u,n,mu)
        ! choose X for change to link
        rv=urv()
        idx=ceiling(2*Nxset*rv)
        UTlink=matmul(Xset(:,:,idx),u(:,:,n,mu))
        u(:,:,n,mu)=UTlink
        ! calc local test link action
        SlT=linkAction(u,n,mu)
        ! test for replacement
        DS=SlT-Sl0
        pkeep=exp(-DS)
c        print *,m,pkeep
        rv=urv()
        nkeep=nkeep+1
        if (rv > pkeep) then
        ! revert back to initial link
          u(:,:,n,mu)=U0link
          nkeep=nkeep-1
        end if
      end do
      print *,"kept",nkeep,"of",Nm
      return
      end subroutine metropolisSU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function linkAction(u,n,mu)
      use options
      implicit none
      complex(prc) u(3,3,Nv,3)
      integer n,mu
      real(prc) linkAction
      complex(prc) A(3,3),B(3,3)
      integer nu1,nu2
   
      nu1=2
      nu2=3
      if (mu.eq.2) then
        nu1=1
      elseif (mu.eq.3) then
        nu1=1
        nu2=2
      endif
      A=staples(u,n,mu,nu1)+staples(u,n,mu,nu2)
      B=-matmul(u(:,:,n,mu),A)

      linkAction=real((B(1,1))+B(2,2)+B(3,3),prc)+4*one
      linkAction=linkAction*gbeta/Nv
      return
      end function linkAction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function staples(u,n,mu,nu)
      use indices
      implicit none
      complex(prc) u(3,3,Nv,3)
      integer n,mu,nu
      complex(prc),dimension(3,3) :: staples
      complex(prc),dimension(3,3) :: TMP1,TMP2,TMP3,TMP4
      integer npmu,npnu,nmnu,npmumnu

      npmu=iu(n,mu)
      npnu=iu(n,nu)
      nmnu=id(n,nu)
      npmumnu=id(npmu,nu)

      TMP1=u(:,:,npmu,nu)
      TMP2=conjg(transpose(u(:,:,npnu,mu)))
      TMP3=conjg(transpose(u(:,:,n,nu)))
      TMP4=matmul(TMP1,TMP2)
      staples=matmul(TMP4,TMP3)

      TMP1=conjg(transpose(u(:,:,npmumnu,nu)))
      TMP2=conjg(transpose(u(:,:,nmnu,mu)))
      TMP3=u(:,:,nmnu,nu)
      TMP4=matmul(TMP1,TMP2)
      staples=staples+matmul(TMP4,TMP4)

      return
      end function staples
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setGellManns()
      implicit none

      GellManns=czero
      GellManns(1,2,1)=one
      GellManns(2,1,1)=one
      GellManns(1,2,2)=-zi
      GellManns(2,1,2)=zi
      GellManns(1,1,3)=one
      GellManns(2,2,3)=-one
      GellManns(1,3,4)=one
      GellManns(3,1,4)=one
      GellManns(1,3,5)=-zi
      GellManns(3,1,5)=zi
      GellManns(2,3,6)=one
      GellManns(3,2,6)=one
      GellManns(2,3,7)=-zi
      GellManns(3,2,7)=zi
      GellManns(1,1,8)=one/rt3
      GellManns(2,2,8)=one/rt3
      GellManns(3,3,8)=-two/rt3
      return
      end subroutine setGellManns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function matrixSU2(eps)
      use rvmodule
      use paulimodule
      implicit none
      real(prc) eps
      complex(prc) matrixSU2(2,2)
      real(prc) x0,xvec(3)
      real(prc) rv0,rv(3)
      real(prc) lr
      integer i

      rv0=urv()-half
      x0=rv0/abs(rv0)*sqrt(one-eps*eps)
      do i=1,3
        rv(i)=urv()-half
      end do
      lr=norm2(rv)
      do i=1,3
        xvec(i)=eps*rv(i)/lr
      end do

      matrixSU2=czero
      matrixSU2(1,1)=x0
      matrixSU2(2,2)=x0
      do i=1,3
        matrixSU2=matrixSU2+zi*xvec(i)*sigPauli(:,:,i)
      end do

c      print *,"det2",det2(matrixSU2)

      return
      end function matrixSU2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function matrixSU3(eps)
      implicit none
      real(prc) eps
      complex(prc),dimension(3,3) :: matrixSU3
      complex(prc),dimension(3,3) :: R,S,T,M
      complex(prc),dimension(2,2) :: TMP
      integer i,j,i1,j1

      TMP=matrixSU2(eps)
      R=czero
      R(1:2,1:2)=TMP
      R(3,3)=cone
c      print *,"R",det3(R)
      TMP=matrixSU2(eps)
      S=czero
      S(1,1)=TMP(1,1)
      S(1,3)=TMP(1,2)
      S(3,1)=TMP(2,1)
      S(3,3)=TMP(2,2)
      S(2,2)=cone
c      print *,"S",det3(S)
      TMP=matrixSU2(eps)
      T=czero
      T(2:3,2:3)=TMP
      T(1,1)=cone
c      print *,"T",det3(T)
      
      M=matmul(S,T)
      matrixSU3=matmul(R,M)
 
c      print *,"SU3",det3(matrixSU3)
     
      return
      end function matrixSU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine reconSU3(S)
      implicit none
!     adjust matrix S to ensure SU3 properties
!       (    u    )
!     U=(    v    )  , det(U)=1
!       ( u* x v* )
      complex(prc) S(3,3),u(3),v(3)
      complex(prc) m
      real(prc) lu,lv

      u=S(1,:)
      lu=sqrt(dot_product(u,u))
      u=u/lu
      v=S(2,:)
      m=dot_product(v,u)
      v=v-m*u
      lv=sqrt(dot_product(v,v))
      v=v/lv
      S(1,:)=u
      S(2,:)=v
      S(3,:)=cross(conjg(u),conjg(v))

      return
      end subroutine reconSU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function cross(a,b)
      implicit none
      complex(prc) a(3),b(3)
      complex(prc) cross(3)
  
      cross(1)=a(2)*b(3)-a(3)*b(2)
      cross(2)=a(3)*b(1)-a(1)*b(3)
      cross(3)=a(1)*b(2)-a(2)*b(1)

      return
      end function cross
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine createXset()
      implicit none
      complex(prc) X(3,3)
      real(prc) eps
      integer i

      eps=0.1
      do i=1,Nxset 
        X=matrixSU3(eps)
        Xset(:,:,i)=X
        Xset(:,:,NXset+i)=conjg(transpose(X)) 
      end do
      return
      end subroutine createXset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SG_SU3(u,SG)
      use arraysizes
      use numbers
      use options
      implicit none
      complex(prc) u(3,3,Nv,3)
      real(prc) SG
      integer mu,nu,i,t
      complex(prc) :: plaq(3,3)

      SG=zero
      do mu=1,3
        do nu=1,mu-1
          do i=1,Nv
            plaq=plaquette(u,i,mu,nu)
            do t=1,3
              SG=SG+real(one-plaq(t,t))
            end do   
          end do
        end do
      end do
      SG=gbeta*SG/3
      return
      end subroutine SG_SU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function plaquette(uSU3,n,mu,nu)
      use indices
      implicit none
      complex(prc),intent(in) :: uSU3(3,3,Nv,3)
      integer :: n,mu,nu
      complex(prc),dimension(3,3) :: plaquette
      complex(prc),dimension(3,3) :: TMP,TMP2
      integer npmu,npnu     

      npmu=iu(n,mu)
      npnu=iu(n,nu)

      TMP=uSU3(:,:,n,mu)
      TMP2=uSU3(:,:,npmu,nu)
      plaquette=matmul(TMP,TMP2)
      TMP=conjg(transpose(uSU3(:,:,npnu,mu)))
      TMP2=matmul(plaquette,TMP)
      TMP=conjg(transpose(uSU3(:,:,n,nu)))
      plaquette=matmul(TMP2,TMP)
      
      return
      end function plaquette
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writeGaugeFieldSU3(fname)
      implicit none
      character*1024 fname

      open(unit=10,file=fname,status='unknown')
      write(10,*) Nv
      write(10,*) thetaSU3
      write(10,*) uSU3
      close(10)
      return
      end subroutine writeGaugeFieldSU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gaugeTransformSU3(u,alpha)
!     local gauge transform of U, by alpha
      use pacc
      use arraysizes
      use indices
      use numbers
      implicit none
      complex(prc) u(Nv,3)
      real(prc) alpha(Nv)
      integer i,mu

      do mu=1,3
        do i=1,Nv
          u(i,mu)=exp(cmplx(zero,alpha(i)))*u(i,mu)*
     &                 exp(-cmplx(zero,alpha(iu(i,mu))))
        enddo
      enddo
      return
      end subroutine gaugeTransformSU3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module gaugefieldSU3
