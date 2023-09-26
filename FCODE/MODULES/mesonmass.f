      module mesonmassmodule
      use pacc
      use arraysizes
      use numbers
      use indices
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalMesonMass()
!     calculate condensate 
      use options
      use ratfuncs
      use gammas
      use gaugefield
      use gaugemodule
      use IOmodule
      use basicdiracopsmod
      use domainwallmod
      use overlapmoduledev
      implicit none
      integer,parameter :: lmax=1
      complex(prc) :: deltaA(Nv,4),deltaDDW(Nv,4,Ls)
      complex(prc) :: DR(Nv,4),DRDDW(Nv,4,Ls)
      complex(prc) :: DR4(Nv,4,4)
      complex(prc) :: A1(4,4),A2(4,4)
      complex(prc) :: vec(4)
      complex dp,trace(Nt,lmax)
      complex trav(Nt)
      type(sgnratfunc) :: SRF
      integer xt,x1,x2,l,ix,iy
      integer idxT,idxT0,idirac,g5idx
      character(len=80) fname 

      call setHTcoeffs(20,SRF)
      trace=czero
      trav=czero
      open(unit=121,file='meson.dat',status='unknown')
!     choose time nt, and get slice at (x1,x2,nt)
c      xt=Nt/2
c      do xt=1,Nt
      do l=1,lmax
        print *,l,"of",lmax
!       get theta field
        call readConvertedThetaFileName(l,fname)
        call readConvertedThirringGaugeField(fname,theta)
        call coef(u,theta)
        print *,"read gauge field ",l*10

!       get indices at t=T
        call ia(1,1,1,idxT0)
!       loop over dirac indices
        do idirac=1,4
          deltaA=zero
          deltaA(idxT0,idirac)=cone
!          call IDW(deltaA,DR,u,.false.,baremass,czero)
          call IKDDW4(deltaA,DR,u,.false.,baremass)
!          call IDOverlap(deltaA,DR,u,.false.,baremass,SRF)
c          DR4(:,:,idirac)=one/(one-baremass)*DR
          DR4(:,:,idirac)=DR-deltaA
        end do
!            call IKDDW(deltaDDW,DRDDW,u,.false.,baremass,czero)

        do xt=1,Nt
          do ix=1,Ns
            do iy=1,Ns
          call ia(ix,iy,xt,idxT)
          A1=czero
          do idirac=1,4
c            g5idx=gamin(5,idirac)
c            A1(g5idx,1:4)=gamval(5,idirac)*DR4(idxT,1:4,idirac)
            A1(idirac,1:4)=DR4(idxT,1:4,idirac)
          end do ! idirac
c          print *,"A1"
c          print *,real(A1(1,1:4))
c          print *,real(A1(2,1:4))
c          print *,real(A1(3,1:4))
c          print *,real(A1(4,1:4))
c          A2=matmul(transpose(conjg(A1)),A1)
          A2=matmul((conjg(A1)),A1)
c          print *,"A2"
c          print *,real(A2(1,1:4))
c          print *,real(A2(2,1:4))
c          print *,real(A2(3,1:4))
c          print *,real(A2(4,1:4))

          do idirac=1,4
            trace(xt,l)=trace(xt,l)-A2(idirac,idirac)
          end do ! idirac

            end do
          end do

c          do idirac=1,4
c            dp=dot_product(conjg(A1(idirac,1:4)),A1(1:4,idirac))
c            trace(xt,l)=trace(xt,l)-dp
c          end do ! idirac
        trav(xt)=trav(xt)+trace(xt,l)
        end do ! xt
      end do ! l

      print *,"trace:",trace
      trav=trav/lmax
      do xt=1,Nt
        write(121,*) real(trace(xt,1)),aimag(trace(xt,1)),
     &                         real(trav(xt)),aimag(trav(xt))
      end do
      close(121)

      print *,real(trace(1:Nt,1))
      print *,real(trace(1:Nt,lmax))

      print *,"Ns",Ns

      return
      end subroutine evalMesonMass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module mesonmassmodule
