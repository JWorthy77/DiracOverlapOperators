      module mesonmodule
      use pacc
      use arraysizes
      use numbers
      use options
      use indices
      use ratfuncs
      use gaugefield
      use gaugemodule
      use overlapmoduleprod
      use domainwallmod
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcColDm1(ix,iy,it,id,DAGGER,mass,SRF,DR)
      implicit none
      integer ix,iy,it,id
      logical DAGGER
      real(prc) :: mass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc),intent(out) :: DR(Nv,4)
      integer iv
      complex(prc) :: R(Nv,4)

      call ia(ix,iy,it,iv)
      R=zero ; R(iv,id)=one
      if (dwkernel.eq.1) then
        call IKDDW4(R,DR,u,DAGGER,mass)
      elseif (dwkernel.eq.2) then
        call IDOLop3(R,DR,u,DAGGER,mass,SRF)
      endif

      return
      end subroutine calcColDm1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcOLColumnPropagators(ix1,iy1,it1,DAGGER, 
     &                            mass,SRF,props)
!     <psi(ix1,iy1,it1)barpsi(ix2,iy2,it2)> for all column i
      implicit none
      integer ix1,iy1,it1
      logical DAGGER
      real(prc) :: mass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) :: DR(Nv,4)
      complex(prc) :: props(Nv,4,4)
      integer id,iv
      
      do id=1,4
        call calcColDm1(ix1,iy1,it1,id,DAGGER,mass,SRF,DR)
        do iv=1,Nv
          props(iv,id,1)=DR(iv,1)
          props(iv,id,2)=DR(iv,2)
          props(iv,id,3)=DR(iv,3)
          props(iv,id,4)=DR(iv,4)
        end do
      end do
      return
      end subroutine calcOLColumnPropagators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcOLColumnDaggerPropagators(pp,ppdg)
      implicit none
      complex(prc) :: pp(Nv,4,4),ppdg(Nv,4,4)
      complex(prc) :: ppr(4,4)
      integer i,j,k

      do i=1,Nv
        do j=1,4
          do k=1,4
            ppr(j,k)=conjg(pp(i,k,j))
          end do
        end do
        ppdg(i,:,:)=ppr
      end do

      return
      end subroutine calcOLColumnDaggerPropagators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcOLPropagator(ix1,iy1,it1,ix2,iy2,it2,DAGGER, 
     &                            mass,SRF,prop)
!     <psi(ix1,iy1,it1)barpsi(ix2,iy2,it2)>
      implicit none
      integer ix1,iy1,it1,ix2,iy2,it2
      logical DAGGER
      real(prc) :: mass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) :: DR(Nv,4)
      complex(prc) :: prop(4,4)
      integer id,iv
      
      do id=1,4
        call calcColDm1(ix1,iy1,it1,id,DAGGER,mass,SRF,DR)
        call ia(ix2,iy2,it2,iv)
        prop(id,1)=DR(iv,1)
        prop(id,2)=DR(iv,2)
        prop(id,3)=DR(iv,3)
        prop(id,4)=DR(iv,4)
      end do
      return
      end subroutine calcOLPropagator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcTwoPointCorrelator(ix1,iy1,it1,ix2,iy2,it2,DAGGER, 
     &                            mass,SRF,tpc)
!     <barpsi(ix2,iy2,it2)psi(ix1,iy1,it1)>
      implicit none
      integer ix1,iy1,it1,ix2,iy2,it2
      logical DAGGER
      real(prc) :: mass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) tpc
      complex(prc) :: prop(4,4)

      call calcOLPropagator(ix1,iy1,it1,ix2,iy2,it2,DAGGER, 
     &                            mass,SRF,prop)
      tpc=prop(1,1)+prop(2,2)+prop(3,3)+prop(4,4)

      return
      end subroutine calcTwoPointCorrelator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcMesonCorrelator(ix1,iy1,it1,ix2,iy2,it2,DAGGER, 
     &                               mass,SRF,mc)
!     <barpsi(iv1)psi(iv1)barpsi(iv2)psi(iv2)>
!  =  S(iv1,iv1).S(iv2,iv2) - S(iv1,iv2).S(iv2,iv1)
      implicit none
      integer ix1,iy1,it1,ix2,iy2,it2
      logical DAGGER
      real(prc) :: mass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) mc
      complex(prc) tpc

      mc=one
      call calcTwoPointCorrelator(ix1,iy1,it1,ix2,iy2,it2,DAGGER, 
     &                            mass,SRF,tpc)
      mc=mc*tpc;
      print *,"tpc:",tpc
      call calcTwoPointCorrelator(ix2,iy2,it2,ix1,iy1,it1,DAGGER, 
     &                            mass,SRF,tpc)
      print *,"tpc:",tpc
      mc=mc*tpc
      print *,"mc:",mc

      return
      end subroutine calcMesonCorrelator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcQuenchedMesonPropagator(N)
      implicit none
      integer N
      integer ix,iy,it,id
      integer ix2,iy2,it2
      complex(prc) :: DR(Nv,4)
      type(sgnratfunc) :: SRF
      real(prc) :: mass
      complex(prc) :: prop(4,4)
      complex(prc) :: mc
      integer iloop

      dwkernel=2 ! 2-Wilson
      GAUGETYPE=1
      mass=0.05
      gbeta=2

      do iloop=1,N

      if (GAUGETYPE.eq.1) then
        call makeQuenchedCosineThirringField()
      elseif (GAUGETYPE.eq.2) then
        call makeQuenchedGaussianThirringField()
      endif
      
      ix=1 ; iy=1 ; it=1 ; id=1
      call setHTcoeffs(24,SRF)
      call calcColDm1(ix,iy,it,id,.false.,mass,SRF,DR)
      print *,DR

      ix2=1 ; iy2=1 ; it2=2 ; 
      call calcOLPropagator(ix,iy,it,ix2,iy2,it2,.false.,mass,SRF,
     &                      prop)
      print *,prop

      call calcMesonCorrelator(ix,iy,it,ix2,iy2,it2,.false., 
     &                               mass,SRF,mc)
      print *,"mc:",mc
      open(unit=11,file='meson.dat',access='append',
     &           status='unknown',form='formatted')
      write(11,*) gbeta,ix2,iy2,it2,mc
      close(11)

      end do

      return
      end subroutine calcQuenchedMesonPropagator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcAxialWardRHS(iv,mass,SRF,rhs)
      use gammas
      use indices
      use ratfuncs 
      implicit none
      integer iv
      real(prc) :: mass
      type(sgnratfunc),intent(in) :: SRF
      real(prc) rhs
      complex(prc) rhspart
      integer ix,iy,it
      complex(prc),dimension(Nv,4,4) :: props,propsdg
      complex(prc),dimension(4,4) :: ppl,ppr,pp,tl,tr,tt
      integer i

      call getLatticeCoords(iv,ix,iy,it)

      call calcOLColumnPropagators(ix,iy,it,.false.,mass,SRF,props)
      call calcOLColumnDaggerPropagators(props,propsdg)
      rhs=0
      do i=1,Nv
c        print *,i
        ppl = props(i,:,:)
c        tl=ppl
c        print *,"ppl",ppl
c        call matrixGmu(ppl,5) ! g5.ppl
c        print *,"ppl",ppl
        ppr = propsdg(i,:,:)
c        tr=ppr
c        print *,"ppr",ppr
c        call matrixGmu(ppl,5) ! g5.ppr
c        print *,"ppr",ppr
        pp =  matmul(ppl,ppr) 
c        tt =  matmul(tl,tr) 
c        print *,"pp",pp
        rhspart=pp(1,1)+pp(2,2)+pp(3,3)+pp(4,4)
        print *,i,rhspart
c        rhspart=tt(1,1)+tt(2,2)+tt(3,3)+tt(4,4)
c        print *,i,rhspart
        rhs=rhs+real(rhspart)
c        print *,"pp",pp
      end do
      print *,"rhs",rhs

      return
      end subroutine calcAxialWardRHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcQuenchedAxialWard()
      implicit none
      integer Nloop
      integer iv
      type(sgnratfunc) :: SRF
c      real(prc) :: mass
      real(prc) :: rhs
      integer iloop,Nsrf
      logical ZOLO,FEXISTS
      real(prc) zmin,zmax

      print *,"read input file"
      open(unit=31,file='qaxwardopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) Nloop,dwkernel,GAUGETYPE,MTYPE,baremass,gbeta,Nsrf,
     &ZOLO,zmin,zmax
        close(31)
        print *,"calcQuenchedAxialWard"
        print *,"Nloop:",Nloop
        print *,"kernel:",dwkernel
        print *,"gaugetype:",GAUGETYPE
        print *,"masstype:",MTYPE
        print *,"baremass:",baremass
        print *,"gbeta:",gbeta
        if (ZOLO) then
          print *,"Zolotarev",zmin,zmax
        else
          print *,"HT"
        end if
        print *,"Nsrf:",Nsrf
      else
        print *,"overlap options file not found"
        stop
      endif
      close(31)
      if (ZOLO) then
        call setZoloCoeffs(Nsrf,SRF,zmin,zmax)
      else
        call setHTcoeffs(Nsrf,SRF)
      end if

      do iloop=1,Nloop

      if (GAUGETYPE.eq.1) then
        call makeQuenchedCosineThirringField()
      elseif (GAUGETYPE.eq.2) then
        call makeQuenchedGaussianThirringField()
      endif
      
      iv=1
      call calcAxialWardRHS(iv,baremass,SRF,rhs)

      open(unit=11,file='awrhs.dat',access='append',
     &           status='unknown',form='formatted')
      write(11,"(2F8.4,3I4,F14.8)") gbeta,baremass,Ns,iv,Nsrf,rhs
      close(11)

      end do
      return
      end subroutine calcQuenchedAxialWard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcAxialWard(fnum,iv)
      use IOmodule
      implicit none
      integer Nloop
      integer fnum,iv
      type(sgnratfunc) :: SRF
c      real(prc) :: mass
      real(prc) :: rhs
      integer iloop,Nsrf
      logical ZOLO,FEXISTS
      real(prc) zmin,zmax
      character(len=80) fname 

      print *,"read input file"
      open(unit=31,file='axwardopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) Nloop,dwkernel,GAUGETYPE,MTYPE,baremass,gbeta,Nsrf,
     &ZOLO,zmin,zmax
        close(31)
        print *,"calcAxialWard"
        print *,"Nloop:",Nloop
        print *,"kernel:",dwkernel
        print *,"gaugetype:",GAUGETYPE
        print *,"masstype:",MTYPE
        print *,"baremass:",baremass
        print *,"gbeta:",gbeta
        if (ZOLO) then
          print *,"Zolotarev",zmin,zmax
        else
          print *,"HT"
        end if
        print *,"Nsrf:",Nsrf
      else
        print *,"overlap options file not found"
        stop
      endif
      close(31)
      if (ZOLO) then
        call setZoloCoeffs(Nsrf,SRF,zmin,zmax)
      else
        call setHTcoeffs(Nsrf,SRF)
      end if

!      do iloop=1,Nloop

      call readThetaFileName(fnum,fname)
      call readMPIConFile(fname,theta)
      call coef(u,theta)

!      if (GAUGETYPE.eq.1) then
!        call makeQuenchedCosineThirringField()
!      elseif (GAUGETYPE.eq.2) then
!        call makeQuenchedGaussianThirringField()
!      endif
      
      call calcAxialWardRHS(iv,baremass,SRF,rhs)

      open(unit=11,file='awrhs.dat',access='append',
     &           status='unknown',form='formatted')
      write(11,"(2F8.4,3I5,F14.8)") gbeta,baremass,fnum,iv,Nsrf,rhs
      close(11)

!      end do
      return
      end subroutine calcAxialWard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module mesonmodule
