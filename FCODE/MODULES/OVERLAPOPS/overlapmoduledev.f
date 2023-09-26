      module overlapmoduledev
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use axbmodule1
      use basicdiracopsmod
      use VOLmodule
      use VOLNKmodule
      use VOLarchive
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLop(R,DR,u,DAGGER,mass,SRF)
      use options
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      type(sgnratfunc),intent(in) :: SRF

      if (dwkernel.eq.1) then
        call DOLS(R,DR,u,DAGGER,mass,SRF)
      elseif (dwkernel.eq.2) then
        call DOverlap(R,DR,u,DAGGER,mass,SRF)
      elseif (dwkernel.eq.3) then
        call DOverlap(R,DR,u,DAGGER,mass,SRF)
      endif

      return
      end subroutine DOLop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLop(R,DR,u,DAGGER,mass,SRF)
      use options
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      type(sgnratfunc),intent(in) :: SRF

      if (dwkernel.eq.1) then
        call IDOLS(R,DR,u,DAGGER,mass,SRF)
      elseif (dwkernel.eq.2) then
        call IDOverlap(R,DR,u,DAGGER,mass,SRF)
      endif

      return
      end subroutine IDOLop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLop2(R,DR,u,DAGGER,mass,SRF)
      use options
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Vptr ! VOLGfac, VOLGpf, VOLGrem, VOLMpf
      procedure(),pointer :: Mptr ! DdagDpC, SdagSpC, H2pC
      procedure(),pointer :: Dptr ! DWilson, DShamir, DHermWilson

      Vptr=>VOLMpf
      if (dwkernel.eq.1) then
        Mptr=>SdagSpC
        Dptr=>DShamir
      elseif (dwkernel.eq.2) then
        Mptr=>DdagDpC
        Dptr=>DWilson
      endif
      call DOL(R,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)

      return
      end subroutine DOLop2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLop2(R,DR,u,DAGGER,mass,SRF)
      use options
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Vptr ! VOLGfac, VOLGpf, VOLGrem, VOLMpf
      procedure(),pointer :: Mptr ! DdagDpC, SdagSpC, H2pC
      procedure(),pointer :: Dptr ! DWilson, DShamir, DHermWilson

      Vptr=>VOLMpf
      if (dwkernel.eq.1) then
        Mptr=>SdagSpC
        Dptr=>DShamir
      elseif (dwkernel.eq.2) then
        Mptr=>DdagDpC
        Dptr=>DWilson
      endif
      call IDOL(R,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)

      return
      end subroutine IDOLop2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLop4(R,DR,u,DAGGER,mass,SRF1,SRF2,SRF3)
      use options
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      type(sgnratfunc),intent(in) :: SRF1,SRF2,SRF3
      procedure(),pointer :: Vptr ! VOLGfac, VOLGpf, VOLGrem, VOLMpf
      procedure(),pointer :: Mptr ! DdagDpC, SdagSpC, H2pC
      procedure(),pointer :: Dptr ! DWilson, DShamir, DHermWilson

      Vptr=>VOLMpf
      if (dwkernel.eq.1) then
        Mptr=>SdagSpC
        Dptr=>DShamir
      elseif (dwkernel.eq.2) then
        Mptr=>DdagDpC
        Dptr=>DWilson
      endif
      call IDOLMSvar(R,DR,u,DAGGER,mass,SRF1,SRF2,SRF3,Vptr,Mptr,Dptr)

      return
      end subroutine IDOLop4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOverlap(R,DR,u,DAGGER,mass,SRF)
      use options
      use gammas
      implicit none
!     calculates DR = DOl*R 
!     DOl = (1+m)/2+(1-m)/2.V(Dw)
!     DOl3 = (1-im.g3)/2+(1+im.g3)/2.V(Dw)
!     DOl3^dag = (1+im.g3)/2 + V(Dw)^dag.(1-im.g3)/2
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Sptr => NULL()
      complex(prc) :: TMP(Nv,4),VR(Nv,4),DL(Nv,4)

      if (RFTYPE.eq.1) then
        Sptr => VOLpf
      elseif (RFTYPE.eq.2) then
        Sptr => VOLfac
      elseif (RFTYPE.eq.3) then
        Sptr => VOLpoly
      endif

      if (MTYPE.eq.1) then ! DOL = (1+m) + (1-m).V

        call Sptr(R,DR,u,DAGGER,-MDW,SRF)
        DR=(one+mass)/two*R + (one-mass)/two*DR

      elseif (MTYPE.eq.3) then ! DOL = (1+i.m.g3) + (1-i.m.g3).V
!       this is not the form corresponding to the domain wall formulation MTYPE=3
        if (.not.DAGGER) then
          call Sptr(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          DR=zi*mass*(R-VR)/two
          call mGmu(DR,4)
          DR=DL+DR
        elseif (DAGGER) then
          call Sptr(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          TMP=R
          call mGmu(TMP,4)
          call Sptr(TMP,VR,u,DAGGER,-MDW,SRF)
          DR=-zi*mass*(TMP-VR)/two
          DR=DL+DR
        endif

      elseif (MTYPE.eq.4) then 

        if (.not.DAGGER) then ! DOL=(1+im.g3) + V(1-im.g3)
          TMP=R
          call mGmu(TMP,4)
          DL=R+zi*mass*TMP
          VR=R-zi*mass*TMP
          call Sptr(VR,DR,u,DAGGER,-MDW,SRF)
          DR=(DL+DR)/two
        elseif (DAGGER) then ! DOLdag=(1-im.g3) + (1+im.g3)Vdag
          TMP=R
          call mGmu(TMP,4)
          DL=R-zi*mass*TMP
          call Sptr(R,VR,u,DAGGER,-MDW,SRF)
          TMP=VR
          call mGmu(TMP,4)
          DR=VR+zi*mass*TMP
          DR=(DL+DR)/two
        endif

      elseif (MTYPE.eq.5) then

        if (.not.DAGGER) then ! DOL=(1-im.g3) + V(1+im.g3)
          TMP=R
          call mGmu(TMP,4)
          DL=R-zi*mass*TMP
          VR=R+zi*mass*TMP
          call Sptr(VR,DR,u,DAGGER,-MDW,SRF)
          DR=(DL+DR)/two
        elseif (DAGGER) then ! DOLdag=(1+im.g3) + (1-im.g3)Vdag
          TMP=R
          call mGmu(TMP,4)
          DL=R+zi*mass*TMP
          call Sptr(R,VR,u,DAGGER,-MDW,SRF)
          TMP=VR
          call mGmu(TMP,4)
          DR=VR-zi*mass*TMP
          DR=(DL+DR)/two
       endif


      endif

      return
      end subroutine DOverlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOverlap(R,DR,u,DAGGER,mass,SRF)
      use options
      use axbmodule2
      implicit none
!     solves DOL.DR = R 
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()
      type(ioptions)  iopts

      iopts%SRF=SRF
      iopts%mass=mass
      Mptr => DOverlap
      call IMnonsymOpts(R,DR,u,DAGGER,Mptr,iopts)

      return
      end subroutine IDOverlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLS(R,DR,u,DAGGER,mass,SRF)
      use options
      use ratfuncs
      use axbmodule1
      use gammas
!     approximate Voverlap for Shamir kernel using partial fraction rational functions
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Sptr => NULL()
      complex(prc) :: TMP(Nv,4),VR(Nv,4),DL(Nv,4),IMG3R(Nv,4)

      if (VERBOSE.gt.2) print *,'DOLS'
      Sptr=>VOLSpf
      if (MTYPE.eq.1) then

        call Sptr(R,TMP,u,DAGGER,-MDW,SRF)
        DR=(half+mass/two)*R+(half-mass/two)*TMP

      elseif (MTYPE.eq.3) then ! DOL = (1+i.m.g3) + (1-i.m.g3).V

        if (.not.DAGGER) then
          call Sptr(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          DR=zi*mass*(R-VR)/two
          call mGmu(DR,4)
          DR=DL+DR
        elseif (DAGGER) then
          call Sptr(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          TMP=R
          call mGmu(TMP,4)
          call Sptr(TMP,VR,u,DAGGER,-MDW,SRF)
          DR=-zi*mass*(TMP-VR)/two
          DR=DL+DR
        endif

      elseif (MTYPE.eq.4) then 

        if (.not.DAGGER) then ! DOL=(1+im.g3) + V(1-im.g3)
          TMP=R
          call mGmu(TMP,4)
          DL=R+zi*mass*TMP
          VR=R-zi*mass*TMP
          call Sptr(VR,DR,u,DAGGER,-MDW,SRF)
          DR=(DL+DR)/two
        elseif (DAGGER) then ! DOLdag=(1-im.g3) + (1+im.g3)Vdag
          TMP=R
          call mGmu(TMP,4)
          DL=R-zi*mass*TMP
          call Sptr(R,VR,u,DAGGER,-MDW,SRF)
          TMP=VR
          call mGmu(TMP,4)
          DR=VR+zi*mass*TMP
          DR=(DL+DR)/two
        endif

      else

        print *,"MTYPE",MTYPE,"not supported"

      endif
      return
      end subroutine DOLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLS(R,DR,u,DAGGER,mass,SRF)
      use options
      use axbmodule2
      implicit none
!     solves DOL.DR = R 
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()
      type(ioptions)  iopts

      if (VERBOSE.eq.2) then ; print *,"IDOLS" ; endif
      iopts%SRF=SRF
      iopts%mass=mass
      Mptr => DOLS
      call IMnonsymOpts(R,DR,u,DAGGER,Mptr,iopts)

      return
      end subroutine IDOLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOL(R,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      use options
      use gammas
      implicit none
!     calculates DR = DOl*R 
!     DOl = (1+m)/2+(1-m)/2.V(Dw)
!     DOl3 = (1-im.g3)/2+(1+im.g3)/2.V(Dw)
!     DOl3^dag = (1+im.g3)/2 + V(Dw)^dag.(1-im.g3)/2
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer,intent(in) :: Vptr ! VOLGfac, VOLGpf, VOLGrem, VOLMpf
      procedure(),pointer,intent(in) :: Mptr ! DdagDpC, SdagSpC, H2pC
      procedure(),pointer,intent(in) :: Dptr ! DWilson, DShamir, DHermWilson
      complex(prc) :: TMP(Nv,4),VR(Nv,4),DL(Nv,4)

      if (MTYPE.eq.1) then

        call Vptr(R,DR,u,DAGGER,-MDW,SRF,Mptr,Dptr)
        DR=(one+mass)/two*R + (one-mass)/two*DR

      elseif (MTYPE.eq.3) then

        if (.not.DAGGER) then
          call Vptr(R,VR,u,DAGGER,-MDW,SRF,Mptr,Dptr)
          DL=(R+VR)/two
          DR=zi*mass*(R-VR)/two
          call mGmu(DR,4)
          DR=DL+DR
        elseif (DAGGER) then
          call Vptr(R,VR,u,DAGGER,-MDW,SRF,Mptr,Dptr)
          DL=(R+VR)/two
          TMP=R
          call mGmu(TMP,4)
          call Vptr(TMP,VR,u,DAGGER,-MDW,SRF,Mptr,Dptr)
          DR=-zi*mass*(TMP-VR)/two
          DR=DL+DR
        endif

      else
        print *,"MTYPE not implemented"
        stop
      endif

      return
      end subroutine DOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOL(RR,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr) 
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Vptr,Mptr,Dptr
      complex(prc) TMP(Nv,4)

      if (.not. DAGGER) then
        call DOL(RR,TMP,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        call IDOLdagDOL(TMP,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      elseif (DAGGER) then
        call IDOLdagDOL(RR,TMP,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        call DOL(TMP,DR,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      end if

      return
      end subroutine IDOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLdagDOL(RR,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Vptr,Mptr,Dptr
      complex(prc) add

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call DOL(RR,x1,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      call DOL(x1,x2,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call DOL(p,x1,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        call DOL(x1,x2,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        alphan=sum(conjg(r)*r)
        alphad=sum(conjg(p)*x2)
        alpha=alphan/alphad

        DR=DR+alpha*p
        r=r-alpha*x2
        betan=sum(conjg(r)*r)
c        print *,i,betan
        if (betan.lt.resid) goto 8
        beta=betan/alphan
        p=r+beta*p
      end do

8     continue
!      print *,itercg,niterc,betan
      return
      end subroutine IDOLdagDOL     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLMSvar(RR,DR,u,DAGGER,mass,SRF1,SRF2,SRF3,
     &                                              Vptr,Mptr,Dptr) 
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF1,SRF2,SRF3
      procedure(),pointer :: Vptr,Mptr,Dptr
      complex(prc) TMP(Nv,4)

      if (.not. DAGGER) then
        call DOL(RR,TMP,u,.not.DAGGER,mass,SRF3,Vptr,Mptr,Dptr)
        call IOLdOLMSvar(TMP,DR,u,DAGGER,mass,SRF1,SRF2,SRF3,
     &                                               Vptr,Mptr,Dptr)
      elseif (DAGGER) then
        call IOLdOLMSvar(RR,TMP,u,DAGGER,mass,SRF1,SRF2,SRF3,
     &                                               Vptr,Mptr,Dptr)
        call DOL(TMP,DR,u,.not.DAGGER,mass,SRF3,Vptr,Mptr,Dptr)
      end if

      return
      end subroutine IDOLMSvar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IOLdOLMSvar(RR,DR,u,DAGGER,mass,SRF1,SRF2,SRF3,
     &                                           Vptr,Mptr,Dptr)
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc),target :: SRF1,SRF2,SRF3
      procedure(),pointer :: Vptr,Mptr,Dptr
      complex(prc) add

      type(sgnratfunc),pointer :: SRF
      integer ptr
      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      SRF=>SRF1
      ptr=1
      DR=RR
      call DOL(RR,x1,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      call DOL(x1,x2,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call DOL(p,x1,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        call DOL(x1,x2,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        alphan=sum(conjg(r)*r)
        alphad=sum(conjg(p)*x2)
        alpha=alphan/alphad

        DR=DR+alpha*p
        r=r-alpha*x2
        betan=sum(conjg(r)*r)
        print *,i,betan
        if (betan.lt.resid) then
          if (ptr.eq.1) then
            SRF=>SRF2
            ptr=2
          elseif (ptr.eq.2) then
            SRF=>SRF3
            ptr=3
          else
            goto 8
          endif
        endif
        beta=betan/alphan
        p=r+beta*p
      end do

8     continue
      print *,itercg,niterc,betan
      return
      end subroutine IOLdOLMSvar    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module overlapmoduledev
