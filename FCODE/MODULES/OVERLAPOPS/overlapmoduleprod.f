      module overlapmoduleprod ! code used for main results in thesis
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use axbmodule1
      use basicdiracopsmod
      use VOLmodule
      implicit none
      logical,parameter :: VB_OL=.false.
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLop3(R,DR,u,DAGGER,mass,SRF)
      use options
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      type(sgnratfunc),intent(in) :: SRF

      if (dwkernel.eq.1) then
        call DOLMS(R,DR,u,DAGGER,mass,SRF)
      elseif (dwkernel.eq.2) then
        call DOLMW(R,DR,u,DAGGER,mass,SRF)
      endif

      return
      end subroutine DOLop3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLop3(R,DR,u,DAGGER,mass,SRF)
      use options
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      type(sgnratfunc),intent(in) :: SRF
      if(VB_OL)then ; print *,"IDOLop3" ; endif
      if (dwkernel.eq.1) then
        call IDOLMS(R,DR,u,DAGGER,mass,SRF)
      elseif (dwkernel.eq.2) then
        if(VB_OL)then ; print *,"IDOLMW" ; endif
        call IDOLMW(R,DR,u,DAGGER,mass,SRF)
      endif

      return
      end subroutine IDOLop3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLMW(R,DR,u,DAGGER,mass,SRF)
      use options
      use gammas
      implicit none
!     calculates DR = DOl*R  for Wilson kernel using multishift cg
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      complex(prc) :: TMP(Nv,4),VR(Nv,4),DL(Nv,4)

      if (VB_OL) then ; print *,"DOLMW" ; end if;

      if (MTYPE.eq.1) then ! DOL = (1+m) + (1-m).V

        call VOLMWpf(R,DR,u,DAGGER,-MDW,SRF)
        DR=(one+mass)/two*R + (one-mass)/two*DR

      elseif (MTYPE.eq.3) then ! DOL = (1+i.m.g3) + (1-i.m.g3).V
!       this is not the form corresponding to the domain wall formulation MTYPE=3
        if (.not.DAGGER) then
          call VOLMWpf(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          DR=zi*mass*(R-VR)/two
          call mGmu(DR,4)
          DR=DL+DR
        elseif (DAGGER) then
          call VOLMWpf(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          TMP=R
          call mGmu(TMP,4)
          call VOLMWpf(TMP,VR,u,DAGGER,-MDW,SRF)
          DR=-zi*mass*(TMP-VR)/two
          DR=DL+DR
        endif

      elseif (MTYPE.eq.4) then 

        if (.not.DAGGER) then ! DOL=(1+im.g3) + V(1-im.g3)
          TMP=R
          call mGmu(TMP,4)
          DL=R+zi*mass*TMP
          VR=R-zi*mass*TMP
          call VOLMWpf(VR,DR,u,DAGGER,-MDW,SRF)
          DR=(DL+DR)/two
        elseif (DAGGER) then ! DOLdag=(1-im.g3) + (1+im.g3)Vdag
          TMP=R
          call mGmu(TMP,4)
          DL=R-zi*mass*TMP
          call VOLMWpf(R,VR,u,DAGGER,-MDW,SRF)
          TMP=VR
          call mGmu(TMP,4)
          DR=VR+zi*mass*TMP
          DR=(DL+DR)/two
        endif

      endif

      return
      end subroutine DOLMW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLMW(R,DR,u,DAGGER,mass,SRF)
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
      if(VB_OL)then ; print *,"IDOLMW" ; endif
      iopts%SRF=SRF
      iopts%mass=mass
      Mptr => DOLMW
      call IMnonsymOpts(R,DR,u,DAGGER,Mptr,iopts)

      return
      end subroutine IDOLMW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLMS(R,DR,u,DAGGER,mass,SRF)
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
      complex(prc) :: TMP(Nv,4),VR(Nv,4),DL(Nv,4),IMG3R(Nv,4)

      if (MTYPE.eq.1) then

        call VOLMSpf(R,TMP,u,DAGGER,-MDW,SRF)
        DR=(half+mass/two)*R+(half-mass/two)*TMP

      elseif (MTYPE.eq.3) then ! DOL = (1+i.m.g3) + (1-i.m.g3).V

        if (.not.DAGGER) then
          call VOLMSpf(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          DR=zi*mass*(R-VR)/two
          call mGmu(DR,4)
          DR=DL+DR
        elseif (DAGGER) then
          call VOLMSpf(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          TMP=R
          call mGmu(TMP,4)
          call VOLMSpf(TMP,VR,u,DAGGER,-MDW,SRF)
          DR=-zi*mass*(TMP-VR)/two
          DR=DL+DR
        endif

      elseif (MTYPE.eq.4) then 

        if (.not.DAGGER) then ! DOL=(1+im.g3) + V(1-im.g3)
          TMP=R
          call mGmu(TMP,4)
          DL=R+zi*mass*TMP
          VR=R-zi*mass*TMP
          call VOLMSpf(VR,DR,u,DAGGER,-MDW,SRF)
          DR=(DL+DR)/two
        elseif (DAGGER) then ! DOLdag=(1-im.g3) + (1+im.g3)Vdag
          TMP=R
          call mGmu(TMP,4)
          DL=R-zi*mass*TMP
          call VOLMSpf(R,VR,u,DAGGER,-MDW,SRF)
          TMP=VR
          call mGmu(TMP,4)
          DR=VR+zi*mass*TMP
          DR=(DL+DR)/two
        endif

      else

        print *,"MTYPE",MTYPE,"not supported"

      endif
      return
      end subroutine DOLMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLMS(R,DR,u,DAGGER,mass,SRF)
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
      Mptr => DOLMS
      call IMnonsymOpts(R,DR,u,DAGGER,Mptr,iopts)

      return
      end subroutine IDOLMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module overlapmoduleprod
