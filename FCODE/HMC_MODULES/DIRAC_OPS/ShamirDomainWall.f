!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module ShamirDomWall
      use arraysizes
      use options
      use WilsonDirac
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_Shamir(R,DR,u,DAGGER,mass)
      use gammas
      use ratfuncs
      implicit none
c     calculates DR = DDW*R where DDW is the Shamir domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TR(Nv,4)
      integer s
      real(prc) as

c     diagonal blocks
      do s=1,Ls
        call DWilson(R(:,:,s),DR(:,:,s),u,DAGGER,-MDW)
      end do
      DR=DR+R

c     projection blocks, P_-.R  = R(:,3:4,s), P_+.R = R(:,1:2,s)
      if (.not.DAGGER) then
        do s=1,Ls-1
          DR(:,3:4,s)=DR(:,3:4,s)-R(:,3:4,s+1) 
          DR(:,1:2,s+1)=DR(:,1:2,s+1)-R(:,1:2,s)
        end do
      else
        do s=1,Ls-1
          DR(:,1:2,s)=DR(:,1:2,s)-R(:,1:2,s+1) 
          DR(:,3:4,s+1)=DR(:,3:4,s+1)-R(:,3:4,s)
        end do
      end if

c     mass terms
 
      if (MTYPE.eq.1) then
        if (.not.DAGGER) then
          DR(:,1:2,1)=DR(:,1:2,1)+mass*R(:,1:2,Ls)
          DR(:,3:4,Ls)=DR(:,3:4,Ls)+mass*R(:,3:4,1)
        else
          DR(:,3:4,1)=DR(:,3:4,1)+mass*R(:,3:4,Ls)
          DR(:,1:2,Ls)=DR(:,1:2,Ls)+mass*R(:,1:2,1)
        end if
      elseif (MTYPE.eq.4) then
        if (.not.DAGGER) then
          TR=zi*R(:,:,Ls)
          call mGmu(TR,4)
          DR(:,1:2,1)=DR(:,1:2,1)+mass*TR(:,1:2)
          TR=zi*R(:,:,1)
          call mGmu(TR,4)
          DR(:,3:4,Ls)=DR(:,3:4,Ls)+mass*TR(:,3:4)
        else
          TR=zi*R(:,:,Ls)
          call mGmu(TR,4)
          DR(:,3:4,1)=DR(:,3:4,1)-mass*TR(:,3:4)
          TR=zi*R(:,:,1)
          call mGmu(TR,4)
          DR(:,1:2,Ls)=DR(:,1:2,Ls)-mass*TR(:,1:2)
        endif
      elseif (MTYPE.eq.2) then
        if (.not.DAGGER) then
          DR(:,1:2,1)=DR(:,1:2,1)+zi*mass*R(:,1:2,Ls)
          DR(:,3:4,Ls)=DR(:,3:4,Ls)-zi*mass*R(:,3:4,1)
        else
          DR(:,3:4,1)=DR(:,3:4,1)+zi*mass*R(:,3:4,Ls)
          DR(:,1:2,Ls)=DR(:,1:2,Ls)-zi*mass*R(:,1:2,1)
        endif
      elseif (MTYPE.eq.3) then
        if (.not.DAGGER) then
          DR(:,1:2,1)=DR(:,1:2,1)-zi*mass*R(:,1:2,Ls)
          DR(:,3:4,Ls)=DR(:,3:4,Ls)+zi*mass*R(:,3:4,1)
        else
          DR(:,3:4,1)=DR(:,3:4,1)-zi*mass*R(:,3:4,Ls)
          DR(:,1:2,Ls)=DR(:,1:2,Ls)+zi*mass*R(:,1:2,1)
        endif
      endif

      return
      end subroutine DDW_Shamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ShamirDomainWallDerivs(dSdA,eta,nu,DAG)
      use numbers
      use gammas
      use indices
      use WilsonDirac
      implicit none
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls),intent(in) ::  eta,nu
      logical,intent(in) :: DAG
      real(prc) :: tmp(Nv,3)
      complex(prc),dimension(Nv,4) ::  eta_l,nu_l
      integer l
      integer pm1

      dSdA=0
      do l=1,Ls  ! diagonal terms
        eta_l=eta(:,:,l)
        nu_l=nu(:,:,l)
        call WilsonDerivs(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+tmp
      end do

      return
      end subroutine ShamirDomainWallDerivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module ShamirDomWall
