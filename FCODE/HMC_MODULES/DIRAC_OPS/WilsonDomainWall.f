!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module dwcoeffs
      use pacc
      use options
      implicit none
      real(zprc) omega(Ls)
      end module dwcoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module WilsonDomWall
      use arraysizes
      use options
      use WilsonDirac
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_Wilson(R,DR,u,DAGGER,mass)
      use gammas
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
c     with Wilson kernel
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TR(Nv,4),TR2(Nv,4)
      integer s,gi
      complex(prc) zkappa

c     diagonal blocks
      do s=1,Ls
        call DWilson(R(:,:,s),DR(:,:,s),u,DAGGER,-MDW)
        DR(:,:,s)=DR(:,:,s)+R(:,:,s)
      end do

      gi=4 ! use gamma3
      if (.not. DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call Pminus(R(:,:,s+1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
          DR(:,:,s)=DR(:,:,s)+TR2
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call Pplus(R(:,:,s),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
          DR(:,:,s+1)=DR(:,:,s+1)+TR2
        end do
      elseif (DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s+1),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,s+1)
          call Pplus(TR,TR2,gi)
          DR(:,:,s)=DR(:,:,s)+TR2
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,s)
          call Pminus(TR,TR2,gi)
          DR(:,:,s+1)=DR(:,:,s+1)+TR2
        end do
      end if

c     mass terms
      if (MTYPE.eq.1) then

        if (.not.DAGGER) then
          call Pminus(R(:,:,1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,1)
          call Pplus(TR,TR2,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)-mass*TR2

        if (.not.DAGGER) then
          call Pplus(R(:,:,Ls),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,Ls)
          call Pminus(TR,TR2,gi)
        endif
        DR(:,:,1)=DR(:,:,1)-mass*TR2

      elseif ((MTYPE.eq.2).or.(MTYPE.eq.3)) then

        if (.not.DAGGER) then
          zkappa=cmplx(0,-mass)
          call Pminus(R(:,:,1),TR,gi)
          call mGmu(TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
        elseif (DAGGER) then
          zkappa=cmplx(0,mass)
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,1)
          call mGmu(TR,gi)
          call Pplus(TR,TR2,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)+zkappa*TR2

        if (.not.DAGGER) then
          zkappa=cmplx(0,-mass)
          call Pplus(R(:,:,Ls),TR,gi)
          call mGmu(TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
        elseif(DAGGER) then
          zkappa=cmplx(0,mass)
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,Ls)
          call mGmu(TR,gi)
          call Pminus(TR,TR2,gi)
        endif
        DR(:,:,1)=DR(:,:,1)+zkappa*TR2

      endif

      return
      end subroutine DDW_Wilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_OWilson(R,DR,u,DAGGER,mass)
      use dwcoeffs
      implicit none
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TR(Nv,4),TR2(Nv,4),TR3(Nv,4)
      integer s,gi
      complex(prc) zkappa

c     diagonal blocks
      do s=1,Ls
        call DWilson(R(:,:,s),TR,u,DAGGER,-MDW)
        DR(:,:,s)=omega(s)*TR+R(:,:,s)
      end do

      gi=4
      if (.not. DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call Pminus(R(:,:,s+1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          DR(:,:,s)=DR(:,:,s)+omega(s)*TR2-TR
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call Pplus(R(:,:,s),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          DR(:,:,s+1)=DR(:,:,s+1)+omega(s+1)*TR2-TR
        end do
      elseif (DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s+1),TR,u,DAGGER,-MDW)
          TR2=omega(s+1)*TR-R(:,:,s+1)
          call Pplus(TR2,TR3,gi)
          DR(:,:,s)=DR(:,:,s)+TR3
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s),TR,u,DAGGER,-MDW)
          TR2=omega(s)*TR-R(:,:,s)
          call Pminus(TR2,TR3,gi)
          DR(:,:,s+1)=DR(:,:,s+1)+TR3
        end do
      end if

c     mass terms
      if (MTYPE.eq.1) then

        if (.not.DAGGER) then
          call Pminus(R(:,:,1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR3=omega(Ls)*TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW)
          TR2=omega(1)*TR-R(:,:,1)
          call Pplus(TR2,TR3,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)-mass*TR3

        if (.not.DAGGER) then
          call Pplus(R(:,:,Ls),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR3=omega(1)*TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW)
          TR2=omega(Ls)*TR-R(:,:,Ls)
          call Pminus(TR2,TR3,gi)
        endif
        DR(:,:,1)=DR(:,:,1)-mass*TR3

      elseif ((MTYPE.eq.2).or.(MTYPE.eq.3)) then

        if (.not.DAGGER) then
          zkappa=cmplx(0,-mass)
          call Pminus(R(:,:,1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR3=omega(Ls)*TR2-TR
        elseif (DAGGER) then
          zkappa=cmplx(0,-mass)
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW)
          TR2=omega(1)*TR-R(:,:,1)
          call Pplus(TR2,TR3,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)+zkappa*TR3

        if (.not.DAGGER) then
          zkappa=cmplx(0,mass)
          call Pplus(R(:,:,Ls),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR3=omega(1)*TR2-TR
        elseif(DAGGER) then
          zkappa=cmplx(0,mass)
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW)
          TR2=omega(Ls)*TR-R(:,:,Ls)
          call Pminus(TR2,TR3,gi)
        endif
        DR(:,:,1)=DR(:,:,1)+zkappa*TR3

      endif

      return
      end subroutine DDW_OWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pplus(R,DR,gi)
      use gammas
      implicit none
c     calculates DR = (1+gamma)/2 R (gi should be 4 typically)
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      integer gi
      integer v,d,di

      do v=1,Nv
        do d=1,4
          di=gamin(gi,d)
          DR(v,d)=(R(v,d)+gamval(gi,d)*R(v,di))/two
        enddo
      enddo
      
      return
      end subroutine Pplus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pminus(R,DR,gi)
      use gammas
      implicit none
c     calculates DR = (1-gamma)/2 R (gi should be 4 typically)
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      integer gi
      integer v,d,di

      do v=1,Nv
        do d=1,4
          di=gamin(gi,d)
          DR(v,d)=(R(v,d)-gamval(gi,d)*R(v,di))/two
        enddo
      enddo
      
      return
      end subroutine Pminus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDomainWallDerivs(dSdA,eta,nu,DAG)
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

      do l=1,Ls-1  ! upper diagonal
        if (.not.DAG) then
          eta_l=eta(:,:,l)
        else
          call Pplus(eta(:,:,l),eta_l,4)
        endif
        if (.not.DAG) then
          call Pminus(nu(:,:,l+1),nu_l,4)
        else
          nu_l=nu(:,:,l+1)
        endif
        call WilsonDerivs(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+tmp
      end do

      do l=2,Ls  ! lower diagonal
        if (.not.DAG) then
          eta_l=eta(:,:,l)
        else
          call Pminus(eta(:,:,l),eta_l,4)
        endif
        if (.not.DAG) then
          call Pplus(nu(:,:,l-1),nu_l,4)
        else
          nu_l=nu(:,:,l-1)
        endif
        call WilsonDerivs(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+tmp
      end do

      return
      end subroutine WilsonDomainWallDerivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module WilsonDomWall
