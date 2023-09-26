      module dwcoeffs
      use pacc
      use options
      implicit none
      real(zprc) omega(Ls)
      end module dwcoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module domainwallmod
      use arraysizes
      use options
      use basicdiracopsmod
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_Shamir(R,DR,u,DAGGER,mass)
      use gammas
      use ratfuncs
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TR(Nv,4)
      integer s
      real(prc) as

c     diagonal blocks
c      do s=1,Ls
c        call DWilson(R(:,:,s),DR(:,:,s),u,DAGGER,-MDW,cone)
c      end do
      as=one
      do s=1,Ls
        call DWilson(R(:,:,s),DR(:,:,s),u,DAGGER,-MDW,czero)
      end do
      DR=as*DR+R

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_OShamir(R,DR,u,DAGGER,mass)
      use gammas
      use ratfuncs
      use dwcoeffs
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TR(Nv,4)
      integer s
      real(prc) as

c     diagonal blocks
      as=one
      do s=1,Ls
        call DWilson(R(:,:,s),DR(:,:,s),u,DAGGER,-MDW,czero)
        DR(:,:,s)=omega(s)*DR(:,:,s)+R(:,:,s)
      end do

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
      end subroutine DDW_OShamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_Wilson(R,DR,u,DAGGER,mass)
      use gammas
      use ratfuncs
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
        call DWilson(R(:,:,s),DR(:,:,s),u,DAGGER,-MDW,cone)
      end do

      gi=4 ! use gamma3
      if (.not. DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call Pminus(R(:,:,s+1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,-cone)
          DR(:,:,s)=DR(:,:,s)+TR2
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call Pplus(R(:,:,s),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,-cone)
          DR(:,:,s+1)=DR(:,:,s+1)+TR2
        end do
      elseif (DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s+1),TR,u,DAGGER,-MDW,-cone)
          call Pplus(TR,TR2,gi)
          DR(:,:,s)=DR(:,:,s)+TR2
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s),TR,u,DAGGER,-MDW,-cone)
          call Pminus(TR,TR2,gi)
          DR(:,:,s+1)=DR(:,:,s+1)+TR2
        end do
      end if

c     mass terms
      if (MTYPE.eq.1) then

        if (.not.DAGGER) then
          call Pminus(R(:,:,1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,-cone)
        elseif (DAGGER) then
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW,-cone)
          call Pplus(TR,TR2,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)-mass*TR2

        if (.not.DAGGER) then
          call Pplus(R(:,:,Ls),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,-cone)
        elseif (DAGGER) then
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW,-cone)
          call Pminus(TR,TR2,gi)
        endif
        DR(:,:,1)=DR(:,:,1)-mass*TR2

      elseif ((MTYPE.eq.2).or.(MTYPE.eq.3)) then

        if (.not.DAGGER) then
          zkappa=cmplx(0,-mass)
          call Pminus(R(:,:,1),TR,gi)
          call mGmu(TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,-cone)
        elseif (DAGGER) then
          zkappa=cmplx(0,mass)
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW,-cone)
          call mGmu(TR,gi)
          call Pplus(TR,TR2,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)+zkappa*TR2

        if (.not.DAGGER) then
          zkappa=cmplx(0,-mass)
          call Pplus(R(:,:,Ls),TR,gi)
          call mGmu(TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,-cone)
        elseif(DAGGER) then
          zkappa=cmplx(0,mass)
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW,-cone)
          call mGmu(TR,gi)
          call Pminus(TR,TR2,gi)
        endif
        DR(:,:,1)=DR(:,:,1)+zkappa*TR2

      endif

      return
      end subroutine DDW_Wilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_OWilson(R,DR,u,DAGGER,mass)
      use gammas
      use ratfuncs
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

c      omega=1e0

c     diagonal blocks
      do s=1,Ls
        call DWilson(R(:,:,s),TR,u,DAGGER,-MDW,czero)
        DR(:,:,s)=omega(s)*TR+R(:,:,s)
      end do

      gi=4
      if (.not. DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call Pminus(R(:,:,s+1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,czero)
          DR(:,:,s)=DR(:,:,s)+omega(s)*TR2-TR
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call Pplus(R(:,:,s),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,czero)
          DR(:,:,s+1)=DR(:,:,s+1)+omega(s+1)*TR2-TR
        end do
      elseif (DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s+1),TR,u,DAGGER,-MDW,czero)
          TR2=omega(s+1)*TR-R(:,:,s+1)
c          TR2=omega(s)*TR-R(:,:,s+1)
          call Pplus(TR2,TR3,gi)
          DR(:,:,s)=DR(:,:,s)+TR3
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s),TR,u,DAGGER,-MDW,czero)
          TR2=omega(s)*TR-R(:,:,s)
c          TR2=omega(s+1)*TR-R(:,:,s)
          call Pminus(TR2,TR3,gi)
          DR(:,:,s+1)=DR(:,:,s+1)+TR3
        end do
      end if

c     mass terms
      if (MTYPE.eq.1) then

        if (.not.DAGGER) then
          call Pminus(R(:,:,1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,czero)
          TR3=omega(Ls)*TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW,czero)
          TR2=omega(1)*TR-R(:,:,1)
          call Pplus(TR2,TR3,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)-mass*TR3

        if (.not.DAGGER) then
          call Pplus(R(:,:,Ls),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,czero)
          TR3=omega(1)*TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW,czero)
          TR2=omega(Ls)*TR-R(:,:,Ls)
          call Pminus(TR2,TR3,gi)
        endif
        DR(:,:,1)=DR(:,:,1)-mass*TR3

      elseif ((MTYPE.eq.2).or.(MTYPE.eq.3)) then

        if (.not.DAGGER) then
          zkappa=cmplx(0,-mass)
          call Pminus(R(:,:,1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,czero)
          TR3=omega(Ls)*TR2-TR
        elseif (DAGGER) then
          zkappa=cmplx(0,-mass)
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW,czero)
          TR2=omega(1)*TR-R(:,:,1)
          call Pplus(TR2,TR3,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)+zkappa*TR3

        if (.not.DAGGER) then
          zkappa=cmplx(0,mass)
          call Pplus(R(:,:,Ls),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW,czero)
          TR3=omega(1)*TR2-TR
        elseif(DAGGER) then
          zkappa=cmplx(0,mass)
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW,czero)
          TR2=omega(Ls)*TR-R(:,:,Ls)
          call Pminus(TR2,TR3,gi)
        endif
        DR(:,:,1)=DR(:,:,1)+zkappa*TR3

      endif

      return
      end subroutine DDW_OWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW(R,DR,u,DAGGER,mass)
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass

      if (DWkernel.eq.1) then
        call DDW_Shamir(R,DR,u,DAGGER,mass)
      elseif (DWkernel.eq.2) then
        call DDW_Wilson(R,DR,u,DAGGER,mass)
      elseif (DWkernel.eq.3) then
        call DDW_OWilson(R,DR,u,DAGGER,mass)
      elseif (DWkernel.eq.4) then
        call DDW_OShamir(R,DR,u,DAGGER,mass)
      else
        print *,"DWkernel not set properly"
        stop
      endif

      return
      end subroutine DDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDDW(RR,DR,u,DAGGER,mass)
!     solve DDW.DR = RR
      use options
      implicit none
      complex(prc) RR(Nv,4,Ls),DR(Nv,4,Ls)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) TMP(Nv,4,Ls)
      procedure(),pointer :: Dptr=>NULL()

      if (DWkernel.eq.1) then
        Dptr => DDW_Shamir
      elseif (DWkernel.eq.2) then
        Dptr => DDW_Wilson
      elseif (DWkernel.eq.3) then
        Dptr => DDW_OWilson
      elseif (DWkernel.eq.4) then
        Dptr => DDW_OShamir
      else
        print *, "Domain Wall kernel not set properly"
        stop
      endif

      if (.not. DAGGER) then
        call Dptr(RR,TMP,u,.true.,mass)
        call IMdagM_DWkernel(TMP,DR,u,.false.,mass,Dptr) ! rem dagger is redundant here
      elseif (DAGGER) then
        call IMdagM_DWkernel(RR,TMP,u,.true.,mass,Dptr) ! rem dagger is redundant here
        call Dptr(TMP,DR,u,.false.,mass)
      end if

      return
      end subroutine IDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine IMdagM_DW(RR,DR,u,DAGGER,mass)
!     solve MdagM.DR = RR for M=DDW
c      implicit none
c      integer,parameter :: kferm = 4*Nv*Ls
c      integer,parameter :: niterc=Nv*Ls
c      complex(prc) RR(kferm), DR(kferm)
c      complex(prc) u(Nv,3)
c      logical DAGGER
c      real(prc) mass

c      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
c      integer itercg
c      integer nx,i
c      real(prc) beta,betan,betad,alpha,alphad,alphan

c      itercg=0

c     initialise
c      DR=RR
c      call DDW_Shamir(RR,x1,u,.false.,mass)
c      call DDW_Shamir(x1,x2,u,.true.,mass)
c      r=RR-x2
c      p=r
c      betan=sum(conjg(r)*r)
c      if (betan.lt.resid) goto 8
c      
c      do i=1,niterc
c        itercg=itercg+1
c        call DDW_Shamir(p,x1,u,.false.,mass)
c        call DDW_Shamir(x1,x2,u,.true.,mass)
c        alphan=sum(conjg(r)*r)
c        alphad=sum(conjg(p)*x2)
c        alpha=alphan/alphad
c
c        DR=DR+alpha*p
c        r=r-alpha*x2
c        betan=sum(conjg(r)*r)
c        if (betan.lt.resid) goto 8
c        beta=betan/alphan
c        p=r+beta*p
c      end do
c
c8     continue
c      print *,itercg,niterc,betan
c      return
c      end subroutine IMdagM_DW    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagM_DWkernel(RR,DR,u,DAGGER,mass,Dptr)
!     solve MdagM.DR = RR for M=DDW
      use countmod
      implicit none
      integer,parameter :: kferm = 4*Nv*Ls
      integer,parameter :: niterc=Nv*Ls
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      procedure(),pointer :: Dptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Dptr(RR,x1,u,.false.,mass)
      call Dptr(x1,x2,u,.true.,mass)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Dptr(p,x1,u,.false.,mass)
        call Dptr(x1,x2,u,.true.,mass)
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
!      print *,itercg,niterc,betan
      oc_idx=oc_idx+1
      outer_count=outer_count+itercg

      return
      end subroutine IMdagM_DWkernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_calcPhi(R,DR,u,DAGGER)
      use options
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      complex(prc) :: TMP1(Nv,4,Ls),TMP2(Nv,4,Ls)

      call DDW_Shamir(R,TMP1,u,.false.,baremass)
!      call DDW(R,TMP1,u,.false.,baremass)
      DR(:,1:2)=TMP1(:,1:2,Ls)
      DR(:,3:4)=TMP1(:,3:4,Ls)

      end subroutine DDW_calcPhi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDDW_calcPhi(R,DR,u,DAGGER)
      use options
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      complex(prc) :: TMP1(Nv,4,Ls),TMP2(Nv,4,Ls)

      call IDDW(R,TMP1,u,.false.,baremass)
      DR(:,1:2)=TMP1(:,1:2,Ls)
      DR(:,3:4)=TMP1(:,3:4,Ls)

      end subroutine IDDW_calcPhi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine PermM(R,DR,DAGGER,gi)
      use options
      implicit none
c     calculates DR(:,:,l) = Pminus.R(:,:,l)+Pplus.R(:,:,l+1)
c     if DAGGER=.true. then
c     calculates DR(:,:,l) = Pminus.R(:,:,l)+Pplus.R(:,:,l-1)
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      logical DAGGER
      integer gi
      complex(prc) :: TMP(Nv,4)
      integer l

      do l=1,Ls
        call Pminus(R(:,:,l),DR(:,:,l),gi)
      enddo
      do l=1,Ls
        call Pplus(R(:,:,l),TMP,gi)
        if (.not.DAGGER) then
          if (l.eq.1) then
            DR(:,:,Ls)=DR(:,:,Ls)+TMP
          else
            DR(:,:,l-1)=DR(:,:,l-1)+TMP
          endif
        else
          if (l.eq.Ls) then
            DR(:,:,1)=DR(:,:,1)+TMP
          else
            DR(:,:,l+1)=DR(:,:,l+1)+TMP
          endif
        endif
      enddo

      return
      end subroutine PermM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine KDDW(R,DR,u,DAGGER,mass)
      implicit none
c     calculates DR = Pdag.IDDW(1).DDW(m).P.R where P is the permutation matrix
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TMP(Nv,4,Ls)
      integer gi
      procedure(),pointer :: Dptr => NULL()
      integer MTMP

      if (DWkernel.eq.1) then
        Dptr=>DDW_Shamir
        print *,"KDDW with Shamir"
      elseif (DWkernel.eq.2) then
        Dptr=>DDW_Wilson
        print *,"KDDW with Wilson"
      elseif (DWkernel.eq.3) then
        Dptr=>DDW_OWilson
        print *,"KDDW with OWilson"
      elseif (DWkernel.eq.4) then
        Dptr=>DDW_OShamir
        print *,"KDDW with OShamir"
      else
        print *,"DW kernel not set properly"
        stop
      endif
      print *,"MTYPE:",MTYPE

      MTMP=MTYPE
      gi=4
      if (.not.DAGGER) then
        call PermM(R,TMP,.false.,gi)
        call Dptr(TMP,DR,u,.false.,mass)
        MTYPE=1
        call IDDW(DR,TMP,u,.false.,one)
        call PermM(TMP,DR,.true.,gi)
      elseif(DAGGER) then
        call PermM(R,TMP,.false.,gi)
        MTYPE=1
        call IDDW(TMP,DR,u,.true.,one)
        MTYPE=MTMP
        call Dptr(DR,TMP,u,.true.,mass)
        call PermM(TMP,DR,.true.,gi)
      endif
      MTYPE=MTMP

      return
      end subroutine KDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IKDDW(R,DR,u,DAGGER,mass)
      implicit none
c     calculates DR = Pdag.IDDW(m).DDW(1).P.R where P is the permutation matrix
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TMP(Nv,4,Ls)
      procedure(),pointer :: Dptr => NULL()
      integer gi
      integer MTMP

      if (DWkernel.eq.1) then
        Dptr=>DDW_Shamir
        print *,"KDDW with Shamir"
      elseif (DWkernel.eq.2) then
        Dptr=>DDW_Wilson
        print *,"KDDW with Wilson"
      elseif (DWkernel.eq.3) then
        Dptr=>DDW_OWilson
        print *,"KDDW with OWilson"
      elseif (DWkernel.eq.4) then
        Dptr=>DDW_OShamir
        print *,"KDDW with OShamir"
      else
        print *,"DW kernel not set properly"
        stop
      endif
      MTMP=MTYPE
      gi=4
      if (.not.DAGGER) then
        call PermM(R,TMP,DAGGER,gi)
        MTYPE=1
c        if (DWkernel.eq.1) then
c          call DDW_Shamir(TMP,DR,u,DAGGER,one)
c        elseif (DWkernel.eq.2) then
c          call DDW_Wilson(TMP,DR,u,DAGGER,one)
c        endif
        call Dptr(TMP,DR,u,DAGGER,one)
        MTYPE=MTMP
        call IDDW(DR,TMP,u,DAGGER,mass)
        call PermM(TMP,DR,.not.DAGGER,gi)
      elseif (DAGGER) then
        call PermM(R,TMP,.not.DAGGER,gi)
        call IDDW(TMP,DR,u,DAGGER,mass)
        MTYPE=1
        call Dptr(DR,TMP,u,DAGGER,one)
        call PermM(TMP,DR,DAGGER,gi)
        MTYPE=MTMP
      endif

      return
      end subroutine IKDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine KDDW4(R,DR,u,DAGGER,mass)
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: R5(Nv,4,Ls)
      complex(prc) :: DR5(Nv,4,Ls)
      complex(prc) :: TMP(Nv,4,Ls)
      integer gi

      R5=czero
      R5(:,:,1)=R
      call KDDW(R5,DR5,u,DAGGER,mass)
      DR=DR5(:,:,1)
      
      return
      end subroutine KDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IKDDW4(R,DR,u,DAGGER,mass)
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: R5(Nv,4,Ls)
      complex(prc) :: DR5(Nv,4,Ls)
      complex(prc) :: TMP(Nv,4,Ls)
      integer gi

      R5=czero
      R5(:,:,1)=R
      call IKDDW(R5,DR5,u,DAGGER,mass)
      DR=DR5(:,:,1)
      
      return
      end subroutine IKDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine calcKDWDeterminant(u,DAGGER,mass)
c      use options
c      implicit none
c      complex(prc),intent(in) :: u(Nv,3)
c      logical DAGGER
c      real(prc) mass
c      complex(prc) :: R(4*Nv)
c      complex(prc) :: COL(4*Nv)
c      integer i,Nv4
c      integer ipiv(4*Nv),err
c      complex(prc) :: DenseOL(4*Nv,4*Nv)
c      complex(prc) DdD(4*Nv,4*Nv),det
c
c      DenseOL=czero
c      Nv4=4*Nv
c      do i=1,Nv4
c        R=czero
c        R(i)=cone
c        call KDDW4(R,COL,u,DAGGER,mass)
c        DenseOL(1:Nv4,i)=COL(1:Nv4)
c      end do
c      print *,DenseOL
c      DdD=matmul(DenseOL,conjg(transpose(DenseOL)))
c      print *,"calc K determinant"
c      call  cgetrf(Nv4,Nv4,DdD,Nv4,ipiv,err)
c      print *,"err",err
c
c      det=1
c      do i=1,Nv4
c        det=det*DdD(i,i)
c        print *,DdD(i,i)
c      end do
c      print *,"det K",det
c
c      return
c      end subroutine calcKDWDeterminant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module domainwallmod
