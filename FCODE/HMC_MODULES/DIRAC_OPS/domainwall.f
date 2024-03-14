!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      use WilsonDirac
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_Wilson(R,DR,u,DAGGER,mass)
      use gammas
!      use ratfuncs
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
!      use gammas
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW(R,DR,u,DAGGER,mass)
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass

      if (DWkernel.eq.2) then
        call DDW_Wilson(R,DR,u,DAGGER,mass)
      elseif (DWkernel.eq.3) then
        call DDW_OWilson(R,DR,u,DAGGER,mass)
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

      if (DWkernel.eq.2) then
        Dptr => DDW_Wilson
      elseif (DWkernel.eq.3) then
        Dptr => DDW_OWilson
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

      call DDW_Wilson(R,TMP1,u,.false.,baremass)
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

      if (DWkernel.eq.2) then
        Dptr=>DDW_Wilson
        print *,"KDDW with Wilson"
      elseif (DWkernel.eq.2) then
        Dptr=>DDW_OWilson
        print *,"KDDW with OWilson"
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

      if (DWkernel.eq.2) then
        Dptr=>DDW_Wilson
        print *,"KDDW with Wilson"
      elseif (DWkernel.eq.3) then
        Dptr=>DDW_OWilson
        print *,"KDDW with OWilson"
      else
        print *,"DW kernel not set properly"
        stop
      endif
      MTMP=MTYPE
      gi=4
      if (.not.DAGGER) then
        call PermM(R,TMP,DAGGER,gi)
        MTYPE=1
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
      end module domainwallmod
