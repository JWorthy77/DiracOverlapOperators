      module weylmod
      use pacc
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DWeyl(R,DR,u,DAGGER,mass,add)
      use paulimodule
      use indices
      implicit none
c     calculates DR = D*R
      complex(prc),intent(in) :: R(Nv,2)
      complex(prc),intent(out) :: DR(Nv,2)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      real(prc) mult
      complex(prc) cmass
      integer i,idirac
      integer mu,igork

      mult=one
      if (DAGGER) then
        mult=-one
      end if

      cmass=mass+add
      if (DAGGER) then
        cmass=conjg(mass+add)
      end if

!     mass term 
      DR=cmass*R

!     Dirac term (anti-Hermitian)
      do mu=1,3
        do idirac=1,2
          igork=pauliin(mu,idirac)
          do i=1,Nv
            DR(i,idirac)=DR(i,idirac)
     &        +mult*paulival(mu,idirac)*
     &         (u(i,mu)*R(iu(i,mu),igork)
     &         -conjg(u(id(i,mu),mu))*R(id(i,mu),igork))/two
          enddo
        enddo
      enddo

!     Wilson term (Hermitian)
      do mu=1,3
        do idirac=1,2
          do i=1,Nv
            DR(i,idirac)=DR(i,idirac)
     &        -( u(i,mu)*R(iu(i,mu),idirac) - two*R(i,idirac)
     &         +conjg(u(id(i,mu),mu))*R(id(i,mu),idirac))/two
          enddo
        enddo
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDWeyl(RR,DR,u,DAGGER,mass,add)
!     solve Dw.DR = RR or Dw(dag).DR = RR
      implicit none
      complex(prc),intent(in) :: RR(Nv,2)
      complex(prc),intent(out) :: DR(Nv,2)
      complex(prc),intent(in) ::  u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,2)

      if (.not. DAGGER) then
        call DWeyl(RR,TMP,u,.true.,mass,add)
        call IDWdagDW(TMP,DR,u,.false.,mass,add)
      elseif (DAGGER) then
        call DWeyl(RR,TMP,u,.false.,mass,add)
        call IDWdagDW(TMP,DR,u,.true.,mass,add)
      end if

      return
      end subroutine IDWeyl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDWdagDW(RR,DR,u,DAGGER,mass,add)
      implicit none
      integer,parameter :: kferm = 2*Nv
      integer,parameter :: niterc=Nv*10
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
      call DWeyl(RR,x1,u,DAGGER,mass,add)
      call DWeyl(x1,x2,u,.not.DAGGER,mass,add)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call DWeyl(p,x1,u,DAGGER,mass,add)
        call DWeyl(x1,x2,u,.not.DAGGER,mass,add)
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
      end subroutine IDWdagDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WdagW(RR,DR,u,DAGGER,mass,add)
      implicit none
      complex(prc),intent(in) :: RR(Nv,2)
      complex(prc),intent(out) :: DR(Nv,2)
      complex(prc),intent(in) ::  u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,2)

      call DWeyl(RR,TMP,u,DAGGER,mass,add)
      call DWeyl(TMP,DR,u,.not.DAGGER,mass,add)

      return
      end subroutine WdagW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WdagWpC(RR,DR,u,DAGGER,mass,add) 
      implicit none
      complex(prc),intent(in) :: RR(Nv,2)
      complex(prc),intent(out) :: DR(Nv,2)
      complex(prc),intent(in) ::  u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,2)

      call DWeyl(RR,TMP,u,DAGGER,mass,czero)
      call DWeyl(TMP,DR,u,.NOT.DAGGER,mass,czero)
      DR=DR+add*RR

      return
      end subroutine WdagWpC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testweyl()
      use rvmodule
      use gaugefield
      implicit none
      complex(prc) R(Nv,2),DR(Nv,2),TMP(Nv,2),TMP2(Nv,2)
    
      call setRVs(2*Nv,R)

      call DWeyl(R,TMP2,u,.TRUE.,zero,czero)
      call DWeyl(TMP2,TMP,u,.FALSE.,zero,czero)
      call WdagW(R,TMP2,u,.TRUE.,zero,czero)
      print *,"WdagW-Wd.W:",maxval(abs(TMP2-TMP))
      call IDWdagDW(TMP,DR,u,.true.,zero,czero)
      print *,"WdagW.IWdagW (true):",maxval(abs(DR-R))

      call DWeyl(R,TMP,u,.false.,zero,czero)
      call IDWeyl(TMP,DR,u,.false.,zero,czero)
      print *,"W.IW (false):",maxval(abs(DR-R))

      call DWeyl(R,TMP,u,.true.,zero,czero)
      call IDWeyl(TMP,DR,u,.true.,zero,czero)
      print *,"W.IW (true):",maxval(abs(DR-R))

      call DWeyl(R,TMP,u,.true.,zero,czero)
      call DWeyl(TMP,TMP2,u,.false.,zero,czero)
      call IDWeyl(TMP2,TMP,u,.false.,zero,czero)
      call IDWeyl(TMP,DR,u,.true.,zero,czero)
      print *,"Wdag.W.IW.IWdag:",maxval(abs(DR-R))

      call DWeyl(R,TMP,u,.true.,zero,czero)
      call DWeyl(TMP,TMP2,u,.false.,zero,czero)
      call WdagWpC(R,TMP,u,.true.,zero,czero)
      print *,"WdagWpC:",maxval(abs(TMP2-TMP))
      return
      end subroutine testweyl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module weylmod

