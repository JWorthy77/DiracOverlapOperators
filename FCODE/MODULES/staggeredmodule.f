      module staggeredmodule
#ifdef STAGGERED
      use pacc
      use arraysizes
      use numbers
      use options
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DStagGN(R,DR,sigma,DAGGER) ! from Hands/Kocic/Kogut
      use indices
      implicit none
c     calculates DR = DStagGN*R
      real(prc),intent(in) :: R(Nv)
      real(prc),intent(out) :: DR(Nv)
      real(prc),intent(in) :: sigma(Nv)
      logical DAGGER
      real(prc) mult
      integer i,i1,i2,i3,nn(3),idx
      integer mu
      integer id111,id011,id101,id110,id001,id010,id100,id000
      real(prc) :: s111,s011,s101,s110,s001,s010,s100,s000
      real(prc) :: stot

      mult=one
      if (DAGGER) then
        mult=-one
      end if

      DR = zero
      nn(1)=1
      do i=1,Nv
        idx=i
        call getLatticeCoords(idx,i1,i2,i3)
c        print *,i1,i2,i3
        nn(2)=(-1)**i1
        nn(3)=(-1)**(i1+i2)
        do mu=1,NDT
          DR(i)=DR(i)+half*mult*nn(mu)*(R(iu(i,mu))-R(id(i,mu)))
        enddo
        id111=i
        id011=id(i,1)
        id101=id(i,2)
        id110=id(i,3)
        id001=id(id011,2)
        id010=id(id011,3)
        id100=id(id101,3)
        id000=id(id100,1)

        s111=sigma(id111)
        s011=sigma(id011)
        s101=sigma(id101)
        s110=sigma(id110)
        s001=sigma(id001)
        s010=sigma(id010)
        s100=sigma(id100)
        s000=sigma(id000)

        stot=s111+s011+s101+s110+s001+s010+s100+s000

        DR(i)=DR(i)+stot/8*R(i)
      enddo
      
      return
      end subroutine DStagGN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DStagGNderiv(R,DR,sigma,DAGGER,i) 
      use indices
      implicit none
c     calculates DR = dDStagGN/dsigma_i*R
      real(prc),intent(in) :: R(Nv)
      real(prc),intent(out) :: DR(Nv)
      real(prc),intent(in) :: sigma(Nv)
      logical DAGGER
      integer i
      real(prc) mult
      integer id111,id011,id101,id110,id001,id010,id100,id000

      DR = zero
      id111=i
      id011=id(i,1)
      id101=id(i,2)
      id110=id(i,3)
      id001=id(id011,2)
      id010=id(id011,3)
      id100=id(id101,3)
      id000=id(id100,1)

      DR(id111)=R(id111)/8
      DR(id011)=R(id011)/8
      DR(id101)=R(id101)/8
      DR(id110)=R(id110)/8
      DR(id001)=R(id001)/8
      DR(id010)=R(id010)/8
      DR(id100)=R(id100)/8
      DR(id000)=R(id000)/8
      
      return
      end subroutine DStagGNderiv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcStagGNCondensate(sigma,pbp) 
      use axbmoduleReal
      implicit none
c     calculates <phi.phi> = tr[DStagGN^-1]
      real(prc),intent(in) :: sigma(Nv)
      real(prc),intent(out) :: pbp
      real(prc) :: PNT(Nv),DR(Nv)
      procedure(),pointer :: Mptr => NULL()
      integer i

      pbp=zero
      Mptr => DStagGN
      do i=1,Nv
        PNT=zero
        PNT(i)=one
        call IMnonsymReal(PNT,DR,sigma,.false.,Mptr)
        pbp=pbp+DR(i)
      end do
      pbp=pbp/Nv
      print *,"Exact",pbp

      return
      end subroutine calcStagGNCondensate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcStagGNCondensateNoisy(sigma,pbp) 
      use axbmoduleReal
      use rvmodule
      implicit none
c     calculates <chi.chi> = tr[DStagGN^-1]
      real(prc),intent(in) :: sigma(Nv)
      real(prc),intent(out) :: pbp
      real(prc) :: eta(Nv),DR(Nv)
      procedure(),pointer :: Mptr => NULL()
      integer i,Nnoise

      Nnoise=20
      pbp=zero
      Mptr => DStagGN
      do i=1,Nnoise
        call setGRVs(Nv,eta)
        call IMnonsymReal(eta,DR,sigma,.false.,Mptr)
        pbp=pbp+sum(eta*DR)/Nv
      end do
      pbp=pbp/Nnoise
      print *,"Noisy",pbp

      return
      end subroutine calcStagGNCondensateNoisy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DStagThirring(R,DR,u,DAGGER,mass,add) ! from Debbio/Hands
      use indices
      implicit none
c     calculates DR = DWilson*R
      complex(prc),intent(in) :: R(Nv)
      complex(prc),intent(out) :: DR(Nv)
      complex(prc),intent(in) :: u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      real(prc) mult
      complex(prc) cmass
      integer i,i1,i2,i3,nn(3)
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
      DR = cmass*R

      if (.not.DAGGER) then

        nn(1)=1
        do i=1,Nv
          call ia(i1,i2,i3,i)
          nn(2)=(-1)**i1
          nn(3)=(-1)**(i1+i2)
          do mu=1,NDT
            DR(i)=DR(i)+nn(mu)*u(i,mu)*R(iu(i,mu))
          enddo
        enddo
      
      else ! dagger

        nn(1)=1
        do i=1,Nv
          call ia(i1,i2,i3,i)
          nn(2)=(-1)**i1
          nn(3)=(-1)**(i1+i2)
          do mu=1,NDT
            DR(i)=DR(i)+nn(mu)*conjg(u(i,mu))*R(id(i,mu))
          enddo
        enddo

      endif

      return
      end subroutine DStagThirring
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testStagGN
c      use gaugefield
      use rvmodule
      use axbmodule1
      use axbmoduleReal
      implicit none
      real(prc) :: R(Nv)
      real(prc) :: TMP1(Nv),TMP2(Nv)
      real(prc) :: TMP3(Nv),TMP4(Nv)
      real(prc) :: sigma(Nv)
      procedure(),pointer :: Mptr => NULL()

c      R=one
      call setRealRVs(Nv,R)
      sigma=one

      call DStagGN(R,TMP1,sigma,.false.)
      call DStagGN(R,TMP3,sigma,.true.)

      Mptr => DStagGN
      call IMnonsymReal(TMP1,TMP2,sigma,.false.,Mptr)
      print *,"max err:",maxval(abs(TMP2-R))
      call IMnonsymReal(TMP3,TMP4,sigma,.true.,Mptr)
      print *,"max err:",maxval(abs(TMP4-R))
     
      return
      end subroutine testStagGN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
      end module staggeredmodule
