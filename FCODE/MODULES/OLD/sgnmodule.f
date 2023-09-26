      module sgnmodule
      use pacc
      use arraysizes
      use ratfuncs
      use wilsonmodule
      implicit none

      private :: D_dagD

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine G5DW2pC(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = (G5DW^2+C).R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) :: TMP(Nv,4)

      call G5DW(R,TMP,u,DAGGER,mass,czero)
      call G5DW(TMP,DR,u,DAGGER,mass,czero)
      DR=DR+add*R
      
      return
      end subroutine G5DW2pC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IG5DW2pC(R,DR,u,DAGGER,mass,add)
      use axbmodule1
      implicit none
c     calculates DR = (G5DW^2+C).R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr => NULL()

      Mptr => G5DW2pC
      call IMnonsym(R,DR,u,DAGGER,mass,add,Mptr)

      return
      end subroutine IG5DW2pC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine D_dagD(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DWilson*R
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,4)

      call G5DW2pC(R,TMP,u,DAGGER,mass,add)
      call G5DW2pC(TMP,DR,u,.not.DAGGER,mass,add)

      return
      end subroutine D_dagD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SGNfactor(RR,S,u,DAGGER,dwmass,SRF)
!     approximate SGN(g5.Dw(-M)) using factor form
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,nn,nd
c      procedure(),pointer :: Mptr => NULL()

c      print *,'SGNfactor'
c      call printSRatFunc(SRF) 

      nn=size(SRF%frf%num%zeros)
      nd=size(SRF%frf%denom%zeros)
c      print *,nn,nd

      if (DAGGER) then
        call G5DW(RR,TMP1,u,DAGGER,dwmass,czero)
      else
        TMP1=RR
      end if
c      print *,'SGNfactor B'
      do j=1,nn
!       numerator
        add = -SRF%frf%num%zeros(j)
c        print *,j,add
        call G5DW2pC(TMP1,TMP2,u,DAGGER,dwmass,add)
!       denominator
c      print *,'SGNfactor BA'
        add = -SRF%frf%denom%zeros(j)
        call IG5DW2pC(TMP2,TMP1,u,DAGGER,dwmass,add)
      end do
c      print *,'SGNfactor C'
      do j=nn+1,nd
!       denominator
        add = -SRF%frf%denom%zeros(j)
        call IG5DW2pC(TMP1,TMP2,u,DAGGER,dwmass,add)
        TMP1=TMP2
      end do
      if (.not.DAGGER) then
        call G5DW(TMP1,S,u,DAGGER,dwmass,czero)
      else
        S=TMP1
      endif
      S=SRF%mult*S

c      print *,'SGNfactor Z'
      return
      end subroutine SGNfactor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ISGNfactor(RR,S,u,DAGGER,dwmass,SRF)
!     approximate SGN(g5.Dw(-M)) using factor form
      use options
      use axbmodule2
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Mptr => NULL()
      type(ioptions) :: iopts

      print *,'ISGNfactor'
      Mptr => SGNfactor
      iopts%SRF=SRF
      iopts%mass=dwmass
      call IMnonsymOpts(RR,S,u,DAGGER,Mptr,iopts)

      return
      end subroutine ISGNfactor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testSGN(ZOLO)
      use numbers
      use options
      use gaugefield
      use axbmodule1
      use arpackmodule
      implicit none
      logical ZOLO
      type(sgnratfunc) :: SRF
      complex(prc) psi(Nv,4),S(Nv,4)
      complex(prc) R(Nv,4),DR(Nv,4),TMP(Nv,4)
      real(prc) mass,MS
      complex(prc) mtest(Nv,4,200)
      procedure(),pointer :: Sptr => NULL()
      integer j,jmin,jmax
      real ev
      real(prc) lmin,lmax

      print *,'test basic SGN operators'

      R=cone
      call G5DW(R,TMP,u,.false.,-MDW,czero)
      call IG5DW(TMP,DR,u,.false.,-MDW,czero)
      print *,maxval(abs(R-DR))
      R=cone
      Sptr=>G5DW2pC
      call G5DW2pC(R,TMP,u,.false.,-MDW,cone)
      call IG5DW2pC(TMP,DR,u,.false.,-MDW,cone)
c      Sptr=>G5DW2pC
c      call IMnonsym(TMP,DR,u,.false.,-MDW,cone,Sptr)
      print *,maxval(abs(R-DR))

      R=czero
      R(1,1)=cone
      call setHTcoeffs(20,SRF)
      call SGNfactor(R,TMP,u,.false.,-MDW,SRF)
      call SGNfactor(TMP,DR,u,.false.,-MDW,SRF)
      print *,'SGN^2',maxval(abs(DR-R))

      call setHTcoeffs(20,SRF)
      call SGNfactor(R,TMP,u,.false.,-MDW,SRF)
      call ISGNfactor(TMP,DR,u,.false.,-MDW,SRF)
      print *,'SGN.ISGN',maxval(abs(R-DR))

      return
      end subroutine testSGN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testSGNConv(ZOLO)
      use numbers
      use options
      use gaugefield
      use axbmodule1
      use arpackmodule
      implicit none
      logical ZOLO
      type(sgnratfunc) :: SRF
      complex(prc) psi(Nv,4),S(Nv,4)
      complex(prc) R(Nv,4),DR(Nv,4),TMP(Nv,4)
      real(prc) mass,MS
      complex(prc) mtest(Nv,4,200)
      procedure(),pointer :: Sptr => NULL()
      integer j,jmin,jmax
      real ev(1)
      real(prc) lmin,lmax

      print *,'check convergence of SGN operators'

      if (.not.ZOLO) then
        open(unit=10,file='onesHT.dat',status='unknown',
     &                                               form='formatted')
      elseif (ZOLO) then
        call calceigs('SM',1,ev,1,Nv)
        print *,'ev min:',ev
        lmin=ev(1)
        call calceigs('LM',1,ev,1,Nv)
        print *,'ev max:',ev
        lmax=ev(1)
        open(unit=10,file='onesZolo.dat',status='unknown',
     &                                               form='formatted')
      end if


      stop

      Sptr => SGNfactor
      jmin=4
      jmax=4

!     psi=1
      psi=cone
!     psi=delta
      psi=czero
      psi(1,1)=cone
      if (.not.ZOLO) then
        call setHTcoeffs(jmin-1,SRF)
      elseif(ZOLO) then
        call setZoloCoeffs(jmin-1,SRF,lmin,lmax)
      endif
      call Sptr(psi,mtest(:,:,jmin-1),u,.false.,-MDW,SRF)
      do j=jmin,jmax
        if (.not.ZOLO) then
          call setHTcoeffs(j,SRF)
        elseif(ZOLO) then
          call setZoloCoeffs(j,SRF,lmin,lmax)
        endif
        call Sptr(psi,mtest(:,:,j),u,.false.,-MDW,SRF)
        print *,j,maxval(abs(mtest(:,:,j)-mtest(:,:,j-1)))
        write(10,*) j,maxval(abs(mtest(:,:,j)-mtest(:,:,j-1)))
      end do

      close(10)

      return
      end subroutine testSGNConv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine compareSGN()
      use numbers
      use options
      use gaugefield
      use axbmodule1
      use arpackmodule
      implicit none
      type(sgnratfunc) :: RFHT,RFZ
      complex(prc) R(Nv,4),DRHT(Nv,4),DRZ(Nv,4)
      procedure(),pointer :: Sptr => NULL()
      real ev(1)
      real(prc) lmin,lmax

      print *,'compare HT and Zolo SGN operators'

      R=cone
      Sptr => SGNfactor
      call setHTcoeffs(10,RFHT)
c      call printSRatFunc(RFHT)
c      stop 
      call calceigs('SM',1,ev,1,Nv)
      print *,'ev min:',ev
      lmin=ev(1)
c      lmin=0.1
      call calceigs('LM',1,ev,1,Nv)
      print *,'ev max:',ev
      lmax=ev(1)
c      lmax=10
      call setZoloCoeffs(5,RFZ,0.5*lmin,2*lmax)
c      call printSRatFunc(RFZ)
c      stop 
      call Sptr(R,DRHT,u,.false.,-MDW,RFHT)
c      print *,'SGN HT done'
c      stop
      call Sptr(R,DRZ,u,.false.,-MDW,RFZ)
      print *,'lmin',lmin
      print *,'lmax',lmax
      print *,'Z and HT difference:',maxval(abs(DRHT-DRZ))
      return
      end subroutine compareSGN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module sgnmodule
