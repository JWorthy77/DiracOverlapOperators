      module chirality
!     contains routines to test the chirality of the dirac operators
      use arraysizes
      use numbers
      use options
      use rvmodule
      use ratfuncs
      use arpackmodule
      use gaugefield
      use overlapmoduledev
      use shamirmodule
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDWchirality()
      implicit none
      logical ZOLO
      integer Nht
      type(sgnratfunc) :: SRF
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) DR1(Nv,4),DR2(Nv,4)
      procedure(),pointer :: Sptr => NULL()
      complex(prc) mtest(Nv,4,200)
      real ev(1)
      real(prc) lmin,lmax
      integer j,jmin,jmax,shft

      print *,'test convergence of VOL operators'
      if (ZOLO) then
        call calceigs('SM',1,ev,1)
        print *,'ev min:',ev
        lmin=ev(1)
        call calceigs('LM',1,ev,1)
        print *,'ev max:',ev
        lmax=ev(1)
      end if

      Sptr => DOverlap
      jmin=6
      jmax=40
!     psi=delta
!      R=czero
!      R(1,1)=cone
      call setRVs(Nv*4,R)
      if (.not.ZOLO) then
        shft=1
        call setHTcoeffs(jmin-1,SRF)
      elseif(ZOLO) then
        shft=2
        call setZoloCoeffs((jmin-2)/2,SRF,lmin,lmax)
      endif
      call Sptr(R,mtest(:,:,jmin-shft),u,.false.,SRF)
      do j=jmin,jmax,shft
        if (.not.ZOLO) then
          call setHTcoeffs(j,SRF)
        elseif(ZOLO) then
          call setZoloCoeffs(j/2,SRF,lmin,lmax)
        endif
        call Sptr(R,mtest(:,:,j),u,.false.,SRF)
      end do

      if (.not.ZOLO) then
        open(unit=11,file='convHT.dat',status='unknown',
     &                                               form='formatted')
      elseif (ZOLO) then
        open(unit=11,file='convZolo.dat',status='unknown',
     &                                               form='formatted')
      end if

      do j=jmin,jmax,shft
        print *,j,maxval(abs(mtest(:,:,j)-mtest(:,:,j-shft)))
        write(11,*) j,maxval(abs(mtest(:,:,j)-mtest(:,:,j-shft)))
      end do

      return
      end subroutine testDWchirality
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module chirality
