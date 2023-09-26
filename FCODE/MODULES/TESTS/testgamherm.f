      module testgamhermmod
      use arraysizes
      use numbers
      use options
      use gaugefield
      use rvmodule
      use ratfuncs
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testGammaHermiticity()
!     test gamma hermiticity of Dirac operators
      use gammas
      implicit none
      complex(prc),dimension(4*Nv) :: R,DR
      complex(prc),dimension(4*Nv,Ls) :: R5,DR5
      real(prc) :: err
      type(sgnratfunc) :: HRF
      type(sgnratfunc) :: ZRF

      call setRVs(Nv*4,R)
      R5=czero ; R5(:,1)=R
      call setRVs(Nv*4*Ls,R5)
      call setHTcoeffs(10,HRF)
c      call setZoloCoeffs(10,ZRF,lmin,lmax)
      
      if (VERBOSE.eq.1) then
        print *,'Gamma Hermiticity: g.D.g=D^dagger'
      endif

      goto 1000

      ! test Wilson operator - independent of mass formulation
      if (VERBOSE.eq.1) then ; print *,'Wilson Dirac operator' ; endif
      call testGH_DOL(R,err,1) 

      ! test Wilson operator - independent of mass formulation
      if (VERBOSE.eq.1) then ; print *,'Shamir (H) Dirac operator';endif
      call testGH_DOL(R,err,2) 

      ! test V(Wilson) partial fraction - independent of mass formulation
      if (VERBOSE.eq.1) then ; print *,'VOLpf(Wilson) operator' ; endif
      call testGH_DOL(R,err,3,HRF) 

      ! test V(Wilson) multiplicative - independent of mass formulation
      if (VERBOSE.eq.1) then ; print *,'VOLfac(Wilson) operator' ; endif
      call testGH_DOL(R,err,4,HRF) 

      ! test DOL(Wilson) - M1
      MTYPE=1
      if (RFTYPE.eq.1) then
        if (VERBOSE.eq.1) then ; print *,'DOL(Wilson) VOLpf M1' ; endif
      elseif (RFTYPE.eq.2) then
        if (VERBOSE.eq.1) then ; print *,'DOL(Wilson) VOLfac M1' ; endif
      endif
      call testGH_DOL(R,err,5,HRF) 

      ! test DOL(Wilson) - M3
      MTYPE=3
      if (RFTYPE.eq.1) then
        if (VERBOSE.eq.1) then ; print *,'DOL(Wilson) VOLpf M3' ; endif
      elseif (RFTYPE.eq.2) then
        if (VERBOSE.eq.1) then ; print *,'DOL(Wilson) VOLfac M3' ; endif
      endif
      call testGH_DOL(R,err,5,HRF) 

      ! test DOL(Shamir) - M1
      MTYPE=1
c      if (RFTYPE.eq.1) then
        if (VERBOSE.eq.1) then ; print *,'DOL(Shamir) VOLpf M1' ; endif
c      elseif (RFTYPE.eq.2) then
c        if (VERBOSE.eq.1) then ; print *,'DOL(Shamir) VOLfac M1' ; endif
c      endif
      call testGH_DOL(R,err,6,HRF) 

      ! test DOL(Shamir) - M3
      MTYPE=3
c      if (RFTYPE.eq.1) then
        if (VERBOSE.eq.1) then ; print *,'DOL(Shamir) VOLpf M3' ; endif
c      elseif (RFTYPE.eq.2) then
c        if (VERBOSE.eq.1) then ; print *,'DOL(Shamir) VOLfac M3' ; endif
c      endif
      call testGH_DOL(R,err,7,HRF) 


1000  continue

      ! test DDW(Shamir) - M1
      MTYPE=1 ; DWkernel=1
      if (VERBOSE.eq.1) then ; print *,'DDW(Shamir) M1' ; endif
      call testGH_DDW(R,R5,err,1) 

      ! test DDW(Wilson) - M1
      MTYPE=1 ; DWkernel=2
      if (VERBOSE.eq.1) then ; print *,'DDW(Wilson) M1' ; endif
      call testGH_DDW(R,R5,err,1) 

      return
      end subroutine testGammaHermiticity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testGH_DOL(R,err,GHTEST,SRF)
!     test gamma hermiticity of Overlap Dirac operators
      use gammas
      use basicdiracopsmod
      use overlapmoduledev
c      use shamirmodule
      implicit none
      complex(prc),dimension(4*Nv),intent(in) :: R
      real(prc),intent(out) :: err
      integer,intent(in) :: GHTEST
      type(sgnratfunc),intent(in),optional :: SRF
      complex(prc),dimension(4*Nv) :: TMP,DRL,DRR,DR,FIX

!     G3.D.G3
      TMP=R
      call mGmu(TMP,4)
      if (GHTEST.eq.1) then
        call DWilson(TMP,DRL,u,.false.,baremass,czero)
      elseif (GHTEST.eq.2) then
        call DShamir(TMP,DRL,u,.false.,baremass,czero)
      elseif (GHTEST.eq.3) then
        call VOLpf(TMP,DRL,u,.false.,-MDW,SRF)
      elseif (GHTEST.eq.4) then
        call VOLfac(TMP,DRL,u,.false.,-MDW,SRF)
      elseif (GHTEST.eq.5) then
        call Doverlap(TMP,DRL,u,.false.,baremass,SRF)
      elseif (GHTEST.eq.6) then
        call DOLS(TMP,DRL,u,.false.,baremass,SRF)
      endif
      call mGmu(DRL,4)

!     D^dagger
      if (GHTEST.eq.1) then
        call DWilson(R,DRR,u,.true.,baremass,czero)
      elseif (GHTEST.eq.2) then
        call DShamir(R,DRR,u,.true.,baremass,czero)
      elseif (GHTEST.eq.3) then
        call VOLpf(R,DRR,u,.true.,-MDW,SRF)
      elseif (GHTEST.eq.4) then
        call VOLfac(R,DRR,u,.true.,-MDW,SRF)
      elseif (GHTEST.eq.5) then
        call Doverlap(R,DRR,u,.true.,baremass,SRF)
        if (MTYPE.eq.3) then
          call VOLfac(R,TMP,u,.false.,-MDW,SRF)
          call VOLfac(R,DR,u,.true.,-MDW,SRF)
          FIX=zi*baremass*(2*R-TMP-DR)/2
          call mGmu(FIX,4)
        endif
      elseif (GHTEST.eq.6) then
        call DOLS(R,DRR,u,.true.,baremass,SRF)
        if (MTYPE.eq.3) then
          call VOLSpf(R,TMP,u,.false.,-MDW,SRF)
          call VOLSpf(R,DR,u,.true.,-MDW,SRF)
          FIX=zi*baremass*(2*R-TMP-DR)/2
          call mGmu(FIX,4)
        endif
      endif
    
!     G3.D.G3 - D^dagger
      DR=DRL-DRR
      if (MTYPE.eq.3) then
!       G3.D.G3 - D^dagger - im.G3.(2-V-Vdag)/2
        DR=DR-FIX
      end if
      err=maxval(abs(DR))
      print *,"max err:",err
      return
      end subroutine testGH_DOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testGH_DDW(R,R5,err,GHTEST)
!     test gamma hermiticity of Domain Wall Dirac operators
      use gammas
      use domainwallmod
      implicit none
      complex(prc),dimension(4*Nv),intent(in) :: R
      complex(prc),dimension(4*Nv,Ls),intent(in) :: R5
      real(prc),intent(out) :: err
      integer,intent(in) :: GHTEST
      complex(prc),dimension(4*Nv) :: TMP,DRL,DRR,DR,FIX
      complex(prc),dimension(4*Nv,Ls) :: TMP5,DRL5,DRR5,DR5,FIX5

!     G3.D.G3
!      TMP=R
      TMP5=R5
!      call mGmu(TMP,4)
      call mGmu5(TMP5,4)
      if (GHTEST.eq.1) then
        call DDW(TMP5,DRL5,u,.false.,baremass)
!      elseif (GHTEST.eq.3) then
!        call KDDW4(TMP,DRL,u,.false.,baremass)
      endif
!      call mGmu(DRL,4)
      call mGmu5(DRL5,4)

!     D^dagger
      if (GHTEST.eq.1) then
        call DDW(R5,DRR5,u,.true.,baremass)
!      elseif (GHTEST.eq.3) then
!        call KDDW4(R,DRR,u,.true.,baremass)
      endif
    
!     G3.D.G3 - D^dagger
!      DR=DRL-DRR
      DR5=DRL5-DRR5
      if (MTYPE.eq.3) then
!       G3.D.G3 - D^dagger - im.G3.(2-V-Vdag)/2
        DR=DR-FIX
      end if
      err=maxval(abs(DR5))
      print *,"max err:",err
      return
      end subroutine testGH_DDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testgamhermmod
