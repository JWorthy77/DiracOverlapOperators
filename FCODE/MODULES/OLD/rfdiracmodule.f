      module rfdiracmodule
      use arraysizes
      use gammas
      use options
      use ratfuncs
      use wilsonmodule
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DSPF(R,DR,u,DAGGER,SRF)
      implicit none
c     calculates DR = DSPF*R where DSPF is the partial fraction formulation
c     of the sign function of G5.DW
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      type(sgnratfunc) :: SRF
      integer s
      real(prc) p,q,m1,md

      associate(pf => SRF%pfrf%pf,qf=>SRF%pfrf%denom%zeros,
     &           front=>SRF%pfrf%front%coeffs,mult=>SRF%mult)

      md=one
      if(DAGGER) then
        md=-one
      endif

      do s=1,Npf
        p=mult*pf(s)
        q=qf(s)
        call G5DW(R(:,:,2*s-1),DR(:,:,2*s-1),u,DAGGER,-MDW,czero)
        DR(:,:,2*s-1)=DR(:,:,2*s-1)/p + R(:,:,2*s) + md*R(:,:,Ls)
        call G5DW(R(:,:,2*s),DR(:,:,2*s),u,DAGGER,-MDW,czero)
        DR(:,:,2*s)=p*DR(:,:,2*s)/q + R(:,:,2*s-1)
      end do

      DR(:,:,Ls)=czero
      if (allocated(SRF%pfrf%front%coeffs)) then
        if (size(front).ne.1) then
          print *,'DSPF only ready for constant at front of partial frac
     &tion'
          stop
        else
          call G5DW(R(:,:,Ls),DR(:,:,Ls),u,DAGGER,-MDW,czero)
          DR(:,:,Ls)=front(1)*mult*DR(:,:,Ls)
        end if
      end if
      do s=1,Npf
        DR(:,:,Ls)=DR(:,:,Ls) - md*R(:,:,2*s-1)
      end do

      end associate

      return
      end subroutine DSPF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DSPF2(R,DR,u,DAGGER,SRF)
      implicit none
c     calculates DR = DSPF*R where DSPF is the partial fraction formulation
c     of the sign function of V
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      type(sgnratfunc) :: SRF
      complex(prc) :: TMP(Nv,4)
      integer s
      real(prc) p,q,m1,md

      associate(pf => SRF%pfrf%pf,qf=>SRF%pfrf%denom%zeros,
     &           front=>SRF%pfrf%front%coeffs,mult=>SRF%mult)

      md=one
      if(DAGGER) then
        md=-one
      endif

      do s=1,Npf
        p=mult*pf(s)
        q=qf(s)
        call DWilson(R(:,:,2*s-1),DR(:,:,2*s-1),u,DAGGER,-MDW,
     &                                                        czero)
        TMP=R(:,:,2*s)+md*R(:,:,Ls)
        DR(:,:,2*s-1)=DR(:,:,2*s-1)/p + TMP
        call DWilson(R(:,:,2*s),DR(:,:,2*s),u,.not.DAGGER,-MDW,czero)
        TMP=R(:,:,2*s-1)
        DR(:,:,2*s)=p*DR(:,:,2*s)/q + TMP
      end do

      DR(:,:,Ls)=czero
      if (allocated(SRF%pfrf%front%coeffs)) then
        if (size(front).ne.1) then
          print *,'DSPF only ready for constant at front of partial frac
     &tion'
          stop
        else
          call DWilson(R(:,:,Ls),TMP,u,.not.DAGGER,-MDW,czero)
          DR(:,:,Ls)=front(1)*mult*TMP
        end if
      end if
      do s=1,Npf
        TMP=-md*R(:,:,2*s-1)
        DR(:,:,Ls)=DR(:,:,Ls) + TMP
      end do

      end associate

      return
      end subroutine DSPF2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DPF(R,DR,u,DAGGER,SRF)
      implicit none
c     calculates DR = DPF*R where DPF is the partial fraction formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      type(sgnratfunc) :: SRF
      complex(prc) :: TMP(Nv,4)
      integer s
      real(prc) p,q,m1,md

      associate(pf => SRF%pfrf%pf,qf=>SRF%pfrf%denom%zeros,
     &           front=>SRF%pfrf%front%coeffs,mult=>SRF%mult)

      md=one
      if(DAGGER) then
        md=-one
      endif

      do s=1,Npf
        p=(half-baremass/two)*mult*pf(s)
        q=qf(s)
        call G5DW(R(:,:,2*s-1),DR(:,:,2*s-1),u,DAGGER,-MDW,czero)
        DR(:,:,2*s-1) = DR(:,:,2*s-1)/p + 
     &                   R(:,:,2*s) + md*R(:,:,Ls)
        call G5DW(R(:,:,2*s),DR(:,:,2*s),u,DAGGER,-MDW,czero)
        DR(:,:,2*s) = p*DR(:,:,2*s)/q + R(:,:,2*s-1)
      end do

      DR(:,:,Ls)=czero
      if (allocated(SRF%pfrf%front%coeffs)) then
        if (size(front).ne.1) then
          print *,'DSPF only ready for constant at front of partial frac
     &tion'
          stop
        else
          call G5DW(R(:,:,Ls),DR(:,:,Ls),u,DAGGER,-MDW,czero)
          DR(:,:,Ls) = (half-baremass/two)*front(1)*mult*DR(:,:,Ls)
        end if
      end if
      do s=1,Npf
        DR(:,:,Ls) = DR(:,:,Ls) - md*R(:,:,2*s-1)
      end do
      TMP=R(:,:,Ls)
      call mGmu(TMP,5)
      DR(:,:,Ls)=DR(:,:,Ls)+(half+baremass/two)*TMP
      end associate

      return
      end subroutine DPF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DPF2(R,DR,u,DAGGER,SRF)
      implicit none
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      type(sgnratfunc) :: SRF
      complex(prc) :: TMP(Nv,4)
      integer s
      real(prc) p,q,m1,md

      associate(pf => SRF%pfrf%pf,qf=>SRF%pfrf%denom%zeros,
     &           front=>SRF%pfrf%front%coeffs,mult=>SRF%mult)

      md=one
      if(DAGGER) then
        md=-one
      endif

      do s=1,Npf
        p=(half-baremass/two)*mult*pf(s)
        q=qf(s)
        call DWilson(R(:,:,2*s-1),DR(:,:,2*s-1),u,.not.DAGGER,-MDW,
     &                                                        czero)
        TMP=R(:,:,2*s)+md*R(:,:,Ls)
        DR(:,:,2*s-1)=DR(:,:,2*s-1)/p + TMP
        call DWilson(R(:,:,2*s),DR(:,:,2*s),u,DAGGER,-MDW,czero)
        TMP=R(:,:,2*s-1)
        DR(:,:,2*s)=p*DR(:,:,2*s)/q + TMP
      end do

      DR(:,:,Ls)=czero
      if (allocated(SRF%pfrf%front%coeffs)) then
        if (size(front).ne.1) then
          print *,'DSPF only ready for constant at front of partial frac
     &tion'
          stop
        else
          call DWilson(R(:,:,Ls),DR(:,:,Ls),u,DAGGER,-MDW,czero)
          DR(:,:,Ls)=(half-baremass/two)*front(1)*mult*DR(:,:,Ls)
        end if
      end if
      do s=1,Npf
        TMP=-md*R(:,:,2*s-1)
        DR(:,:,Ls)=DR(:,:,Ls) + TMP
      end do
      DR(:,:,Ls)=DR(:,:,Ls)+(half+baremass/two)*R(:,:,Ls)

      end associate

      return
      end subroutine DPF2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDM5(RR,DR,u,DAGGER,Mptr,SRF)
      use axbmodule3
      implicit none
      complex(prc) RR(Nv,4,Ls),DR(Nv,4,Ls)
      complex(prc) u(Nv,3)
      logical DAGGER
      procedure(),pointer :: Mptr
      type(sgnratfunc) :: SRF
      complex(prc) TMP(Nv,4,Ls)
      type(ioptions) iopts

      iopts%SRF=SRF
      if (.not.DAGGER) then
        call Mptr(RR,TMP,u,.not.DAGGER,SRF)
        call IM5dagM5Opts(TMP,DR,u,DAGGER,Mptr,iopts)
      else
        call IM5dagM5Opts(TMP,DR,u,DAGGER,Mptr,iopts)
        call Mptr(RR,TMP,u,.not.DAGGER,SRF)
      end if

      return
      end subroutine IDM5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testSPF()
      use numbers
      use gammas
      use rvmodule
      use gaugefield
      use options
      use sgnmodule
      implicit none
      complex(prc) R(Nv,4,Ls),DR(Nv,4,Ls)
      complex(prc) TMP(Nv,4,Ls)
      complex(prc) R3(Nv,4),DR3(Nv,4)
      procedure(),pointer :: Mptr => NULL()
      type(sgnratfunc) :: SRF
      integer Nht
      
      call setRVs(Nv*4*Ls,R)
      Nht=2*Npf
c      Nht=2*Npf+1

      print *,"test partial fraction inversion"
      print *,'Npf:',Npf,'Nht:',Nht,'Ls:',Ls

      call setHTcoeffs(Nht,SRF)

      Mptr => DSPF
      call Mptr(R,TMP,u,.false.,SRF)
      call IDM5(TMP,DR,u,.false.,Mptr,SRF)
      print *,'DSPF.IDSPF',maxval(abs(R-DR))

      Mptr => DSPF2
      call Mptr(R,TMP,u,.false.,SRF)
      call IDM5(TMP,DR,u,.false.,Mptr,SRF)
      print *,'DSPF2.IDSPF2',maxval(abs(R-DR))

      R=czero
      call setRVs(Nv*4,R(:,:,Ls))
      R3=R(:,:,Ls)
      Mptr => DSPF
      call IDM5(R,DR,u,.false.,Mptr,SRF)
      Mptr => DSPF2
      call IDM5(R,TMP,u,.false.,Mptr,SRF)
      call mGmu(TMP(:,:,Ls),5)
      print *,'IDSPF(:,:,Ls)-IDSPF2(:,:,Ls)',
     &          maxval(abs(DR(:,:,Ls)-TMP(:,:,Ls)))
      call ISGNfactor(R3,DR3,u,.false.,-MDW,SRF)
      print *,'IDSPF(:,:,Ls)-IDSGN',maxval(abs(DR(:,:,Ls)-DR3))
      print *,'IDSPF2(:,:,Ls)-IDSGN',maxval(abs(TMP(:,:,Ls)-DR3))

      return
      end subroutine testSPF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testPF()
      use numbers
      use gammas
      use rvmodule
      use gaugefield
      use options
      use sgnmodule
      use overlapmodule
      use axbmodule3
      implicit none
      complex(prc) R(Nv,4,Ls),DR(Nv,4,Ls)
      complex(prc) TMP(Nv,4,Ls)
      complex(prc) R3(Nv,4),DR3(Nv,4)
      procedure(),pointer :: Mptr => NULL()
      type(sgnratfunc) :: SRF
      integer Nht
      
      call setRVs(Nv*4*Ls,R)
c      Nht=2*Npf
      Nht=2*Npf+1

      print *,"test partial fraction inversion"
      print *,'Npf:',Npf,'Nht:',Nht,'Ls:',Ls

      call setHTcoeffs(Nht,SRF)

      Mptr => DPF
      call Mptr(R,TMP,u,.false.,SRF)
      call IDM5(TMP,DR,u,.false.,Mptr,SRF)
      print *,'DPF.IDPF',maxval(abs(R-DR))

      R=czero
      call setRVs(Nv*4,R(:,:,Ls))
      R3=R(:,:,Ls)
      call mGmu(R(:,:,Ls),5)
      call IDM5(R,DR,u,.false.,Mptr,SRF)
      call IDOverlap(R3,DR3,u,.false.,baremass,SRF)
      print *,'IDPF(:,:,Ls)-IDOL',maxval(abs(DR(:,:,Ls)-DR3))

      stop

      Mptr => DPF2
      call Mptr(R,TMP,u,.false.,SRF)
      call IDM5(TMP,DR,u,.false.,Mptr,SRF)
      print *,'DPF2.IDPF',maxval(abs(R-DR))

      R=czero
      call setRVs(Nv*4,R(:,:,Ls))
      R3=R(:,:,Ls)
      call IDM5(R,DR,u,.false.,Mptr,SRF)
      call IDOverlap(R3,DR3,u,.false.,baremass,SRF)
      print *,'IDPF2(:,:,Ls)-IDOL',maxval(abs(DR(:,:,Ls)-DR3))
      
      return
      end subroutine testPF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module rfdiracmodule
