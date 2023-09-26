      module diracmcmodule
      use arraysizes
      use gammas
      use options
      use ratfuncs
      use basicdiracopsmod
      implicit none
      integer,parameter :: Lv=2*Ls-1
c      integer,parameter :: Lv=Ls
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DHMCSPF(R,DR,u,DAGGER,SRF)
      implicit none
c     calculates DR = DSPF*R where DSPF is the partial fraction formulation
c     of the sign function of G5.DW
      complex(prc),intent(in) :: R(Nv,4,Lv)
      complex(prc),intent(out) :: DR(Nv,4,Lv)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      type(sgnratfunc) :: SRF
      complex(prc) :: TMP(Nv,4)
      integer s,vs
      real(prc) p,q,m1,md

      associate(pf => SRF%pfrf%pf,qf=>SRF%pfrf%denom%zeros,
     &           front=>SRF%pfrf%front%coeffs,mult=>SRF%mult)

      vs=Ls-1

      md=one
      if(DAGGER) then
        md=-one
      endif

      do s=1,Npf
        p=mult*pf(s)
        q=qf(s)
        call G5DW(R(:,:,2*s-1),DR(:,:,2*s-1),u,DAGGER,-MDW,czero)
        DR(:,:,2*s-1)=DR(:,:,2*s-1)/p + R(:,:,2*s) + md*R(:,:,Lv)
        call G5DW(R(:,:,2*s),DR(:,:,2*s),u,DAGGER,-MDW,czero)
        DR(:,:,2*s)=p*DR(:,:,2*s)/q + R(:,:,2*s-1)
      end do
      do s=1,Npf
        p=mult*pf(s)
        q=qf(s)
        call G5DW(R(:,:,vs+2*s-1),DR(:,:,vs+2*s-1),u,
     &                                      .not.DAGGER,-MDW,czero)
        DR(:,:,vs+2*s-1)=DR(:,:,vs+2*s-1)/p + 
     &                         R(:,:,vs+2*s) - md*R(:,:,Lv)
        call G5DW(R(:,:,vs+2*s),DR(:,:,vs+2*s),u,
     &                         .not.DAGGER,-MDW,czero)
        DR(:,:,vs+2*s)=p*DR(:,:,vs+2*s)/q + R(:,:,vs+2*s-1)
      end do

      DR(:,:,Lv)=czero
      if (allocated(SRF%pfrf%front%coeffs)) then
        if (size(front).ne.1) then
          print *,'DHMCSPF only ready for constant at front of partial f
     &raction'
          stop
        else
          call G5DW(R(:,:,Lv),TMP,u,DAGGER,-MDW,czero)
          DR(:,:,Lv)=front(1)*mult*TMP
          call G5DW(R(:,:,Lv),TMP,u,.not.DAGGER,-MDW,czero)
          DR(:,:,Lv)=DR(:,:,Lv)+front(1)*mult*TMP
        end if
      end if
      do s=1,Npf
        DR(:,:,Lv)=DR(:,:,Lv) - md*R(:,:,2*s-1)
      end do
      do s=Npf+1,2*Npf
        DR(:,:,Lv)=DR(:,:,Lv) + md*R(:,:,2*s-1)
      end do

      end associate

      return
      end subroutine DHMCSPF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DHMCSPF2(R,DR,u,DAGGER,SRF)
      implicit none
c     calculates DR = DSPF*R where DSPF is the partial fraction formulation
c     of the sign function of V
      complex(prc),intent(in) :: R(Nv,4,Lv)
      complex(prc),intent(out) :: DR(Nv,4,Lv)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      type(sgnratfunc) :: SRF
      complex(prc) :: TMP(Nv,4)
      integer s,vs
      real(prc) p,q,m1,md

      associate(pf => SRF%pfrf%pf,qf=>SRF%pfrf%denom%zeros,
     &           front=>SRF%pfrf%front%coeffs,mult=>SRF%mult)

      vs=Ls-1

      md=one
      if(DAGGER) then
        md=-one
      endif

      do s=1,Npf
        p=mult*pf(s)
        q=qf(s)
        call DWilson(R(:,:,2*s-1),DR(:,:,2*s-1),u,DAGGER,-MDW,
     &                                                        czero)
        TMP=R(:,:,2*s)+md*R(:,:,Lv)
        DR(:,:,2*s-1)=DR(:,:,2*s-1)/p + TMP
        call DWilson(R(:,:,2*s),DR(:,:,2*s),u,.not.DAGGER,-MDW,czero)
        TMP=R(:,:,2*s-1)
        DR(:,:,2*s)=p*DR(:,:,2*s)/q + TMP
      end do
      do s=1,Npf
        p=mult*pf(s)
        q=qf(s)
        call DWilson(R(:,:,vs+2*s-1),DR(:,:,vs+2*s-1),u,.not.DAGGER,
     &                                                    -MDW,czero)
        TMP=R(:,:,vs+2*s)-md*R(:,:,Lv)
        DR(:,:,vs+2*s-1)=DR(:,:,vs+2*s-1)/p + TMP
        call DWilson(R(:,:,vs+2*s),DR(:,:,vs+2*s),u,DAGGER,-MDW,czero)
        TMP=R(:,:,vs+2*s-1)
        DR(:,:,vs+2*s)=p*DR(:,:,vs+2*s)/q + TMP
      end do

      DR(:,:,Lv)=czero
      if (allocated(SRF%pfrf%front%coeffs)) then
        if (size(front).ne.1) then
          print *,'DSPF only ready for constant at front of partial frac
     &tion'
          stop
        else
          call DWilson(R(:,:,Lv),TMP,u,.not.DAGGER,-MDW,czero)
          DR(:,:,Lv)=front(1)*mult*TMP
          call DWilson(R(:,:,Lv),TMP,u,DAGGER,-MDW,czero)
          DR(:,:,Lv)=DR(:,:,Lv)+front(1)*mult*TMP
        end if
      end if
      do s=1,Npf
        TMP=-md*R(:,:,2*s-1)
        DR(:,:,Lv)=DR(:,:,Lv) + TMP
      end do
      do s=Npf+1,2*Npf
        TMP=md*R(:,:,2*s-1)
        DR(:,:,Lv)=DR(:,:,Lv) + TMP
      end do

      end associate

      return
      end subroutine DHMCSPF2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DHMCPF(R,DR,u,DAGGER,SRF)
      implicit none
c     calculates DR = DPF*R where DPF is the partial fraction formulation
      complex(prc),intent(in) :: R(Nv,4,Lv)
      complex(prc),intent(out) :: DR(Nv,4,Lv)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      type(sgnratfunc) :: SRF
      complex(prc) :: TMP(Nv,4)
      integer s,vs
      real(prc) p,q,m1,md

      associate(pf => SRF%pfrf%pf,qf=>SRF%pfrf%denom%zeros,
     &           front=>SRF%pfrf%front%coeffs,mult=>SRF%mult)

      vs=Ls-1

      md=one
      if(DAGGER) then
        md=-one
      endif

      do s=1,Npf
        p=(half-baremass/two)*(half+baremass/two)*mult*pf(s)
        q=qf(s)
        call G5DW(R(:,:,2*s-1),DR(:,:,2*s-1),u,DAGGER,-MDW,czero)
        DR(:,:,2*s-1) = DR(:,:,2*s-1)/p + 
     &                   R(:,:,2*s) + md*R(:,:,Lv)
        call G5DW(R(:,:,2*s),DR(:,:,2*s),u,DAGGER,-MDW,czero)
        DR(:,:,2*s) = p*DR(:,:,2*s)/q + R(:,:,2*s-1)
      end do
      do s=1,Npf
        p=(half-baremass/two)*(half+baremass/two)*mult*pf(s)
        q=qf(s)
        call G5DW(R(:,:,vs+2*s-1),DR(:,:,vs+2*s-1),u,.not.DAGGER,
     &                                                  -MDW,czero)
        DR(:,:,vs+2*s-1) = DR(:,:,vs+2*s-1)/p + 
     &                   R(:,:,vs+2*s) + md*R(:,:,Lv)
        call G5DW(R(:,:,vs+2*s),DR(:,:,vs+2*s),u,.not.DAGGER,-MDW,czero)
        DR(:,:,vs+2*s) = p*DR(:,:,vs+2*s)/q + R(:,:,vs+2*s-1)
      end do


      DR(:,:,Lv)=czero
      if (allocated(SRF%pfrf%front%coeffs)) then
        if (size(front).ne.1) then
          print *,'DSPF only ready for constant at front of partial frac
     &tion'
          stop
        else
          call G5DW(R(:,:,Lv),DR(:,:,Lv),u,DAGGER,-MDW,czero)
          DR(:,:,Lv) = (half-baremass/two)*front(1)*mult*DR(:,:,Lv)
        end if
      end if
      do s=1,2*Npf
        DR(:,:,Lv) = DR(:,:,Lv) - md*R(:,:,2*s-1)
      end do
      TMP=R(:,:,Lv)
      call mGmu(TMP,5)
      DR(:,:,Lv)=DR(:,:,Lv)+(half+baremass/two)*(half+baremass/two)*TMP
      end associate

      return
      end subroutine DHMCPF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DHMCPF2(R,DR,u,DAGGER,SRF)
      implicit none
      complex(prc),intent(in) :: R(Nv,4,Lv)
      complex(prc),intent(out) :: DR(Nv,4,Lv)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      type(sgnratfunc) :: SRF
      complex(prc) :: TMP(Nv,4)
      integer s,vs
      real(prc) p,q,m1,md

      associate(pf => SRF%pfrf%pf,qf=>SRF%pfrf%denom%zeros,
     &           front=>SRF%pfrf%front%coeffs,mult=>SRF%mult)

      vs=Ls-1

      md=one
      if(DAGGER) then
        md=-one
      endif

      do s=1,Npf
        p=(half-baremass/two)*(half+baremass/two)*mult*pf(s)
        q=qf(s)
        call DWilson(R(:,:,2*s-1),DR(:,:,2*s-1),u,.not.DAGGER,-MDW,
     &                                                        czero)
        TMP=R(:,:,2*s)+md*R(:,:,Lv)
        DR(:,:,2*s-1)=DR(:,:,2*s-1)/p + TMP
        call DWilson(R(:,:,2*s),DR(:,:,2*s),u,DAGGER,-MDW,czero)
        TMP=R(:,:,2*s-1)
        DR(:,:,2*s)=p*DR(:,:,2*s)/q + TMP
      end do
      do s=1,Npf
        p=(half-baremass/two)*(half+baremass/two)*mult*pf(s)
        q=qf(s)
        call DWilson(R(:,:,vs+2*s-1),DR(:,:,vs+2*s-1),u,DAGGER,-MDW,
     &                                                        czero)
        TMP=R(:,:,vs+2*s)+md*R(:,:,Lv)
        DR(:,:,vs+2*s-1)=DR(:,:,vs+2*s-1)/p + TMP
        call DWilson(R(:,:,vs+2*s),DR(:,:,vs+2*s),u,.not.DAGGER,
     &                                                    -MDW,czero)
        TMP=R(:,:,vs+2*s-1)
        DR(:,:,vs+2*s)=p*DR(:,:,vs+2*s)/q + TMP
      end do

      DR(:,:,Lv)=czero
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
      do s=1,2*Npf
        TMP=-md*R(:,:,2*s-1)
        DR(:,:,Lv)=DR(:,:,Lv) + TMP
      end do
      DR(:,:,Lv)=DR(:,:,Lv)+(half+baremass/two)*R(:,:,Lv)

      end associate

      return
      end subroutine DHMCPF2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDM5(RR,DR,u,DAGGER,Mptr,SRF)
      use axbmodule3
      implicit none
      complex(prc) RR(Nv,4,Lv),DR(Nv,4,Lv)
      complex(prc) u(Nv,3)
      logical DAGGER
      procedure(),pointer :: Mptr
      type(sgnratfunc) :: SRF
      complex(prc) TMP(Nv,4,Lv)
      type(ioptions) iopts

      iopts%SRF=SRF
      if (.not.DAGGER) then
        call Mptr(RR,TMP,u,.not.DAGGER,SRF)
        call IM5dagM5OptsVS(4*Nv*Lv,TMP,DR,u,DAGGER,Mptr,iopts)
      else
        call IM5dagM5OptsVS(4*Nv*Lv,TMP,DR,u,DAGGER,Mptr,iopts)
        call Mptr(RR,TMP,u,.not.DAGGER,SRF)
      end if

      return
      end subroutine IDM5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testHMCSPF()
      use numbers
      use gammas
      use rvmodule
      use gaugefield
      use options
c      use sgnmodule
      implicit none
      complex(prc) R(Nv,4,Lv),DR(Nv,4,Lv)
      complex(prc) TMP(Nv,4,Lv)
      complex(prc) R3(Nv,4),DR3(Nv,4),TMP3(Nv,4)
      procedure(),pointer :: Mptr => NULL()
      type(sgnratfunc) :: SRF
      integer Nht
      
      call setRVs(Nv*4*Lv,R)
c      Nht=2*Npf
      Nht=2*Npf+1

      print *,"test partial fraction inversion"
      print *,'Npf:',Npf,'Nht:',Nht,'Ls:',Ls

      call setHTcoeffs(Nht,SRF)
c      call setZoloCoeffs(Nht,0.1,5.0,SRF)

c      call setRVs(Nv*4,R3)
c      call SGNfactor(R3,TMP3,u,.false.,SRF)
c      call SGNfactor(R3,DR3,u,.true.,SRF)
c      print *,'error:',maxval(abs(TMP3-DR3))

c      Mptr => DHMCSPF
c      call Mptr(R,TMP,u,.false.,SRF)
c      call IDM5(TMP,DR,u,.false.,Mptr,SRF)
c      print *,'DSPF.IDSPF',maxval(abs(R-DR))

c      R=czero
c      call setRVs(Nv*4,R(:,:,Lv))
c      call IDM5(R,DR,u,.false.,Mptr,SRF)
c      R3=DR(:,:,Lv)
c      call SGNfactor(R3,TMP3,u,.false.,SRF)
c      call SGNfactor(R3,DR3,u,.true.,SRF)
c      print *,'error:',maxval(abs(R(:,:,Lv)-TMP3-DR3))

      Mptr => DHMCSPF2
      call Mptr(R,TMP,u,.false.,SRF)
      call IDM5(TMP,DR,u,.false.,Mptr,SRF)
      print *,'DHMCSPF2.IDHMCSPF2',maxval(abs(R-DR))

c      R=czero
c      call setRVs(Nv*4,R(:,:,Ls))
c      R3=R(:,:,Ls)
c      Mptr => DHMCSPF
c      call IDM5(R,DR,u,.false.,Mptr,SRF)
c      Mptr => DSPF2
c      call IDM5(R,TMP,u,.false.,Mptr,SRF)
c      call mGmu(TMP(:,:,Ls),5)
c      print *,'IDSPF(:,:,Ls)-IDSPF2(:,:,Ls)',
c     &          maxval(abs(DR(:,:,Ls)-TMP(:,:,Ls)))
c      call ISGNfactor(R3,DR3,u,.false.,SRF)
c      print *,'IDSPF(:,:,Ls)-IDSGN',maxval(abs(DR(:,:,Ls)-DR3))
c      print *,'IDSPF2(:,:,Ls)-IDSGN',maxval(abs(TMP(:,:,Ls)-DR3))

      return
      end subroutine testHMCSPF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testHMCPF()
      use numbers
      use gammas
      use rvmodule
      use gaugefield
      use options
!      use sgnmodule
      use overlapmoduledev
      use axbmodule3
      implicit none
      complex(prc) R(Nv,4,Lv),DR(Nv,4,Lv)
      complex(prc) TMP(Nv,4,Lv)
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

      Mptr => DHMCPF
      call Mptr(R,TMP,u,.false.,SRF)
      call IDM5(TMP,DR,u,.false.,Mptr,SRF)
      print *,'DHMCPF.IDHMCPF',maxval(abs(R-DR))

      stop

c      R=czero
c      call setRVs(Nv*4,R(:,:,Ls))
c      R3=R(:,:,Ls)
c      call mGmu(R(:,:,Ls),5)
c      call IDM5(R,DR,u,.false.,Mptr,SRF)
c      call IDOverlap(R3,DR3,u,.false.,SRF)
c      print *,'IDPF(:,:,Ls)-IDOL',maxval(abs(DR(:,:,Ls)-DR3))

c      stop

c      Mptr => DPF2
c      call Mptr(R,TMP,u,.false.,SRF)
c      call IDM5(TMP,DR,u,.false.,Mptr,SRF)
c      print *,'DPF2.IDPF',maxval(abs(R-DR))

c      R=czero
c      call setRVs(Nv*4,R(:,:,Ls))
c      R3=R(:,:,Ls)
c      call IDM5(R,DR,u,.false.,Mptr,SRF)
c      call IDOverlap(R3,DR3,u,.false.,SRF)
c      print *,'IDPF2(:,:,Ls)-IDOL',maxval(abs(DR(:,:,Ls)-DR3))
      
      return
      end subroutine testHMCPF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module diracmcmodule
