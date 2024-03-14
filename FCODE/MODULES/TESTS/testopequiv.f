      module testDequivmod
      use pacc
      use arraysizes
      use numbers
      use options
      use basicdiracopsmod
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOperatorEquivalence()
      implicit none

!      call equivShamir()
!      call equivVOLW()
!      call equivVOLS()
c      call equivDOL()
      call equivDOL_DDW()

      return
      end subroutine testOperatorEquivalence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine equivShamir()
      use overlapmoduledev
      use polycoeffsmod
      use ratfuncs
      use gaugefield
      use rvmodule
      use spectra
      use kernelspectrarange
      implicit none
      complex(prc),dimension(Nv,4) :: R,DR1,DR2
      real(prc) :: err,lmin,lmax

      print *,"Check equivalence of Shamir polynomial formulation"
      print *,"Calc Wilson eigenvalues"
      dwkernel=1
      call estimateKernelExtrema(10,1,lmax)
      call estimateKernelExtrema(10,-1,lmin)
      print *,"lmin:",lmin
      print *,"lmax:",lmax
      call setRVs(Nv*4,R)
      call setPolyCoeffs(pfcs,1,14,real(0.0,prc),real(0.0,prc)) ! Shamir - 1/(2+x) [1,5]
      if(VERBOSE.eq.1)then ; print *,'Shamir-ShamirPoly' ; endif
      call DShamir(R,DR1,u,.false.,baremass,czero)
      print *,maxval(abs(DR1))
      call DShamirPoly(R,DR2,u,.false.,baremass,czero)
      print *,maxval(abs(DR2))
      err=maxval(abs(DR1-DR2))
      print *,'err(!DAGGER) ',err
      call DShamir(R,DR1,u,.true.,baremass,czero)
      call DShamirPoly(R,DR2,u,.true.,baremass,czero)
      err=maxval(abs(DR1-DR2))
      print *,'err(DAGGER) ',err

      return
      end subroutine equivShamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine equivVOLW()
      use overlapmoduledev
      use polycoeffsmod
      use ratfuncs
      use gaugefield
      use rvmodule
      use spectra
      implicit none
      complex(prc),dimension(Nv,4) :: R,DR1,DR2,DR3,DR4
      real(prc) :: err1,err2,err3,lmin,lmax
      type(sgnratfunc) :: HT,PF,EPF
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: Dptr => NULL()

      print *,"Check equivalence of VOLW multishift form"
      print *,"Calc Wilson eigenvalues"
!      call calcDEigs(10,lmin,lmax,2,.false.)
      print *,"lmin:",lmin
      print *,"lmax:",lmax
      call setRVs(Nv*4,R)
!      call setPolyRFCoeffs(PF,2,10,real(0.0,prc),real(0.0,prc)) ! 2:  x^-0.5 [1,5] 
c      call setPolyRFCoeffs(PF,3,10,real(0.0,prc),real(0.0,prc)) ! 2:  x^-0.5 [1,25] 
c      call setEPolyRFCoeffs(EPF,3,-3,8,real(0.0,prc),real(0.0,prc)) ! 3: x^-0.5 [1,25]
c      print *,PF%pfrf%front%coeffs
      call setHTcoeffs(20,HT) 

      if(VERBOSE.eq.1)then ; print *,'VOLpf-VOLpoly' ; endif
      call VOLpf(R,DR1,u,.false.,baremass,HT)
      print *,maxval(abs(DR1))
      Mptr => DdagD
      Dptr => DWilson
      call VOLMpf(R,DR2,u,.false.,baremass,HT,Mptr,Dptr)
      print *,maxval(abs(DR2))
c      call VOLpoly(R,DR2,u,.false.,baremass,PF)
c      print *,maxval(abs(DR3))
c      call VOLepoly(R,DR3,u,.false.,baremass,EPF)
c      print *,maxval(abs(DR4))
      err1=maxval(abs(DR1-DR2))
c      err2=maxval(abs(DR1-DR3))
c      err3=maxval(abs(DR1-DR4))
      print *,'err(!DAGGER - Mpf) ',err1
c      print *,'err(!DAGGER - poly) ',err2
c      print *,'err(!DAGGER - epoly) ',err3

!      call VOLpf(R,DR1,u,.true.,baremass,HT)
!      call VOLpoly(R,DR2,u,.true.,baremass,PF)
!      err=maxval(abs(DR1-DR2))
!      print *,'err(DAGGER) ',err

      return
      end subroutine equivVOLW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine equivVOLS()
      use overlapmoduledev
      use polycoeffsmod
      use ratfuncs
      use gaugefield
      use rvmodule
      use spectra
      implicit none
      complex(prc),dimension(Nv,4) :: R,DR1,DR2,DR3
      real(prc) :: err1,err2,lmin,lmax
      type(sgnratfunc) :: HT,PF,EPF

      print *,"Check equivalence of VOLS polynomial formulation"
      print *,"Calc Wilson eigenvalues"
!      call calcDEigs(10,lmin,lmax,2,.false.)
      print *,"lmin:",lmin
      print *,"lmax:",lmax
      call setRVs(Nv*4,R)
!      call setPolyRFCoeffs(PF,2,10,real(0.0,prc),real(0.0,prc)) ! 2:  x^-0.5 [1,5] 
      call setPolyRFCoeffs(PF,3,10,real(0.0,prc),real(0.0,prc)) ! 3:  x^-0.5 [1,25] 
      call setPolyCoeffs(pfcs,1,10,real(0.0,prc),real(0.0,prc)) ! Shamir - 1/(2+x) [1,5]
      call setEPolyRFCoeffs(EPF,3,-3,8,real(0.0,prc),real(0.0,prc)) ! 3: x^-0.5 [1,25]
      print *,PF%pfrf%front%coeffs
      print *,EPF%epoly%cfs
      call setHTcoeffs(20,HT) 

      if(VERBOSE.eq.1)then ; print *,'VOLSpf-VOLSpoly' ; endif
c      call VOLSpf(R,DR1,u,.false.,baremass,HT)
      print *,maxval(abs(DR1))
c      call VOLSpoly(R,DR2,u,.false.,baremass,PF)
      print *,maxval(abs(DR2))
      call VOLSepoly(R,DR3,u,.false.,baremass,EPF)
      print *,maxval(abs(DR3))
      err1=maxval(abs(DR1-DR2))
      err2=maxval(abs(DR1-DR3))
      print *,'err(!DAGGER) poly',err1
      print *,'err(!DAGGER) epoly',err2

c      call VOLSpf(R,DR1,u,.true.,baremass,HT)
c      call VOLSpoly(R,DR2,u,.true.,baremass,PF)
c      err=maxval(abs(DR1-DR2))
c      print *,'err(DAGGER) ',err

      return
      end subroutine equivVOLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine equivDOL()
      use gammas
      use overlapmoduledev
      use domainwallmod
      use rvmodule
      use gaugefield
      implicit none
      complex(prc),dimension(Nv,4) :: R,DR1,DR2
      real(prc) :: err
      type(sgnratfunc) :: SRF
      integer i1,i2

      call setRVs(Nv*4,R)
      call setHTcoeffs(20,SRF)
!      call setZoloCoeffs(10,ZRF,lmin,lmax)

      print *,""
      print *,"TEST EQUIVALENCE OF OVERLAP FORMULATIONS"
      print *,""

      MTYPE=1
      dwkernel=2
      call Doverlap(R,DR1,u,.false.,baremass,SRF)
      call DOLop(R,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'Doverlap-DOL:',err

      dwkernel=1
      call DOLS(R,DR1,u,.false.,baremass,SRF)
      call DOLop(R,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOLS-DOL:',err

      return

      end subroutine equivDOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine equivDOL_DDW()
      use gammas
      use overlapmoduledev
      use dwcoeffs
      use zolomodule
      use domainwallmod
      use rvmodule
      use gaugefield
      use kernelspectrarange
      implicit none
      complex(prc),dimension(Nv,4) :: RR,DR1,DR2,TMP
      real(prc) :: err
      type(sgnratfunc) :: SRF
      type(zolotarev) :: zolo
      integer i1,i2
      real(prc) lmin,lmax
      complex(prc) :: R5(Nv,4,Ls)
      complex(prc) :: DR5(Nv,4,Ls)
      complex(prc) :: TMP5(Nv,4,Ls)

      call getMinMaxHEigs(lmin,lmax)
      call setZoloCoeffs(Ls,SRF,lmin,lmax)
      call setZolo(lmin,lmax,Ls,zolo)
      call getRoots(zolo)
      omega=one/zolo%roots
      
      print *,zolo%roots
      call setRVs(Nv*4,RR)
!      call setHTcoeffs(Ls,SRF)

      baremass=0.05
      print *,"omega:",omega

      MTYPE=1
      dwkernel=3
      call KDDW4(RR,DR1,u,.false.,baremass)
      dwkernel=2
      call DOLop(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOL-KDDW4 ',err

      MTYPE=1
      dwkernel=3
      call IKDDW4(RR,DR1,u,.false.,baremass)
      dwkernel=2
      call IDOLop(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'IDOL-IKDDW4 ',err

      stop

      print *,""
      print *,"TEST EQUIVALENCE OF OVERLAP AND DOMAIN WALL FORMULATIONS"
      print *,""

      MTYPE=1
      do i2=1,2
        dwkernel=i2 ! 1:Shamir 2:Wilson
        print *,""
        print *,"MTYPE:",MTYPE,"DWKERNEL:",dwkernel
        print *,""
        call KDDW4(RR,DR1,u,.false.,baremass)
        call DOLop(RR,DR2,u,.false.,baremass,SRF)
        err=maxval(abs(DR1-DR2))
        print *,'DOL-KDDW4:',err
        call KDDW4(RR,DR1,u,.true.,baremass)
        call DOLop(RR,DR2,u,.true.,baremass,SRF)
        err=maxval(abs(DR1-DR2))
        print *,'DOLdag-KDDW4dag ',err
      end do

      return

      dwkernel=1
      MTYPE=3
      call KDDW4(RR,DR1,u,.false.,baremass)
      MTYPE=4
      call DOLop(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOL-KDDW4 M4/M3 Shamir:',err
      MTYPE=3
      call KDDW4(RR,DR1,u,.true.,baremass)
      MTYPE=4
      call DOLop(RR,DR2,u,.true.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOLdag-KDDW4dag M4/M3 Shamir:',err

      dwkernel=2
      MTYPE=3
      call KDDW4(RR,DR1,u,.false.,baremass)
      MTYPE=4
      call DOLop(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOL-KDDW4 M4/M3 Wilson:',err
      MTYPE=3
      call KDDW4(RR,DR1,u,.true.,baremass)
      MTYPE=4
      call DOLop(RR,DR2,u,.true.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOLdag-KDDW4dag M4/M3 Wilson:',err

      return

ccc   move the following to gamhermtest

!     test g5.KDDW4.g5-KDDW4dag = 0 
      DR1=RR
      call mGmu(DR1,4)
      call KDDW4(DR1,DR2,u,.false.,baremass)
      call mGmu(DR2,4)
      call KDDW4(RR,DR1,u,.true.,baremass)
      err=maxval(abs(DR1-DR2))
      print *,'err(KDDW4 g3 hermiticity)',err 

!     test g3.DOL.g3-DOLdag = 0 
      DR1=RR
      call mGmu(DR1,4)
      call Doverlap(DR1,DR2,u,.false.,baremass,SRF)
      call mGmu(DR2,4)
      call Doverlap(RR,DR1,u,.true.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'err(DOL g3 hermiticity)',err 

!     test g5.DOL.g5-DOLdag = 0 
      DR1=RR
      call mGmu(DR1,5)
      call Doverlap(DR1,DR2,u,.false.,baremass,SRF)
      call mGmu(DR2,5)
      call Doverlap(RR,DR1,u,.true.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'err(DOL g5 hermiticity)',err 


      return

      if(VERBOSE.eq.1)then ; print *,'DOL-KDDW4-Shamir-M1-HT' ; endif
      DWkernel=1
      MTYPE=1
      call KDDW4(RR,DR1,u,.false.,baremass)
      call DOLS(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'err(!DAGGER) ',err
      call KDDW4(RR,DR1,u,.true.,baremass)
      call DOLS(RR,DR2,u,.true.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'err(DAGGER) ',err

      return
      end subroutine equivDOL_DDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine outputKDDW4()
      use gammas
      use overlapmoduledev
      use dwcoeffs
      use zolomodule
      use domainwallmod
      use rvmodule
      use gaugefield
      use condnoisytoolsmod
      implicit none
      complex(prc),dimension(Nv,4) :: RR,DR1,DR2,TMP,DR
      real(prc) :: err
      type(sgnratfunc) :: SRF
      type(zolotarev) :: zolo
      integer i1,i2,j,l,d,Nnoise,Ninstance,i
      real(prc) lmin,lmax
      complex(prc) :: R5(Nv,4,Ls)
      complex(prc) :: DR5(Nv,4,Ls)
      complex(prc) :: TMP5(Nv,4,Ls)
      complex(prc) pbp,pbptot,pbpi

      print *,"Test Output KDDW4"

      lmin=1d-1  
      lmax=5.0  
      call setZoloCoeffs(Ls,SRF,lmin,lmax)
      call setZolo(lmin,lmax,Ls,zolo)
      call getRoots(zolo)
      omega=one/zolo%roots
      RR=one
      baremass=0.01
      print *,"omega:",omega
      MTYPE=3
      dwkernel=3
      R5=czero
      R5(:,:,1)=RR

      print *,""
      print *,"Tests:"
      print *,""

      if (.false.) then
        call KDDW(R5,DR5,u,.false.,baremass)
        print *,"T1:"
        do l=1,Ls
         print *,"l:",l
          do d=1,4
           print *,"d:",d
            do j=1,Nv
              print *,DR5(j,d,l)
            end do
          end do
        end do 
      end if

      if (.false.) then
        call KDDW4(RR,DR,u,.false.,baremass)
        print *,"KDDW4:"
          do d=1,4
           print *,"d:",d
            do j=1,Nv
              print *,DR(j,d)
            end do
          end do
      end if

      if (.false.) then
        call IKDDW4(RR,DR,u,.false.,baremass)
        print *,"IKDDW4:"
          do d=1,4
           print *,"d:",d
            do j=1,Nv
              print *,DR(j,d)
            end do
          end do
      end if

      if (.false.) then ! single noisy calculation
        pbptot=0
        Nnoise=10
          call setGRVs(3*Nv,theta)
          call coef(u,theta)
          do j=1,Nnoise
            dwkernel=3
            call evalCondNoisy_KDDW4(u,pbp)
            pbptot=pbptot+pbp
!          dwkernel=2
!          call evalCondNoisy_OL(u,pbp,SRF)
        end do
        pbptot=pbptot/Nnoise
        open(unit=10,file="instanceKDDW4.dat",status='unknown',
     &                                             access='append')
        print *,"pbp instance:",pbp
        write(10,*) "pbp instance:",pbp
        close(10)
      end if


      if (.true.) then ! measurement
        open(unit=11,file="condKDDW4.dat",status='unknown',
     &                                             access='append')
        pbptot=0
        Ninstance=10
        Nnoise=10
        do i=1,Ninstance
          call setGRVs(3*Nv,theta)
          theta=theta/3.0
          call coef(u,theta)
          pbpi=0
          do j=1,Nnoise
            dwkernel=3
            call evalCondNoisy_KDDW4(u,pbp)
            pbpi=pbpi+pbp
          end do
          pbpi=pbpi/Nnoise
          pbptot=pbptot+pbpi
          print *,"pbp instance:",pbpi
          write(11,*) "pbp instance:",pbpi
        end do
        pbptot=pbptot/Ninstance
        print *,"pbp average:",pbptot
        write(11,*) "pbp average:",pbptot
        close(11)
      end if

      goto 910
!      call DDW_OWilson(R5,TMP5,u,.false.,baremass)
      dwkernel=3
      call IDDW(R5,TMP5,u,.false.,baremass)
      call DDW_OWilson(TMP5,DR5,u,.false.,baremass)
      err=maxval(abs(R5-DR5))
      print *,"err",err


      stop

        call PermM(R5,TMP5,.false.,4)
        call DDW_OWilson(TMP5,DR5,u,.false.,baremass)
        MTYPE=1
        call IDDW(DR5,TMP5,u,.false.,one)
        call PermM(TMP5,DR5,.true.,4)
      TMP=DR5(:,:,1)
910   continue

      return
      end subroutine outputKDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testDequivmod
