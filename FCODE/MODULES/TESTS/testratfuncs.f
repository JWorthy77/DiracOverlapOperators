      module testratfuncsmod
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testRatFuncs
      implicit none

      call testZoloCoeffs()
!      call testRationalFunctions()
!      call checkRationalFunctionCoeffs()

      return
      end subroutine testRatFuncs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testZoloCoeffs
      use zolomodule
      use ratfuncs
      implicit none
      type(zolotarev) :: zolo
      type(sgnratfunc) :: SRF

      call setZolo(1d-3,10d0,21,zolo)
      call getRoots(zolo)
      call setZoloCoeffs(21,SRF,1d-3,10d0)
      
      print *,"roots: ",zolo%roots
      end subroutine testZoloCoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testRationalFunctions()
      use numbers
      use ratfuncs
      use zolomodule
      implicit none
      real(prc) lmin,lmax
      integer Nr
      type(zolotarev) :: zolo

c      call testDividePoly
c      call testPartFrac
c      call testPartialFraction()

      Nr=256
c      lmin=real(1.0,prc)
c      lmax=real(100.0,prc)
c      call setZolo(lmin,lmax,Nr,zolo)
c      call getRoots(zolo)
c      call testZolo(Nr)
      call testHTFunctions(Nr,one/1000,1000*one)
c      call testZoloFunctions(Nr,one/2,2*one)

c      call testHTFunctions(21)
c      call testZoloFunctions(36,one,1000*one)
c      call testZoloFunctions(27,one,1000*one)
c      call testZoloFunctions(21,one/100000,one/10)
c      call testZoloFunctions(22,one/100000,one/10)

      return
      end subroutine testRationalFunctions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine checkRationalFunctionCoeffs()
      use numbers
      use ratfuncs
      use zolomodule
      implicit none
      type(sgnratfunc) :: HT,ZOLO,HTodd,HTeven
      real(prc) lmin,lmax

      call setHTcoeffs(40,HT)
      call setZoloCoeffs(40,ZOLO,0.0001*one,0.1*one)
      print *,"HT:"
      call printSRatFunc(HT)
      print *,"ZOLO:"
      call printSRatFunc(ZOLO)
      
      call setHTcoeffs(256,HTeven)
      call setHTcoeffs(11,HTodd)
      print *,"HTeven:"
      call printSRatFunc(HTeven)
      print *,"HTodd:"
      call printSRatFunc(HTodd)
      return
      end subroutine checkRationalFunctionCoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testratfuncsmod
