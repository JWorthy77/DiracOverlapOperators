      module overlapspectrarange
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use kernelspectrarange
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcOverlapRange(fnum)
      use gaugefield
      use gaugemodule
      use options
      use IOmodule
      use statsmod
      use countmod
      implicit none
      integer fnum
      integer icf,idx
      character(len=80) fname 
      character(len=5) fnumc
      real(prc) lmax,lmin,lmaxav,lmaxsd,lminav,lminsd,cnav,cnsd
      real(prc),allocatable,dimension(:) :: vlmax,vlmin,vcn
      integer fstart,fstop,fskip,Nf
      logical evalMax,evalMin
      integer Nmax,Nmin,MAXMIN
      logical FEXISTS,ZOLO
      integer Nsrf
      real(prc) zmin,zmax
      type(sgnratfunc) SRF

      fname="overlapopts.txt"
      if(.not.QUENCHED) then
        fnumc=itoa(fnum)
        print *,fnum,fnumc
        fname="overlapopts.txt"//trim(fnumc)
        print *,fname
      endif
      open(unit=31,file=fname,status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) fstart,fstop,fskip,evalMax,Nmax,evalMin,Nmin,dwkernel
     &,GAUGETYPE,baremass,Nsrf,ZOLO,zmin,zmax,gbeta
        close(31)
        print *,"calcOverlapRange"
        print *,"fstart,fstop,fskip:",fstart,fstop,fskip
        print *,"evalMax,evalMin:",evalMax,evalMin
        print *,"Nmax,Nmin:",Nmax,Nmin
        print *,"kernel:",dwkernel
        print *,"gaugetype:",GAUGETYPE
        print *,"baremass:",baremass
        print *,"beta:",gbeta
        print *,"quenched:",QUENCHED
        if (ZOLO) then
          print *,"Zolotarev",zmin,zmax
        else
          print *,"HT"
        end if
        print *,"Ls:",Nsrf
      else
        print *,"kernel options file not found"
        stop
      endif
      close(31)

      if (ZOLO) then
        call setZoloCoeffs(Nsrf,SRF,zmin,zmax)
      else
        call setHTcoeffs(Nsrf,SRF)
      end if

      lmax=0
      lmin=0
      Nf=(fstop-fstart)/fskip+1
      allocate(vlmax(Nf),vlmin(Nf),vcn(Nf))
      idx=1
      do icf=fstart,fstop,fskip
        if (QUENCHED) then
          if (GAUGETYPE.eq.1) then
            call makeQuenchedCosineThirringField()
          elseif (GAUGETYPE.eq.2) then
            call makeQuenchedGaussianThirringField()
          endif
        else
          call readThetaFromFile(icf,theta)
          call coef(u,theta)
        endif
  
        if (evalMax) then
          MAXMIN=1
          call estimateOverlapExtrema(Nmax,MAXMIN,lmax,SRF)
        endif

        if (evalMin) then
          MAXMIN=-1
          call estimateOverlapExtrema(Nmin,MAXMIN,lmin,SRF)
        endif

        vlmax(idx)=lmax
        vlmin(idx)=lmin
        vcn(idx)=lmax/lmin
        idx=idx+1

        open(unit=14,file='Extrema.dat',access='append',
     &           status='unknown',form='formatted')
        write(14,'(I4,F6.2,2E12.4)') icf,MDW,lmax,lmin
        close(14)
      end do

      call calcVarReal(Nf,vlmax,lmaxav,lmaxsd) 
      call calcVarReal(Nf,vlmin,lminav,lminsd) 
      call calcVarReal(Nf,vcn,cnav,cnsd) 

      if (.false.) then
     
      open(unit=11,file='OverlapExtrema.dat',access='append',
     &           status='unknown',form='formatted')
        write(11,*) "Nterms,max,min,sd,err"
        write(11,*) "Max:"
        write(11,*) Nf,lmaxav,maxval(vlmax),minval(vlmax),lmaxsd,lmaxsd/
     &                                                 sqrt(real(Nf))
        write(11,*) "Min:"
        write(11,*) Nf,lminav,maxval(vlmin),minval(vlmin),lminsd,lminsd/
     &                                                 sqrt(real(Nf))
        write(11,*) "Condition number:"
        write(11,*) Nf,cnav,maxval(vcn),minval(vcn),cnsd,cnsd/
     &                                                 sqrt(real(Nf))
      close(11)

      else

      open(unit=11,file='DOLinvCount.dat',access='append',
     &           status='unknown',form='formatted')
      print *,"outer count av:",oc_idx,real(outer_count)/real(oc_idx)
      print *,"inner count av:",ic_idx,real(inner_count)/real(ic_idx)
      write(11,'(I10,F12.2,I10,F12.2)') oc_idx,
     &    real(outer_count)/real(oc_idx),ic_idx,
     &    real(inner_count)/real(ic_idx)
      close(11)
      end if

      return
      end subroutine calcOverlapRange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine overlaprangetask(fid,Nmax,emax,Nmin,emin,SRF)
      use gaugefield
      use gaugemodule
      use IOmodule
      use rrparams
      use ritzmod
      use rrspectrum
      use options
      implicit none
      integer fid,Nmin,Nmax
      real(prc) emax,emin
      type(sgnratfunc) :: SRF
      real(prc) lambda(nev)
  
      print *,"overlaprangetask,fid:",fid
      if(nev.ne.1)then ; print *,"no evs not 1" ; stop ; endif
      if (QUENCHED) then
        if (GAUGETYPE.eq.1) then
          call makeQuenchedCosineThirringField()
        elseif (GAUGETYPE.eq.2) then
          call makeQuenchedGaussianThirringField()
        endif
      else
        call readThetaFromFile(fid,theta)
        call coef(u,theta)
      endif
      if (Nmax > 0) call estimateOverlapExtrema(Nmax,1,emax,SRF)
c      if (Nmin > 0) call estimateOverlapExtrema(Nmin,-1,emin,SRF)
      KERNEL=.false. ! this is an rrparams option
      if (Nmin > 0) then
        call spectrum(-MDW,SRF,lambda)
        emin=lambda(1)
      endif

      return
      end subroutine overlaprangetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine estimateOverlapExtrema(Nmax,MAXMIN,eig,SRF)
      use options
      use rvmodule
      use gaugefield
      use gammas
      use overlapmoduledev
      use overlapmoduleprod
      use domainwallmod
      implicit none
      integer Nmax,MAXMIN
      real(prc) eig
      type(sgnratfunc) SRF
      complex(prc),dimension(Nv,Ndc) :: R,TMP,DR
      integer i

      call setRVs(Nv*4,R)
      call normalise(R)
      open(unit=81,file='DetailsEigs.dat',access='append',
     &           status='unknown',form='formatted')
      do i=1,Nmax
        if (MAXMIN.eq.1) then
          if (dwkernel.eq.1) then
            call KDDW4(R,DR,u,.false.,baremass)
            call mGmu(DR,4)
          elseif(dwkernel.eq.2) then
            call DOLop3(R,DR,u,.false.,baremass,SRF)
            call mGmu(DR,4)
          else
            print *,"only dwkernel 1 and 2 options implemented"
            stop
          endif
        elseif (MAXMIN.eq.-1) then
          if (dwkernel.eq.1) then
            call IKDDW4(R,DR,u,.false.,baremass)
            call mGmu(DR,4)
          elseif(dwkernel.eq.2) then
            call IDOLop3(R,DR,u,.false.,baremass,SRF)
            call mGmu(DR,4)
          else
            print *,"only dwkernel 1 and 2 options implemented"
            stop
          endif
        endif
        R=DR
        eig=mag(R)
        R=R/eig
        print *,i,eig
        write(81,*) i,eig
      end do
      close(81)
      if (MAXMIN.eq.-1) then
        eig=one/eig
      endif

      return
      end subroutine estimateOverlapExtrema
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module overlapspectrarange
