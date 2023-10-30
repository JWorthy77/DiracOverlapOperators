!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module parallelmod
#ifdef PARALLEL
      use basicparallelmod
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initialisetasks()
      use rvmod
      use iomodule
      implicit none
      integer lseed
      character(len=80) dfname 
      character(len=3) cnum
      integer i,nc

!     initialise problem 
      lseed=137+rank ! note this makes the noisy cases with nproc>2 not
                     ! guaranteed to be identical since the distribution 
                     ! of the files is not fixed
      call setrn(lseed)

      cnum=itoa(rank)
      nc=len(trim(cnum))
      do i=nc+1,3
        cnum="0"//cnum
      enddo
      dfname="DetailsNKspectraHT"//cnum//".dat"
      open(unit=88,file=dfname,access='append',
     &           status='unknown',form='formatted')

      return
      end subroutine initialisetasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initialiseolspectratasks(ntasks,fstart,fstop,fskip,
     &                                       Nmax,Nmin,Nsrf,SRF)
      use rvmod
      use iomodule
      use ratfuncs
      use options
      implicit none
      integer ntasks,fstart,fstop,fskip,Nmax,Nmin,Nsrf
      logical evalMax,evalMin
      type(sgnratfunc) SRF
      integer lseed
      character(len=80) dfname 
      character(len=3) cnum
      integer i,nc,Nterms
      logical FEXISTS,ZOLO
      real(prc) zmin,zmax
 
!     initialise problem
      lseed=137
      call setrn(lseed)

      cnum=itoa(rank)
      nc=len(trim(cnum))
      do i=nc+1,3
        cnum="0"//cnum
      enddo

      dfname="DetailsNKspectraHT"//cnum//".dat"
      open(unit=88,file=dfname,access='append',
     &           status='unknown',form='formatted')

      open(unit=31,file='overlapopts.txt',status='unknown')
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
        if (ZOLO) then
          print *,"Zolotarev",zmin,zmax
        else
          print *,"HT"
        end if
        print *,"Ls:",Nsrf
      else
        print *,"overlap options file not found"
        stop
      endif
      close(31)

      if (ZOLO) then
        call setZoloCoeffs(Nsrf,SRF,zmin,zmax)
      else
        call setHTcoeffs(Nsrf,SRF)
      end if

      ntasks=(fstop-fstart)/fskip+1

      return
      end subroutine initialiseolspectratasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initialiseoverlaprangetasks(ntasks,fstart,fstop,fskip,
     &                                       Nmax,Nmin,Nsrf,SRF)
      use rvmod
      use iomodule
      use ratfuncs
      use options
      implicit none
      integer ntasks,fstart,fstop,fskip,Nmax,Nmin,Nsrf
      logical evalMax,evalMin
      type(sgnratfunc) SRF
      integer lseed
      character(len=80) dfname 
      character(len=3) cnum
      integer i,nc,Nterms
      logical FEXISTS,ZOLO
      real(prc) zmin,zmax
 
!     initialise problem
      lseed=137
      call setrn(lseed)

      cnum=itoa(rank)
      nc=len(trim(cnum))
      do i=nc+1,3
        cnum="0"//cnum
      enddo

      dfname="overlapDetails"//cnum//".dat"
      open(unit=88,file=dfname,access='append',
     &           status='unknown',form='formatted')

      open(unit=31,file='overlapopts.txt',status='unknown')
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
        if (ZOLO) then
          print *,"Zolotarev",zmin,zmax
        else
          print *,"HT"
        end if
        print *,"Ls:",Nsrf
      else
        print *,"overlap options file not found"
        stop
      endif
      close(31)

      if (ZOLO) then
        call setZoloCoeffs(Nsrf,SRF,zmin,zmax)
      else
        call setHTcoeffs(Nsrf,SRF)
      end if

      ntasks=(fstop-fstart)/fskip+1

      return
      end subroutine initialiseoverlaprangetasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initialisecondensatetasks(ntasks,fstart,fstop,fskip,
     &                                                 nnoise,SRF)
      use rvmod
      use iomodule
      use ratfuncs
      use options
      implicit none
      integer ntasks,fstart,fstop,fskip,nnoise
      type(sgnratfunc) SRF
      integer lseed
      character(len=80) dfname 
      character(len=3) cnum
      integer i,nc,Nterms
      logical FEXISTS,ZOLO
      real(prc) lmin,lmax
 

!     initialise problem
!      lseed=137
!      call setrn(lseed)

      cnum=itoa(rank)
      nc=len(trim(cnum))
      do i=nc+1,3
        cnum="0"//cnum
      enddo

      dfname="condDetails"//cnum//".dat"
      open(unit=88,file=dfname,access='append',
     &           status='unknown',form='formatted')


      open(unit=31,file='noisycondopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) fstart,fstop,fskip,nnoise
        close(31)
        if (rank .eq. 0) then
          print *,"opts1:",fstart,fstop,fskip,nnoise
        end if
      else
        print *,"Condensate options file not found"
        stop
      endif
      close(31)

      ntasks=(fstop-fstart)/fskip+1

      open(unit=31,file='olopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) dwkernel,MTYPE,baremass,Nterms,ZOLO,lmin,lmax
        close(31)
        if (rank .eq. 0) then
          print *,"opts2:",dwkernel,MTYPE,baremass,Nterms,ZOLO,lmin,lmax
        end if
      else
        print *,"Overlap options file not found"
        stop
      endif
      close(31)
      if (ZOLO) then
        call setZoloCoeffs(Nterms,SRF,lmin,lmax)
      else
        call setHTcoeffs(Nterms,SRF)
      endif

      return
      end subroutine initialisecondensatetasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initialisequenchedcondensatetasks(ntasks,nnoise,
     &                               strength,SRF)
      use rvmod
      use iomodule
      use ratfuncs
      use options
      implicit none
      integer ntasks,nnoise
      real(prc) strength
      type(sgnratfunc) SRF
      integer lseed
      character(len=80) dfname 
      character(len=3) cnum
      integer i,nc,Nterms
      logical FEXISTS,ZOLO
      real(prc) lmin,lmax
 

!     initialise problem
      lseed=137
      call setrn(lseed)

      cnum=itoa(rank)
      nc=len(trim(cnum))
      do i=nc+1,3
        cnum="0"//cnum
      enddo

c      dfname="condDetails"//cnum//".dat"
c      open(unit=88,file=dfname,access='append',
c     &           status='unknown',form='formatted')


      open(unit=31,file='noisyqcondopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) ntasks,nnoise,strength
        close(31)
        if (rank .eq. 0) then
          print *,"opts1:",ntasks,nnoise,strength
        end if
      else
        print *,"Condensate options file not found"
        stop
      endif
      close(31)

      open(unit=31,file='olopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) dwkernel,MTYPE,baremass,Nterms,ZOLO,lmin,lmax
        close(31)
        if (rank .eq. 0) then
          print *,"opts2:",dwkernel,MTYPE,baremass,Nterms,ZOLO,lmin,lmax
        end if
      else
        print *,"Overlap options file not found"
        stop
      endif
      close(31)
      if (ZOLO) then
        call setZoloCoeffs(Nterms,SRF,lmin,lmax)
      else
        call setHTcoeffs(Nterms,SRF)
      endif

      return
      end subroutine initialisequenchedcondensatetasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine finalisecondensatetasks(ntasks,pbps)
      use options
      use statsmod
      implicit none
      integer ntasks
      complex(prc),dimension(:) :: pbps
      complex(prc) pbp,sd
      character(len=80) scmd 

      call calcVar(ntasks,pbps,pbp,sd)
      open(unit=24,file='out.dat',status='unknown')
      write(24,*) "Condensate measurements"
      write(24,*) "pbp:",real(pbp),pbp
      write(24,*) "sd:",real(sd),sd
      write(24,*) "err:",real(sd)/sqrt(real(ntasks))
      write(24,*) "kernel:",dwkernel
      write(24,*) "MTYPE:",MTYPE
      write(24,*) "nfields:",ntasks
      close(24)     

c      scmd="rm confiledir.txt"
c      call system(scmd)
c      scmd="rm machine*"
c      call system(scmd)
c      scmd="rm olopts.txt"
c      call system(scmd)
c      scmd="rm noisycondopts.txt"
c      call system(scmd)
c      scmd="rm out.dat"
c      call system(scmd)
c      scmd="rm condDetails000.dat"
c      call system(scmd)
c      scmd="rm condDetails001.dat"
c      call system(scmd)
c      scmd="rm out.dat"
c      call system(scmd)
c      scmd="rm scr.run"
c      call system(scmd)

      return
      end subroutine finalisecondensatetasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine finaliseoverlaprangetasks(ntasks,emaxs,emins)
      use options
      use statsmod
      implicit none
      integer ntasks
      real(prc),dimension(:) :: emaxs,emins
      real(prc) emax,emin,sdmax,sdmin
      real(prc) denom
      integer i
      real(prc) Lmin,Lmax,Lmins(ntasks),Lmaxs(ntasks),Lsdmax,Lsdmin

      call calcVarReal(ntasks,emaxs,emax,sdmax)
      call calcVarReal(ntasks,emins,emin,sdmin)
      open(unit=24,file='OLrange.dat',status='unknown')
      write(24,*) "Overlap Range"
      write(24,*) "min/max:",emin,emax
      write(24,*) "sd:",sdmin,sdmax
      denom=sqrt(real(ntasks))
      write(24,*) "err:",sdmin/denom,sdmax/denom
      write(24,*) "kernel:",dwkernel
      write(24,*) "MTYPE:",MTYPE
      write(24,*) "nfields:",ntasks
      do i=1,ntasks
        Lmins(i)=emins(i)/sqrt(1-emins(i)*emins(i))
        Lmaxs(i)=emaxs(i)/sqrt(1-emaxs(i)*emaxs(i))
      enddo
      call calcVarReal(ntasks,Lmaxs,Lmax,Lsdmax)
      call calcVarReal(ntasks,Lmins,Lmin,Lsdmin)
      write(24,*) "Lmin/Lmax:",Lmin,Lmax,sqrt(Lmin),sqrt(Lmax)
      write(24,*) "L*Lmin/L*Lmax:",Ns*Lmin,Ns*Lmax,Ns*sqrt(Lmin),
     &                             Ns*sqrt(Lmax)
      close(24)     

      return
      end subroutine finaliseoverlaprangetasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dotask(taskid,lambda)
      use spectra
      use rrparams
      implicit none
      integer taskid
      real*8 lambda(nev)
      integer fid

      fid=500+10*taskid
c      print *,'rank',rank,'doing task',taskid
c      call spectratask(fid,lambda)     

      return
      end subroutine dotask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dooverlapspectratask(taskid,lambda,fstart,fskip,SRF)
      use rrparams
      use spectra
      use ratfuncs
      implicit none
      integer taskid,ntasks,fstart,fskip
      real*8 lambda(nev)
      type(sgnratfunc) SRF
      integer fid

      fid=fstart+(taskid-1)*fskip
      print *,'rank',rank,'doing task',taskid,"fid:",fid
      call spectratask(fid,lambda,SRF)

      return
      end subroutine dooverlapspectratask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine docondensatetask(taskid,pbp,ntasks,fstart,fstop,fskip,
     &                                               nnoise,SRF)
      use condensatemodule
      implicit none
      integer ntasks,fstart,fstop,fskip,nnoise
      type(sgnratfunc) SRF
      integer taskid
      complex(prc) pbp
      integer fid

      fid=fstart+(taskid-1)*fskip
      print *,'rank',rank,'doing task',taskid
      call condensatetask(fid,pbp,nnoise,SRF)     

      return
      end subroutine docondensatetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine doquenchedcondensatetask(taskid,pbp,nnoise,
     &                                               strength,SRF)
      use condensatemodule
      implicit none
      integer nnoise
      real(prc) strength
      type(sgnratfunc) SRF
      integer taskid
      complex(prc) pbp
      integer fid

      print *,'rank',rank,'doing quenched condensate task',taskid
      call qcondensatetask(pbp,nnoise,strength,SRF)     
      print *,'rank',rank,'done quenched condensate task',taskid

      return
      end subroutine doquenchedcondensatetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine docondensateminitask(taskid,pbp,ntasks,
     &                                   fstart,fstop,fskip,nnoise,SRF)
      use condensatemodule
      implicit none
      integer ntasks,fstart,fstop,fskip,nnoise
      type(sgnratfunc) SRF
      integer taskid
      complex(prc) pbp
      integer fid,overtaskid,didx

      print *,"docondensateminitask"

      overtaskid=taskid/(4*nnoise)
      if (mod(taskid,nnoise*4).eq.0) then
        overtaskid=overtaskid-1
      endif
      fid=fstart+overtaskid*fskip
      didx=mod(taskid,4)
      if (didx.eq.0) then
        didx=4
      end if
      print *,'rank',rank,'mini task',taskid,"fid:",fid,"didx:",didx
      call condensatesubtask(fid,pbp,didx,SRF)     
      print *,'rank',rank,'mini task done'

      return
      end subroutine docondensateminitask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dooverlapspectrarangetask(taskid,emax,emin,
     &    fstart,fskip,Nmax,Nmin,SRF)
      use overlapspectrarange
      implicit none
      integer taskid,ntasks,fstart,fskip,Nmax,Nmin
      real(prc) emax,emin
      type(sgnratfunc) SRF
      complex(prc) pbp
      integer fid

      fid=fstart+(taskid-1)*fskip
      print *,'rank',rank,'doing task',taskid,"fid:",fid
      call overlaprangetask(fid,Nmax,emax,Nmin,emin,SRF)

      return
      end subroutine dooverlapspectrarangetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writetask(taskid,lambda)
      use spectra
      use rrparams
      implicit none
      integer taskid
      real*8 lambda(nev)
      integer fid

      fid=500+10*taskid
      open(unit=11,file='NKspectraHT.dat',access='append',
     &           status='unknown',form='formatted')
      write(11,*) fid,lambda
      close(11)

      return
      end subroutine writetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writeolspectratask(taskid,lambda,fstart,fskip)
      use spectra
      use rrparams
      implicit none
      integer taskid,fstart,fskip
      real*8 lambda(nev)
      integer fid

      fid=fstart+fskip*taskid
      open(unit=11,file='NKspectraHT.dat',access='append',
     &           status='unknown',form='formatted')
      write(11,*) fid,lambda
      close(11)

      return
      end subroutine writeolspectratask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writeoverlaprangetask(taskid,emax,emin,fstart,fskip)
      use pacc
      implicit none
      integer taskid
      real(prc) emax,emin
      integer fstart,fskip,fid

      fid=fstart+(taskid-1)*fskip
      open(unit=11,file='OLSpectraRange.dat',access='append',
     &           status='unknown',form='formatted')
      write(11,*) fid,emax,emin
      close(11)

      return
      end subroutine writeoverlaprangetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writecondensatetask(taskid,pbp,fstart,fskip)
c      use spectra
c      use rrparams
      use iomodule
      implicit none
      integer taskid
      complex(prc) pbp
      integer fstart,fskip,fid,nc,i
      character(len=80) cfname 
      character(len=4) cnum

      fid=fstart+(taskid-1)*fskip
      cnum=itoa(fid)
      nc=len(trim(cnum))
      do i=nc+1,4
        cnum="0"//cnum
      enddo
      cfname="../cond"//cnum//".dat"

      open(unit=11,file=cfname,access='append',
     &           status='unknown',form='formatted')
      write(11,*) fid,real(pbp,prc),dimag(pbp)
      close(11)

      return
      end subroutine writecondensatetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writeqcondensatetask(taskid,pbp,nnoise)
      use options
      implicit none
      integer taskid,nnoise
      complex(prc) pbp

      open(unit=11,file='cond.dat',access='append',
     &           status='unknown',form='formatted')
      write(11,*) taskid,nnoise,real(pbp,prc),dimag(pbp)
      close(11)

      return
      end subroutine writeqcondensatetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writecondensateminitask(taskid,pbp,fstart,fskip,nnoise)
      use spectra
      use rrparams
      implicit none
      integer taskid
      complex(prc) pbp
      integer fstart,fskip,fid,nnoise
      integer overtask,didx,inoise
  
      overtask=taskid/(nnoise*4)
      if (mod(taskid,nnoise*4).eq.0) then
        overtask=overtask-1
      endif
      fid=fstart+overtask*fskip

      didx=mod(taskid,4)
      if (didx.eq.0) then
        didx=4
      endif

      overtask=mod(taskid,nnoise*4)      
      if (mod(taskid,nnoise*4).eq.0) then
        overtask=nnoise*4
      endif
      inoise=(overtask-1)/4+1

      open(unit=11,file='cond.dat',access='append',
     &           status='unknown',form='formatted')
      write(11,'(3I6,2E19.10)') fid,inoise,didx,real(pbp,prc),dimag(pbp)
      close(11)

      return
      end subroutine writecondensateminitask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine parallelSequence(ntasks,SWITCH)
      use mpi
      use pacc
      use rrparams
      use ratfuncs
      use options
      implicit none
      integer ntasks
      integer SWITCH ! 1 is spectra calcs, 2 is condensate
      integer i,p,taskid,next_proc,tag
      integer task(nprocs),itaskid
      logical INUSE(nprocs),OPEN_PROCESSES,UNASSIGNED_TASKS
      logical AVAILABLE_PROCS
      logical SENDINFO,RECEIVEINFO
      integer sourceproc
      integer sendloc,recvloc,ierr
c      real*8 senddouble,recvdouble
      logical sendbool,recvbool
      integer status(MPI_STATUS_SIZE),request
      integer closedproc(nprocs-1)
      real*8 lambda(nev),emm(2)
      complex(prc) pbp
      complex(prc),allocatable,dimension(:) :: pbps
      real(prc),allocatable,dimension(:) :: emaxs,emins
      real(prc) emax,emin,strength
      logical,parameter :: VB_PS=.false.
      integer fstart,fstop,fskip,nnoise,Nmax,Nmin,Nsrf
      type(sgnratfunc) SRF

!      print *,"parallelSequence"
!      print *,rank,"nprocs",nprocs
      if (nprocs.eq.1) then
        print *,"parallel sequence: only 1 processor - stopping"
        stop
      endif

      if (SWITCH.eq.1) then
c        call initialisetasks()
        call initialiseolspectratasks(ntasks,fstart,fstop,fskip,
     &                                       Nmax,Nmin,Nsrf,SRF)
      elseif (SWITCH.eq.2) then ! each task covers each aux field
        call initialisecondensatetasks(ntasks,fstart,fstop,fskip,
     &                                                 nnoise,SRF)
        allocate(pbps(ntasks))
      elseif (SWITCH.eq.3) then ! each task covers each noisy iteration
        call initialisecondensatetasks(ntasks,fstart,fstop,fskip,
     &                                                 nnoise,SRF)
        ntasks=nnoise*ntasks*4
        allocate(pbps(ntasks))
      elseif (SWITCH.eq.4) then
        call initialiseoverlaprangetasks(ntasks,fstart,fstop,fskip,
     &                                       Nmax,Nmin,Nsrf,SRF)
        allocate(emaxs(ntasks),emins(ntasks))
      elseif (SWITCH.eq.5) then ! quenched condensate - each task covers each aux field
        call initialisequenchedcondensatetasks(ntasks,nnoise,
     &                                               strength,SRF)
        allocate(pbps(ntasks))

      elseif (SWITCH.eq.0) then
        print *,rank,"TEST",ntasks
      else
        print *,"parallelSequence SWITCH not set"
        stop
      endif

      sendloc=-1
      recvloc=-1
      closedproc=0

      UNASSIGNED_TASKS=.true.
      if (rank.ne.0) then

        SENDINFO=.false.
        do while (UNASSIGNED_TASKS)

          ! send availablity to master
          sendloc=rank
          call MPI_SEND(sendloc,1,MPI_INTEGER,0,10,MPI_COMM_WORLD,ierr)
          if (VB_PS) print *,'rank:',rank,'sent rank availability:',
     &                                                          sendloc

          ! indicate if any information from last completed task to be sent
          sendbool=SENDINFO
          call MPI_SEND(sendbool,1,MPI_LOGICAL,0,9,MPI_COMM_WORLD,ierr)
          if (VB_PS) print *,'rank:',rank,'send info?',sendbool

          ! send task info
          if (SENDINFO) then
            itaskid=recvloc ! taskid
            call MPI_SEND(itaskid,1,MPI_INTEGER,0,8,MPI_COMM_WORLD,ierr)
c            print *,'rank',rank,'sending info for taskid:',sendloc
            if (SWITCH.eq.1) then
              call MPI_SEND(lambda,nev,MPI_DOUBLE_PRECISION,0,8,
     &                                              MPI_COMM_WORLD,ierr)
            elseif ((SWITCH.eq.2).or.(SWITCH.eq.3).or.(SWITCH.eq.5))then
              call MPI_SEND(pbp,1,MPI_DOUBLE_COMPLEX,0,8,
     &                                              MPI_COMM_WORLD,ierr)
            elseif (SWITCH.eq.4) then
              emm(1)=emax
              emm(2)=emin
              call MPI_SEND(emm,2,MPI_DOUBLE_PRECISION,0,8,
     &                                              MPI_COMM_WORLD,ierr)
            endif
          endif

          ! receive taskid
          call MPI_RECV(recvloc,1,MPI_INT,0,10+rank,
     &                                 MPI_COMM_WORLD,status,ierr)
          print *,'rank:',rank,'recieved taskid:',recvloc

          ! carry out task
          if (recvloc.eq.-1) then
            ! cleanup
            UNASSIGNED_TASKS=.false.
            print *,'rank',rank,'no more tasks'
          else
            ! carry out primary task
            print *,'rank',rank,'doing task',recvloc
            if (SWITCH.eq.1) then
c              call dotask(recvloc,lambda)
              call dooverlapspectratask(taskid,lambda,fstart,fskip,SRF)
            elseif (SWITCH.eq.2) then
              call docondensatetask(recvloc,pbp,ntasks,fstart,fstop,
     &                                               fskip,nnoise,SRF)
            elseif (SWITCH.eq.3) then
              call docondensateminitask(recvloc,pbp,ntasks,fstart,fstop,
     &                                               fskip,nnoise,SRF)
            elseif (SWITCH.eq.4) then

              call dooverlapspectrarangetask(recvloc,emax,emin,
     &                                  fstart,fskip,Nmax,Nmin,SRF)
            elseif (SWITCH.eq.5) then
              call doquenchedcondensatetask(recvloc,pbp,nnoise,
     &                                                strength,SRF)
            endif
            SENDINFO=.true. 
c            senddouble=dble(recvloc+0.5) ! info
          endif

        enddo

      elseif (rank.eq.0) then

        OPEN_PROCESSES=.true.
        taskid=1
        do while (OPEN_PROCESSES)
          print *,'master allocate task:',taskid

          ! get available process
          call MPI_RECV(recvloc,1,MPI_INT,MPI_ANY_SOURCE,10,
     &                                 MPI_COMM_WORLD,status,ierr)
          if (VB_PS) print *,'master avail proc:',recvloc

          ! is there information to collect?
          sourceproc=recvloc
          call MPI_RECV(recvbool,1,MPI_LOGICAL,sourceproc,9,
     &                                 MPI_COMM_WORLD,status,ierr)
          if (VB_PS) print *,'master collect info? from:',recvbool

          ! collect task information
          RECEIVEINFO=recvbool
          if (RECEIVEINFO) then
            call MPI_RECV(itaskid,1,MPI_INT,sourceproc,8,
     &                                 MPI_COMM_WORLD,status,ierr)
            if(VB_PS) then
              print *,'master receiving from:',
     &                             sourceproc,'taskid/fid:',recvloc
            end if
            if (SWITCH.eq.1) then
              call MPI_RECV(lambda,nev,MPI_DOUBLE,sourceproc,8,
     &                                 MPI_COMM_WORLD,status,ierr)
c              call writetask(itaskid,lambda)
              call writeolspectratask(itaskid,lambda,fstart,fskip)
            elseif (SWITCH.eq.2) then
              call MPI_RECV(pbp,1,MPI_DOUBLE_COMPLEX,sourceproc,8,
     &                                 MPI_COMM_WORLD,status,ierr)
              call writecondensatetask(itaskid,pbp,fstart,fskip)
              pbps(itaskid)=pbp
            elseif (SWITCH.eq.3) then
              call MPI_RECV(pbp,1,MPI_DOUBLE_COMPLEX,sourceproc,8,
     &                                 MPI_COMM_WORLD,status,ierr)
              call writecondensateminitask(itaskid,pbp,
     &                                    fstart,fskip,nnoise)
              pbps(itaskid)=pbp
            elseif (SWITCH.eq.4) then
              call MPI_RECV(emm,2,MPI_DOUBLE_PRECISION,sourceproc,8,
     &                                 MPI_COMM_WORLD,status,ierr)
              emax=emm(1)
              emin=emm(2)
              call writeoverlaprangetask(itaskid,emax,emin,fstart,fskip)
              emaxs(itaskid)=emax
              emins(itaskid)=emin
            elseif (SWITCH.eq.5) then
              call MPI_RECV(pbp,1,MPI_DOUBLE_COMPLEX,sourceproc,8,
     &                                 MPI_COMM_WORLD,status,ierr)
              call writeqcondensatetask(itaskid,pbp,nnoise)
              pbps(itaskid)=pbp
            endif
            print *,'rank',rank,'received from:',
     &                             sourceproc,'eigenvalues:',lambda
          endif


          ! assign task
          if (taskid.gt.ntasks) then
            sendloc=-1
            closedproc(sourceproc)=1
          else
            sendloc=taskid
            taskid=taskid+1
          endif
          call MPI_SEND(sendloc,1,MPI_INT,recvloc,10+recvloc,
     &                                 MPI_COMM_WORLD,ierr)
          if (VB_PS) print *,'master sent task:',sendloc,'to',recvloc
          if (sum(closedproc).eq.nprocs-1) then
            OPEN_PROCESSES=.false.
          endif
        end do

        if ((SWITCH.eq.2).or.(SWITCH.eq.3).or.(SWITCH.eq.5))then
          call finalisecondensatetasks(ntasks,pbps)
        else if (SWITCH.eq.4)then
          call finaliseoverlaprangetasks(ntasks,emaxs,emins)
        endif

      endif

      return
      end subroutine parallelSequence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
      end module parallelmod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
