      module localitymodule
!     contains routines to test the locality of the overlap operator
      use arraysizes
      use numbers
      use options
      use rvmodule
      use ratfuncs
      use arpackmodule
      use gaugefield
      use overlapmoduledev
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testLocality(SWITCH,FILETYPE,Ng) ! test routine
      use indices
      use gaugefield
      use gaugemodule
      use domainwallmod
      use IOmodule
      implicit none
      integer SWITCH ! 1=HT Wilson,2=HT Shamir,3=Zolo Wilson,4=Zolo Shamir,5=DW,6=I HT Wilson
      integer FILETYPE ! 1 = Free field, 2 = HMC constructed, 3 = from converted file, 4 = from binary con file
      integer Ng
      type(sgnratfunc) :: SRF
      complex(prc) psi(Nv,4),DOL(Nv,4)
      complex(prc) psiDW(Nv,4,Ls),DDWpsi(Nv,4,Ls)
      procedure(),pointer :: Sptr => NULL()
      integer j,jdx
      real ev(1)
      real(prc) lmin,lmax
      real(prc) f(2*Ns+Nt,Ng),avf(2*Ns+Nt),sdf(2*Ns+Nt),err(2*Ns+Nt)
      real(prc) div(2*Ns+Nt)
      real(prc) fn
      integer rmax,dr
      integer nr(2*Ns+Nt)
      integer ig
      character(len=80) :: fname

      print *,"testLocality"
      print *,"SWITCH:",SWITCH
      print *,"FILETYPE:",FILETYPE
      print *,"Ng:",Ng

!     psi=delta
      psi=(zero,zero)
      psi(1,2)=(one,zero)
      if (SWITCH.eq.1) then

        jdx=40
        open(unit=12,file='localHTWilson.dat',status='unknown',
     &                                               form='formatted')
        Sptr => DOverlap

      elseif (SWITCH.eq.2) then
        jdx=40
        open(unit=12,file='localHTShamir.dat',status='unknown',
     &                                               form='formatted')
        Sptr => DOLS
      
      elseif (SWITCH.eq.3) then
        jdx=40
        open(unit=12,file='localZoloWilson.dat',status='unknown',
     &                                               form='formatted')
        Sptr => DOverlap

      elseif (SWITCH.eq.4) then
        jdx=40
        open(unit=12,file='localZoloShamir.dat',status='unknown',
     &                                               form='formatted')
        Sptr => DOLS

      elseif (SWITCH.eq.5) then
        open(unit=11,file='localDomWallAll.dat',status='unknown',
     &                                               form='formatted')
        open(unit=12,file='localDomWall.dat',status='unknown',
     &                                               form='formatted')
        Sptr => KDDW
        psiDW=(zero,zero)
        psiDW(1,1,1)=(one,zero)

      end if


      f=zero
      do ig=1,Ng
        print *,'ig:',ig
        if (FILETYPE.eq.1) then
          call makeGaugeField(.true.)
        elseif (FILETYPE.eq.2) then
          call makeGaugeField(.false.)
        elseif (FILETYPE.eq.3) then
c          ThirringFileDir=
c     &    "/home/jude/2019/Thirring/Sunbird/converted/b_0.33/"
          call readConvertedThetaFileName(ig,fname)
          call readConvertedThirringGaugeField(fname,theta)
          call coef(u,theta)
        elseif (FILETYPE.eq.4) then
          call readThetaFileName(100,fname)
          call readMPIConFile(fname,theta)
          call coef(u,theta)
        endif

        if (SWITCH.ne.5) then
          if ((SWITCH.eq.1).or.(SWITCH.eq.2).or.(SWITCH.eq.6)) then
            call setHTcoeffs(jdx,SRF)
          elseif((SWITCH.eq.3).or.(SWITCH.eq.4)) then
!            call calceigs('SM',1,ev,2,Nv)
            print *,"in locality module - update calceigs"
            stop
            print *,'ev min:',ev
            lmin=ev(1)
!            call calceigs('LM',1,ev,2,Nv)
            print *,"in locality module - update calceigs"
            stop
            print *,'ev max:',ev
            lmax=ev(1)
            call setZoloCoeffs(jdx,SRF,lmin,lmax)
          endif
          call Sptr(psi,DOL,u,.false.,baremass,SRF)
        else
          call Sptr(psiDW,DDWpsi,u,.false.,baremass)
          DOL=DDWpsi(:,:,1)
        endif

        open(unit=21,file='DOLgd.dat',status='unknown',form='formatted')
        write(21,*) DOL
        close(21)

!       now look at locality
        print *,"calc locality"
        rmax=0
        nr=0
        do j=1,Nv
          dr=dist(j,1,1,1)
          nr(dr+1)=nr(dr+1)+1
          rmax=max(rmax,dr)
          fn=dot_product(DOL(j,:),DOL(j,:))**0.5
          f(dr+1,ig)=max(f(dr+1,ig),fn)
        end do
        print *,"rmax:",rmax
        do j=0,rmax
          print *,j,f(j+1,ig),nr(j+1)
          write(11,*) 95+5*ig,j,f(j+1,ig),nr(j+1)
        end do
        flush(11)

      end do
      
      do j=0,rmax
        avf(j+1)=sum(f(j+1,:))/Ng
        sdf(j+1)=(sum((f(j+1,:)-avf(j+1))**2)/(Ng-1))**0.5
        err(j+1)=sdf(j+1)/real(Ng)**0.5
        if (j.gt.0) then
          div(j+1)=avf(j+1)/avf(j)
          write(12,*) j,avf(j+1),sdf(j+1),err(j+1),nr(j+1),div(j+1)
        else
          div(j+1)=0
          write(12,*) j,avf(j+1),sdf(j+1),err(j+1),nr(j+1)
        end if
        print *,j,avf(j+1),sdf(j+1),err(j+1),nr(j+1),div(j+1)
      end do

      close(11)
      close(12)

      return
      end subroutine testLocality
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalLocality(DOL,ix,iy,it,Ng,ig,rmax,nr,cnt,f)
      use indices
      implicit none
      complex(prc) DOL(Nv,4)
      integer Ng,ig,rmax,cnt
      integer nr(2*Ns+Nt)
      real(prc) f(2*Ns+Nt,64*Ng) ! note 64*Ng is dependent on outer loop
      integer ix,iy,it
      real(prc) fn
      integer j,dr

      print *,"calc locality"
      rmax=0 ; nr=0
      do j=1,Nv
        dr=dist(j,ix,iy,it)
        if (dr < 2) then
          print *,j,dr
        endif
        nr(dr+1)=nr(dr+1)+1
        rmax=max(rmax,dr)
        fn=dot_product(DOL(j,:),DOL(j,:))**0.5
        f(dr+1,cnt)=max(f(dr+1,cnt),fn)
      end do
      print *,rmax
      do j=0,rmax
        print *,j,f(j+1,cnt),nr(j+1)
        write(11,*) cnt,ig,ix,iy,it,j,f(j+1,cnt),nr(j+1)
      end do
      flush(11)
      cnt=cnt+1
      
      return
      end subroutine evalLocality
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalLocalityStats(rmax,Ng,f,nr,outname)
      implicit none
      integer rmax,Ng
      real(prc) f(2*Ns+Nt,Ng)
      integer nr(2*Ns+Nt)
      character(len=80) :: outname
      real(prc) avf(2*Ns+Nt),sdf(2*Ns+Nt),err(2*Ns+Nt)
      real(prc) div(2*Ns+Nt)
      integer j

      open(unit=12,file=outname,status='unknown',form='formatted')

      do j=0,rmax
        avf(j+1)=sum(f(j+1,:))/Ng
        sdf(j+1)=(sum((f(j+1,:)-avf(j+1))**2)/(Ng-1))**0.5
        err(j+1)=sdf(j+1)/real(Ng)**0.5
        if (j.gt.0) then
          div(j+1)=avf(j+1)/avf(j)
          write(12,*) j,avf(j+1),sdf(j+1),err(j+1),nr(j+1),div(j+1)
        else
          div(j+1)=0
          write(12,*) j,avf(j+1),sdf(j+1),err(j+1),nr(j+1)
        end if
        print *,j,avf(j+1),sdf(j+1),err(j+1),nr(j+1),div(j+1)
      end do

      close(12)

      return
      end subroutine evalLocalityStats
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testConvertedLocality(SWITCH,Ng)
      use indices
      use gaugefield
      use gaugemodule
      use domainwallmod
      use IOmodule
      implicit none
      integer SWITCH ! 1=HT Wilson,2=HT Shamir, 5=DW
      integer Ng
      type(sgnratfunc) :: SRF
      complex(prc) psi(Nv,4),DOL(Nv,4)
      complex(prc) psiDW(Nv,4,Ls),DDWpsi(Nv,4,Ls)
      procedure(),pointer :: Sptr => NULL()
      integer j,jdx
      real ev(1)
      real(prc) lmin,lmax
      real(prc) f(2*Ns+Nt,64*Ng),avf(2*Ns+Nt),sdf(2*Ns+Nt),err(2*Ns+Nt)
      real(prc) div(2*Ns+Nt)
      real(prc) fn
      integer rmax,dr
      integer nr(2*Ns+Nt)
      integer ig,ix,iy,it,idx,cnt
      character(len=80) :: fname

      if (SWITCH.eq.1) then

        jdx=48
        open(unit=12,file='localHTWilson.dat',status='unknown',
     &                                               form='formatted')
        Sptr => DOverlap
        call setHTcoeffs(jdx,SRF)

      elseif (SWITCH.eq.2) then

        jdx=40
        open(unit=12,file='localHTShamir.dat',status='unknown',
     &                                               form='formatted')
        Sptr => DOLS
        call setHTcoeffs(jdx,SRF)
      
      elseif (SWITCH.eq.5) then

        open(unit=11,file='localDomWallAll.dat',status='unknown',
     &                                               form='formatted')
        open(unit=12,file='localDomWall.dat',status='unknown',
     &                                               form='formatted')
        OLTYPE=2  ! 1=direct, 2=indirect
        DWkernel=1 ! 1=Shamir, 2=Wilson
        MTYPE=3 
        Sptr => KDDW
        GAUGETYPE=2
        baremass=0.005
        
      end if
c      ThirringFileDir=
c     &     "/home/jude/2019/Thirring/Sunbird/converted/b_0.33/"


      f=zero
      cnt=1
      do ig=1,12
        print *,'ig:',ig
        call readConvertedThetaFileName(ig,fname)
        call readConvertedThirringGaugeField(fname,theta)
        call coef(u,theta)
      
        do ix=1,Ns,8
          do iy=1,Ns,8
            do it=1,Nt,8

              call ia(ix,iy,it,idx)
!             psi=delta
              psi=(zero,zero)
              psi(idx,1)=(one,zero)

              print *,"calc D.psi",ix,iy,it,idx

              if (SWITCH.ne.5) then
                call Sptr(psi,DOL,u,.false.,baremass,SRF)
              else
                call KDDW4(psi,DOL,u,.false.,baremass)
              endif


!             now look at locality
              print *,"calc locality"
              rmax=0
              nr=0
              do j=1,Nv
                dr=dist(j,ix,iy,it)
                if (dr < 2) then
                  print *,j,dr
                endif
                nr(dr+1)=nr(dr+1)+1
                rmax=max(rmax,dr)
                fn=dot_product(DOL(j,:),DOL(j,:))**0.5
                f(dr+1,cnt)=max(f(dr+1,cnt),fn)
              end do
              print *,rmax
              do j=0,rmax
                print *,j,f(j+1,cnt),nr(j+1)
                write(11,*) cnt,ig,ix,iy,it,j,f(j+1,cnt),nr(j+1)
              end do
              flush(11)
              cnt=cnt+1
            enddo
          enddo
        enddo

      end do
      
      do j=0,rmax
        avf(j+1)=sum(f(j+1,:))/cnt
        sdf(j+1)=(sum((f(j+1,:)-avf(j+1))**2)/(cnt-1))**0.5
        err(j+1)=sdf(j+1)/real(cnt)**0.5
        if (j.gt.0) then
          div(j+1)=avf(j+1)/avf(j)
          write(12,*) j,avf(j+1),sdf(j+1),err(j+1),nr(j+1),div(j+1)
        else
          div(j+1)=0
          write(12,*) j,avf(j+1),sdf(j+1),err(j+1),nr(j+1)
        end if
        print *,j,avf(j+1),sdf(j+1),err(j+1),nr(j+1),div(j+1)
      end do

      close(11)
      close(12)

      return
      end subroutine testConvertedLocality
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcLocality() ! PhD results routine
c      use indices
      use gaugefield
      use gaugemodule
      use domainwallmod
      use IOmodule
      implicit none
      integer,parameter :: Ng=11
      complex(prc) psi(Nv,4),DOL(Nv,4)
      real(prc) f(2*Ns+Nt,64*Ng) ! note 8*Ng is dependent on loop
      integer rmax,dr
      integer nr(2*Ns+Nt)
      integer ig,ix,iy,it,idx,cnt
      integer j,fid
      character(len=80) :: fname,outname,outnamehist

!     this routine is tailored for the 20x20 Sunbird aux fields 500-1000
!     locality is calculated using the domain wall formulation

      outname='localityDomWall.dat'
      outnamehist='localityDomWallHist.dat'
      open(unit=11,file=outnamehist,status='unknown',form='formatted')
      OLTYPE=2  ! 1=direct, 2=indirect
      DWkernel=1 ! 1=Shamir, 2=Wilson
      MTYPE=3
      MDW=one 

      f=zero ; cnt=1 
      do ig=1,Ng
        fid=500+(ig-1)*50 ! 500-1000
        call readThetaFileName(fid,fname)
        call readMPIConFile(fname,theta)
        call coef(u,theta)
      
        do ix=1,Ns,5
          do iy=1,Ns,5
            do it=1,Nt,5 ! there are 64 of these - used in f array size

              call ia(ix,iy,it,idx)
              psi=(zero,zero)
              psi(idx,1)=(one,zero)

              print *,"calc D.psi",ix,iy,it,idx

              call KDDW4(psi,DOL,u,.false.,baremass)
              call evalLocality(DOL,ix,iy,it,Ng,ig,rmax,nr,cnt,f) ! writes to unit 11

            enddo
          enddo
        enddo

      end do
      close(11)
      
      call evalLocalityStats(rmax,Ng,f,nr,outname)

      return
      end subroutine calcLocality
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcQuenchedLocality() ! PhD results routine 
      use gaugefield
      use gaugemodule
      use overlapmoduleprod
      use IOmodule
      implicit none
!      integer,parameter :: Ng=11
      complex(prc) psi(Nv,4),DOL(Nv,4)
      real(prc),allocatable,dimension(:,:) :: f ! f(2*Ns+Nt,64*Ng) ! note 8*Ng is dependent on loop
      integer rmax,dr
      integer nr(2*Ns+Nt),Ng
      integer ig,ix,iy,it,idx,cnt
      integer j,fid
      character(len=80) :: fname,outname,outnamehist
      logical FEXISTS
      type(sgnratfunc) :: SRF

!     this routine is tailored for 16x16 quenched cases
      MTYPE=3
      MDW=one 
      GAUGETYPE=2
 
      open(unit=31,file='quenchedlocality.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) dwkernel,baremass,Ng,gbeta
        close(31)
        print *,dwkernel,baremass,Ng,gbeta
      else
        print *,"quenched locality options file not found"
        stop
      endif
      close(31)
      
      allocate(f(2*Ns+Nt,64*Ng)) ! note 8*Ng is dependent on loop


      if (dwkernel.eq.1) then
        print *,"oh dear"
        stop
      elseif (dwkernel.eq.2) then
        call setZoloCoeffs(24,SRF,0.001d0,10d0)
      endif

      outname='locality.dat'
      outnamehist='localityHist.dat'
      open(unit=11,file=outnamehist,status='unknown',form='formatted')

      f=zero ; cnt=1 
      do ig=1,Ng
        if (GAUGETYPE.eq.1) then
          call makeQuenchedCosineThirringField()
        elseif (GAUGETYPE.eq.2) then
          call makeQuenchedGaussianThirringField()
        endif
        call coef(u,theta)
      
        do ix=1,Ns,Ns/2
          do iy=1,Ns,Ns/2
            do it=1,Nt,Nt/2 ! there are 8 of these - used in f array size

              call ia(ix,iy,it,idx)
              psi=(zero,zero)
              psi(idx,1)=(one,zero)

              print *,"calc D.psi",ix,iy,it,idx
!              call KDDW4(psi,DOL,u,.false.,baremass)
              call DOLop3(psi,DOL,u,.false.,baremass,SRF)
              call evalLocality(DOL,ix,iy,it,Ng,ig,rmax,nr,cnt,f) ! writes to unit 11

            enddo
          enddo
        enddo

      end do
      close(11)
      
      call evalLocalityStats(rmax,Ng,f,nr,outname)

      return
      end subroutine calcQuenchedLocality
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module localitymodule
