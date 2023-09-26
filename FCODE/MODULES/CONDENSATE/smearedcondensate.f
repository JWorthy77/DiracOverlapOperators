      module smearedcondensate
      use pacc
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalSmearedCondensateDomWall()
!     calculate overlap condensate with direct (partial fraction) 
!     or domain wall formulation with either Wilson or Shamir kernel
      use rvmodule
      use gammas
      use options
      use ratfuncs
      use overlapmoduledev
      use gaugefield
      use domainwallmod
      use statsmod
      use indices
      use IOmodule
      implicit none
      integer,parameter :: Nl=4
      integer SWITCH
      complex(prc) pbp(Nl,Nl)
      complex(prc),dimension(Nv,4,Ls) :: delta,DR
      complex(prc),dimension(Nv,4,Ls) :: delta1,delta2,deltaN,deltaNm1
      integer j,idr,ix,iy,it
      integer li,lj,tidx,Nvol,Np
      complex(prc) :: mean,stdev
      real(prc) :: rval
      character*128 OutFile,mkcmd
      character(len=80) fname 

      print *,"Eval domain wall smeared condensate"
!     done in Simon's code
      return

c#ifndef THIRRINGCASEDIR
c#define THIRRINGCASEDIR '8x8/Ls20/Bp32/mp01/g3/'
c#endif
c      ThirringCase=THIRRINGCASEDIR
c      print *,'Case',ThirringCase
c      ThirringFileDir='/home/jude/2020/GaugeFields/'//ThirringCase
c      OutFile=trim('/home/jude/2020/RUNS/SMEARED/'//ThirringCase)
c      mkcmd=trim('mkdir -p '//Outfile)
c      call system(mkcmd)
c      Outfile=trim(Outfile)//'smearedShamirLs40.dat'

      Nvol=0
      pbp=czero
      do tidx=10,100,10
        call readConvertedThirringGaugeField(fname,theta)
        call coef(u,theta)
        print *,"eval condensate for gauge field",tidx

      SWITCH=1
      if (SWITCH.eq.1) then ! point method
!      do j=1,Nv,4*Ns
      Np=0
      do ix=1,Ns,Ns/2
      do iy=1,Ns,Ns/2
      do it=1,Nt,Nt/2
        call ia(ix,iy,it,j)
        Np=Np+1
        Nvol=Nvol+1
        do idr=1,4
            Np=Np+1
            print *,idr,j,Nv
            if (idr < 3) then
              do li=1,Nl
                delta=czero
                delta(j,idr,li)=cone
                call IDDW(delta,DR,u,.false.,baremass)
                if (MTYPE.eq.3) then
                  call mGmu5(DR,4)
                endif
                do lj=1,Nl
                  pbp(li,lj)=pbp(li,lj)+DR(j,idr,Ls+1-lj)
                end do
              end do
            elseif (idr > 2) then
              do li=1,Nl
                delta=czero
                delta(j,idr,Ls+1-li)=cone
                call IDDW(delta,DR,u,.false.,baremass)
                if (MTYPE.eq.3) then
                  call mGmu5(DR,4)
                endif
                do lj=1,Nl
                  pbp(li,lj)=pbp(li,lj)+DR(j,idr,lj)
                end do
              end do
            endif

        end do
      end do
      end do
      end do
      
      elseif (SWITCH.eq.2) then ! noisy estimator
        print *,"noisy estimator not implemented"
        stop
      endif

      end do ! Thirring file loop

      if (MTYPE .eq. 3) then
        do li=1,Nl
          do lj=1,Nl
            pbp(lj,li)=-zi*pbp(lj,li)
          end do
        end do
      end if
 
c      call calcVar(Np,pbpVec(1:Np),mean,stdev)
c      print *,"pbp",pbp/Np,"sd",stdev,"err",(stdev*stdev/Np)**half
c      pbp=pbp/Np
c      open(unit=12,file='smearedShamir.dat',status='unknown',
c     &     access='append',form='formatted')
c      write(12,*) Ns,Ls,pbp1,pbp12,pbp21,pbp2
c      close(12)

      open(unit=13,file=Outfile,status='unknown',
     &     access='append',form='formatted')

c      print *,pbp1,pbp12,pbp21,pbp2
      write(13,fmt="(i3,1x)",advance='no') Ls
      do li=1,Nl
        do lj=1,Nl
!          rval=realpart(pbp(li,lj))/Nvol
          rval=real(pbp(li,lj),prc)/Nvol
          write(13,fmt="(1e12.4,1x)",advance='no') rval
          print *,li,lj,rval
        end do
      end do
      write(13,fmt="(1x)")
      close(13)
      print *,"Eval domain wall smeared condensate DONE"
      return
      end subroutine evalSmearedCondensateDomWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module smearedcondensate


