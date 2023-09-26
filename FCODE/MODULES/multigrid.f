      module multigridmod
      use pacc
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine MGLs(Nmax)
      use gaugefield
      use overlapmoduledev
      use options
      use rvmodule
      use condensatemodule
      implicit none
      integer Nmax
      integer Nht
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) BASE(Nv,4)
      type(sgnratfunc) :: SRF
      procedure(),pointer :: IMptr => NULL()
      complex(prc) pbp,pbpbase
      real err
      integer N


      call setRVs(Nv*4,R)
      N=4
      call setHTcoeffs(N,SRF)
      call DOverlap(R,BASE,u,.false.,baremass,SRF)
      OLTYPE=1
c      IMptr => IDOLS
      IMptr=>IDOverlap
c      call evalCondensateOLOld(u,baremass,1,pbpbase,SRF,IMptr)
      do while(N<=Nmax)    
        N=2*N
        call setHTcoeffs(N,SRF)
        call DOverlap(R,DR,u,.false.,baremass,SRF)
        err=maxval(abs(DR-BASE))
        print *,"err:",err
        BASE=DR
c        call evalCondensateOLOld(u,baremass,1,pbp,SRF,IMptr)
        print *,"pbp",pbp,pbp-pbpbase
        pbpbase=pbp
      end do



      

      return
      end subroutine MGLs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DW(Ncx,Ncy,Nct,dx,dy,dt,R,DR,u,DAGGER,mass)
      use gammas
      use indices
      implicit none
c     calculates DR = DWilson*R
      integer Ncx,Ncy,Nct
      real(prc) dx,dy,dt
      complex(prc),intent(in) :: R(Ncx,Ncy,Nct,4)
      complex(prc),intent(out) :: DR(Ncx,Ncy,Nct,4)
      complex(prc),intent(in) :: u(Ncx,Ncy,Nct,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      real(prc) mult
      complex(prc) cmass
      integer ix,iy,it,idirac
      integer mu,igork
      real evol

      mult=one
      if (DAGGER) then
        mult=-one
      end if
      evol=dx*dy*dt

c      cmass=mass+add
c      if (DAGGER) then
c        cmass=conjg(mass+add)
c      end if
!     mass term 
      DR = mass*R

!     Dirac term (anti-Hermitian)
      mu=1
      do idirac=1,4
        igork=gamin(mu,idirac)
        do ix=1,Ncx
          do iy=1,Ncy
            do it=1,Nct
              DR(ix,iy,it,idirac)=DR(ix,iy,it,idirac)+
     &            mult*gamval(mu,idirac)*
     &       ( u(ix,iy,it,mu)*R(np1(ix,Ncx),iy,it,igork)
     &     -conjg(u(nm1(ix,Ncx),iy,it,mu))*
     &            R(nm1(ix,Ncx),iy,it,igork) )/ (two*dx)
            end do
          enddo
        enddo
      enddo
      mu=2
      do idirac=1,4
        igork=gamin(mu,idirac)
        do ix=1,Ncx
          do iy=1,Ncy
            do it=1,Nct
              DR(ix,iy,it,idirac)=DR(ix,iy,it,idirac)+
     &            mult*gamval(mu,idirac)*
     &       ( u(ix,iy,it,mu)*R(ix,np1(iy,Ncy),it,igork)
     &     -conjg(u(ix,nm1(iy,Ncy),it,mu))*
     &            R(ix,nm1(iy,Ncy),it,igork) )/ (two*dy)
            end do
          enddo
        enddo
      enddo
      mu=3
      do idirac=1,4
        igork=gamin(mu,idirac)
        do ix=1,Ncx
          do iy=1,Ncy
            do it=1,Nct
              DR(ix,iy,it,idirac)=DR(ix,iy,it,idirac)+
     &            mult*gamval(mu,idirac)*
     &       ( u(ix,iy,it,mu)*R(ix,iy,np1(it,Nct),igork)
     &     -conjg(u(ix,iy,nm1(it,Nct),mu))*
     &            R(ix,iy,nm1(it,Nct),igork) )/ (two*dt)
            end do
          enddo
        enddo
      enddo

!     Wilson term (Hermitian)
      mu=1
      do idirac=1,4
        do ix=1,Ncx
          do iy=1,Ncy
            do it=1,Nct
              DR(ix,iy,it,idirac)=DR(ix,iy,it,idirac)
     &        -( u(ix,iy,it,mu)*R(np1(ix,Ncx),iy,it,idirac) 
     &               - two*R(ix,iy,it,idirac)
     &         +conjg(u(nm1(ix,Ncx),iy,it,mu))*
     &                R(nm1(ix,Ncx),iy,it,idirac))/(two*dx*dx)
            enddo
          enddo
        enddo
      enddo
      mu=2
      do idirac=1,4
        do ix=1,Ncx
          do iy=1,Ncy
            do it=1,Nct
              DR(ix,iy,it,idirac)=DR(ix,iy,it,idirac)
     &        -( u(ix,iy,it,mu)*R(ix,np1(iy,Ncy),it,idirac) 
     &               - two*R(ix,iy,it,idirac)
     &         +conjg(u(ix,nm1(iy,Ncy),it,mu))*
     &                R(ix,nm1(iy,Ncy),it,idirac))/(two*dy*dy)
            enddo
          enddo
        enddo
      enddo
      mu=3
      do idirac=1,4
        do ix=1,Ncx
          do iy=1,Ncy
            do it=1,Nct
              DR(ix,iy,it,idirac)=DR(ix,iy,it,idirac)
     &        -( u(ix,iy,it,mu)*R(ix,iy,np1(it,Nct),idirac) 
     &               - two*R(ix,iy,it,idirac)
     &         +conjg(u(ix,iy,nm1(it,Nct),mu))*
     &                R(ix,iy,nm1(it,Nct),idirac))/(two*dt*dt)
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine DW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine restrictPhi(Ncx,Ncy,Nct,phiU,phiL,uU,uL)
      use gammas
      use indices
      implicit none
c     calculates DR = DWilson*R
      integer Ncx,Ncy,Nct
      real(prc) dx,dy,dt
      complex(prc),intent(in) :: phiU(Ncx,Ncy,Nct,4)
      complex(prc),intent(out) :: phiL(Ncx/2,Ncy/2,Nct/2,4)
      complex(prc),intent(in) :: uU(Ncx,Ncy,Nct,3)
      complex(prc),intent(out) :: uL(Ncx/2,Ncy/2,Nct/2,3)
      integer ix,iy,it,idr,iux,iuy,iut

      do idr=1,4
        do ix=1,Ncx/2
          iux=2*ix-1
          do iy=1,Ncy/2
            iuy=2*iy-1
            do it=1,Nct/2
              iut=2*it-1
              phiL(ix,iy,it,idr)=(8*phiU(iux,iuy,iut,idr) +
     &                           4*phiU(np1(iux,Ncx),iuy,iut,idr) +
     &                           4*phiU(nm1(iux,Ncx),iuy,iut,idr) +
     &                           4*phiU(iux,np1(iuy,Ncy),iut,idr) +
     &                           4*phiU(iux,nm1(iuy,Ncy),iut,idr) +
     &                           4*phiU(iux,iuy,np1(iut,Nct),idr) +
     &                           4*phiU(iux,iuy,nm1(iut,Nct),idr) +
     &                       2*phiU(np1(iux,Ncx),np1(iuy,Ncy),iut,idr) +
     &                       2*phiU(np1(iux,Ncx),nm1(iuy,Ncy),iut,idr) +
     &                       2*phiU(nm1(iux,Ncx),np1(iuy,Ncy),iut,idr) +
     &                       2*phiU(nm1(iux,Ncx),nm1(iuy,Ncy),iut,idr) +
     &                       2*phiU(np1(iux,Ncx),iuy,np1(iut,Nct),idr) +
     &                       2*phiU(np1(iux,Ncx),iuy,nm1(iut,Nct),idr) +
     &                       2*phiU(nm1(iux,Ncx),iuy,np1(iut,Nct),idr) +
     &                       2*phiU(nm1(iux,Ncx),iuy,nm1(iut,Nct),idr) +
     &                       2*phiU(iux,np1(iuy,Ncy),np1(iut,Nct),idr) +
     &                       2*phiU(iux,np1(iuy,Ncy),nm1(iut,Nct),idr) +
     &                       2*phiU(iux,nm1(iuy,Ncy),np1(iut,Nct),idr) +
     &                       2*phiU(iux,nm1(iuy,Ncy),nm1(iut,Nct),idr) +
     &                phiU(np1(iux,Ncx),np1(iuy,Ncy),np1(iut,Nct),idr) +
     &                phiU(np1(iux,Ncx),np1(iuy,Ncy),nm1(iut,Nct),idr) +
     &                phiU(np1(iux,Ncx),nm1(iuy,Ncy),np1(iut,Nct),idr) +
     &                phiU(np1(iux,Ncx),nm1(iuy,Ncy),nm1(iut,Nct),idr) +
     &                phiU(nm1(iux,Ncx),np1(iuy,Ncy),np1(iut,Nct),idr) +
     &                phiU(nm1(iux,Ncx),np1(iuy,Ncy),nm1(iut,Nct),idr) +
     &                phiU(nm1(iux,Ncx),nm1(iuy,Ncy),np1(iut,Nct),idr) +
     &                phiU(nm1(iux,Ncx),nm1(iuy,Ncy),nm1(iut,Nct),idr))
     &                                                /64
              
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine restrictPhi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine prolongErr(Ncx,Ncy,Nct,errU,errL)
      use gammas
      use indices
      implicit none
      integer Ncx,Ncy,Nct
      real(prc) dx,dy,dt
      complex(prc),intent(out) :: errU(Ncx,Ncy,Nct,4)
      complex(prc),intent(in) :: errL(Ncx/2,Ncy/2,Nct/2,4)
      integer ix,iy,it,idr,iux,iuy,iut,xp1,yp1,tp1

      do idr=1,4
        do ix=1,Ncx/2
          iux=2*ix-1
          do iy=1,Ncy/2
            iuy=2*iy-1
            do it=1,Nct/2
              iut=2*it-1

              errU(iux,iuy,iut,idr)=errL(iux,iuy,iut,idr)
 
              xp1=1
              if (ix.eq.Ncx/2) then
                xp1=-iux
              endif
              yp1=1
              if (iy.eq.Ncy/2) then
                yp1=-iuy
              endif
              tp1=1
              if (it.eq.Nct/2) then
                tp1=-iut
              endif

              errU(iux+1,iuy,iut,idr)=(errL(iux,iuy,iut,idr) 
     &                                +errL(iux+xp1,iuy,iut,idr))/two 
              errU(iux,iuy+1,iut,idr)=(errL(iux,iuy,iut,idr) 
     &                                +errL(iux,iuy+yp1,iut,idr))/two 
              errU(iux,iuy,iut+1,idr)=(errL(iux,iuy,iut,idr) 
     &                                +errL(iux,iuy,iut+tp1,idr))/two 

              errU(iux+1,iuy+1,iut,idr)=(errL(iux,iuy,iut,idr) 
     &                                  +errL(iux+xp1,iuy,iut,idr)
     &                                  +errL(iux,iuy+yp1,iut,idr)
     &                                  +errL(iux+xp1,iuy+yp1,iut,idr))
     &                                       /(two*two) 
              errU(iux+1,iuy,iut+1,idr)=(errL(iux,iuy,iut,idr) 
     &                                  +errL(iux+xp1,iuy,iut,idr)
     &                                  +errL(iux,iuy,iut+tp1,idr)
     &                                  +errL(iux+xp1,iuy,iut+tp1,idr))
     &                                       /(two*two) 
              errU(iux,iuy+1,iut+1,idr)=(errL(iux,iuy,iut,idr) 
     &                                  +errL(iux,iuy,iut+tp1,idr)
     &                                  +errL(iux,iuy+yp1,iut,idr)
     &                                  +errL(iux,iuy+yp1,iut+tp1,idr))
     &                                       /(two*two) 

              errU(iux,iuy+1,iut+1,idr+1)=(errL(iux,iuy,iut,idr) 
     &                              +errL(iux+xp1,iuy,iut,idr)
     &                              +errL(iux,iuy+yp1,iut,idr)
     &                              +errL(iux,iuy,iut+tp1,idr)
     &                              +errL(iux+xp1,iuy+yp1,iut,idr)
     &                              +errL(iux+xp1,iuy,iut+tp1,idr)
     &                              +errL(iux,iuy+yp1,iut+tp1,idr)
     &                              +errL(iux+xp1,iuy+yp1,iut+tp1,idr))
     &                                       /(two*two*two) 
            end do
          end do
        end do
      end do

      return
      end subroutine prolongErr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine restrictU(Ncx,Ncy,Nct,phiU,phiL,uU,uL)
      use gammas
      use indices
      implicit none
c     calculates DR = DWilson*R
      integer Ncx,Ncy,Nct
      real(prc) dx,dy,dt
      complex(prc),intent(in) :: phiU(Ncx,Ncy,Nct,4)
      complex(prc),intent(out) :: phiL(Ncx/2,Ncy/2,Nct/2,4)
      complex(prc),intent(in) :: uU(Ncx,Ncy,Nct,3)
      complex(prc),intent(out) :: uL(Ncx/2,Ncy/2,Nct/2,3)
      integer ix,iy,it,idr,iux,iuy,iut


c     X component
      idr=1
      do ix=1,Ncx/2
        iux=2*ix-1
        do iy=1,Ncy/2
          iuy=2*iy-1
          do it=1,Nct/2
            iut=2*it-1
            uL(ix,iy,it,1)=(4*uU(iux,iuy,iut,idr) 
     &                    + 4*uU(iux+1,iuy,iut,idr) 
     &                    + 2*uU(iux,np1(iuy,Ncy),iut,idr) 
     &                    + 2*uU(iux+1,np1(iuy,Ncy),iut,idr) 
     &                    + 2*uU(iux,nm1(iuy,Ncy),iut,idr) 
     &                    + 2*uU(iux+1,nm1(iuy,Ncy),iut,idr) 
     &                    + 2*uU(iux,iuy,np1(iut,Nct),idr) 
     &                    + 2*uU(iux+1,iuy,np1(iut,Nct),idr) 
     &                    + 2*uU(iux,iuy,nm1(iut,Nct),idr) 
     &                    + 2*uU(iux+1,iuy,nm1(iut,Nct),idr) 
     &               + uU(iux,np1(iuy,Ncy),np1(iut,Nct),idr) 
     &               + uU(iux+1,np1(iuy,Ncy),np1(iut,Nct),idr) 
     &               + uU(iux,np1(iuy,Ncy),nm1(iut,Nct),idr) 
     &               + uU(iux+1,np1(iuy,Ncy),nm1(iut,Nct),idr) 
     &               + uU(iux,nm1(iuy,Ncy),np1(iut,Nct),idr) 
     &               + uU(iux+1,nm1(iuy,Ncy),np1(iut,Nct),idr) 
     &               + uU(iux,nm1(iuy,Ncy),nm1(iut,Nct),idr) 
     &               + uU(iux+1,nm1(iuy,Ncy),nm1(iut,Nct),idr))/32 
          enddo
        enddo
      enddo

c     Y component
      idr=2
      do ix=1,Ncx/2
        iux=2*ix-1
        do iy=1,Ncy/2
          iuy=2*iy-1
          do it=1,Nct/2
            iut=2*it-1
            uL(ix,iy,it,2)=(4*uU(iux,iuy,iut,idr) 
     &                    + 4*uU(iux,iuy+1,iut,idr) 
     &                    + 2*uU(np1(iux,Ncx),iuy,iut,idr) 
     &                    + 2*uU(np1(iux,Ncx),iuy+1,iut,idr) 
     &                    + 2*uU(nm1(iux,Ncx),iuy,iut,idr) 
     &                    + 2*uU(nm1(iux,Ncx),iuy+1,iut,idr) 
     &                    + 2*uU(iux,iuy,np1(iut,Nct),idr) 
     &                    + 2*uU(iux,iuy+1,np1(iut,Nct),idr) 
     &                    + 2*uU(iux,iuy,nm1(iut,Nct),idr) 
     &                    + 2*uU(iux,iuy+1,nm1(iut,Nct),idr) 
     &               + uU(np1(iux,Ncx),iuy,np1(iut,Nct),idr) 
     &               + uU(np1(iux,Ncx),iuy+1,np1(iut,Nct),idr) 
     &               + uU(np1(iux,Ncx),iuy,nm1(iut,Nct),idr) 
     &               + uU(np1(iux,Ncx),iuy+1,nm1(iut,Nct),idr) 
     &               + uU(nm1(iux,Ncx),iuy,np1(iut,Nct),idr) 
     &               + uU(nm1(iux,Ncx),iuy+1,np1(iut,Nct),idr) 
     &               + uU(nm1(iux,Ncx),iuy,nm1(iut,Nct),idr) 
     &               + uU(nm1(iux,Ncx),iuy+1,nm1(iut,Nct),idr))/32 

          enddo
        enddo
      enddo

c     T component
      idr=3
      do ix=1,Ncx/2
        iux=2*ix-1
        do iy=1,Ncy/2
          iuy=2*iy-1
          do it=1,Nct/2
            iut=2*it-1
            uL(ix,iy,it,3)=(4*uU(iux,iuy,iut,idr) 
     &                    + 4*uU(iux,iuy,iut+1,idr) 
     &                    + 2*uU(np1(iux,Ncx),iuy,iut,idr) 
     &                    + 2*uU(np1(iux,Ncx),iuy,iut+1,idr) 
     &                    + 2*uU(nm1(iux,Ncx),iuy,iut,idr) 
     &                    + 2*uU(nm1(iux,Ncx),iuy,iut+1,idr) 
     &                    + 2*uU(iux,np1(iuy,Ncy),iut,idr) 
     &                    + 2*uU(iux,np1(iuy,Ncy),iut+1,idr) 
     &                    + 2*uU(iux,nm1(iuy,Ncy),iut,idr) 
     &                    + 2*uU(iux,nm1(iuy,Ncy),iut+1,idr) 
     &               + uU(np1(iux,Ncx),np1(iuy,Ncy),iut,idr) 
     &               + uU(np1(iux,Ncx),np1(iuy,Ncy),iut+1,idr) 
     &               + uU(np1(iux,Ncx),nm1(iuy,Ncy),iut,idr) 
     &               + uU(np1(iux,Ncx),nm1(iuy,Ncy),iut+1,idr) 
     &               + uU(nm1(iux,Ncx),np1(iuy,Ncy),iut,idr) 
     &               + uU(nm1(iux,Ncx),np1(iuy,Ncy),iut+1,idr) 
     &               + uU(nm1(iux,Ncx),nm1(iuy,Ncy),iut,idr) 
     &               + uU(nm1(iux,Ncx),nm1(iuy,Ncy),iut+1,idr))/32 

          enddo
        enddo
      enddo

      return
      end subroutine restrictU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testMultigrid
      use gaugefield
      use basicdiracopsmod
      use rvmodule
      use options
      implicit none
      complex(prc) R(Nv,4),DR1(Nv,4),DR2(Nv,4)

      call setRVs(4*Nv,R)
      call DW(Ns,Ns,Nt,one,one,one,R,DR1,u,.false.,baremass)
      call DWilson(R,DR2,u,.false.,baremass,czero)
      print *,maxval(abs(DR1-DR2))
      print *,maxval(abs(DR1))

      return
      end subroutine testMultigrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module multigridmod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
