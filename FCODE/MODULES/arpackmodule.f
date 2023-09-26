      module arpackmodule
#ifdef USEARPACK      
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calceigs(which,nevs,peval,kernel,Nv)
c      use arraysizes
      implicit none

      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=250000, maxnev=160, maxncv=400, 
     &                                                   ldv=maxn)

      integer Nv

c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Complex*16 
     &                  ax(maxn), d(maxncv), 
     &                  v(ldv,maxncv), workd(3*maxn), 
     &                  workev(2*maxncv), resid(maxn), 
     &                  workl(3*maxncv*maxncv+5*maxncv)
      Double precision  
     &                  rwork(maxncv), rd(maxncv,3)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1
      character         which*2  ! 'SM'=smallest 'LM'=largest eigenvalues
      integer           ido, n, nx, nev, ncv, lworkl, info, ierr,
     &                  j, ishfts, maxitr, mode1, nconv
      integer           kernel ! 1=Wilson,2=Shamir
      Complex*16 
     &                  sigma
      Double precision 
     &                  tol
      logical           rvec
      integer           nevs
      real              peval(nevs) !eigenvalues
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision 
     &                  dznrm2 , dlapy2 
      external          dznrm2 , zaxpy , dlapy2  
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting mcaupd = 1                     |
c     %-------------------------------------------------%
c
c
c\SCCS Information: @(#) 
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer  logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
c      include 'debug.h' - contains the above content

      integer i

c      return

      ndigit = -3
      logfil = 6
      mcaitr = 0 
      mcapps = 0
      mcaupd = 0 ! 1 ! this is the output setting
      mcaup2 = 0
      mceigh = 0
      mceupd = 0
c
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
c
c      nx    = 10 
      n     = 4*Nv 

      do i=1,maxncv
        do j=1,ldv
          v(ldv,maxncv)=0
        enddo 
      enddo 
c
c     %-----------------------------------------------%
c     |                                               | 
c     | Specifications for ARPACK usage are set       | 
c     | below:                                        |
c     |                                               |
c     |    1) NEV = 4  asks for 4 eigenvalues to be   |  
c     |       computed.                               | 
c     |                                               |
c     |    2) NCV = 20 sets the length of the Arnoldi |
c     |       factorization                           |
c     |                                               |
c     |    3) This is a standard problem              |
c     |         (indicated by bmat  = 'I')            |
c     |                                               |
c     |    4) Ask for the NEV eigenvalues of          |
c     |       largest magnitude                       |
c     |         (indicated by which = 'LM')           |
c     |       See documentation in ZNAUPD  for the     |
c     |       other options SM, LR, SR, LI, SI.       | 
c     |                                               |
c     | Note: NEV and NCV must satisfy the following  |
c     | conditions:                                   |
c     |              NEV <= MAXNEV                    |
c     |          NEV + 2 <= NCV <= MAXNCV             |
c     |                                               |
c     %-----------------------------------------------%
c
      nev   = nevs
c      ncv   = 30 
      ncv   = 60 
      bmat  = 'I'
c
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling ZNAUPD                     |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from ZNAUPD . (see usage below)                 |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to ZNAUPD .                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID).  | 
c     |                                                     |
c     | The work array WORKL is used in ZNAUPD  as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl  = 3*ncv**2+5*ncv 
      tol    = 1d-9 
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting IPARAM(1) = 1).             |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | ZNAUPD .                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 20000
c      maxitr = 2000
c      mode1 = 2
      mode1 = 1
c
      iparam(1) = ishfts
c                
      iparam(3) = maxitr
c                  
      iparam(7) = mode1
c
c     %------------------------------------------------%
c     | M A I N   L O O P (Reverse Communication Loop) | 
c     %------------------------------------------------%
c
 10   continue
c   
c        %---------------------------------------------%
c        | Repeatedly call the routine ZNAUPD  and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%

         call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv,
     &                 v, ldv, iparam, ipntr, workd, workl, lworkl,
     &                 rwork,info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                                           |
c           |                y <--- A*x                 |
c           |                                           |
c           | The user should supply his/her own        |
c           | matrix vector multiplication routine here |
c           | that takes workd(ipntr(1)) as the input   |
c           | vector x , and returns the resulting      |
c           | matrix-vector product y = A*x in the      |
c           | array workd(ipntr(2)).                    | 
c           %-------------------------------------------%
c
c            call av (n, workd(ipntr(1)), workd(ipntr(2)))
            if (kernel == 1) then
              call GWils (n, workd(ipntr(1)), workd(ipntr(2)))
            else if (kernel == 2) then
              call Shamir (n, workd(ipntr(1)), workd(ipntr(2)))
            else if (kernel == 3) then
              call DW1 (n, workd(ipntr(1)), workd(ipntr(2)))
            else
              print *,'ARPACK kernel not specified'
              stop
            end if
c 
c           %-----------------------------------------%
c           | L O O P   B A C K to call ZNAUPD  again. |
c           %-----------------------------------------%
c
            go to 10
c
         endif
c 
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message, check the |
c        | documentation in ZNAUPD   |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd'
         print *, ' '
         stop
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using ZNEUPD .                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        |                                           |
c        | The routine ZNEUPD  now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1.)                                   |
c        |                                           |
c        %-------------------------------------------%
c           
         rvec = .false.
c

         call zneupd  (rvec, 'A', select, D, V, ldv, sigma,
     &        workev, bmat, n, which, nev, tol, resid, ncv,
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork, ierr)

      open(unit=11,file='DWilsonEigs.dat',access='append',
     &            status='unknown',form='formatted')
         if (kernel.ne.3) then
           print *,"SQRT EIGS ",which
           do j=1,nev
             print *,j,sqrt(D(j))
             write(11,*) which,j,D(j),abs(D(j)),sqrt(abs(D(j)))
           end do
           do j=1,nev
           if (which .eq. 'LM') then
c            print *,'D(nev):', D(nev)
             peval(j) = sqrt(abs(D(j)))
           elseif (which .eq. 'SM') then
c            print *,'D(nev):', D(nev)
             peval(j) = sqrt(abs(D(j)))
           else
c            print *,'D(nev):', D(nev)
             peval(j) = abs(D(j))
           endif
           end do
         else
c      open(unit=11,file='DWilsonEigs.dat',access='append',
c     &            status='unknown',form='formatted')
           print *,"EIGS ",which
           do j=1,nev
             print *,j,D(j),abs(D(j))
             write(11,*) which,j,D(j),abs(D(j))
             peval(j) = abs(D(j))
           end do
c      close(11)
         endif
      close(11)
c
c        %-----------------------------------------------%
c        | Eigenvalues are returned in the one           |
c        | dimensional array D and the corresponding     |
c        | eigenvectors are returned in the first        |
c        | NCONV (=IPARAM(5)) columns of the two         |
c        | dimensional array V if requested.  Otherwise, |
c        | an orthogonal basis for the invariant         |
c        | subspace corresponding to the eigenvalues in  |
c        | D is returned in V.                           |
c        %-----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of ZNEUPD . |
c            %------------------------------------%
c
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
c
         else
c
             nconv =  iparam(5)
             print *,"nconv:",nconv
             do 20 j=1, nconv
c
c                %---------------------------%
c                | Compute the residual norm |
c                |                           |
c                |   ||  A*x - lambda*x ||   |
c                |                           |
c                | for the NCONV accurately  |
c                | computed eigenvalues and  |
c                | eigenvectors.  (iparam(5) |
c                | indicates how many are    |
c                | accurate to the requested |
c                | tolerance)                |
c                %---------------------------%
c
c                 call av(n, v(1,j), ax)
                 if (kernel == 1) then
                   call GWils(n, v(1,j), ax)
                 else if (kernel == 2) then
                   call Shamir(n, v(1,j), ax)
                 else if (kernel == 3) then
                   call DW1(n, v(1,j), ax)
                 else
                   print *,'ARPACK kernel not set'
                   stop
                 end if

                 call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
                 rd(j,1) = dble (d(j))
                 rd(j,2) = dimag (d(j))
                 rd(j,3) = dznrm2 (n, ax, 1)
                 rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
 20          continue

             if (.false.) then 
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
             call dmout (6, nconv, 3, rd, maxncv, -6,
     &            'Ritz values (Real, Imag) and relative residuals')
             end if
         end if
c
c        %-------------------------------------------%
c        | Print additional convergence information. |
c        %-------------------------------------------%
c
       if (.true.) then 
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit',
     &                ' Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      
c
         print *, ' '
         print *, '_NSIMP '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', 
     &              nconv 
         print *, ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
      end if
c
      end if
c
c     %---------------------------%
c     | Done with program znsimp . |
c     %---------------------------%
c
 9000 continue
c
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DW1(n, v, w)
      use arraysizes
      use gaugefield
      use numbers
      use options
      use basicdiracopsmod
      implicit none
      integer n, j, lo
      Complex*16 v(n), w(n)

      call DWilson(v,w,u,.false.,-MDW,czero)
c      call DWilson(v,w,u,.false.,zero,czero)
c      w=v

      return
      end subroutine DW1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine GWILS(n, v, w)
      use arraysizes
      use gaugefield
      use numbers
      use options
      use basicdiracopsmod
      implicit none
      integer n, j, lo
      Complex*16 v(n), w(n), t(n)

      call DWilson(v,w,u,.false.,-MDW,czero)
c      call G5(w,v)
      call DWilson(w,t,u,.true.,-MDW,czero)
c      call G5(w,v)

      w=t+0*v

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SHAMIR(n, v, w)
      use arraysizes
      use gaugefield
      use numbers
      use options
      use basicdiracopsmod
      implicit none
      integer n, j, lo
      Complex*16 v(n), w(n)
      procedure(),pointer :: Mptr => NULL()

      call DShamir(v,w,u,.false.,-MDW,czero)      
      call DShamir(w,v,u,.true.,-MDW,czero)      
      w=v
      return

      call DWilson(v,w,u,.false.,-MDW,czero)      
c      call IDWDWD2(w,v,u)
      call IDW(w,v,u,.false.,-MDW,ctwo)
      call IDW(v,w,u,.true.,-MDW,ctwo)
      call DWilson(w,v,u,.true.,-MDW,czero)
      w=v

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine IDWDDWp2(n, v, w)
      use gaugefield
      use axbmodule1
      use options
      use basicdiracopsmod
      implicit none
      integer n, j, lo
      Complex*16 v(n), w(n)
      procedure(),pointer :: Mptr => NULL()

      Mptr => DdagDpC      
      call IM(v,w,u,.false.,-MDW,ctwo,Mptr)

      return
      end subroutine IDWDDWp2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine G5(v, w)
      use arraysizes
      implicit none
      integer n, j, lo, idirac, g5i
      Complex*16   v(Nv,4), w(Nv,4)

      print *,"G5 routine not implemented"
      stop

c      do j=1,kvol
c        do idirac=1,4
c          w(j,idirac)=2*gamval(g5i,idirac)*v(j,gamin(g5i,idirac))
c        enddo
c      enddo
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine testEigs()
      use gammas
      use gaugefield
      use gaugemodule
      implicit none
      real ev(1)
     
      call setGammas
#ifdef TWODIMENSIONS
      print *,"testEigs not ready for TWODIMENSIONS"
      stop
#else
      call setIndices
#endif
      theta=zero
      call makeGaugeField(.true.)

      call calceigs('SM',1,ev,1,Nv)
      print *,'ev min:',ev(1)
      call calceigs('LM',1,ev,1,Nv)
      print *,'ev max:',ev
      call calceigs('SM',1,ev,2,Nv)
      print *,'ev min:',ev
      call calceigs('LM',1,ev,2,Nv)
      print *,'ev max:',ev
      call calceigs('SM',1,ev,3,Nv)
      print *,'ev min:',ev
      call calceigs('LM',1,ev,3,Nv)
      print *,'ev max:',ev

      return
      end subroutine testEigs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine IDWDWD2(RR,DR,u)
      use pacc
      use arraysizes
      use numbers
      use options
      use basicdiracopsmod
      implicit none
      integer,parameter :: kferm=4*Nv
      integer,parameter :: niterc=10*Nv
      
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      complex(prc) cbetan

      integer nx,i,itercg
      real(prc) beta,betad,betan,alpha,alphad,alphan
c
c      resid=1d-8
      itercg=0

      do nx=1,kferm
        DR(nx)=cmplx(one,zero,prc)
      enddo

c
c      do 1 nx=1,niterc
      do 1 nx=1,1000
c      print *,"it:",nx

      itercg=itercg+1
      if(nx.gt.1) goto 51
c
c   initialise p=x, r=Phi(na)
c
      do 2 i=1,kferm
      p(i)=DR(i)
      r(i)=RR(i)
2     continue
      betad=1.0
      alpha=1.0
51    alphad=0.0

c      print *,"p:"
c      call ptosCfermion(p)
c      print *,"U:"
c      call ptosCboson(u)
  
c
c  x1=Mp
c
      call DWilson(p,x1,u,.false.,-MDW,ctwo)
c      print *,"x1",x1
c      call ptosCfermion(x1)
c      stop  
c
      if(nx.eq.1) goto 201
c
c   alpha=(r,r)/(p,(Mdagger)Mp)
c 
      alphad=0.0
      do 32 i=1,kferm
      alphad=alphad+conjg(x1(i))*x1(i)
32    continue
      alpha=alphan/alphad
c      
c   x=x+alpha*p
c
c      print *,"alphad: ",alphad
c      print *,"alphan: ",alphan
c      print *,"alpha: ",alpha
      do 4 i=1,kferm
      DR(i)=DR(i)+alpha*p(i)
4     continue
c      print *,"x:"
c      call ptosCfermion(x)
201   continue
c     
c   x2=(Mdagger)x1=(Mdagger)Mp
c
      call DWilson(x1,x2,u,.true.,-MDW,ctwo)

c      print *,"x2:"
c      call ptosCfermion(x2)  
  
c
c   r=r-alpha*(Mdagger)Mp
c
      do 6 i=1,kferm
      r(i)=r(i)-alpha*x2(i)
6     continue
c      print *,"r:"
c      call ptosCfermion(r)
c
c   beta=(r_k+1,r_k+1)/(r_k,r_k)
c
      betan=0.0 
      cbetan=(0.0,0.0) 
      do 61 i=1,kferm
      betan=betan+conjg(r(i))*r(i) 
      cbetan=cbetan+conjg(r(i))*r(i) 
61    continue 
      beta=betan/betad
      betad=betan
      alphan=betan
c      print *,"beta: ",betan
c      print *,"betad: ",betad
c      print *,"beta: ",beta
c
      if(nx.eq.1) beta=0.0
c
c   p=r+beta*p
c
      do 7 i=1,kferm
      p(i)=r(i)+beta*p(i)
7     continue

c      print *,"p:"
c      call ptosCfermion(p)
c      if (nx .eq. 2) stop

c      write(52,*) nx,betan,resid
c      print *,nx,betan,resid
c      print *,"x:"
c      call ptosCfermion(x)
      
      if(betan.lt.resid*1000000) goto 8
1     continue
8     continue
c      print *,nx,betan,resid
c      print *,"x:"
c      call ptosCfermion(x)
c      stop
      return
      end
#endif      
      end module arpackmodule


