      module rrspectrum
      use arraysizes
      use numbers
      use ratfuncs
      use rrparams
      use ritzmod
      use diagritzmod
      use eigenritzmod
      implicit none

      contains

      subroutine spectrum(mass,SRF,lambda)
      use gaugefield
      use basicdiracopsmod
      use overlapmoduledev
      use options
      implicit none
      REAL*8 mass
      type(sgnratfunc),intent(in),optional :: SRF
      REAL*8 lambda(nev)

      COMPLEX*16 ur(nev,Ns,Ns,Nt,4)
      COMPLEX*16 x(Ns,Ns,Nt,4)
      COMPLEX*16 yt(Ns,Ns,Nt,4)
      COMPLEX*16 y(Ns,Ns,Nt,4)
      COMPLEX*16 p(Ns,Ns,Nt,4)
      COMPLEX*16 t(Ns,Ns,Nt,4)
      COMPLEX*16 overc
      REAL*8 over,muritz,pnorm,gnorm,gstart,w
      INTEGER ks,iev,i,iter_cum,iflag
      INTEGER is,ic,i1,i2,i3,i4

      common /ritzev/ ur
      common /hwoper/ iter_cum

      procedure(),pointer :: Vptr => NULL()
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: Dptr => NULL()

      if (dwkernel.eq.1) then
        SPECTRA_KERNEL=2
      elseif (dwkernel.eq.2) then
        SPECTRA_KERNEL=1
      else
        print *,"dwkernel not set"
        stop
      endif

!     KERNEL must be set as well
      if (SPECTRA_KERNEL.eq.1) then
        Vptr => VOLMpf
        Mptr => DdagDpC
        Dptr => DWilson
        print *,"Calculate Eigenvalues for Wilson Overlap Operator"
      elseif (SPECTRA_KERNEL.eq.2) then
        Vptr => VOLSpf
!        Mptr => SdagSpC ! do not assign these
!        Dptr => DShamir ! do not assign these
        print *,"Calculate Eigenvalues for Shamir Overlap Operator"
      endif

      iter_cum=0
      iflag=0

      print *,"Calc Spectrum"
  
      call startritz()

      DO ks=1,nks
        DO iev=1,nev
          flush(88)
          x=ur(iev,:,:,:,:)
          CALL GS(x,iev)  

          over=sum(x*dconjg(x))

          if(over.lt.1.e-18) then
            write(99,*) 'zero vector over at start',over,ks,iev
          endif
          over=1.0/sqrt(over)
          x = over*x
          ur(iev,:,:,:,:) = x

          if (KERNEL.and.(SPECTRA_KERNEL.eq.1)) then
            CALL DWilson(x,yt,u,.false.,mass,czero) 
            CALL DWilson(yt,y,u,.true.,mass,czero) 
          elseif (KERNEL.and.(SPECTRA_KERNEL.eq.2)) then
            CALL DShamir(x,yt,u,.false.,mass,czero) 
            CALL DShamir(yt,y,u,.true.,mass,czero) 
          elseif (.not.KERNEL) then
            call VOLNK(x,y,u,.false.,mass,SRF,Vptr,Mptr,Dptr)
          endif
c          print *,"VOLNK 1"
          muritz=sum(dconjg(x)*y)
          p=y-muritz*x

          CALL GS(p,iev)  
          pnorm=sum(p*dconjg(p))
          pnorm=dsqrt(pnorm)
          if(pnorm.lt.1.e-18) then
            write(99,*) 'zero search at start',pnorm,ks,iev
          endif
          gnorm=pnorm 
          gstart=gnorm
          DO i=2,np
c            print *,"ks:",ks,"iev:",iev,"i:",i
c            print *,"gnorm:",gnorm
            CALL ritz(x,y,p,iev,muritz,pnorm,gnorm,mass,SRF)
            w=w2*muritz
            if(w1.gt.w) w=w1
            if(i.gt.10) then
c              print *,"gnorm:",gnorm,"w:",w,"t2:",gstart*cv_fac
              if(gnorm.lt.w) go to 2
              if(gnorm.lt.gstart*cv_fac) go to 2
            endif
            if(mod(i,10).eq.0) then
              CALL GS(x,iev)  
              over=sum(x*dconjg(x))
              if(over.lt.1.e-18) then
                write(99,*) 'zero vector over',over,ks,iev,i
              endif
              over=sqrt(over)
              x = (1.d0/over)*x
              p = over *p
              if (KERNEL.and.(SPECTRA_KERNEL.eq.1)) then
                CALL DWilson(x,yt,u,.false.,mass,czero) 
                CALL DWilson(yt,y,u,.true.,mass,czero) 
              elseif (KERNEL.and.(SPECTRA_KERNEL.eq.2)) then
                CALL DShamir(x,yt,u,.false.,mass,czero) 
                CALL DShamir(yt,y,u,.true.,mass,czero) 
              elseif (.not.KERNEL) then
                call VOLNK(x,y,u,.false.,mass,SRF,Vptr,Mptr,Dptr)
              endif
c              print *,"VOLNK 2"
              muritz=sum(dconjg(x)*y)
              gnorm=0.d0
              t=y-muritz*x
              gnorm=sum(t*dconjg(t))
              gnorm=dsqrt(gnorm)
              overc=sum(dconjg(x)*p)
              p=p-overc*x
              overc=sum(dconjg(y)*p)
              overc=(overc/(gnorm*gnorm))-1.d0
              p=p-overc*t
              CALL GS(p,iev)  
              pnorm=0.d0
              pnorm=sum(p*dconjg(p))
              pnorm=dsqrt(pnorm)
            endif
            if(pnorm.lt.1.e-18) then
              write(99,*) 'zero search in Ritz',pnorm,ks,iev,i
            endif
            print *,"muritz:",muritz
          ENDDO

 2        continue
          if(i.gt.np) then
            write(88,*) 'No of Ritz steps=', np
          else
            write(88,*) 'No of Ritz steps=', i
          endif
          ur(iev,:,:,:,:)=x

        ENDDO
        CALL diagritz(lambda,iflag,ks,mass,SRF)

        if(iflag.eq.1) go to 1

      ENDDO

 1    continue
      if(ks.gt.nks) then
         write(88,*) 'No of KS steps= ', nks
      else
         write(88,*) 'No of KS steps= ', ks
      endif
      write(88,*) 'mass=', mass, ';   No of total action of Hw=', 
     &     iter_cum

      CALL eigenritz(lambda,mass,SRF)
         

      return
      end

      end module rrspectrum



