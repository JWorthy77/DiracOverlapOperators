      module diagritzmod
      use arraysizes
      use rrparams
      use mathtoolsmod
      use gaugefield
      use ratfuncs
      use basicdiracopsmod
      use overlapmoduledev
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE diagritz(lambda,iflag,ks,mass,SRF)
      implicit none
      type(sgnratfunc),intent(in),optional :: SRF
      COMPLEX*16 ur(nev,Ns,Ns,Nt,4)
      COMPLEX*16 up(nev,Ns,Ns,Nt,4)
      COMPLEX*16 x(Ns,Ns,Nt,4)
      COMPLEX*16 t(Ns,Ns,Nt,4)
      COMPLEX*16 y(Ns,Ns,Nt,4)
      COMPLEX*16 yt(Ns,Ns,Nt,4)
      COMPLEX*16 hm(nev,nev),p(nev),ee(nev)
      REAL*8 d(nev),ep(nev)
      REAL*8 lambda(nev)
      COMPLEX*16 tempc
      REAL*8 temp,mass
      REAL*8 diff,w
      INTEGER ks,i,j,k,iflag
      INTEGER is,ic,i1,i2,i3,i4

      common /ritzev/ ur

      procedure(),pointer :: Vptr => NULL()
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: Dptr => NULL()

      if (SPECTRA_KERNEL.eq.1) then
        Vptr => VOLMpf
        Mptr => DdagDpC
        Dptr => DWilson 
      elseif (SPECTRA_KERNEL.eq.2) then
        Vptr => VOLSpf
c        Mptr => SdagSpC
c        Dptr => DShamir 
      endif

      DO j=1,nev
        x=ur(j,:,:,:,:)
        if (KERNEL.and.(SPECTRA_KERNEL.eq.1)) then
          CALL DWilson(x,yt,u,.false.,mass,czero)      
          CALL DWilson(yt,y,u,.true.,mass,czero)      
        elseif (KERNEL.and.(SPECTRA_KERNEL.eq.2)) then
          CALL DShamir(x,yt,u,.false.,mass,czero) 
          CALL DShamir(yt,y,u,.true.,mass,czero) 
        elseif (.not.KERNEL) then
          call VOLNK(x,y,u,.false.,mass,SRF,Vptr,Mptr,Dptr)
        endif
        DO i=1,j
          hm(i,j)=sum(dconjg(ur(i,:,:,:,:))*y)
          hm(j,i)=dconjg(hm(i,j))
        ENDDO
      ENDDO

      call TRIDIAG (hm,nev,d,ee,ep,p)
      call TQLI(d,ep,nev,hm)      
      
      DO i=1,nev-1
         DO j=i+1,nev
            if(d(j).lt.d(i)) then
               temp=d(i)
               d(i)=d(j)
               d(j)=temp
               DO k=1,nev
                  tempc=hm(k,i)
                  hm(k,i)=hm(k,j)
                  hm(k,j)=tempc
               ENDDO
            endif
         ENDDO
      ENDDO
      DO i=1,nev
         d(i)=dsqrt(dabs(d(i)))
      ENDDO
      iflag=1
      DO i=1,nev
         if(ks.gt.1) then
            diff=dabs(DABS(d(i))-lambda(i))
            w=10.d0*w2*dabs(lambda(i))
	    if(10.d0*w1.gt.w) w=10.d0*w1
	    if(diff.gt.w) iflag=0
         else
            iflag=0
         endif
         lambda(i)=DABS(d(i))
         write(88,*) 'hov2', i,lambda(i)
      ENDDO

      DO i=1,nev
        up(i,:,:,:,:)=0.d0
        DO j=1,nev
          up(i,:,:,:,:)=up(i,:,:,:,:)+dconjg(hm(j,i))*ur(j,:,:,:,:)
        ENDDO
      ENDDO

      DO i=1,nev
        ur(i,:,:,:,:)=up(i,:,:,:,:)
      ENDDO

      RETURN
      END
      
      end module diagritzmod
