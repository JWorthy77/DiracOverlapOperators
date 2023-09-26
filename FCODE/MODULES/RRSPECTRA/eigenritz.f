      module eigenritzmod
      use arraysizes
      use rrparams
      use mathtoolsmod
      use gaugefield
      use basicdiracopsmod
      use overlapmoduledev
      contains

      SUBROUTINE eigenritz(lambda,mass,SRF)
      implicit none
      type(sgnratfunc),intent(in),optional :: SRF

      COMPLEX*16 ur(nev,Ns,Ns,Nt,4)
      COMPLEX*16 up(nev,Ns,Ns,Nt,4)
      COMPLEX*16 x(Ns,Ns,Nt,4)
      COMPLEX*16 y(Ns,Ns,Nt,4)
      COMPLEX*16 yt(Ns,Ns,Nt,4)
      COMPLEX*16 hm(nev,nev),p(nev),ee(nev)
      REAL*8 lambda(nev),ep(nev)

      COMPLEX*16 tempc
      REAL*8 chiral,hsum,mass,temp
      INTEGER i,j,k
      INTEGER i1,i2,i3,i4,id

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
          hm(i,j)=0.d0
          hm(i,j)=hm(i,j)+sum(dconjg(ur(i,:,:,:,:))*y)
          hm(j,i)=hm(i,j)
        ENDDO
      ENDDO

      call TRIDIAG (hm,nev,lambda,ee,ep,p)
      call TQLI(lambda,ep,nev,hm)      

      DO i=1,nev-1
         DO j=i+1,nev
            if(lambda(j).lt.lambda(i)) then
               temp=lambda(i)
               lambda(i)=lambda(j)
               lambda(j)=temp
               DO k=1,nev
                  tempc=hm(k,i)
                  hm(k,i)=hm(k,j)
                  hm(k,j)=tempc
               ENDDO
            endif
         ENDDO
      ENDDO
      DO i=1,nev
         lambda(i)=dsqrt(dabs(lambda(i)))
        write(25,*) mass,lambda(i)
      ENDDO

      DO i=1,nev
        up(i,:,:,:,:)=0.d0
        DO j=1,nev
          up(i,:,:,:,:)=up(i,:,:,:,:)+dconjg(hm(j,i))*ur(j,:,:,:,:)
        ENDDO
      ENDDO

      ur=up

      do i=1,nev
        do j=1,Nt
          hsum=0.d0
          do id=1,4
            do i1=1,Ns
              do i2=1,Ns
                hsum=hsum+dconjg(ur(i,i1,i2,j,id))
     &                          *ur(i,i1,i2,j,id)
              enddo
            enddo
          enddo
          write(27,*) hsum
        enddo
      enddo

      RETURN

      END

      end module eigenritzmod
