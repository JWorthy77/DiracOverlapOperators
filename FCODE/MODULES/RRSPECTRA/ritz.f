      module ritzmod
      use arraysizes
      use numbers
      use rrparams
      use gaugefield
      use ratfuncs
      use basicdiracopsmod
      use overlapmoduledev
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE startritz()
      use rvmod
      implicit none

      INTEGER iv,is,i1,i2,i3
      COMPLEX*16 ur(nev,Ns,Ns,Nt,4)
      complex*16 ai
c      real*8 sranf

      common /ritzev/ ur

      ai=(0.d0,1.d0)

      do iv=1,nev
        do i1=1,Ns
          do i2=1,Ns
            do i3=1,Nt
              do is=1,4
                ur(iv,i1,i2,i3,is) = sranf()+ai*sranf()
               enddo
             enddo
           enddo
         enddo
       enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ritz(x,y,p,iev,muritz,pnorm,gnorm,mass,SRF)
      implicit none
      type(sgnratfunc),intent(in) :: SRF
      COMPLEX*16 ur(nev,Ns,Ns,Nt,4)
      COMPLEX*16 x(Ns,Ns,Nt,4)
      COMPLEX*16 y(Ns,Ns,Nt,4)
      COMPLEX*16 p(Ns,Ns,Nt,4)
      COMPLEX*16 t(Ns,Ns,Nt,4)
      COMPLEX*16 z(Ns,Ns,Nt,4)
      COMPLEX*16 zt(Ns,Ns,Nt,4)

      COMPLEX*16 e1c
      REAL*8 muritz,pnorm,gnorm,tnorm,e1,mass
      REAL*8 s1,s2,s3,a,cosdel,sindel,costh,sinth,betar

      INTEGER iev
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

      if (KERNEL.and.(SPECTRA_KERNEL.eq.1)) then
        call DWilson(p,zt,u,.false.,mass,czero)
        call DWilson(zt,z,u,.true.,mass,czero)
      elseif (KERNEL.and.(SPECTRA_KERNEL.eq.2)) then
        CALL DShamir(p,zt,u,.false.,mass,czero) 
        CALL DShamir(zt,z,u,.true.,mass,czero) 
      elseif (.not.KERNEL) then
        call VOLNK(p,z,u,.false.,mass,SRF,Vptr,Mptr,Dptr)
      endif

      e1=sum(dconjg(p)*z)
      s1=0.5d0*(muritz+(e1/(pnorm*pnorm)))
      s2=0.5d0*(muritz-(e1/(pnorm*pnorm)))
      s3=(gnorm*gnorm)/pnorm
      if(dabs(s2).ge.s3) then
         a=dabs(s2)*dsqrt(1.d0+(s3/s2)**2)
      else
         a=s3*dsqrt(1.d0+(s2/s3)**2)
      endif
      cosdel=s2/a
      sindel=s3/a
      if(cosdel.gt.0) then
         sinth=-sqrt(0.5d0*(1.d0+cosdel))
         costh=-0.5d0*sindel/sinth
      else
         costh=sqrt(0.5d0*(1.d0-cosdel))
         sinth=-0.5d0*sindel/costh
      endif

      muritz=muritz-2.d0*a*sinth*sinth

      sinth=sinth/pnorm

      x = costh*x + sinth*p
      y = costh*y + sinth*z

      z = y - muritz*x

      tnorm = sum(z*dconjg(z))
      tnorm = dsqrt(tnorm)

      betar=(costh*tnorm*tnorm)/(gnorm*gnorm)
      gnorm=tnorm
      if((betar*pnorm*costh).gt.(20.d0*gnorm))betar=0.0

      e1c = sum(dconjg(x)*p)
      p = p-e1c*x
      p = z + betar*p
      
      CALL GS(p,iev)  
      
      pnorm=sum(p*dconjg(p))
      pnorm=DSQRT(pnorm)

      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GS(x,iv)
      implicit none

      COMPLEX*16 ur(nev,Ns,Ns,Nt,4)
      COMPLEX*16 x(Ns,Ns,Nt,4)

      COMPLEX*16 overc
      INTEGER iv,i

      common /ritzev/ ur

      DO i=1,iv-1
        overc=sum(dconjg(ur(i,:,:,:,:))*x)
        x=x-overc*ur(i,:,:,:,:)
      ENDDO

      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module ritzmod
