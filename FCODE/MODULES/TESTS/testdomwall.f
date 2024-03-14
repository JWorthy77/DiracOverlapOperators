      module testdomwallmod
      use pacc
      use arraysizes
      use numbers
      use options
      use rvmodule
      use ratfuncs
      use gaugefield
      use domainwallmod
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDomainWall()
      implicit none
      print *,""
      print *,"TEST DOMAIN WALL ROUTINES"
      print *,""
!      call testDomainWallHT()
      call testDomainWallZolo()
!      call testDomainWallWilson()
!      call testDomainWallShamir()
      return
      end subroutine testDomainWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDomainWallHT()
      use gammas
c      use dwcoeffs
      implicit none
      complex(prc),dimension(Nv,4,Ls) :: R5,DR5,TMP5 
      complex(prc),dimension(Nv,4) :: R,DR,TMP 
      integer gi,i1,i2

      call setRVs(Nv*4,R)
      call setRVs(Nv*4*Ls,R5)
      R5=czero
      R5(:,:,1)=R
      call PermM(R5,TMP5,.false.,4)
      call PermM(TMP5,DR5,.true.,4)
      print *,maxval(abs(TMP5-R5))
      print *,"PermM.PermMdag:",maxval(abs(DR5-R5))

c      omega=one
c      omega(1)=half

      MTYPE=1
      dwkernel=2
      TMP=R
      call mGmu(TMP,4)
      call KDDW4(TMP,DR,u,.false.,one/10)
      call mGmu(DR,4)
      call KDDW4(R,TMP,u,.true.,one/10)
      print *,"KDDW4 M1 G3 hermiticity:",maxval(abs(DR-TMP))
      TMP5=R5
      call mGmu5(TMP5,4)
      call KDDW(TMP5,DR5,u,.false.,one/10)
      call mGmu5(DR5,4)
      call KDDW(R5,TMP5,u,.true.,one/10)
      print *,"KDDW M1 G3 hermiticity:",maxval(abs(DR5-TMP5))

      return

      do i1=1,3,2
        MTYPE=i1
        do i2=1,2
          dwkernel=i2 ! 1:Shamir 2:Wilson
          print *,""
          print *,"MTYPE:",MTYPE,"DWKERNEL:",dwkernel
          print *,""
  
      call setRVs(Nv*4*Ls,R5)

      call DDW(R5,TMP5,u,.false.,baremass)
      call IDDW(TMP5,DR5,u,.false.,baremass)
      print *,"DDW.IDDW:",maxval(abs(DR5-R5))
      call DDW(R5,TMP5,u,.true.,baremass)
      call IDDW(TMP5,DR5,u,.true.,baremass)
      print *,"DDWdag.IDDWdag:",maxval(abs(DR5-R5))

      call KDDW(R5,TMP5,u,.false.,baremass)
      call IKDDW(TMP5,DR5,u,.false.,baremass)
      print *,"KDDW.IKDDW:",maxval(abs(DR5-R5))
      call KDDW(R5,TMP5,u,.true.,baremass)
      call IKDDW(TMP5,DR5,u,.true.,baremass)
      print *,"KDDWdag.IKDDWdag:",maxval(abs(DR5-R5))

        enddo
      enddo

      return

c      RDW=czero
c      call setRVs(Nv*4,RDW(:,:,1))
c      RDW(:,:,1)=cone
c      ROLS(:,:)=RDW(:,:,1)
c      ROLW(:,:)=RDW(:,:,1)
c      call KDDW(RDW,DRDW,u,.false.,baremass)
c
c      call setHTcoeffs(Ls,SRF)
c      call DOLS(ROLS,DROLS,u,.false.,baremass,SRF)
c      call DOverlap(ROLW,DROLW,u,.false.,baremass,SRF)
c      print *,"maxval Shamir:"
c      print *,maxval(abs(DROLS(:,:)))
c      print *,"difference Shamir:"
c      print *,maxval(abs(DROLS(:,:)-DRDW(:,:,1)))
c      print *,"difference Wilson:"
c      print *,maxval(abs(DROLW(:,:)-DRDW(:,:,1)))

    
        
c      call setRVs(Nv*4,R4)
c      call Action_DOL(R4,SB,u,baremass,SRF)

c      stop

c      MTYPE=1
c      DR4A=czero
cc      call DWilson(R4,DR4A,u,.false.,baremass,czero)
cc      call DOLS(R4,DR4A,u,.false.,baremass,SRF)
c      call DOverlap(R4,DR4A,u,.false.,baremass,SRF)
cc      call Action_DOL(R4,SA,u,baremass,SRF)
c      MTYPE=3
c      DR4B=czero
cc      call DWilson(R4,DR4B,u,.false.,baremass,czero)
cc      call DOLS(R4,DR4B,u,.false.,baremass,SRF)
c      call DOverlap(R4,DR4B,u,.false.,baremass,SRF)
cc      call Action_DOL(R4,SB,u,baremass,SRF)

c      print *,maxval(abs(DR4B-DR4A))
cc      print *,SA,SB,SA-SB

      return
      end subroutine testDomainWallHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDomainWallZolo()
      use zolomodule
      use dwcoeffs
      use overlapmoduledev
      use gammas
      implicit none
      complex(prc),dimension(Nv,4) :: R4,DR1,DR2,TMP 
      complex(prc),dimension(Nv,4,Ls) :: R5,DR5A,DR5B 
      type(zolotarev) :: zolo
      type(SGNratfunc) :: zolorf

      call setRVs(Nv*4*Ls,R5)
      call setRVs(Nv*4,R4)
!     set roots for domain wall coefficients
      call setZolo(1d-1,2d0,Ls,zolo)
      call getRoots(zolo)
      omega=zolo%roots
      print *,"omega:",omega
!     set coefficients for direct evaluation
      call setZoloCoeffs(2*Ls,zolorf,1d-1,2d0)

      

c      MTYPE=1
c      dwkernel=3
c      DR1=R4
c      call mGmu(DR1,4)
c      call KDDW4(DR1,DR2,u,.false.,one/10)
c      call mGmu(DR2,4)
c      call KDDW4(R4,DR1,u,.true.,one/10)
c      print *,"KDDW4 M1 G3 hermiticity:",maxval(abs(DR1-DR2))
c      MTYPE=1
c      dwkernel=3
c      DR5A=R5
c      call mGmu5(DR5A,4)
c      call DDW(DR5A,DR5B,u,.false.,one/10)
c      call mGmu5(DR5B,4)
c      call DDW(DR5A,DR5A,u,.true.,one/10)
c      print *,"DDW M1 G3 hermiticity:",maxval(abs(DR5A-DR5B))

      MTYPE=1
      dwkernel=2
      call DOLop(R4,DR1,u,.false.,zero,zolorf)
      print *,"DOL done"
      dwkernel=3
      call KDDW4(R4,DR2,u,.false.,zero)
      print *,"TEST equiv M1:",maxval(abs(DR1-DR2))

      MTYPE=1
      dwkernel=2
      call DOLop(R4,DR1,u,.false.,one/10,zolorf)
      print *,"DOL done"
      dwkernel=3
      call KDDW4(R4,DR2,u,.false.,one/10)
      print *,"TEST equiv M1:",maxval(abs(DR1-DR2))

c      MTYPE=1
c      dwkernel=2
c      call IDOLop(DR1,DR2,u,.false.,one/10,zolorf)
c      print *,"TEST:",maxval(abs(R4-DR2))

c      dwkernel=3
c      call IKDDW4(DR1,DR2,u,.false.,one/10)
c      print *,"TEST:",maxval(abs(R4-DR2))

c      MTYPE=1
c      dwkernel=2
c      call DOLop(R4,DR1,u,.true.,one/10,zolorf)
c      print *,"DOL done"
c      dwkernel=3
c      call KDDW4(R4,DR2,u,.true.,one/10)
c      print *,"TEST dag equiv M1:",maxval(abs(DR1-DR2))

      MTYPE=5
      dwkernel=2
      call DOLop(R4,DR1,u,.false.,zero,zolorf)
      print *,"DOL done"
      MTYPE=3
      dwkernel=3
      call KDDW4(R4,DR2,u,.false.,zero)
      print *,"TEST equiv M3/M5:",maxval(abs(DR1-DR2))


      MTYPE=5
      dwkernel=2
      call DOLop(R4,DR1,u,.true.,one/10,zolorf)
      print *,"DOL done"
      MTYPE=3
      dwkernel=3
      call KDDW4(R4,DR2,u,.true.,one/10)
      print *,"TEST dag equiv M3/M5:",maxval(abs(DR1-DR2))


      dwkernel=3
      call KDDW4(R4,DR1,u,.false.,one/10)
      call IKDDW4(DR1,DR2,u,.false.,one/10)
      print *,"TEST Zolo KDDW4.IKDDW4:",maxval(abs(R4-DR2))

      dwkernel=2
      call KDDW4(R4,DR1,u,.false.,one/10)
      call IKDDW4(DR1,DR2,u,.false.,one/10)
      print *,"TEST HT KDDW4.IKDDW4:",maxval(abs(R4-DR2))


      return


      MTYPE=1
      dwkernel=3
      call DDW(R5,DR5A,u,.false.,one/10)
      dwkernel=2
      call DDW(R5,DR5B,u,.false.,one/10)
      print *,"TEST:",maxval(abs(DR5A-DR5B))

      dwkernel=3
      call DDW(R5,DR5A,u,.true.,one/10)
      dwkernel=2
      call DDW(R5,DR5B,u,.true.,one/10)
      print *,"TESTdag:",maxval(abs(DR5A-DR5B))

      dwkernel=3
      call KDDW4(R4,DR1,u,.false.,one/10)
      dwkernel=2
      call KDDW4(R4,DR2,u,.false.,one/10)
      print *,"TEST K:",maxval(abs(DR1-DR2))

      dwkernel=3
      call KDDW4(R4,DR1,u,.true.,one/10)
      dwkernel=2
      call KDDW4(R4,DR2,u,.true.,one/10)
      print *,"TEST K dag:",maxval(abs(DR1-DR2))

      dwkernel=2
      call IDDW(R5,DR5B,u,.false.,one/10)
      dwkernel=3
      call IDDW(R5,DR5A,u,.false.,one/10)
      print *,"TEST:",maxval(abs(DR5A-DR5B))



c      MTYPE=1
c      dwkernel=3
c      call KDDW4(R4,DR1,u,.false.,one/10)
c      dwkernel=2
c      call DOLop(R4,DR2,u,.false.,one/10,zolorf)
c      print *,"KDDW4-DOL M1 Zolo:",maxval(abs(DR1-DR2))
      return
      end subroutine testDomainWallZolo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDomainWallWilson()
      use zolomodule
      use dwcoeffs
      use overlapmoduledev
      use gammas
      implicit none
      complex(prc),dimension(Nv,4) :: R4,DR1,DR2,TMP 
      complex(prc),dimension(Nv,4,Ls) :: R5,DR5A,DR5B 
      type(zolotarev) :: zolo
      type(SGNratfunc) :: zolorf

      MTYPE=3
      MDW=one
      baremass=0.0

      call setGRVs(3*Nv,theta)
      theta=theta*0.5
      call coef(u,theta)

      call setRVs(Nv*4*Ls,R5)
      call setRVs(Nv*4,R4)
!     set roots for domain wall coefficients
      call setZolo(1d-2,10d0,Ls,zolo)
      call getRoots(zolo)
      omega=zolo%roots
      print *,"omega: ",omega

      dwkernel=3
      call KDDW4(R4,DR1,u,.false.,baremass)
      dwkernel=2
      call KDDW4(R4,DR2,u,.false.,baremass)

      print *,"KDDW4 Z-HT: M3:",maxval(abs(DR1-DR2))

      return
      end subroutine testDomainWallWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDomainWallShamir()
      use zolomodule
      use dwcoeffs
      use overlapmoduleprod
      use gammas
      implicit none
      complex(prc),dimension(Nv,4) :: R4,DR1,DR2,TMP 
      complex(prc),dimension(Nv,4,Ls) :: R5,DR5A,DR5B 
      type(zolotarev) :: zolo
      type(sgnratfunc) :: HTcoeffs
      type(sgnratfunc) :: Zcoeffs


      MTYPE=3
      MDW=one
      baremass=0.0

      call setGRVs(3*Nv,theta)
      theta=theta*1.5
      call coef(u,theta)

      call setRVs(Nv*4*Ls,R5)
      call setRVs(Nv*4,R4)
!     set roots for domain wall coefficients
      call setZolo(1d-2,10d0,Ls,zolo)
      call getRoots(zolo)
      omega=zolo%roots
      print *,"omega: ",omega

      dwkernel=1
c      call KDDW4(R4,DR1,u,.false.,baremass)
      call setHTcoeffs(30,HTcoeffs)
      call DOLMS(R4,DR1,u,.false.,baremass,HTcoeffs)
      dwkernel=4
c      call KDDW4(R4,DR2,u,.false.,baremass)
      call setZoloCoeffs(30,Zcoeffs,1d-2,10d0)
      call DOLMS(R4,DR2,u,.false.,baremass,Zcoeffs)

      print *,"KDDW4 Z-HT: S M3:",maxval(abs(DR1-DR2))

      return
      end subroutine testDomainWallShamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testdomwallmod
