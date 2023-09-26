      module VOLmodule
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use axbmodule1
      use basicdiracopsmod
      implicit none
      logical,parameter :: VB_VOL=.false.
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLGpf(RR,S,u,DAGGER,dwmass,SRF,Mptr,Dptr)
      use options
!     approximate Voverlap using partial fraction rational functions
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Mptr
      procedure(),pointer :: Dptr
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      real(prc) mult

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      nf=nn-nd+1
      S=zero
      if (nf.ge.1) then ! it should only ever be 0 or 1
        do j=0,nf-1
          mult=front(j+1)
!         T = m.DWilson^j.RR
          TMP1=mult*RR
          do p=1,j
            call Mptr(TMP1,TMP2,u,DAGGER,dwmass,czero)
            TMP1=TMP2
          end do
          S=S+TMP1
        end do
      end if

      do j=1,nd
         add = -denom(j)
         call IM(RR,TMP1,u,DAGGER,dwmass,add,Mptr)
         TMP2 = pf(j)*TMP1
         S=S+TMP2
      end do

      call Dptr(S,TMP1,u,DAGGER,dwmass,czero)
      S=SRF%mult*TMP1

      end associate
      return
      end subroutine VOLGpf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLMpf(RR,S,u,DAGGER,dwmass,SRF,Mptr,Dptr)
      use options
!     approximate Voverlap using partial fraction rational functions
!     with a multishift cg method
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Mptr
      procedure(),pointer :: Dptr
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      real(prc) mult

      if(VB_VOL)then ; print *,"VOLMpf" ; endif

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      if(VB_VOL)then ; print *,"nn",nn,"nd",nd ; endif
      nf=nn-nd+1
      S=zero
      if (nf.ge.1) then ! it should only ever be 0 or 1
        do j=0,nf-1
          mult=front(j+1)
!         T = m.DWilson^j.RR
          TMP1=mult*RR
          do p=1,j
            call Mptr(TMP1,TMP2,u,DAGGER,dwmass,czero)
            TMP1=TMP2
          end do
          S=S+TMP1
        end do
      end if

c      print *,"MSCG"
      call MSCG(RR,S,u,DAGGER,dwmass,SRF,Mptr)

      call Dptr(S,TMP1,u,DAGGER,dwmass,czero)
      S=SRF%mult*TMP1

      end associate
      return
      end subroutine VOLMpf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLMWpf(RR,S,u,DAGGER,dwmass,SRF) ! problems with Intel compiler with pointy version
      use options
!     approximate Wilson Voverlap using partial fraction rational functions with a multishift cg method
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      real(prc) mult
      if(VB_VOL)then ; print *,"VOLMWpf" ; endif

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      if(VB_VOL)then ; print *,"nn",nn,"nd",nd ; endif
      nf=nn-nd+1
      S=zero
      if (nf.ge.1) then ! it should only ever be 0 or 1
        S=front(1)*RR
      end if
      call MSCGW(RR,S,u,DAGGER,dwmass,SRF) ! S=S+...
      call DWilson(S,TMP1,u,DAGGER,dwmass,czero)
      S=SRF%mult*TMP1

      end associate
      return
      end subroutine VOLMWpf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLMSpf(RR,S,u,DAGGER,dwmass,SRF) ! problems with Intel compiler with pointy version
      use options
!     approximate Shamir Voverlap using partial fraction rational functions with a multishift cg method
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      real(prc) mult

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      nf=nn-nd+1
      S=zero
      if (nf.ge.1) then ! it should only ever be 0 or 1
        S=front(1)*RR
      end if
      call MSCGS(RR,S,u,DAGGER,dwmass,SRF) ! S=S+...
      call DShamir(S,TMP1,u,DAGGER,dwmass,czero)
      S=SRF%mult*TMP1

      end associate
      return
      end subroutine VOLMSpf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine MSCG(rhs,DR,u,DAGGER,dwmass,SRF,Mptr)
      use options
      use countmod
!     a multishift cg method solves all (Mptr+sigma_i)Z_i=rhs
      implicit none
      complex(prc) rhs(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Mptr
      complex(prc) Zm1ns(Nv,4)
      integer ns
      integer k,o,it,maxit,nsr
      complex(prc),dimension(:,:,:),allocatable :: Z,Zp1,P,Pp1
      real(prc),dimension(:),allocatable :: gamm1,gam,gamp1,sig,ocount
      complex(prc) R(Nv,4),Rp1(Nv,4),B(Nv,4)
      real(prc) alpham1,alpha,beta,betap1,dnm,conv,cns,drr,drp1rp1

      if (VB_VOL) then ; print *,"MSCG" ; end if

      maxit=10*Nv
      associate(sigma => SRF%pfrf%denom%zeros,
     &          mult => SRF%pfrf%pf)

      ns=size(sigma)
      nsr=ns
      allocate(Z(Nv,4,ns),P(Nv,4,ns))
      allocate(gamm1(ns),gam(ns),gamp1(ns),sig(ns),ocount(ns))
      Z=0
      sig=-sigma
      ocount=0
      do k=1,nsr
        P(:,:,k)=rhs
      end do
      R=rhs
      alpham1=one
      beta=zero
      gamm1=one
      gam=one
      gamp1=one
      dnm=dopr(rhs,rhs)
      do it=1,maxit
        call Mptr(P(:,:,1),B,u,DAGGER,dwmass,czero)
        B=B+sig(1)*P(:,:,1)
        drr=dopr(R,R)
        alpha=drr/dopr(P(:,:,1),B)
        Rp1=R-alpha*B
        drp1rp1=dopr(Rp1,Rp1)
        betap1=drp1rp1/drr
        Z(:,:,1)=Z(:,:,1)+alpha*P(:,:,1)
        P(:,:,1)=Rp1+betap1*P(:,:,1)
        Zm1ns=Z(:,:,nsr)
        do o=2,nsr
          gamp1(o)=gam(o)*gamm1(o)*alpham1/
     &      ( alpha*beta*(gamm1(o)-gam(o)) +
     &          gamm1(o)*alpham1*(one+alpha*(sig(o)-sig(1))) )
          Z(:,:,o)=Z(:,:,o)+alpha*gamp1(o)/gam(o)*P(:,:,o)
          P(:,:,o)=gamp1(o)*Rp1+betap1*(gamp1(o)/gam(o))**2*P(:,:,o)
        end do
        beta=betap1
        conv=sqrt(drp1rp1/dnm)
       
c        print *,maxval(abs(P))
c        print *,maxval(abs(Z))
c        print *,gamp1
c        print *,gamp1/gam
c        print *,"alpha:",alpha,"beta:",beta
c        print *,it," of maxit:",maxit,"conv:",conv
        cns=maxval(abs(Z(:,:,nsr)-Zm1ns))/maxval(abs(Z(:,:,nsr)))
        if (cns.lt.1d-18) then
          ocount(nsr)=it
          nsr=nsr-1
c          print *,"ns:",ns
          if (nsr.eq.0) then
            goto 501
          endif
        endif
        if (conv.lt.1d-10) then
          goto 501
        endif
        alpham1=alpha
        R=Rp1
        gamm1=gam
        gam=gamp1
      end do

501   continue
      if (VB_VOL) then 
        print *,it," of maxit:",maxit,"conv:",conv 
      end if
c      print *,"ocount:",ocount
c      ns=size(sigma)
!      ic_idx=ic_idx+1
!      inner_count(ic_idx)=it
      do o=1,ns
        DR=DR+mult(o)*Z(:,:,o)
      end do

c      print *,maxval(abs(Z)),minval(abs(Z))
c      print *,maxval(abs(DR)),minval(abs(DR))
c      print *,mult
c      stop     

      deallocate(gamm1,gam,gamp1,sig,ocount)
      deallocate(Z,P)
      end associate
      return
      end subroutine MSCG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine MSCGW(rhs,DR,u,DAGGER,dwmass,SRF)
      use options
      use countmod
!     a multishift cg method solves all (DdagD+sigma_i)Z_i=rhs
      implicit none
      complex(prc) rhs(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
c      procedure(),pointer :: Mptr
      complex(prc) Zm1ns(Nv,4)
      integer ns
      integer k,o,it,maxit,nsr
      complex(prc),dimension(:,:,:),allocatable :: Z,Zp1,P,Pp1
      real(prc),dimension(:),allocatable :: gamm1,gam,gamp1,sig,ocount
      complex(prc) R(Nv,4),Rp1(Nv,4),B(Nv,4)
      real(prc) alpham1,alpha,beta,betap1,dnm,conv,cns,drr,drp1rp1

      if(VB_VOL)then ; print *,"MSCGW" ; endif
      maxit=10*Nv
      associate(sigma => SRF%pfrf%denom%zeros,
     &          mult => SRF%pfrf%pf)

      if(VB_VOL)then ; print *,"sigma:",sigma ; endif
      if(VB_VOL)then ; print *,"mult:",mult ; endif
      ns=size(sigma)
      nsr=ns
      allocate(Z(Nv,4,ns),P(Nv,4,ns))
      allocate(gamm1(ns),gam(ns),gamp1(ns),sig(ns),ocount(ns))
      Z=0
      sig=-sigma
      ocount=0
      do k=1,nsr
        P(:,:,k)=rhs
      end do
      R=rhs
      alpham1=one
      beta=zero
      gamm1=one
      gam=one
      gamp1=one
      dnm=dopr(rhs,rhs)
      if(VB_VOL)then ; print *,"dnm:",dnm ; endif
      do it=1,maxit
        call DdagDpC(P(:,:,1),B,u,DAGGER,dwmass,czero)
        B=B+sig(1)*P(:,:,1)
        drr=dopr(R,R)
        alpha=drr/dopr(P(:,:,1),B)
        Rp1=R-alpha*B
        drp1rp1=dopr(Rp1,Rp1)
        if(VB_VOL)then ; print *,"drp1rp1:",drp1rp1 ; endif
        betap1=drp1rp1/drr
        Z(:,:,1)=Z(:,:,1)+alpha*P(:,:,1)
        P(:,:,1)=Rp1+betap1*P(:,:,1)
        Zm1ns=Z(:,:,nsr)
        do o=2,nsr
          gamp1(o)=gam(o)*gamm1(o)*alpham1/
     &      ( alpha*beta*(gamm1(o)-gam(o)) +
     &          gamm1(o)*alpham1*(one+alpha*(sig(o)-sig(1))) )
          Z(:,:,o)=Z(:,:,o)+alpha*gamp1(o)/gam(o)*P(:,:,o)
          P(:,:,o)=gamp1(o)*Rp1+betap1*(gamp1(o)/gam(o))**2*P(:,:,o)
        end do
        beta=betap1
        conv=sqrt(drp1rp1/dnm)
        if(VB_VOL)then ; print *,"beta:",beta,"conv:",conv ; endif
       
        cns=maxval(abs(Z(:,:,nsr)-Zm1ns))/maxval(abs(Z(:,:,nsr)))
        if (cns.lt.1d-18) then
          ocount(nsr)=it
          nsr=nsr-1
          if (nsr.eq.0) then
            goto 501
          endif
        endif
        if (conv.lt.1d-10) then
          goto 501
        endif
        alpham1=alpha
        R=Rp1
        gamm1=gam
        gam=gamp1
      end do

501   continue
      do o=1,ns
        DR=DR+mult(o)*Z(:,:,o)
      end do

      ic_idx=ic_idx+1
      inner_count=inner_count+it

      deallocate(gamm1,gam,gamp1,sig,ocount)
      deallocate(Z,P)
      end associate
      return
      end subroutine MSCGW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine MSCGS(rhs,DR,u,DAGGER,dwmass,SRF)
      use options
!     a multishift cg method solves all (SdagS+sigma_i)Z_i=rhs
      implicit none
      complex(prc) rhs(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Mptr
      complex(prc) Zm1ns(Nv,4)
      integer ns
      integer k,o,it,maxit,nsr
      complex(prc),dimension(:,:,:),allocatable :: Z,Zp1,P,Pp1
      real(prc),dimension(:),allocatable :: gamm1,gam,gamp1,sig,ocount
      complex(prc) R(Nv,4),Rp1(Nv,4),B(Nv,4)
      real(prc) alpham1,alpha,beta,betap1,dnm,conv,cns,drr,drp1rp1

      maxit=10*Nv
      associate(sigma => SRF%pfrf%denom%zeros,
     &          mult => SRF%pfrf%pf)

      ns=size(sigma)
      nsr=ns
      allocate(Z(Nv,4,ns),P(Nv,4,ns))
      allocate(gamm1(ns),gam(ns),gamp1(ns),sig(ns),ocount(ns))
      Z=0
      sig=-sigma
      ocount=0
      do k=1,nsr
        P(:,:,k)=rhs
      end do
      R=rhs
      alpham1=one
      beta=zero
      gamm1=one
      gam=one
      gamp1=one
      dnm=dopr(rhs,rhs)
      do it=1,maxit
        call SdagSpC(P(:,:,1),B,u,DAGGER,dwmass,czero)
        B=B+sig(1)*P(:,:,1)
        drr=dopr(R,R)
        alpha=drr/dopr(P(:,:,1),B)
        Rp1=R-alpha*B
        drp1rp1=dopr(Rp1,Rp1)
        betap1=drp1rp1/drr
        Z(:,:,1)=Z(:,:,1)+alpha*P(:,:,1)
        P(:,:,1)=Rp1+betap1*P(:,:,1)
        Zm1ns=Z(:,:,nsr)
        do o=2,nsr
          gamp1(o)=gam(o)*gamm1(o)*alpham1/
     &      ( alpha*beta*(gamm1(o)-gam(o)) +
     &          gamm1(o)*alpham1*(one+alpha*(sig(o)-sig(1))) )
          Z(:,:,o)=Z(:,:,o)+alpha*gamp1(o)/gam(o)*P(:,:,o)
          P(:,:,o)=gamp1(o)*Rp1+betap1*(gamp1(o)/gam(o))**2*P(:,:,o)
        end do
        beta=betap1
        conv=sqrt(drp1rp1/dnm)
       
        cns=maxval(abs(Z(:,:,nsr)-Zm1ns))/maxval(abs(Z(:,:,nsr)))
        if (cns.lt.1d-18) then
          ocount(nsr)=it
          nsr=nsr-1
          if (nsr.eq.0) then
            goto 501
          endif
        endif
        if (conv.lt.1d-10) then
          goto 501
        endif
        alpham1=alpha
        R=Rp1
        gamm1=gam
        gam=gamp1
      end do

501   continue
      do o=1,ns
        DR=DR+mult(o)*Z(:,:,o)
      end do

      deallocate(gamm1,gam,gamp1,sig,ocount)
      deallocate(Z,P)
      end associate
      return
      end subroutine MSCGS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function dopr(R1,R2) result(val)
      implicit none
      complex(prc) R1(Nv,4),R2(Nv,4)
      real(prc) val
      complex(prc) cval
      integer v,d
     
      cval=czero
      do d=1,4
        do v=1,Nv
          cval=cval+conjg(R1(v,d))*R2(v,d)
        end do
      end do
c      print *,"cval:",cval
      val=real(cval,prc)
      return 
      end function dopr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLpf(RR,S,u,DAGGER,dwmass,SRF)
      use options
!     approximate Voverlap using partial fraction rational functions
!     
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      procedure(),pointer :: Mptr => NULL()
      real(prc) mult

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      nf=nn-nd+1
      S=zero
      if (nf.ge.1) then ! it should only ever be 0 or 1
        do j=0,nf-1
          mult=front(j+1)
!         T = m.DWilson^j.RR
          TMP1=mult*RR
          do p=1,j
            call DdagD(TMP1,TMP2,u,DAGGER,dwmass,czero)
            TMP1=TMP2
          end do
          S=S+TMP1
        end do
      end if

      Mptr => DdagDpC
      do j=1,nd
         add = -denom(j)
         call IM(RR,TMP1,u,DAGGER,dwmass,add,Mptr)
         TMP2 = pf(j)*TMP1
         S=S+TMP2
      end do

      call DWilson(S,TMP1,u,DAGGER,dwmass,czero)
      S=SRF%mult*TMP1

      end associate
      return
      end subroutine VOLpf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLSpf(RR,S,u,DAGGER,dwmass,SRF)
      use options
      use ratfuncs
      use axbmodule1
!     approximate Voverlap using partial fraction rational functions
!     
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      complex(prc) TMP3(Nv,4),TMP4(Nv,4)
      integer j,p,nn,nd,nf
c      procedure(),pointer :: Mptr => NULL()
      real(prc) mult
      complex(prc) cmult

c      print *,"VOLSpf"

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      nf=nn-nd+1
      S=zero
      if (nf.eq.1) then ! it should only ever be 0 or 1
!       b0.RR
        S=front(1)*RR
      end if

c      Mptr => SdagSpC
      do j=1,nd
         add = -denom(j)
c         print *,j,"add:",add
c         call IM(RR,TMP2,u,DAGGER,dwmass,add,Mptr)
c         call IM(RR,TMP4,u,.not.DAGGER,dwmass,add,Mptr)
         call ISdagSpC(RR,TMP1,u,DAGGER,dwmass,add)
c         call ISdagSpC(RR,TMP3,u,.not.DAGGER,dwmass,add)
c         print *,j,"IMdag-ISdag:",maxval(abs(TMP1-TMP2))
c         print *,j,"IM-IS:",maxval(abs(TMP3-TMP4))
c         print *,j,"IMdag-IM:",maxval(abs(TMP2-TMP4))
c         print *,j,"ISdag-IS:",maxval(abs(TMP1-TMP3))
         TMP2 = pf(j)*TMP1
         S=S+TMP2
      end do

      call DShamir(S,TMP1,u,DAGGER,dwmass,czero)
      S=SRF%mult*TMP1

      end associate
      return
      end subroutine VOLSpf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLscaled(RR,S,u,DAGGER,dwmass,SRF,alpha)
      use options
!     approximate scaled Voverlap using partial fraction rational functions
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      real(prc) alpha
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      procedure(),pointer :: Mptr => NULL()
      real(prc) mult

      print *,"VOLscaled"

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      nf=nn-nd+1
      S=zero
      if (nf.eq.1) then ! it should only ever be 0 or 1
        mult=front(1) ! not affected by alpha
        S=mult*RR
      end if

      Mptr => scaledDdagDpC
      do j=1,nd
         add = -denom(j)
         call IMscaled(RR,TMP1,u,DAGGER,dwmass,add,Mptr,alpha)
         TMP2 = pf(j)*TMP1
         S=S+TMP2
      end do

      call DWilson(S,TMP1,u,DAGGER,dwmass,czero)
      S=SRF%mult*alpha*TMP1

      end associate
      return
      end subroutine VOLscaled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module VOLmodule
