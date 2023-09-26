      module mathtoolsmod
      use pacc
      contains
C **********************************************************************
      SUBROUTINE TRIDIAG (A,NP,D,E,EP,P)
c***********************************************************************
	implicit none
	integer I,J,K,L,NP
c adaptation to hermitian case of minimal storage hoseholder in
c numerical recipes. Nov. 22. 1995.
	complex*16 A(NP,NP),E(NP),P(NP),F,G
	real*8 D(NP), H, SCALE, HH, EP(NP),FA
      IF(NP.GT.1)THEN
        DO 18 I=NP,2,-1
          L=I-1 
          H=0.d0
          SCALE=0.d0
          IF(L.GT.1)THEN
            DO 11 K=1,L 
              SCALE=SCALE+dabs(dreal(A(I,K)))+dabs(dimag(A(I,K)))
11          CONTINUE
            IF(SCALE.eq.0.d0)THEN ! skip transformation
              E(I)=A(I,L) 
            ELSE
              DO 12 K=1,L 
                A(I,K)=A(I,K)/SCALE ! use scaled A's for transformation
                H=H+(dreal(A(I,K)))**2+(dimag(A(I,K)))**2 ! form numrec 
c  sigma in H
12            CONTINUE
              F=A(I,L)
	      FA=cdabs(A(I,L))
		if(FA.ne.0.d0) then
              G=-F*(dsqrt(H)/FA) ! G=-e^iphase(F) * sqrt(sigma)
	         else
	      G=-dsqrt(H)
	         end if
              E(I)=SCALE*G ! stores correct off diagonal value 
c in final location
              H=H-dreal(dconjg(F)*G) ! makes H into 1/2 |u|^2
              A(I,L)=F-G ! makes A(I,L) into last entry in u^\dagger:
c As a result, now the row A(I,j) j=1,2,...L, has u^\dagger in it.
              F=(0.d0,0.d0) ! To accumulate parameter K= u^\dagger p
              DO 15 J=1,L ! J.le.L=I-1
                A(J,I)=dconjg(A(I,J))/H ! stores u/H in column I.
                G=(0.d0,0.d0) ! form the J-th element of u^\dagger A in G
                DO 13 K=1,J ! K.le.J.le.L=I-1
                  G=G+A(I,K)*dconjg(A(J,K)) ! J.le.K.lt.I so rows have not been 
c altered yet.
13              CONTINUE
                IF(L.GT.J)THEN ! protects next do-loop
                  DO 14 K=J+1,L 
                    G=G+A(I,K)*A(K,J) ! continue bulding J-th element of
c u^\dagger A in G. Now J.lt.K and even columns haven't been altered yet.
14                CONTINUE
                ENDIF 
                E(J)=G/H ! stores p^\dagger =u^\dagger A / H in E for J=1 to L.
                F=F+E(J)*dconjg(A(I,J)) ! p^\dagger u.
15            CONTINUE
              HH=dreal(F)/(H+H) ! parameter K in numrec.: p^\dagger u /2H
              DO 17 J=1,L ! Forms q^\dagger= p^\dagger - K u^\dagger in E
                F=A(I,J)
                G=E(J)-HH*F 
                E(J)=G
                DO 16 K=1,J 
                  A(J,K)=A(J,K)-dconjg(F)*E(K)-dconjg(G)*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF 
          ELSE
            E(I)=A(I,L) 
          ENDIF 
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.d0 
      E(1)=(0.d0,0.d0) 
      DO 23 I=1,NP 
        L=I-1 
        IF(D(I).ne.0.d0)THEN
          DO 21 J=1,L 
            G=(0.d0,0.d0)
            DO 19 K=1,L 
              G=G+A(I,K)*A(K,J) ! Forms u^\Dagger * Q
19          CONTINUE
            DO 20 K=1,L 
              A(K,J)=A(K,J)-G*A(K,I) ! Puts in the factor u/H 
20          CONTINUE
21        CONTINUE
        ENDIF 
        D(I)=dreal(A(I,I)) 
        A(I,I)=(1.d0, 0.d0)
        IF(L.GE.1)THEN
          DO 22 J=1,L 
            A(I,J)=(0.d0,0.d0) 
            A(J,I)=(0.d0,0.d0) 
22        CONTINUE
        ENDIF 
23    CONTINUE
c construct diagonal matrices needed to fix accumulated unitary
c transformations so that triangular matrix be real. 
	P(1)=(1.d0,0.d0)
	DO 100 I=2,NP
	EP(I)=cdabs(E(I))! to make triangular matrix real
	if(EP(I).eq.0.d0) then
	P(I)=(1.d0,0.d0)
	else
	P(I)=EP(I)*P(I-1)/E(I)
	end if
100	CONTINUE
c done
c fix accumulated unitary trnasformations
	DO 102 I=1,NP
	F=dconjg(P(I))
	DO 101 J=1,NP
	A(J,I)=A(J,I)*F
101	CONTINUE
102	CONTINUE
c done
      RETURN
      END 

C****************************************************************
      SUBROUTINE TQLI(D,E,NP,Z) 
C****************************************************************
      implicit none
      integer I, K, L, M , NP, ITER
      real*8 D(NP),E(NP),G,R,S,P,B,C,F
      real*8 DD
      complex*16 Z(NP,NP), FC
      IF (NP.GT.1) THEN
        DO 11 I=2,NP 
          E(I-1)=E(I) 
11      CONTINUE
        E(NP)=0.d0
        DO 15 L=1,NP 
          ITER=0
1	continue
	if(L.lt.NP) then
         DO 12 M=L,NP-1 
             DD=DABS(D(M))+DABS(D(M+1))
             IF (DABS(E(M))+DD.EQ.DD) GO TO 2 
12        CONTINUE
	end if
          M=NP 
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30) write(99,*) 'TOO MANY ITERATIONS' 
            ITER=ITER+1 
            G=(D(L+1)-D(L))/(2.d0*E(L)) 
            R=DSQRT(G**2+1.d0) 
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.d0
            C=1.d0
            P=0.d0
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(DABS(F).GE.DABS(G))THEN
                C=G/F 
                R=SQRT(C**2+1.d0) 
                E(I+1)=F*R
                S=1.d0/R
                C=C*S 
              ELSE
                S=F/G 
                R=SQRT(S**2+1.d0) 
                E(I+1)=G*R
                C=1.d0/R
                S=S*C 
              ENDIF 
              G=D(I+1)-P
              R=(D(I)-G)*S+2.d0*C*B 
              P=S*R 
              D(I+1)=G+P
              G=C*R-B 
              DO 13 K=1,NP 
                FC=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*FC 
                Z(K,I)=C*Z(K,I)-S*FC 
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P 
            E(L)=G
            E(M)=0.d0 
            GO TO 1 
          ENDIF 
15      CONTINUE
      ENDIF 
      RETURN
      END 


C****************************************************************
      SUBROUTINE TQLI1(D,E,NP,N,Z) 
C****************************************************************
      implicit none
      integer I, K, L, M , NP, ITER, N
      real*8 D(N),E(N),G,R,S,P,B,C,F
      real*8 DD
      real*8 Z(N,N), FC
      IF (NP.GT.1) THEN
        DO 11 I=2,NP 
          E(I-1)=E(I) 
11      CONTINUE
        E(NP)=0.d0
        DO 15 L=1,NP 
          ITER=0
1	continue
	if(L.lt.NP) then
         DO 12 M=L,NP-1 
             DD=DABS(D(M))+DABS(D(M+1))
             IF (DABS(E(M))+DD.EQ.DD) GO TO 2 
12        CONTINUE
	end if
          M=NP 
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30) write(99,*) 'TOO MANY ITERATIONS' 
            ITER=ITER+1 
            G=(D(L+1)-D(L))/(2.d0*E(L)) 
            R=DSQRT(G**2+1.d0) 
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.d0
            C=1.d0
            P=0.d0
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(DABS(F).GE.DABS(G))THEN
                C=G/F 
                R=SQRT(C**2+1.d0) 
                E(I+1)=F*R
                S=1.d0/R
                C=C*S 
              ELSE
                S=F/G 
                R=SQRT(S**2+1.d0) 
                E(I+1)=G*R
                C=1.d0/R
                S=S*C 
              ENDIF 
              G=D(I+1)-P
              R=(D(I)-G)*S+2.d0*C*B 
              P=S*R 
              D(I+1)=G+P
              G=C*R-B 
              DO 13 K=1,NP 
                FC=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*FC 
                Z(K,I)=C*Z(K,I)-S*FC 
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P 
            E(L)=G
            E(M)=0.d0 
            GO TO 1 
          ENDIF 
15      CONTINUE
      ENDIF 
      RETURN
      END 

      end module mathtoolsmod
