c     **** Subroutine for Fast Fourier Transform (FFT) ****
	SUBROUTINE fft(a,n,inv)
      COMPLEX,DIMENSION(:) :: a
      INTEGER :: n,inv
      COMPLEX u,w,t
      REAL,PARAMETER :: Pi=3.14159265359
      m=INT(log(1.*n)/log(2.))
      nv2=n/2
      nm1=n-1
      j=1
      DO 7 i=1,nm1
      IF(i.ge.j) go to 5
      t=a(j)
      a(j)=a(i)
      a(i)=t
 5    k=nv2
 6    IF(k.ge.j) go to 7
      j=j-k
      k=k/2
      go to 6
 7    j=j+k
      DO 20 l=1,m
      le=2**l
      le1=le/2
      u=CMPLX(1.0,0.0)
      IF(inv.eq.0) THEN
           w=CMPLX(COS(Pi/le1),-SIN(Pi/le1))
      else
           w=CMPLX(COS(Pi/le1),SIN(Pi/le1))
      ENDIF
      DO 20 j=1,le1
      DO 10 i=j,n,le
      ip=i+le1
      t=a(ip)*u
      a(ip)=a(i)-t
 10   a(i)=a(i)+t
 20   u=u*w
      IF(inv.eq.0) go to 30
      DO 40 i=1,n
 40   a(i)=a(i)/n
 30   ENDSUBROUTINE fft