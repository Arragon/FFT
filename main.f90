PROGRAM ComputeFFT
    IMPLICIT NONE
    !INTERFACE
    !    SUBROUTINE fft(a,n,inv)
    !        COMPLEX,DIMENSION(:) :: a
    !        INTEGER,INTENT(IN) :: n,inv
    !    ENDSUBROUTINE fft
    !ENDINTERFACE

    INTEGER :: count = 0
    INTEGER :: Maxm = 32768
    INTEGER :: ios,ii,x,inv
    REAL :: dt,df,f,fst,fen,data_read
    REAL,DIMENSION(:),ALLOCATABLE :: data1,data2
    COMPLEX,DIMENSION(:),ALLOCATABLE :: four1
    REAL,PARAMETER    :: C=2.997925e+8,Pi=3.14159265359

    ALLOCATE (data1(0:Maxm),data2(0:Maxm))
    ALLOCATE (four1(0:Maxm))
    OPEN(1,file="D:\Files\qfdtd90_test\qfdtd_ader_forward_difference\qfdtd90_new\es.csv")
    DO
        READ(1,*,IOSTAT=ios) ii,data_read
        !write(*,*) II,data_read
        IF(ios /= 0) EXIT
        data1(count) = data_read
        four1(count) = data1(count)
        count = count + 1
    ENDDO
    CLOSE(1)
    DO x=count,Maxm
        four1(x)=(0.0,0.0)
    ENDDO
    
    do x=1, maxm
      if(four1(X).ne.0) print *,four1(X)
    enddo

    OPEN(10,file='spectra.csv')
    inv=0
    CALL fft(four1,Maxm,inv)
    DO x=0,Maxm
        data2(x) = CABS(four1(x))
    ENDDO
    dt = 1.1824031E-13
    fst = 0
    fen = 10
    df=1./(dt*Maxm)
    write(*,*) df
    x=0
    DO
        f=x*df
        IF(f > (fen*1.e+9)) EXIT
        IF((f >= (fst*1.e+9)).and.(f <= (fen*1.e+9))) THEN
            write(*,*) f
            !pause
            WRITE(10,*)f*1.e-9,",",CABS(four1(x))/MAXVAL(data2)
            print *, x*df, CABS(four1(x)), MAXVAL(data2)
            !pause
        ENDIF
        x=x+1
    ENDDO
    CLOSE(10)
    DEALLOCATE (data1,data2)
    DEALLOCATE (four1)
  contains
  	SUBROUTINE fft(a,n,inv)
      COMPLEX,DIMENSION(:) :: a
      INTEGER :: n,inv,I,j,k,l,m,nv2,nm1,LE,LE1,IP
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
  END PROGRAM ComputeFFT
  



