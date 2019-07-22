      subroutine pssli(logL,np,n)

      implicit double precision (a-h,o-z)
      DIMENSION x1(4000,10),theta1(4000),
     *work1(4000),y1(4000),
     *beta1(10),bt1(10),fact(0:130)

      double precision logL,fact
      integer y1,i0,i1,np,n,n0,m,mpar,link1,maxy1
      external fpss

      COMMON/param/x1,theta1,work1,y1,
     *beta1,bt1,m,mpar,omega1,rho1,link1,maxy1

      data zero/0.0d0/

      call mati(x1,beta1,work1,4000,10,1,n,np+1)
      call fac(fact,130)
      i0=0
      do 10 i=1,n
      if (link1.eq.0) then 
      theta1(i)=work1(i)
      else if (link1.eq.1) then 
      theta1(i)=dexp(work1(i))
      end if
   10 continue
      i0=1
   20 if (y1(i0).eq.(-1)) then
      i0=i0+1
      go to 20
      end if

      n0 = n
   30 if (y1(n0).eq.(-1)) then
      n0=n0-1
      go to 30
      end if

      logL=0

      logL=-theta1(i0)+y1(i0)*dlog(theta1(i0))-dlog(fact(y1(i0)))
C     logL=-theta1(i0)+y1(i0)*dlog(theta1(i0))

      if (i0.eq.n0) return
      i = i0+1
   40 if (i.le.n0) then
      i1 = i
   50 if (y1(i1).eq.(-1)) then
      i1=i1+1
      go to 50
      end if

C  i0 is the most recent (past) observation time
C  i1 is the next observation time
  
      prob=fpss(i0,y1(i0),i1,y1(i1),theta1,rho1,fact)
      logL = logL+dlog(prob)
      
C      logL = logL+dlog(prob*fact(y1(i1)))
      ! logL = logL+ dlog(prob*fact(y(i1))) 
      ! moltiplico prob per fact(y(i1)) per prevenire -Inf,
      ! inoltre questo sostanzialmente e` simile a togliere il fact() 
      ! a denominatore in una Poisson
      
      i0=i1
      i=i0+1
      go to 40
      end if
      return
      end

