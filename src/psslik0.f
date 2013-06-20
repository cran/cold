      subroutine psslik0(logL,beta,np,x,y,theta,work,n,fact,link)

      implicit double precision (a-h,o-z)
      DIMENSION x(n,np-1),beta(np-1),y(n),theta(n),work(n), 
     *fact(0:200)

      double precision logL
      integer y,n,link,np,i,i0,i1,n0

      data zero/0.0d0/

      call matp(x,beta,work,n,np-1,1)
      i0=0
      do 10 i=1,n
        if (link.eq.0) then 
	theta(i)=work(i)
        else if (link.eq.1) then 
	theta(i)=dexp(work(i))
       end if
   10 continue
      i0=1
   20   if (y(i0).eq.(-1)) then
      i0=i0+1
      go to 20
      end if	

      n0 = n
   30   if (y(n0).eq.(-1)) then
      n0=n0-1
      go to 30
      end if

      logL=0
      logL=-theta(i0)+y(i0)*dlog(theta(i0))-dlog(fact(y(i0)))

      if (i0.eq.n0) return
      i = i0+1
   40 if (i.le.n0) then
      i1 = i
   50   if (y(i1).eq.(-1)) then
      i1=i1+1
      go to 50
      end if

C  i0 is the most recent (past) observation time
C  i1 is the next observation time
  
      prob=-theta(i1)+y(i1)*dlog(theta(i1))
      logL = logL+prob-dlog(fact(y(i1)))

      i0=i1
      i=i0+1
      go to 40
      end if
      return
      end

