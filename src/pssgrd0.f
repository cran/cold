      subroutine pssgrd0(grad,beta,npar,x,y,theta,work,n,link)

      implicit double precision (a-h,o-z)
      DIMENSION x(n,npar-1),beta(npar-1),y(n),theta(n),work(n),
     *grad(npar)
      integer y,n,p,i,i0,i1,npar,link,n0,k

      data zero/0.0d0/, one/1.0d0/
      p=npar-1
      call matp(x,beta,work,n,p,1)

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

      do 40 k=1,p
        if (link.eq.0) then 
      gl=one
        else if (link.eq.1) then 
      gl=theta(i0)
        else 
      gl=zero
       end if
      grad(k)=(-one+y(i0)/theta(i0))*gl*x(i0,k)
   40 continue

      grad(npar)=zero
      if (i0.eq.n0) return
        i = i0+1
   50 if (i.le.n0) then
      i1 = i
   60   if (y(i1).eq.(-1)) then
      i1=i1+1
      go to 60
      end if

C  i0 is the most recent (past) observation time
C  i1 is the next observation time

      do 70 ip=1,p
      if (link.eq.0) then
      gl1=one
      else if (link.eq.1) then
      gl1=theta(i1)
      else
      gl1=zero
      end if

      a=(-one+y(i1)/theta(i1))*gl1*x(i1,ip)
      grad(ip)=grad(ip)+a
   70 continue
   

C  derivatives with respect to rho
      grad(npar)=zero

      i0=i1
      i=i0+1

      go to 50
      end if
      return
      end
