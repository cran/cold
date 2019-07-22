      subroutine pssgrd(grad,beta,rho,npar,x,y,theta,work,n,f,link)

      implicit double precision (a-h,o-z)
      DIMENSION x(n,npar-1),beta(npar-1),y(n),theta(n),work(n),
     *grad(npar),f(0:200)
      double precision nu
      integer y,n,p,r,rmax,i,i0,i1,npar,link,n0,k,k0,k1
      external fpss
c      G(r,k)= rm**r * (1-rm)**(k-r)
c      H(r)= expnu * nu**r /f(r)
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
   20 if (y(i0).eq.(-1)) then
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

C  derivatives with respect to beta
      rm=rho**(i1-i0)
      k1=y(i1)
      k0=y(i0)
      prob=fpss(i0,k0,i1,k1,theta,rho,f)
      a=-prob
      if (k1.gt.0) then
      a=a+fpss(i0,k0,i1,k1-1,theta,rho,f)
      end if
      do 70 ip=1,p
      if (link.eq.0) then
      gl0=one
      gl1=one
      else if (link.eq.1) then
      gl0=theta(i0)
      gl1=theta(i1)
      else
      gl0=zero
      gl1=zero
      end if
      b=gl1*x(i1,ip)-rm*gl0*x(i0,ip)
      grad(ip)=grad(ip)+a*b/prob
   70 continue

C  derivatives with respect to rho
      rmax=min(k1,k0)
      nu=theta(i1)-rm*theta(i0)
      expnu=dexp(-nu)
      dprm=zero
      do 80 r=0,rmax
      comb=f(k0)/(f(r)*f(k0-r))
      g1= rm**(r-1) * (1-rm)**(k0-r)
      g2= rm**r * (1-rm)**(k0-r-1)
      g3= rm**r * (1-rm)**(k0-r)
      h1= expnu * nu**(k1-r) /f(k1-r)
      dg=r*g1-(k0-r)*g2
      dh=-h1
      if (k1.gt.r) then 
      h2= expnu * nu**(k1-r-1) /f(k1-r-1)
      dh=dh+h2
      end if
      dh=dh*(-theta(i0))
      dprm=dprm+comb*(dg*h1+g3*dh)
   80 continue
      grad(npar)=grad(npar)+(i1-i0)*rm*dprm/(prob*rho)
      i0=i1  
      i=i0+1

      go to 50
      end if
      return
      end
