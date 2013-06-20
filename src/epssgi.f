      subroutine pssgi(grad,npar,n)

      implicit double precision (a-h,o-z)
      DIMENSION x1(4000,10),theta1(4000),
     *work1(4000),y1(4000),
     *beta1(10),bt1(10),grad(11),fact1(0:130)
      double precision nu
      integer y1,i0,i1,npar,n,n0,k,m,mpar,link1,p,r,rmax,
     *maxy1
      external fpss

      COMMON/param/x1,theta1,work1,y1,
     *beta1,bt1,m,mpar,omega1,rho1,link1,maxy1

      G(r,k)=rm**r*(1-rm)**(k-r)
      H(r)=expnu*nu**r/fact1(r)
      data zero/0.0d0/, one/1.0d0/

      p=npar-1
      call mati(x1,beta1,work1,4000,10,1,n,npar+1)
      call fac(fact1,130)
      do 10 i=1,n
        if (link1.eq.0) then 
	theta1(i)=work1(i)
        else if (link1.eq.1) then 
	theta1(i)=dexp(work1(i))
       end if
   10 continue
      i0=1
   20   if (y1(i0).eq.(-1)) then
      i0=i0+1
      go to 20
      end if	

      n0 = n
   30   if (y1(n0).eq.(-1)) then
      n0=n0-1
      go to 30
      end if

      do 40 k=1,p
        if (link1.eq.0) then 
	gl=one
        else if (link1.eq.1) then 
	gl=theta1(i0)
        else 
	gl=zero
       end if
      grad(k)=(-one+y1(i0)/theta1(i0))*gl*x1(i0,k)
   40 continue

      grad(npar)=zero
      if (i0.eq.n0) return
        i = i0+1
   50 if (i.le.n0) then
      i1 = i
   60   if (y1(i1).eq.(-1)) then
      i1=i1+1
      go to 60
      end if

C  i0 is the most recent (past) observation time
C  i1 is the next observation time

C  derivatives with respect to beta
      rm=rho1**(i1-i0)
      k1=y1(i1)
      k0=y1(i0)
      prob=fpss(i0,k0,i1,k1,theta1,rho1,fact1)
      a=-prob
      if (k1.gt.0) then
	 a=a+fpss(i0,k0,i1,k1-1,theta1,rho1,fact1)
      end if
      do 70 ip=1,p
      if (link1.eq.0) then
	gl0=one
	gl1=one
      else if (link1.eq.1) then
	gl0=theta1(i0)
	gl1=theta1(i1)
      else
	gl=zero
      end if
      b=gl1*x1(i1,ip)-rm*gl0*x1(i0,ip)
      grad(ip)=grad(ip)+a*b/prob
   70 continue

C  derivatives with respect to rho
      rmax=min(k1,k0)
      nu=theta1(i1)-rm*theta1(i0)
      expnu=dexp(-nu)
      dprm=zero
      do 80 r=0,rmax
	comb=fact1(k0)/(fact1(r)*fact1(k0-r))
	dg=r*G(r-1,k0-1)-(k0-r)*G(r,k0-1)
	dh=-H(k1-r)
	if (k1.gt.r) then 
		dh=dh+H(k1-r-1)
	end if
	dh=dh*(-theta1(i0))
	dprm=dprm+comb*(dg*H(k1-r)+G(r,k0)*dh)
   80 continue
      grad(npar)=grad(npar)+(i1-i0)*rm*dprm/(prob*rho1)
      i0=i1
      i=i0+1

      go to 50
      end if
      return
      end
