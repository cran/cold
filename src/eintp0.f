      subroutine intp0(logL,bt2,beta2,omega,npar,link,
     *maxy,x2,y2,theta2,work2,n,li,ls,epsabs,epsrel,key,limit)
      implicit double precision (a-h,o-z)    
      External fp0

      DOUBLE PRECISION logL,li,ls
      INTEGER n,npar,y2,key,limit,neval,ier,iord,last,
     *m,mpar,y1,k,k2,link,link1,maxy,maxy1

      DIMENSION x1(4000,10),theta1(4000),
     *work1(4000),y1(4000),beta1(10),bt1(10),
     *x2(n,npar-1),beta2(npar-1),theta2(n),
     *y2(n),work2(n),bt2(npar-1),
     *alist(limit),blist(limit),elist(limit),iord(limit),
     *rlist(limit)


      COMMON/param/x1,theta1,work1,y1,
     *beta1,bt1,m,mpar,omega1,link1,maxy1


      do 10 k=1,(npar-1)
      bt1(k)=bt2(k)
      beta1(k)=beta2(k)
   10 continue 
      do 30 k2=1,n 
          do 40 k=1,(npar-1)
           x1(k2,k)=x2(k2,k)
   40 continue 
      y1(k2)=y2(k2)
      theta1(k2)=theta2(k2)
      work1(k2)=work2(k2)
   30 continue 
      m=n
      mpar=npar 
      omega1=omega
      link1=link
      maxy1=maxy
      a=li*dexp(omega/(2.0d0))
      b=ls*dexp(omega/(2.0d0))

      CALL dqager(fp0,a,b,epsabs,epsrel,key,limit,result,abserr,
     *neval,ier,alist,blist,rlist,elist,iord,last,1)
       logL=result
      RETURN
      END
      
