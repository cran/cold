C CCCC
C  function to calculate the probability  for poisson model 
CCCCC

      FUNCTION fpss(t0,k0,t1,k1,theta,rho,f)
      implicit double precision (a-h,o-z)
      DOUBLE PRECISION nu, fpss
      INTEGER t0,t1,r,k1,k0, maxr
      DIMENSION theta(t1), f(0:130)
C      DIMENSION theta(t1), f(0:200)
      
      nu= theta(t1)-rho*theta(t0)
C     if (nu.le.0.0) return
      if(nu.le.0.0) then
      call dblepr('Parameters values unfeasible',28,nu,1)
c      return  !coloquei comentário  a 11-07-19 p ver se não dá erro
      end if

      maxr=min(k1,k0)
      prob=0.0d0
      rhom=rho**(t1-t0)
      do 10 r =0, maxr
         a=f(k0)*rhom**r*(1-rhom)**(k0-r)*dexp(-nu)*nu**(k1-r)
         b=f(r)*f(k0-r)*f(k1-r)
      prob=prob+a/b
   10 continue
      fpss=prob
      return
      end
