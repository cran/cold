
      FUNCTION fp(v)
      DOUBLE PRECISION logL,fp,z,v,
     *bt1,beta1,rho1,x1,theta1,work1,
     *omega1
      INTEGER m,mpar,y1,maxy1,link1

      DIMENSION x1(4000,10),theta1(4000),
     *work1(4000),y1(4000),beta1(10),bt1(10)

      COMMON/param1/x1,theta1,work1,y1,
     *beta1,bt1,m,mpar,omega1,rho1,link1,maxy1

      beta1(1) = v +bt1(1)
  
      CALL pssli(logL,mpar,m)
      z= logL
      fp= dexp(z-(v**2)/(2*dexp(omega1)))
      return
      end
