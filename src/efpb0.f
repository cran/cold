
      FUNCTION fpb0(v,i)
      DOUBLE PRECISION logL,fpb0,z,v,
     *bt1,beta1,x1,theta1,work1,omega1,
     *d1,grad
      INTEGER m,mpar,y1,i,maxy1,link1

      DIMENSION x1(4000,10),theta1(4000),
     *work1(4000),y1(4000),
     *beta1(10),bt1(10),grad(10)

      COMMON/param/x1,theta1,work1,y1,
     *beta1,bt1,m,mpar,omega1,link1,maxy1

      beta1(1) = v +bt1(1)

      CALL pssli0(logL,mpar,m)
      z= logL
      CALL pssgi0(grad,mpar,m)
      d1=grad(i)
      fpb0= dexp(z-(v**2)/(2*dexp(omega1)))*d1

      return
      end
