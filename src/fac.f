

      subroutine fac(fact,m)
C      ! computes factorial
      implicit double precision (a-h,o-z)
      dimension fact(m+1)   
      fact(1)=1
      do 10 i=2,(m+1)
      fact(i)=fact(i-1)*(i-1)
   10 continue 
      return
      end
