      subroutine sumF(x, n, s)
      integer n, i
      double precision x(n), s
      s=0
      do 10, i=1, n
        s = s + x(i)
   10 continue
      return
      end
