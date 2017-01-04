	subroutine choldcorigin(a,n,p)
	
	
c
c     a real version of the Numerical Recipes 
c	routine for the Cholesky decomposition
c

c
c	declare variables
c
      integer np3,n
      double precision a(n,n), p(n)
      
      integer i, j, k
      double precision sum
c
      np3=n
      do i=1, np3
         do j=1, np3

            sum = a(i,j)
            do k=i-1, 1, -1
               sum = sum - a(i,k)*a(j,k)
            end do
            if(i.eq.j)then
              if(sum.le.0.0d0)then
			  write(*,*)i,'choldc failed'
			  
c			  stop 
		    endif
              p(i) = sqrt(sum)
            else
    	        a(j,i) = sum/p(i)
            endif

          end do
      end do
c
      return
      end
c
