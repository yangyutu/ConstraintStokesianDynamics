	subroutine choldc(a,p)
	
	common / np_3 / np3
c
c     a real version of the Numerical Recipes 
c	routine for the Cholesky decomposition
c

c
c	declare variables
c
      integer np3
      double precision a(np3,np3), p(np3)

      integer i, j, k
      double precision sum
c
      do i=1, np3
         do j=1, np3

            sum = a(i,j)
            do k=i-1, 1, -1
               sum = sum - a(i,k)*a(j,k)
            end do
            if(i.eq.j)then
              if(sum.le.0.0d0)then
			  write(*,*)i
			  
c			  stop 'choldc failed'
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
