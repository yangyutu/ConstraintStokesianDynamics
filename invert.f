	SUBROUTINE INVERT(a,np)
c

c	This routine replaces the matrix a(np,np) by its inverse
c	Subroutine taken from Numerical Recipes
c

c
	INTEGER np, i, j
	double precision d, a(np, np)

	integer, allocatable :: indx(:)

	double precision, allocatable :: y(:, :)

	allocate ( indx(np), y(np, np) )
	
	do i=1,np 
	  do j=1,np
	     y(i,j)=0d0
	  enddo
	  y(i,i)=1d0
	enddo

	call ludcmp(a,np,indx,d) 

	do j=1,np 
	  call lubksb(a,np,indx,y(1,j))
	enddo
	do i=1,np
	  do j=1,np
	    a(i,j)=y(i,j)
	  end do
	end do

	deallocate ( indx, y )
	
	return
	END
c 
