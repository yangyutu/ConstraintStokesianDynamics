	SUBROUTINE lubksb(a,np,indx,b)
c

c
c	See Numerical Recipes for description on this subroutine
c

c
      INTEGER n,np,indx(np)
	double precision b(np)
	INTEGER i,ii,j,ll
	double precision sum

	double precision a(np, np) 

	n=np
	ii=0 
	do i=1,n
	  ll=indx(i)
	  sum=b(ll)
	  b(ll)=b(i)
	  if (ii.ne.0.0)then
	    do j=ii,i-1
	      sum=sum-a(i,j)*b(j)
		enddo
	  else if (sum.ne.0.0) then
	    ii=i 
	  endif
	  b(i)=sum
	enddo
	do i=n,1,-1 
	  sum=b(i)
	  do j=i+1,n
		sum=sum-a(i,j)*b(j)
	  enddo
	  b(i)=sum/a(i,i) 
	enddo

	return 
	END
c
