	SUBROUTINE random(CholDiag, D, randisp)

	  common / np_3 / np3
	  common / seed / iDummy

c    ******************************************************************
c    ** subroutine to calculate the random displacements             **
c    ******************************************************************

	  integer j, k, iDummy, np3
	  double precision CholDiag(np3), D(np3,np3), 
     +				   randisp(np3), stdNormDev
	  double precision gasdev

	  do j=1,np3
	   randisp(j) = 0.0d0
	  end do

        do j=1,np3
         stdNormDev = gasdev(iDummy)
c         ...diagonal term ...
         randisp(j) = randisp(j) + stdNormDev*(CholDiag(j))
c         ... other terms ...
          do k=j+1,np3
            randisp(k) = randisp(k) + stdNormDev*(D(k,j))
          end do
	  end do

	  return
	  end
