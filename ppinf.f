	SUBROUTINE ppinf(mabinf, n1, n2, r, nxyz)

        common / num_par / np
	  common / np_3 / np3

	  common / boxlen / boxlen

	  common / radius / a
	  common / delta_fn / delta
        common /constraint/ Kconstr,kconstr3,nxyzbin
c    *****************************************************************************
c    ** subroutine to calculate the far-field particle-particle mobility tensor **
c    ** by the multipole expansion of the Green's function                      **
c    *****************************************************************************

	  integer	n1, n2, i, j, alpha, beta, k, l, np, np3,
     +  Kconstr,kconstr3,nxyzbin

	  integer	nxyz(3, nxyzbin)

	  double precision		r(np3), a, mabinf(6,6), dist, delta(6,6), 
     +						stklet, dipole, rij(3), boxlen, boxlenz			   

c    **	calculate particle-particle separation   **
	  
	  mabinf = 0.0

		

	  do i=1, 3

	   rij(i) = r(nxyz(i, n1)) - r(nxyz(i, n2))
c	   rij(i) = rij(i) - anint (rij(i)/boxlen) * boxlen

	  enddo

	  dist=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)

c    ** run loops over particles and particle coordinates. i,j denote particle **
c    ** indices. alpha and beta denote particle coordinates                    **

	  do i=1,2
	    if(i.eq.1)then
	      k=0
	    else
	      k=3
	    endif
	    do j=1,2
		  if(j.eq.1)then
	        l=0
		  else
		    l=3
		  endif

	      do alpha=1,3
	        do beta=1,3
		      if(i.ne.j)then
			  			  					    
			    stklet=delta(alpha,beta)/dist+rij(alpha)*rij(beta)
     +                      /dist**3
			    dipole=delta(alpha,beta)/dist**3-3*rij(alpha)
     +        *rij(beta)/dist**5

			  else
			    stklet=0.0
			    dipole=0.0

			  endif


			mabinf(k+alpha,l+beta)=delta(i,j)*delta(alpha,beta)
     +		+3.0*a/4.0*(stklet+dipole*2.0/3.0*a**2)
		if( i .gt. 10) then
		write(*,*) 'h'
		end if
		    end do
		  end do
	    end do
	  end do

		

		
		call invert (mabinf, 6)

	  return
	  end
