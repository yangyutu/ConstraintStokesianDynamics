	SUBROUTINE pwff(pwinf, p1, p2, r,rimage, nxyz)

        common / num_tot / n
	  common / n_3 / n3
	  common / num_par / np
	  common / num_wall / nw
	  common / nw_3 / nw3
         common / np_3 / np3
	  common / radius / a
	  common / delta_fn / delta
	  	  
	  common / boxlen_z / boxlenz
	  common / boxlen_xy / boxlen
	   common /constraint/ Kconstr,kconstr3,nxyzbin
         common /up_wall/ upperwall
        common / twowall/ twowallflag
        
c
c
	integer		np, nw, nw3, i, j, k, l, n, n3, pw, n1,p1,p2

	integer   Kconstr,kconstr3,nxyzbin,twowallflag
       integer nxyz(3,nxyzbin)

	double precision  a, r(np3), rijd, delta(6,6), stklet, dipole, 
     +				  rij(3), pwinf(3, 3), boxlenz, boxlen, 
     +				  mw2b(nw3, nw3),rimage(np3),imagerij(3), imagedist,
     +  wallgreen,upperwall,hh,pwinf2(3,3)

c	double precision, allocatable :: mpw(:, :)
	
c	allocate ( mpw(nw3+3, nw3+3) )
c
c	the particle-wall far-field resistances are added in these loops.
c	the procedure is as follows: a mobility tensor is calculated
c	for the free particle and all the wall particles. the mobility
c	tensor is inverted to get the resistance tensor. for calculating
c	particle velocities, only a 3x3 tensor for the single particle
c	is required
c
	pwinf = 0.0
        hh=r(nxyz(3,p1))/a        
	  do i = 1, 3
	    do j = 1, 3      
c       now adding the reflection wall terms
        pwinf(i,j)=pwinf(i,j)+delta(i,j)
     +   -1.0/16.0*(9.0/hh-2.0/hh**3+1.0/hh**5)*
     +   (delta(i,j)-delta(i,3)*delta(j,3))-
     + 1.0/8.0*(9.0/hh-4.0/hh**3+1.0/hh**5)*delta(i,3)*delta(j,3)
		 end do
          end do
c	end do
c
c
          call invert (pwinf, 3)
          
      if(twowallflag .eq. 1) then
          pwinf2=0.0
      hh=(upperwall-r(nxyz(3,p1)))/a        
	  do i = 1, 3
	    do j = 1, 3     
c       now adding the reflection wall terms

        pwinf2(i,j)=pwinf2(i,j)+delta(i,j)
     +   -1.0/16.0*(9.0/hh-2.0/hh**3+1.0/hh**5)*
     +   (delta(i,j)-delta(i,3)*delta(j,3))-
     + 1.0/8.0*(9.0/hh-4.0/hh**3+1.0/hh**5)*delta(i,3)*delta(j,3)

		 end do
          end do
       call invert (pwinf2, 3) 
       pwinf=pwinf+pwinf2
          
        end if



c	call invert (pwinf, 3)
c
c	
	
c
c	deallocate ( mpw )
	return
	end
