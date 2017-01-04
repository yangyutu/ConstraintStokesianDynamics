	SUBROUTINE grndrm(rinfi, rmove, rgrnd, t,r, nxyz)

	  common / num_tot / n
	  common / n_3 / n3
	  common / num_par / np
	  common / np_3 / np3
	  
	  common / nw_3 / nw3

	  common / radius / a

	  common / boxlen_xy / boxlen
         common /constraint/ Kconstr,kconstr3,nxyzbin
c    ******************************************************************
c    ** This subroutine sets up the grand resistance matrix as       **
c    ** defined by Brenner et al. It is 3n*3n matrix.                **
c    ******************************************************************

	  integer	n, np, np3, n3, nw3, i, j, k, l, p1, p2, alpha, beta

	  integer   Kconstr,kconstr3,nxyzbin
       integer nxyz(3,nxyzbin)

	  double precision      r(np3), a, rmove(np3), t, ravg, rsq, 
     +						rexact(6, 6), rabinf(6, 6), pw(3, 3), 
     +						pwinf(3, 3)

	  double precision		rx, ry, rz, rsep, boxlen

	  double precision      rinfi(np3, np3), rgrnd(np3, np3), 
     +			 rimage(np3),rinfi2(np3,np3),
     + Gact(np3,np3) 

	  

c
c	calculate average displacement for each particle
c
	  ravg=0.d0
	  
	    do p1=1,np
	  rimage(nxyz(1,p1))=r(nxyz(1,p1))
	  rimage(nxyz(2,p1))=r(nxyz(2,p1))
	  rimage(nxyz(3,p1))=-r(nxyz(3,p1))
	  end do
	  
c	reinitialize resistance tensor to zero
c	  
		rinfi = 0.d0
c
c	compute the far-field interactions by summing up for each wall 
c	"particle"
c
c			call greenact(Gact,r,nxyz)
c			call invert (Gact, np3)
c		call farfldres(rinfi,  r, rimage, nxyz)
		
          call farfldres2(rinfi, r, rimage, nxyz)
c
c	change coordinates for rmove 
c
		do i=1,np3
		  
		  rmove(i)=r(i)

		end do
c
c	  endif
c
c	  
	  do i = 1, np3
	    do j = 1, np3
			
		  rgrnd(i,j) = rinfi(i,j)

		end do
		    
	  end do

c	call invert(rgrnd,np3)
c	do i = 1, np
c	write(*,*)i
c	write(*,*)
c	write(*,*)rgrnd(nxyz(1,i):nxyz(3,i),nxyz(1,i):nxyz(3,i))
c	pause
c	end do
c
c	add particle-particle lubrication forces to the far-field 
c	resistance tensor
c
	
	  do p1 = 1, np-1

	    do p2 = p1+1, np

c    ******************************************************************
c	call subroutines for calculating the 6X6 exact and far-field   **
c	particle-particle mobility tensors. invert them to get the     **
c	exact and far-field particle-particle resistance matrices.     **
c    ******************************************************************
			
		  rx = r( nxyz(1, p1) ) - r( nxyz(1, p2) )
		  ry = r( nxyz(2, p1) ) - r( nxyz(2, p2) )
		  rz = r( nxyz(3, p1) ) - r( nxyz(3, p2) )
		  
c		  rx = rx - anint(rx / boxlen) * boxlen
c		  ry = ry - anint(ry / boxlen) * boxlen
		  
		  rsep = sqrt(rx**2 + ry**2 + rz**2) - 2.0 * a

		  if(rsep .lt. 6.0*a)then 	
				
		    call ppexct(rexact, p1, p2, t, r, nxyz)

		    call ppinf(rabinf, p1, p2, r, nxyz)
c       note: in the ppinf, i have not calculated the particle particle far field interaction
c     because in the pp exct, i have not calculate the wall contribution to the particle particle
c     exact interactions, so i do not need to substract them to avoid double countinng issue
c    ******************************************************************************
c	begin loops over corordinates. i,j=1,3 is for particle 1 i,j=4,6 is for    **
c	particle 2 in the two body mobility tensor								   **
c    ******************************************************************************
		  
		    do i=1,6
		      do j=1,6

			    if(i.gt.3)then
			      k=p2
			      alpha=3
			    else
				  k=p1
				  alpha=0
			    endif

			    if(j.gt.3)then
				  l=p2
				  beta=3
			    else
				  l=p1
				  beta=0
			    endif

			    rgrnd(nxyz(i-alpha,k),nxyz(j-beta,l))=
     +		    rgrnd(nxyz(i-alpha,k),nxyz(j-beta,l))+rexact(i,j)
     +									             -rabinf(i,j)
	  
			  end do
		    end do

		  endif
		
	    end do

	  end do
c
c	call invert(rgrnd,np3)
c	do i = 1, np
c	write(*,*)i
c	write(*,*)
c	write(*,*)rgrnd(nxyz(1,i):nxyz(3,i),nxyz(1,i):nxyz(3,i))
c	pause
c	end do
c
c	add particle-wall lubrication forces. the particle-wall exact solution
c	is calculated, and added to the grand resitance tensor as
c	before. 
c
	 do p1 = 1, np	

	   call pwexact (pw, p1, r, nxyz)
	   call pwff(pwinf,p1,p1,r,rimage,nxyz)
	   
	   do i = 1, 3
		 do j = 1, 3
   
	   
	   			rgrnd( nxyz(i, p1), nxyz(j, p1) ) = 
     +		rgrnd( nxyz(i, p1), nxyz(j, p1) ) + pw(i,j)-pwinf(i,j)
        end do
        end do
      end do
      
c      for the far-field two-body resistnaces, a mobility tensor is
c	constructed for each particle using image method.
c	this is inverted to get the far-field two body resistance tensor, which
c	is added to the grand resistance tensor to get the complete hydrodynamic
c	forces (Nott & Brady, Singh & Nott).   
        
     
       
c	call invert(rgrnd,np3)
c	do i = 1, np
c	write(*,*)i
c	write(*,*)
c	write(*,'(3f7.3)')rgrnd(nxyz(1,i):nxyz(3,i),nxyz(1,i):nxyz(3,i))
c	pause
c	end do	   
c	
	
	return
	end