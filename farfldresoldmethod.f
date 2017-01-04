	SUBROUTINE farfldres(rinfi, r,rimage, nxyz)

	  common / num_par / np
	  common / np_3 / np3

	  common / selftm_vol / ms, vol
	  
	   common / delta_fn / delta
	   	  common / radius / a
        common /constraint/ Kconstr,kconstr3,nxyzbin
c    ******************************************************************
c     This subroutine sets up a 3n*3n resistance tensor and returns 
c     a 3np*3np subset for calculating particle velocities. 
c    ******************************************************************
	
	integer  n, np, n3, np3, nw3, i, j, k, l
        
	integer   Kconstr,kconstr3,nxyzbin
       integer nxyz(3,nxyzbin),p1,p2

	double precision  mr(3, 3), mk(3, 3), ms(3, 3), vol,  
     +				  rinfi(np3, np3), r(np3),rabinf(6,6),
     +   delta(6,6),rij(3),a,dist,stklet,dipole,imagerij(3),
     + imagestklet,imagedipole, rimage(np3),imagedist,wallgreen,
     + tempsum
     


	  rinfi = 0.0
	  
	
	  
	  do p1=1, np
	  
	      do p2=1, np
	      
	      
	       do i=1, 3

	   rij(i) = r(nxyz(i, p1)) - r(nxyz(i, p2))
c	   rij(i) = rij(i) - anint (rij(i)/boxlen) * boxlen
            imagerij(i)=r(nxyz(i,p1))-rimage(nxyz(i,p2))
	  enddo

	  dist=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
	  imagedist=sqrt(imagerij(1)**2+imagerij(2)**2+imagerij(3)**2)
	      
        if(dist .le. 2*a .and. p1 .ne. p2) then
        dist=2*a+1
        end if
	      
	  		    do i=1,3
		      do j=1,3

            
               stklet=delta(i,j)/dist+rij(i)*rij(j)
     +                      /dist**3
			    dipole=delta(i,j)/dist**3-3*rij(i)
     +        *rij(j)/dist**5
     
	      imagestklet=delta(i,j)/imagedist+
     +      imagerij(i)*imagerij(j)
     +                      /imagedist**3
     
             imagedipole=delta(i,j)/imagedist**3-3*imagerij(i)
     +        *imagerij(j)/imagedist**5
     
            
     
	      if (p1 .eq. p2) then
	      
	      
	      
c	      rinfi(nxyz(i,p1),nxyz(j,p2))=
c     +		    rinfi(nxyz(i,p1),nxyz(j,p2))+delta(i,j)
c     +      -3.0*a/4.0*(imagestklet+imagedipole*2.0/3.0*a**2)
     
            tempsum=-3.0*a/4.0*(imagestklet+imagedipole*2.0/3.0*a**2)
     
           call calwalgnsub(i,j,r,rimage,nxyz,p1,p2,
     + imagerij,imagedist,wallgreen)
            
c             write(*,*) tempsum, wallgreen
             
             rinfi(nxyz(i,p1),nxyz(j,p2))=
     +		    rinfi(nxyz(i,p1),nxyz(j,p2))+delta(i,j)
     +      +wallgreen
           
	      
	      else
	 
			   
c			    rinfi(nxyz(i,p1),nxyz(j,p2))=
c     +		    rinfi(nxyz(i,p1),nxyz(j,p2))
c     +		+3.0*a/4.0*(stklet+dipole*2.0/3.0*a**2)
c     +      -3.0*a/4.0*(imagestklet+imagedipole*2.0/3.0*a**2)
     
          tempsum=-3.0*a/4.0*(imagestklet+imagedipole*2.0/3.0*a**2)
          
             call calwalgnsub(i,j,r,rimage,nxyz,p1,p2,
     + imagerij,imagedist,wallgreen)
     
c        write(*,*) tempsum, wallgreen
        
        
               rinfi(nxyz(i,p1),nxyz(j,p2))=
     +		    rinfi(nxyz(i,p1),nxyz(j,p2))
     +		+3.0*a/4.0*(stklet+dipole*2.0/3.0*a**2)
     +      +wallgreen
	  
	  
	      end if
			  end do
		    end do

	    end do
	    end do
	    
	  
	  call invert (rinfi, np3)

	  

	  return
	  end
	