	SUBROUTINE farfldres2(rinfi, r,rimage, nxyz)

	  common / num_par / np
	  common / np_3 / np3

	  common / selftm_vol / ms, vol
	  
	   common / delta_fn / delta
	   	  common / radius / a
        common /constraint/ Kconstr,kconstr3,nxyzbin
        common /up_wall/ upperwall
        common / twowall/ twowallflag
c    ******************************************************************
c     This subroutine sets up a 3n*3n resistance tensor and returns 
c     a 3np*3np subset for calculating particle velocities. 
c    ******************************************************************
	
	integer  n, np, n3, np3, nw3, i, j, k, l
        
	integer   Kconstr,kconstr3,nxyzbin
       integer nxyz(3,nxyzbin),p1,p2, twowallflag

	double precision  mr(3, 3), mk(3, 3), ms(3, 3), vol,  
     +				  rinfi(np3, np3), r(np3),rabinf(6,6),
     +   delta(6,6),rij(3),a,dist,stklet,dipole,imagerij(3),
     + imagestklet,imagedipole, rimage(np3),imagedist,wallgreen,
     + tempsum,e(3),h,real_R, norm_R,e_real(3),mabuf(3,3),
     + maauf(3,3) ,hh, upperwall,rimage2(np3)
     


	  rinfi = 0.0
	  
	
	  
	  do p1=1, np
	  
	      do p2=1, np
	      
	      
	       do i=1, 3

	   rij(i) = r(nxyz(i, p1)) - r(nxyz(i, p2))
            imagerij(i)=r(nxyz(i,p1))-rimage(nxyz(i,p2))
	  enddo

	  dist=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
	  imagedist=sqrt(imagerij(1)**2+imagerij(2)**2+imagerij(3)**2)
	      
        if(dist .le. 2*a .and. p1 .ne. p2) then
        dist=2*a+1
        end if
	      
	  h=r(nxyz(3,p2))/imagerij(3)
	  e(1)=imagerij(1)/imagedist
	  e(2)=imagerij(2)/imagedist
	  e(3)=imagerij(3)/imagedist
	  
	  e_real(1)=rij(1)/dist
	  e_real(2)=rij(2)/dist
	  e_real(3)=rij(3)/dist
	  
	  if(p1 .ne. p2) then
	  
	  Norm_R=imagedist/a
	  real_R=dist/a
	  		    do i=1,3
		      do j=1,3

c       first construct the MUF component for real particle p1, p2
c      which is from paper Dynamic simulation of hydrodynamically interacting particles
        mabuF(i,j)=(1.5/real_R-1.0/real_R**3.0)*e_real(i)*e_real(j)+
     +  (0.75/real_R+0.5/real_R**3.0)*(delta(i,j)-e_real(i)*e_real(j)) 
c      then we construct the extra MUF component for the particle p1, p2 above wall
c      which is from paper simulation of hydrodynamically interacting particles near a no-slip wall      
        mabUF(i,j)=mabUF(i,j)
     +   -1.0/4.0*(3.0*(1.0+2.0*h*(1-h)*e(3)**2.0)/Norm_R+
     +2.0*(1.0-3.0*e(3)**2.0)/Norm_R**3.0-2.0*(1.0-5.0*e(3)**2.0)/
     + Norm_R**5.0)*delta(i,j)
     + -1.0/4.0*(3.0*(1-6.0*h*(1-h)*e(3)**2)/Norm_R-
     + 6.0*(1.0-5.0*e(3)**2)/Norm_R**3+
     + 10.0*(1.0-7.0*e(3)**2)/Norm_R**5.0)*e(i)*e(j)
     + +1.0/2.0*e(3)*(3.0*h*(1.0-6.0*(1-h)*e(3)**2)/Norm_R-
     + 6.0*(1.0-5.0*e(3)**2)/Norm_R**3+
     + 10.0*(2.0-7.0*e(3)**2)/Norm_R**5.0)*e(i)*delta(j,3)
     + +0.5*e(3)*(3.0*h/Norm_R-10.0/Norm_R**5.0)*delta(i,3)*e(j)
     + -(3.0*h**2.0*e(3)**2/Norm_R+3.0*e(3)**2.0/Norm_R**3+
     + (2.0-15.0*e(3)**2)/Norm_R**5.0)*delta(i,3)*delta(j,3)
             end do
                  end do
                  
        do i=1,3
        do j=1,3
        rinfi(nxyz(i,p1),nxyz(j,p2))=
     +	rinfi(nxyz(i,p1),nxyz(j,p2))+mabUF(i,j)
          
        end do
        end do
                      
       
        else
             
c       p1=p2
        do i=1,3
        do j=1,3
        maaUF(i,j)=delta(i,j)
        end do
        end do
        
                hh=r(nxyz(3,p1))/a
        
c       now adding the reflection wall terms
        do i=1,3
        do j=1,3
        maaUF(i,j)=maaUF(i,j)-1.0/16.0*(9.0/hh-2.0/hh**3+1.0/hh**5)*
     +   (delta(i,j)-delta(i,3)*delta(j,3))-
     + 1.0/8.0*(9.0/hh-4.0/hh**3+1.0/hh**5)*delta(i,3)*delta(j,3)
        end do
        end do
         
         do i=1,3
        do j=1,3
        rinfi(nxyz(i,p1),nxyz(j,p2))=
     +	rinfi(nxyz(i,p1),nxyz(j,p2))+maaUF(i,j)
         
        end do
        end do   

        end if
        
        end do
        end do
        
        
        if(twowallflag .eq. 1) then
            
            maaUF=0.0
            mabUF=0.0
         	    do p1=1,np
	  rimage2(nxyz(1,p1))=r(nxyz(1,p1))
	  rimage2(nxyz(2,p1))=r(nxyz(2,p1))
	  rimage2(nxyz(3,p1))=2.0*upperwall-r(nxyz(3,p1))
              end do              
         	  do p1=1, np	  
	      do p2=1, np	      
	       do i=1, 3

	   rij(i) = r(nxyz(i, p1)) - r(nxyz(i, p2))
            imagerij(i)=r(nxyz(i,p1))-rimage2(nxyz(i,p2))
	  enddo

	  dist=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
	  imagedist=sqrt(imagerij(1)**2+imagerij(2)**2+imagerij(3)**2)
	      
        if(dist .le. 2*a .and. p1 .ne. p2) then
        dist=2*a+1
        end if
c note here that h is shifted, and e(3) is changing the sign     
	  h=(r(nxyz(3,p2))-upperwall)/imagerij(3)
	  e(1)=imagerij(1)/imagedist
	  e(2)=imagerij(2)/imagedist
	  e(3)=imagerij(3)/imagedist
	  
	  e_real(1)=rij(1)/dist
	  e_real(2)=rij(2)/dist
	  e_real(3)=rij(3)/dist
	  
	  if(p1 .ne. p2) then
	  
	  Norm_R=imagedist/a
	  real_R=dist/a
	  		    do i=1,3
		      do j=1,3

c       first construct the MUF component for real particle p1, p2
c      which is from paper Dynamic simulation of hydrodynamically interacting particles
c        mabuF(i,j)=(1.5/real_R-1.0/real_R**3.0)*e_real(i)*e_real(j)+
c     +  (0.75/real_R+0.5/real_R**3.0)*(delta(i,j)-e_real(i)*e_real(j)) 
c      then we construct the extra MUF component for the particle p1, p2 above wall
c      which is from paper simulation of hydrodynamically interacting particles near a no-slip wall      
        mabUF(i,j)=
     +   -1.0/4.0*(3.0*(1.0+2.0*h*(1-h)*e(3)**2.0)/Norm_R+
     +2.0*(1.0-3.0*e(3)**2.0)/Norm_R**3.0-2.0*(1.0-5.0*e(3)**2.0)/
     + Norm_R**5.0)*delta(i,j)
     + -1.0/4.0*(3.0*(1-6.0*h*(1-h)*e(3)**2)/Norm_R-
     + 6.0*(1.0-5.0*e(3)**2)/Norm_R**3+
     + 10.0*(1.0-7.0*e(3)**2)/Norm_R**5.0)*e(i)*e(j)
     + +1.0/2.0*e(3)*(3.0*h*(1.0-6.0*(1-h)*e(3)**2)/Norm_R-
     + 6.0*(1.0-5.0*e(3)**2)/Norm_R**3+
     + 10.0*(2.0-7.0*e(3)**2)/Norm_R**5.0)*e(i)*delta(j,3)
     + +0.5*e(3)*(3.0*h/Norm_R-10.0/Norm_R**5.0)*delta(i,3)*e(j)
     + -(3.0*h**2.0*e(3)**2/Norm_R+3.0*e(3)**2.0/Norm_R**3+
     + (2.0-15.0*e(3)**2)/Norm_R**5.0)*delta(i,3)*delta(j,3)
             end do
                  end do
                  
        do i=1,3
        do j=1,3
        rinfi(nxyz(i,p1),nxyz(j,p2))=
     +	rinfi(nxyz(i,p1),nxyz(j,p2))+mabUF(i,j)
          
        end do
        end do
                      
       
        else
             
c       p1=p2
c        do i=1,3
c        do j=1,3
c        maaUF(i,j)=delta(i,j)
c        end do
c        end do
        
                hh=(upperwall-r(nxyz(3,p1)))/a
        
c       now adding the reflection wall terms
        do i=1,3
        do j=1,3
        maaUF(i,j)=-1.0/16.0*(9.0/hh-2.0/hh**3+1.0/hh**5)*
     +   (delta(i,j)-delta(i,3)*delta(j,3))-
     + 1.0/8.0*(9.0/hh-4.0/hh**3+1.0/hh**5)*delta(i,3)*delta(j,3)
        end do
        end do
         
         do i=1,3
        do j=1,3
        rinfi(nxyz(i,p1),nxyz(j,p2))=
     +	rinfi(nxyz(i,p1),nxyz(j,p2))+maaUF(i,j)
         
        end do
        end do   

        end if
        
        end do
        end do     
        
      end if
        
  
c       here 

	  
        
c                open(131,file='mobility.txt')
c        open(132,file='constrain.txt')
c        do i=1,np3
c        do j=1,np3	
c        write(131,*) rinfi(i,j)
c        end do
c        end do
c        close(131)
	  call invert (rinfi, np3)

	  return
	  end
	