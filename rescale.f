      subroutine rescale(r, nxyz,fixdist,distcheck,rescalecheck,
     +  preheight) 		
      		

	  common / num_par / np
	  common / np_3 / np3
	  common / radius / a 		
       common /constraint/ Kconstr,kconstr3,nxyzbin
       
      	  integer, parameter :: m = 2, maxbin=5000  
	  integer    i, j, np, np3, n, n3,mu, nv,kconstr,kconstr3,nxyzbin
	  
	  integer	 nxyz(3, nxyzbin) ,centralindex,rescalecheck
	  double precision r(np3),fixdist,distcheck(maxbin),
     +   direct_vectortemp(3),direct_vectornorm, direct_vector(3), 
     +  center(3), preheight,newcenter(3)
       
        center=0.0
        do i=1,np
        center(1)=center(1)+r(nxyz(1,i))
        center(2)=center(2)+r(nxyz(2,i)) 
        center(3)=center(3)+r(nxyz(3,i)) 
        end do       
        center(1)=center(1)/dble(np)
         center(2)=center(2)/dble(np)
          center(3)=center(3)/dble(np)

		centralindex=int(np/2)+1
		
		r(nxyz(1,centralindex))=center(1)
        r(nxyz(2,centralindex))=center(2)
        r(nxyz(3,centralindex))=center(3)
   
		 rescalecheck=1
		 direct_vectortemp(1)=r(nxyz(1,np))-r(nxyz(1,1))
		 direct_vectortemp(2)=r(nxyz(2,np))-r(nxyz(2,1))   
		 direct_vectortemp(3)=r(nxyz(3,np))-r(nxyz(3,1))  
		 
		 direct_vectornorm=sqrt(direct_vectortemp(1)**2+
     +	direct_vectortemp(2)**2+direct_vectortemp(3)**2)
		 direct_vector(1)=direct_vectortemp(1)/direct_vectornorm
		 direct_vector(2)=direct_vectortemp(2)/direct_vectornorm
        direct_vector(3)=direct_vectortemp(3)/direct_vectornorm
        
        do j=1,np
        
        
        if(j .ne. centralindex) then
        r(nxyz(1,j))=r(nxyz(1,centralindex))+
     +  (j-centralindex)*fixdist*direct_vector(1)
                 r(nxyz(2,j))=r(nxyz(2,centralindex))+
     +  (j-centralindex)*fixdist*direct_vector(2)
                 r(nxyz(3,j))=r(nxyz(3,centralindex))+
     +  (j-centralindex)*fixdist*direct_vector(3)
        
         end if
         
         end do
         
                 newcenter=0.0
        do i=1,np
        newcenter(1)=newcenter(1)+r(nxyz(1,i))
        newcenter(2)=newcenter(2)+r(nxyz(2,i)) 
        newcenter(3)=newcenter(3)+r(nxyz(3,i)) 
        end do       
        newcenter(1)=newcenter(1)/dble(np)
         newcenter(2)=newcenter(2)/dble(np)
          newcenter(3)=newcenter(3)/dble(np)
         
          do j=1,np
        
        

        r(nxyz(1,j))=r(nxyz(1,j))+center(1)-newcenter(1)
        r(nxyz(2,j))=r(nxyz(2,j))+center(2)-newcenter(2)
        r(nxyz(3,j))=r(nxyz(3,j))+center(3)-newcenter(3)


         
         end do
         
         end