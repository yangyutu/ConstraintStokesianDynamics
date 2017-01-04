        subroutine calrodop(r,nxyz,center,ori_vector,
     +    angle_azi,angle_polar)
        
        common / num_par / np
	  common / np_3 / np3
        common /constraint/ Kconstr,kconstr3,nxyzbin
        
        
         integer i,j,np,np3,Kconstr,kconstr3,nxyzbin
       integer nxyz(3,nxyzbin)
       
       double precision  r(np3), center(3), ori_vector(3),
     +  direct_vectortemp(3), direct_vectornorm ,angle_azi,
     + angle_polar 
       
       
       
       direct_vectortemp(1)=r(nxyz(1,np))-r(nxyz(1,1))
		 direct_vectortemp(2)=r(nxyz(2,np))-r(nxyz(2,1))
		 direct_vectortemp(3)=r(nxyz(3,np))-r(nxyz(3,1))   
		 
		 direct_vectornorm=sqrt(direct_vectortemp(1)**2+
     +	direct_vectortemp(2)**2+direct_vectortemp(3)**2)
		 ori_vector(1)=direct_vectortemp(1)/direct_vectornorm
		 ori_vector(2)=direct_vectortemp(2)/direct_vectornorm
		 ori_vector(3)=direct_vectortemp(3)/direct_vectornorm
		 center=0.0
		 do i=1,np
		 
		 center(1)=r(nxyz(1,i))+center(1)
		 center(2)=center(2)+r(nxyz(2,i))
		 center(3)=center(3)+r(nxyz(3,i))
		 end do
		 
		 center(1)=center(1)/dble(np)
		 center(2)=center(2)/dble(np)
		 center(3)=center(3)/dble(np)
		 
		 call vector2angle(angle_azi,angle_polar,
     +	 ori_vector(1),ori_vector(2),ori_vector(3))
		 
		 end
		 