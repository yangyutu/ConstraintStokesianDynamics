      subroutine calvelocity(r, nxyz,r_rod, u, U_rod,
     + ori_vector1,ori_vector2,ori_vector3  ) 		
      		

	  common / num_par / np
	  common / np_3 / np3
	  common / radius / a 		
       common /constraint/ Kconstr,kconstr3,nxyzbin
       common/ permu/ eps
      	  integer, parameter :: m = 2, maxbin=5000  
	  integer    i, j, np, np3, n, n3,mu, nv,kconstr,kconstr3,nxyzbin
	  
	  integer	 nxyz(3, nxyzbin) ,centralindex,rescalecheck,step
	  double precision r(np3),a,
     +   direct_vectortemp(3),direct_vectornorm, direct_vector(3), 
     +  center(3),vcenter(3),vrelative(np3),v_angu(np),u(np3),
     +  D(np3,np3),F(np3),fac3,totalF(3),v_para,v_vert,
     +  F_para, F_vert,tempr,mobility_para,mobility_vert,mobility_anglu,
     +diff_para,diff_vert,Diff_anglu,U_rod(6),omega(3),mass(3),r_rod(3),
     + omega_para,omega_vert1,omega_vert2 ,
     + ori_vector1(3),ori_vector2(3),ori_vector3(3) ,distance  
        
        double precision eps(3,3,3) 
     
        double precision,parameter:: kb=1.3806488e-23
c       first calculat the velocity

     
        
       mass(1)=r_rod(1)
       mass(2)=r_rod(2)
       mass(3)=r_rod(3)
       
          
c       now calculate of the tranlation velocity and totalforce

        vcenter=0
        do i=1,np
        vcenter(1)=vcenter(1)+u(nxyz(1,i)) 
        vcenter(2)=vcenter(2)+u(nxyz(2,i))
        vcenter(3)=vcenter(3)+u(nxyz(3,i))     
        end do
        
        vcenter(1)=vcenter(1)/dble(np)
        vcenter(2)=vcenter(2)/dble(np)
        vcenter(3)=vcenter(3)/dble(np)
        
       U_rod(1)=vcenter(1)
       U_rod(2)=vcenter(2)
       U_rod(3)=vcenter(3)
       
        
c       now calculate the relative velocity of the each bead

        do i=1,np
        
        vrelative(nxyz(1,i))=u(nxyz(1,i))-vcenter(1)
        vrelative(nxyz(2,i))=u(nxyz(2,i))-vcenter(2)
        vrelative(nxyz(3,i))=u(nxyz(3,i))-vcenter(3)           
        
        end do
        
                distance=(r(nxyz(1,1))-mass(1))**2+
     +   (r(nxyz(2,1))-mass(2))**2+(r(nxyz(3,1))-mass(3))**2
                
                 omega=0.0
           do i=1,3
                do j=1,3
                    do k=1,3
        omega(i)=omega(i)+eps(i,j,k)*(r(nxyz(j,1))-mass(j))*
     +   vrelative(nxyz(k,1))/distance
                    end do
                end do
           end do
           
           
            omega_para=omega(1)*ori_vector1(1)+omega(2)*ori_vector1(2)
     + + omega(3)*ori_vector1(3)
         
        omega_vert1=omega(1)*ori_vector2(1)+omega(2)*ori_vector2(2)
     + + omega(3)*ori_vector2(3)
        
         omega_vert2=omega(1)*ori_vector3(1)+omega(2)*ori_vector3(2)
     + + omega(3)*ori_vector3(3)
         
         U_rod(4)=omega_vert1
         U_rod(5)=omega_vert2
        
        
         
        end          
                       