      subroutine caldiff(ori_vector1,ori_vector2,ori_vector3, F_net,
     + U_rodtemp,r_rod,r,nxyz,tempr)
      
      common / num_par / np
      common / np_3 / np3
      common / radius / a
            common /constraint/ Kconstr,kconstr3,nxyzbin
      integer    i, j, np,np3,Kconstr,kconstr3,nxyzbin,nxyz(3,nxyzbin)
      double precision a, r(np3), F_net(np3), rij(3), sep, boxlen, 
     +   tempr,U_rod(6), massx,massy,massz, r_rod(3),F_total(3),
     + F_general(6),totaltorque(3), ori_vector1(3), ori_vector2(3), 
     + ori_vector3(3),M(6),Diff(6),U_rodtemp(6) 
      
      double precision, parameter :: kb = 1.380658E-23, pi = 3.1416 
      massx=r_rod(1)
      massy=r_rod(2)
      massz=r_rod(3)
      U_rod=U_rodtemp
      
      F_total=0
c         unit of F is nN      
      do i=1,np
      F_total(1)=F_total(1)+F_net(nxyz(1,i))
      F_total(2)=F_total(2)+F_net(nxyz(2,i))
      F_total(3)=F_total(3)+F_net(nxyz(3,i))
      end do

c       now i need to calculate the parallel force, perpendicular force
        F_general(1)=F_total(1)*ori_vector1(1)+F_total(2)*ori_vector1(2)
     + +F_total(3)*ori_vector1(3)   

        F_general(2)=F_total(1)*ori_vector2(1)+F_total(2)*ori_vector2(2)
     + +F_total(3)*ori_vector2(3)
        
        F_general(3)=F_total(1)*ori_vector3(1)+F_total(2)*ori_vector3(2)
     + +F_total(3)*ori_vector3(3) 
        
        U_rod(1)=U_rodtemp(1)*ori_vector1(1)+U_rodtemp(2)*ori_vector1(2)
     + +U_rodtemp(3)*ori_vector1(3)
      U_rod(2)=U_rodtemp(1)*ori_vector2(1)+U_rodtemp(2)*ori_vector2(2)
     + +U_rodtemp(3)*ori_vector2(3)
      U_rod(3)=U_rodtemp(1)*ori_vector3(1)+U_rodtemp(2)*ori_vector3(2)
     + +U_rodtemp(3)*ori_vector3(3)
        
      totaltorque=0.0
        do j=1,np

        totaltorque(1)=totaltorque(1)+(r(nxyz(2,j))-massy)*
     =   F_net(nxyz(3,j))-(r(nxyz(3,j))-massz)*F_net(nxyz(2,j))
        totaltorque(2)=totaltorque(2)+(r(nxyz(3,j))-massz)*
     +   F_net(nxyz(1,j))-(r(nxyz(1,j))-massx)*F_net(nxyz(3,j))
        totaltorque(3)=totaltorque(3)+(r(nxyz(1,j))-massx)*
     +   F_net(nxyz(2,j))-(r(nxyz(2,j))-massy)*F_net(nxyz(1,j))
    
        end do
        
             	  F_general(4)=totaltorque(1)*ori_vector2(1)+
     +   totaltorque(2)*ori_vector2(2)
     + +totaltorque(3)*ori_vector2(3)  
        F_general(5)=totaltorque(1)*ori_vector3(1)+
     +   totaltorque(2)*ori_vector3(2)
     + +totaltorque(3)*ori_vector3(3)
        
c         velocity is nm/s    , rad/s    
        do jj=1,5     	         
     	  M(jj)=U_rod(jj)/F_general(jj)     	  
     	  Diff(jj)=kb*Tempr*M(jj)   
        end do
              write(146,'(i8,E12.5, 3f14.5,5E14.5)') np, r_rod(3)/a,
     +	  ori_vector1(1),ori_vector1(2),
     +	   ori_vector1(3), 
     +	  Diff(1),Diff(2),Diff(3), Diff(4)*1e18,
     +  Diff(5)*1e18
  
        
        
      end
      