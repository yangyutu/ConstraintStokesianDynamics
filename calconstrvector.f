      subroutine calconstrvector(granddelta,nxyz,r,constrvector,fixdist)
      	  common / num_par / np
	  common / np_3 / np3
        common /constraint/ Kconstr,kconstr3,nxyzbin
        
      integer i,j,np,np3,Kconstr,kconstr3,nxyzbin
       integer nxyz(3,nxyzbin),mu2,indextemp
      
      double precision granddelta(nxyzbin*3,nxyzbin*3), r(np3), 
     + constrvector(kconstr,np3),fixdist 
      
      
      constrvector=0.0
      do mu=1,Kconstr
      
      do j=1,np
      
      if(mu .eq. 1) then
      
      do i=1,3
      
      constrvector(mu,nxyz(i,j))=(
     + r(nxyz(i,np))-r(nxyz(i,1)))/((np-1)*fixdist)*
     + (Granddelta(j,np)-granddelta(j,1)) 
       
c      write(*,*)  constrvector(mu,nxyz(i,j)),mu,
c     + nxyz(i,j) ,mu,j
      end do
      
      else if( mu .le. ((np-2)+1)) then
      
      indextemp=mu
      
      constrvector(mu,nxyz(1,j))=
     + Granddelta(j,indextemp)-Granddelta(1,j)
     + -(indextemp-1)*Granddelta(np,j)/dble(np-1)+
     + (indextemp-1)*granddelta(1,j)/dble(np-1) 
          
     
      else if (mu .le. (2*(np-2)+1)) then
      
        indextemp=mu-(np-2)
      
            constrvector(mu,nxyz(2,j))=
     + Granddelta(j,indextemp)-Granddelta(1,j)
     + -(indextemp-1)*Granddelta(np,j)/dble(np-1)+
     + (indextemp-1)*granddelta(1,j)/dble(np-1) 
            
      
     
      else 
      
            indextemp=mu-(2*np-4)
      
 
                       constrvector(mu,nxyz(3,j))=
     + Granddelta(j,indextemp)-Granddelta(1,j)
     + -(indextemp-1)*Granddelta(np,j)/dble(np-1)+
     + (indextemp-1)*granddelta(1,j)/dble(np-1) 
     
      
      end if
      
      end do 
      
      end do
      
      
      end