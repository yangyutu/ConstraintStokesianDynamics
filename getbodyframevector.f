      subroutine getbodyframevector(ori_vector1,ori_vector2,ori_vector3)
      
        common/ permu/ eps
        common /np_rod/ nrod,nrod3,npinrod
        
          integer, parameter :: m = 2, maxbin=5000,nmax=3000,nmax2=300 
        integer nrod,nrod3, npinrod
       double precision ori_vector1(3),ori_vector2(3),
     +  ori_vector3(3)
      
      
      double precision eps(3,3,3),ez(3)
      
      double precision i,j,k,ii,norm
      
      ez(1)=0
      ez(2)=0
      ez(3)=1
      ori_vector2=0
      ori_vector3=0
         
      
      
    
 
      
      
      norm=0
      
      
      
          
                norm=0.0
      do i=1,3
          norm=norm+ori_vector1(i)**2
      end do
      norm=sqrt(norm)
      do i=1,3
      ori_vector1(i)=ori_vector1(i)/norm
      end do
          
          
      
      if(ori_vector1(1) .eq. 0 .and. ori_vector1(2) .eq.0) then
        ori_vector2(1)=0
        ori_vector2(2)=1
        ori_vector2(3)=0
      else
      do i=1,3
          do j=1,3
              do k=1,3
      ori_vector2(i)=ori_vector2(i)+eps(i,j,k)*
     +  ez(j)*ori_vector1(k)
      
              end do
          end do
      end do
      norm=0.0
      do i=1,3
          norm=norm+ori_vector2(i)**2
      end do
      norm=sqrt(norm)
      do i=1,3
      ori_vector2(i)=ori_vector2(i)/norm
      end do
      
      end if
      do i=1,3
          do j=1,3
              do k=1,3
      ori_vector3(i)=ori_vector3(i)+
     +  eps(i,j,k)*ori_vector2(j)*ori_vector1(k)
              end do
          end do
      end do
      
     
      
    
      
      
          
      end