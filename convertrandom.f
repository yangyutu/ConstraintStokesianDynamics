      subroutine convertrandom2(ori_vector1, ori_vector2, ori_vector3,
     + ori_vector1mid,ori_vector2mid,ori_vector3mid, F_random,nxyz  )
      
      common / num_par / np
      common / np_3 / np3
      common /constraint/ Kconstr,kconstr3,nxyzbin
      
      integer    i, j,k, np, np3, n, n3,mu, nv,kconstr,kconstr3,nxyzbin
	  
	  integer	 nxyz(3, nxyzbin)
        
      double precision ori_vector1(3), ori_vector2(3), ori_vector3(3),
     + ori_vector1mid(3),ori_vector2mid(3),ori_vector3mid(3),tran(3,3),
     + F_random2(np3), F_random_body1(np3),F_random_body2(np3),
     + F_random(np3)
      
      
c     first transform the F_random in lab frame into the initial body frame
c     then transform to intermediate body frame
c     finally transform back to lab frame
      
          do k=1,3
              tran(1,k)=ori_vector1(k)
              tran(2,k)=ori_vector2(k)
              tran(3,k)=ori_vector3(k)
          end do
          
      F_random_body1=0.0
      do i=1,np
          do j=1,3
              do k=1,3
      F_random_body1(nxyz(j,i))=F_random_body1(nxyz(j,i))+
     + tran(j,k)*F_random(nxyz(k,i))
              end do
          end do
      end do
      
      tran=0
          do k=1,3
              tran(1,1)=tran(1,1)+ori_vector1mid(k)*ori_vector1(k)
              tran(1,2)=tran(1,2)+ori_vector1mid(k)*ori_vector2(k)
              tran(1,3)=tran(1,3)+ori_vector1mid(k)*ori_vector3(k)
              tran(2,1)=tran(2,1)+ori_vector2mid(k)*ori_vector1(k)
              tran(2,2)=tran(2,2)+ori_vector2mid(k)*ori_vector2(k)
              tran(2,3)=tran(2,3)+ori_vector2mid(k)*ori_vector3(k)
              tran(3,1)=tran(3,1)+ori_vector3mid(k)*ori_vector1(k)
              tran(3,2)=tran(3,2)+ori_vector3mid(k)*ori_vector2(k)
              tran(3,3)=tran(3,3)+ori_vector3mid(k)*ori_vector3(k)      
          end do
      
      F_random_body2=0.0
      do i=1,np
          do j=1,3
              do k=1,3
      F_random_body2(nxyz(j,i))=F_random_body2(nxyz(j,i))+
     + tran(j,k)*F_random_body1(nxyz(k,i))
              end do
          end do
      end do
      
      
       do k=1,3
              tran(k,1)=ori_vector1mid(k)
              tran(k,2)=ori_vector2mid(k)
              tran(k,3)=ori_vector3mid(k)
       end do
       
             F_random2=0.0
      do i=1,np
          do j=1,3
              do k=1,3
      F_random2(nxyz(j,i))=F_random2(nxyz(j,i))+
     + tran(j,k)*F_random_body2(nxyz(k,i))
              end do
          end do
      end do
      
      end
      