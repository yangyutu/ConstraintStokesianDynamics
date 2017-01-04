      subroutine permutator()
      
      common/ permu/ eps
      
      double precision eps(3,3,3)
      integer i,j ,k
      
      eps=0
      do i=1,3
          do j=1,3
              do k=1,3
                  if(i .eq. 1 .and. j .eq. 2 .and. k .eq. 3) then
                      eps(i,j,k)=1
                  end if
                  if(i .eq. 2 .and. j .eq. 3 .and. k .eq. 1) then
                      eps(i,j,k)=1
                  end if
                  if(i .eq. 3 .and. j .eq. 1 .and. k .eq. 2) then
                      eps(i,j,k)=1
                  end if
                    if(i .eq. 1 .and. j .eq. 3 .and. k .eq. 2) then
                      eps(i,j,k)=-1
                  end if
                  if(i .eq. 3 .and. j .eq. 2 .and. k .eq. 1) then
                      eps(i,j,k)=-1
                  end if
                  if(i .eq. 2 .and. j .eq. 1 .and. k .eq. 3) then
                      eps(i,j,k)=-1
                  end if
                  
              end do
          end do
      end do
      end
              