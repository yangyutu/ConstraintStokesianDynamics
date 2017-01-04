      subroutine vector2angle(angle_azi,angle_polar,
     + ori_vectorx,ori_vectory, ori_vectorz)
      
      integer j
      
      
      double precision angle, ori_vectorx,ori_vectory,
     +  angletemp,angle_azi,angle_polar,ori_vectorz,rtemp
     
     
     
       
       angle_azi=atan2(real(ori_vectory),real(ori_vectorx))
       
c       if(ori_vectory .lt. 0) then 
c       angle=angletemp+2*3.1415926
c       else
        rtemp=sqrt(ori_vectory**2+ori_vectorx**2)
       angle_polar=atan2(real(rtemp),real(ori_vectorz))
c       end if
       
       

       end 
       
        
      
      
      