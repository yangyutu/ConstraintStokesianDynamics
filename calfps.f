      subroutine calfps(F_ps, r, nxyz, granddelta,fixdist)
      
      	  common / num_par / np
	  common / np_3 / np3
        common /constraint/ Kconstr,kconstr3,nxyzbin
               	  common / f_pot / lambda, fcm, dg,tempr
      integer i,j,np,jj,np3,Kconstr,kconstr3,nxyzbin,mu,nv
       integer	 nxyz(3, nxyzbin)
      
       double precision r(np3), G_matrix(kconstr,kconstr),
     + constrvector(kconstr,np3),r_plus(np3),h,r_minus(np3),
     + G_matrix_plus(kconstr,kconstr),detG_plus,detG_minus,
     + p(kconstr),granddelta(nxyzbin*3,nxyzbin*3),F_ps(np3),
     + lambda,fcm,dg,tempr,fixdist,G_matrix_minus(kconstr,kconstr)
      double precision, parameter :: kb = 1.380658E-23
      
      h=0.01
      F_ps=0
      do jj=1,np3
c           only the first bead and the last bead position has relation with G       
      if( jj .le. 3 .or. jj .ge. (np3-3)) then
      r_plus=r
      
      r_plus(jj)=r_plus(jj)+h
      
      call calconstrvector(granddelta,nxyz,r_plus,constrvector,fixdist)
      
      
       G_matrix_plus=0.0
        do mu=1,kconstr
        do nv=1,kconstr
        
        do j=1,np3	
        	
        G_matrix_plus(mu,nv)=G_matrix_plus(mu,nv)+constrvector(mu,j)*
     + constrvector(nv,j)
        end do
        end do
        end do
        
        
        call choldcorigin(G_matrix_plus, kconstr, p)
        
        detG_plus=p(1)
        do i=2,kconstr
        detG_plus=p(i)*detG_plus
        end do
        
        detG_plus=sqrt(detG_plus)
        
             r_minus=r
      
      r_minus(jj)=r_minus(jj)-h
      
      call calconstrvector(granddelta,nxyz,r_minus,constrvector,fixdist)
      
      
       G_matrix_minus=0.0
        do mu=1,kconstr
        do nv=1,kconstr
        
        do j=1,np3	
        	
        G_matrix_minus(mu,nv)=G_matrix_minus(mu,nv)+constrvector(mu,j)*
     + constrvector(nv,j)
        end do
        end do
        end do 
        
       call choldcorigin(G_matrix_minus, kconstr,  p)
        
        detG_minus=p(1)
        do i=2,kconstr
        detG_minus=p(i)*detG_minus
        end do
        
        detG_minus=sqrt(detG_minus)
        
        F_ps(jj)=1.0e18*kb*tempr*(1/2.0/h)*
     +  (log(detG_plus)-log(detG_minus))
     
        end if
c           here kT/nm is 10e18 nN        
        end do
        
        end 