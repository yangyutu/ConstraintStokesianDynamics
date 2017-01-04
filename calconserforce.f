	SUBROUTINE calconserF(F, r, nxyz,kappa,tempr,Bpw, ori_vector1)

	  common / num_par / np, nw
	  common / np_3 / np3, nw3

	  
	  common / boxlen / boxlenx,boxleny
	  
	  common / radius / a

	  common / grav_pot / Fgrav, Fgrav2, Fhw
c	  common / dl / cutoff, op, deps,kappa,Bpp 
	  common/forcethresh/ ppthresh,wallthresh
	      common /up_wall/ upperwall
        common / twowall/ twowallflag
          common/ asym/ asymflag
           common/active/ activeforce
c    ******************************************************************
c    ** evlauates the different force components                     **
c    ******************************************************************

	  integer, parameter :: rmax = 1
	  
	  integer    i, j, np, np3, n, n3, nw, nw3
	  
	  integer	 nxyz(3, np),asymflag,twowallflag

	  double precision		 a, r(np3), F(np3), rij(3), sep, boxlen, 
     +			 rijsep, Fss, rijtemp(3), brush_length, 
     +			 Fpp, boxlenz, sep2, Fgrav, Fsw, Fpw, Fhw, 
     +			 cutoff, op, deps, pi, rw(nw3),tempr,kappa,Bpp,Bpw,
     +  wallthresh, ppthresh,boxlenx,boxleny,Fao,lambda,fcm,dg, Fgrav2,
     + upperwall,activeforce,ori_vector1(3)
     
        double precision,parameter:: kb=1.3806503e-23	
	  PI = 3.1415927	  
	  F = 0.d0
C		DIFFUSERS-WALL PARTICLE FORCES
	  do i=1,np
		  if( r(nxyz(3,i)) .le. a ) then
c			write(*,*) i, j, (rijsep-2*a)
		    Fpp = 0.01

		 else		
		  Fpp = kb*tempr*Bpw*kappa*1.0e18*exp(-kappa*(r(nxyz(3,i))-a)/a)
c     		Fao=pi*op*kb*tempr*(((rijsep/2.0)**2)-((a*(1+deps))**2))*1.0e18
            if(twowallflag .eq. 1 ) then
           Fpp = Fpp-kb*tempr*Bpw*kappa*1.0e18* 
     +  exp(-kappa*(upperwall-r(nxyz(3,i))-a)/a) 
            end if     		
c     		    Fao=0.0
c            Fpp=Fpp+Fao
		  endif
		  F(nxyz(3,i))=F(nxyz(3,i))+Fpp
c		gravity!!!!

		F(nxyz(3,i))=F(nxyz(3,i))-Fgrav
          
          
c         here i average the total active force onto each bead 
              
          F(nxyz(1,i))=F(nxyz(1,i))+activeforce*
     +     ori_vector1(1)/dble(np)/a*kb*tempr*1e18
      F(nxyz(2,i))=F(nxyz(2,i))+activeforce*ori_vector1(2)
     +     /dble(np)/a*kb*tempr*1e18
         
       F(nxyz(3,i))=F(nxyz(3,i))+activeforce*ori_vector1(3)
     +     /dble(np)/a*kb*tempr*1e18
c         change the force unit to nN    
          
          
          
          if(asymflag .eq. 1) then
          
          if(mod(np,2) .eq. 1) then
              if( i .lt. ((dble(np)+1)/2.0)) then
              F(nxyz(3,i))=F(nxyz(3,i))-(Fgrav2-Fgrav)
              elseif( i .eq. ((np+1)/2)) then
              F(nxyz(3,i))=F(nxyz(3,i))-0.5*(Fgrav2-Fgrav)
              end if
      else
            if( i .lt. ((dble(np)+1.0)/2.0)) then
              F(nxyz(3,i))=F(nxyz(3,i))-(Fgrav2-Fgrav)
      endif
      end if
          end if
          
          

	  end do


	

	  return
	  end
