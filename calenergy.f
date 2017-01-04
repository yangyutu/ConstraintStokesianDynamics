        	SUBROUTINE calenergy(r, nxyz,kappa,tempr,Bpw,energy)

	  common / num_par / np, nw
	  common / np_3 / np3, nw3

	  
	  common / boxlen / boxlenx,boxleny
	  
	  common / radius / a

	  common / grav_pot / Fgrav, Fhw
 
	  common/forcethresh/ ppthresh,wallthresh
	  

c    ******************************************************************
c    ** evlauates the different force components                     **
c    ******************************************************************

	  integer, parameter :: rmax = 1
	  
	  integer    i, j, np, np3, n, n3, nw, nw3
	  
	  integer	 nxyz(3, np)

	  double precision		 a, r(np3), F(np3), rij(3), sep, boxlen, 
     +			 rijsep, Fss, rijtemp(3), brush_length, 
     +			 Fpp, boxlenz, sep2, Fgrav, Fsw, Fpw, Fhw, 
     +			 cutoff, op, deps, pi, rw(nw3),tempr,kappa,Bpp,bpw,
     +  wallthresh, ppthresh,boxlenx,boxleny,Fao,lambda,fcm,dg,energy
     
        double precision,parameter:: kb=1.3806503e-23
	
	  PI = 3.1415927
	  
	  F = 0.d0
	  
	  energy=0.0

c	  Fsw = 0.0
	  
c	  r(1050)=100*a
c
c		DIFFUSER particle-particle forces
cc       Here according to the Newton's third law, the sum of internal force
c       should be zeros
        F=0.0
        energy=0.0


C		DIFFUSERS-WALL PARTICLE FORCES

	  do i=1,np
	  

		  if( r(nxyz(3,i)) .le. a ) then

c			write(*,*) i, j, (rijsep-2*a)
		    Fpp = 0.01
		    energy=energy+1e5

		 else		
c		  Fpp = kb*tempr*Bpw*kappa*1.0e18*exp(-kappa*(r(nxyz(3,i))-a)/a)
c     		Fao=pi*op*kb*tempr*(((rijsep/2.0)**2)-((a*(1+deps))**2))*1.0e18
     		
     		  energy=energy+Bpw*a*exp(-kappa*(r(nxyz(3,i))-a)/a)
     		
     		    Fao=0.0

            Fpp=Fpp+Fao

		  endif

		  F(nxyz(3,i))=F(nxyz(3,i))+Fpp

	   

 
c		gravity!!!!

		F(nxyz(3,i))=F(nxyz(3,i))-Fgrav
		energy=energy+Fgrav*r(nxyz(3,i))*1.0e-18/kb/tempr

	  end do


	

	  return
	  end
