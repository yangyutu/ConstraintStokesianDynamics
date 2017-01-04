	SUBROUTINE netforces(F, r, nxyz,granddelta, H, rgrnd,dt,fixdist,
     +  constrvector,F_conser,F_random,F_ps)

	  common / num_par / np
	  common / np_3 / np3

	  
	  common / boxlen / boxlen
	  
	  common / radius / a

	  common / grav_pot / Fhw

	  common / dl / kappa, cutoff, bpp
	  
        	  common / f_pot / lambda, fcm, dg,tempr
  
        common /constraint/ Kconstr,kconstr3,nxyzbin
c    ******************************************************************
c    ** evlauates the different force components                     **
c    ******************************************************************

	  integer, parameter :: rmax = 1
	  
	  integer    i, j, np, np3, n, n3,mu, nv,kconstr,kconstr3,nxyzbin
	  
	  integer	 nxyz(3, nxyzbin)

	  double precision		 a, r(np3), F(np3), rij(3), sep, boxlen, 
     +			 rijsep, Fss, rijtemp(3), brush_length, 
     +			 Fpp, boxlenz, sep2, Fgrav, Fsw, Fpw, Fhw,tempr, 
     +			 kappa, cutoff, bpp,lambda, fcm, dg
	 double precision		 Exi, Eyi, Ezi, Exj, Eyj, Ezj, Fo, 
     +						 F1, F2, F3, felx, fely, felz

	  double precision		 RT, EX1, EY1, EX2, EY2, EX3, EY3,
     +						 EX4, EY4, RX, RY, STEP, EMAGI, EMAGJ,EMAG
     
        double precision dE2x, dE2y, Fdepx, Fdepy, rdep,
     +   granddelta(nxyzbin*3,nxyzbin*3),  
     + constrvector(kconstr,np3), H(np3,np3),rgrnd(np3,np3),
     + randomforce_unproj(np3), choldiag(np3), randisp(np3),
     + rgrndtemp(np3,np3),dt,G_matrix(kconstr,kconstr),
     + P_vector(kconstr),randomforce_hardcom(kconstr),
     + tempsum,randomforce_proj(np3),
     + H_hatmatrix(kconstr,kconstr),Q_vector(kconstr),
     + lam_constr(kconstr),fixdist,testsum(100),F_conser(np3),
     +   F_random(np3),F_unconstr(np3),F_ps(np3)
     

        double precision, parameter :: kb = 1.380658E-23, pi = 3.1416

	  F = 0.d0
	  
        F_unconstr=F_conser+F_random+F_ps
        
c       third calculate the corrective pseudo force (we later calculate it)

c       call calfps(F_ps, r, nxyz, granddelta,fixdist)
        
c        F_ps=0.0
               
                
c       fourth, calculate the contraint force lam_constr
c       here is a test for the solver
c        randomforce_proj=1
c        open(131,file='mobility.txt')
c        open(132,file='constrain.txt')
c        do i=1,np3
c        do j=1,np3	
c        write(131,*) rgrnd(i,j)
c        end do
c        end do
c        do mu=1,kconstr
c        do i=1,np3
c        write(132,*) constrvector(mu,i)
c        end do
c        end do
        call calconstrvector(granddelta,nxyz,r,constrvector,fixdist)
                H_hatmatrix=0.0
        do mu=1,kconstr
        
        do nv=1,kconstr
        
        do i=1,np3
        do j=1,np3	
  	
        H_hatmatrix(mu,nv)=H_hatmatrix(mu,nv)+constrvector(mu,i)*
     + H(i,j)*constrvector(nv,j)
c       the last two terms are the total constraint force in the i direction     
        end do
        end do
        end do
        end do
        
        call invert(H_hatmatrix,kconstr)
        
        Q_vector=0.0
        
        do mu=1,kconstr
        
        do i=1,np3
        do j=1,np3
        
        Q_vector(mu)=Q_vector(mu)+constrvector(mu,i)*H(i,j)*
     +     F_unconstr(j)
        end do
        end do
        end do
        
        lam_constr=0.0
        
        
        do nv=1,kconstr
        do mu=1,kconstr
        lam_constr(nv)=lam_constr(nv)+
     +  H_hatmatrix(nv,mu)*Q_vector(mu)
        
        end do
        
        end do
        
c       now calculate the total force        
        do j=1,np3
        
        tempsum=0
        do mu=1,kconstr
        
        tempsum=tempsum-lam_constr(mu)*constrvector(mu,j)
        
        end do
        
        F(j)=F_unconstr(j)+tempsum
        
        end do
   
        
        
             
                
        
		 

	  return
	  end
