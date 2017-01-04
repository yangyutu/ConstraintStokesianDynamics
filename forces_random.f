	SUBROUTINE forces_random(F, r, nxyz,granddelta,  rgrnd,dt,fixdist,
     +  constrvector,randomforce_unproj)

	  common / num_par / np
	  common / np_3 / np3

	  
	  common / boxlen / boxlen
	  
	  common / radius / a


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
     + tempsum,randomforce_proj(np3),F_ps(np3),
     + H_hatmatrix(kconstr,kconstr),Q_vector(kconstr),
     + lam_constr(kconstr),fixdist,testsum(100),F_conser(np3)
     

        double precision, parameter :: kb = 1.380658E-23, pi = 3.1416

	  F = 0.d0
	  
        rgrndtemp=rgrnd
	  Fsw = 0.0
c        open(131,file='mobility.txt')
c        open(132,file='constrain.txt')
c        do i=1,np3
c        do j=1,np3	
c        write(131,*) rgrnd(i,j)
c        end do
c        end do
c        close(131)
c       first calculate the constraint force directional vector
        call calconstrvector(granddelta,nxyz,r,constrvector,fixdist)

c       second generate the unprojected random force vector
         call choldc(rgrndtemp, CholDiag)
	    call random(CholDiag, rgrndtemp, randisp)
	    do i=1,np3
        randomforce_unproj(i)=randisp(i)
         end do
         
c       now calculate the projected random force        
c        randomforce_unproj=1
        G_matrix=0.0
        do mu=1,kconstr
        do nv=1,kconstr
        
        do j=1,np3	
        	
        G_matrix(mu,nv)=G_matrix(mu,nv)+constrvector(mu,j)*
     + constrvector(nv,j)
        end do
        end do
        end do
        
        call invert(G_matrix,kconstr)
        P_vector=0.0
        
        do mu=1,kconstr
        do j=1,np3
        P_vector(mu)=P_vector(mu)+randomforce_unproj(j)*
     +  constrvector(mu,j) 
        end do
        end do
        
        
        
        randomforce_hardcom=0.0
        do mu=1,kconstr
        
        do j=1,kconstr
        randomforce_hardcom(mu)=randomforce_hardcom(mu)+
     +  G_matrix(mu,j)*P_vector(j)
        
        end do
        
        end do
        
        
        do j=1,np3
        tempsum=0.0
        
        do mu=1,kconstr
        
        tempsum=tempsum+randomforce_hardcom(mu)*constrvector(mu,j)
        
        end do
        
        randomforce_proj(j)=randomforce_unproj(j)-tempsum
        
        end do
        
        
c       third calculate the corrective pseudo force (we later calculate it)

       
        
c        F_ps=0.0
          testsum=0.0
          do mu=1,kconstr
          do i=1,np3
         testsum(mu)=testsum(mu)+randomforce_proj(i)*constrvector(mu,i)
         end do
         end do
            

        F=randomforce_proj
        
        
                
        
		 

	  return
	  end
