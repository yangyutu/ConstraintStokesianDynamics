	SUBROUTINE ppexct(rexact, n1, n2, t, r, nxyz)

	  common / num_par / np
	  common / num_tot / n
	  common / num_wall / nw
	  common / np_3 / np3

	  common / boxlen_xy / boxlen
	  common / boxlen_z / boxlenz

	  common / radius / a
	  common / delta_fn / delta

	  common / p_p_lub / fx, fy, g1, g2, g3
        common /constraint/ Kconstr,kconstr3,nxyzbin
        common /lub/ cylinderflag
c    *************************************************************************
c    **  subroutine for calculating the exact two particle mobility tensor. **
c    **  (analytical solution as given by Batchelor, 1976)                  **
c    *************************************************************************
 
	  integer	n1, n2, i, j, k, l, alpha, beta, m2, m, nw, n, n3,np3,np

	 integer   Kconstr,kconstr3,nxyzbin
       integer nxyz(3,nxyzbin),cylinderflag

	  double precision a, r(np3), rij(3), delta(6,6), X(2,2), Y(2,2), 
     +			       boxlen, rijsq, rijsep, rexact(6,6), s, m1, 
     +			       fx(0:11), fy(0:11), g1, g2, g3, t, one_4ssq,
     +			       two_pow_m, boxlenz,mabexact(6,6),AA(2,2), BB(2,2)

c    ** calculate separation between the two particles. define dimensionless **
c    **	separation 'u'              
c       here is the documented value for A(2,2) and B(2,2) from Batchelor
           AA(1,1)=0.775
           AA(2,2)=0.775
           AA(1,2)=0.775
           AA(2,1)=0.775
           
           AA(1,1)=AA(1,1)+1e-4
           AA(2,2)=AA(2,2)+1e-4
           AA(1,2)=AA(1,2)-1e-4
           AA(2,1)=AA(2,1)-1e-4
           
           BB(1,1)=0.891
           BB(1,2)=0.490
           BB(2,1)=0.490
           BB(2,2)=0.891
           
	
	rexact = 0.d0
	mabexact=0.0
c
       if(abs(n1-n2) .gt. 1) then
	do i=1, 3

	  rij(i) = r(nxyz(i,n2)) - r(nxyz(i,n1))
	  

c		rij(i) = rij(i) - anint (rij(i)/boxlenz) * boxlenz

c	  endif

	enddo

	rijsq = rij(1)**2 + rij(2)**2 + rij(3)**2
	rijsep = sqrt ( rijsq )
	s = rijsep / a

	if(s.le.2)then

c	  if((n1.lt.nw).or.(n2.lt.nw))write(*,*)t,n1,n2,s

	  s = 2.0 + 1d-8
	  rijsep = s * a
	  rijsq = rijsep ** 2

	endif
c
c
c
			one_4ssq = 1.0 - 4.0 / s**2

			X(1,1) = g1 / one_4ssq - g2 * log(one_4ssq)
     +			   - g3 * one_4ssq * log(one_4ssq) + fx(0) - g1

			do m2 = 1, 5
			  
			  m = 2 * m2

			  two_pow_m = 1.0 / real(2.0 ** m)

			  if(m.eq.2)then
				m1=-2d0
			  else
				m1=real(m-2)
			  endif

			  X(1,1) = X(1,1) + (two_pow_m*two_pow_m*fx(m)
     +				 - g1 - 2d0/real(m) * g2 + 4d0/real(m)/m1 * g3) 
     +				 * (2d0/s) ** m

			end do

			X(2,2)=X(1,1)
c
c
			X(1,2) = -2d0/s * g1 / one_4ssq 
     +			   - g2 * log((s+2.0)/(s-2.0))
     +			   - g3 * one_4ssq * log((s+2)/(s-2)) - 4d0 * g3 / s


			do m2=1,6
			  
			  m = 2 * m2 - 1

			  two_pow_m = 1.0 / (2.0 ** m)

			  if(m.eq.2)then
				m1=-2d0
			  else
				m1=real(m-2)
			  endif

			  X(1,2) = X(1,2) - ( two_pow_m*two_pow_m*fx(m)
     +				 - g1 - 2d0/real(m) * g2 + 4d0/real(m)/m1*g3 ) 
     +				 * (2d0/s)**m

			end do

			X(2,1)=X(1,2)
c
c
			Y(1,1) = -0.167d0 * log(one_4ssq) + fy(0)

			do m2=1,5
			  
			  m = 2 * m2

			  if(m.eq.2)then
				m1=-2d0
			  else
				m1=real(m-2)
			  endif

			  Y(1,1) = Y(1,1) + (1d0/real(2**m)*1d0/real(2**m)*fy(m)
     +				 - 2d0/real(m) * 0.167) * (2d0/s) ** m

			end do

			Y(2,2) = Y(1,1)
c
c
			Y(1,2) = -0.167 * log((s+2)/(s-2))

			do m2=1,6
			  
			  m=2d0*m2-1d0

			  if(m.eq.2)then
				m1=-2d0
			  else
				m1=real(m-2)
			  endif

			  Y(1,2) = Y(1,2) - ( 1d0/real(2**m)*1d0/real(2**m)*fy(m)
     +				 - 2d0/m*0.167 ) * (2d0/s)**m

			end do

			Y(2,1)=Y(1,2)
c
c    **	begin loops over particles and coordinates. i,j are particle indices. **
c    **   alpha, beta are particle coordinates                                  **
c

			  do i=1,2
				if(i.eq.1)then
				  k=0
				else
				  k=3
				endif
				
				do j=1,2
				  if(j.eq.1)then
					l=0
				  else
					l=3
				  endif
				  
				  do alpha=1,3
					do beta=1,3

c    **   calculate the exact particle-particle mobility tensor   **

					  rexact(k+alpha,l+beta)=rexact(k+alpha,l+beta)+
     +				  X(i,j)*rij(alpha)*rij(beta)/rijsq
     +				 +Y(i,j)*(delta(alpha,beta)
     +				 -rij(alpha)*rij(beta)/rijsq)

c    **	end loops over particle coordinates   **

					end do
				  end do

c    **	end loops over particle indices   **

				end do
			  end do  
			  
			  
		else
c       for the neighboring particles
        	do i=1, 3

	  rij(i) = r(nxyz(i,n2)) - r(nxyz(i,n1))

	 enddo
	 
	 s=2.0+1.0e-3
	 
	 one_4ssq = 1.0 - 4.0 / s**2
	 
	 X(1,1)=0.125*1.0/one_4ssq-0.225*log(one_4ssq)-
     + 0.027*one_4ssq*log(one_4ssq)+1.0-0.125+0.9954-1.0
     + +0.25*0.125
     
        X(1,1)=12500+0.307
       X(2,2)=X(1,1)
       X(1,2)=-2.0/s*0.125/one_4ssq-0.225*log((s+2.0)/(s-2.0))-
     + 0.027*one_4ssq*log((s+2.0)/(s-2.0))-4.0*0.027/s+0.3502-
     +0.25*0.125-2.0*0.225*log(2.0)-2.0*0.027
        
       X(1,2)=-12500+0.307
       X(2,1)=X(1,2)
       
       
              if(cylinderflag .eq. 1) then
               X(1,1)=12500+0.307
       X(2,2)=X(1,1)    
              X(1,2)=-12500+0.307
       X(2,1)=X(1,2)
c         the correction will lead to A=0.814       
       else
       X(1,1)=12500+0.3225
       X(2,2)=X(1,1)    
              X(1,2)=-12500+0.3225
       X(2,1)=X(1,2)
c        this correction will lead to A=0.775, just the two sphere model
       end if
       
       Y(1,1)=-0.1667*log((s-2.0))+0.9983
       Y(1,2)=0.1667*log((s-2.0))-0.2737
        Y(2,2)=Y(1,1)
        Y(2,1)=Y(1,2)
	rijsq = rij(1)**2 + rij(2)**2 + rij(3)**2
	rijsep = sqrt ( rijsq )
	
				  do i=1,2
				if(i.eq.1)then
				  k=0
				else
				  k=3
				endif
				
				do j=1,2
				  if(j.eq.1)then
					l=0
				  else
					l=3
				  endif
				  
				  do alpha=1,3
					do beta=1,3

c    **   calculate the exact particle-particle mobility tensor   **

					  rexact(k+alpha,l+beta)=rexact(k+alpha,l+beta)+
     +				 X(i,j)*rij(alpha)*rij(beta)/rijsq
     +				 +Y(i,j)*(delta(alpha,beta)
     +				 -rij(alpha)*rij(beta)/rijsq)

c    **	end loops over particle coordinates   **

					end do
				  end do

c    **	end loops over particle indices   **

				end do
			  end do  

			  
			  end if

	  return
	  end