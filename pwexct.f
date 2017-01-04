	subroutine pwexact(pw, n1, r, nxyz)
	
	  common / num_tot / n
	  common / n_3 / n3
        common / num_par / np
        common / np_3 / np3
	  common / boxlen_z / boxlenz

	  common / radius / a
         common /constraint/ Kconstr,kconstr3,nxyzbin
	  common / brsh / brush_length
                common /up_wall/ upperwall
        common / twowall/ twowallflag
c
c	subroutine to calculate the exact particle-wall mobility tensor
c	declare variables
c
	integer		n1, i, j, n, n3, twowallflag

	       
	integer   Kconstr,kconstr3,nxyzbin
       integer nxyz(3,nxyzbin)

	double precision  a, r(np3), PW(3,3), v, Aw, Bw, boxlenz, 
     +				  brush_length, upperwall
c
	PW = 0.d0

	v = ( r(nxyz(3,n1)) - a ) / a

	if(v.lt.0.d0) v = 1.0d-8

c
c	calculate the A and B constants in the expression for the mobility tensor
c
c	Aw = (368.0*v**3 + 559.0*v**2 + 81.0*v) 
c     +   / (368.0*v**3 + 779.0*v**2 + 250.0*v)
     
      Aw=(28.8158*v**2+13.1198*v+0.232)/
     + (28.9502*v**2+28.3834*v+1.0)
           

	Bw = (6.0*v**2+2.0*v) / (6.0*v**2+9.0*v+2.0)
c
c	start loops over particle coordinates
c
	do i=1,3
	  do j=1,3

	    if(i.eq.j)then
		  if(i.eq.3)then
		    PW(i,j)= 1.d0 / Bw
		  else
		    PW(i,j)=1.d0 / Aw
		  endif
		else   
	      PW(i,j)=0d0
		endif

	  end do
      end do
c
c
c
      if(twowallflag .eq. 1) then
          
      v = (upperwall- r(nxyz(3,n1)) - a ) / a

	if(v.lt.0.d0) v = 1.0d-8

c
c	calculate the A and B constants in the expression for the mobility tensor
c
c	Aw = (368.0*v**3 + 559.0*v**2 + 81.0*v) 
c     +   / (368.0*v**3 + 779.0*v**2 + 250.0*v)
     
      Aw=(28.8158*v**2+13.1198*v+0.232)/
     + (28.9502*v**2+28.3834*v+1.0)
           

	Bw = (6.0*v**2+2.0*v) / (6.0*v**2+9.0*v+2.0)
c
c	start loops over particle coordinates
c
	do i=1,3
	  do j=1,3

	    if(i.eq.j)then
		  if(i.eq.3)then
		    PW(i,j)= pw(i,j)+1.d0 / Bw
		  else
		    PW(i,j)=pw(i,j)+ 1.d0 / Aw
		  endif
	    
		endif

	  end do
      end do
      
      end if
      
	return
	end
c