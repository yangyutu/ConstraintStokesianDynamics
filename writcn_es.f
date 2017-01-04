	
      SUBROUTINE writ_head (nstep, iprint, a, tempr, phi, dt, np)

	integer np, iprint, nstep
	double precision tempr, phi, dt, a

	write(20,*)nstep*dt
	write(20,*)iprint*dt
	write(20,*)np
	write(20,*)2*a
	write(20,*)5
	write(20,*)'SD simulations for polymer-coated particles near'
      write(20,*)'a hard interface with complete hydrodynamics'
	write(20,*)'area fraction (based on core radius) =',phi
	write(20,*)'Temperature (C) =',tempr
	write(20,*)' '

	return
	end
	  
	  
	  
	  
	  SUBROUTINE writcn ( t, r, nxyz )

        common / num_par / np
	  common / np_3 / np3
         common / radius / a

c    *******************************************************************
c    ** subroutine to write out the configuration to unit 20          **
c    *******************************************************************

	  integer	n, np3, np, i 
	  
	  integer	nxyz(3, np)

        double precision		t, r(np3),a


c   ********************************************************************
	

	  do i = 1, np

		write(20,10) i,r(nxyz(1,i))/a,r(nxyz(2,i))/a, r(nxyz(3,i))/a

	  end do			

10	format (i7,3f14.5)


        return
        end
