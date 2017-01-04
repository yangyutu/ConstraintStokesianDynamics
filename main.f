	PROGRAM SBD

	  common / num_par / np
	  common / np_3 / np3

	  common / boxlen / boxlen
	  common / radius / a



	  common / p_p_lub / fx, fy, g1, g2, g3

	  common / seed / iDummy

	 common / grav_pot /  Fgrav,Fgrav2, Fhw


	  common / delta_fn / delta

	  common / dl / kappa, cutoff, bpp
	  
        common / f_pot / lambda, fcm, dg,tempr

        COMMON / ORDER_PAR / CONN6AVG, RMIN,PSI6,RGMEAN
        
        COMMON / EXP_PARA / expbox, VAR
        
        common /constraint/ Kconstr,kconstr3,nxyzbin
        
        common /lub/ cylinderflag
        
             common /up_wall/ upperwall
        common / twowall/ twowallflag
        common/ asym/ asymflag
        common/active/ activeforce
c    ******************************************************************
c    ** This program runs the actual simulation and writes           **
c    ** coordinates into an output file.                             **
c    ******************************************************************
	  
	  integer, parameter ::  maxbin=5000
	  
	  integer          i, j, k, l, np, nw, n3, np3, n, nw3,mu,nv
	  
	  integer		   nstep, iprint, istart, iDummy,nxyzbin

	  integer, allocatable :: nxyz(:, :)

	  double precision, parameter :: sqpi = 1.7725,m = 2, 
     +								 twopi = 6.2831853, tol=1d-5

	  double precision, allocatable :: F(:), D(:,:), rgrnd(:,:),
     +					   choldiag(:), randisp(:), 
     +					   rmove(:), r0(:), u0(:), u(:),u_final(:),
     +					   ud(:), rinfi(:,:), mwall(:,:),
     +					   mw2b(:,:), r(:), conhist(:), rghist(:),
     +              psihist(:),granddelta(:,:),F_ps(:),
     + constrvector(:,:),velocity_origin(:),velocity_pseudo(:),
     + F_conser(:),F_net(:),F_random(:),distcheck(:),
     + randomforce_unproj(:),F_net2(:),F_random2(:),u2(:)   

	  double precision t, delta(6,6), a, tempr, dt,
     +				   boxlen, zi, ms(3,3), brush_length, dt_m,
     +				   m_2, fac1, fac2, istart_dt, t_min_isdt, vol, 
     +				   tpi_bxln, fx(0:11), fy(0:11), g1, g2, g3, 
     +				   boxlenz, tpi_bxlnz, phi, deltat, Fgrav, Fhw

        character        par_in*60, par_out*60, check, 
     +	cphi*60, cconc*60, dss_out*60,RGHISTFILE*30, psihistfile*30,
     + conhistfile*30,filetemp*3


	  double precision p, rcut, kcut, rr(3, -20:20), kk(3, -20:20)

	  integer rmax, kmax, bin, rgbin,conbin, psibin,conbinmax,
     +  rgbinmax, psibinmax,Kconstr,kconstr3,rescalecheck,cyclenum,
     + cycleindex,hydroflag,cylinderflag,twowallflag,asymflag,
     +   caldiffflag
		
	  double precision dx, dy, dz, dxt, dyt, dzt, kappa, cutoff, bpp,
     +   drbin, rlower, rupper, ghist(maxbin), nideal, rxij,ryij,rzij,
     +    rij,lambda, fcm, dg,delrg,delpsi,delcon,rgmin,rgmean,rmin,
     + rescaleflag,direct_vectortemp(3),direct_vector(3),center(3),
     + direct_vectornorm,ori_vector1(3),ori_vector2(3),ori_vector3(3),
     + ori_vector1mid(3),ori_vector2mid(3), ori_vector3mid(3),angle,
     + ori_vectorxtemp, ori_vectorytemp, ori_vectornorm,ori_vectorztemp,
     + twodimflag, angle_azi, angle_polar,testsum(maxbin),
     + testsum2(20,20),toler,preheight,Fac3,rhof,rhop,en_total,bpw,grav,
     + r_rod(3),r_rodmid(3), U0_rod(6),Ud_rod(6),rtemp(3),Udrift(6),
     + U_rodmid(6),centerindex,Udrift2(6),U_rodfinal(6),
     + upperwall,grav2, Fgrav2,activeforce,F_net3(6),U_rodfinal2(6)
     
        double precision ophmax,expbox(2),conn6avg,psi6,fixdist,vis,norm
        
        double precision, parameter ::kb=1.3806488e-23, pi=3.14159265359
     

	open(1,file='run.txt')
	read(1,*)
	read(1,*)np
	read(1,*)
	read(1,*)nstep
	read(1,*)
	read(1,*)iprint
	read(1,*)
	read(1,*)istart
	read(1,*)
	read(1,*)a
	read(1,*)
	read(1,*)tempr
	read(1,*)
	read(1,*)dt
	read(1,*)
	read(1,*)t
	read(1,*) 
	read(1,*) idummy
      read(1,*) 
      read(1,*) fixdist
      read(1,*)
      read(1,*) rescaleflag
      read(1,*)
      read(1,*) cyclenum
      read(1,*) 
      read(1,*) ori_vectorxtemp, ori_vectorytemp, ori_vectorztemp
      read(1,*)
      read(1,*) preheight
      read(1,*) 
      read(1,*) toler
      read(1,*)
	read(1,*)rhof
	read(1,*)
	read(1,*)rhop
	read(1,*)
      read(1,*) Grav
      read(1,*) 
	read(1,*)kappa
	read(1,*)
	read(1,*) Bpw
      read(1,*)
      read(1,*) hydroflag
      read(1,*)
      read(1,*) cylinderflag
      read(1,*)
      read(1,*) twowallflag
      read(1,*)
      read(1,*) upperwall
      read(1,*)
      read(1,*) asymflag
      read(1,*)
      read(1,*) grav2
      read(1,*) 
      read(1,*) activeforce
      read(1,*)
      read(1,*) caldiffflag
      
      if(caldiffflag .eq. 1) then
          open(146,file='cal_diff.txt')
      end if
      

	  
	  fixdist=2.0*a+fixdist*a
	  preheight=preheight*a
	  toler=toler*a
	  dt=dt/1000.0
	  Kconstr=3*np-5
	  kconstr3=3*kconstr
        nxyzbin=max(kconstr,np)
       call permutator()

        
	  cutoff = cutoff * a
	  expbox = expbox * a
	  dg=dg*a
	  	  	rhof=rhof*1000.0
	      rhop=rhop*1000.0
	

	Fgrav=1e9*4.0/3.0*3.1415926*(rhop-rhof)*((a*1.0e-9)**3.0)*9.8066
      Fgrav=Grav*kb*tempr/a*1e18
	  Fgrav2=Grav2*kb*tempr/a*1e18

  

	  par_out = 'sd_xyz.txt'

	  dss_out = 'sd_dss.txt'
	  
	  open(144,file='distcheck.txt')
	  open(145,file='rodop.txt')



        
  	
	  np3 = np * 3

	  allocate ( F(np3), D(np3, np3), rgrnd(np3, np3), 
     +			 choldiag(np3), randisp(np3), rmove(np3), 
     +			 r0(np3), u0(np3), u(np3), ud(np3), 
     +			 rinfi(np3, np3), r(np3), nxyz(3, nxyzbin),
     + granddelta(nxyzbin*3,nxyzbin*3),F_ps(np3),
     +  constrvector(kconstr,np3),velocity_origin(np3),
     + velocity_pseudo(np3),F_conser(np3),F_random(np3),
     +  F_net(np3),distcheck(np),u_final(np3),randomforce_unproj(np3),
     + F_net2(np3),F_random2(np3),u2(np3)   )


         conhist=0
	  ghist=0
	  rghist=0
	  psihist=0
	  
	  
	  do i = 1, nxyzbin
		
	    nxyz(1, i) = 3 * (i-1) + 1
	    
		nxyz(2, i) = nxyz(1, i) + 1

		nxyz(3, i) = nxyz(2, i) + 1

	  end do  

	  boxlen = a * (4.0/3.0 * pi * np / phi) ** (1.0/3.0)
	  
	  boxlen=dg

	  vol = boxlen ** 3

c	  do i=1,np
	  
c	  r(nxyz(1,i))=(i-1)*fixdist
c	  r(nxyz(2,i))=0
c	  r(nxyz(3,i))=0
c	  end do
	  
c	  r(1)=0
c	  r(2)=0
c	  r(3)=0
c	  r(4)=fixdist
c	  r(5)=0
c	  r(6)=0
c	  r(7)=0.5*fixdist
c	  r(8)=0.5*fixdist*sqrt(3.0)
c	  r(9)=0
	  

	  do i=1,6
	   do j=1,6
	     delta(i,j)=0d0
	 	 if(i.eq.j)delta(i,j)=1d0
	   end do
	  end do
	  
	  do i=1,nxyzbin*3
	  do j=1,nxyzbin*3
	  Granddelta(i,j)=0d0
	  if(i .eq. j) Granddelta(i,j)=1d0
	  end do
	  end do
c



c
c	compute the self term (needs to be calculated only once)
c
	fx(0)=1.d0
	fx(1)=3.d0
	fx(2)=9.d0
	fx(3)=19.d0
	fx(4)=93.d0
	fx(5)=387.d0
	fx(6)=1197.d0
	fx(7)=5331.d0
	fx(8)=19821.d0
	fx(9)=76115.d0
	fx(10)=320173.d0
	fx(11)=1178451.d0

	fy(0)=1.d0
	fy(1)=1.5d0
	fy(2)=2.25d0
	fy(3)=7.375d0
	fy(4)=29.063d0
	fy(5)=70.594d0
	fy(6)=230.391d0
	fy(7)=700.336d0
	fy(8)=2226.629d0
	fy(9)=8694.131d0
	fy(10)=32889.478d0
	fy(11)=130304.1382d0

	g1=0.25
	g2=0.225
	g3=0.0268
c
c		call selftm(zi, ms)

          Fhw = 0.417
          
          
 
                  
c
	  dt_m = dt / real(m)

	  m_2 = real(m) / 2.d0

	  fac1 = 5.3052E+7 / a
        
        vis=1e-3
        
        fac1=6.0*pi*vis*a*1e-9

	  fac2 = 38.2749 * sqrt((tempr) / a / dt)
	  
	  fac3=5.3052e+10 /a

	  istart_dt = real(istart) * dt
	  	
c    **	Open file for writing simulation output   **

        do cycleindex=1, cyclenum
        
        if(cycleindex .lt. 10) then
        write(filetemp,'(i1)') cycleindex
        else if(cycleindex .lt. 100) then
        write(filetemp, '(i2)') cycleindex
        else 
        write(filetemp, '(i3)') cycleindex
        end if
        
         par_out = trim('sd_xyz')//trim(filetemp)//'.txt'
         
	  dss_out = trim('sd_dss')//trim(filetemp)//'.txt'
	  
	  open(144,file='distcheck.txt')
	  open(145,file=trim('sd_rodop')//trim(filetemp)//'.txt')
      
	  open ( unit = 20, file = par_out )

	  open ( unit = 30, file = dss_out )
	  
	  ori_vectornorm=sqrt(ori_vectorxtemp**2+ori_vectorytemp**2+
     +  ori_vectorztemp**2 )
	  
	  ori_vector1(1)=ori_vectorxtemp/ori_vectornorm
	  ori_vector1(2)=ori_vectorytemp/ori_vectornorm
        ori_vector1(3)=ori_vectorztemp/ori_vectornorm
        
        r_rod(1)=0
        r_rod(2)=0
        r_rod(3)=preheight
        
        do i=1,np
	  centerindex=(np+1)/2.0	  
	  r(nxyz(1,i))=(i-centerindex)*
     +	  fixdist*ori_vector1(1)+
     +  r_rod(1)
	  r(nxyz(2,i))=(i-centerindex)*
     +  fixdist*ori_vector1(2)+
     +   r_rod(2)
	  r(nxyz(3,i))=(i-centerindex)*
     +  fixdist*ori_vector1(3)+
     +   r_rod(3)
        enddo
        
        
        
      call getbodyframevector(ori_vector1,ori_vector2,ori_vector3)
c
c	  if(check.eq.'n') then

c          call writ_head(nstep, iprint, a, tempr, phi, dt, np)
          call writcn( t, r, nxyz )

c        endif
c
c	initialize arrays
c
	  rmove = 0.0
c
c    **   Now start the dynamics   **
c 
	  dxt = 0.0
	  dyt = 0.0
	  dzt = 0.0
	  t=0
		write(*,*) 'begin', cycleindex


	  do 50 l=1, nstep
	  
	    t=t+dt

	    t_min_isdt = t - istart_dt

c    ******************************************************************
c    ** call subroutines for calculating diffusion coefficients,     ** 
c    ** forces and random displacement terms                         **
c    ******************************************************************
	
          if(hydroflag .eq. 1) then
	    call grndrm(rinfi, rmove, rgrnd, t, r, nxyz)
c          write(*,*) l
          else
              rgrnd=0.0
          do i=1,np3
              rgrnd(i,i)=1
          end do
          end if
          

		do j=1,np3
		  
		  if(j.le.np3)then

		    r0(j) = r(j)
		    u0(j) = 0.0
		    u(j) = 0.0

		  endif

	      do k=1,np3

		    D(j,k)=rgrnd(j,k)
		  
		  end do
	    
		end do

          if(l .eq. 454) then
              write(*,*)
              end if
c
		call invert(D, np3)
          rgrnd=rgrnd*fac1
          D=D/fac1
         

        
        F_conser=0.0
	    call forces_random(F_random, r, nxyz,granddelta, 
     +     rgrnd,dt,fixdist,constrvector,randomforce_unproj )
          
          F_random=F_random*sqrt(2.0*kb*tempr/dt)*1e9
c         the random force is in nN          
c      F_random=0.0
      call calconserF(F_conser, r, nxyz,kappa,tempr,Bpw,ori_vector1)
c      F_conser=0.0
     
      
      call calfps(F_ps, r, nxyz, granddelta,fixdist)
      
c      F_ps=0.0
c      F_random=0.0
c      F_conser=0.0
c      F_conser(3)=1.0e-5
c      F_conser(6)=-1.0e-5
      call netforces(F_net, r, nxyz,granddelta, D, 
     +    rgrnd,dt,fixdist,
     +    constrvector,F_conser,F_random, F_ps) 
      
      
      	do j = 1, np3
          u0(j)=0.0
	      do k = 1, np3
	      
       u0(j)=u0(j)+D(j,k)*(F_net(k))
      
c         now the velocity is in nm/s      
c      velocity_pseudo(j)=velocity_pseudo(j)+D(j,k) *fac1
		  end do
c		  u0(j) = velocity_origin(j)!+velocity_pseudo(j)	  
c		  r(j) = r0(j) + u0(j) * dt/m
          end do	
          
         call calvelocity(r, nxyz, r_rod, u0, U0_rod,
     +    ori_vector1,ori_vector2,ori_vector3) 
      
         if(caldiffflag .eq. 1) then
c             F_net3=F_conser+F_random+F_ps
      call caldiff(ori_vector1,ori_vector2,ori_vector3, F_net,
     + U0_rod,r_rod,r,nxyz,tempr)
       end if
             u0_rod(4)=0
             u0_rod(3)=0
          rtemp(1)=u0_rod(1)*dt/m
        rtemp(2)=u0_rod(2)*dt/m
        rtemp(3)=u0_rod(3)*dt/m
c      change the unit to nm       
       
          thetatemp=u0_rod(4)*dt/m
        phitemp=u0_rod(5)*dt/m
        
c        thetatemp=0
c        phitemp=0
       
      
c     here is to transform body frame to lab frame
       r_rodmid(1)=r_rod(1)+rtemp(1)
          r_rodmid(2)=r_rod(2)+rtemp(2)         
            r_rodmid(3)=r_rod(3)+rtemp(3)
            
          ori_vector1mid(1)=ori_vector1(1)+phitemp*ori_vector2(1)+
     +  thetatemp*ori_vector3(1)
        ori_vector1mid(2)=ori_vector1(2)+phitemp*ori_vector2(2)+
     +  thetatemp*ori_vector3(2)
        ori_vector1mid(3)=ori_vector1(3)+phitemp*ori_vector2(3)+
     +  thetatemp*ori_vector3(3)
        


	 call getbodyframevector(ori_vector1mid,ori_vector2mid,ori_vector3mid)	
c       here i perform the calculation of the mass center translation velocity and mass centered rotation velocity        			
c       here i perform the test: weather the constraining force is vertial to the velocity
    	
		 	  
	  do i=1,np
	  
	  centerindex=(np+1)/2.0	  
	  r(nxyz(1,i))=(i-centerindex)*
     +	  fixdist*ori_vector1mid(1)+
     +  r_rodmid(1)
	  r(nxyz(2,i))=(i-centerindex)*
     +  fixdist*ori_vector1mid(2)+
     +   r_rodmid(2)
	  r(nxyz(3,i))=(i-centerindex)*
     +  fixdist*ori_vector1mid(3)+
     +   r_rodmid(3)
	  end do
           
           
        testsum=0.0
        do j=1,kconstr
        
        do i=1,np3
        
        testsum(j)=testsum(j)+constrvector(j,i)*velocity_origin(i)
        
        end do
        
        end do		
		
c		call calvelocity(r0,nxyz,D,F,fac3,tempr,l)	
		
c
c    **	call subroutines for calculating D and F based on midpoint coordinates  ** 
c               
             if(hydroflag .eq. 1) then
	    call grndrm(rinfi, rmove, rgrnd, t, r, nxyz)
          else
              rgrnd=0.0
          do i=1,np3
              rgrnd(i,i)=1
          end do
          end if
c         call convertrandom(F_random2, r, nxyz,granddelta, 
c     +    rgrnd,dt,fixdist,constrvector,randomforce_unproj )
c          F_random2=F_random2*sqrt(2.0*kb*tempr/dt)*1e9
          D=rgrnd
		call invert(D, np3)

          rgrnd=rgrnd*fac1
          D=D/fac1
        	call calconserF(F_conser, r, nxyz,kappa,tempr,Bpw,ori_vector1)
c          F_conser=0.0
      call calfps(F_ps, r, nxyz, granddelta,fixdist)
          F_ps=0.0
      call netforces(F_net, r, nxyz,granddelta, D, 
     +    rgrnd,dt,fixdist,
     +    constrvector,F_conser,F_random,F_ps)   
           call netforces(F_net2, r, nxyz,granddelta, D, 
     +    rgrnd,dt,fixdist,
     +    constrvector,F_conser,F_random2,F_ps)   
c
c    **	calculate new particle coordinates   **
c
	  u=0
c       u2=0
	    do j = 1, np3

	      do k = 1, np3

		    u(j) = u(j) + D(j,k)*( F_net(k))
c              u2(j) = u2(j) + D(j,k)*( F_net2(k))
c          now the velocity is in nm/s              
	
            end do
            ud(j)=m/2.0*(u(j)-u0(j))
            
c		  r(j) = r0(j) + (u0(j)+ud(j)) * dt

		end do
		
c          call calvelocity(r, nxyz, r_rod, ud, Udrift,
c     +     ori_vector1,ori_vector2,ori_vector3) 
c         here ud is for ori_vector1mid or ori_vector1          
c          call calvelocity(r, nxyz, r_rodmid, ud, Udrift2,
c     +     ori_vector1mid,ori_vector2mid,ori_vector3mid) 
c          call calvelocity(r, nxyz, r_rodmid, u2, U_rodmid,
c     +     ori_vector1mid,ori_vector2mid,ori_vector3mid)
c            U_rodfinal2=(U_rodmid-U0_rod)*1.0+U0_rod        
              u_final=u0+ud      
          call calvelocity(r, nxyz, r_rod, u_final, U_rodfinal,
     +     ori_vector1,ori_vector2,ori_vector3)            
              U_rodfinal(4)=0
              U_rodfinal(3)=0
c          do i=1,6
c       Udrift2(i)=m_2*(u_rodmid(i)-u0_rod(i))
c      end do
      
      rtemp(1)=(U_rodfinal(1))*dt
        rtemp(2)=(U_rodfinal(2))*dt
        rtemp(3)=(U_rodfinal(3))*dt
c      change the unit to nm       
       
c        Dis_random=sqrt(2.0*Dr*dt)*gasdev(iDummy)
c        Dis_drift= torque_vert2(j)*Dr*dt
      thetatemp=(U_rodfinal(4))*dt
        phitemp=(U_rodfinal(5))*dt
		 r_rod(1)=r_rod(1)+rtemp(1)
          r_rod(2)=r_rod(2)+rtemp(2)          
            r_rod(3)=r_rod(3)+rtemp(3)
c             thetatemp=0
c        phitemp=0
       
c          phitemp=0
c          thetatemp=0
       
       ori_vector1(1)=ori_vector1(1)+phitemp*ori_vector2(1)+
     +  thetatemp*ori_vector3(1)
        ori_vector1(2)=ori_vector1(2)+phitemp*ori_vector2(2)+
     +  thetatemp*ori_vector3(2)
        ori_vector1(3)=ori_vector1(3)+phitemp*ori_vector2(3)+
     +  thetatemp*ori_vector3(3)
        
c      now normalize it
        
      call getbodyframevector(ori_vector1,ori_vector2,ori_vector3)	
        
        
        
c        ori_vector1(1,j)=1
c       ori_vector1(2,j)=0
c       ori_vector1(3,j)=0
        phi=atan2(real(ori_vector1(2)),real(ori_vector1(1)))
       dtemp=sqrt(ori_vector1(2)**2+ori_vector1(1)**2)
       theta=atan2(real(dtemp),real(ori_vector1(3)))
       
        do i=1,np
	  
	  centerindex=(np+1)/2.0	  
	  r(nxyz(1,i))=(i-centerindex)*
     +	  fixdist*ori_vector1(1)+
     +  r_rod(1)
	  r(nxyz(2,i))=(i-centerindex)*
     +  fixdist*ori_vector1(2)+
     +   r_rod(2)
	  r(nxyz(3,i))=(i-centerindex)*
     +  fixdist*ori_vector1(3)+
     +   r_rod(3)
      end do
      
       
		   testsum=0.0
        do j=1,kconstr
        
        do i=1,np3
        
        testsum(j)=testsum(j)+constrvector(j,i)*(u0(j)+ud(j))
        
        end do
        
        end do	
		
		

	
		if( (l.gt.istart) .and. (mod(l,iprint).eq.0)) then
		
		

      	call writcn(t_min_isdt, r, nxyz )
      	
      	write(145,'(f17.4,11f12.5)') l*dt, r_rod(1)/a,r_rod(2)/a,
     + 	r_rod(3)/a,ori_vector1(1),ori_vector1(2),ori_vector1(3),
     + phi,theta,phi*180.0/3.1415926,
     + theta*180.0/3.1415926, min(r(nxyz(3,1)),r(nxyz(3,np)))/a
     
       
        
        end if
        
        
        

50	  continue

        close(20)
        close(30)
        close(145)
        end do
c
        
      	
      	
           
           
           

	  deallocate ( F, D, rgrnd, choldiag, randisp, rmove, 
     +			   r0, u0, u, ud, rinfi, r, nxyz )

	  stop
	  end
c
c
