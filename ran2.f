	FUNCTION ran2(idummy)

c    ******************************************************************
c    ** Function to generate uniform random deviates between 0 and 1 **
c    **					  (Numerical Recipes).                     **
c    ******************************************************************

        integer idummy,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        double precision ran2, AM, EPS, RNMX
        parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     &	  	  IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     &		 IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        integer idum2,j,k,iv(NTAB),iy
        save iv,iy,idum2
        data idum2/123456789/, iv/NTAB*0/, iy/0/

        if (idummy.le.0) then
           idummy=max(-idummy,1)
           idum2=idummy
           do j=NTAB+8,1,-1
              k=idummy/IQ1
              idummy=IA1*(idummy-k*IQ1)-k*IR1
              if (idummy.lt.0) idummy=idummy+IM1
              if (j.le.NTAB) iv(j) = idummy
           end do
           iy=iv(1)
        endif
        k=idummy/IQ1
        idummy=IA1*(idummy-k*IQ1)-k*IR1
        if (idummy.lt.0) idummy=idummy+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2 = idum2+IM2
        j=1+iy/NDIV
        iy=iv(j)-idum2
        iv(j)=idummy
        if (iy.lt.1) iy=iy+IMM1
        ran2=min(AM*iy,RNMX)

        return
        end