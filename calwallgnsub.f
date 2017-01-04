      subroutine calwalgnsub(i,j,r,rimage,nxyz,p1,p2,
     + imagerij,imagedist,finalresult)
      
      
      common / num_par / np
	common / np_3 / np3
	common / delta_fn / delta
	common / radius / a
      
      
      integer  n, np, n3, np3, nw3, i, j, k, l

	integer nxyz(3, np),p1,p2

      double precision r(np3),rimage(np3),imagedist,imagerij(3),a,
     + delta(6,6),green, lapxgreen, lapygreen, finalresult, 
     + imagerijdxplus(3),imagerijdxminus(3),imagedistdxplus,
     + imagedistdxminus,imagerijdyplus(3),imagerijdyminus(3),
     + imagedistdyplus,
     + imagedistdyminus,sumwallgnfuns,y3,y3tempplus,y3tempminus,
     + firstterm,secondterm,thirdterm,testgreen,testgreen2,
     + wallgnfunone,wallgnfunthree,wallgnfuntwo



      y3=r(nxyz(3,p2))
      
c     first calculate the green function itself
      green=0.0
      lapxgreen=0.0
      lapygreen=0.0
      testgreen=0.0
      testgreen2=0.0
      tempsum=0.0
      
      h=0.01

      green=sumwallgnfuns(i,j,imagerij,imagedist,y3,delta)
     
c      now calculate the lapacian of green fucntion to x=vector(x_alpha)

       do k=1,3
       imagerijdxplus(1)=r(nxyz(1,p1))+h*delta(1,k)-rimage(nxyz(1,p2))
       imagerijdxplus(2)=r(nxyz(2,p1))+h*delta(2,k)-rimage(nxyz(2,p2))
      imagerijdxplus(3)=r(nxyz(3,p1))+h*delta(3,k)-rimage(nxyz(3,p2))
      
      imagedistdxplus=sqrt(imagerijdxplus(1)**2+imagerijdxplus(2)**2+
     + imagerijdxplus(3)**2 )
      
       imagerijdxminus(1)=r(nxyz(1,p1))-h*delta(1,k)-rimage(nxyz(1,p2))
       imagerijdxminus(2)=r(nxyz(2,p1))-h*delta(2,k)-rimage(nxyz(2,p2))
      imagerijdxminus(3)=r(nxyz(3,p1))-h*delta(3,k)-rimage(nxyz(3,p2))
      
      imagedistdxminus=sqrt(imagerijdxminus(1)**2+imagerijdxminus(2)**2+
     + imagerijdxminus(3)**2 )

      lapxgreen=lapxgreen+1.0/(h**2)*(
     + sumwallgnfuns(i,j,imagerijdxplus,imagedistdxplus,y3,delta)
     + -2.0*sumwallgnfuns(i,j,imagerij,imagedist,y3,delta)
     + +sumwallgnfuns(i,j,imagerijdxminus,imagedistdxminus,y3,delta))
     
c      testgreen=1.0/(h**2)*(
c     + wallgnfunone(i,j,imagerijdxplus,imagedistdxplus,y3,delta)
c     + -2.0*wallgnfunone(i,j,imagerij,imagedist,y3,delta)
c     + +wallgnfunone(i,j,imagerijdxminus,imagedistdxminus,y3,delta))
     
c      tempsum=tempsum+testgreen
      end do
      
c      now calculate the lapacian of green fucntion to y=vector(x_beta)
c       how i choose the sign of the h, when increase the origin particle in 1 and 2 direction
c       the rimage in 1 and 2 direction increase too. However, when i increase the original 
c       particle in 3 or z direction, the rimage will decrease in 3 or z diection      
c       tempsum=0.0
     
       do k=1,3
      imagerijdyplus(1)=r(nxyz(1,p1))-(rimage(nxyz(1,p2))+h*delta(1,k))
      imagerijdyplus(2)=r(nxyz(2,p1))-(rimage(nxyz(2,p2))+h*delta(2,k))
      imagerijdyplus(3)=r(nxyz(3,p1))-(rimage(nxyz(3,p2))-h*delta(3,k))
      
      imagedistdyplus=sqrt(imagerijdyplus(1)**2+imagerijdyplus(2)**2+
     + imagerijdyplus(3)**2)
      
       imagerijdyminus(1)=r(nxyz(1,p1))+h*delta(1,k)-rimage(nxyz(1,p2))
       imagerijdyminus(2)=r(nxyz(2,p1))+h*delta(2,k)-rimage(nxyz(2,p2))
      imagerijdyminus(3)=r(nxyz(3,p1))-h*delta(3,k)-rimage(nxyz(3,p2))
      
      imagedistdyminus=sqrt(imagerijdyminus(1)**2+imagerijdyminus(2)**2+
     + imagerijdyminus(3)**2 )
     
      y3tempplus=y3+h*delta(3,k)
      y3tempminus=y3-h*delta(3,k)
      
      lapygreen=lapygreen+1.0/(h**2)*(
     +sumwallgnfuns(i,j,imagerijdyplus,imagedistdyplus,y3tempplus,delta)
     + -2.0*sumwallgnfuns(i,j,imagerij,imagedist,y3,delta)
     + +sumwallgnfuns(i,j,imagerijdyminus,imagedistdyminus,
     + y3tempminus,delta))
     
     
         testgreen3=1.0/(h**2)*(
     + wallgnfunthree(i,j,imagerijdyplus,
     +  imagedistdyplus,y3tempplus,delta)
     + -2.0*wallgnfunthree(i,j,imagerij,imagedist,y3,delta)
     + +wallgnfunthree(i,j,imagerijdyminus,imagedistdyminus,
     + y3tempminus,delta))
     
        testgreen2=1.0/(h**2)*(
     + wallgnfuntwo(i,j,imagerijdyplus,
     +  imagedistdyplus,y3tempplus,delta)
     + -2.0*wallgnfuntwo(i,j,imagerij,imagedist,y3,delta)
     + +wallgnfuntwo(i,j,imagerijdyminus,imagedistdyminus,
     + y3tempminus,delta))
     
      testgreen1=1.0/(h**2)*(
     + wallgnfunone(i,j,imagerijdyplus,
     +  imagedistdyplus,y3tempplus,delta)
     + -2.0*wallgnfunone(i,j,imagerij,imagedist,y3,delta)
     + +wallgnfunone(i,j,imagerijdyminus,imagedistdyminus,
     + y3tempminus,delta))
     
      tempsum=tempsum+testgreen
      end do
      
c       note 3.0a/4.0 comes from 1/8 pi mu*6pi mu a
c       note a^3/8.0 comes from 1/8 pi mu*6pi mu a* 1/6*a^2      
      firstterm=3.0*a/4.0*green
      secondterm=a**3/8.0*lapxgreen
      thirdterm=a**3/8.0*lapygreen
      
      finalresult=firstterm+secondterm+thirdterm
      

      
      end