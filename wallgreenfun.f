      function wallgnfunone(i,j,rij,rdist,y3,delta)
      
      
      double precision wallgnfunone,rij(3),rdist,y3
      integer i,j
       double precision delta(6,6)
      
      
      wallgnfunone=-(delta(i,j)/rdist+rij(i)*rij(j)/rdist**3)
     
      return
      
      end
      
      function wallgnfuntwo(i,j,rij,rdist,y3,delta)
      
      
      double precision wallgnfuntwo,rij(3),rdist,y3
      integer i,j
      double precision delta(6,6),line1,line2,line3,line4
      
      
      wallgnfuntwo=y3**2*(delta(i,j)*2.0/rdist**3-
     +  6.0/rdist**5*rij(i)*rij(j))
     + -2*y3**2*(2.0/rdist**3*delta(i,j)*delta(j,3)-
     + 6.0/rdist**5*rij(i)*rij(j)*delta(j,3))
     
c      wallgnfuntwo=wallgnfuntwo
     
      line1=y3**2*delta(i,j)*2.0/rdist**3
      line2=-y3**2*6.0/rdist**5*rij(i)*rij(j)
      line3=-2*y3**2*2.0/rdist**3*delta(i,j)*delta(j,3)
      line4=2*y3**2*6.0/rdist**5*rij(i)*rij(j)*delta(j,3)
      return
      
      end   
      
      
      function wallgnfunthree(i,j,rij,rdist,y3,delta)
      
      
      double precision wallgnfunthree,rij(3),rdist,y3
      integer i,j
      double precision delta(6,6),line1,line2,line3,line4,
     + line5,line6,line7,tempsum
      
c      wallgnfunthree=2.0/rdist**3*delta(i,j)*y3*rij(3)-
c     + 2.0*y3/rdist**3*(delta(j,3)*rij(i)+delta(i,3)*rij(j))+
c     + 6.0/rdist**5*rij(i)*rij(j)*y3*rij(3)-
c     + 4.0/rdist**3*delta(i,j)*delta(j,3)*y3*rij(3)+
c     + 4.0*(y3)/rdist**3*(delta(3,j)*rij(i)+
c     + delta(j,3)*delta(i,3)*rij(j))-
c     + 12.0/rdist**5*delta(3,j)*rij(i)*rij(j)*y3*rij(3)
     
c      wallgnfunthree=-2.0*y3*(delta(i,j)-2*delta(3,i)*delta(3,j))*
c     + (-delta(3,j)*rij(i)/rdist**3+delta(3,i)*rij(j)/rdist**3+
c     + rij(3)*delta(i,j)/rdist**3-3.0*rij(3)*rij(i)*rij(j)/rdist**5)
     
     
      wallgnfunthree=2.0*y3*
     + (-delta(i,j)*rij(3)/rdist**3+3.0*rij(3)*rij(i)*rij(j)/rdist**5+
     +  delta(i,3)*rij(j)/rdist**3-delta(3,j)*rij(i)/rdist**3)-
     + 4.0*y3*delta(j,3)*(-delta(i,3)*rij(3)/rdist**3+
     + 3.0*rij(i)*rij(3)**2/rdist**5+
     + delta(i,3)*rij(3)/rdist**3-rij(i)/rdist**3)
     
     
  

c       according to the reference of Blake, this should be -2h*...
c       However, in Brady's paper, this is 2h     
       line1=2.0/rdist**3*delta(i,j)*y3*rij(3)
       line2=-2.0*y3/rdist**3*(delta(j,3)*rij(i)+delta(i,3)*rij(j))
       line3=6.0/rdist**5*rij(i)*rij(j)*y3*rij(3)
       line4=-4.0/rdist**3*delta(i,j)*delta(j,3)*y3*rij(3)
       line5=4.0*(y3)/rdist**3*(delta(3,j)*rij(i))
       line6=4.0*(y3)/rdist**3*delta(j,3)*delta(i,3)*rij(j)
       line7=-12.0/rdist**5*delta(3,j)*rij(i)*rij(j)*y3*rij(3)
       
       tempsum=line1+line2+line3+line4+line5+line6+line7
     
      return
      
      end
      
      
      function sumwallgnfuns(i,j,rij,rdist,y3,delta)
      double precision wallgnfunthree,wallgnfunone,wallgnfuntwo,
     + rij(3),rdist,y3, sumwallgnfuns
      integer i,j
      double precision delta(6,6)
      
      sumwallgnfuns=wallgnfunone(i,j,rij,rdist,y3,delta)+
     + wallgnfuntwo(i,j,rij,rdist,y3,delta)+
     + wallgnfunthree(i,j,rij,rdist,y3,delta)
      
      
        return
        
        end