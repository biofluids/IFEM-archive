      subroutine calsts_el(disp,stress,mgk,ik)
c
c***********************************************************
c
c*** update the stress status of each Gauss point 
c    for elastic problem
c      
cc*** input parameters:
c    Lgp
c    Lmgp
c    disp
c    stress: of LAST time step at this Gauss point
c    ik:     the # of current gauss point
c
c*** output parameters:
c    stress: the Cauchy stress matrix (transpose) of CURRENT time step
c    pstn:   of CURRENT time step
c    ys_k:   of CURRENT
c
c
c
      implicit none
      include 'parameter.h'
c
	! I&O parameters
c
      integer mgk,ik
      integer Lgp(maxGP),Lmgp(mnsch,maxGP),
     &        lnods(maxNode,maxElem),
     &        Lmnp(mnsch,maxNumnp),
     &        Lnp(maxNumnp)
      real*8 shpk(mnsch,maxGP),shpkdx(mnsch,maxGP),
     &       shpkdy(mnsch,maxGP),shpn(mnsch,maxNumnp),
     &       disp(2,maxNumnp)
      real*8 stress(3,3)
      
      	! local arrays
      real*8 du_dx(3,3),stn_sym(3,3),sts(3,3)

        ! local vars
      integer jp,ir,i,j,k,mLoop

      real*8 tmp,
     &       d11,d12,d21,d22,d33,stn1,stn2,stn3


	! common blocks
      real*8 el_E,el_V,el_G,el_K,
     &       pl_k0,pl_H,pl_EP
c
      common /shapeK/shpk,shpkdx,shpkdy,shpn
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /elastic/ el_E, el_V, el_G, el_K      
      common /plastic/ pl_k0, pl_H, pl_EP

      
      	        ! cal du_dx, where,
		!
		!
		!
                !			    du 
                !                             i  
		!    du_dx(i,j)= -----------------------
                !                           dx
		!                             j

         do i = 1, 3
            do j = 1, 3
               du_dx(i,j) = 0.
            enddo
         enddo

         mLoop = Lgp(ik)
       
         do ir = 1, mLoop
            jp = Lmgp(ir,ik)
            du_dx(1,1) = du_dx(1,1)
     &                 + shpkdx(ir,ik)*disp(1,jp)
            du_dx(1,2) = du_dx(1,2) 
     &                 + shpkdy(ir,ik)*disp(1,jp)
            du_dx(2,1) = du_dx(2,1) 
     &                 + shpkdx(ir,ik)*disp(2,jp)
            du_dx(2,2) = du_dx(2,2) 
     &                 + shpkdy(ir,ik)*disp(2,jp)
         enddo
         du_dx(3,3) = 0.
        
         ! predicting ( trial stress )

      do i=1,3
         do j=1,3
            stn_sym(i,j)=0.
         end do
      end do

      do i=1,2
         do j=1,2
                ! symmetric part
            stn_sym(i,j)=0.5* ( du_dx(i,j)+du_dx(j,i) )
         end do
      end do      
      stn_sym(3,3)=du_dx(3,3)
         
      do i=1,3
         do j=1,3
            sts(i,j)=0.
         end do
      end do
      

      tmp=el_E*(1.-el_V) / ( (1.+el_V)*(1.-2.*el_V) )
      d11=tmp
      d12=tmp*el_V / (1.-el_V)
      d21=d12
      d22=tmp
      d33=tmp*(1.-2.*el_V )/ (2.*(1.-el_V) )
      stn1=stn_sym(1,1)
      stn2=stn_sym(2,2)
      stn3=stn_sym(1,2)
      
      sts(1,1)=d11*stn1+d12*stn2
      sts(2,2)=d21*stn1+d22*stn2
      sts(1,2)=d33*(2.*stn3)
      sts(2,1)=sts(1,2)
      sts(3,3)=el_V*( sts(1,1)+sts(2,2) )
 

c     !!! stress( , ) is the transpose of tau( , ) !!!

      do i=1,3
         do j=1,3
         
            stress(i,j)=sts(j,i)
            
         end do
      end do

      end
