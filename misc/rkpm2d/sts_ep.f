c
      subroutine calsts_ep(disp,d_disp,
     &                     sts_CHY,sts_PK1,
     &                     pstn,ys_k,effsts,
     &                     mgk,ik)
c
c*****************************************************************
c
c*** update the stress status of each Gauss point 
c    for large deformation rate-independent
c    elasto-plastic problem (radial return mapping)
c      
c*** input parameters:
c    Lgp
c    Lmgp
c    d_disp
c    sts_CHY: of LAST time step at this Gauss point ( Cauchy stress )
c    pstn:   of LAST time step at this Gauss point
c    ys_k:   yield shear stress of LAST time step at this Gauss point
c    ik:     the # of current gauss point
c
c*** output parameters:
c    sts_CHY: the Cauchy stress matrix of CURRENT time step
c    sts_PK1: the PK1 stress
c    pstn:   of CURRENT time step
c    ys_k:   of CURRENT
c      
c
c
c*******************************************************************
c
      implicit none
      include 'parameter.h'

	! I&O parameters
      integer mgk,ik
      real*8 disp(2,maxNumnp),d_disp(2,maxNumnp)
      real*8 sts_CHY(3,3),sts_PK1(3,3)

      real*8 pstn,ys_k,effsts
      
      	! local arrays
      real*8 sts_last(3,3),sts(3,3),devsts(3,3)
      real*8 d_sts_trial(3,3),el_sts_trial(3,3),
     &       sts_trial(3,3),devsts_trial(3,3)
      real*8 ddu_dx(3,3),d_stn_sym(3,3),d_stn_anti(3,3)
      real*8 unit_normal(3,3)		! unit normal of yield surface
      					! i.e. unit normal of devsts tensor
                                        !
        ! local vars
      real*8 F_matx(3,3),Finv_matx(3,3),detF
      real*8 Fhf_matx(3,3),Fhfinv_matx(3,3),detFhf,
     &       detTmp
      real*8 pstn_last,ys_k_last
      real*8 sqrt23,
     &       vlength_devtrial,g_val,grad_g
      real*8 h_prime, ys_k_prime
      real*8 d_lmda, d_h
      real*8 pstn_iter
      real*8 ys_k_iter
      real*8 epis_yieldcond,epis_detF
      real*8 trace_ststr,trace_sts
c      
      integer jp,ir,i,j,k
      integer mLoop
      integer iiter,maxIIter
      integer jiter,maxJiter
c
      real*8  tmpLength,tmp
      real*8  d11,d12,d21,d22,d33,dstn1,dstn2,dstn3
c
      real*8  tmp_scale
      integer n_three
      real*8  unit_matx(3,3),Q_matx(3,3),Qt_matx(3,3),
     &        tmp1_matx(3,3),tmp2_matx(3,3),tmp3_matx(3,3),
     &        tmp4_matx(3,3),tmp5_matx(3,3),tmp6_matx(3,3)
      integer i_rotation_trial_sts_method
      integer i_smalldef		! =1 : small deformation
      					!  0 : large deformation

      real*8  shpk_dsmallx,shpk_dsmally

      real*8  tmpmean,tmpsts(3,3),tmpratio
      logical do_radial_return

	! common blocks
      real*8 el_E,el_V,el_G,el_K
      real*8 pl_k0,pl_H,pl_EP  
c
      integer lnods(maxNode,maxElem),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
c
      real*8 shpk(mnsch,maxGP),shpkdx(mnsch,maxGP),
     &       shpkdy(mnsch,maxGP),shpn(mnsch,maxNumnp)
c
      common /elastic/ el_E, el_V, el_G, el_K
      common /plastic/ pl_k0, pl_H, pl_EP
      common /shapeK/ shpk,shpkdx,shpkdy,shpn
      common /connect/ Lnp,Lmnp,Lgp,Lmgp,lnods
c
c-----------------------------------------------------------------
c
      i_rotation_trial_sts_method=1
      					! to determine the method to cal the rotational part
      					! of trial stress
                                        ! == 1: use Hughes method(1980)
       					!    2: use the forward euler 
    					    
      i_smalldef=0			
      					! == 1: small deformation
                                        !    0: large deformation
                                            
      sqrt23=sqrt(2./3.)
      n_three=3
      
      do i=1,3
         do j=1,3
            sts_last(i,j)=sts_CHY(i,j)
         end do
      end do
      
      
      pstn_last=pstn
      ys_k_last=ys_k
      
      	        ! cal ddu_dx, where,
		!
		!
		!
                !			d(delta_u  )
                !                                i  
		!    ddu_dx(i,j)= -----------------------                
                !                           dx
		!                             j
		!
		!		  d(delta_u  )      d X  
                !                           i          k
                !               = -------------- -----------
                !                      dX           d x
                !                        k             j
                !
                !

         mLoop = Lgp(ik)
         
		! cal the deformation gradient
         do i=1,3
            do j=1,3       
               F_matx(i,j)=0.
               Fhf_matx(i,j)=0.
            enddo
         enddo
         do ir = 1, mLoop
            jp = Lmgp(ir,ik)
            F_matx(1,1) = F_matx(1,1)
     &                 + shpkdx(ir,ik)*disp(1,jp)
            F_matx(1,2) = F_matx(1,2) 
     &                 + shpkdy(ir,ik)*disp(1,jp)
            F_matx(2,1) = F_matx(2,1) 
     &                 + shpkdx(ir,ik)*disp(2,jp)
            F_matx(2,2) = F_matx(2,2) 
     &                 + shpkdy(ir,ik)*disp(2,jp)

            Fhf_matx(1,1) = Fhf_matx(1,1)
     &                 + shpkdx(ir,ik)*(disp(1,jp)-d_disp(1,jp)*0.5)
            Fhf_matx(1,2) = Fhf_matx(1,2) 
     &                 + shpkdy(ir,ik)*(disp(1,jp)-d_disp(1,jp)*0.5)
            Fhf_matx(2,1) = Fhf_matx(2,1) 
     &                 + shpkdx(ir,ik)*(disp(2,jp)-d_disp(2,jp)*0.5)
            Fhf_matx(2,2) = Fhf_matx(2,2) 
     &                 + shpkdy(ir,ik)*(disp(2,jp)-d_disp(2,jp)*0.5)
          
         enddo
         do i=1,3
            F_matx(i,i)=F_matx(i,i)+1.
            Fhf_matx(i,i)=Fhf_matx(i,i)+1.
         enddo
         
		!!!!! reset F_matx to be unit matx only
         if ( i_smalldef.eq.1) then       
            do i=1,3
               do j=1,3
                  F_matx(i,j)=0.
                  Fhf_matx(i,j)=0.
               enddo
               F_matx(i,i)=1.
               Fhf_matx(i,j)=1.
            enddo
         endif
         
         
         	! cal the F^(-1)
         call invers_matx33(F_matx,Finv_matx,detF)
         call invers_matx33(Fhf_matx,Fhfinv_matx,detFhf)
         
         do i = 1, 3
            do j = 1, 3
               ddu_dx(i,j) = 0.
            enddo
         enddo

         if ( i_smalldef.eq.1) then	  
         		! small deformation
            do ir = 1, mLoop
               jp = Lmgp(ir,ik)
               ddu_dx(1,1) = ddu_dx(1,1)
     &                    + shpkdx(ir,ik)*d_disp(1,jp)
               ddu_dx(1,2) = ddu_dx(1,2)
     &                    + shpkdy(ir,ik)*d_disp(1,jp)
               ddu_dx(2,1) = ddu_dx(2,1)
     &                    + shpkdx(ir,ik)*d_disp(2,jp)
               ddu_dx(2,2) = ddu_dx(2,2)
     &                    + shpkdy(ir,ik)*d_disp(2,jp)
            enddo
            ddu_dx(3,3) = 0.
         else   
         		! large deformation
            do ir = 1, mLoop
               jp = Lmgp(ir,ik)
               
               shpk_dsmallx=shpkdx(ir,ik)*Fhfinv_matx(1,1)
     &                     +shpkdy(ir,ik)*Fhfinv_matx(2,1)  
     
               shpk_dsmally=shpkdx(ir,ik)*Fhfinv_matx(1,2)
     &                     +shpkdy(ir,ik)*Fhfinv_matx(2,2)  
     
               
               ddu_dx(1,1) = ddu_dx(1,1)
     &                    + shpk_dsmallx*d_disp(1,jp)
               ddu_dx(1,2) = ddu_dx(1,2)
     &                    + shpk_dsmally*d_disp(1,jp)
               ddu_dx(2,1) = ddu_dx(2,1)
     &                    + shpk_dsmallx*d_disp(2,jp)
               ddu_dx(2,2) = ddu_dx(2,2)
     &                    + shpk_dsmally*d_disp(2,jp)
            enddo
            ddu_dx(3,3) = 0.            
         endif		! end to cal ddu_dx for small deformation
                        ! or large deformation
	                        
        
         ! predicting ( trial stress )

      do i=1,3
         do j=1,3
            d_stn_sym(i,j)  = 0.0
            d_stn_anti(i,j) = 0.0
         end do
      end do

      do i=1,2
         do j=1,2
                ! symmetric part
            d_stn_sym(i,j)=0.5* ( ddu_dx(i,j)+ddu_dx(j,i) )
                ! anti-symmetric part
            d_stn_anti(i,j)=0.5* ( ddu_dx(i,j)-ddu_dx(j,i) )
         end do
      end do
c      
      d_stn_sym(3,3)   = ddu_dx(3,3)
      d_stn_anti(3,3)  = 0.0
c      
      if ( i_smalldef.eq.1) then
         do i=1,3
            do j=1,3
               d_stn_anti(i,j) = 0.0
            enddo
         enddo
      endif
c      
c         
c    !** elastic stress increment for trial stress
c
c
      do i=1,3
         do j=1,3
            el_sts_trial(i,j) = 0.0
            d_sts_trial(i,j)  = 0.0
         enddo
      enddo
      
      	! for plain stain
c      E0=el_E/(1-el_V**2)
c      V0=el_V/(1-el_V)
c      D0=E0/(1-V0*V0)
c      
c      D1111=D0
c      D1122=D0*V0
c      D2211=D0*V0
c      D2222=D1111
c      D1212=D0*(1-V0)/2
c      
c      el_sts_trial(1,1)= D1111*d_stn_sym(1,1)+D1122*d_stn_sym(2,2)
c      el_sts_trial(2,2)= D2211*d_stn_sym(1,1)+D2222*d_stn_sym(2,2)
c      el_sts_trial(1,2)=2*D1212*d_stn_sym(1,2)
c      el_sts_trial(2,1)=el_sts_trial(1,2)
c      el_sts_trial(3,3)=el_V*(el_sts_trial(1,1)+el_sts_trial(2,2) )
	
	   ! for plain strain     	
c       vkapa=3.-4.*el_V
c       tmp=(3.-vkapa)/(2.*(vkapa-1.) )
c       trace=d_stn_sym(1,1)+d_stn_sym(2,2)
c       el_sts_trial(1,1)=2.*el_G*(d_stn_sym(1,1)+tmp*trace)
c       el_sts_trial(1,2)=2.*el_G*d_stn_sym(1,2)
c       el_sts_trial(2,1)=el_sts_trial(1,2)
c       el_sts_trial(2,2)=2.*el_G*(d_stn_sym(2,2)+tmp*trace)
c       el_sts_trial(3,3)=el_V*(el_sts_trial(1,1)+el_sts_trial(2,2) )

 
c 	 vlmda=(2.*el_G*el_v)/(1.-2.*el_V)
c        trace=d_stn_sym(1,1)+d_stn_sym(2,2)
c        el_sts_trial(1,1)=2.*el_G*d_stn_sym(1,1)+vlmda*trace
c        el_sts_trial(1,2)=2.*el_G*d_stn_sym(1,2)
c        el_sts_trial(2,1)=el_sts_trial(1,2)
c        el_sts_trial(2,2)=2.*el_G*d_stn_sym(2,2)+vlmda*trace
c        el_sts_trial(3,3)=2.*el_G*d_stn_sym(3,3)+vlmda*trace

      tmp=el_E*(1.-el_V) / ( (1.+el_V)*(1.-2.*el_V) )
      d11=tmp
      d12=tmp*el_V / (1.-el_V)
      d21=d12
      d22=tmp
      d33=tmp*(1.-2.*el_V )/ (2.*(1.-el_V) )
      dstn1=d_stn_sym(1,1)
      dstn2=d_stn_sym(2,2)
      dstn3=d_stn_sym(1,2)
      
      el_sts_trial(1,1)=d11*dstn1+d12*dstn2
      el_sts_trial(2,2)=d21*dstn1+d22*dstn2
      el_sts_trial(1,2)=d33*(2.*dstn3)
      el_sts_trial(2,1)=el_sts_trial(1,2)
      el_sts_trial(3,3)=el_V*( el_sts_trial(1,1)+el_sts_trial(2,2) )
       
		! eval the rotational part of trial stress
       
      if  (i_rotation_trial_sts_method.eq.1) then
c       
c    ! eval the trial stress via the method of Hughes(1980) 
c
         do i=1,3
            do j=1,3
               unit_matx(i,j) = 0.0
               Q_matx(i,j)    = 0.0
               Qt_matx(i,j)   = 0.0
               tmp1_matx(i,j) = 0.0
               tmp2_matx(i,j) = 0.0
               tmp3_matx(i,j) = 0.0
               tmp4_matx(i,j) = 0.0
               tmp5_matx(i,j) = 0.0
               tmp6_matx(i,j) = 0.0
             end do
 
             unit_matx(i,i)=1.
         end do
c
         tmp_scale=0.5
         call scalemul_matx(d_stn_anti,tmp_scale,tmp1_matx,
     &                     n_three,n_three)
                                          ! tmp1 = 1/2 * W
 
         call minus_matx(unit_matx,tmp1_matx,tmp2_matx,n_three,n_three)
                                          ! tmp2= I - 1/2 * W
 
         call invers_matx33(tmp2_matx,tmp3_matx,detTmp)
                                          ! tmp3= inv (I-1/2 * W)
 
         call mul_matx(tmp3_matx,d_stn_anti,tmp4_matx,
     &                 n_three,n_three,n_three)
                                                   ! tmp4= inv(i-1/2 * W) * W
 
         call add_matx(tmp4_matx,unit_matx,Q_matx,n_three,n_three)
                                                   ! Q_matx
 
           ! cal the trial stress
         call trans_matx(Q_matx, Qt_matx, n_three, n_three)
 
         call mul_matx(Q_matx, sts_last, tmp5_matx,
     &                 n_three,n_three,n_three)
                                                   ! tmp5= Q*tao
 
         call mul_matx(tmp5_matx, Qt_matx, tmp6_matx,
     &                 n_three,n_three,n_three)
                                                   ! tmp6=Q* tao* Qt
 
         do i=1,3
            do j=1,3
               sts_trial(i,j)=el_sts_trial(i,j)+tmp6_matx(i,j)
            enddo
         enddo
      
      elseif (i_rotation_trial_sts_method.eq.2) then
	      ! rotation part of Jaumann rate increment
          
         if ( 1.eq.0) then     	! consider Jaumann rate increment
         
            do i=1,3
               do j=1,3
                  d_sts_trial(i,j)=el_sts_trial(i,j)
                  do k=1,3
                     d_sts_trial(i,j)=d_sts_trial(i,j)
     &                     +sts_last(i,k)*d_stn_anti(j,k)
     &                     +sts_last(j,k)*d_stn_anti(i,k)
                  end do
               end do
            end do
         
         else
            do i=1,3
               do j=1,3
                 d_sts_trial(i,j)=el_sts_trial(i,j)
               end do
           end do
           
         endif

 
           ! get the trial stress
          do i=1,3
            do j=1,3
               sts_trial(i,j)=sts_last(i,j)+d_sts_trial(i,j)
            end do
          end do
      else
          write(*,*) 'bad value of i_rotation_trial_sts_method'
          write(*,*) 'i_rotation_trial_sts_method=',
     &               i_rotation_trial_sts_method
      endif
          
      
      trace_ststr=0.
      do i=1,3
         trace_ststr=trace_ststr+sts_trial(i,i)
      enddo      
      do i=1,3
         do j=1,3
            devsts_trial(i,j)=sts_trial(i,j)
         enddo
         devsts_trial(i,i)=sts_trial(i,i)-trace_ststr/3.
      enddo
      
      vlength_devtrial=0.  
      do i=1,3
        do j=1,3
           vlength_devtrial=vlength_devtrial+devsts_trial(i,j)**2
        end do
      end do
      vlength_devtrial=sqrt(vlength_devtrial)
      
      
        ! radial return algorithm
      g_val=-sqrt23*ys_k+vlength_devtrial
      iiter=0
      maxIIter=5
      
                 ! the tolerance of the iteration to eval the d_lmda
      epis_yieldcond=1.e-3

  
      if ( g_val/pl_k0 .le. epis_yieldcond )  then 
             	! elastic loading or unloading, still inside the yield surface
         
         do_radial_return=.false.
                
         do i=1,3
            do j=1,3
               sts(i,j)=sts_trial(i,j)
            end do
         end do
         
c         if (g_val/pl_k0 .ge. -0.5 ) then
c         
c                     ! add elastic volume change to update the total stress values
c            do i=1,3
c               do j=1,3
c                  sts(i,j)=devsts_trial(i,j)
c               enddo
c            enddo                     
c                     
c                     
c            trace_dstn=0.
c            do i=1,3
c               trace_dstn=trace_dstn+d_stn_sym(i,i)
c            end do
c           
c            do i=1,3
c               sts(i,i)=sts(i,i)+trace_laststs/3.+el_K*trace_dstn
c            end do         
c
c        endif
         

      else 	
      		! outside the yield surface, need use radial return mapping  
                
         do_radial_return=.true.
                
          iiter=iiter+1   
          
c          write(*,'(1x,a,i5,A,4e11.3)' ) 
c     &              'iiter=',iiter, 
c     &              ' sts_trial= ',sts_trial(1,1),sts_trial(1,2),
c     &              sts_trial(2,2)

          if (iiter.gt.maxIIter) then
             write(*,*) 'Iiter exceed limit. Divergent.stop'
             stop
          endif
      
           ! compute the unit normal  ( assume the back stress is zero)
         do i=1,3
            do j=1,3
               unit_normal(i,j)=devsts_trial(i,j)/vlength_devtrial
            end do
         end do   
        
         
           ! find the d_lmda by iteration
           
           ! h_prime is the derivative of back stress module vs. plastic strain
           ! yield_shrsts_prime is the derivative of the yield_k vs. plastic strain
           !   ( assume isotropic linear hardening  )
         h_prime=0.
c         ys_k_prime=pl_H
         ys_k_prime=pl_EP
           
           ! init value of iteration
         d_lmda=0.	  
         pstn_iter=pstn_last
         ys_k_iter=ys_k_last
         
         jiter=0
         maxjIter=10
         do while ( abs(g_val)/pl_k0 .gt. epis_yieldcond )
         
            jiter=jiter+1
            
            if ( jiter.gt.maxjIter) then
               write(*,*) 'jIter exceed limit. Divergent. stop'
               write(*,*) 'maxjIter=',maxjIter
               stop
            endif
            
         	! gradient of the g function
            grad_g=-2.*el_G*
     &        (1. + (h_prime+ys_k_prime)/(3.*el_G) )
      
            d_lmda=d_lmda-g_val/grad_g

            pstn_iter=pstn_last+sqrt23*d_lmda

c            ys_k_iter=pl_k0+pstn_iter*pl_H

			! linear harding
	    ys_k_iter=pl_k0+pstn_iter*pl_EP
c            ys_k_iter=ys_k_last+sqrt23*d_lmda*pl_EP

            d_h=0.	! assume isotropic hardening
            
            g_val=-sqrt23*ys_k_iter+vlength_devtrial
     &            -(2.*el_G*d_lmda+sqrt23*d_h)


         end do
         
         pstn=pstn_iter
         ys_k=ys_k_iter
         


         
           ! get the updated stresses
         do i=1,3
            do j=1,3
                sts(i,j)=sts_trial(i,j)
     &         -d_lmda*(2.*el_G*unit_normal(i,j))
            end do
         end do

      end if	! if use radial return

      if ( 1.eq.0) then     	! consider Jaumann rate increment
         
            do i=1,3
               do j=1,3
                  do k=1,3
                     sts(i,j)=sts(i,j)
     &                     +sts_last(i,k)*d_stn_anti(j,k)
     &                     +sts_last(j,k)*d_stn_anti(i,k)
                  end do
               end do
            end do
                        
       endif  

		! to test if the yield consistent condition satisfied         
                
      if ( do_radial_return  ) then    
          
         tmpmean=(sts(1,1)+sts(2,2)+sts(3,3))/3.
         do i=1,3
            do j=1,3
               tmpsts(i,j)=sts(i,j)
            enddo
            tmpsts(i,i)=tmpsts(i,i)-tmpmean
         enddo
 
         tmpLength=0.
         do i=1,3
            do j=1,3
               tmpLength=tmpLength+tmpsts(i,j)**2
            enddo
         enddo
         tmpLength=sqrt(tmpLength)
         tmpratio=tmpLength/(sqrt23*ys_K)
         if ( abs( tmpratio -1. ). gt. 1.e-2 ) then
             write(*,*)
     &        'abs( tmpLength/(sqrt23*ys_K) -1. ). gt. 1.e-2'
             write(*,*) 'tmpratio=',tmpratio
             write(*,*) 'stop'
             stop
         endif
 
      endif


      do i=1,3
        do j=1,3
            sts_CHY(i,j)=sts(i,j)
         end do
      end do
c
      call mul_matx(Finv_matx,sts_CHY,sts_PK1,n_three,n_three,n_three)
c
      epis_detF=1.e-6
      if ( detF.lt.epis_detF ) then
         write(*,*) 'detF < epis_detF in sts_ep.f'
         write(*,'(1x,a,i5,a,e12.4)') 'ik=',ik,' detF=',detF
         stop
      endif
     
      call scalemul_matx(sts_PK1,detF,sts_PK1,n_three,n_three)


      trace_sts=sts(1,1)+sts(2,2)+sts(3,3)
      do i=1,3
         do j=1,3
            devsts(i,j)=sts(i,j)
         enddo
         devsts(i,i)=sts(i,i)-trace_sts/3.
      enddo

      effsts=0.
      do i=1,3
         do j=1,3
            effsts=effsts+devsts(i,j)*devsts(i,j)
         enddo
      enddo
      effsts=sqrt(1.5*effsts)
      
      end
