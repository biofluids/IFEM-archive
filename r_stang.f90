subroutine r_stang
  use run_variables, only: ntsbout,its
  use solid_variables
  use r_common
!  use solid_fem_surface_traction
  implicit none

  real*8 :: x(3,nis)
  real*8 :: y(3,nis)
  real*8 toc(3,3)  !...Cauchy Green Deformation Tensor at t=0
  real*8 xto(3,3)  !...F    = dx/dX
  real*8 xot(3,3)  !...F^-1 = dX/dx
  real*8 xj(3,3)   !...
  real*8 xji(3,3)  !...
  real*8 rs(3)     !...Shapefunction value for each spatial direction at quadrature point (isoparametric)
  real*8 toxj(3,3)
  real*8 toxji(3,3)
  real*8 xmj(3),xmi(3),dxmj(3,6),ddxmj(3,6,6),obc(6,6),ocuu(6,6),ocup(6)
  real*8 xfrtem(6,6),tem(6),ten(6),ttm(6)

  integer :: ine,in,n,k,nu1,mu1,ip,jp

  real*8 :: w,wp,tot_vol_init,tot_vol_curr
  integer :: nos,ntem
  integer :: ixq,iyq,izq,isd
  real*8 :: det   !...Jacobian Determinante of Deformation Gradient
  real*8 :: todet !...Jacobian Determinante of Deformation Gradient at t=0
  real*8 :: wto
  real*8 :: ocpp

  write(*,*) " calculate internal + inertial forces (r_stang)"

  call r_common_allocate(ne_solid,nup,nis)


  predrf(1:3*nn_solid) = 0.0d0

  tot_vol_init = 0.0d0
  tot_vol_curr = 0.0d0

  element: do ine=1,ne_solid
 !...position
     do jp=1,nump
        xfp(jp,ine)=0.0d0

        xkup(1:3*nis,jp,ine)=0.0d0

        do ip=1,nump
           xkpp(ip,jp,ine)=0.0d0
	    enddo
	 enddo
     
     do nos=1,nis
        ntem=solid_fem_con(ine,nos) !...connectivity
        x(1:3,nos)=solid_coor_init(1:3,ntem)
        y(1:3,nos)=solid_coor_curr(1:3,ntem)
     enddo
    !...gauss integration
    !...update 06.03.2002, Axel G.: can handle tetrahedral elements as well
     int_x: do ixq=1,nint
	    select case (nis)
	       case (8); rs(1)=xg(ixq,nint)
	       case (4); rs(1)=xg_tetra(ixq,nint)
	    end select
        !rs(1)=xg(ixq,nint)
	    int_y: do iyq=1,nint
		   select case (nis)
	          case (8); rs(2)=xg(iyq,nint)
	          case (4); rs(2)=xg_tetra(iyq,nint)
	       end select
		   !rs(2)=xg(iyq,nint)
		   int_z: do izq=1,nint
		      select case (nis)
	             case (8); rs(3)=xg(izq,nint)
	             case (4); rs(3)=xg_tetra(izq,nint)
	          end select
			  !rs(3)=xg(izq,nint)
!     isoparametric interpolation
              call r_element(rs)
!     y-(r,s)
              call r_jacob(y,xj,xji,det)
!     x-(r,s)
              call r_jacob(x,toxj,toxji,todet)
!     derivative about ox and x
              call r_bdpd_init(toxji)
			  call r_bdpd_curr(xji)
!     deformation gradient
              call r_stoxc(xto,xot,xj,xji,toxj,toxji,toc)
!     material j
              call r_smaterj(wto,toc,xmi,xmj,dxmj,ddxmj)
!     discretized pressure
              call r_spress(rs,ine)
!     continuous pressure
              call r_sbpress(dxmj,ddxmj,xmj)
!     material c
              call r_sboc(obc,ocpp,ocuu,ocup,xmj,dxmj,ddxmj)
!     strain
              call r_sstrain(toc,xto,ixq,iyq,izq,ine)
!     stress
              call r_spiola(ocpp,xmj,dxmj,xto)
              
!     calculate cauchy stress for output
              !if (mod(its,ntsbout).eq.0) then
              call r_scauchy(det,todet,xto,ixq,iyq,izq,ine)
			  !endif

              select case (nis)
	          case (8)
                 wp = wgt(ixq,nint) * wgt(iyq,nint) * wgt(izq,nint)
	             w = wp * todet
				 tot_vol_curr = tot_vol_curr + wp * det
	          case (4)
                 wp = wgt_tetra(ixq,nint) * wgt_tetra(iyq,nint) * wgt_tetra(izq,nint)
	             w = wp * todet / 6  !...integral in tetraheder V_tetra = 1/6 V_cube
                 tot_vol_curr = tot_vol_curr + wp * det / 6
              end select

	          !write(*,'("wp=",f9.6," det=",f9.6," todet=",f9.6)') wp,det,todet
	


	          tot_vol_init = tot_vol_init + w
              call r_sstif(ocpp,ocuu,ocup,ine,w,toxj,ixq,iyq,izq)



		   enddo int_z
		enddo int_y
	 enddo int_x   
  enddo element
  write(*,'("  total solid volume (init) = ",f12.6)') tot_vol_init
  write(*,'("  total solid volume (curr) = ",f12.6)') tot_vol_curr




 !
 !     pressure condensation, inverse kpp
 !
  if (1 == 0) then
  do ine=1,ne_solid
     do ip=1,nump
        tem(ip)=xfp(ip,ine)
        do jp=1,nump
           xfrtem(ip,jp)=xkpp(ip,jp,ine)
 	    enddo
	 enddo
	 !write(*,*) xkpp
	 !stop

     call gaussj(xfrtem,nump,6,tem,1,1)

     do ip=1,nump
        ttm(ip)=0.0d0
        do nos=1,nis
           do isd=1,3
              nu1 = (isd-1)*nn_solid + solid_fem_con(ine,nos)
              mu1 = (isd-1)*nis + nos
              ttm(ip) = ttm(ip) + xkup(mu1,ip,ine)*du(isd,solid_fem_con(ine,nos))
              predrf(nu1) = predrf(nu1) + xkup(mu1,ip,ine)*tem(ip)
			  write(*,*) xkup(mu1,ip,ine)*tem(ip)
		   enddo
	    enddo
     enddo

!...  storage

     do ip=1,nump
        ten(ip)=0.0d0
        do jp=1,nump
           ten(ip)=ten(ip)+xfrtem(ip,jp)*ttm(jp)
	    enddo
        pre(ip,ine)=-tem(ip)-ten(ip)
	 enddo
  enddo

  endif


 !...calculate pressure and stress for the structure output
  if (mod(its,ntsbout).eq.0) then
   
   write(*,*) "  calculate stresses for output"
   do in=1,nn_solid
    
     solid_stress(1:6,in)=0.0d0
     solid_strain(1:6,in)=0.0d0
     
     ntem=0
     pave(in)=0
     if (nsd_solid .ne. 0) then  !do not calculate if it is a point
        do ine=1,ne_solid
           do nos=1,nis
              if (solid_fem_con(ine,nos) .eq. in) then
                 ntem=ntem+1
                 pave(in)=pave(in)+pre(1,ine)
                 do n=1,6
                    solid_stress(n,in)=solid_stress(n,in)+cstr(n,ine,1,1,1)
                    solid_strain(n,in)=solid_strain(n,in)+ge(n,ine,1,1,1)
         	     enddo
                 go to 541
              endif
		   enddo
541	    enddo
        do k=1,6
           solid_stress(k,in)=solid_stress(k,in)/ntem
           solid_strain(k,in)=solid_strain(k,in)/ntem
	    enddo
	    pave(in)=pave(in)/ntem
	 endif
  
   enddo
  
  endif


  call r_common_deallocate

  write(*,*) " done                                 (r_stang)"
  return
 end subroutine r_stang
