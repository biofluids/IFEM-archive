subroutine r_stang
  use solid_variables
  use r_common
  implicit none

  real*8 x(3,9),toc(3,3),xto(3,3),xot(3,3),xj(3,3),xji(3,3),rs(3),toxj(3,3),toxji(3,3)
  real*8 xmj(3),xmi(3),dxmj(3,6),ddxmj(3,6,6),obc(6,6),ocuu(6,6),ocup(6)
  real*8 xfrtem(6,6),tem(6),ten(6),ttm(6)

  integer :: ine,i,j,m,k,nu1,mu1

  real*8 :: w,wp,tot_vol,body_force
  integer :: nos,ntem
  integer :: ixq,iyq,izq,isd
  real*8 :: det,todet,wto,ocpp


  write(*,*) "start r_stang"

  do i=1,3*nn_solid
     predrf2(i) = predrf(i)
     predrf(i) = 0.0d0
     drf2(i) = drf(i)
  enddo
    
  tot_vol=0.0
  body_force=0.0
  element: do ine=1,ne_solid
 !...position
     do j=1,nump
        xfp(j,ine)=0.0d0
        do i=1,nis
           xkup(i,j,ine)=0.0d0
           xkup(i+nis,j,ine)=0.0d0
	       xkup(i+2*nis,j,ine)=0.0d0
	    enddo
        do i=1,nump
           xkpp(i,j,ine)=0.0d0
	    enddo
	 enddo
     
     do nos=1,nis
        ntem=nea(ine,nos) !...connectivity
        do isd=1,3
           x(isd,nos)=solid_coor_init(isd,ntem)
		   y(isd,nos)=solid_coor_curr(isd,ntem)
           !y(isd,nos)=coor(ntem,isd)+dis(isd,ntem)
        enddo
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
	          !write(*,*) rs
	          !write(*,*) "nis = ",nis
!     isoparametric interpolation
              call r_element(rs)
!     y-(r,s)
              call r_jacob(y,xj,xji,det)
!     x-(r,s)
              call r_jacob(x,toxj,toxji,todet)
	          !write(*,*) det,todet
!     derivative about ox
              call r_bdpd(toxji)
!     deformation gradient
              call r_stoxc(xto,xot,xj,xji,toxj,toxji,toc)
	          !write(*,*) xto,xot
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
              call r_spiola(ocpp,xmj,dxmj)
!     discretized pressure
              select case (nis)
	          case (8)
                 wp = wgt(ixq,nint) * wgt(iyq,nint) * wgt(izq,nint)
	             w = wp * todet
	          case (4)
                 wp = wgt_tetra(ixq,nint) * wgt_tetra(iyq,nint) * wgt_tetra(izq,nint)
	             w = wp * todet / 6  !...integral in tetraheder V_tetra = 1/6 V_cube
              end select

	          !write(*,'("wp=",f9.6," det=",f9.6," todet=",f9.6)') wp,det,todet
	
	          tot_vol=tot_vol+w
              call r_sstif(ocpp,ocuu,ocup,ine,w,toxj,body_force)
 
              call r_scauchy(det,todet,xto,ixq,iyq,izq,ine)

		   enddo int_z
		enddo int_y
	 enddo int_x   
  enddo element
  write(*,'(" total solid volume = ",f12.6)') tot_vol
  !write(*,*) 'body force=',body_force

 !
 !     pressure condensation, inverse kpp
 !
  do ine=1,ne_solid
     do i=1,nump
        tem(i)=xfp(i,ine)
        do j=1,nump
           xfrtem(i,j)=xkpp(i,j,ine)
 	    enddo
	 enddo

     call gaussj(xfrtem,nump,6,tem,1,1)

     do i=1,nump
        ttm(i)=0.0d0
        do k=1,nis
           do m=1,3
              nu1 = (m-1)*nn_solid + nea(ine,k)
              mu1 = (m-1)*nis + k
              ttm(i) = ttm(i) + xkup(mu1,i,ine)*du(m,nea(ine,k))
              predrf(nu1) = predrf(nu1) + xkup(mu1,i,ine)*tem(i)
		   enddo
	    enddo
     enddo
!
!     storage
!
     do i=1,nump
        ten(i)=0.0d0
        do j=1,nump
           ten(i)=ten(i)+xfrtem(i,j)*ttm(j)
	    enddo
        pre(i,ine)=-tem(i)-ten(i)
	 enddo
  enddo

  write(*,*) "end r_stang"
  return
 end subroutine r_stang
