!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module delta nonuniform 
!
! Lucy Zhang, Axel Gerstenberger
! NWU, 04/22/2003
!
! contains:
!   - variables          !...all variables related to the delta function
!   - delta_initialize   !...calculate domain of influence for each solid node
!   - delta_exchange     !...performs exchange of information in both directions (fluid <--> solid)

module delta_nonuniform
  implicit none
  save

public
  !...use these parameters to define the direction of information flow
  integer,parameter :: delta_exchange_fluid_to_solid = 1
  integer,parameter :: delta_exchange_solid_to_fluid = 2  

  integer :: maxconn !...used to define connectivity matrix size

  integer :: ndelta !...defines type of delta function used -> right now, only "1" (RKPM) is available

  real(8),allocatable :: shrknode(:,:)      !...shape function for each node, contains the weights
  integer,allocatable :: cnn(:,:),ncnn(:)  !...connectivity arrays for domain of incluence for each solid node

 !...private subroutines
  private :: getinf  !,correct3d

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine delta_initialize(nn_solids,x_solids,xna,ien,dwjp)
! Subroutine rkpm_delta
! Lucy Zhang
! 11/06/02
! Northwestern University

! This subroutine calculate the delta function using RKPM which is used for both
! ninterpolation and distribution of the velocities and forces respectively
! between fluids and solids domain.
  use fluid_variables
  implicit none

 !...solids variables
  integer,intent(in)  :: nn_solids
  real(8),intent(in)  :: x_solids(nsd,nn_solids)

 !...fluids variables
  real(8) :: xna(nsd,nn),xn(nsd,nen)
  integer,intent(in)  :: ien(nen,ne)
  real(8),intent(out) :: dwjp(nn)
  real(8) :: adist(nsd,nn)

 !...local variables
  real(8) :: x(nsd), y(nsd), a(nsd)
  real(8) :: xr(nsd,nsd), cf(nsd,nsd) 
  real(8) :: b(4), bd(nsd,4)
  real(8) :: shp, det
  real* 8 :: xmax(nsd),vol
  real(8) :: coef,avginf
  integer :: iq
  integer :: maxinf,mininf,totinf
  integer :: ie,inl,isd,nnum,node
  integer :: inf(maxconn),ninf
  integer :: i,n,error_id


  
  write(*,*) "*** Calculating RKPM delta function ***"

  if (allocated(shrknode)) then
     deallocate(shrknode)
  end if
  if (allocated(cnn)) then
     deallocate(cnn)
  end if
  if (allocated(ncnn)) then
     deallocate(ncnn)
  end if


  allocate(shrknode(maxconn,nn_solids),stat=error_id)
  allocate(cnn(maxconn,nn_solids)     ,stat=error_id)
  allocate(ncnn(nn_solids)            ,stat=error_id)
  coef = 0.7d0
  write(*,*) '***RKPM coefficient equals***', coef
  maxinf = 0
  mininf = 9999
  avginf = 0
  cnn(:,:)=0
  ncnn(:)=0
  shrknode(:,:)=0.0d0
  
 !...Calculate nodal weights
  dwjp(:) = 0.0
  adist(:,:) = 0.0
  totinf = 0
  do ie = 1,ne
     do inl=1,nen
        do isd=1,nsd
           nnum = ien(inl,ie)
           xn(isd,inl) = xna(isd,nnum)
        enddo
     enddo

	do isd=1,nsd
		xmax(isd) = coef*(maxval(xn(isd,1:nen)) - minval(xn(isd,1:nen)))
	enddo

     do inl = 1,nen
        node = ien(inl,ie) 
        adist(1:nsd,node) = max(adist(1:nsd,node),xmax(1:nsd))
     enddo

 !...Calculate volume
	if (nsd == 2) then
	  if (nen == 3) then
			include "vol2d3n.fi"
	  elseif (nen == 4) then
			include "vol2d4n.fi"
	  endif
	elseif (nsd == 3) then
      if (nen == 4) then
			include "vol3d4n.fi"
      elseif (nen == 8) then 
			include "vol3d8n.fi"
	  endif
    endif

     do inl = 1,nen
        nnum = ien(inl,ie)
        dwjp(nnum) = dwjp(nnum) + vol/nen
     enddo
  enddo

! Calculate the RKPM shape function for the solids points
  do i = 1, nn_solids
     x(1:nsd)=x_solids(1:nsd,i) !get solids point coordinate
     ninf=0
     inf(:)=0
! get a list of influence nodes from the fluids grid
     call getinf(inf,ninf,x,xna,adist,nn,nsd,maxconn)
    cnn(1:ninf,i)=inf(1:ninf)
     ncnn(i)=ninf

     if (ninf > maxinf) then
			maxinf = ninf
     elseif (ninf < mininf) then 
			mininf = ninf
	 endif
     totinf = totinf + ninf
! calculate the correction function
     if (nsd==2) then 
		call correct2d(b,bd,x,xna,adist,dwjp,nn,1,inf,ninf,maxconn)
     elseif (nsd==3) then 
		call correct3d(b,bd,x,xna,adist,dwjp,nn,1,inf,ninf,maxconn)
	 endif

     do n = 1, ninf
        nnum = inf(n)
        do isd = 1,nsd
           y(isd) = xna(isd,nnum)
           a(isd) = adist(isd,nnum)
        enddo

		 if (nsd==2) then 
	        call RKPMshape2d(shp,b,bd,x,y,a,dwjp(nnum))
		 elseif (nsd==3) then 
			call RKPMshape3d(shp,b,bd,x,y,a,dwjp(nnum))
		 endif
        shrknode(n,i)=shp
     enddo

  enddo
  avginf = totinf/nn_solids
  write(6,'("  Maximum Influence Nodes = ",i7)') maxinf
  write(6,'("  Minimum Influence Nodes = ",i7)') mininf
  write(6,'("  Average Influence Nodes = ",f7.2)') avginf
!  write(6,'("Come up man! =",i7 )')  avaginf
  return
end subroutine delta_initialize

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c This subroutine finds the influence points of point x
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getinf(inf,ninf,x,xna,adist,nn,nsd,maxconn)
  implicit none

  integer :: ninf,nn,nsd,maxconn
  real(8) x(nsd), xna(nsd,nn), adist(nsd,nn)
  real(8) r(nsd)
  integer inf(maxconn)
  integer i

!cccccccccccccccccc
!   x = the coordinate of the point to be calculated for
!   xna = the coordinate of all points
!   r = the distance between point x and other points in the system
!   inf = a collection of all the influence points
!   ninf = total number of influence points
!   adist = the radial distance of the influence domain
!cccccccccccccccccc
 !  open(unit=400, file='interface_RKPM.dat',status='unknown')
  ninf = 0
  do i = 1,nn
     r(1:nsd) = x(1:nsd) - xna(1:nsd,i)
	 if (nsd==3) then
		if ((abs(r(1))<=2*adist(1,i)).and.(abs(r(2))<=2*adist(2,i)).and.(abs(r(3))<=2*adist(3,i))) then
			ninf = ninf + 1
			inf(ninf) = i
		endif
	 elseif (nsd==2) then
		if ((abs(r(1))<=2*adist(1,i)).and.(abs(r(2))<=2*adist(2,i))) then
			ninf = ninf + 1
			inf(ninf) = i
                     !   write(400,*) i  
                     !   write(*,*) '****I am in getinf ****'      
		endif
	 endif
   enddo
  !write(*,*) 'ninf=  ', inf(1)
  if (ninf > maxconn) then
     write (*,*) "Too many influence nodes!"
     write (*,*) ninf
  elseif (ninf.lt.4) then
     write (*,*) "Not enough influence nodes!"
     write (*,*) ninf
  endif

  return
end subroutine getinf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!SUBROUTINE DELTA
! Lucy Zhang
! 11/06/02
!
! Northwestern University
! This subroutine calculate the delta function which is used for both
! interpolation and distribution of the velocities and forces respectively
! between fluids and solids domain.
!
! There are 3 options for calculating
! 1. RKPM - cubic spline for non-uniform spacing
! 2. RKPM - cubic spline for uniform spacing
! 3. Original delta function for uniform spacing 

subroutine delta_exchange(data_solids,nn_solids,data_fluids,nn_fluids,ndelta,dv,nsd,ibuf)
  implicit none
  integer,intent(in) :: ibuf,ndelta
  integer nsd
 !...solids variables
  integer,intent(in)   :: nn_solids
  real(8),intent(inout) :: data_solids(nsd,nn_solids)

 !...fluids variables
  integer,intent(in)   :: nn_fluids
  real(8),intent(inout) :: data_fluids(nsd,nn_fluids)
  real(8),intent(in)    :: dv(nn_fluids)

 !...local variables
  integer :: inn,icnn,pt
  real(8)  :: tot_vel(nn_fluids),tot_vel_fluid,vol_inf 
  real(8) force_solid
  real(8) force_fluid
  tot_vel_fluid = 0
  tot_vel(1:nn_fluids)=0.0
  force_solid=0
  force_fluid=0
  if (ndelta == 1) then                  !c    If non-uniform grid 

     if (ibuf == delta_exchange_fluid_to_solid) then  !velocity interpolation
        data_solids(:,:)=0
        do inn=1,nn_solids
           vol_inf=0.0d0
           do icnn=1,ncnn(inn)
              pt=cnn(icnn,inn)
              data_solids(1:nsd,inn) = data_solids(1:nsd,inn) + data_fluids(1:nsd,pt) * shrknode(icnn,inn)
              tot_vel(pt)=data_fluids(1,pt)
              vol_inf = vol_inf + dv(pt)
           enddo
        enddo

     elseif (ibuf == delta_exchange_solid_to_fluid) then !force distribution
        write(*,*) '*** Distributing Forces onto Fluid ***'
        data_fluids(:,:)=0
        do inn=1,nn_solids
           do icnn=1,ncnn(inn)
              pt=cnn(icnn,inn)
              data_fluids(1:nsd,pt) = data_fluids(1:nsd,pt) + data_solids(1:nsd,inn) * shrknode(icnn,inn)  
		   enddo    
        enddo
       write(*,*) 'sum of solid', sum(data_solids(1,:))
       write(*,*) 'sum of fluid', sum(data_fluids(1,:))
       do inn=1,nn_solids
          if (data_solids(1,inn)>0) then
             force_solid=force_solid+data_solids(1,inn)
          end if
       end do
       do inn=1,nn_fluids
          if (data_fluids(1,inn)>0) then
             force_fluid=force_fluid+data_fluids(1,inn)
          end if
       end do
       write(*,*) 'positve force of solid', force_solid
       write(*,*) 'positve force of fluid', force_fluid



     endif
  else


  endif ! end of option for delta function


  return
end subroutine delta_exchange

end module delta_nonuniform
