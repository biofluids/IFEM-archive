!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine do interpolation for velocity (ibuf==1) or make distribution for FSI force (ibuf==2)
! input:
! data_solids: source (velocity) for 1 & result (FSI force) for 2
! data_fluids: result (velocity) for 2 & source (FSI force) for 1
! nn_fluids: number of nodes in fluid domain
! nn_solids: number of nodes in solid domain
! nsd: spatial degree freedom
! ibuf: control parameter
! ne: number of elements in fluid domain
! nen: number of nodes per element in fluid domain
! ne_solid: number of elements in solid domain
! nen_solid: number of nodes per element in solid domain
! xyz_solid: xyz coordinates for solid points
! ien_solid_m(ne_solid,nen_solid) !conectivity matrix for solid from hypo
! ien_solid(nen_solid,ne_solid) ! connectivity matrix for FEM interpolation
! xyz_fluid(nsd, nn_fluids) !xyz coordinates for fluid points
! ien(nen,ne) !conectivity matrix for fluid
! infdomain: influence domain
! output:
! data_solids for ibuf==1
! data_fluids for ibuf==2
! Xingshi Wang 2008
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine data_exchange_FEM(data_solids,nn_solids,data_fluids,nn_fluids,dv,nsd,ibuf,ne,nen,&
                              ne_solid,nen_solid,xyz_solid,ien_solid_m,xyz_fluid,ien,infdomain)
  implicit none
  integer,intent(in) :: ibuf
  integer nsd
  integer ne
  integer nen
  integer nen_solid
  integer ne_solid

 !...solids variables
  integer,intent(in)   :: nn_solids
  real(8),intent(inout) :: data_solids(nsd,nn_solids)
  real(8) xyz_solid(nsd,nn_solids) ! xyz coordinates for solid points
  integer ien_solid_m(ne_solid,nen_solid) !conectivity matrix for solid from hypo
  integer ien_solid(nen_solid,ne_solid) ! connectivity matrix for FEM interpolation

 !...fluids variables
  integer,intent(in)   :: nn_fluids
  real(8),intent(inout) :: data_fluids(nsd,nn_fluids)
  real(8),intent(in)    :: dv(nn_fluids)
  real(8) xyz_fluid(nsd, nn_fluids) !xyz coordinates for fluid points
  integer ien(nen,ne) !conectivity matrix for fluid
  
!===========================================================
  integer infdomain(nn_solids) ! influence domain matrix
  
!********************************
! distribution loop variables

  integer index_el(nen)
  real(8) xyz_el(nsd,nen)
  !distribution coefficients for force
  real(8) dis_solid(nen)
  
!***********************************
  
 ! local varialbes
  integer inn ! loop variable for global nodals
  real(8) x(nsd) ! searching point for getinf_el_3d
  integer finf ! index of the searched element for getinf_el_3d
  integer maxconn ! constant for getinf_el_3d, the maximun number we can bear for the step of search
  integer inf_index(nen) ! the global indexes of the nodals formed the element we got



  integer ii ! loop variable for the nodals in the element
  real(8) inf_xyz(nsd,nen) ! xyz coordinates for the influence nodals



  real(8) sh(nen) ! shape function for each element in fluid domain
  
  ! check the sum of force both in fluid and solid
  real(8) force_solid(2)
  real(8) force_fluid(2)

  
  
  integer isd  ! loop variable for nsd
  integer ie  ! loop varible for ne
  do ii=1,nen_solid
     do ie=1,ne_solid
        ien_solid(ii,ie)=ien_solid_m(ie,ii)
     end do
  end do

  if (ibuf ==1) then  !velocity interpolation
        write(*,*) '***Velocity interpolation using FEM method***'

        data_solids(:,:)=0

        open (unit=40, file='interface.dat', status='unknown')
        open (unit=41, file='exchange_solid.dat', status='unknown')
        open (unit=42, file='sh.dat', status='unknown')
           do inn=1,nn_solids

           x(1:nsd)=xyz_solid(1:nsd,inn)

           do ii=1,nen
              inf_index(ii)=ien(ii,infdomain(inn))
              inf_xyz(1:nsd,ii)=xyz_fluid(1:nsd,inf_index(ii))
           end do
!           write(40,*) inf_xyz(1,:), inf_xyz(2,:)
           call sh_exchange(x,inf_xyz,nsd,nen,sh)

          do ii=1,nen
             data_solids(1:nsd,inn) = data_solids(1:nsd,inn) +&
             data_fluids(1:nsd,inf_index(ii)) * sh(ii)


          end do
!         write(42,*) sh(:)
!         write(41,*) data_solids(1,inn), data_solids(2,inn)
        end do





     else if (ibuf ==2) then !force distribution
        open(unit=43, file='dis_coe.dat', status='unknown')
        open(unit=44, file='force_fluid.dat', status='unknown')
        open(unit=45, file='force_solid.dat', status='unknown')
         write(*,*) '*** Distributing Forces onto Fluid using FEM interpolation***'
        ! pre-step pick the element near the solid domain
        data_fluids(:,:)=0
! check if we lost any information when we do the distribution when needed
!        force_solid(:)=0
!        force_fluid(:)=0
        

        
!        write(*,*) 'I am in the inverse change'
       
        do inn=1, nn_solids
           
        
           do ii=1,nen
              index_el(ii)=ien(ii,infdomain(inn))
              do isd=1,nsd
               xyz_el(isd,ii)=xyz_fluid(isd,index_el(ii))
              end do
           end do
           x(:)=xyz_solid(:,inn)

           dis_solid(:)=0
           call sh_exchange(x,xyz_el,nsd,nen,dis_solid)
           do ii =1,nen
           data_fluids(:,index_el(ii))=data_fluids(:,index_el(ii))+data_solids(:,inn)*dis_solid(ii)
           end do
        end do
        do inn=1,nn_solids
     !      do isd=1,nsd
        !      if (data_solids(isd,inn)>0) then
      !           force_solid(isd)=force_solid(isd)+data_solids(isd,inn)
        !      end if
     !      end do
        end do 
      !    write(*,*) 'sum of fsi for solid', force_solid(1) , force_solid(2)
        do inn=1,nn_fluids
 !          write(44,*) data_fluids(1,inn), data_fluids(2,inn)
       !    do isd=1,nsd
       !    if (data_fluids(isd,inn)>0) then
       !       force_fluid(isd)=force_fluid(isd)+data_fluids(isd,inn)
       !    end if
       !   end do
        end do
        ! write(*,*) 'sum of fsi for fluid', force_fluid(1) , force_fluid(2)
        ! write(*,*) 'difference of force= ', force_solid(1)-force_fluid(1)






     end if
return
end
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the shape function for both 2-D and 3-D case
! 3-D case has not been finished yet
! The shape function here also can be used as the cooeffiecints for the force distribution
! Except triangle, other shape function need be fixed manually,
! since the digital error can not be avoided
! Xingshi Wang 2008
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine sh_exchange(x,xe,nsd,nen,sh)
integer nen
integer nsd
real(8) x(nsd)
real(8) xe(nsd,nen) ! xyz coordinates of the element nodals
real(8) sh(nen)
real(16) sh1(nen)
real(16) me(nen,nen)
real(8) me1(nen,nen)
real(16) invme(nen,nen)
real(8) det
real(16) xp(nen)
integer indx(nen)
integer i ! loop variable
integer j
        do i=1,nen
        sh1(i)=0
        do j=1,nen
        me(i,j)=0
        end do
        end do
        
        if (nsd .eq. 2) then  ! 2-D case
           if(nen .eq. 3) then ! 2-D triangle
           do i=1,nen
           me1(i,1)=1.0
           me1(i,2:(nsd+1))=xe(1:nsd,i)
           end do
           call determinant(me1,nen,nen,det)
          ! write(*,*) det
           sh(1)=(xe(1,2)*xe(2,3)-xe(1,3)*xe(2,2)+&
           (xe(2,2)-xe(2,3))*x(1)+&
           (xe(1,3)-xe(1,2))*x(2))/(det)
           
           sh(2)=(xe(1,3)*xe(2,1)-xe(1,1)*xe(2,3)+&
           (xe(2,3)-xe(2,1))*x(1)+&
           (xe(1,1)-xe(1,3))*x(2))/(det)

           sh(3)=(xe(1,1)*xe(2,2)-xe(1,2)*xe(2,1)+&
           (xe(2,1)-xe(2,2))*x(1)+&
           (xe(1,2)-xe(1,1))*x(2))/(det)
          
           sh(3)=1-sh(1)-sh(2)
          
           else if (nen .eq. 4) then ! 2-D quadrilateral


           do i=1,nen
           me(i,1)=1.0
           me(i,2)=xe(1,i)
           me(i,3)=xe(2,i)
           me(i,4)=xe(1,i)*xe(2,i)
           end do
           xp(1)=1
           xp(2)=x(1)
           xp(3)=x(2)
           xp(4)=x(1)*x(2)

           call MIGS(me,nen,invme,indx) ! get inverse

           do j=1,nen
           do i=1,nen
 !          write(*,*) 'invme', invme(i,j)
           sh1(j)=sh1(j)+xp(i)*invme(i,j)
           sh(j)=sh1(j)
           end do
           !write(*,*) 'shape1', sh1(j)
           end do
           end if
        else
           ! 3-D case




        end if
        do i=1,nen
          if (sh(i)<0) then
            call fix_sh(sh,nen)
            
          

        write(*,*) '****errors in shape function, running fix munually***'
          do j=1,nen
!         write(*,*) 'xe=', xe(1,j), xe(2,j)
         write(*,*)'sh= ', sh(j)
!         write(*,*) 'xp= ', xp(j)
         write(*,*) 'sum of shape function should be 1', sum(sh)
          end do
          go to 1000
          end if
       
        end do
        
1000    continue
return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Manually fix the digital error in shape function
! Right now assume we only have one negetive shapefunction
! Take the two smallest number (one is the negetive value)
! sum(these two terms)=1-sum(rest terms)
! further assume that these two equal to each other
! then fix them
! input:
! nen: number of nodes per element
! sh: shape function of this element
! output:
! sh: the fixed shape funtion
! (all larger tham zero and less than one; sum of them equals one)
! Xingshi Wang 2008
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fix_sh(sh,nen)
real(8) sh(nen)
integer nen

integer i
real(8) tempsh(nen)
real(8) fixsh1
real(8) fixsh2
real(8) tempsum
integer index1
integer index2

   fixsh1=1.0
   fixsh2=1.0
   tempsum=0

   tempsh(:)=sh(:)
   do i=1,nen
      if (tempsh(i)< fixsh1) then
      fixsh1=tempsh(i)
      index1=i
      end if
   end do
   do i=1,nen
      if (tempsh(i)< fixsh2 .and. tempsh(i)> fixsh1) then
      fixsh2=tempsh(i)
      index2=i
      end if
   end do
   
   do i=1,nen
   if ((i .ne. index1) .and. (i .ne. index2)) then
   tempsum=tempsum+tempsh(i)
   end if
   end do
   tempsum=1-tempsum
   fixsh1=tempsum/2.0
   fixsh2=tempsum/2.0
   sh(index1)=fixsh1
   sh(index2)=fixsh2

   return
   
end
