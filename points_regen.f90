!=======================================================!
!        Regenerate interface points                    !
!=======================================================!

subroutine points_regen(I_fluid,inter_ele_nn,x_center,I_fluid_center,inter_ele,ne_inter,Ic_inter,corr_Ip,x,x_inter,ien,nn_inter_regen,x_inter_regen,hg,infdomain_inter)

  use fluid_variables, only:nsd,nn,ne,nen
  use interface_variables
  
  real(8) I_fluid(nn)
  integer inter_ele_nn(ne)
  real(8) I_fluid_center(ne)      !initial indicator for element center
  real(8) x_center(nsd,ne)            !coordinates for element center
  integer inter_ele(ne)           !interfacial elements index
  integer ne_inter                !# of interfacial elements
  real(8) Ic_inter                !constant interfacial indicator
  real(8) corr_Ip(maxmatrix)
  real(8) x(nsd,nn)                !fluid coordinates
!  integer nn_inter                 !original # of interface points
  real(8) x_inter(nsd,maxmatrix)   !original coordinats of interface points
  integer ien(nen,ne)
  integer nn_inter_regen           !new # of interface points
  real(8) x_inter_regen(nsd,maxmatrix) !new coordinates of interface points
  real(8) hg(ne)
  integer infdomain_inter(maxmatrix)

  integer i,j,isd,inl,node,ie,icount,jcount
  real(8) x_fluid(nsd,nen)         !global coordinates of the nodes of each element
  integer nn_sub                   !# of candidate points per fluid element
  real(8) x_loc_can(nsd,5000)       !local coordinats of candidate points per element           
  real(8) x_glo_can(nsd,5000)       !global coordinats of candidat points per element

  real(8) sh(nen)
  real(8) temp(nsd)
  real(8) xtemp(nsd)
  real(8) xlocan(nsd)
  real(8) distance
  real(8) hs
  real(8) dx(nsd), Sp, Ip_temp1,Ip_temp2,xIp_temp1
  integer nn_cr                    !minimum # of points each element should contain
  integer nn_local                 !# of points of each element
  real(8) x_inter_regen_local(nsd,50) !#coordinates of local regen points
  real(8) err, signf
  integer flag_loc(nen),flag_sum
  nn_inter_regen=0
  x_inter_regen(:,:)=0.0
  nn_cr=6
!  write(*,*)'infdomain_inter',infdomain_inter(1:nn_inter)
!write(*,*)'inter_ele_nn=',inter_ele_nn(1:ne_inter)
  do ie=1,ne_inter  ! loop over interfacial elements
!  do ie=62,62
     flag_loc(1:nen)=0
     do inl=1,nen
	node=ien(inl,inter_ele(ie))
	x_fluid(1:nsd,inl)=x(1:nsd,node)
	if(I_fluid(node).ge.Ic_inter) then
	  flag_loc(inl)=1
	end if
     end do   !assign global coordinates to each element node

     flag_sum=0
     do inl=1,nen
	flag_sum=flag_sum+flag_loc(inl)
     end do

     if((flag_sum.eq.0) .or. (flag_sum.eq.nen)) then
!	write(*,*)ie,' not an interfacial element'
!	goto 999
     end if


     if (nsd==2 .and. nen==4) then     !2d4nodes
!      nn_sub=3                        ! 4*4 candidate points
      nn_local=0
      icount=2
      do while((nn_local.eq.0).and.(icount.le.2))
         if(icount==1) then
	    signf=1.0
	 else if(icount==2) then
            signf=-1.0
         end if
      nn_sub=4
!      x_inter_regen_local(:,:)=0.0
      do while((nn_local.le.nn_cr).and.(nn_sub.le.30)) ! begin regeneration loop
	nn_local=0                ! reset nn_local
	x_inter_regen_local(:,:)=0.0 ! reset local coordinates
	do i=1,nn_sub
	   do j=1,nn_sub
		x_loc_can(1,nn_sub*(i-1)+j)=2.0/nn_sub*i-1.0-1.0/nn_sub
		x_loc_can(2,nn_sub*(i-1)+j)=2.0/nn_sub*j-1.0-1.0/nn_sub
	   end do
	end do  ! assign local coordinates for candidate points
	x_glo_can(:,:) = 0.0
	do i=1,nn_sub**2  ! loop over the candidate points
	   sh(1)=0.25*(1-x_loc_can(1,i))*(1-x_loc_can(2,i))
	   sh(2)=0.25*(1+x_loc_can(1,i))*(1-x_loc_can(2,i))
	   sh(3)=0.25*(1+x_loc_can(1,i))*(1+x_loc_can(2,i))
	   sh(4)=0.25*(1-x_loc_can(1,i))*(1+x_loc_can(2,i))
!!!!!!!!!!!calculate the shape function for candidate points!!!!!!!!!!!!!!!
	   temp(:) = 0.0
	   do inl=1,nen
		temp(1:nsd)=temp(1:nsd)+sh(inl)*x_fluid(1:nsd,inl)
	   end do
	   xlocan(:)=temp(:)
!!!!!!!!!!calculate the global coordinates for the candidate poinits!!!!!!!!
	   Ip_temp1 = 0.0     
!	   hs=hg(inter_ele(ie))
	   do j=1,ne
	        hs=hg(j)
!		if (I_fluid_center(j).gt.1.0e-6) then
		   dx(:)=abs(temp(:)-x_center(:,j))
		   call B_Spline(dx,hs,nsd,Sp)
		   Ip_temp1=Ip_temp1+I_fluid_center(j)*Sp
!		end if
	   end do  ! get initial Ip for candidate points
	   Ip_temp2 = 0.0
	   do j=1,nn_inter
		hs=hg(infdomain_inter(j))/1.0
		dx(:)=abs(temp(:)-x_inter(:,j))
		call B_Spline(dx,hs,nsd,Sp)
		Ip_temp2=Ip_temp2+corr_Ip(j)*Sp
	   end do ! get corrected Ip for candidate points

	   Ip_temp1=Ip_temp1+Ip_temp2
!           write(*,*)'Ip_ini_can=',Ip_temp1
!	   if( (signf*(Ip_temp1-Ic_inter) .le. 0.5*Ic_inter).and.(signf*(Ip_temp1-Ic_inter).ge.0.0)) then
!           if(((Ip_temp1-Ic_inter).le.(0.25*Ic_inter)).and.((Ic_inter-Ip_temp1).le.(0.25*Ic_inter))) then
         if(((Ip_temp1-Ic_inter).le.0.15).and.((Ic_inter-Ip_temp1).le.0.15)) then

		call point_projection(Ic_inter,temp,x_center,I_fluid_center,x_inter,corr_Ip,hs,Ip_temp1,err,hg,infdomain_inter)
!		write(*,*)'err_regen=',err
		distance=sqrt((xlocan(1)-temp(1))**2+(xlocan(2)-temp(2))**2)
		if(err .lt. 1.0e-8) then
		  if(distance.le.(hs/nn_sub/2)) then
		  nn_local=nn_local+1
		  x_inter_regen_local(1:nsd,nn_local)=temp(1:nsd)
		  else
!		   write(*,*)'distance is too big'
		  end if
		else
!		  write(*,*)'err is too big. newton fails. err=',err
!		  write(*,*)'ne_inter=',ie
		end if
!		nn_inter_regen=nn_inter_regen+1
!		x_inter_regen(1:nsd,nn_inter_regen)=temp(1:nsd)
	   end if
	end do ! end of loop over the candidate points
	
	nn_sub=nn_sub+4
      end do ! end of do while

!write(*,*)'ie=',ie,'nn_local=',nn_local,'icount=',icount
       icount=icount+1
      end do ! end of do while 2
     end if    ! end of 2d4nodes


      nn_inter_regen=nn_inter_regen+nn_local
      x_inter_regen(1:nsd,(nn_inter_regen-nn_local+1):nn_inter_regen)= &
                                x_inter_regen_local(1:nsd,1:nn_local)

999 continue
   end do ! end of loop over the interfacial elements

write(*,*)'nn_inter_regen=',nn_inter_regen
write(*,*)'Ic_inter=',Ic_inter

end subroutine points_regen




















