!remove points that are too close to each other

subroutine points_removal(x_inter,nn_inter)

  use interface_variables,only:maxmatrix
  use fluid_variables, only:nsd

  integer nn_inter_temp
  real(8) temp
  real(8) x_inter_temp(nsd,maxmatrix)
  real(8) x_inter(nsd,maxmatrix)
  integer icount,j
  nn_inter_temp=0
  do icount=1,nn_inter
!	if( (x_inter(1,icount).gt.(0.5-0.05*max_hg)) .or. &
!	    (x_inter(1,icount).lt.(-0.5+0.05*max_hg)).or. &
!	    (x_inter(2,icount).gt.(0.5-0.05*max_hg)) .or. &
!	    (x_inter(2,icount).lt.(-0.5+0.05*max_hg)).or. &
!	    (x_inter(3,icount).gt.(1.0-0.05*max_hg)) .or. &
!	    (x_inter(3,icount).lt.(0.0+0.05*max_hg)) ) then
!         if((x_inter(2,icount).lt.0.0)) then
	 temp=(x_inter(1,icount)+0.5)**2+(x_inter(2,icount)-1.0)**2
	 temp=sqrt(temp)
	 if(temp.lt.0.5) then
            goto 789
	end if
     nn_inter_temp=nn_inter_temp+1
     x_inter_temp(1:nsd,nn_inter_temp)=x_inter(1:nsd,icount)
789 continue
  end do
  nn_inter=nn_inter_temp
  x_inter(1:nsd,1:nn_inter)=x_inter_temp(1:nsd,1:nn_inter)

end subroutine points_removal
