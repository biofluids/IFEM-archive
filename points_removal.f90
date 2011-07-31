!remove points that are too close to each other

subroutine points_removal(x_inter)

  use interface_variables
  use fluid_variables, only:nsd

  integer nn_inter_temp
  real(8) temp
  real(8) x_inter_temp(nsd,maxmatrix)
  real(8) x_inter(nsd,maxmatrix)
  integer icount,j
  nn_inter_temp=0
  do icount=1,nn_inter
     do j=1,nn_inter_temp
        temp=(x_inter_temp(1,j)-x_inter(1,icount))**2+(x_inter_temp(2,j)-x_inter(2,icount))**2
        temp=sqrt(temp)
        if(temp.lt.max_hg/10.0) then
          goto 789
        end if
     end do
        if( (x_inter(1,icount).gt.(0.5)) .or. &
            (x_inter(1,icount).lt.(-0.5)).or. &
            (x_inter(2,icount).gt.(0.5)) .or. &
            (x_inter(2,icount).lt.(-0.5))) then
          goto 789
        end if

     nn_inter_temp=nn_inter_temp+1
     x_inter_temp(1:nsd,nn_inter_temp)=x_inter(1:nsd,icount)
789 continue
  end do
  nn_inter=nn_inter_temp
  x_inter(1:nsd,1:nn_inter)=x_inter_temp(1:nsd,1:nn_inter)

end subroutine points_removal
