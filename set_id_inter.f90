!!!!!!!!!set id for normal!!!!!!!!!!!!!!!!!!!


subroutine set_id_inter(p_var,ids)

  use fluid_variables, only:nn,nsd

  real(8) p_var(nsd,nn)
  integer ids(nn)

  integer i

  do i=1,nn
     if(ids(i) == 1) then
	p_var(1:nsd,i) = 0.0
     end if
  end do

end subroutine set_id_inter
