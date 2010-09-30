!!!!!!!!!set id for normal!!!!!!!!!!!!!!!!!!!


subroutine set_id_curv(p_var,ids)

  use fluid_variables, only:nn

  real(8) p_var(nn)
  integer ids(nn)

  integer i

  do i=1,nn
     if(ids(i) == 1) then
	p_var(i) = 0.0
     end if
  end do

end subroutine set_id_curv
