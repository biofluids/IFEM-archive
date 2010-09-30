!==================================================

!reset I_fluid_center using I_fluid which is solved before

!===================================================

subroutine reset_Icenter(I_fluid,I_fluid_center,ien,Ic_inter)

  use fluid_variables,only:nn,ne,nen
  use interface_variables

  real(8) I_fluid(nn)
  real(8) I_fluid_center(ne)
  integer ien(nen,ne)
  real(8) Ic_inter

  integer flag(nen),node,inl,ie,flag_sum

  do ie=1,ne
     flag(1:nen)=0
     do inl=1,nen
	node=ien(inl,ie)
	if(I_fluid(node).ge.Ic_inter)then
	   flag(inl)=1
	end if
     end do
     flag_sum=0
     do inl=1,nen
	flag_sum=flag_sum+flag(inl)
     end do
     if (flag_sum.eq.0) then
	I_fluid_center(ie)=0
     else if(flag_sum.eq.nen)then
	I_fluid_center(ie)=1
     else
	I_fluid_center(ie)=0.5
     end if
  end do

  endsubroutine reset_Icenter
