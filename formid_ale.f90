subroutine formid_ale(ids,rngface,ien)
	use fluid_variables
	implicit none
	integer :: ids(nsd,nn),rngface(neface,ne),ien(nen,ne)
	integer :: ieface,inface,inl,iec,irng

ids(:,:)=1

do ieface=1,neface
	do inface=1,nnface
	   inl=mapping(ieface,inface,etype)
		do iec=1,ne
			irng=rngface(ieface,iec)
			if (irng.ne.0) then
				ids(:,ien(inl,iec))=0
			end if
		end do
	end do
end do

end subroutine formid_ale
