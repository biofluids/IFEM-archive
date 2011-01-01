subroutine formid_ale(ids,rngface,ien)
	use fluid_variables
	implicit none
	integer :: ids(nsd,nn),rngface(neface,ne),ien(nen,ne)
	integer :: ieface,inface,inl,iec,irng
        integer :: check

check=0
ids(:,:)=1

do ieface=1,neface
	do inface=1,nnface
	   inl=mapping(ieface,inface,etype)
		do iec=1,ne
			irng=rngface(ieface,iec)
			if (irng.ne.0) then
				check=check+1
			end if
		end do
		do iec=1,ne
			irng=rngface(ieface,iec)
			if (irng.ne.0) then
				ids(:,ien(inl,iec))=0
			!----------------------------
			! Adding for left and right wall only
				if ((irng .eq. 1) .or. (irng .eq. 3)) then
					if (check .le. 1) then
						ids(:,ien(inl,iec)) = 1
						ids(3,ien(inl,iec)) = 0
					end if
				end if
			!----------------------------
			end if
		end do
	end do
end do

end subroutine formid_ale
