subroutine error(string, param, fatal)
     	implicit none
      	character string*(*)
      	integer :: param, strlen
      	logical :: fatal

       	if (string.ne.' ') then
        	strlen = len(string)
            	do while (string(strlen:strlen).eq.' ')
               		strlen = strlen - 1
            	end do

            	if (param.ne.-999) then
               		write(unit=0, fmt='(a$)') string
               		write(unit=0, fmt=*) param
            	else
               		write(unit=0, fmt='(a)') string(1:strlen)
            	end if
         end if

      	return
end
