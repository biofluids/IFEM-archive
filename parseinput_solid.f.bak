      subroutine parseinput_solid
      implicit none
      
      include 'r_common'
      include 'main_common'
      include 'ibd0_nonzerocfddiv_vars.fh'
      
           
      	
	character*32 key, keyaux,test
	character*8 onoff
	logical fctrl, getkey, isatty
	logical enough
	data enough /.false./
	  
c       command loop

	do while (.not.enough)
	   enough = .not.getkey(key)
	   if ((key.eq.'abort').or.(key.eq.'exit')) then
	      call exit(1)
	   else if ((key.eq.'done').or.(key.eq.'quit')) then
	      enough = .true.
c       ---------------------------------------------------------------------c
c---- mesh information
	   else if (key.eq.'xshift') then !..................shift on x axis
	      call getreal(xshift)
	   else if (key.eq.'zshift') then !..................shift on z axis
	      call getreal(zshift)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   else
	      
	      if ((.not.enough).and.(key.ne.' ').and.(key(1:1).ne.'#'))
	1	call error("parseinput: strange keyword "//key, -999, .false.)
	   end if
	   
	end do
      
      
      
      return
      end
      