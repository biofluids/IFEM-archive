

      subroutine ReadALine(ifr,ifw,linestr)
c*** 
c*** input parameters:
c    ifr:     the input file handle
c    ifw:     the output(.ECHO) file handle
c*** output parameter:
c    linestr: one line proper input from the input file
c
c*** NOTE:
c    1. the linestr should be long enough, e.g. the lenght > 80 in most cases
c    2. the comments in the input file will be written to the output file, 
c       but the empty line will be ignored.

  
      implicit none       
      integer ifr,ifw
      character*(*) linestr
      
      
      	!** local var
      logical readOK
      integer length,i
      character*1 c, c_TAB
      
      c_TAB=char(9)
      
      readOK=.false.
      
      do while ( .not. readOK )
      	 read(ifr,'(A)') linestr
      
         length=len(linestr)
      
         do 100 i=1,length
            c=linestr(i:i)
            if (    (c.eq.'#') .or. (c.eq.'!') 
     &         .or. (c.eq.'c') .or. (c.eq.'C') ) then
         			! comment statement, 
				! read another line from file
               write(ifw,'(1x,A80)')  linestr
c               exit
               goto 101
                                
            elseif ( (c.eq.' ').or. (c.eq.c_TAB) ) then
c                cycle
               goto 100		! space, ignore it.
         
            else
               readOK=.true.	! this line is proper input
c               exit
               goto 101
            endif
100      continue

101      continue 

      enddo	! endwhile
      
      end
