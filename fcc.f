      logical function bconvert(noeof)
      integer noeof

      if (noeof.eq.0) then 
         bconvert=.false.
      else
         bconvert=.true.
      endif
c      write(*,*) 'bconvert=',bconvert
      return
      end
