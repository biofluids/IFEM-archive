ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine fclear(f,m)
       real* 8 f(m)
       do i=1,m
       f(i) = 0.0
       enddo
       return
       end
