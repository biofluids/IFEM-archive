      program convert
      parameter (nn=64,ne=27,nsd=3,nen=8,neface=6)
      integer i,m,ie,n,ien(nen,ne),ien2(nen,ne),rng(neface,ne)
      character*10 letter
      real x(nsd,nn),alpha,x3
      open(unit=10, file='mxyz.dat')
      open(unit=11, file='mien.dat')
      open(unit=15, file='mrng.dat')
      open(unit=16, file='mxyz')
      open(unit=17, file='mien')
      open(unit=18, file='mrng')
12    format(f16.8,f16.8,f16.8)
      call readx(x,nsd,nn)

      do i=1,nn
         write(10,12) i,x(1,i), x(2,i), x(3,i)
      enddo
      
 13   format(i9,i9,i9,i9)
 14   format(i9,i9,i9,i9,i9,i9)
      call readien(ien2,nen,ne)
      do i=1,ne
         write(11,13) ien2(1,i),ien2(2,i),ien2(3,i),ien2(4,i)
      enddo
      call readrng(rng,neface,ne)
      do i=1,ne
	write(15,14) rng(1,i),rng(2,i),rng(3,i),rng(4,i),rng(5,i),rng(6,i)
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine writex(xn,nsd,nn)
        real* 8 xn(nsd,nn),maxrecl
        maxrecl=nsd*nn
        call ewd_open("mxyz", ifp)
        call ewd_write(ifp, xn, nsd*nn*8)
        call ewd_close(ifp)
        return
        end
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine writeien(ien,nen,ne)
        integer ien(nen,ne),maxrecl
        maxrecl=nen*ne
        call ewd_open("mien", ifp)
        call ewd_write(ifp, ien, nen*ne*8)
        call ewd_close(ifp)
        return
        end
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine writerng(rng,neface,ne)
        integer rng(neface,ne),maxrecl
        maxrecl=neface*ne
        call ewd_open("mrng", ifp)
        call ewd_write(ifp, rng, neface*ne*8)
        call ewd_close(ifp)
        return
        end

        subroutine readx(xn,nsd,nn)
        real* 8 xn(nsd,nn),maxrecl
        maxrecl=nsd*nn
        call ewd_open("mxyz", ifp)
        call ewd_read(ifp, xn, nsd*nn*8)
        call ewd_close(ifp)
        return
        end
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine readien(ien,nen,ne)
        integer ien(nen,ne),maxrecl
        maxrecl=nen*ne
        call ewd_open("mien", ifp)
        call ewd_read(ifp, ien, nen*ne*8)
        call ewd_close(ifp)
        return
        end
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine readrng(rng,neface,ne)
        integer rng(neface,ne),maxrecl
        maxrecl=neface*ne
        call ewd_open("mrng", ifp)
        call ewd_read(ifp, rng, neface*ne*8)
        call ewd_close(ifp)
        return
        end





