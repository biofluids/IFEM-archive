ccccccccccccccccccc
c  tecplotpost.f  c
c  G. Wagner      c
ccccccccccccccccccc

      subroutine tecplotpost
      include "global.h"
      include "malloc.h"

      integer ien(nen,nec)  
      real* 8 xn(nsd,nnc), dn(ndf,nnc), dd(ndf,nnc)
      integer plotmax, plotstep
      pointer (ienptr,ien),(xnptr, xn)
      pointer (dnptr,dn),(dfnptr,dfn),(ddptr,dd)
      character*24 title
      character*9 file
      character*16 et

      ienptr = malloc(nen*nec*isize)
      xnptr = malloc(nnc*nsd*fsize)
      dnptr = malloc(nnc*ndf*fsize)
      dfnptr = malloc(nnc*fsize)
      ddptr = malloc(nnc*(ndf)*fsize)

      ibase = ichar('0')        !! integer value for char '0'

      call error("read ien",-999,.false.)
      call readien(ien)

      call error("read x and rng",-999,.false.)
      call readx(xn)

      call error("write tecout.plt",-999,.false.)
      title = '3D Tet Mesh'
      call headerout(title,nsd,ndf)

      plotmax = 2
      plotstep = 1

      file = "data.0000"

      do i = 0, plotmax, plotstep
        file = "data.0000"

        i4 = i/1000
        i3 = (i-i4*1000)/100
        i2 = (i-i4*1000-i3*100)/10
        i1 = (i-i4*1000-i3*100-i2*10)/1
        i4 = i4 + ibase
        i3 = i3 + ibase
        i2 = i2 + ibase
        i1 = i1 + ibase
        file (6:6) = char(i4)
        file (7:7) = char(i3)
        file (8:8) = char(i2)
        file (9:9) = char(i1)

        write (*,*) file

        call postin(dn,dd,file)
        call postout(i,xn,dd,nsd,ndf,nn,ne,nen)
        if (i == 0) call meshout(ien,nen,nec)
      end do
      
      end
      


