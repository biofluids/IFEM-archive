      integer function calcbinxyz(xcoord,ycoord,zcoord,lenbin1,
     $     lenbin2,lenbin3,nbin1,nbin2,nbin3)
      implicit real*8 (a-h,o-z)

      ibin1 = xcoord / lenbin1
      ibin2 = ycoord / lenbin2
      ibin3 = zcoord / lenbin3

      calcbinxyz = ibin1 + nbin1*(ibin2 + nbin2 * ibin3)

      return
      end


