      integer function calcperbinxyz(xcoord,ycoord,zcoord,
     $     lenbin1,lenbin2,lenbin3,nbin1,nbin2,nbin3)
      implicit real*8 (a-h,o-z)

      integer calcbinxyz
      integer    floor


      nn_kludge_mod  = nbin1*lenbin1*nbin2*lenbin2*nbin3*lenbin3
      rnn_kludge_mod  = nn_kludge_mod

      ncoord1 = int(xcoord + rnn_kludge_mod)
      ncoord2 = int(ycoord + rnn_kludge_mod)
      ncoord3 = int(zcoord + rnn_kludge_mod)

c      nx = mod(ncoord1, nbin1*lenbin1)
c      ny = mod(ncoord2, nbin2*lenbin2)
c      nz = mod(ncoord3, nbin3*lenbin3)

c      xxcoord = nx
c      yycoord = ny
c      zzcoord = nz

c      calcperbinxyz =calcbinxyz(xxcoord,yycoord,zzcoord,lenbin1,
c     $     lenbin2,lenbin3,nbin1,nbin2,nbin3)

      return
      end
