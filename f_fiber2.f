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
      nx = mod(ncoord1, nbin1*lenbin1)
      ny = mod(ncoord2, nbin2*lenbin2)
      nz = mod(ncoord3, nbin3*lenbin3)

      xxcoord = nx
      yycoord = ny
      zzcoord = nz

      calcperbinxyz =calcbinxyz(xxcoord,yycoord,zzcoord,lenbin1,
     $     lenbin2,lenbin3,nbin1,nbin2,nbin3)

      return
      end
