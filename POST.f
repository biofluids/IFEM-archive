        parameter(nn=2072841,ne=2040000,ndf=5,nsd=3,nen=8)
        parameter(neface = 6, nnface = 4)

        real* 8   d(ndf,nn)
        character*16 dummy
        integer ien(nen,ne), rng(neface,ne),nd(nn)
        integer ien2(nen/2,ne),rng2(neface/2,ne)
        real* 8 out(2,nn)
        integer hacemap(neface,nnface)
        data hacemap /
     &  1, 1, 2, 3, 4, 5,
     &  4, 2, 3, 4, 1, 6,
     &  3, 6, 7, 8, 5, 7,
     &  2, 5, 6, 7, 8, 8/

        character*9 datain
        integer i1,i2,i3,i4
        ibase = ichar('0')        !! integer value for char '0'

        call ewd_open("mien", ifd)
        call ewd_read(ifd, ien2, nen*ne*4)
        call ewd_close(ifd)
        ierr = ieg2cray(1,nen*ne,ien2,0,ien)

        call ewd_open("mrng", ifd)
        call ewd_read(ifd, rng2, neface*ne*4)
        call ewd_close(ifd)
        ierr = ieg2cray(1,neface*ne,rng2,0,rng)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  	  call ewd_open("DATA",id)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do in=1,nn
        nd(in) = 1
        enddo

      do ie=1,ne
        do ic=1,neface
        iflag = 0
        if(rng(ic,ie).ne.0) iflag=1
        do in=1,nnface
        inl = hacemap(ic,in)
        node = ien(inl,ie)
        if(iflag.eq.1) nd(node) = 0
        enddo
        enddo
        enddo
        write(6,*)'pass'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        write(6,*) 'last record registered= ?'
        read(5,*) nrec 
        write(6,*) 'increment= ?'
        read(5,*) irec 

        datain  = "data.0000"

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=0,nrec,irec   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        i4 = i/1000
        i3 = (i-i4*1000)/100
        i2 = (i-i4*1000-i3*100)/10
        i1 = (i-i4*1000-i3*100-i2*10)/1
        i4 = i4 + ibase
        i3 = i3 + ibase
        i2 = i2 + ibase
        i1 = i1 + ibase
        datain (6:6) = char(i4)
        datain (7:7) = char(i3)
        datain (8:8) = char(i2)
        datain (9:9) = char(i1)

        call ewd_open(datain,ifd)
        call ewd_read(ifd,d,ndf*nn*8)
        call ewd_close(ifd)
        ierr = ieg2cray(8,ndf*nn,d,0,d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do in=1,nn
        d(5,in) = nd(in)*d(5,in)
        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	 do in=1,nn
	 out(1,in) = d(4,in)
	 out(2,in) = d(5,in)
	 enddo

	ierr = cray2ieg(8,2*nn,out,0,out)
	call ewd_write(id,out,2*nn*8)
	write(6,*) i

      enddo
	call ewd_close(id)
      stop
      end
