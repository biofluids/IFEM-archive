      subroutine  rdpc(name_file,klok,dnext_pt,
     $     dlfreeel_number,dlfreeel_head,dlfreeel_tail,
     $     dlptlocal_number,dlptlocal_head,dlptlocal_tail,
     $     pt_iptexp,coord_pt,attrib_fcu,
     $     num_fiber, num_point)
      implicit real*8 (a-h,o-z)
      include 'r_common' 
      include 'main_common'            

      integer dnext_pt( mn_pt_alloc:mx_pt_alloc )
      integer dlfreeel_number
      integer dlfreeel_head
      integer dlfreeel_tail
      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail
      integer pt_iptexp( mn_pt_alloc:mx_pt_alloc )

      dimension coord_pt( ix:iz, mn_pt_alloc:mx_pt_alloc )
      dimension attrib_fcu( mn_pt_alloc:mx_pt_alloc,
     $     mn_attr_fcu:mx_attr_fcu )

      integer num_point(max_point)

      character*11  fileroot
      character*8 name_file

      parameter (ifileunit = 15)

      if (n_ibmfem .eq. 0) then
         open(ifileunit, file=name_file,status='old')         
         read(ifileunit,*) nptfile, ng1
      endif

      nptibm=nptfile

c     assign fem nnd
c++++

      if (n_ibmfem .eq. 1) then
         open(ifileunit, file=name_file,status='old')         
         read(ifileunit,*) nptibm, ng1
         nptfem = nnd+ntether
         nptfile = nptfem +nptibm
         write(*,*) 'nptibm=',nptibm
         write(*,*) 'nptfem=',nptfem
c         nptfile=nnd+ntether
      endif
c++++


      if (nptfile .gt. n_pt_exp) then
         stop
      endif

      call alloclist(nptfile,mn_pt_alloc,mx_pt_alloc,
     $     null_pt,dnext_pt,dlfreeel_number,dlfreeel_head,
     $     dlfreeel_tail,dlptlocal_number, dlptlocal_head,
     $     dlptlocal_tail)

      ipt = dlptlocal_head

c
c     fem to ibm coordinate transformation
c
      if (n_ibmfem .eq. 1) then
         do 301 i=1,nnd
c            coor(i,1)=coor(i,1)+xnmag*xshift
            coor(i,1)=coor(i,1)+xshift
            coord_pt(1,i)=coor(i,1)
            coord_pt(2,i)=xlay(i)
c            coor(i,2)=coor(i,2)+xnmag*zshift
            coor(i,2)=coor(i,2)+zshift
            coord_pt(3,i)=coor(i,2)
 301     continue
         njc=0
         
         do 37 i=1,numgb
ccccccccccccccccccccccccccccccccccccccccccc
c     Fixed 1,2
ccccccccccccccccccccccccccccccccccccccccccc
            if (ndirgb(i) .eq. 111111) then
               do 38 j=1,numdir(i)
                  nl=nodegb(i,j)
c
                  if (nxt(1,nl) .ne. 0) then
                     njc=njc+1                  
                     coord_pt(2,nnd+njc)=coord_pt(2,nl)
                     coord_pt(3,nnd+njc)=coord_pt(3,nl)
                     coord_pt(1,nnd+njc)=coord_pt(1,nl)
     $                    +xtedis*nxt(1,nl)
                  endif
c     
                  if (nxt(3,nl) .ne. 0) then
                     njc=njc+1
                     coord_pt(2,nnd+njc)=coord_pt(2,nl)
                     coord_pt(1,nnd+njc)=coord_pt(1,nl)
                     coord_pt(3,nnd+njc)=coord_pt(3,nl)
     $                    +xtedis*nxt(3,nl)
                  endif
 38            continue
            endif
 37      continue

c
c     assign pt_iptexp(ipt)
c
c++++ modified on Oct.22
c         do 101 npt =1,nptfile
c            pt_iptexp(ipt)=npt
 101     continue
      endif

         nptfilea=nptfile
         do npt=1, nptfem
            pt_iptexp(ipt)=npt
            ipt = dnext_pt(ipt)
         enddo

         num_temp = 0
         num_fiber  = 0
         do npt = nptfem+1, nptfile
            read(ifileunit,*) ntss,(coord_pt(idim,ipt),idim=ix,iz),
     $           (attrib_fcu(ipt,iattrfcu),iattrfcu=mn_attr_fcu,
     $           mx_attr_fcu)
            pt_iptexp(ipt)=npt
            num_temp = num_temp +1
            if(attrib_fcu(ipt,1) .eq. -1) then
               num_fiber = num_fiber +1
               num_point(num_fiber) = num_temp
               num_temp = 0
            endif
            ipt = dnext_pt(ipt)
            enddo
c+++++++



c
c     immersed boundary read_in
c

      if (n_ibmfem .ne. 1) then
         num_temp = 0

         write(*,*) nptfile
            
         nptfilea=nptfile
         do 100 npt = 1, nptfile
            
            read(ifileunit,*) ntss,(coord_pt(idim,ipt),idim=ix,iz),
     $           (attrib_fcu(ipt,iattrfcu),iattrfcu=mn_attr_fcu,
     $           mx_attr_fcu)
            pt_iptexp(ipt)=npt
            num_temp = num_temp +1
            if(attrib_fcu(ipt,1) .eq. -1) then
               num_fiber = num_fiber +1
               num_point(num_fiber) = num_temp
               num_temp = 0
            endif

            ipt = dnext_pt(ipt)
            
            if (mod(npt, 50000) .eq. 0) then
            endif
            
 100     continue
      endif
 210  format(i6,1x,3f13.8,1x,i4,2f13.8)
c      write(*,*) 'ipt', ipt
      return
      end







