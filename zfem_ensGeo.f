c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c This subroutine generates Ensight6 geometry file
c Lucy Zhang
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine  zfem_ensGeo(klok,td,dnext_pt,
     $     dlptlocal_number,
     $     dlptlocal_head, 
     $     dlptlocal_tail,
     $     pt_iptexp,
     $     coord_pt, attrib_fcu,ien,xn)

      implicit real*8 (a-h,o-z)
      
      include "r_common"
      include "main_common"
	include "global.h"

      integer dnext_pt( mn_pt_alloc:mx_pt_alloc )
      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail

      integer    klok
      real*8       td

      integer pt_iptexp( mn_pt_alloc:mx_pt_alloc )

      real*8 coord_pt( ix:iz, mn_pt_alloc:mx_pt_alloc )

      real*8 attrib_fcu( mn_pt_alloc:mx_pt_alloc,
     $     mn_attr_fcu:mx_attr_fcu )

      integer    nptfile
      integer    ng

      real*8       rsentinel_plus_fiber_end
      parameter (rsentinel_plus_fiber_end = -1.0d0)

      integer    n_fiber_alloc
      integer    mn_fiber_alloc
      integer    mx_fiber_alloc
      parameter (n_fiber_alloc  = max_point)
      parameter (mn_fiber_alloc = 0)
      parameter (mx_fiber_alloc = mn_fiber_alloc + (n_fiber_alloc-1))

      integer    nfiber_file
      integer    nf
      integer    iif
      integer    npf

      character*5 fileroot
      character*12 file_name

      integer    i_file_unit
      parameter (i_file_unit = 15)

      integer    idim
      integer    iattrfcu

      integer    ipt
      integer    npt

      real*8       time_file
c%%%%%%%%%%%%%%%%%%%%%%
c used for ensight
c%%%%%%%%%%%%%%%%%%%%%
      integer sizeX, sizeY, sizeZ, i, j, k, NUM, nskip

	real*8  xn(nsd,nn)
	integer ien(nen,ne)

 200  format(a5   )
 201  format(a4,i1)
 202  format(a3,i2)
 203  format(a2,i3)
 204  format(a1,i4)
 205  format(   i5)

      if (klok .eq. 0) then
         write(fileroot, 200) '00000'
      elseif (klok .lt. 10) then
         write(fileroot, 201) '0000',klok
      elseif (klok .lt. 100) then
         write(fileroot, 202) '000' ,klok
      elseif (klok .lt. 1000) then
         write(fileroot, 203) '00'  ,klok
      elseif (klok .lt. 10000) then
         write(fileroot, 204) '0'   ,klok
      elseif (klok .lt. 100000) then    
         write(fileroot, 205)   ''  ,klok
      else
         write(0,*) 'klok .ge. 100000: modify subroutine createfileroot'
         call exit(1)
      endif


      write(file_name,'(A7, A5)')  'fem.geo', fileroot
      write(6,*) 'writing... ', file_name 

      open(i_file_unit, file=file_name, form='formatted')
      write(i_file_unit, *) 'This is the ensight format geometry file'
      write(i_file_unit, *) 'This is the ensitht format geometry file'
      

      npf = 1
      nf  = 1
      iif  = mn_fiber_alloc

      ipt = dlptlocal_head

       do npt = 1, dlptlocal_number
         
         npoint_fiber(iif) = npf
         
         if (attrib_fcu(ipt,i_attr_rest0) .ne. 
     $        rsentinel_plus_fiber_end) then
            npf = npf + 1
         else
            npf = 1
            nf = nf + 1
            iif = iif + 1
         endif

         ipt = dnext_pt(ipt)
          
	 enddo
      
       nfiber_file = nf - 1

c%%%%%%%%%%%%%%%%%%%%
c using quadrangle element
c%%%%%%%%%%%%%%%%%%%%
       write(i_file_unit, *) 'node id given'
       write(i_file_unit, *) 'element id given'
       write(i_file_unit, *) 'coordinates'
       write(i_file_unit, '(I8)') nnd+nn

cc--> node id, x, y, x
cc-->
	!write structure coordinates
       do i=1, nnd
         write(i_file_unit, 101) i,(coor(i,1)+dis(1,i)),
     +	(coor(i,2)+dis(2,i)),
     +	(coor(i,3)+dis(3,i))
c	    write(i_file_unit, 101) i, coord_pt(1,i), coord_pt(2,i),
c     +			coord_pt(3,i)

	 enddo
101	format(i8,3e12.5)

	!write fluids coordinates
	do i=1,nn
	  write(i_file_unit,101) i+nnd,xn(1,i),xn(2,i),xn(3,i)
	enddo


	!write structure part - connectivity
      write(i_file_unit, *) 'part 1'
      write(i_file_unit, *) ' Structure Model'
	if (nsd_solid .eq. 0) then
		write(i_file_unit,'(a7)') '  point'
		write(i_file_unit,'(i8)') nnd
		do i=1,nnd
			write(i_file_unit,'(2i8)') i,i
		enddo
	elseif (nsd_solid .eq.3) then
		write(i_file_unit,'(a7)') '  hexa8'
	    write(i_file_unit, '(i8)')  numel   ! number of elements
          do i=1, numel
		  write(i_file_unit,'(9i8)') i, (nea(i,j),j=1,nis) !element connectivity
		enddo
	endif

	!write fluids part element connectivity
      write(i_file_unit, *) 'part 2'
      write(i_file_unit, *) ' Fluid Model'
	if (nen .eq. 4) then
	   write(i_file_unit, *) ' tetra4'
	   write(i_file_unit, '(i8)')  ne   ! number of elements
	   do i=1, ne
		 write(i_file_unit,'(5i8)') i, (ien(j,i)+nnd,j=1,nen) !element connectivity
	   enddo
	elseif (nen .eq. 8) then
	   write(i_file_unit, *) 'hexa8'  !element type
         write(i_file_unit, '(I8)')  ne   ! number of elements
         do i=1, ne
           write(i_file_unit,'(9i8)') i, (ien(j,i)+nnd,j=1,nen) !element connectivity
	   enddo
	endif


      close(i_file_unit)

      return
      end  
