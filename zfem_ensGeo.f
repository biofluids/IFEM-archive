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
      integer part, part2, sizeX, sizeY, sizeZ, i, j, k, NUM, nskip
      data part/1/, part2/2/


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
      write(6,*) 'file name is ', file_name 

      open(i_file_unit, file=file_name, form='formatted')
      write(i_file_unit, *) 'This is the ensight format geometry file'
      write(i_file_unit, *) 'This is the ensitht format geometry file'
      

      npf = 1
      nf  = 1
      iif  = mn_fiber_alloc

      ipt = dlptlocal_head

       do 100 npt = 1, dlptlocal_number
         
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
          
 100   continue
      
       nfiber_file = nf - 1

c%%%%%%%%%%%%%%%%%%%%
c using quadrangle element
c%%%%%%%%%%%%%%%%%%%%
       write(i_file_unit, *) 'node id given'
       write(i_file_unit, *) 'element id given'
       write(i_file_unit, *) 'coordinates'
       write(i_file_unit, '(I10)') nnda

cc--> node id, x, y, x
cc-->
101	format(i8,3e12.5)
	!write structure coordinates
       do i=1, nnda  
         write(i_file_unit, 101) i,(coor(i,1)+dis(1,i))*unit_length,
     +	coord_pt(2,i)*unit_length,(coor(i,2)+dis(2,i))*unit_length
	 enddo

	!write fluids coordinates
	do i=1,nn
	  write(i_file_unit,101) i+nnda,xn(1,i),xn(2,i),xn(3,i)
	enddo

	!write structure part - node numbers
102	format('part',i2)
      write(i_file_unit, 102) part
      write(i_file_unit, *) 'Structure Model'
	write(i_file_unit,*) 'point'
	write(i_file_unit,*) nnda

	do i=1,nnda
	write(i_file_unit,103) i,i
	enddo
103	format(2i8)

	!write fluids part - node numbers
      write(i_file_unit, 102) part2
      write(i_file_unit, *) 'Fluid Model'
	write(i_file_unit,*) 'point'
	write(i_file_unit,*) nn
	do i=1,nn
	write(i_file_unit,103) i+nnda,i+nnda
	enddo

	!write structure part element connectivity
      write(i_file_unit, *) 'quad4'  !element type
      write(i_file_unit, '(I10)')  numela   ! number of elements
      do i=1, numela
         write(i_file_unit, 105) i, (nea(i,j),j=1,nis) !element connectivity
	enddo
105	format(5i10)

	!write fluids part element connectivity
      write(i_file_unit, *) 'tetra4'  !element type
      write(i_file_unit, '(I10)')  ne   ! number of elements
      do i=1, ne
        write(i_file_unit, 105) i, (ien(j,i)+nnda,j=1,nen) !element connectivity
	enddo

      close(i_file_unit)

      return
      end  
