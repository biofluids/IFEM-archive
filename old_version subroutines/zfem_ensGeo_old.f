c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c modified from io6.f(wrxf() subroutine) to 
c generate ensight geometry format file
c last modified data: Sep. 28, 2001
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine  zfem_ensGeo(klok,td,dnext_pt,
     $     dlptlocal_number,
     $     dlptlocal_head, 
     $     dlptlocal_tail,
     $     pt_iptexp,
     $     coord_pt, attrib_fcu)

      implicit real*8 (a-h,o-z)
      
      include "r_common"
      include "main_common"

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
       write(i_file_unit, *) 'part'
       write(i_file_unit, '(I10)') part
       write(i_file_unit, *) 'Structure Model'
       write(i_file_unit, *) 'coordinates'
       write(i_file_unit, '(I10)') nnda

cc--> node id
cc-->
       do i=1, nnda
         write(i_file_unit, '(I10)') i
	 enddo

cc--> x,y,x components
cc-->
       do i=1,nnda
	    write(*, '(e12.5)') 
     $         (coor(i,1)+dis(1,i))*unit_length
          write(i_file_unit, '(e12.5)') 
     $         (coor(i,1)+dis(1,i))*unit_length
	 enddo

       do i=1,nnda
          write(i_file_unit, '(e12.5)') coord_pt(2,i)*unit_length
	 enddo
      
	do i=1,nnda
          write(i_file_unit, '(e12.5)') 
     $        (coor(i,2)+dis(2,i))*unit_length
	enddo

cc--> element type
cc-->
       write(i_file_unit, *) 'quad4'

cc--> number of elements
cc-->
c       write(i_file_unit, '(I10)')  numela*4
       write(i_file_unit, '(I10)')  numela
 
cc--> element id
cc-->      
c       do i=1, numela*4
       do i=1, numela
         write(i_file_unit, '(I10)') i
	 enddo

cc--> element connectivities refer to the coordinate array index, 
cc--> not node id  !!!!!!!
c       do i=1, numela*4
       do i=1, numela
         write(i_file_unit, '(4I10)') (nea(i,j),j=1,nis)
	 enddo

       sizeX = n_lat_exp1-1
       sizeY = n_lat_exp2-1
       sizeZ = n_lat_exp3-1

       nskip = 2
       ntempZ = (mxk-mnk)/nskip+1
       ntempY = sizeY+1
       ntempX = (sizeX-nbou-3-nbou)/nskip+1

       write(i_file_unit, *) 'part'
       write(i_file_unit, '(I10)') part2
       write(i_file_unit, *) 'fluid field'
       write(i_file_unit, *) 'block iblanked'
       write(i_file_unit, '(3I10)') ntempX, ntempY, ntempZ  
       

       do k = mnk,mxk, nskip
          do j =  0,sizeY
             do i= 3+nbou,sizeX-nbou, nskip
               write(i_file_unit ,'(e12.5)') i*unit_length
		   enddo
		enddo
	 enddo

        do k = mnk,mxk, nskip
          do j =  0,sizeY
             do i= 3+nbou,sizeX-nbou, nskip
               write(i_file_unit,'(e12.5)') -j*unit_length
	       enddo
		enddo
	enddo

        do k = mnk,mxk, nskip
          do j =  0,sizeY
             do i= 3+nbou,sizeX-nbou, nskip
                write(i_file_unit, '(e12.5)') k*unit_length
		   enddo
		enddo
	enddo

       j=1
       do i = 1, ntempX*ntempY*ntempZ 
               write(i_file_unit, '(I10)') j
	 enddo

      close(i_file_unit)

      return
      end  




