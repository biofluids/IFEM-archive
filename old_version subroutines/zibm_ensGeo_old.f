c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c modified from io6.f(wrxf() subroutine) to 
c generate ensight geometry format file
c last modified data: Sep. 28, 2001
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine  zibm_ensGeo(klok,td,dnext_pt,
     $     dlptlocal_number,
     $     dlptlocal_head, 
     $     dlptlocal_tail,
     $     pt_iptexp,
     $     coord_pt, attrib_fcu)

      implicit real*8 (a-h,o-z)
      include 'main_common'

      integer dnext_pt( mn_pt_alloc:mx_pt_alloc )
      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail

      integer pt_iptexp( mn_pt_alloc:mx_pt_alloc )

      dimension coord_pt( ix:iz, mn_pt_alloc:mx_pt_alloc )
      dimension attrib_fcu( mn_pt_alloc:mx_pt_alloc,
     $     mn_attr_fcu:mx_attr_fcu )

      parameter (rsentinel_plus_fiber_end = -1.0d0)

      parameter (n_fiber_alloc  = max_point)
      parameter (mn_fiber_alloc = 0)
      parameter (mx_fiber_alloc = mn_fiber_alloc + (n_fiber_alloc-1))

      dimension  npoint_fiber(mn_fiber_alloc:mx_fiber_alloc)

      character*5 fileroot
      character*12 file_name

      parameter (i_file_unit = 15)

c%%%%%%%%%%%%%%%%%%%%%%
c used for ensight
c%%%%%%%%%%%%%%%%%%%%%
      integer part, part2, size, i, j, k, NUM, nskip
      data part/1/, part2/2/, size/16/

c      call  createfileroot(klok,  fileroot)

c      integer    klok
c      character*5  fileroot

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

      write(file_name,'(A7, A5)')  'ibm.geo', fileroot
      write(6,*) 'file name is ', file_name 

      open(i_file_unit, file=file_name, form='formatted')
      write(i_file_unit, *) 'This is the ensight format geometry file'
      write(i_file_unit, *) 'This is the ensight format geometry file'

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
c test using point element
c%%%%%%%%%%%%%%%%%%%%
       write(i_file_unit, *) 'node id given'
       write(i_file_unit, *) 'element id given'
       write(i_file_unit, *) 'part'
       write(i_file_unit, '(I10)') part
       write(i_file_unit, *) 'Blood Vessel Model'
       write(i_file_unit, *) 'coordinates'
       write(i_file_unit, '(I10)') dlptlocal_tail

       do i=1, dlptlocal_tail
         write(i_file_unit, '(I10)') i
	 enddo

       do j=ix, iz
          do i=1, dlptlocal_tail
             write(i_file_unit, '(e12.5)') coord_pt(j, i)
		enddo
	 enddo

       write(i_file_unit, *) 'point'
       write(i_file_unit, '(I10)')  dlptlocal_tail
       
       do i=1, dlptlocal_tail
         write(i_file_unit, '(I10)') i
	 enddo
       do i=1, dlptlocal_tail
         write(i_file_unit, '(I10)') i
	 enddo

       NUM = n_lat_exp1-1
      
       write(i_file_unit, *) 'part'
       write(i_file_unit, '(I10)') part2
       write(i_file_unit, *) 'fluid field'
       write(i_file_unit, *) 'block iblanked'
       write(i_file_unit, '(3I10)') NUM/nskipx+1, 
     $      NUM/nskipy+1, NUM/nskipz+1

       do k = 0,NUM, nskipz
         do j =  0,NUM, nskipy
            do i =  0,NUM, nskipx
               write(i_file_unit ,'(e12.5)') i*1.0
		  enddo
	   enddo
	enddo

       do k = 0,NUM, nskipz 
         do j =  0,NUM, nskipy 
            do i =  0,NUM, nskipx
               write(i_file_unit,'(e12.5)') j*1.0
		  enddo
	   enddo
	 enddo

       do k = 0,NUM, nskipz 
         do j =  0,NUM, nskipy 
            do i =  0,NUM, nskipx
               write(i_file_unit, '(e12.5)') k*1.0
		  enddo
	   enddo
	 enddo
       j=1
       do i = 1, (NUM/nskipx+1)*(NUM/nskipx+1)*(NUM/nskipx+1)    
          write(i_file_unit, '(I10)') j
	 enddo
      
      close(i_file_unit)

      return
      end  


