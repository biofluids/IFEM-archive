c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c modified form io11.f file to generate 
c ensight fluid field file
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zfem_ensFluid(u,v,w,p,klok,mnlatwr1,nskipwr1,
     $     mxlatwr1,mnlatwr2,nskipwr2,mxlatwr2,mnlatwr3,
     $     nskipwr3, mxlatwr3,dmnac1, dmnac2, dmnac3,dmxac1,
     $     dmxac2,dmxac3,dmnlat1,dmnlat2,dmnlat3,dmxlat1,
     $     dmxlat2, dmxlat3)

      implicit real*8 (a-h,o-z)
      include 'main_common'
      include 'r_common'

      integer dmnac1
      integer dmnac2
      integer dmnac3
      integer dmxac1
      integer dmxac2
      integer dmxac3
      integer dmnlat1
      integer dmnlat2
      integer dmnlat3
      integer dmxlat1
      integer dmxlat2
      integer dmxlat3

      dimension  u( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  v( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  w( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  p( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  vort( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)

      character*5  fileroot
      character*12 name_file1
      character*12 name_file2
      character*12 name_file3

      parameter (ifileunit = 15)

c%%%%%%%%%%%%%%%%%%%%%%%
      integer part
      data part/2/


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

      write(name_file1,'(A7, A5)')  'fem.vel', fileroot
      write(name_file2,'(A7, A5)')  'fem.pre', fileroot
      write(name_file3,'(A7, A5)')  'fem.vor', fileroot

      write(*,*) 'file name is ', name_file1
      write(*,*) 'file name is ', name_file2
      write(*,*) 'file name is ', name_file3
 
      nskip = 2
      open(ifileunit, file=name_file1, form='formatted')
      write(ifileunit, '(A)') 'fluid field: velocity'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'block'
      
       do 2000 k = mnk,mxk, nskip
          do 2001 j =  0,dmxlat2
             do 2002 i= 3+nbou,dmxlat1-nbou, nskip
               write(ifileunit,'(E12.5)') u(i,0,k)*unit_velocity
 2002        continue
 2001     continue
 2000  continue

       do 2005 k = mnk,mxk, nskip
          do 2006 j =  0,dmxlat2
             do 2007 i= 3+nbou,dmxlat1-nbou, nskip
               write(ifileunit,'(E12.5)') v(i,0,k)*unit_velocity
 2007       continue
 2006    continue
 2005 continue

       do 2010 k = mnk,mxk, nskip
          do 2011 j =  0,dmxlat2
             do 2012 i= 3+nbou,dmxlat1-nbou, nskip
               write(ifileunit,'(E12.5)') w(i,0,k)*unit_velocity
 2012        continue
 2011     continue
 2010  continue
      close(ifileunit)


      open(ifileunit, file=name_file2, form='formatted')
      write(ifileunit, '(A)') 'fluid field: pressure'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'block'
      
       do 2100 k = mnk,mxk, nskip
          do 2101 j =  0,dmxlat2
             do 2102 i= 3+nbou,dmxlat1-nbou, nskip
               write(ifileunit,'(E12.5)') 
     $               p(i,0,k)*unit_pressure
 2102        continue
 2101     continue
 2100  continue
      close(ifileunit)

      j=0
c++++++++
c     calculate vorticity
c++++++++
         do 400 k = mnk, mxk
            do 403 i= dmnlat1+3, dmxlat1
               deltau = (u(i,j,k+1) - u(i,j,k-1))*unit_velocity
               deltaz = 2.0d0*unit_length
               deltaw = (w(i+1,j,k) - w(i-1,j,k))*unit_velocity
               deltax = 2.0d0*unit_length
               vort(i,j,k) = deltau/deltaz - deltaw/deltax
 403        continue
 400     continue

      open(ifileunit, file=name_file3, form='formatted')
      write(ifileunit, '(A)') 'fluid field: vorticity'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'block'
      
       do 2200 k = mnk,mxk, nskip
          do 2201 j =  0,dmxlat2
             do 2202 i= 3+nbou,dmxlat1-nbou, nskip
               write(ifileunit,'(E12.5)') 
     $               vort(i,0,k)
 2202        continue
 2201     continue
 2200  continue
      close(ifileunit)



      return
      end




