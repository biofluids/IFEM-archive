c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c modified form io11.f file to generate 
c ensight fluid field file
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zibm_ensFluid_old(u,v,w,p,klok,mnlatwr1,nskipwr1,
     $     mxlatwr1,mnlatwr2,nskipwr2,mxlatwr2,mnlatwr3,
     $     nskipwr3, mxlatwr3,dmnac1, dmnac2, dmnac3,dmxac1,
     $     dmxac2,dmxac3,dmnlat1,dmnlat2,dmnlat3,dmxlat1,
     $     dmxlat2, dmxlat3)

      implicit real*8 (a-h,o-z)
c     include 'common' 
      include 'main_common'

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

      character*5  fileroot
      character*12 name_file1
      character*12 name_file2

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

      write(name_file1,'(A7, A5)')  'ibm.vel', fileroot
      write(name_file2,'(A7, A5)')  'ibm.pre', fileroot

      write(*,*) 'file name is ', name_file1
      write(*,*) 'file name is ', name_file2
 

      open(ifileunit, file=name_file1, form='formatted')
      write(ifileunit, '(A)') 'fluid field: velocity'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'block'

      do 2000 k = mnlatwr3, mxlatwr3, nskipz
         do 2001 j = mnlatwr2, mxlatwr2, nskipy
            do 2002 i = mnlatwr1, mxlatwr1, nskipx
               write(ifileunit,'(E12.5)') u(i,j,k)*unit_velocity
 2002        continue
 2001     continue
 2000  continue

      do 2005 k = mnlatwr3, mxlatwr3, nskipz
         do 2006 j = mnlatwr2, mxlatwr2, nskipy 
            do 2007 i = mnlatwr1, mxlatwr1, nskipx
               write(ifileunit,'(E12.5)') v(i,j,k)*unit_velocity
 2007       continue
 2006    continue
 2005 continue

      do 2010 k = mnlatwr3, mxlatwr3, nskipz
         do 2011 j = mnlatwr2, mxlatwr2, nskipy 
            do 2012 i = mnlatwr1, mxlatwr1, nskipx
               write(ifileunit,'(E12.5)') w(i,j,k)*unit_velocity
 2012        continue
 2011     continue
 2010  continue
      close(ifileunit)


      open(ifileunit, file=name_file2, form='formatted')
      write(ifileunit, '(A)') 'fluid field: pressure'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'block'

       do 2100 k = mnlatwr3, mxlatwr3,nskipz
         do 2101 j = mnlatwr2, mxlatwr2,nskipy
            do 2102 i = mnlatwr1, mxlatwr1, nskipx
               write(ifileunit,'(E12.5)') p(i,j,k)*unit_pressure
 2102        continue
 2101     continue
 2100  continue
      close(ifileunit)

      return
      end




