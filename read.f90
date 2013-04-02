subroutine Read_Real(var,count)
      implicit none

      INTEGER :: count
      real(8) :: var(count)
      INTEGER :: file_in,echo_out
      integer :: i,j,k,linelen
      COMMON/filename/ file_in,echo_out
      CHARACTER(len=132) :: line

!  Read lines in until reaching a non-blank line that is not commented out
      i = 0
 10   CONTINUE
      READ(file_in,'(a132)',err=911) line
      linelen = 133
 11   CONTINUE
      linelen = linelen - 1
      IF (linelen .GT. 0) then
          if(line(linelen:linelen) .EQ. ' ') GOTO 11
      endif
      
      IF (echo_out .GT. 0)WRITE(echo_out,15)(line(k:k),k=1,linelen)
 15 FORMAT(132a1)
        
      IF (INDEX('!@#%/',line(1:1)) .NE. 0) GOTO 10
      IF (linelen .LE. 0) GOTO 10
        
!  Convert the string beyond the equals sign into the var(count) array

      j = INDEX(line,'=') + 1
 20   CONTINUE
      READ(line(j:linelen),*,err=40) var(i+1)
      i = i + 1
      IF (i .LT. count) THEN
 30     CONTINUE
          k = INDEX(line(j:linelen),' ')
          j = j + k
          IF (k .EQ. 1) GOTO 30
          IF (k .EQ. 0) GOTO 10
          GOTO 20
       ENDIF

 40   continue
      do while (i .lt. count)
        var(i+1) = var(i)
        i = i+1
      enddo


!  Now distribute the array var(count) to all nodes

      RETURN

!  Read Error routine

 911  CONTINUE
      print*,'Read_Real error - ',i,count
 !     STOP'Read_Real error'

end subroutine Read_Real


!********************************************************************
subroutine Read_Int(ivar,count)
      implicit none
      INTEGER :: count,ivar(count),file_in,echo_out
      COMMON/filename/ file_in,echo_out
      CHARACTER(len=132) :: line
      integer :: i,j,k,linelen

! Read lines in until reaching a non-blank line that is not commented out

        i = 0
   10   CONTINUE
      READ(file_in,'(a132)',err=911) line
      linelen = 133
   11 CONTINUE
      linelen = linelen - 1
      IF (linelen > 0) then
          if(line(linelen:linelen) == ' ') GOTO 11
      endif
      IF (echo_out > 0)WRITE(echo_out,15)(line(k:k),k=1,linelen)
   15 FORMAT(132a1)
 
      IF (INDEX('!@#%/',line(1:1)) .NE. 0) GOTO 10
      IF (linelen .LE. 0) GOTO 10

! Convert the string beyond the equals sign into the ivar(count) array

      j = INDEX(line,'=') + 1
   20 CONTINUE
      READ(line(j:linelen),*,err=10) ivar(i+1)
      i = i + 1
      IF (i .LT. count) THEN
   30    CONTINUE
         k = INDEX(line(j:linelen),' ')
         j = j + k
         IF (k .EQ. 1) GOTO 30
         IF (k .EQ. 0) GOTO 10
         GOTO 20
      ENDIF
! Now distribute the array ivar(count) to all nodes


      RETURN

! Read Error routine

  911 CONTINUE
      print*,'Read_Int error - ',i,count
!      STOP'Read_Int error'

end subroutine Read_Int
