c*** MyForLib.FOR
c    MISCELLAneous.FOR + SCR.FOR + STRING.FOR

c-----------------------------------------------------------------------
c***  MISCELLAneous.FOR
c...  subroutine AppExt(fhead,ext,fname)
c...  subroutine beep(n)
c-----------------------------------------------------------------------
c***  SCR.FOR
c...  subroutine ScrText(text,x,y,height,RotAngle)
c...  subroutine ScrDrawLine(x1,y1,x2,y2)
c...  subroutine ScrOutputStr(str)
c...  subroutine ScrSetLayerColor(LayerName,icolor)
c...  subroutine ScrSetCurLayer(LayerName)
c...  subroutine ScrSetNewLayer(LayerName)
c-----------------------------------------------------------------------
c*** STRING.FOR

c...  subroutine AppStr(str1,str2,str3)
c...  subroutine DelHeadStr(str,LenDel)
c...  subroutine DelHeadSpace(str)
c...  subroutine ToUpper(str)
c...  subroutine UnSpace(str,istart,iend)
c-----------------------------------------------------------------------

c      subroutine GenInt2Str(int,str)
c      subroutine Int2Str(int,str)


c*** MISCELLAneous.FOR
c*** the subs with miscellaneous functions
c*** by XUE GANG
c*** 1994.11
c-----------------------------------------------------------------------
c...  subroutine checklim(nowval,limval,str)
c...  subroutine AppExt(fhead,ext,fname)
c...  subroutine beep(n)
c-----------------------------------------------------------------------

        subroutine checklim(nowval,limval,str)
        character *(*) str

        if ( nowval .gt. limval ) then
                write(*,*) 'Exceed the limit of ',str
                write(*,*) 'limval is :',limval,
     .                     '  nowval is :',nowval
                stop
        end if
        end


      subroutine AppExt(fhead,ext,fname)
c***  this subroutine is to APPend the EXTension to the fhead with result
c     in fname ( also insert a '.' between fhead and ext )
c*** input parameters :
c     fhead
c     ext
c*** output parameter :
c     fname           : fname = fhead.ext

      character *(*) fhead,ext,fname

        ! local vars
      integer lenfname,lenhead,lenext
      integer i,j

                ! initialize fname first
      lenfname=len(fname)
      do i=1,lenfname
            fname(i:i)=' '
      end do

      lenhead=len(fhead)
      j=0
      do i=1,lenhead
              if ( fhead(i:i) .ne. ' ') then
                      j=j+1
                      fname(j:j)=fhead(i:i)
              end if
      end do
      j=j+1
      fname(j:j)='.'
      lenext=len(ext)
      do i=1,lenext
              if ( (ext(i:i).ne.' ') .and. (ext(i:i).ne.'.') ) then
                      j=j+1
                      fname(j:j)=ext(i:i)
              end if
      end do
      end   !** ends AppExt


      subroutine beep(n)
c*** this subroutine beep n times
      character*1 beepch
      beepch=char(7)
      write(*,*) (beepch,i=1,n)
      end !** ends beep




c*** STRING.FOR
c*** the subs to manipulate string
c*** by XUE GANG
c*** 1994.11
c-----------------------------------------------------------------------
c...  subroutine AppStr(str1,str2,str3)
c...  subroutine DelHeadStr(str,LenDel)
c...  subroutine DelHeadSpace(str)
c...  subroutine ToUpper(str)
c...  subroutine UnSpace(str,istart,iend)
c-----------------------------------------------------------------------


      subroutine AppStr(str1,str2,str3)
c***  this subroutine is to append str1+str2 ==> str3
c     ( DELETE ALL the heading and tailing space )
c*** input parameters :
c     str1
c     str2
c*** output parameter :
c     str3


      implicit real*8 (a-h,o-z)
      character *(*) str1,str2,str3

            !** preset
      len3=len(str3)
      do i=1,len3
            str3(i:i)=' '
      end do

            !** count the non-space character number of str1
      call UnSpace(str1,istart,iend)
      j=0
      do i=istart,iend        !** copy str1 to str3
                              !** (without heading and tailing spaces)
            j=j+1
            if (j.le.len3 ) then
                  str3(j:j)=str1(i:i)
            else
                  write(*,*) '##Warning in sub AppStr'
                  write(*,*) '  len3 is too small'
                  call beep(1)
c                  exit
                  goto 1001
            end if
      end do    !** i
1001  continue

            !** count the non-space character number of str2
      call UnSpace(str2,istart,iend)
      do i=istart,iend        !** copy str1 to str3
                              !** (without heading and tailing spaces)
            j=j+1
            if (j.le.len3 ) then
                  str3(j:j)=str2(i:i)
            else
                  write(*,*) '##Warning in sub AppStr'
                  write(*,*) '  len3 is too small'
                  call beep(1)
c                  exit
                  goto 1002
            end if
      end do    !** i
1002  continue

      end   !** ends AppStr


      subroutine DelHeadStr(str,LenDel)
c*** this subroutine is to delete the heading length characters.
c*** input parameter :
c    str             : the string to be dealed with
c    LenDel          : the length of the str wanted to be deleted
c*** output parameter:
c    str             : with the length heading characters deleted


      implicit real*8 (a-h,o-z)
      character *(*) str

      length=len(str)
      istart=min(length,LenDel)+1
      imove=istart-1
      do j=istart,length
            str(j-imove:j-imove)=str(j:j)
      end do
      do j=length-imove+1,length
            str(j:j)=' '
      end do
      end   !** ends DelHeadStr

      subroutine DelHeadSpace(str)
c*** this subroutine is to delete all the heading spaces in the string str.
c*** input parameter :
c    str             : the string to be dealed with
c*** output parameter:
c    str             : with all the heading spaces deleted


      implicit real*8 (a-h,o-z)
      character *(*) str

      length=len(str)
      do iNonSp=1,length
            if ( str(iNonSp:iNonSp).ne.' ') then
c                  exit
                  goto 1001
            end if
      end do
1001  continue      

            !** now str(iNonSp:iNonSp) is the first non-space character
      if ( iNonSp.le.length ) then  !** there is non-space char
            imove=iNonSp-1
            do j=iNonSp,length
                  str(j-imove:j-imove)=str(j:j)
            end do
            do j=length-imove+1,length
                  str(j:j)=' '
            end do
      end if
      end   !** ends DelHeadSpace

      subroutine ToUpper(str)
c*** this subroutine is to  transform the string str to Uppercase .
c*** input parameter :
c    str             :
c*** output parameter:
c    str             : with all alphebit UPPERcase

      implicit real*8 (a-h,o-z)
      character *(*) str

      length=len(str)
      do i=1,length
            if ( ( ichar(str(i:i)).ge.ichar('a') )
     .          .and. (ichar(str(i:i)).le.ichar('z') )  ) then
                  str(i:i)=char( ichar( str(i:i) )
     .                     -ichar('a')+ichar('A') )
            end if
      end do !** i
      end   !** ends ToUpper


      subroutine UnSpace(str,istart,iend)
c*** this subroutine is to UNcount the heading and tailing spaces in the string
c    str
c*** input parameter  :
c    str
c*** output parameters:
c    istart           : the first non-space character index
c    iend             : the last non-space character index


      implicit real*8 (a-h,o-z)
      character *(*) str

      istart=1
      iend=len(str)
            !** uncount the heading spaces
      do while ( (str(istart:istart).eq.' ').and.(istart.le.iend) )
            istart=istart+1
      end do
            !** uncount the tailing spaces
      do while ( (str(iend:iend).eq.' ') .and. (iend.ge.istart) )
            iend=iend-1
      end do

      end   !** ends UnSpace
c
      subroutine GenInt2Str(int,str)
c
c*** this subroutine is to convert an integer to a string 
c    ( in General case where the FORTRAN compiler doesn't support interal file)
c*** note: the str should be at least 20-char long

      implicit none
            
      integer int
      character *(*) str

      character*40 tmpstr
      integer itmpf
      integer i,j,istart,iend

      itmpf=98      
      
      open(itmpf,status='scratch')
      write(itmpf,'(1x,i19)') int
      rewind(itmpf)   
      read(itmpf,'(a)') tmpstr
      close(itmpf)

      call UnSpace(tmpstr,istart,iend)
      str=''
      j=1
      do i=istart,iend
         str(j:j)=tmpstr(i:i)
         j=j+1
      enddo
      
      end	!ends GenInt2Str
      
      subroutine Int2Str(int,str)
c
c*** this subroutine is to convert an integer to a string 
c    ( the compilerr must support interal file)
c*** note: the str should be at least 20-char long

      implicit none
            
      integer int
      character *(*) str

      character*40 tmpstr
      integer i,j,istart,iend

      write(tmpstr,'(1x,i19)') int
      call UnSpace(tmpstr,istart,iend)
      str=''
      j=1
      do i=istart,iend
         str(j:j)=tmpstr(i:i)
         j=j+1
      enddo
      
      end	!ends Int2Str
c
c
      subroutine myflush(iunit)
c** A subroutine to write all pending output to a file attached
c   to logical unit "iunit" ( integer)

      implicit none

      character*65 filename
      integer iunit
      
      return
      
c      inquire(iunit,err=10,NAME=filename)
c      close(iunit)
c      open(iunit,file=filename,status='old',access='APPEND',err=20)
c      
c      return
c
c10    write(*,*) 'Cannot inquire the file in subroutine flush'
c      stop
c
c20    write(*,*) 'file handling error in subroutine flush.'
c      stop
      
      end	!ends myflush     
c
