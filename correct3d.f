      subroutine correct3d(b,bd,cpt,cjp,dcjp,
     &                   dwjp,nep,iInter)
c
c......3-D correct function......
c
      implicit none
c
      real*8 b(*),bd(3,*),cpt(3)
      real*8 cjp(3,nep),dcjp(3,nep),dwjp(nep)
      integer nep,iInter
c
      if (iInter .eq. 1) then
	 call correct3dl(b,bd,cpt,cjp,dcjp,dwjp,nep)
      elseif(iInter .eq. 11)  then
	 call correct3dtl(b,bd,cpt,cjp,dcjp,dwjp,nep)
      else
	 print *, 'wrong iInter'
	 stop
      endif
c
c
      return
      end
c
c
