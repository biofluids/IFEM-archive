      subroutine correct_w(b,bd,bdd,wx,wxd,wxdd,wy,wyd,wydd,
     &                     wxx,wxxd,wxxdd,wxy,wxyd,wxydd,
     &                     wyy,wyyd,wyydd,
     &                     cpt,cjp,dcjp,dwjp,nep,iInter)
c
      implicit double precision (a-h,o-z)
c
c      implicit none
c
      integer nep,iInter
      real*8 b(*),bd(2,*),bdd(3,*),cpt(2)
      real*8 wx(*),wxd(2,*),wxdd(3,*)
      real*8 wy(*),wyd(2,*),wydd(3,*)
      real*8 wxx(*),wxxd(2,*),wxxdd(3,*)
      real*8 wxy(*),wxyd(2,*),wxydd(3,*)
      real*8 wyy(*),wyyd(2,*),wyydd(3,*)
      real*8 cjp(2,nep),dcjp(2,nep),dwjp(nep)
c
      if (iInter .eq. 1) then
	 call correct_Law(b,bd,wx,wxd,wy,wyd,
     &                    cpt,cjp,dcjp,dwjp,nep)
      elseif(iInter .eq. 11) then
	 call correct_blw(b,bd,bdd,wx,wxd,wxdd,
     &                   wy,wyd,wydd,wxy,wxyd,wxydd,
     &                   cpt,cjp,dcjp,dwjp,nep)
      elseif(iInter .eq. 2)  then
	 call correct_qdw(b,bd,bdd,wx,wxd,wxdd,
     &                   wy,wyd,wydd,wxx,wxxd,wxxdd,
     &                   wxy,wxyd,wxydd,wyy,wyyd,wyydd,
     &                   cpt,cjp,dcjp,dwjp,nep)
      else
	 call correct_Lbw(b,bd,bdd,wx,wxd,wxdd,
     &                    wy,wyd,wydd,
     &                    cpt,cjp,dcjp,dwjp,nep)
      endif
c
c
      return
      end
c

