      subroutine getun(un,dn,d2,shrknode,cnn2,ncnn2,hn,hm2)

      implicit none
      include "global.h"

      real* 8 un(ndf,nnc),dn(ndf,nnc),d2(ndf,nn_on2)
      real* 8 shrknode(maxconn,nnc)
      integer cnn2(maxconn,nnc),ncnn2(nnc)
      real* 8 hn(nnc),hm2(nn_on2)
      
      integer inn,inl,isd,node

      call grab_all2(d2,dn,ndf,hn,hm2)

      un(:,:) = 0

      do inn = 1,nnc
        do inl = 1,ncnn2(inn)
          node = cnn2(inl,inn)
          un(:,inn) = un(:,inn) + shrknode(inl,inn) * d2(:,node)
        end do
      end do

      return
      end

