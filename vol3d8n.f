      subroutine vol3d8n(vol,xn,xq)

      include "global.h"
      real* 8 vol, xn(nsdpad,nenpad)

      do iq = 1,nquad
      
c  shape function values
        sq(0,1,iq) = (1 - xq(1,iq))
     &       * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
        sq(0,2,iq) = (1 + xq(1,iq))
     &       * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
        sq(0,3,iq) = (1 + xq(1,iq))
     &       * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
        sq(0,4,iq) = (1 - xq(1,iq))
     &       * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
        sq(0,5,iq) = (1 - xq(1,iq))
     &       * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
        sq(0,6,iq) = (1 + xq(1,iq))
     &       * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
        sq(0,7,iq) = (1 + xq(1,iq))
     &       * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
        sq(0,8,iq) = (1 - xq(1,iq))
     &       * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
        
c  local first derivatives
        sq(1,1,iq) = - (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
        sq(1,2,iq) = + (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
        sq(1,3,iq) = + (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
        sq(1,4,iq) = - (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
        sq(1,5,iq) = - (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
        sq(1,6,iq) = + (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
        sq(1,7,iq) = + (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
        sq(1,8,iq) = - (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
        
        sq(2,1,iq) = - (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
        sq(2,2,iq) = - (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
        sq(2,3,iq) = + (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
        sq(2,4,iq) = + (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
        sq(2,5,iq) = - (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
        sq(2,6,iq) = - (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
        sq(2,7,iq) = + (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
        sq(2,8,iq) = + (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
        
        sq(3,1,iq) = - (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
        sq(3,2,iq) = - (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
        sq(3,3,iq) = - (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
        sq(3,4,iq) = - (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
        sq(3,5,iq) = + (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
        sq(3,6,iq) = + (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
        sq(3,7,iq) = + (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
        sq(3,8,iq) = + (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
        


      
