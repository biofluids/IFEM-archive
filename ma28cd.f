
      subroutine ma28cd(n, a, licn, icn, ikeep, rhs, w, mtype)
c
c this subroutine uses the factors from ma28a/ad or ma28b/bd to
c     solve a system of equations without iterative refinement.
c the parameters are ...
c n   integer  order of matrix  not altered by subroutine.
c a   real/double precision array  length licn.  the same array as
c     was used in the most recent call to ma28a/ad or ma28b/bd.
c licn  integer  length of arrays a and icn.  not altered by
c     subroutine.
c icn    integer array of length licn.  same array as output from
c     ma28a/ad.  unchanged by ma28c/cd.
c ikeep  integer array of length 5*n.  same array as output from
c     ma28a/ad.  unchanged by ma28c/cd.
c rhs    real/double precision array  length n.  on entry, it holds the
c     right hand side.  on exit, the solution vector.
c w      real/double precision array  length n. used as workspace by
c     ma30c/cd.
c mtype  integer  used to tell ma30c/cd to solve the direct equation
c     (mtype=1) or its transpose (mtype.ne.1).
c
      real* 8 a(licn), rhs(n), w(n), resid, mresid, eps, rmin
      integer idisp(2)
      integer icn(licn), ikeep(n,5)
      logical abort1, abort2
c common block variables.
c unless otherwise stated common block variables are as in ma28a/ad.
c     those variables referenced by ma28c/cd are mentioned below.
c resid  real/double precision  variable returns maximum residual of
c     equations where pivot was zero.
c mresid  real/double precision variable used by ma28c/cd to
c     communicate between ma28f/fd and ma30h/hd.
c idisp  integer array  length 2  the same as that used by ma28a/ad.
c     it is unchanged by ma28b/bd.
c
c further information on common block variables can be found in block
c     data or ma28a/ad.
      common /ma28fd/ eps, rmin, resid, irncp, icncp, minirn, minicn,
     * irank, abort1, abort2
      common /ma28gd/ idisp
      common /ma30hd/ mresid
c
c this call performs the solution of the set of equations.
      call ma30cd(n, icn, a, licn, ikeep, ikeep(1,4), ikeep(1,5),
     * idisp, ikeep(1,2), ikeep(1,3), rhs, w, mtype)
c transfer common block information.
      resid = mresid
      return
      end
      subroutine ma30cd(n, icn, a, licn, lenr, lenrl, lenoff, idisp, ip,
     * iq, x, w, mtype)
c ma30c/cd uses the factors produced by ma30a/ad or ma30b/bd to solve
c     ax=b or a transpose x=b when the matrix p1*a*q1 (paq) is block
c     lower triangular (including the case of only one diagonal
c     block).
c
c we now describe the argument list for ma30c/cd.
c n  is an integer variable set to the order of the matrix. it is not
c     altered by the subroutine.
c icn is an integer array of length licn. entries idisp(1) to
c     idisp(2) should be unchanged since the last call to ma30a/ad. if
c     the matrix has more than one diagonal block, then column indices
c     corresponding to non-zeros in sub-diagonal blocks of paq must
c     appear in positions 1 to idisp(1)-1. for the same row those
c     entries must be contiguous, with those in row i preceding those
c     in row i+1 (i=1,...,n-1) and no wasted space between rows.
c     entries may be in any order within each row. it is not altered
c     by ma30c/cd.
c a  is a real/double precision array of length licn.  entries
c     idisp(1) to idisp(2) should be unchanged since the last call to
c     ma30a/ad or ma30b/bd.  if the matrix has more than one diagonal
c     block, then the values of the non-zeros in sub-diagonal blocks
c     must be in positions 1 to idisp(1)-1 in the order given by icn.
c     it is not altered by ma30c/cd.
c licn  is an integer variable set to the size of arrays icn and a.
c     it is not altered by ma30c/cd.
c lenr,lenrl are integer arrays of length n which should be
c     unchanged since the last call to ma30a/ad. they are not altered
c     by ma30c/cd.
c lenoff  is an integer array of length n. if the matrix paq (or
c     p1*a*q1) has more than one diagonal block, then lenoff(i),
c     i=1,...,n should be set to the number of non-zeros in row i of
c     the matrix paq which are in sub-diagonal blocks.  if there is
c     only one diagonal block then lenoff(1) may be set to -1, in
c     which case the other entries of lenoff are never accessed. it is
c     not altered by ma30c/cd.
c idisp  is an integer array of length 2 which should be unchanged
c     since the last call to ma30a/ad. it is not altered by ma30c/cd.
c ip,iq are integer arrays of length n which should be unchanged
c     since the last call to ma30a/ad. they are not altered by
c     ma30c/cd.
c x is a real/double precision array of length n. it must be set by
c     the user to the values of the right hand side vector b for the
c     equations being solved.  on exit from ma30c/cd it will be equal
c     to the solution x required.
c w  is a real/double precision array of length n which is used as
c     workspace by ma30c/cd.
c mtype is an integer variable which must be set by the user. if
c     mtype=1, then the solution to the system ax=b is returned; any
c     other value for mtype will return the solution to the system a
c     transpose x=b. it is not altered by ma30c/cd.
c
      real* 8 a(licn), x(n), w(n), wii, wi, resid, zero
      logical neg, nobloc
      integer idisp(2)
      integer icn(licn), lenr(n), lenrl(n), lenoff(n), ip(n), iq(n)
c see block data for comments on variables in common.
      common /ma30hd/ resid
      data zero /0.0/
c
c the final value of resid is the maximum residual for an inconsistent
c     set of equations.
      resid = zero
c nobloc is .true. if subroutine block has been used previously and
c     is .false. otherwise.  the value .false. means that lenoff
c     will not be subsequently accessed.
      nobloc = lenoff(1).lt.0
      if (mtype.ne.1) go to 140
c
c we now solve   a * x = b.
c neg is used to indicate when the last row in a block has been
c     reached.  it is then set to true whereafter backsubstitution is
c     performed on the block.
      neg = .false.
c ip(n) is negated so that the last row of the last block can be
c     recognised.  it is reset to its positive value on exit.
      ip(n) = -ip(n)
c preorder vector ... w(i) = x(ip(i))
      do 10 ii=1,n
        i = ip(ii)
        i = iabs(i)
        w(ii) = x(i)
   10 continue
c lt holds the position of the first non-zero in the current row of the
c     off-diagonal blocks.
      lt = 1
c ifirst holds the index of the first row in the current block.
      ifirst = 1
c iblock holds the position of the first non-zero in the current row
c     of the lu decomposition of the diagonal blocks.
      iblock = idisp(1)
c if i is not the last row of a block, then a pass through this loop
c     adds the inner product of row i of the off-diagonal blocks and w
c     to w and performs forward elimination using row i of the lu
c     decomposition.   if i is the last row of a block then, after
c     performing these aforementioned operations, backsubstitution is
c     performed using the rows of the block.
      do 120 i=1,n
        wi = w(i)
        if (nobloc) go to 30
        if (lenoff(i).eq.0) go to 30
c operations using lower triangular blocks.
c ltend is the end of row i in the off-diagonal blocks.
        ltend = lt + lenoff(i) - 1
        do 20 jj=lt,ltend
          j = icn(jj)
          wi = wi - a(jj)*w(j)
   20   continue
c lt is set the beginning of the next off-diagonal row.
        lt = ltend + 1
c set neg to .true. if we are on the last row of the block.
   30   if (ip(i).lt.0) neg = .true.
        if (lenrl(i).eq.0) go to 50
c forward elimination phase.
c iend is the end of the l part of row i in the lu decomposition.
        iend = iblock + lenrl(i) - 1
        do 40 jj=iblock,iend
          j = icn(jj)
          wi = wi + a(jj)*w(j)
   40   continue
c iblock is adjusted to point to the start of the next row.
   50   iblock = iblock + lenr(i)
        w(i) = wi
        if (.not.neg) go to 120
c back substitution phase.
c j1 is position in a/icn after end of block beginning in row ifirst
c     and ending in row i.
        j1 = iblock
c are there any singularities in this block?  if not, continue with
c     the backsubstitution.
        ib = i
        if (iq(i).gt.0) go to 70
        do 60 iii=ifirst,i
          ib = i - iii + ifirst
          if (iq(ib).gt.0) go to 70
          j1 = j1 - lenr(ib)
          resid = dmax1(resid,dabs(w(ib)))
          w(ib) = zero
   60   continue
c entire block is singular.
        go to 110
c each pass through this loop performs the back-substitution
c     operations for a single row, starting at the end of the block and
c     working through it in reverse order.
   70   do 100 iii=ifirst,ib
          ii = ib - iii + ifirst
c j2 is end of row ii.
          j2 = j1 - 1
c j1 is beginning of row ii.
          j1 = j1 - lenr(ii)
c jpiv is the position of the pivot in row ii.
          jpiv = j1 + lenrl(ii)
          jpivp1 = jpiv + 1
c jump if row  ii of u has no non-zeros.
          if (j2.lt.jpivp1) go to 90
          wii = w(ii)
          do 80 jj=jpivp1,j2
            j = icn(jj)
            wii = wii - a(jj)*w(j)
   80     continue
          w(ii) = wii
   90     w(ii) = w(ii)/a(jpiv)
  100   continue
  110   ifirst = i + 1
        neg = .false.
  120 continue
c
c reorder solution vector ... x(i) = w(iqinverse(i))
      do 130 ii=1,n
        i = iq(ii)
        i = iabs(i)
        x(i) = w(ii)
  130 continue
      ip(n) = -ip(n)
      go to 320
c
c
c we now solve   atranspose * x = b.
c preorder vector ... w(i)=x(iq(i))
  140 do 150 ii=1,n
        i = iq(ii)
        i = iabs(i)
        w(ii) = x(i)
  150 continue
c lj1 points to the beginning the current row in the off-diagonal
c     blocks.
      lj1 = idisp(1)
c iblock is initialized to point to the beginning of the block after
c     the last one ]
      iblock = idisp(2) + 1
c ilast is the last row in the current block.
      ilast = n
c iblend points to the position after the last non-zero in the
c     current block.
      iblend = iblock
c each pass through this loop operates with one diagonal block and
c     the off-diagonal part of the matrix corresponding to the rows
c     of this block.  the blocks are taken in reverse order and the
c     number of times the loop is entered is min(n,no. blocks+1).
      do 290 numblk=1,n
        if (ilast.eq.0) go to 300
        iblock = iblock - lenr(ilast)
c this loop finds the index of the first row in the current block..
c     it is first and iblock is set to the position of the beginning
c     of this first row.
        do 160 k=1,n
          ii = ilast - k
          if (ii.eq.0) go to 170
          if (ip(ii).lt.0) go to 170
          iblock = iblock - lenr(ii)
  160   continue
  170   ifirst = ii + 1
c j1 points to the position of the beginning of row i (lt part) or pivot
        j1 = iblock
c forward elimination.
c each pass through this loop performs the operations for one row of the
c     block.  if the corresponding entry of w is zero then the
c     operations can be avoided.
        do 210 i=ifirst,ilast
          if (w(i).eq.zero) go to 200
c jump if row i singular.
          if (iq(i).lt.0) go to 220
c j2 first points to the pivot in row i and then is made to point to the
c     first non-zero in the u transpose part of the row.
          j2 = j1 + lenrl(i)
          wi = w(i)/a(j2)
          if (lenr(i)-lenrl(i).eq.1) go to 190
          j2 = j2 + 1
c j3 points to the end of row i.
          j3 = j1 + lenr(i) - 1
          do 180 jj=j2,j3
            j = icn(jj)
            w(j) = w(j) - a(jj)*wi
  180     continue
  190     w(i) = wi
  200     j1 = j1 + lenr(i)
  210   continue
        go to 240
c deals with rest of block which is singular.
  220   do 230 ii=i,ilast
          resid = dmax1(resid,dabs(w(ii)))
          w(ii) = zero
  230   continue
c back substitution.
c this loop does the back substitution on the rows of the block in
c     the reverse order doing it simultaneously on the l transpose part
c     of the diagonal blocks and the off-diagonal blocks.
  240   j1 = iblend
        do 280 iback=ifirst,ilast
          i = ilast - iback + ifirst
c j1 points to the beginning of row i.
          j1 = j1 - lenr(i)
          if (lenrl(i).eq.0) go to 260
c j2 points to the end of the l transpose part of row i.
          j2 = j1 + lenrl(i) - 1
          do 250 jj=j1,j2
            j = icn(jj)
            w(j) = w(j) + a(jj)*w(i)
  250     continue
  260     if (nobloc) go to 280
c operations using lower triangular blocks.
          if (lenoff(i).eq.0) go to 280
c lj2 points to the end of row i of the off-diagonal blocks.
          lj2 = lj1 - 1
c lj1 points to the beginning of row i of the off-diagonal blocks.
          lj1 = lj1 - lenoff(i)
          do 270 jj=lj1,lj2
            j = icn(jj)
            w(j) = w(j) - a(jj)*w(i)
  270     continue
  280   continue
        iblend = j1
        ilast = ifirst - 1
  290 continue
c reorder solution vector ... x(i)=w(ipinverse(i))
  300 do 310 ii=1,n
        i = ip(ii)
        i = iabs(i)
        x(i) = w(ii)
  310 continue
c
  320 return
      end
