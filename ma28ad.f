
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias ma28ad ma28bd ma28cd
c###### calls   ma30    mc20    mc22    mc23    mc24
      subroutine ma28ad(n, nz, a, licn, irn, lirn, icn, u, ikeep, iw, w,
     * iflag)
c this subroutine performs the lu factorization of a.
c
c the parameters are as follows.....
c n     order of matrix  not altered by subroutine.
c nz    number of non-zeros in input matrix  not altered by subroutine.
c a is a  real array  length licn.  holds non-zeros of matrix on entry
c     and non-zeros of factors on exit.  reordered by mc20a/ad and
c     mc23a/ad and altered by ma30a/ad.
c licn  integer  length of arrays a and icn.  not altered by subroutine.
c irn   integer array of length lirn.  holds row indices on input.
c     used as workspace by ma30a/ad to hold column orientation of
c     matrix.
c lirn  integer  length of array irn. not altered by the subroutine.
c icn   integer array of length licn.  holds column indices on entry
c     and column indices of decomposed matrix on exit. reordered by
c     mc20a/ad and mc23a/ad and altered by ma30a/ad.
c u     real variable  set by user to control bias towards numeric or
c     sparsity pivoting.  u=1.0 gives partial pivoting while u=0. does
c     not check multipliers at all.  values of u greater than one are
c     treated as one while negative values are treated as zero.  not
c     altered by subroutine.
c ikeep  integer array of length 5*n  used as workspace by ma28a/ad
c     (see later comments).  it is not required to be set on entry
c     and, on exit, it contains information about the decomposition.
c     it should be preserved between this call and subsequent calls
c     to ma28b/bd or ma28c/cd.
c     ikeep(i,1),i=1,n  holds the total length of the part of row i
c     in the diagonal block.
c     row ikeep(i,2),i=1,n  of the input matrix is the ith row in
c     pivot order.
c     column ikeep(i,3),i=1,n  of the input matrix is the ith column
c     in pivot order.
c     ikeep(i,4),i=1,n  holds the length of the part of row i in
c     the l part of the l/u decomposition.
c     ikeep(i,5),i=1,n  holds the length of the part of row i in the
c     off-diagonal blocks.  if there is only one diagonal block,
c     ikeep(1,5) will be set to -1.
c iw    integer array of length 8*n.  if the option nsrch.le.n is
c     used, then the length of array iw can be reduced to 7*n.
c w real array  length n.  used by mc24a/ad both as workspace and to
c     return growth estimate in w(1).  the use of this array by ma28a/ad
c     is thus optional depending on common block logical variable grow.
c iflag  integer variable  used as error flag by routine.  a positive
c     or zero value on exit indicates success.  possible negative
c     values are -1 through -14.
c
      integer n, nz, licn, lirn, iflag
      integer irn(lirn), icn(licn), ikeep(n,5), iw(n,8)
      real* 8 a(licn), u, w(n)
c
c common and private variables.
c     common block ma28f/fd is used merely
c     to communicate with common block ma30f/fd  so that the user
c     need not declare this common block in his main program.
c the common block variables are as follows ...
c lp,mp  integer  default value 6 (line printer).  unit number
c     for error messages and duplicate element warning resp.
c nlp,mlp  integer  unit number for messages from ma30a/ad and
c     mc23a/ad resp.  set by ma28a/ad to value of lp.
c lblock  logical  default value true.  if true mc23a/ad is used
c     to first permute the matrix to block lower triangular form.
c grow    logical  default value true.  if true then an estimate
c     of the increase in size of matrix elements during l/u
c     decomposition is given by mc24a/ad.
c eps,rmin,resid  real/double precision variables not referenced
c     by ma28a/ad.
c irncp,icncp  integer  set to number of compresses on arrays irn and
c     icn/a respectively.
c minirn,minicn  integer  minimum length of arrays irn and icn/a
c     respectively, for success on future runs.
c irank  integer   estimated rank of matrix.
c mirncp,micncp,mirank,mirn,micn integer variables.  used to
c     communicate between ma30f/fd and ma28f/fd values of abovenamed
c     variables with somewhat similar names.
c abort1,abort2  logical variables with default value true.  if false
c     then decomposition will be performed even if the matrix is
c     structurally or numerically singular respectively.
c aborta,abortb  logical variables used to communicate values of
c     abort1 and abort2 to ma30a/ad.
c abort  logical  used to communicate value of abort1 to mc23a/ad.
c abort3  logical variable not referenced by ma28a/ad.
c idisp   integer array  length 2.  used to communicate information
c     on decomposition between this call to ma28a/ad and subsequent
c     calls to ma28b/bd and ma28c/cd.  on exit, idisp(1) and
c     idisp(2) indicate position in arrays a and icn of the
c     first and last elements in the l/u decomposition of the
c     diagonal blocks, respectively.
c numnz  integer  structural rank of matrix.
c num    integer  number of diagonal blocks.
c large  integer  size of largest diagonal block.
c
c see block data for further comments on common block variables.
c see code for comments on private variables.
c
      real* 8 tol, themax, big, dxmax, errmax, dres, cgce,
     * tol1, big1, upriv, rmin, eps, resid, zero
      integer idisp(2)
      logical grow, lblock, abort, abort1, abort2, abort3, aborta,
     * abortb, lbig, lbig1
      common /ma28ed/ lp, mp, lblock, grow
      common /ma28fd/ eps, rmin, resid, irncp, icncp, minirn, minicn,
     * irank, abort1, abort2
      common /ma28gd/ idisp
      common /ma28hd/ tol, themax, big, dxmax, errmax, dres, cgce,
     * ndrop, maxit, noiter, nsrch, istart, lbig
      common /ma30id/ tol1, big1, ndrop1, nsrch1, lbig1
      common /ma30ed/ nlp, aborta, abortb, abort3
      common /ma30fd/ mirncp, micncp, mirank, mirn, micn
      common /mc23bd/ mlp, numnz, num, large, abort
      common /lpivot/ lpiv(10),lnpiv(10),mapiv,manpiv,iavpiv,
     *                ianpiv,kountl
c
c some  initialization and transfer of information between
c     common blocks (see earlier comments).
      data zero /0.0/
      iflag = 0
      aborta = abort1
      abortb = abort2
      abort = abort1
      mlp = lp
      nlp = lp
      tol1 = tol
      lbig1 = lbig
      nsrch1 = nsrch
c upriv private copy of u is used in case it is outside
c     range  zero to one  and  is thus altered by ma30a/ad.
      upriv = u
c simple data check on input variables and array dimensions.
      if (n.gt.0) go to 10
      iflag = -8
      if (lp.ne.0) write (lp,99999) n
      go to 210
   10 if (nz.gt.0) go to 20
      iflag = -9
      if (lp.ne.0) write (lp,99998) nz
      go to 210
   20 if (licn.ge.nz) go to 30
      iflag = -10
      if (lp.ne.0) write (lp,99997) licn
      go to 210
   30 if (lirn.ge.nz) go to 40
      iflag = -11
      if (lp.ne.0) write (lp,99996) lirn
      go to 210
c
c data check to see if all indices lie between 1 and n.
   40 do 50 i=1,nz
        if (irn(i).gt.0 .and. irn(i).le.n .and. icn(i).gt.0 .and.
     *   icn(i).le.n) go to 50
        if (iflag.eq.0 .and. lp.ne.0) write (lp,99995)
        iflag = -12
        if (lp.ne.0) write (lp,99994) i, a(i), irn(i), icn(i)
   50 continue
      if (iflag.lt.0) go to 220
c
c sort matrix into row order.
      call mc20ad(n, nz, a, icn, iw, irn, 0)
c part of ikeep is used here as a work-array.  ikeep(i,2) is
c     the last row to have a non-zero in column i.  ikeep(i,3)
c     is the off-set of column i from the start of the row.
      do 60 i=1,n
        ikeep(i,2) = 0
        ikeep(i,1) = 0
   60 continue
c
c check for duplicate elements .. summing any such entries and
c     printing a warning message on unit mp.
c move is equal to the number of duplicate elements found.
      move = 0
c the loop also calculates the largest element in the matrix, themax.
      themax = zero
c j1 is position in arrays of first non-zero in row.
      j1 = iw(1,1)
      do 130 i=1,n
        iend = nz + 1
        if (i.ne.n) iend = iw(i+1,1)
        length = iend - j1
        if (length.eq.0) go to 130
        j2 = iend - 1
        newj1 = j1 - move
        do 120 jj=j1,j2
          j = icn(jj)
          themax = dmax1(themax,dabs(a(jj)))
          if (ikeep(j,2).eq.i) go to 110
c first time column has ocurred in current row.
          ikeep(j,2) = i
          ikeep(j,3) = jj - move - newj1
          if (move.eq.0) go to 120
c shift necessary because of  previous duplicate element.
          newpos = jj - move
          a(newpos) = a(jj)
          icn(newpos) = icn(jj)
          go to 120
c duplicate element.
  110     move = move + 1
          length = length - 1
          jay = ikeep(j,3) + newj1
          if (mp.ne.0) write (mp,99993) i, j, a(jj)
          a(jay) = a(jay) + a(jj)
          themax = dmax1(themax,dabs(a(jay)))
  120   continue
        ikeep(i,1) = length
        j1 = iend
  130 continue
c
c knum is actual number of non-zeros in matrix with any multiple
c     entries counted only once.
      knum = nz - move
      if (.not.lblock) go to 140
c
c perform block triangularisation.
      call mc23ad(n, icn, a, licn, ikeep, idisp, ikeep(1,2),
     *ikeep(1,3), ikeep(1,5), iw(1,3), iw)
      if (idisp(1).gt.0) go to 170
      iflag = -7
      if (idisp(1).eq.-1) iflag = -1
      if (lp.ne.0) write (lp,99992)
      go to 210
c
c block triangularization not requested.
c move structure to end of data arrays in preparation for
c     ma30a/ad.
c also set lenoff(1) to -1 and set permutation arrays.
  140 do 150 i=1,knum
        ii = knum - i + 1
        newpos = licn - i + 1
        icn(newpos) = icn(ii)
        a(newpos) = a(ii)
  150 continue
      idisp(1) = 1
      idisp(2) = licn - knum + 1
      do 160 i=1,n
        ikeep(i,2) = i
        ikeep(i,3) = i
  160 continue
      ikeep(1,5) = -1
  170 if (lbig) big1 = themax
      if (nsrch.le.n) go to 180
c
c perform l/u decomosition on diagonal blocks.
      call ma30ad(n, icn, a, licn, ikeep, ikeep(1,4), idisp,
     *ikeep(1,2), ikeep(1,3), irn, lirn, iw(1,2), iw(1,3), iw(1,4),
     *iw(1,5), iw(1,6), iw(1,7), iw(1,8), iw, upriv, iflag)
      go to 190
c this call if used if nsrch has been set less than or equal n.
c     in this case, two integer work arrays of length can be saved.
  180 call ma30ad(n, icn, a, licn, ikeep, ikeep(1,4), idisp,
     * ikeep(1,2), ikeep(1,3), irn, lirn, iw(1,2), iw(1,3), iw(1,4),
     * iw(1,5), iw, iw, iw(1,6), iw, upriv, iflag)
c
c transfer common block information.
  190 minirn = max0(mirn,nz)
      minicn = max0(micn,nz)
      irncp = mirncp
      icncp = micncp
      irank = mirank
      ndrop = ndrop1
      if (lbig) big = big1
      if (iflag.ge.0) go to 200
      if (lp.ne.0) write (lp,99991)
      go to 210
c
c reorder off-diagonal blocks according to pivot permutation.
  200 i1 = idisp(1) - 1
      if (i1.ne.0) call mc22ad(n, icn, a, i1, ikeep(1,5), ikeep(1,2),
     * ikeep(1,3), iw, irn)
      i1 = idisp(1)
      iend = licn - i1 + 1
c
c optionally calculate element growth estimate.
      if (grow) call mc24ad(n, icn, a(i1), iend, ikeep, ikeep(1,4), w)
c increment growth estimate by original maximum element.
      if (grow) w(1) = w(1) + themax
      if (grow .and. n.gt.1) w(2) = themax
c set flag if the only error is due to duplicate elements.
      if (iflag.ge.0 .and. move.ne.0) iflag = -14
      go to 220
  210 if (lp.ne.0) write (lp,99990)
  220 return
99999 format (36x, 17hn out of range = , i10)
99998 format (36x, 18hnz non positive = , i10)
99997 format (36x, 17hlicn too small = , i10)
99996 format (36x, 17hlirn too small = , i10)
99995 format (54h error return from ma28a/ad because indices found out ,
     * 8hof range)
99994 format (1x, i6, 22hth element with value , 1pd22.14, 9h is out o,
     * 21hf range with indices , i8, 2h ,, i8)
99993 format (31h duplicate element in position , i8, 2h ,, i8,
     * 12h with value , 1pd22.14)
99992 format (36x, 26herror return from mc23a/ad)
99991 format (36x, 26herror return from ma30a/ad)
99990 format (36h+error return from ma28a/ad because )
      end
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias ma30ad
      subroutine ma30ad(nn, icn, a, licn, lenr, lenrl, idisp, ip, iq,
     * irn, lirn, lenc, ifirst, lastr, nextr, lastc, nextc, iptr, ipc,
     * u, iflag)
c if  the user requires a more convenient data interface then the ma28
c     package should be used.  the ma28 subroutines call the ma30
c     subroutines after checking the user's input data and optionally
c     using mc23a/ad to permute the matrix to block triangular form.
c this package of subroutines (ma30a/ad, ma30b/bd, ma30c/cd and
c     ma30d/dd) performs operations pertinent to the solution of a
c     general sparse n by n system of linear equations (i.e. solve
c     ax=b). structually singular matrices are permitted including
c     those with row or columns consisting entirely of zeros (i.e.
c     including rectangular matrices).  it is assumed that the
c     non-zeros of the matrix a do not differ widely in size.  if
c     necessary a prior call of the scaling subroutine mc19a/ad may be
c     made.
c a discussion of the design of these subroutines is given by duff and
c     reid (acm trans math software 5 pp 18-35,1979 (css 48)) while
c     fuller details of the implementation are given in duff (harwell
c     report aere-r 8730,1977).  the additional pivoting option in
c     ma30a/ad and the use of drop tolerances (see common block
c     ma30i/id) were added to the package after joint work with reid,
c     schaumburg, wasniewski and zlatev (duff, reid, schaumburg,
c     wasniewski and zlatev, harwell report css 135, 1983).
c
c ma30a/ad performs the lu decomposition of the diagonal blocks of the
c     permutation paq of a sparse matrix a, where input permutations
c     p1 and q1 are used to define the diagonal blocks.  there may be
c     non-zeros in the off-diagonal blocks but they are unaffected by
c     ma30a/ad. p and p1 differ only within blocks as do q and q1. the
c     permutations p1 and q1 may be found by calling mc23a/ad or the
c     matrix may be treated as a single block by using p1=q1=i. the
c     matrix non-zeros should be held compactly by rows, although it
c     should be noted that the user can supply the matrix by columns
c     to get the lu decomposition of a transpose.
c
c the parameters are...
c this description should also be consulted for further information on
c     most of the parameters of ma30b/bd and ma30c/cd.
c
c n  is an integer variable which must be set by the user to the order
c     of the matrix.  it is not altered by ma30a/ad.
c icn is an integer array of length licn. positions idisp(2) to
c     licn must be set by the user to contain the column indices of
c     the non-zeros in the diagonal blocks of p1*a*q1. those belonging
c     to a single row must be contiguous but the ordering of column
c     indices with each row is unimportant. the non-zeros of row i
c     precede those of row i+1,i=1,...,n-1 and no wasted space is
c     allowed between the rows.  on output the column indices of the
c     lu decomposition of paq are held in positions idisp(1) to
c     idisp(2), the rows are in pivotal order, and the column indices
c     of the l part of each row are in pivotal order and precede those
c     of u. again there is no wasted space either within a row or
c     between the rows. icn(1) to icn(idisp(1)-1), are neither
c     required nor altered. if mc23a/ad been called, these will hold
c     information about the off-diagonal blocks.
c a is a real/double precision array of length licn whose entries
c     idisp(2) to licn must be set by the user to the  values of the
c     non-zero entries of the matrix in the order indicated by  icn.
c     on output a will hold the lu factors of the matrix where again
c     the position in the matrix is determined by the corresponding
c     values in icn. a(1) to a(idisp(1)-1) are neither required nor
c     altered.
c licn  is an integer variable which must be set by the user to the
c     length of arrays icn and a. it must be big enough for a and icn
c     to hold all the non-zeros of l and u and leave some "elbow
c     room".  it is possible to calculate a minimum value for licn by
c     a preliminary run of ma30a/ad. the adequacy of the elbow room
c     can be judged by the size of the common block variable icncp. it
c     is not altered by ma30a/ad.
c lenr  is an integer array of length n.  on input, lenr(i) should
c     equal the number of non-zeros in row i, i=1,...,n of the
c     diagonal blocks of p1*a*q1. on output, lenr(i) will equal the
c     total number of non-zeros in row i of l and row i of u.
c lenrl  is an integer array of length n. on output from ma30a/ad,
c     lenrl(i) will hold the number of non-zeros in row i of l.
c idisp  is an integer array of length 2. the user should set idisp(1)
c     to be the first available position in a/icn for the lu
c     decomposition while idisp(2) is set to the position in a/icn of
c     the first non-zero in the diagonal blocks of p1*a*q1. on output,
c     idisp(1) will be unaltered while idisp(2) will be set to the
c     position in a/icn of the last non-zero of the lu decomposition.
c ip  is an integer array of length n which holds a permutation of
c     the integers 1 to n.  on input to ma30a/ad, the absolute value of
c     ip(i) must be set to the row of a which is row i of p1*a*q1. a
c     negative value for ip(i) indicates that row i is at the end of a
c     diagonal block.  on output from ma30a/ad, ip(i) indicates the row
c     of a which is the i th row in paq. ip(i) will still be negative
c     for the last row of each block (except the last).
c iq is an integer array of length n which again holds a
c     permutation of the integers 1 to n.  on input to ma30a/ad, iq(j)
c     must be set to the column of a which is column j of p1*a*q1. on
c     output from ma30a/ad, the absolute value of iq(j) indicates the
c     column of a which is the j th in paq.  for rows, i say, in which
c     structural or numerical singularity is detected iq(i) is
c     negated.
c irn  is an integer array of length lirn used as workspace by
c     ma30a/ad.
c lirn  is an integer variable. it should be greater than the
c     largest number of non-zeros in a diagonal block of p1*a*q1 but
c     need not be as large as licn. it is the length of array irn and
c     should be large enough to hold the active part of any block,
c     plus some "elbow room", the  a posteriori  adequacy of which can
c     be estimated by examining the size of common block variable
c     irncp.
c lenc,ifirst,lastr,nextr,lastc,nextc are all integer arrays of
c     length n which are used as workspace by ma30a/ad.  if nsrch is
c     set to a value less than or equal to n, then arrays lastc and
c     nextc are not referenced by ma30a/ad and so can be dummied in
c     the call to ma30a/ad.
c iptr,ipc are integer arrays of length n which are used as workspace
c     by ma30a/ad.
c u  is a real/double precision variable which should be set by the
c     user to a value between 0. and 1.0. if less than zero it is
c     reset to zero and if its value is 1.0 or greater it is reset to
c     0.9999 (0.999999999 in d version).  it determines the balance
c     between pivoting for sparsity and for stability, values near
c     zero emphasizing sparsity and values near one emphasizing
c     stability. we recommend u=0.1 as a posible first trial value.
c     the stability can be judged by a later call to mc24a/ad or by
c     setting lbig to .true.
c iflag  is an integer variable. it will have a non-negative value if
c     ma30a/ad is successful. negative values indicate error
c     conditions while positive values indicate that the matrix has
c     been successfully decomposed but is singular. for each non-zero
c     value, an appropriate message is output on unit lp.  possible
c     non-zero values for iflag are ...
c
c -1  the matrix is structually singular with rank given by irank in
c     common block ma30f/fd.
c +1  if, however, the user wants the lu decomposition of a
c     structurally singular matrix and sets the common block variable
c     abort1 to .false., then, in the event of singularity and a
c     successful decomposition, iflag is returned with the value +1
c     and no message is output.
c -2  the matrix is numerically singular (it may also be structually
c     singular) with estimated rank given by irank in common block
c     ma30f/fd.
c +2  the  user can choose to continue the decomposition even when a
c     zero pivot is encountered by setting common block variable
c     abort2 to .false.  if a singularity is encountered, iflag will
c     then return with a value of +2, and no message is output if the
c     decomposition has been completed successfully.
c -3  lirn has not been large enough to continue with the
c     decomposition.  if the stage was zero then common block variable
c     minirn gives the length sufficient to start the decomposition on
c     this block.  for a successful decomposition on this block the
c     user should make lirn slightly (say about n/2) greater than this
c     value.
c -4  licn not large enough to continue with the decomposition.
c -5  the decomposition has been completed but some of the lu factors
c     have been discarded to create enough room in a/icn to continue
c     the decomposition. the variable minicn in common block ma30f/fd
c     then gives the size that licn should be to enable the
c     factorization to be successful.  if the user sets common block
c     variable abort3 to .true., then the subroutine will exit
c     immediately instead of destroying any factors and continuing.
c -6  both licn and lirn are too small. termination has been caused by
c     lack of space in irn (see error iflag= -3), but already some of
c     the lu factors in a/icn have been lost (see error iflag= -5).
c     minicn gives the minimum amount of space required in a/icn for
c     decomposition up to this point.
c
      real* 8 a(licn), u, au, umax, amax, zero, pivrat, pivr,
     * tol, big, anew, aanew, scale
      integer iptr(nn), pivot, pivend, dispc, oldpiv, oldend, pivrow,
     * rowi, ipc(nn), idisp(2), colupd
      integer icn(licn), lenr(nn), lenrl(nn), ip(nn), iq(nn),
     * lenc(nn), irn(lirn), ifirst(nn), lastr(nn), nextr(nn),
     * lastc(nn), nextc(nn)
      logical abort1, abort2, abort3, lbig
c for comments of common block variables see block data subprogram.
      common /ma30ed/ lp, abort1, abort2, abort3
      common /ma30fd/ irncp, icncp, irank, minirn, minicn
      common /ma30id/ tol, big, ndrop, nsrch, lbig
      common /lpivot/ lpiv(10),lnpiv(10),mapiv,manpiv,iavpiv,
     *                ianpiv,kountl
c
      data umax/.999999999/
      data zero /0.0/
      msrch = nsrch
      ndrop = 0
      do 1272 kk=1,10
        lnpiv(kk)=0
        lpiv(kk)=0
 1272 continue
      mapiv = 0
      manpiv = 0
      iavpiv = 0
      ianpiv = 0
      kountl = 0
      minirn = 0
      minicn = idisp(1) - 1
      morei = 0
      irank = nn
      irncp = 0
      icncp = 0
      iflag = 0
c reset u if necessary.
      u = dmin1(u,umax)
c ibeg is the position of the next pivot row after elimination step
c     using it.
      u = dmax1(u,zero)
      ibeg = idisp(1)
c iactiv is the position of the first entry in the active part of a/icn.
      iactiv = idisp(2)
c nzrow is current number of non-zeros in active and unprocessed part
c     of row file icn.
      nzrow = licn - iactiv + 1
      minicn = nzrow + minicn
c
c count the number of diagonal blocks and set up pointers to the
c     beginnings of the rows.
c num is the number of diagonal blocks.
      num = 1
      iptr(1) = iactiv
      if (nn.eq.1) go to 20
      nnm1 = nn - 1
      do 10 i=1,nnm1
        if (ip(i).lt.0) num = num + 1
        iptr(i+1) = iptr(i) + lenr(i)
   10 continue
c ilast is the last row in the previous block.
   20 ilast = 0
c
c ***********************************************
c ****    lu decomposition of block nblock   ****
c ***********************************************
c
c each pass through this loop performs lu decomposition on one
c     of the diagonal blocks.
      do 1000 nblock=1,num
        istart = ilast + 1
        do 30 irows=istart,nn
          if (ip(irows).lt.0) go to 40
   30   continue
        irows = nn
   40   ilast = irows
c n is the number of rows in the current block.
c istart is the index of the first row in the current block.
c ilast is the index of the last row in the current block.
c iactiv is the position of the first entry in the block.
c itop is the position of the last entry in the block.
        n = ilast - istart + 1
        if (n.ne.1) go to 90
c
c code for dealing with 1x1 block.
        lenrl(ilast) = 0
        ising = istart
        if (lenr(ilast).ne.0) go to 50
c block is structurally singular.
        irank = irank - 1
        ising = -ising
        if (iflag.ne.2 .and. iflag.ne.-5) iflag = 1
        if (.not.abort1) go to 80
        idisp(2) = iactiv
        iflag = -1
        if (lp.ne.0) write (lp,99999)
c     return
        go to 1120
   50   scale = dabs(a(iactiv))
        if (scale.eq.zero) go to 60
        if (lbig) big = dmax1(big,scale)
        go to 70
   60   ising = -ising
        irank = irank - 1
        iptr(ilast) = 0
        if (iflag.ne.-5) iflag = 2
        if (.not.abort2) go to 70
        idisp(2) = iactiv
        iflag = -2
        if (lp.ne.0) write (lp,99998)
        go to 1120
   70   a(ibeg) = a(iactiv)
        icn(ibeg) = icn(iactiv)
        iactiv = iactiv + 1
        iptr(istart) = 0
        ibeg = ibeg + 1
        nzrow = nzrow - 1
   80   lastr(istart) = istart
        ipc(istart) = -ising
        go to 1000
c
c non-trivial block.
   90   itop = licn
        if (ilast.ne.nn) itop = iptr(ilast+1) - 1
c
c set up column oriented storage.
        do 100 i=istart,ilast
          lenrl(i) = 0
          lenc(i) = 0
  100   continue
        if (itop-iactiv.lt.lirn) go to 110
        minirn = itop - iactiv + 1
        pivot = istart - 1
        go to 1100
c
c calculate column counts.
  110   do 120 ii=iactiv,itop
          i = icn(ii)
          lenc(i) = lenc(i) + 1
  120   continue
c set up column pointers so that ipc(j) points to position after end
c     of column j in column file.
        ipc(ilast) = lirn + 1
        j1 = istart + 1
        do 130 jj=j1,ilast
          j = ilast - jj + j1 - 1
          ipc(j) = ipc(j+1) - lenc(j+1)
  130   continue
        do 150 indrow=istart,ilast
          j1 = iptr(indrow)
          j2 = j1 + lenr(indrow) - 1
          if (j1.gt.j2) go to 150
          do 140 jj=j1,j2
            j = icn(jj)
            ipos = ipc(j) - 1
            irn(ipos) = indrow
            ipc(j) = ipos
  140     continue
  150   continue
c dispc is the lowest indexed active location in the column file.
        dispc = ipc(istart)
        nzcol = lirn - dispc + 1
        minirn = max0(nzcol,minirn)
        nzmin = 1
c
c initialize array ifirst.  ifirst(i) = +/- k indicates that row/col
c     k has i non-zeros.  if ifirst(i) = 0, there is no row or column
c     with i non zeros.
        do 160 i=1,n
          ifirst(i) = 0
  160   continue
c
c compute ordering of row and column counts.
c first run through columns (from column n to column 1).
        do 180 jj=istart,ilast
          j = ilast - jj + istart
          nz = lenc(j)
          if (nz.ne.0) go to 170
          ipc(j) = 0
          go to 180
  170     if (nsrch.le.nn) go to 180
          isw = ifirst(nz)
          ifirst(nz) = -j
          lastc(j) = 0
          nextc(j) = -isw
          isw1 = iabs(isw)
          if (isw.ne.0) lastc(isw1) = j
  180   continue
c now run through rows (again from n to 1).
        do 210 ii=istart,ilast
          i = ilast - ii + istart
          nz = lenr(i)
          if (nz.ne.0) go to 190
          iptr(i) = 0
          lastr(i) = 0
          go to 210
  190     isw = ifirst(nz)
          ifirst(nz) = i
          if (isw.gt.0) go to 200
          nextr(i) = 0
          lastr(i) = isw
          go to 210
  200     nextr(i) = isw
          lastr(i) = lastr(isw)
          lastr(isw) = i
  210   continue
c
c
c **********************************************
c ****    start of main elimination loop    ****
c **********************************************
        do 980 pivot=istart,ilast
c
c first find the pivot using markowitz criterion with stability
c     control.
c jcost is the markowitz cost of the best pivot so far,.. this
c     pivot is in row ipiv and column jpiv.
          nz2 = nzmin
          jcost = n*n
c
c examine rows/columns in order of ascending count.
          do 340 l=1,2
            pivrat = zero
            isrch = 1
            ll = l
c a pass with l equal to 2 is only performed in the case of singularity.
            do 330 nz=nz2,n
              if (jcost.le.(nz-1)**2) go to 420
              ijfir = ifirst(nz)
              if (ijfir) 230, 220, 240
  220         if (ll.eq.1) nzmin = nz + 1
              go to 330
  230         ll = 2
              ijfir = -ijfir
              go to 290
  240         ll = 2
c scan rows with nz non-zeros.
              do 270 idummy=1,n
                if (jcost.le.(nz-1)**2) go to 420
                if (isrch.gt.msrch) go to 420
                if (ijfir.eq.0) go to 280
c row ijfir is now examined.
                i = ijfir
                ijfir = nextr(i)
c first calculate multiplier threshold level.
                amax = zero
                j1 = iptr(i) + lenrl(i)
                j2 = iptr(i) + lenr(i) - 1
                do 250 jj=j1,j2
                  amax = dmax1(amax,dabs(a(jj)))
  250           continue
                au = amax*u
                isrch = isrch + 1
c scan row for possible pivots
                do 260 jj=j1,j2
                  if (dabs(a(jj)).le.au .and. l.eq.1) go to 260
                  j = icn(jj)
                  kcost = (nz-1)*(lenc(j)-1)
                  if (kcost.gt.jcost) go to 260
                  pivr = zero
                  if (amax.ne.zero) pivr = dabs(a(jj))/amax
                  if (kcost.eq.jcost .and. (pivr.le.pivrat .or.
     *             nsrch.gt.nn+1)) go to 260
c best pivot so far is found.
                  jcost = kcost
                  ijpos = jj
                  ipiv = i
                  jpiv = j
                  if (msrch.gt.nn+1 .and. jcost.le.(nz-1)**2) go to 420
                  pivrat = pivr
  260           continue
  270         continue
c
c columns with nz non-zeros now examined.
  280         ijfir = ifirst(nz)
              ijfir = -lastr(ijfir)
  290         if (jcost.le.nz*(nz-1)) go to 420
              if (msrch.le.nn) go to 330
              do 320 idummy=1,n
                if (ijfir.eq.0) go to 330
                j = ijfir
                ijfir = nextc(ijfir)
                i1 = ipc(j)
                i2 = i1 + nz - 1
c scan column j.
                do 310 ii=i1,i2
                  i = irn(ii)
                  kcost = (nz-1)*(lenr(i)-lenrl(i)-1)
                  if (kcost.ge.jcost) go to 310
c pivot has best markowitz count so far ... now check its
c     suitability on numeric grounds by examining the other non-zeros
c     in its row.
                  j1 = iptr(i) + lenrl(i)
                  j2 = iptr(i) + lenr(i) - 1
c we need a stability check on singleton columns because of possible
c     problems with underdetermined systems.
                  amax = zero
                  do 300 jj=j1,j2
                    amax = dmax1(amax,dabs(a(jj)))
                    if (icn(jj).eq.j) jpos = jj
  300             continue
                  if (dabs(a(jpos)).le.amax*u .and. l.eq.1) go to 310
                  jcost = kcost
                  ipiv = i
                  jpiv = j
                  ijpos = jpos
                  if (amax.ne.zero) pivrat = dabs(a(jpos))/amax
                  if (jcost.le.nz*(nz-1)) go to 420
  310           continue
c
  320         continue
c
  330       continue
c in the event of singularity, we must make sure all rows and columns
c are tested.
            msrch = n
c
c matrix is numerically or structurally singular  ... which it is will
c     be diagnosed later.
            irank = irank - 1
  340     continue
c assign rest of rows and columns to ordering array.
c matrix is structurally singular.
          if (iflag.ne.2 .and. iflag.ne.-5) iflag = 1
          irank = irank - ilast + pivot + 1
          if (.not.abort1) go to 350
          idisp(2) = iactiv
          iflag = -1
          if (lp.ne.0) write (lp,99999)
          go to 1120
  350     k = pivot - 1
          do 390 i=istart,ilast
            if (lastr(i).ne.0) go to 390
            k = k + 1
            lastr(i) = k
            if (lenrl(i).eq.0) go to 380
            minicn = max0(minicn,nzrow+ibeg-1+morei+lenrl(i))
            if (iactiv-ibeg.ge.lenrl(i)) go to 360
            call ma30dd(a, icn, iptr(istart), n, iactiv, itop, .true.)
c check now to see if ma30d/dd has created enough available space.
            if (iactiv-ibeg.ge.lenrl(i)) go to 360
c create more space by destroying previously created lu factors.
            morei = morei + ibeg - idisp(1)
            ibeg = idisp(1)
            if (lp.ne.0) write (lp,99997)
            iflag = -5
            if (abort3) go to 1090
  360       j1 = iptr(i)
            j2 = j1 + lenrl(i) - 1
            iptr(i) = 0
            do 370 jj=j1,j2
              a(ibeg) = a(jj)
              icn(ibeg) = icn(jj)
              icn(jj) = 0
              ibeg = ibeg + 1
  370       continue
            nzrow = nzrow - lenrl(i)
  380       if (k.eq.ilast) go to 400
  390     continue
  400     k = pivot - 1
          do 410 i=istart,ilast
            if (ipc(i).ne.0) go to 410
            k = k + 1
            ipc(i) = k
            if (k.eq.ilast) go to 990
  410     continue
c
c the pivot has now been found in position (ipiv,jpiv) in location
c     ijpos in row file.
c update column and row ordering arrays to correspond with removal
c     of the active part of the matrix.
  420     ising = pivot
          if (a(ijpos).ne.zero) go to 430
c numerical singularity is recorded here.
          ising = -ising
          if (iflag.ne.-5) iflag = 2
          if (.not.abort2) go to 430
          idisp(2) = iactiv
          iflag = -2
          if (lp.ne.0) write (lp,99998)
          go to 1120
  430     oldpiv = iptr(ipiv) + lenrl(ipiv)
          oldend = iptr(ipiv) + lenr(ipiv) - 1
c changes to column ordering.
          if (nsrch.le.nn) go to 460
          colupd = nn + 1
            lenpp = oldend-oldpiv+1
            if (lenpp.lt.4) lpiv(1) = lpiv(1) + 1
            if (lenpp.ge.4 .and. lenpp.le.6) lpiv(2) = lpiv(2) + 1
            if (lenpp.ge.7 .and. lenpp.le.10) lpiv(3) = lpiv(3) + 1
            if (lenpp.ge.11 .and. lenpp.le.15) lpiv(4) = lpiv(4) + 1
            if (lenpp.ge.16 .and. lenpp.le.20) lpiv(5) = lpiv(5) + 1
            if (lenpp.ge.21 .and. lenpp.le.30) lpiv(6) = lpiv(6) + 1
            if (lenpp.ge.31 .and. lenpp.le.50) lpiv(7) = lpiv(7) + 1
            if (lenpp.ge.51 .and. lenpp.le.70) lpiv(8) = lpiv(8) + 1
            if (lenpp.ge.71 .and. lenpp.le.100) lpiv(9) = lpiv(9) + 1
            if (lenpp.ge.101) lpiv(10) = lpiv(10) + 1
            mapiv = max0(mapiv,lenpp)
            iavpiv = iavpiv + lenpp
          do 450 jj=oldpiv,oldend
            j = icn(jj)
            lc = lastc(j)
            nc = nextc(j)
            nextc(j) = -colupd
            if (jj.ne.ijpos) colupd = j
            if (nc.ne.0) lastc(nc) = lc
            if (lc.eq.0) go to 440
            nextc(lc) = nc
            go to 450
  440       nz = lenc(j)
            isw = ifirst(nz)
            if (isw.gt.0) lastr(isw) = -nc
            if (isw.lt.0) ifirst(nz) = -nc
  450     continue
c changes to row ordering.
  460     i1 = ipc(jpiv)
          i2 = i1 + lenc(jpiv) - 1
          do 480 ii=i1,i2
            i = irn(ii)
            lr = lastr(i)
            nr = nextr(i)
            if (nr.ne.0) lastr(nr) = lr
            if (lr.le.0) go to 470
            nextr(lr) = nr
            go to 480
  470       nz = lenr(i) - lenrl(i)
            if (nr.ne.0) ifirst(nz) = nr
            if (nr.eq.0) ifirst(nz) = lr
  480     continue
c
c move pivot to position lenrl+1 in pivot row and move pivot row
c     to the beginning of the available storage.
c the l part and the pivot in the old copy of the pivot row is
c     nullified while, in the strictly upper triangular part, the
c     column indices, j say, are overwritten by the corresponding
c     entry of iq (iq(j)) and iq(j) is set to the negative of the
c     displacement of the column index from the pivot entry.
          if (oldpiv.eq.ijpos) go to 490
          au = a(oldpiv)
          a(oldpiv) = a(ijpos)
          a(ijpos) = au
          icn(ijpos) = icn(oldpiv)
          icn(oldpiv) = jpiv
c check to see if there is space immediately available in a/icn to
c     hold new copy of pivot row.
  490     minicn = max0(minicn,nzrow+ibeg-1+morei+lenr(ipiv))
          if (iactiv-ibeg.ge.lenr(ipiv)) go to 500
          call ma30dd(a, icn, iptr(istart), n, iactiv, itop, .true.)
          oldpiv = iptr(ipiv) + lenrl(ipiv)
          oldend = iptr(ipiv) + lenr(ipiv) - 1
c check now to see if ma30d/dd has created enough available space.
          if (iactiv-ibeg.ge.lenr(ipiv)) go to 500
c create more space by destroying previously created lu factors.
          morei = morei + ibeg - idisp(1)
          ibeg = idisp(1)
          if (lp.ne.0) write (lp,99997)
          iflag = -5
          if (abort3) go to 1090
          if (iactiv-ibeg.ge.lenr(ipiv)) go to 500
c there is still not enough room in a/icn.
          iflag = -4
          go to 1090
c copy pivot row and set up iq array.
  500     ijpos = 0
          j1 = iptr(ipiv)
c
          do 530 jj=j1,oldend
            a(ibeg) = a(jj)
            icn(ibeg) = icn(jj)
            if (ijpos.ne.0) go to 510
            if (icn(jj).eq.jpiv) ijpos = ibeg
            icn(jj) = 0
            go to 520
  510       k = ibeg - ijpos
            j = icn(jj)
            icn(jj) = iq(j)
            iq(j) = -k
  520       ibeg = ibeg + 1
  530     continue
c
          ijp1 = ijpos + 1
          pivend = ibeg - 1
          lenpiv = pivend - ijpos
          nzrow = nzrow - lenrl(ipiv) - 1
          iptr(ipiv) = oldpiv + 1
          if (lenpiv.eq.0) iptr(ipiv) = 0
c
c remove pivot row (including pivot) from column oriented file.
          do 560 jj=ijpos,pivend
            j = icn(jj)
            i1 = ipc(j)
            lenc(j) = lenc(j) - 1
c i2 is last position in new column.
            i2 = ipc(j) + lenc(j) - 1
            if (i2.lt.i1) go to 550
            do 540 ii=i1,i2
              if (irn(ii).ne.ipiv) go to 540
              irn(ii) = irn(i2+1)
              go to 550
  540       continue
  550       irn(i2+1) = 0
  560     continue
          nzcol = nzcol - lenpiv - 1
c
c go down the pivot column and for each row with a non-zero add
c     the appropriate multiple of the pivot row to it.
c we loop on the number of non-zeros in the pivot column since
c     ma30d/dd may change its actual position.
c
          nzpc = lenc(jpiv)
          if (nzpc.eq.0) go to 900
          do 840 iii=1,nzpc
            ii = ipc(jpiv) + iii - 1
            i = irn(ii)
c search row i for non-zero to be eliminated, calculate multiplier,
c     and place it in position lenrl+1 in its row.
c  idrop is the number of non-zero entries dropped from row    i
c        because these fall beneath tolerance level.
c
            idrop = 0
            j1 = iptr(i) + lenrl(i)
            iend = iptr(i) + lenr(i) - 1
            do 570 jj=j1,iend
              if (icn(jj).ne.jpiv) go to 570
c if pivot is zero, rest of column is and so multiplier is zero.
              au = zero
              if (a(ijpos).ne.zero) au = -a(jj)/a(ijpos)
              if (lbig) big = dmax1(big,dabs(au))
              a(jj) = a(j1)
              a(j1) = au
              icn(jj) = icn(j1)
              icn(j1) = jpiv
              lenrl(i) = lenrl(i) + 1
              go to 580
  570       continue
c jump if pivot row is a singleton.
  580       if (lenpiv.eq.0) go to 840
c now perform necessary operations on rest of non-pivot row i.
            rowi = j1 + 1
            iop = 0
c jump if all the pivot row causes fill-in.
            if (rowi.gt.iend) go to 650
c perform operations on current non-zeros in row i.
c innermost loop.
            lenpp = iend-rowi+1
            if (lenpp.lt.4) lnpiv(1) = lnpiv(1) + 1
            if (lenpp.ge.4 .and. lenpp.le.6) lnpiv(2) = lnpiv(2) + 1
            if (lenpp.ge.7 .and. lenpp.le.10) lnpiv(3) = lnpiv(3) + 1
            if (lenpp.ge.11 .and. lenpp.le.15) lnpiv(4) = lnpiv(4) + 1
            if (lenpp.ge.16 .and. lenpp.le.20) lnpiv(5) = lnpiv(5) + 1
            if (lenpp.ge.21 .and. lenpp.le.30) lnpiv(6) = lnpiv(6) + 1
            if (lenpp.ge.31 .and. lenpp.le.50) lnpiv(7) = lnpiv(7) + 1
            if (lenpp.ge.51 .and. lenpp.le.70) lnpiv(8) = lnpiv(8) + 1
            if (lenpp.ge.71 .and. lenpp.le.100) lnpiv(9) = lnpiv(9) + 1
            if (lenpp.ge.101) lnpiv(10) = lnpiv(10) + 1
            manpiv = max0(manpiv,lenpp)
            ianpiv = ianpiv + lenpp
            kountl = kountl + 1
            do 590 jj=rowi,iend
              j = icn(jj)
              if (iq(j).gt.0) go to 590
              iop = iop + 1
              pivrow = ijpos - iq(j)
              a(jj) = a(jj) + au*a(pivrow)
              if (lbig) big = dmax1(dabs(a(jj)),big)
              icn(pivrow) = -icn(pivrow)
              if (dabs(a(jj)).lt.tol) idrop = idrop + 1
  590       continue
c
c  jump if no non-zeros in non-pivot row have been removed
c       because these are beneath the drop-tolerance  tol.
c
            if (idrop.eq.0) go to 650
c
c  run through non-pivot row compressing row so that only
c      non-zeros greater than   tol   are stored.  all non-zeros
c      less than   tol   are also removed from the column structure.
c
            jnew = rowi
            do 630 jj=rowi,iend
              if (dabs(a(jj)).lt.tol) go to 600
              a(jnew) = a(jj)
              icn(jnew) = icn(jj)
              jnew = jnew + 1
              go to 630
c
c  remove non-zero entry from column structure.
c
  600         j = icn(jj)
              i1 = ipc(j)
              i2 = i1 + lenc(j) - 1
              do 610 ii=i1,i2
                if (irn(ii).eq.i) go to 620
  610         continue
  620         irn(ii) = irn(i2)
              irn(i2) = 0
              lenc(j) = lenc(j) - 1
              if (nsrch.le.nn) go to 630
c remove column from column chain and place in update chain.
              if (nextc(j).lt.0) go to 630
c jump if column already in update chain.
              lc = lastc(j)
              nc = nextc(j)
              nextc(j) = -colupd
              colupd = j
              if (nc.ne.0) lastc(nc) = lc
              if (lc.eq.0) go to 622
              nextc(lc) = nc
              go to 630
  622         nz = lenc(j) + 1
              isw = ifirst(nz)
              if (isw.gt.0) lastr(isw) = -nc
              if (isw.lt.0) ifirst(nz) = -nc
  630       continue
            do 640 jj=jnew,iend
              icn(jj) = 0
  640       continue
c the value of idrop might be different from that calculated earlier
c     because, we may now have dropped some non-zeros which were not
c     modified by the pivot row.
            idrop = iend + 1 - jnew
            iend = jnew - 1
            lenr(i) = lenr(i) - idrop
            nzrow = nzrow - idrop
            nzcol = nzcol - idrop
            ndrop = ndrop + idrop
  650       ifill = lenpiv - iop
c jump is if there is no fill-in.
            if (ifill.eq.0) go to 750
c now for the fill-in.
            minicn = max0(minicn,morei+ibeg-1+nzrow+ifill+lenr(i))
c see if there is room for fill-in.
c get maximum space for row i in situ.
            do 660 jdiff=1,ifill
              jnpos = iend + jdiff
              if (jnpos.gt.licn) go to 670
              if (icn(jnpos).ne.0) go to 670
  660       continue
c there is room for all the fill-in after the end of the row so it
c     can be left in situ.
c next available space for fill-in.
            iend = iend + 1
            go to 750
c jmore spaces for fill-in are required in front of row.
  670       jmore = ifill - jdiff + 1
            i1 = iptr(i)
c we now look in front of the row to see if there is space for
c     the rest of the fill-in.
            do 680 jdiff=1,jmore
              jnpos = i1 - jdiff
              if (jnpos.lt.iactiv) go to 690
              if (icn(jnpos).ne.0) go to 700
  680       continue
  690       jnpos = i1 - jmore
            go to 710
c whole row must be moved to the beginning of available storage.
  700       jnpos = iactiv - lenr(i) - ifill
c jump if there is space immediately available for the shifted row.
  710       if (jnpos.ge.ibeg) go to 730
            call ma30dd(a, icn, iptr(istart), n, iactiv, itop, .true.)
            i1 = iptr(i)
            iend = i1 + lenr(i) - 1
            jnpos = iactiv - lenr(i) - ifill
            if (jnpos.ge.ibeg) go to 730
c no space available so try to create some by throwing away previous
c     lu decomposition.
            morei = morei + ibeg - idisp(1) - lenpiv - 1
            if (lp.ne.0) write (lp,99997)
            iflag = -5
            if (abort3) go to 1090
c keep record of current pivot row.
            ibeg = idisp(1)
            icn(ibeg) = jpiv
            a(ibeg) = a(ijpos)
            ijpos = ibeg
            do 720 jj=ijp1,pivend
              ibeg = ibeg + 1
              a(ibeg) = a(jj)
              icn(ibeg) = icn(jj)
  720       continue
            ijp1 = ijpos + 1
            pivend = ibeg
            ibeg = ibeg + 1
            if (jnpos.ge.ibeg) go to 730
c this still does not give enough room.
            iflag = -4
            go to 1090
  730       iactiv = min0(iactiv,jnpos)
c move non-pivot row i.
            iptr(i) = jnpos
            do 740 jj=i1,iend
              a(jnpos) = a(jj)
              icn(jnpos) = icn(jj)
              jnpos = jnpos + 1
              icn(jj) = 0
  740       continue
c first new available space.
            iend = jnpos
  750       nzrow = nzrow + ifill
c innermost fill-in loop which also resets icn.
            idrop = 0
            do 830 jj=ijp1,pivend
              j = icn(jj)
              if (j.lt.0) go to 820
              anew = au*a(jj)
              aanew = dabs(anew)
              if (aanew.ge.tol) go to 760
              idrop = idrop + 1
              ndrop = ndrop + 1
              nzrow = nzrow - 1
              minicn = minicn - 1
              ifill = ifill - 1
              go to 830
  760         if (lbig) big = dmax1(aanew,big)
              a(iend) = anew
              icn(iend) = j
              iend = iend + 1
c
c put new entry in column file.
              minirn = max0(minirn,nzcol+lenc(j)+1)
              jend = ipc(j) + lenc(j)
              jroom = nzpc - iii + 1 + lenc(j)
              if (jend.gt.lirn) go to 770
              if (irn(jend).eq.0) go to 810
  770         if (jroom.lt.dispc) go to 780
c compress column file to obtain space for new copy of column.
              call ma30dd(a, irn, ipc(istart), n, dispc, lirn, .false.)
              if (jroom.lt.dispc) go to 780
              jroom = dispc - 1
              if (jroom.ge.lenc(j)+1) go to 780
c column file is not large enough.
              go to 1100
c copy column to beginning of file.
  780         jbeg = ipc(j)
              jend = ipc(j) + lenc(j) - 1
              jzero = dispc - 1
              dispc = dispc - jroom
              idispc = dispc
              do 790 ii=jbeg,jend
                irn(idispc) = irn(ii)
                irn(ii) = 0
                idispc = idispc + 1
  790         continue
              ipc(j) = dispc
              jend = idispc
              do 800 ii=jend,jzero
                irn(ii) = 0
  800         continue
  810         irn(jend) = i
              nzcol = nzcol + 1
              lenc(j) = lenc(j) + 1
c end of adjustment to column file.
              go to 830
c
  820         icn(jj) = -j
  830       continue
            if (idrop.eq.0) go to 834
            do 832 kdrop=1,idrop
            icn(iend) = 0
            iend = iend + 1
  832       continue
  834       lenr(i) = lenr(i) + ifill
c end of scan of pivot column.
  840     continue
c
c
c remove pivot column from column oriented storage and update row
c     ordering arrays.
          i1 = ipc(jpiv)
          i2 = ipc(jpiv) + lenc(jpiv) - 1
          nzcol = nzcol - lenc(jpiv)
          do 890 ii=i1,i2
            i = irn(ii)
            irn(ii) = 0
            nz = lenr(i) - lenrl(i)
            if (nz.ne.0) go to 850
            lastr(i) = 0
            go to 890
  850       ifir = ifirst(nz)
            ifirst(nz) = i
            if (ifir) 860, 880, 870
  860       lastr(i) = ifir
            nextr(i) = 0
            go to 890
  870       lastr(i) = lastr(ifir)
            nextr(i) = ifir
            lastr(ifir) = i
            go to 890
  880       lastr(i) = 0
            nextr(i) = 0
            nzmin = min0(nzmin,nz)
  890     continue
c restore iq and nullify u part of old pivot row.
c    record the column permutation in lastc(jpiv) and the row
c    permutation in lastr(ipiv).
  900     ipc(jpiv) = -ising
          lastr(ipiv) = pivot
          if (lenpiv.eq.0) go to 980
          nzrow = nzrow - lenpiv
          jval = ijp1
          jzer = iptr(ipiv)
          iptr(ipiv) = 0
          do 910 jcount=1,lenpiv
            j = icn(jval)
            iq(j) = icn(jzer)
            icn(jzer) = 0
            jval = jval + 1
            jzer = jzer + 1
  910     continue
c adjust column ordering arrays.
          if (nsrch.gt.nn) go to 920
          do 916 jj=ijp1,pivend
            j = icn(jj)
            nz = lenc(j)
            if (nz.ne.0) go to 914
            ipc(j) = 0
            go to 916
  914       nzmin = min0(nzmin,nz)
  916     continue
          go to 980
  920     jj = colupd
          do 970 jdummy=1,nn
            j = jj
            if (j.eq.nn+1) go to 980
            jj = -nextc(j)
            nz = lenc(j)
            if (nz.ne.0) go to 924
            ipc(j) = 0
            go to 970
  924       ifir = ifirst(nz)
            lastc(j) = 0
            if (ifir) 930, 940, 950
  930       ifirst(nz) = -j
            ifir = -ifir
            lastc(ifir) = j
            nextc(j) = ifir
            go to 970
  940       ifirst(nz) = -j
            nextc(j) = 0
            go to 960
  950       lc = -lastr(ifir)
            lastr(ifir) = -j
            nextc(j) = lc
            if (lc.ne.0) lastc(lc) = j
  960       nzmin = min0(nzmin,nz)
  970     continue
  980   continue
c ********************************************
c ****    end of main elimination loop    ****
c ********************************************
c
c reset iactiv to point to the beginning of the next block.
  990   if (ilast.ne.nn) iactiv = iptr(ilast+1)
 1000 continue
c
c ********************************************
c ****    end of deomposition of block    ****
c ********************************************
c
c record singularity (if any) in iq array.
      if (irank.eq.nn) go to 1020
      do 1010 i=1,nn
        if (ipc(i).lt.0) go to 1010
        ising = ipc(i)
        iq(ising) = -iq(ising)
        ipc(i) = -ising
 1010 continue
c
c run through lu decomposition changing column indices to that of new
c     order and permuting lenr and lenrl arrays according to pivot
c     permutations.
 1020 istart = idisp(1)
      iend = ibeg - 1
      if (iend.lt.istart) go to 1040
      do 1030 jj=istart,iend
        jold = icn(jj)
        icn(jj) = -ipc(jold)
 1030 continue
 1040 do 1050 ii=1,nn
        i = lastr(ii)
        nextr(i) = lenr(ii)
        iptr(i) = lenrl(ii)
 1050 continue
      do 1060 i=1,nn
        lenrl(i) = iptr(i)
        lenr(i) = nextr(i)
 1060 continue
c
c update permutation arrays ip and iq.
      do 1070 ii=1,nn
        i = lastr(ii)
        j = -ipc(ii)
        nextr(i) = iabs(ip(ii)+0)
        iptr(j) = iabs(iq(ii)+0)
 1070 continue
      do 1080 i=1,nn
        if (ip(i).lt.0) nextr(i) = -nextr(i)
        ip(i) = nextr(i)
        if (iq(i).lt.0) iptr(i) = -iptr(i)
        iq(i) = iptr(i)
 1080 continue
      ip(nn) = iabs(ip(nn)+0)
      idisp(2) = iend
      go to 1120
c
c   ***    error returns    ***
 1090 idisp(2) = iactiv
      if (lp.eq.0) go to 1120
      write (lp,99996)
      go to 1110
 1100 if (iflag.eq.-5) iflag = -6
      if (iflag.ne.-6) iflag = -3
      idisp(2) = iactiv
      if (lp.eq.0) go to 1120
      if (iflag.eq.-3) write (lp,99995)
      if (iflag.eq.-6) write (lp,99994)
 1110 pivot = pivot - istart + 1
      write (lp,99993) pivot, nblock, istart, ilast
      if (pivot.eq.0) write (lp,99992) minirn
c
c
 1120 return
99999 format (54h error return from ma30a/ad because matrix is structur,
     * 13hally singular)
99998 format (54h error return from ma30a/ad because matrix is numerica,
     * 12hlly singular)
99997 format (48h lu decomposition destroyed to create more space)
99996 format (54h error return from ma30a/ad because licn not big enoug,
     * 1hh)
99995 format (54h error return from ma30a/ad because lirn not big enoug,
     * 1hh)
99994 format (51h error return from ma30a/ad lirn and licn too small)
99993 format (10h at stage , i5, 10h in block , i5, 16h with first row ,
     * i5, 14h and last row , i5)
99992 format (34h to continue set lirn to at least , i8)
      end
      subroutine ma30dd(a, icn, iptr, n, iactiv, itop, reals)
c this subroutine performs garbage collection operations on the
c     arrays a, icn and irn.
c iactiv is the first position in arrays a/icn from which the compress
c     starts.  on exit, iactiv equals the position of the first entry
c     in the compressed part of a/icn
c
      real* 8 a(itop)
      logical reals
      integer iptr(n)
      integer icn(itop)
c see block data for comments on variables in common.
      common /ma30fd/ irncp, icncp, irank, minirn, minicn
c
      if (reals) icncp = icncp + 1
      if (.not.reals) irncp = irncp + 1
c set the first non-zero entry in each row to the negative of the
c     row/col number and hold this row/col index in the row/col
c     pointer.  this is so that the beginning of each row/col can
c     be recognized in the subsequent scan.
      do 10 j=1,n
        k = iptr(j)
        if (k.lt.iactiv) go to 10
        iptr(j) = icn(k)
        icn(k) = -j
   10 continue
      kn = itop + 1
      kl = itop - iactiv + 1
c go through arrays in reverse order compressing to the back so
c     that there are no zeros held in positions iactiv to itop in icn.
c     reset first entry of each row/col and pointer array iptr.
      do 30 k=1,kl
        jpos = itop - k + 1
        if (icn(jpos).eq.0) go to 30
        kn = kn - 1
        if (reals) a(kn) = a(jpos)
        if (icn(jpos).ge.0) go to 20
c first non-zero of row/col has been located
        j = -icn(jpos)
        icn(jpos) = iptr(j)
        iptr(j) = kn
   20   icn(kn) = icn(jpos)
   30 continue
      iactiv = kn
      return
      end
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc13d
      subroutine mc13d(n,icn,licn,ip,lenr,ior,ib,num,iw)
      integer ip(n)
      integer icn(licn),lenr(n),ior(n),ib(n),iw(n,3)
      call mc13e(n,icn,licn,ip,lenr,ior,ib,num,iw(1,1),iw(1,2),iw(1,3))
      return
      end
      subroutine mc13e(n,icn,licn,ip,lenr,arp,ib,num,lowl,numb,prev)
      integer stp,dummy
      integer ip(n)
c
c arp(i) is one less than the number of unsearched edges leaving
c     node i.  at the end of the algorithm it is set to a
c     permutation which puts the matrix in block lower
c     triangular form.
c ib(i) is the position in the ordering of the start of the ith
c     block.  ib(n+1-i) holds the node number of the ith node
c     on the stack.
c lowl(i) is the smallest stack position of any node to which a path
c     from node i has been found.  it is set to n+1 when node i
c     is removed from the stack.
c numb(i) is the position of node i in the stack if it is on
c     it, is the permuted order of node i for those nodes
c     whose final position has been found and is otherwise zero.
c prev(i) is the node at the end of the path when node i was
c     placed on the stack.
      integer icn(licn),lenr(n),arp(n),ib(n),lowl(n),numb(n),
     1prev(n)
c
c
c   icnt is the number of nodes whose positions in final ordering have
c     been found.
      icnt=0
c num is the number of blocks that have been found.
      num=0
      nnm1=n+n-1
c
c initialization of arrays.
      do 20 j=1,n
      numb(j)=0
      arp(j)=lenr(j)-1
   20 continue
c
c
      do 120 isn=1,n
c look for a starting node
      if (numb(isn).ne.0) go to 120
      iv=isn
c ist is the number of nodes on the stack ... it is the stack pointer.
      ist=1
c put node iv at beginning of stack.
      lowl(iv)=1
      numb(iv)=1
      ib(n)=iv
c
c the body of this loop puts a new node on the stack or backtracks.
      do 110 dummy=1,nnm1
      i1=arp(iv)
c have all edges leaving node iv been searched.
      if (i1.lt.0) go to 60
      i2=ip(iv)+lenr(iv)-1
      i1=i2-i1
c
c look at edges leaving node iv until one enters a new node or
c     all edges are exhausted.
      do 50 ii=i1,i2
      iw=icn(ii)
c has node iw been on stack already.
      if (numb(iw).eq.0) go to 100
c update value of lowl(iv) if necessary.
  50  lowl(iv)=min0(lowl(iv),lowl(iw))
c
c there are no more edges leaving node iv.
      arp(iv)=-1
c is node iv the root of a block.
   60 if (lowl(iv).lt.numb(iv)) go to 90
c
c order nodes in a block.
      num=num+1
      ist1=n+1-ist
      lcnt=icnt+1
c peel block off the top of the stack starting at the top and
c     working down to the root of the block.
      do 70 stp=ist1,n
      iw=ib(stp)
      lowl(iw)=n+1
      icnt=icnt+1
      numb(iw)=icnt
      if (iw.eq.iv) go to 80
   70 continue
   80 ist=n-stp
      ib(num)=lcnt
c are there any nodes left on the stack.
      if (ist.ne.0) go to 90
c have all the nodes been ordered.
      if (icnt.lt.n) go to 120
      go to 130
c
c backtrack to previous node on path.
   90 iw=iv
      iv=prev(iv)
c update value of lowl(iv) if necessary.
      lowl(iv)=min0(lowl(iv),lowl(iw))
      go to 110
c
c put new node on the stack.
 100  arp(iv)=i2-ii-1
      prev(iw)=iv
      iv=iw
      ist=ist+1
      lowl(iv)=ist
      numb(iv)=ist
      k=n+1-ist
      ib(k)=iv
  110 continue
c
  120 continue
c
c
c put permutation in the required form.
  130 do 140 i=1,n
      ii=numb(i)
 140  arp(ii)=i
      return
      end
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc20ad mc20bd
      subroutine mc20ad(nc,maxa,a,inum,jptr,jnum,jdisp)
c
      integer   inum(maxa),jnum(maxa)
      real* 8 a(maxa),ace,acep
      dimension jptr(nc)
c
c     ******************************************************************
c
      null=-jdisp
c**      clear jptr
      do 60 j=1,nc
   60 jptr(j)=0
c**      count the number of elements in each column.
      do 120 k=1,maxa
      j=jnum(k)+jdisp
      jptr(j)=jptr(j)+1
  120 continue
c**      set the jptr array
      k=1
      do 150 j=1,nc
      kr=k+jptr(j)
      jptr(j)=k
  150 k=kr
c
c**      reorder the elements into column order.  the algorithm is an
c        in-place sort and is of order maxa.
      do 230 i=1,maxa
c        establish the current entry.
      jce=jnum(i)+jdisp
      if(jce.eq.0) go to 230
      ace=a(i)
      ice=inum(i)
c        clear the location vacated.
      jnum(i)=null
c        chain from current entry to store items.
      do 200 j=1,maxa
c        current entry not in correct position.  determine correct
c        position to store entry.
      loc=jptr(jce)
      jptr(jce)=jptr(jce)+1
c        save contents of that location.
      acep=a(loc)
      icep=inum(loc)
      jcep=jnum(loc)
c        store current entry.
      a(loc)=ace
      inum(loc)=ice
      jnum(loc)=null
c        check if next current entry needs to be processed.
      if(jcep.eq.null) go to 230
c        it does.  copy into current entry.
      ace=acep
      ice=icep
  200 jce=jcep+jdisp
c
  230 continue
c
c**      reset jptr vector.
      ja=1
      do 250 j=1,nc
      jb=jptr(j)
      jptr(j)=ja
  250 ja=jb
      return
      end
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc21a
      subroutine mc21a(n,icn,licn,ip,lenr,iperm,numnz,iw)
      integer ip(n)
      integer icn(licn),lenr(n),iperm(n),iw(n,4)
      call mc21b(n,icn,licn,ip,lenr,iperm,numnz,iw(1,1),iw(1,2),iw(1,3),
     1iw(1,4))
      return
      end
      subroutine mc21b(n,icn,licn,ip,lenr,iperm,numnz,pr,arp,cv,out)
      integer ip(n)
c   pr(i) is the previous row to i in the depth first search.
c it is used as a work array in the sorting algorithm.
c   elements (iperm(i),i) i=1, ... n  are non-zero at the end of the
c algorithm unless n assignments have not been made.  in which case
c (iperm(i),i) will be zero for n-numnz entries.
c   cv(i) is the most recent row extension at which column i
c was visited.
c   arp(i) is one less than the number of non-zeros in row i
c which have not been scanned when looking for a cheap assignment.
c   out(i) is one less than the number of non-zeros in row i
c which have not been scanned during one pass through the main loop.
      integer icn(licn),lenr(n),iperm(n),pr(n),cv(n),
     1arp(n),out(n)
c
c   initialization of arrays.
      do 10 i=1,n
      arp(i)=lenr(i)-1
      cv(i)=0
   10 iperm(i)=0
      numnz=0
c
c
c   main loop.
c   each pass round this loop either results in a new assignment
c or gives a row with no assignment.
      do 130 jord=1,n
      j=jord
      pr(j)=-1
      do 100 k=1,jord
c look for a cheap assignment
      in1=arp(j)
      if (in1.lt.0) go to 60
      in2=ip(j)+lenr(j)-1
      in1=in2-in1
      do 50 ii=in1,in2
      i=icn(ii)
      if (iperm(i).eq.0) go to 110
   50 continue
c   no cheap assignment in row.
      arp(j)=-1
c   begin looking for assignment chain starting with row j.
   60 out(j)=lenr(j)-1
c inner loop.  extends chain by one or backtracks.
      do 90 kk=1,jord
      in1=out(j)
      if (in1.lt.0) go to 80
      in2=ip(j)+lenr(j)-1
      in1=in2-in1
c forward scan.
      do 70 ii=in1,in2
      i=icn(ii)
      if (cv(i).eq.jord) go to 70
c   column i has not yet been accessed during this pass.
      j1=j
      j=iperm(i)
      cv(i)=jord
      pr(j)=j1
      out(j1)=in2-ii-1
      go to 100
   70 continue
c
c   backtracking step.
   80 j=pr(j)
      if (j.eq.-1) go to 130
   90 continue
c
  100 continue
c
c   new assignment is made.
  110 iperm(i)=j
      arp(j)=in2-ii-1
      numnz=numnz+1
      do 120 k=1,jord
      j=pr(j)
      if (j.eq.-1) go to 130
      ii=ip(j)+lenr(j)-out(j)-2
      i=icn(ii)
      iperm(i)=j
  120 continue
c
  130 continue
c
c   if matrix is structurally singular, we now complete the
c permutation iperm.
      if (numnz.eq.n) return
      do 140 i=1,n
  140 arp(i)=0
      k=0
      do 160 i=1,n
      if (iperm(i).ne.0) go to 150
      k=k+1
      out(k)=i
      go to 160
  150 j=iperm(i)
      arp(j)=i
  160 continue
      k=0
      do 170 i=1,n
      if (arp(i).ne.0) go to 170
      k=k+1
      ioutk=out(k)
      iperm(ioutk)=i
  170 continue
      return
      end
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc22ad
      subroutine mc22ad(n,icn,a,nz,lenrow,ip,iq,iw,iw1)
      real* 8 a(nz),aval
      integer iw(n,2)
      integer   icn(nz),lenrow(n),ip(n),iq(n),iw1(nz)
      if (nz.le.0) go to 1000
      if (n.le.0) go to 1000
c set start of row i in iw(i,1) and lenrow(i) in iw(i,2)
      iw(1,1)=1
      iw(1,2)=lenrow(1)
      do 10 i=2,n
      iw(i,1)=iw(i-1,1)+lenrow(i-1)
 10   iw(i,2)=lenrow(i)
c permute lenrow according to ip.  set off-sets for new position
c     of row iold in iw(iold,1) and put old row indices in iw1 in
c     positions corresponding to the new position of this row in a/icn.
      jj=1
      do 20 i=1,n
      iold=ip(i)
      iold=iabs(iold)
      length=iw(iold,2)
      lenrow(i)=length
      if (length.eq.0) go to 20
      iw(iold,1)=iw(iold,1)-jj
      j2=jj+length-1
      do 15 j=jj,j2
 15   iw1(j)=iold
      jj=j2+1
 20   continue
c set inverse permutation to iq in iw(.,2).
      do 30 i=1,n
      iold=iq(i)
      iold=iabs(iold)
 30   iw(iold,2)=i
c permute a and icn in place, changing to new column numbers.
c
c ***   main loop   ***
c each pass through this loop places a closed chain of column indices
c     in their new (and final) positions ... this is recorded by
c     setting the iw1 entry to zero so that any which are subsequently
c     encountered during this major scan can be bypassed.
      do 200 i=1,nz
      iold=iw1(i)
      if (iold.eq.0) go to 200
      ipos=i
      jval=icn(i)
c if row iold is in same positions after permutation go to 150.
      if (iw(iold,1).eq.0) go to 150
      aval=a(i)
c **  chain loop  **
c each pass through this loop places one (permuted) column index
c     in its final position  .. viz. ipos.
      do 100 ichain=1,nz
c newpos is the original position in a/icn of the element to be placed
c in position ipos.  it is also the position of the next element in
c     the chain.
      newpos=ipos+iw(iold,1)
c is chain complete ?
      if (newpos.eq.i) go to 130
      a(ipos)=a(newpos)
      jnum=icn(newpos)
      icn(ipos)=iw(jnum,2)
      ipos=newpos
      iold=iw1(ipos)
      iw1(ipos)=0
c **  end of chain loop  **
 100  continue
 130  a(ipos)=aval
 150  icn(ipos)=iw(jval,2)
c ***   end of main loop   ***
 200  continue
c
 1000 return
      end
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc23ad
c###### calls   mc13    mc21
      subroutine mc23ad(n,icn,a,licn,lenr,idisp,ip,iq,lenoff,iw,iw1)
      real* 8 a(licn)
      integer idisp(2),iw1(n,2)
      logical abort
      integer   icn(licn),lenr(n),ip(n),iq(n),lenoff(n),iw(n,5)
      common /mc23bd/ lp,numnz,num,large,abort
c input ... n,icn .. a,icn,lenr ....
c
c set up pointers iw(.,1) to the beginning of the rows and set lenoff
c     equal to lenr.
      iw1(1,1)=1
      lenoff(1)=lenr(1)
      if (n.eq.1) go to 20
      do 10 i=2,n
      lenoff(i)=lenr(i)
   10 iw1(i,1)=iw1(i-1,1)+lenr(i-1)
c idisp(1) points to the first position in a/icn after the
c     off-diagonal blocks and untreated rows.
   20 idisp(1)=iw1(n,1)+lenr(n)
c
c find row permutation ip to make diagonal zero-free.
      call mc21a(n,icn,licn,iw1,lenr,ip,numnz,iw)
c
c possible error return for structurally singular matrices.
      if (numnz.ne.n.and.abort) go to 170
c
c iw1(.,2) and lenr are permutations of iw1(.,1) and lenr/lenoff
c     suitable for entry
c     to mc13d since matrix with these row pointer and length arrays
c     has maximum number of non-zeros on the diagonal.
      do 30 ii=1,n
      i=ip(ii)
      iw1(ii,2)=iw1(i,1)
   30 lenr(ii)=lenoff(i)
c
c find symmetric permutation iq to block lower triangular form.
      call mc13d(n,icn,licn,iw1(1,2),lenr,iq,iw(1,4),num,iw)
c
      if (num.ne.1) go to 60
c
c action taken if matrix is irreducible.
c whole matrix is just moved to the end of the storage.
      do 40 i=1,n
      lenr(i)=lenoff(i)
      ip(i)=i
   40 iq(i)=i
      lenoff(1)=-1
c idisp(1) is the first position after the last element in the
c     off-diagonal blocks and untreated rows.
      nz=idisp(1)-1
      idisp(1)=1
c idisp(2) is the position in a/icn of the first element in the
c     diagonal blocks.
      idisp(2)=licn-nz+1
      large=n
      if (nz.eq.licn) go to 230
      do 50 k=1,nz
      j=nz-k+1
      jj=licn-k+1
      a(jj)=a(j)
   50 icn(jj)=icn(j)
c 230 = return
      go to 230
c
c data structure reordered.
c
c form composite row permutation ... ip(i) = ip(iq(i)).
   60 do 70 ii=1,n
      i=iq(ii)
   70 iw(ii,1)=ip(i)
      do 80 i=1,n
   80 ip(i)=iw(i,1)
c
c run through blocks in reverse order separating diagonal blocks
c     which are moved to the end of the storage.  elements in
c     off-diagonal blocks are left in place unless a compress is
c     necessary.
c
c ibeg indicates the lowest value of j for which icn(j) has been
c     set to zero when element in position j was moved to the
c     diagonal block part of storage.
      ibeg=licn+1
c iend is the position of the first element of those treated rows
c     which are in diagonal blocks.
      iend=licn+1
c large is the dimension of the largest block encountered so far.
      large=0
c
c num is the number of diagonal blocks.
      do 150 k=1,num
      iblock=num-k+1
c i1 is first row (in permuted form) of block iblock.
c i2 is last row (in permuted form) of block iblock.
      i1=iw(iblock,4)
      i2=n
      if (k.ne.1) i2=iw(iblock+1,4)-1
      large=max0(large,i2-i1+1)
c go through the rows of block iblock in the reverse order.
      do 140 ii=i1,i2
      inew=i2-ii+i1
c we now deal with row inew in permuted form (row iold in original
c     matrix).
      iold=ip(inew)
c if there is space to move up diagonal block portion of row go to 110
      if (iend-idisp(1).ge.lenoff(iold)) go to 110
c
c in-line compress.
c moves separated off-diagonal elements and untreated rows to
c     front of storage.
      jnpos=ibeg
      ilend=idisp(1)-1
      if (ilend.lt.ibeg) go to 190
      do 90 j=ibeg,ilend
      if (icn(j).eq.0) go to 90
      icn(jnpos)=icn(j)
      a(jnpos)=a(j)
      jnpos=jnpos+1
   90 continue
      idisp(1)=jnpos
      if (iend-jnpos.lt.lenoff(iold)) go to 190
      ibeg=licn+1
c reset pointers to the beginning of the rows.
      do 100 i=2,n
  100 iw1(i,1)=iw1(i-1,1)+lenoff(i-1)
c
c row iold is now split into diag. and off-diag. parts.
  110 irowb=iw1(iold,1)
      leni=0
      irowe=irowb+lenoff(iold)-1
c backward scan of whole of row iold (in original matrix).
      if (irowe.lt.irowb) go to 130
      do 120 jj=irowb,irowe
      j=irowe-jj+irowb
      jold=icn(j)
c iw(.,2) holds the inverse permutation to iq.
c     ..... it was set to this in mc13d.
      jnew=iw(jold,2)
c if (jnew.lt.i1) then ....
c element is in off-diagonal block and so is left in situ.
      if (jnew.lt.i1) go to 120
c element is in diagonal block and is moved to the end of the storage.
      iend=iend-1
      a(iend)=a(j)
      icn(iend)=jnew
      ibeg=min0(ibeg,j)
      icn(j)=0
      leni=leni+1
  120 continue
c
      lenoff(iold)=lenoff(iold)-leni
  130 lenr(inew)=leni
  140 continue
c
      ip(i2)=-ip(i2)
  150 continue
c resets ip(n) to positive value.
      ip(n)=-ip(n)
c idisp(2) is position of first element in diagonal blocks.
      idisp(2)=iend
c
c this compress is used to move all off-diagonal elements to the
c     front of the storage.
      if (ibeg.gt.licn) go to 230
      jnpos=ibeg
      ilend=idisp(1)-1
      do 160 j=ibeg,ilend
      if (icn(j).eq.0) go to 160
      icn(jnpos)=icn(j)
      a(jnpos)=a(j)
      jnpos=jnpos+1
  160 continue
c idisp(1) is first position after last element of off-diagonal blocks.
      idisp(1)=jnpos
      go to 230
c
c
c error return
  170 if (lp.ne.0) write(lp,180) numnz
  180 format(33x,41h matrix is structurally singular, rank = ,i6)
      idisp(1)=-1
      go to 210
  190 if (lp.ne.0) write(lp,200) n
  200 format(33x,33h licn not big enough increase by ,i6)
      idisp(1)=-2
  210 if (lp.ne.0) write(lp,220)
  220 format(33h+error return from mc23ad because)
c
  230 return
      end
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc24ad
      subroutine mc24ad(n,icn,a,licn,lenr,lenrl,w)
      real* 8 a(licn),w(n),amaxl,wrowl,amaxu,zero
      integer   icn(licn),lenr(n),lenrl(n)
      data zero/0.0/
      amaxl=zero
      do 10 i=1,n
 10   w(i)=zero
      j0=1
      do 100 i=1,n
      if (lenr(i).eq.0) go to 100
      j2=j0+lenr(i)-1
      if (lenrl(i).eq.0) go to 50
c calculation of 1-norm of l.
      j1=j0+lenrl(i)-1
      wrowl=zero
      do 30 jj=j0,j1
 30   wrowl=wrowl+dabs(a(jj))
c amaxl is the maximum norm of columns of l so far found.
      amaxl=dmax1(amaxl,wrowl)
      j0=j1+1
c calculation of norms of columns of u (max-norms).
 50   j0=j0+1
      if (j0.gt.j2) go to 90
      do 80 jj=j0,j2
      j=icn(jj)
 80   w(j)=dmax1(dabs(a(jj)),w(j))
 90   j0=j2+1
 100  continue
c amaxu is set to maximum max-norm of columns of u.
      amaxu=zero
      do 200 i=1,n
 200  amaxu=dmax1(amaxu,w(i))
c grofac is max u max-norm times max l 1-norm.
      w(1)=amaxl*amaxu
      return
      end
