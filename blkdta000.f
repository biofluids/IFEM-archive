
      block data blkdta000
c
c comments on all the common block variables are given here even
c     though some are not initialized by block data.
c lp,mp are used by the subroutine as the unit numbers for its warning
c     and diagnostic messages. default value for both is 6 (for line
c     printer output). the user can either reset them to a different
c     stream number or suppress the output by setting them to zero.
c     while lp directs the output of error diagnostics from the
c     principal subroutines and internally called subroutines, mp
c     controls only the output of a message which warns the user that he
c     has input two or more non-zeros a(i), . . ,a(k) with the same row
c     and column indices.  the action taken in this case is to proceed
c     using a numerical value of a(i)+...+a(k). in the absence of other
c     errors, iflag will equal -14 on exit.
c lblock is a logical variable which controls an option of first
c     preordering the matrix to block lower triangular form (using
c     harwell subroutine mc23a). the preordering is performed if lblock
c     is equal to its default value of .true. if lblock is set to
c     .false. , the option is not invoked and the space allocated to
c     ikeep can be reduced to 4*n+1.
c grow is a logical variable. if it is left at its default value of
c     .true. , then on return from ma28a/ad or ma28b/bd, w(1) will give
c     an estimate (an upper bound) of the increase in size of elements
c     encountered during the decomposition. if the matrix is well
c     scaled, then a high value for w(1), relative to the largest entry
c     in the input matrix, indicates that the lu decomposition may be
c     inaccurate and the user should be wary of his results and perhaps
c     increase u for subsequent runs.  we would like to emphasise that
c     this value only relates to the accuracy of our lu decomposition
c     and gives no indication as to the singularity of the matrix or the
c     accuracy of the solution.  this upper bound can be a significant
c     overestimate particularly if the matrix is badly scaled. if an
c     accurate value for the growth is required, lbig (q.v.) should be
c     set to .true.
c eps,rmin are real variables. if, on entry to ma28b/bd, eps is less
c     than one, then rmin will give the smallest ratio of the pivot to
c     the largest element in the corresponding row of the upper
c     triangular factor thus monitoring the stability of successive
c     factorizations. if rmin becomes very large and w(1) from
c     ma28b/bd is also very large, it may be advisable to perform a
c     new decomposition using ma28a/ad.
c resid is a real variable which on exit from ma28c/cd gives the value
c     of the maximum residual over all the equations unsatisfied because
c     of dependency (zero pivots).
c irncp,icncp are integer variables which monitor the adequacy of "elbow
c     room" in irn and a/icn respectively. if either is quite large (say
c     greater than n/10), it will probably pay to increase the size of
c     the corresponding array for subsequent runs. if either is very low
c     or zero then one can perhaps save storage by reducing the size of
c     the corresponding array.
c minirn,minicn are integer variables which, in the event of a
c     successful return (iflag ge 0 or iflag=-14) give the minimum size
c     of irn and a/icn respectively which would enable a successful run
c     on an identical matrix. on an exit with iflag equal to -5, minicn
c     gives the minimum value of icn for success on subsequent runs on
c     an identical matrix. in the event of failure with iflag= -6, -4,
c     -3, -2, or -1, then minicn and minirn give the minimum value of
c     licn and lirn respectively which would be required for a
c     successful decomposition up to the point at which the failure
c     occurred.
c irank is an integer variable which gives an upper bound on the rank of
c     the matrix.
c abort1 is a logical variable with default value .true.  if abort1 is
c     set to .false.  then ma28a/ad will decompose structurally singular
c     matrices (including rectangular ones).
c abort2 is a logical variable with default value .true.  if abort2 is
c     set to .false. then ma28a/ad will decompose numerically singular
c     matrices.
c idisp is an integer array of length 2. on output from ma28a/ad, the
c     indices of the diagonal blocks of the factors lie in positions
c     idisp(1) to idisp(2) of a/icn. this array must be preserved
c     between a call to ma28a/ad and subsequent calls to ma28b/bd,
c     ma28c/cd or ma28i/id.
c tol is a real variable.  if it is set to a positive value, then any
c     non-zero whose modulus is less than tol will be dropped from the
c     factorization.  the factorization will then require less storage
c     but will be inaccurate.  after a run of ma28a/ad with tol positive
c     it is not possible to use ma28b/bd and the user is recommended to
c     use ma28i/id to obtain the solution.  the default value for tol is
c     0.0.
c themax is a real variable.  on exit from ma28a/ad, it will hold the
c     largest entry of the original matrix.
c big is a real variable. if lbig has been set to .true., big will hold
c     the largest entry encountered during the factorization by ma28a/ad
c     or ma28b/bd.
c dxmax is a real variable. on exit from ma28i/id, dxmax will be set to
c     the largest component of the solution.
c errmax is a real variable.  on exit from ma28i/id, if maxit is
c     positive, errmax will be set to the largest component in the
c     estimate of the error.
c dres is a real variable.  on exit from ma28i/id, if maxit is positive,
c     dres will be set to the largest component of the residual.
c cgce is a real variable. it is used by ma28i/id to check the
c     convergence rate.  if the ratio of successive corrections is
c     not less than cgce then we terminate since the convergence
c     rate is adjudged too slow.
c ndrop is an integer variable. if tol has been set positive, on exit
c     from ma28a/ad, ndrop will hold the number of entries dropped from
c     the data structure.
c maxit is an integer variable. it is the maximum number of iterations
c     performed by ma28i/id. it has a default value of 16.
c noiter is an integer variable. it is set by ma28i/id to the number of
c     iterative refinement iterations actually used.
c nsrch is an integer variable. if nsrch is set to a value less than n,
c     then a different pivot option will be employed by ma28a/ad.  this
c     may result in different fill-in and execution time for ma28a/ad.
c     if nsrch is less than or equal to n, the workspace array iw can be
c     reduced in length.  the default value for nsrch is 32768.
c istart is an integer variable. if istart is set to a value other than
c     zero, then the user must supply an estimate of the solution to
c     ma28i/id.  the default value for istart is zero.
c lbig is a logical variable. if lbig is set to .true., the value of the
c     largest element encountered in the factorization by ma28a/ad or
c     ma28b/bd is returned in big.  setting lbig to .true.  will
c     increase the time for ma28a/ad marginally and that for ma28b/bd
c     by about 20%.  the default value for lbig is .false.
c
      real* 8 eps, rmin, resid, tol, themax, big, dxmax,
     * errmax, dres, cgce
      logical lblock, grow, abort1, abort2, lbig
      common /ma28ed/ lp, mp, lblock, grow
      common /ma28fd/ eps, rmin, resid, irncp, icncp, minirn, minicn,
     * irank, abort1, abort2
      common /ma28gd/ idisp(2)
      common /ma28hd/ tol, themax, big, dxmax, errmax, dres, cgce,
     * ndrop, maxit, noiter, nsrch, istart, lbig
      data eps /1.0e-4/, tol /0.0/, cgce /0.5/
      data maxit /16/
      data lp /6/, mp /6/, nsrch /32768/, istart /0/
      data lblock /.true./, grow /.true./, lbig /.false./
      data abort1 /.true./, abort2 /.true./
      end
