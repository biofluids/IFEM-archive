
      block data blkdta001
c although all common block variables do not have default values,
c     we comment on all the common block variables here.
c
c common block ma30e/ed holds control parameters ....
c     common /ma30ed/ lp, abort1, abort2, abort3
c the integer lp is the unit number to which the error messages are
c     sent. lp has a default value of 6.  this default value can be
c     reset by the user, if desired.  a value of 0 suppresses all
c     messages.
c the logical variables abort1,abort2,abort3 are used to control the
c     conditions under which the subroutine will terminate.
c if abort1 is .true. then the subroutine will exit  immediately on
c     detecting structural singularity.
c if abort2 is .true. then the subroutine will exit immediately on
c     detecting numerical singularity.
c if abort3 is .true. then the subroutine will exit immediately when
c     the available space in a/icn is filled up by the previously
c     decomposed, active, and undecomposed parts of the matrix.
c the default values for abort1,abort2,abort3 are set to .true.,.true.
c     and .false. respectively.
c
c the variables in the common block ma30f/fd are used to provide the
c     user with information on the decomposition.
c     common /ma30fd/ irncp, icncp, irank, minirn, minicn
c irncp and icncp are integer variables used to monitor the adequacy
c     of the allocated space in arrays irn and a/icn respectively, by
c     taking account of the number of data management compresses
c     required on these arrays. if irncp or icncp is fairly large (say
c     greater than n/10), it may be advantageous to increase the size
c     of the corresponding array(s).  irncp and icncp are initialized
c     to zero on entry to ma30a/ad and are incremented each time the
c     compressing routine ma30d/dd is entered.
c icncp is the number of compresses on a/icn.
c irncp is the number of compresses on irn.
c irank is an integer variable which gives an estimate (actually an
c     upper bound) of the rank of the matrix. on an exit with iflag
c     equal to 0, this will be equal to n.
c minirn is an integer variable which, after a successful call to
c     ma30a/ad, indicates the minimum length to which irn can be
c     reduced while still permitting a successful decomposition of the
c     same matrix. if, however, the user were to decrease the length
c     of irn to that size, the number of compresses (irncp) may be
c     very high and quite costly. if lirn is not large enough to begin
c     the decomposition on a diagonal block, minirn will be equal to
c     the value required to continue the decomposition and iflag will
c     be set to -3 or -6. a value of lirn slightly greater than this
c     (say about n/2) will usually provide enough space to complete
c     the decomposition on that block. in the event of any other
c     failure minirn gives the minimum size of irn required for a
c     successful decomposition up to that point.
c minicn is an integer variable which after a successful call to
c     ma30a/ad, indicates the minimum size of licn required to enable
c     a successful decomposition. in the event of failure with iflag=
c     -5, minicn will, if abort3 is left set to .false., indicate the
c     minimum length that would be sufficient to prevent this error in
c     a subsequent run on an identical matrix. again the user may
c     prefer to use a value of icn slightly greater than minicn for
c     subsequent runs to avoid too many conpresses (icncp). in the
c     event of failure with iflag equal to any negative value except
c     -4, minicn will give the minimum length to which licn could be
c     reduced to enable a successful decomposition to the point at
c     which failure occurred.  notice that, on a successful entry
c     idisp(2) gives the amount of space in a/icn required for the
c     decomposition while minicn will usually be slightly greater
c     because of the need for "elbow room".  if the user is very
c     unsure how large to make licn, the variable minicn can be used
c     to provide that information. a preliminary run should be
c     performed with abort3 left set to .false. and licn about 3/2
c     times as big as the number of non-zeros in the original matrix.
c     unless the initial problem is very sparse (when the run will be
c     successful) or fills in extremely badly (giving an error return
c     with iflag equal to -4), an error return with iflag equal to -5
c     should result and minicn will give the amount of space required
c     for a successful decomposition.
c
c common block ma30g/gd is used by the ma30b/bd entry only.
c     common /ma30gd/ eps, rmin
c eps is a real/double precision variable. it is used to test for
c     small pivots. its default value is 1.0e-4 (1.0d-4 in d version).
c     if the user sets eps to any value greater than 1.0, then no
c     check is made on the size of the pivots. although the absence of
c     such a check would fail to warn the user of bad instability, its
c     absence will enable ma30b/bd to run slightly faster. an  a
c     posteriori  check on the stability of the factorization can be
c     obtained from mc24a/ad.
c rmin is a real/double precision variable which gives the user some
c     information about the stability of the decomposition.  at each
c     stage of the lu decomposition the magnitude of the pivot apiv
c     is compared with the largest off-diagonal entry currently in its
c     row (row of u), rowmax say. if the ratio
c                       min (apiv/rowmax)
c     where the minimum is taken over all the rows, is less than eps
c     then rmin is set to this minimum value and iflag is returned
c     with the value +i where i is the row in which this minimum
c     occurs.  if the user sets eps greater than one, then this test
c     is not performed. in this case, and when there are no small
c     pivots rmin will be set equal to eps.
c
c common block ma30h/hd is used by ma30c/cd only.
c     common /ma30hd/ resid
c resid is a real/double precision variable. in the case of singular
c     or rectangular matrices its final value will be equal to the
c     maximum residual for the unsatisfied equations; otherwise its
c     value will be set to zero.
c
c common  block ma30i/id controls the use of drop tolerances, the
c     modified pivot option and the the calculation of the largest
c     entry in the factorization process. this common block was added
c     to the ma30 package in february, 1983.
c     common /ma30id/ tol, big, ndrop, nsrch, lbig
c tol is a real/double precision variable.  if it is set to a positive
c     value, then ma30a/ad will drop from the factors any non-zero
c     whose modulus is less than tol.  the factorization will then
c     require less storage but will be inaccurate.  after a run of
c     ma30a/ad where entries have been dropped, ma30b/bd  should not
c     be called.  the default value for tol is 0.0.
c big is a real/double precision variable.  if lbig has been set to
c     .true., big will be set to the largest entry encountered during
c     the factorization.
c ndrop is an integer variable. if tol has been set positive, on exit
c     from ma30a/ad, ndrop will hold the number of entries dropped
c     from the data structure.
c nsrch is an integer variable. if nsrch is set to a value less than
c     or equal to n, then a different pivot option will be employed by
c     ma30a/ad.  this may result in different fill-in and execution
c     time for ma30a/ad. if nsrch is less than or equal to n, the
c     workspace arrays lastc and nextc are not referenced by ma30a/ad.
c     the default value for nsrch is 32768.
c lbig is a logical variable. if lbig is set to .true., the value of
c     the largest entry encountered in the factorization by ma30a/ad
c     is returned in big.  setting lbig to .true.  will marginally
c     increase the factorization time for ma30a/ad and will increase
c     that for ma30b/bd by about 20%.  the default value for lbig is
c     .false.
c
      real* 8 eps, rmin, tol, big
      logical abort1, abort2, abort3, lbig
      common /ma30ed/ lp, abort1, abort2, abort3
      common /ma30gd/ eps, rmin
      common /ma30id/ tol, big, ndrop, nsrch, lbig
      data eps /1.0e-4/, tol /0.0/, big /0.0/
      data lp /6/, nsrch /32768/
      data lbig /.false./
      data abort1 /.true./, abort2 /.true./, abort3 /.false./
      end
