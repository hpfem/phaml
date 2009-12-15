*---------------------------------------------------------------------!
*                                PHAML                                !
*                                                                     !
* The Parallel Hierarchical Adaptive MultiLevel code for solving      !
* linear elliptic partial differential equations of the form          !
* (PUx)x + (QUy)y + RU = F on 2D polygonal domains with mixed         !
* boundary conditions, and eigenvalue problems where F is lambda*U.   !
*                                                                     !
* PHAML is public domain software.  It was produced as part of work   !
* done by the U.S. Government, and is not subject to copyright in     !
* the United States.                                                  !
*                                                                     !
*     William F. Mitchell                                             !
*     Mathematical and Computational Sciences Division                !
*     National Institute of Standards and Technology                  !
*     william.mitchell@nist.gov                                       !
*     http://math.nist.gov/phaml                                      !
*                                                                     !
*---------------------------------------------------------------------!

* This file contains the Conjugate Gradient and GMRES routines from the
* Template book (see reference below).  The code was obtained from netlib.

* Changes that have been made since obtaining it are flagged with WFM.
* 05/25/07 Commented out second declaration of NDX1 and NDX2 in CGREVCOM
* 05/25/07 In four places, the return from a reverse communication branched
*          to a point inside a block (two if-then and two do) which is not
*          legal.  Made changes to remove that.
* 05/31/07 Fixed a bug in the scaling of the initial residual in GMRES
* 10/18/07 Added printing of convergence

*
* WFM 10/18/07 Added PRCONV to control printing of convergence

      SUBROUTINE CG( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $               PSOLVE, INFO, PRCONV )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      DOUBLE PRECISION   RESID
* WFM 10/18/07 Added PRCONV to control printing of convergence
      LOGICAL            PRCONV
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), B( * ),  WORK( * )
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL           MATVEC, PSOLVE
*     ..
*
*  Purpose
*  =======
*
*  CG solves the linear system Ax = b using the
*  Conjugate Gradient iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*  --Done in CGREVCOM.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*          Set by CGREVCOM.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension ( * ).
*          Workspace for residual, direction vector, etc.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*          Set by CGREVCOM.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*          Set by CGREVCOM.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must 
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*          Set by CGREVCOM()
*  ============================================================
*
*     ..
*     .. Local Scalars ..
*This variable used to communicate requests between CG() and CGREVCOM()
*CG -> CGREVCOM: 1 = init, 
*                2 = use saved state to resume flow.
*CGREVCOM -> CG: -1 = done, return to main, 
*                 1 = matvec using SCLR1/2, NDX1/2 
*                 2 = solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*     Arg/Result indices into WORK[].
      INTEGER          NDX1, NDX2
*     Scalars passed from CGREVCOM to CG.
      DOUBLE PRECISION SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      DOUBLE PRECISION TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         CGREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == R; 2 == Z; 3 == P; 4 == Q; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

      CALL CGREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $              NDX1, NDX2, SCLR1, SCLR2, IJOB)

*     On a return from CGREVCOM() we use the table (CGREVCOM -> CG)
*     to figure out what is reqd.
      IF (IJOB .eq. -1) THEN
*        revcom wants to terminate, so do it.
         GOTO 2
      ELSEIF (IJOB .EQ. 1) THEN
*        call matvec.
         CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
      ELSEIF (IJOB .EQ. 2) THEN
*        call solve.
         CALL PSOLVE(WORK(NDX1), WORK(NDX2))
      ELSEIF (IJOB .EQ. 3) THEN
*        call matvec with X.
         CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
      ELSEIF (IJOB .EQ. 4) THEN
*        do stopping test 2
*        if first time, set INFO so that BNRM2 is computed.
         IF( FTFLG ) INFO = -1
* WFM 10/18/07 added PRCONV
         CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO,
     *                  PRCONV)
         FTFLG = .FALSE.
      ENDIF
*
*     done what revcom asked, set IJOB & go back to it.
      IJOB = 2
      GOTO 1
*
*     come here to terminate
 2    CONTINUE

      RETURN
*
*     End of CG
*
      END
*
      SUBROUTINE CGREVCOM( N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                     NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      DOUBLE PRECISION   RESID
      INTEGER            NDX1, NDX2
      DOUBLE PRECISION   SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), B( * ),  WORK( LDW,* )
*
*     (output) for matvec and solve. These index into WORK[]
* WFM 5/25/07 commented out this line because they are already declared
*      INTEGER NDX1, NDX2
*     ..
*
*  Purpose
*  =======
*
*  CG solves the linear system Ax = b using the
*  Conjugate Gradient iterative method with preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,4).
*          Workspace for residual, direction vector, etc.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter.
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) DOUBLE PRECISION.
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2
*  ============================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER          MAXIT, R, Z, P, Q, NEED1, NEED2
      DOUBLE PRECISION ALPHA, BETA, RHO, RHO1, DDOT, DNRM2,
     $                 TOL
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL

* WFM 5/25/07 the "move 2 continue" fix uses this
      LOGICAL RET2
* end WFM

*
*     saving all.
      SAVE
*     ..
*     .. External Routines ..
      EXTERNAL         DAXPY, DCOPY, DDOT, DNRM2
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
* init.
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL = RESID
*
*     Alias workspace columns.
*
      R = 1
      Z = 2
      P = 3
      Q = 4
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((Q - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((Q - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set initial residual.
*
      CALL DCOPY( N, B, 1, WORK(1,R), 1 )
* WFM 5/25/07 the "move 2 continue" fix uses this
      RET2 = .FALSE.
* end WFM
      IF ( DNRM2( N, X, 1 ).NE.ZERO ) THEN
* WFM 5/25/07 the "move 2 continue" fix uses this
         RET2 = .TRUE.
* end WFM

*********CALL MATVEC( -ONE, X, ONE, WORK(1,R) )
*
*        Set args for revcom return
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN

* WFM 5/25/07 this is the "move 2 continue" fix
*             On return with RLBL .eq. 2 the original code jumps into
*             this if-then block, which Fortran does not allow.  The
*             original code (with ***** in front of it) is:

******
**********************
***** 2       CONTINUE
**********************
******
*****         IF ( DNRM2( N, WORK(1,R), 1 ).LT.TOL ) GO TO 30
*****      ENDIF
******

*             This is the code with 2 continue moved outside the block:

      ENDIF

*****************
 2    CONTINUE
*****************

      IF (RET2) THEN
         IF ( DNRM2( N, WORK(1,R), 1 ).LT.TOL ) GO TO 30
      ENDIF

* end WFM

      ITER = 0
*
   10 CONTINUE
*
*        Perform Preconditioned Conjugate Gradient iteration.
*
         ITER = ITER + 1
*
*        Preconditioner Solve.
*
*********CALL PSOLVE( WORK(1,Z), WORK(1,R) )
*
         NDX1 = ((Z - 1) * LDW) + 1
         NDX2 = ((R - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
*****************
*
         RHO = DDOT( N, WORK(1,R), 1, WORK(1,Z), 1 )
*
*        Compute direction vector P.
*
         IF ( ITER.GT.1 ) THEN
            BETA = RHO / RHO1
            CALL DAXPY( N, BETA, WORK(1,P), 1, WORK(1,Z), 1 )
*
            CALL DCOPY( N, WORK(1,Z), 1, WORK(1,P), 1 )
         ELSE
            CALL DCOPY( N, WORK(1,Z), 1, WORK(1,P), 1 )
         ENDIF
*
*        Compute scalar ALPHA (save A*P to Q).
*
*********CALL MATVEC( ONE, WORK(1,P), ZERO, WORK(1,Q) )
*
         NDX1 = ((P - 1) * LDW) + 1
         NDX2 = ((Q - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ZERO
         RLBL = 4
         IJOB = 1
         RETURN
*
*****************
 4       CONTINUE
*****************
*
         ALPHA =  RHO / DDOT( N, WORK(1,P), 1, WORK(1,Q), 1 )
*
*        Compute current solution vector X.
*
         CALL DAXPY( N, ALPHA, WORK(1,P), 1, X, 1 )
*
*        Compute residual vector R, find norm,
*        then check for tolerance.
*
         CALL DAXPY( N, -ALPHA,  WORK(1,Q), 1, WORK(1,R), 1 )
*
*********RESID = DNRM2( N, WORK(1,R), 1 ) / BNRM2
*********IF ( RESID.LE.TOL ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 5
         IJOB = 4
         RETURN
*
*****************
 5       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
         RHO1 = RHO
*
         GO TO 10
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of CGREVCOM
*
      END
*
* WFM 10/18/07 Added PRCONV to control printing of convergence

      SUBROUTINE GMRES( N, B, X, RESTRT, WORK, LDW, WORK2, LDW2, 
     $                  ITER, RESID, MATVEC, PSOLVE, INFO, PRCONV )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     EiITERkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, RESTRT, LDW, LDW2, ITER, INFO
      DOUBLE PRECISION   RESID
* WFM 10/18/07 Added PRCONV to control printing of convergence
      LOGICAL            PRCONV
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( * ), X( * ), WORK( * ), WORK2( * )
*     ..
*     .. Function Arguments ..
*
      EXTERNAL           MATVEC, PSOLVE
*
*  Purpose
*  =======
*
*  GMRES solves the linear system Ax = b using the
*  Generalized Minimal Residual iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*  RESTRT  (input) INTEGER
*          Restart parameter, <= N. This parameter controls the amount
*          of memory required for matrix WORK2.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,5).
*          Note that if the initial guess is the zero vector, then 
*          storing the initial residual is not necessary.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  WORK2   (workspace) DOUBLE PRECISION array, dimension (LDW2,2*RESTRT+2).
*          This workspace is used for constructing and storing the
*          upper Hessenberg matrix. The two extra columns are used to
*          store the Givens rotation matrices.
*
*  LDW2    (input) INTEGER
*          The leading dimension of the array WORK2.
*          LDW2 >= max(1,RESTRT).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable error tolerance.
*          On ouput, the norm of the residual vector if solution
*          approximated to tolerance, otherwise reset to input
*          tolerance.
*
*  INFO    (output) INTEGER
*          =  0:  successful exit
*          =  1:  maximum number of iterations performed;
*                 convergence not achieved.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product A*x = y.
*          Vector x must remain unchanged. The solution is
*          over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( X, Y )
*
*  ============================================================
*     ..
*     .. Parameters ..
      INTEGER            OFSET
      PARAMETER        ( OFSET = 1000 )
*
*This variable used to communicate requests between GMRES & GMRESREVCOM
*GMRES -> GMRESREVCOM: 1 = init, 
*                      2 = use saved state to resume flow.
*GMRESREVCOM -> GMRES: -1 = done, return to main, 
*                       1 = matvec using SCLR1/2, NDX1/2 
*                       2 = solve using NDX1/2
*                       3 = matvec using WORK2 NDX1 & WORK NDX2
      INTEGER          IJOB
      LOGICAL          FTFLG
*
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from GMRESREVCOM to GMRES
      DOUBLE PRECISION SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      DOUBLE PRECISION TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         GMRESREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ELSE IF ( LDW2.LT.RESTRT ) THEN
         INFO = -4
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*     Note: V, and GIV contain # of history vectors each.
*     To access the i'th vector in V, add i to V*OFSET 1<=i<=RESTRT
*     To access the i'th vector in GIV, add i to GIV*OFSET 1<=i<=RESTRT
*
*     1 == R; 2 == S; 3 == W; 4 == Y; 5 == AV; 6 == H; 7*OFSET+i == V;
*     8*OFSET+i == GIV; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

          CALL GMRESREVCOM(N, B, X, RESTRT, WORK, LDW, WORK2, LDW2, 
     $                     ITER, RESID, INFO, NDX1, NDX2, SCLR1, 
     $                     SCLR2, IJOB)


*         On a return from REVCOM() we use the table
*         to decode IJOB.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec directly with X.
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 2) THEN
*           call solve.
            CALL PSOLVE(WORK(NDX1), WORK(NDX2))
          ELSEIF (IJOB .eq. 3) THEN
*           call matvec.
            CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
         ELSEIF (IJOB .EQ. 4) THEN
*           do stopping test 2
*           if first time, set INFO so that BNRM2 is computed.
            IF( FTFLG ) INFO = -1
* WFM 10/18/07 added PRCONV
            CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO,
     *                     PRCONV)
            FTFLG = .FALSE.
         ENDIF
*
*         done what revcom asked, set IJOB & go back to it.
          IJOB = 2
          GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
      RETURN
*
*     End of GMRES.f
*
      END
*
*
       SUBROUTINE GMRESREVCOM(N, B, X, RESTRT, WORK, LDW, WORK2,
     $                  LDW2, ITER, RESID, INFO, NDX1, NDX2, SCLR1, 
     $                  SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     EiITERkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, RESTRT, LDW, LDW2, ITER, INFO
      DOUBLE PRECISION   RESID
      INTEGER            NDX1, NDX2
      DOUBLE PRECISION   SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( * ), X( * ), WORK( LDW,* ), WORK2( LDW2,* )
*     ..
*
*  Purpose
*  =======
*
*  GMRES solves the linear system Ax = b using the
*  Generalized Minimal Residual iterative method with preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*  RESTRT  (input) INTEGER
*          Restart parameter, <= N. This parameter controls the amount
*          of memory required for matrix WORK2.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,5).
*          Note that if the initial guess is the zero vector, then 
*          storing the initial residual is not necessary.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  WORK2   (workspace) DOUBLE PRECISION array, dimension (LDW2,2*RESTRT+2).
*          This workspace is used for constructing and storing the
*          upper Hessenberg matrix. The two extra columns are used to
*          store the Givens rotation matrices.
*
*  LDW2    (input) INTEGER
*          The leading dimension of the array WORK2.
*          LDW2 >= max(1,RESTRT).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable error tolerance.
*          On ouput, the norm of the residual vector if solution
*          approximated to tolerance, otherwise reset to input
*          tolerance.
*
*  INFO    (output) INTEGER
*          =  0:  successful exit
*          =  1:  maximum number of iterations performed;
*                 convergence not achieved.
*            -5: Erroneous NDX1/NDX2 in INIT call.
*            -6: Erroneous RLBL.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) DOUBLE PRECISION.
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  ============================================================
*
*     .. Parameters ..
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      INTEGER             OFSET
      PARAMETER         ( OFSET = 1000 )
*     ..
*     .. Local Scalars ..
      INTEGER             I, MAXIT, AV, GIV, H, R, S, V, W, Y,
     $                    NEED1, NEED2
      DOUBLE PRECISION    BNRM2, RNORM, TOL, APPROXRES, DDOT, DNRM2
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL

* WFM 5/25/07 the "move 2 continue" fix uses this
      LOGICAL RET2
* end WFM

*
*     saving all.
      SAVE
*
*     ..
*     .. External Routines ..
      EXTERNAL            DAXPY, DCOPY, DDOT, DNRM2
*     ..
*     .. Executable Statements ..
*
* Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 200
      ENDIF
*
* init.
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R   = 1
      S   = 2
      W   = 3
      Y   = 4
      AV  = 5
      V   = 6
*
      H   = 1
      GIV = H + RESTRT
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((W - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((Y - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((AV - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( ( NDX1.GT.V*OFSET ) .AND.
     $           ( NDX1.LE.V*OFSET+RESTRT ) ) THEN
            NEED1 = ((NDX1-V*OFSET - 1) * LDW) + 1
         ELSEIF( ( NDX1.GT.GIV*OFSET ) .AND.
     $           ( NDX1.LE.GIV*OFSET+RESTRT ) ) THEN
            NEED1 = ((NDX1-GIV*OFSET - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 100
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((W - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((Y - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((AV - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( ( NDX2.GT.V*OFSET ) .AND.
     $           ( NDX2.LE.V*OFSET+RESTRT ) ) THEN
            NEED2 = ((NDX2-V*OFSET - 1) * LDW) + 1
         ELSEIF( ( NDX2.GT.GIV*OFSET ) .AND.
     $           ( NDX2.LE.GIV*OFSET+RESTRT ) ) THEN
            NEED2 = ((NDX2-GIV*OFSET - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 100
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set initial residual.
*
      CALL DCOPY( N, B, 1, WORK(1,R), 1 )
* WFM 5/25/07 the "move 2 continue" fix uses this
      RET2 = .FALSE.
* end WFM
      IF ( DNRM2( N, X, 1 ).NE.ZERO ) THEN
* WFM 5/25/07 the "move 2 continue" fix uses this
         RET2 = .TRUE.
* end WFM
*********CALL MATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using X directly
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 1
         RETURN
* WFM 5/25/07 this is the "move 2 continue" fix
*             On return with RLBL .eq. 2 the original code jumps into
*             this if-then block, which Fortran does not allow.  The
*             original code (with ***** in front of it) is:

******
**********************
***** 2       CONTINUE
**********************
******
*****         IF ( DNRM2( N, WORK(1,R), 1 ).LT.TOL ) GO TO 200
*****      ENDIF

*             This is the code with 2 continue moved outside the block:

      ENDIF

*****************
 2    CONTINUE
*****************

      IF (RET2) THEN
         IF ( DNRM2( N, WORK(1,R), 1 ).LT.TOL ) GO TO 200
      ENDIF

* end WFM

      BNRM2 = DNRM2( N, B, 1 )
      IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
*
      ITER = 0
   10 CONTINUE
*
         ITER = ITER + 1
*
*        Construct the first column of V, and initialize S to the
*        elementary vector E1 scaled by RNORM.
*
*********CALL PSOLVE( WORK( 1,V ), WORK( 1,R ) )
*
         NDX1 = ((V - 1) * LDW) + 1
         NDX2 = ((R - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
*****************
*
         RNORM = DNRM2( N, WORK( 1,V ), 1 )
* WFM 5/31/07 should be dividing by RNORM, not multiplying
*  was    CALL DSCAL( N, RNORM, WORK( 1,V ), 1 )
         CALL DSCAL( N, ONE / RNORM, WORK( 1,V ), 1 )
         CALL ELEMVEC( 1, N, RNORM, WORK( 1,S ) )
*
* WFM 5/25/07 RLBL .eq. 4 and RLBL .eq. 5 jump into this DO loop, which
*             Fortran does not allow.  Replace the DO loop with a GOTO loop.
*             The original statement was:
*****         DO 50 I = 1, RESTRT

         I = 1
 501     CONTINUE

* end WFM

************CALL MATVEC( ONE, WORK( 1,V+I-1 ), ZERO, WORK( 1,AV ) )
*
         NDX1 = ((V+I-1 - 1) * LDW) + 1
         NDX2 = ((AV    - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ZERO
         RLBL = 4
         IJOB = 3
         RETURN
*
*****************
 4       CONTINUE
*****************
*
*********CALL PSOLVE( WORK( 1,W ), WORK( 1,AV ) )
*
         NDX1 = ((W  - 1) * LDW) + 1
         NDX2 = ((AV - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 5
         IJOB = 2
         RETURN
*
*****************
 5       CONTINUE
*****************
*
*           Construct I-th column of H so that it is orthnormal to 
*           the previous I-1 columns.
*
            CALL ORTHOH( I, N, WORK2( 1,I+H-1 ), WORK( 1,V ), LDW,
     $                   WORK( 1,W ) )
*
            IF ( I.GT.1 )
*
*              Apply Givens rotations to the I-th column of H. This
*              effectively reduces the Hessenberg matrix to upper
*              triangular form during the RESTRT iterations.
*
     $         CALL APPLYGIVENS( I, WORK2( 1,I+H-1 ), WORK2( 1,GIV ),
     $                           LDW2 )
*
*           Approxiyate residual norm. Check tolerance. If okay, compute
*           final approximation vector X and quit.
*
            RESID = APPROXRES( I, WORK2( 1,I+H-1 ), WORK( 1,S ),
     $                         WORK2( 1,GIV ), LDW2 ) / BNRM2
            IF ( RESID.LE.TOL ) THEN
               CALL UPDATE(I, N, X, WORK2( 1,H ), LDW2, 
     $                     WORK(1,Y), WORK(1,S), WORK( 1,V ), LDW)
               GO TO 200
            ENDIF
* WFM 5/25/07 RLBL .eq. 4 and RLBL .eq. 5 jump into this DO loop, which
*             Fortran does not allow.  Replace the DO loop with a GOTO loop.
*             The original statement was:
*****   50    CONTINUE

            I = I + 1
            IF (I .LE. RESTRT) GO TO 501
   50    CONTINUE

* end WFM

*
*        Compute current solution vector X.
*
         CALL UPDATE(RESTRT, N, X, WORK2( 1,H ), LDW2,
     $               WORK(1,Y), WORK( 1,S ), WORK( 1,V ), LDW )
*
*        Compute residual vector R, find norm,
*        then check for tolerance.
*
         CALL DCOPY( N, B, 1, WORK( 1,R ), 1 )
*********CALL MATVEC( -ONE, X, ONE, WORK( 1,R ) )
*
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = -ONE
         SCLR2 = ONE
         RLBL = 6
         IJOB = 1
         RETURN
*
*****************
 6       CONTINUE
*****************
*
         WORK( I+1,S ) = DNRM2( N, WORK( 1,R ), 1 )
*
*********RESID = WORK( I+1,S ) / BNRM2
*********IF ( RESID.LE.TOL  ) GO TO 200
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 7
         IJOB = 4
         RETURN
*
*****************
 7       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 200
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 100
         ENDIF
*
         GO TO 10
*
  100 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
  200 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1

      RETURN
*
*     End of GMRESREVCOM
*
      END
*
*     =========================================================
      SUBROUTINE ORTHOH( I, N, H, V, LDV, W )
*
      INTEGER            I, N, LDV
      DOUBLE PRECISION   H( * ), W( * ), V( LDV,* )
*
*     Construct the I-th column of the upper Hessenberg matrix H
*     using the Gram-Schmidt process on V and W.
*
      INTEGER            K
      DOUBLE PRECISION   DDOT, DNRM2, ONE
      PARAMETER        ( ONE = 1.0D+0 )
      EXTERNAL           DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*
      DO 10 K = 1, I
         H( K ) = DDOT( N, W, 1, V( 1,K ), 1 )
         CALL DAXPY( N, -H( K ), V( 1,K ), 1, W, 1 )
   10 CONTINUE
      H( I+1 ) = DNRM2( N, W, 1 )
      CALL DCOPY( N, W, 1, V( 1,I+1 ), 1 )
      CALL DSCAL( N, ONE / H( I+1 ), V( 1,I+1 ), 1 )
*
      RETURN
*
      END
*     =========================================================
      SUBROUTINE APPLYGIVENS( I, H, GIVENS, LDG )
*
      INTEGER            I, LDG
      DOUBLE PRECISION   H( * ), GIVENS( LDG,* )
*
*     This routine applies a sequence of I-1 Givens rotations to
*     the I-th column of H. The Givens parameters are stored, so that
*     the first I-2 Givens rotation matrices are known. The I-1st
*     Givens rotation is computed using BLAS 1 routine DROTG. Each
*     rotation is applied to the 2x1 vector [H( J ), H( J+1 )]',
*     which results in H( J+1 ) = 0.
*
      INTEGER            J
*      DOUBLE PRECISION   TEMP
      EXTERNAL           DROTG
*
*     .. Executable Statements ..
*
*     Construct I-1st rotation matrix.
*
*     CALL DROTG( H( I ), H( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
      CALL GETGIV( H( I ), H( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
*
*     Apply 1,...,I-1st rotation matrices to the I-th column of H.
*
      DO 10 J = 1, I-1
         CALL ROTVEC( H( J ), H( J+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
*        TEMP     =  GIVENS( J,1 ) * H( J ) + GIVENS( J,2 ) * H( J+1 ) 
*        H( J+1 ) = -GIVENS( J,2 ) * H( J ) + GIVENS( J,1 ) * H( J+1 )
*        H( J ) = TEMP
 10   CONTINUE
*
      RETURN
*
      END
*
*     ===============================================================
      DOUBLE PRECISION FUNCTION APPROXRES( I, H, S, GIVENS, LDG )
*
      INTEGER            I, LDG
      DOUBLE PRECISION   H( * ), S( * ), GIVENS( LDG,* )
*
*     This function allows the user to approximate the residual
*     using an updating scheme involving Givens rotations. The
*     rotation matrix is formed using [H( I ),H( I+1 )]' with the
*     intent of zeroing H( I+1 ), but here is applied to the 2x1
*     vector [S(I), S(I+1)]'.
*
      INTRINSIC          DABS
      EXTERNAL           DROTG
*
*     .. Executable Statements ..
*
*     CALL DROTG( H( I ), H( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
      CALL GETGIV( H( I ), H( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
      CALL ROTVEC( S( I ), S( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
*
      APPROXRES = DABS( S( I+1 ) )
*
      RETURN
*
      END
*     ===============================================================
      SUBROUTINE UPDATE( I, N, X, H, LDH, Y, S, V, LDV )
*
      INTEGER            N, I, J, LDH, LDV
      DOUBLE PRECISION   X( * ), Y( * ), S( * ), H( LDH,* ), V( LDV,* )
      EXTERNAL           DAXPY, DCOPY, DTRSV
*
*     Solve H*y = s for upper triangualar H.
*
      CALL DCOPY( I, S, 1, Y, 1 )
      CALL DTRSV( 'UPPER', 'NOTRANS', 'NONUNIT', I, H, LDH, Y, 1 )
*
*     Compute current solution vector X.
*
      DO 10 J = 1, I
         CALL DAXPY( N, Y( J ), V( 1,J ), 1, X, 1 )
   10 CONTINUE
*
      RETURN
*
      END
*
*     ===============================================================
      SUBROUTINE GETGIV( A, B, C, S )
*
      DOUBLE PRECISION   A, B, C, S, TEMP, ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*
      IF ( B.EQ.ZERO ) THEN
         C = ONE
         S = ZERO
      ELSE IF ( ABS( B ).GT.ABS( A ) ) THEN
         TEMP = -A / B
         S = ONE / SQRT( ONE + TEMP**2 )
         C = TEMP * S
      ELSE
         TEMP = -B / A
         C = ONE / SQRT( ONE + TEMP**2 )
         S = TEMP * C
      ENDIF
*
      RETURN
*
      END
*
*     ================================================================
      SUBROUTINE ROTVEC( X, Y, C, S )
*
      DOUBLE PRECISION   X, Y, C, S, TEMP

*
      TEMP = C * X - S * Y
      Y    = S * X + C * Y
      X    = TEMP
*
      RETURN
*
      END
*
*     ===============================================================
      SUBROUTINE ELEMVEC( I, N, ALPHA, E )
*
*     Construct the I-th elementary vector E, scaled by ALPHA.
*
      INTEGER            I, J, N
      DOUBLE PRECISION   ALPHA, E( * )
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*
      DO 10 J = 1, N
         E( J ) = ZERO
   10 CONTINUE
      E( I ) = ALPHA
*
      RETURN
*
      END
*
* WFM 10/18/07 added PRCONV
      SUBROUTINE STOPTEST2( N, R, B, BNRM2, RESID, TOL, INFO, PRCONV )
*
*     .. Scalar Arguments ..
      INTEGER            N, INFO
      DOUBLE PRECISION   RESID, TOL, BNRM2
* WFM 10/18/07 added PRCONV
      LOGICAL            PRCONV
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   R( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  Computes the stopping criterion 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  INFO    (output) INTEGER
*          On exit, 1/0 depending on whether stopping criterion
*          was met or not.
*
*  BNRM2   (input/output) DOUBLE PRECISION.
*          On first time entry, will be -1.0.
*          On first time exit will contain norm2(B)
*          On all subsequent entry/exit's unchanged.
*
*  RESID   (output) DOUBLE PRECISION.
*          On exit, the computed stopping measure.
*
*  TOL     (input) DOUBLE PRECISION.
*          On input, the allowable convergence measure.
*
*  R       (input) DOUBLE PRECISION array, dimension N.
*          On entry, the residual.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  BLAS CALLS:   DNRM2
*  ============================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. External Routines ..
      DOUBLE PRECISION DNRM2
      EXTERNAL         DNRM2
* WFM 10/18/07 save the old residual to print convergence rate and iteration count
      INTEGER, SAVE :: ITER
      DOUBLE PRECISION, SAVE :: OLDRES = 0.0D0
*     ..
*     .. Executable Statements ..
*
      IF( INFO.EQ.-1 ) THEN
         BNRM2 = DNRM2( N, B, 1 )
         IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
      ENDIF
*
      RESID = DNRM2( N, R, 1 ) / BNRM2
* WFM 10/18/07 print the convergence history
      IF (PRCONV) THEN
         IF ( INFO.EQ.-1 ) THEN
            write(6,"(A)")
            write(6,"(A)") "L2 norm of linear system residual"
            write(6,"(A)")
            write(6,"(1x,A,A)") "     iteration",
     *      "    absolute residual    relative residual      reduction"
            write(6,"(A)")
            write(6,"(SS,1P,1X,I11,2E18.10E2)") 0,dnrm2(N,R,1),
     *           RESID
            OLDRES = RESID
            ITER = 1
         ELSE
            write(6,"(SS,1P,1X,I11,3E18.10E2)") ITER,dnrm2(N,R,1),
     *           RESID,RESID/OLDRES
            OLDRES = RESID
            ITER = ITER + 1
         ENDIF
      ENDIF
* WFM 10/18/07 end
*
      INFO = 0
      IF ( RESID.LE.TOL )
     $     INFO = 1
*
      RETURN
*
*     End of STOPTEST2
*
      END
