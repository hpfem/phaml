C WFM 6/3/03 Removed sequence numbers in columns 73-80
C            Replaced computed goto with regular if statements
C     6/4/07 Replaced FORTRAN 66 style DIMENSION(1) with assumed size arrays
C
C     ALGORITHM 582, COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.8, NO. 2,
C     JUN., 1982, P. 190.
C     ==============================================================
C
C     GIBBS-POOLE-STOCKMEYER AND GIBBS-KING ALGORITHMS ...
C
      SUBROUTINE   GPSKCA   (N, DEGREE, RSTART, CONNEC, OPTPRO, WRKLEN, 
     1                       PERMUT, WORK, BANDWD, PROFIL, ERROR, SPACE)
C
C     ==================================================================
C     ==================================================================
C     =                                                                =
C     = B A N D W I D T H    OR    P R O F I L E    R E D U C T I O N  =
C     =        FOR A SPARSE AND (STRUCTURALLY) SYMMETRIC MATRIX,       =
C     =                         USING EITHER                           =
C     =                                                                =
C     =   THE GIBBS-POOLE-STOCKMEYER ALGORITHM (BANDWIDTH REDUCTION)   =
C     =                               OR                               =
C     =          THE GIBBS-KING ALGORITHM (PROFILE REDUCTION)          =
C     =                                                                =
C     ==================================================================
C     ==================================================================
C     =     THIS CODE SUPERSEDES TOMS ALGORITHMS 508 AND 509 IN THE    =
C     =     COLLECTED ALGORITHMS OF THE ACM (CALGO).                   =
C     ==================================================================
C     ==================================================================
C
C     -------------------
C     P A R A M E T E R S
C     -------------------
C
      INTEGER     N, RSTART(N), WRKLEN, BANDWD, PROFIL, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(*), PERMUT(N), WORK(WRKLEN)
      INTEGER     DEGREE(N), CONNEC(*), PERMUT(N), WORK(WRKLEN)
C
      LOGICAL     OPTPRO
C
C     ------------------------------------------------------------------
C
C     INPUT PARAMETERS:
C     ----- ----------
C
C         N      -- THE DIMENSION OF THE MATRIX
C
C         DEGREE,
C         RSTART,
C         CONNEC -- DESCRIBE THE STRUCTURE OF THE SPARSE MATRIX.
C                   DEGREE(I) SPECIFIES THE NUMBER OF NON-ZERO
C                   OFF-DIAGONAL ENTRIES IN THE  I-TH  ROW OF THE
C                   SPARSE MATRIX.  THE COLUMN INDICES OF THESE
C                   ENTRIES ARE GIVEN IN CONSECUTIVE LOCATIONS IN
C                   CONNEC, STARTING AT LOCATION  RSTART(I).
C                   IN OTHER WORDS, THE INDICES OF THE NON-ZERO
C                   OFF-DIAGONAL ELEMENTS OF THE  I-TH  ROW ARE FOUND
C                   IN:
C                       CONNEC (RSTART(I)),
C                       CONNEC (RSTART(I) + 1),
C                                . . .
C                       CONNEC (RSTART(I) + DEGREE(I) - 1)
C
C                   DIMENSIONS:
C                       RSTART IS DIMENSION  N  (OR LONGER).
C                       DEGREE IS DIMENSION  N  (OR LONGER).
C                       CONNEC IS DIMENSION ROUGHLY THE NUMBER OF NON-
C                              ZERO ENTRIES IN THE MATRIX.
C
C         OPTPRO -- .TRUE. IF REDUCING THE PROFILE OF THE MATRIX
C                          IS MORE IMPORTANT THAN REDUCING THE
C                          BANDWIDTH
C                   .FALSE. IF BANDWIDTH REDUCTION IS MOST IMPORTANT
C
C         WRKLEN -- THE  ACTUAL  LENGTH OF THE VECTOR  WORK  AS SUPPLIED
C                   BY THE USER.  SEE THE DISCUSSION OF THE WORKSPACE
C                   'WORK'  BELOW FOR TYPICAL STORAGE REQUIREMENTS.
C                   THE VALUE OF  WRKLEN  WILL BE USED TO ENSURE THAT
C                   THE ROUTINE WILL NOT USE MORE STORAGE THAN IS
C                   AVAILABLE.  IF NOT ENOUGH SPACE IS GIVEN IN  WORK
C                   TO PERMIT A SOLUTION TO BE FOUND, THE  ERROR  FLAG
C                   WILL BE SET AND FURTHER COMPUTATION STOPPED.
C
C
C     INPUT AND OUTPUT PARAMETER:
C     ----- --- ------ ---------
C
C         PERMUT -- ON INPUT, AN ALTERNATIVE REORDERING FOR THE
C                   ROWS AND COLUMNS OF THE MATRIX.  PERMUT(I) GIVES
C                   THE POSITION IN WHICH ROW AND COLUMN  I  SHOULD
C                   BE PLACED TO REDUCE THE BANDWIDTH OR THE PROFILE.
C                   IF THE USER HAS NO ALTERNATIVE TO THE NATURAL
C                   ORDERING IMPLICIT IN  DEGREE,  RSTART  AND  CONNEC,
C                   HE SHOULD INITIALIZE  PERMUT  TO BE THE IDENTITY
C                   PERMUTATION  PERMUT(I) = I .
C
C                   ON OUTPUT,  PERMUT  WILL CONTAIN THE PERMUTATION
C                   FOR REORDERING THE ROWS AND COLUMNS WHICH REDUCES
C                   THE BANDWIDTH AND/OR PROFILE.  THE RESULT WILL BE
C                   THE REORDERING FOUND BY 'GPSKCA' OR THE REORDERING
C                   GIVEN BY THE USER IN 'PERMUT', WHICHEVER DOES THE
C                   JOB BETTER.
C
C
C     OUTPUT PARAMETERS:
C     ------ ----------
C
C         WORK   -- A TEMPORARY STORAGE VECTOR, OF LENGTH SOMEWHAT
C                   GREATER THAN  3N.  THE SPACE BEYOND  3N  REQUIRED
C                   IS PROBLEM-DEPENDENT.  ANY PROBLEM CAN BE SOLVED
C                   IN  6N+3  LOCATIONS.
C                   MOST PROBLEMS CAN BE REORDERED WITH  4N  LOCATIONS
C                   IN 'WORK'.  IF SPACE IS NOT A CONSTRAINT, PROVIDE
C                   6N+3  LOCATIONS IN 'WORK'.  OTHERWISE, PROVIDE AS
C                   MUCH MORE THAN  3N  AS IS CONVENIENT AND CHECK THE
C                   ERROR FLAG AND SPACE REQUIRED PARAMETERS (SEE BELOW)
C
C                   ON OUTPUT, THE 1ST  N  LOCATIONS OF WORK WILL BE
C                   A LISTING OF THE ORIGINAL ROW AND COLUMN INDICES AS
C                   THEY APPEAR IN THE COMPUTED REORDERING.
C                   LOCATIONS  N+1, ... , 2N  OF  WORK  WILL CONTAIN
C                   THE NEW POSITIONS FOR THE EQUATIONS IN THE ORDER
C                   FOUND BY GPSKCA.  THUS, THE TWO VECTORS ARE INVERSE
C                   PERMUTATIONS OF EACH OTHER.  IF THE ORDERING
C                   FOUND BY THIS ALGORITHM IS BETTER THAN THE USER-
C                   SUPPLIED ORDER, THE SECOND PERMUTATION VECTOR IS
C                   IDENTICAL TO THE RESULT RETURNED IN  'PERMUT'.
C
C         BANDWD -- THE BANDWIDTH OF THE MATRIX WHEN ROWS AND COLUMNS
C                   ARE REORDERED IN THE ORDERING RETURNED IN  PERMUT.
C
C         PROFIL -- THE PROFILE OF THE MATRIX WHEN ROWS AND COLUMNS ARE
C                   REORDERED IN THE ORDERING RETURNED IN  PERMUT.
C
C         ERROR  -- WILL BE EQUAL TO ZERO IF A NEW NUMBERING COULD BE
C                   FOUND IN THE SPACE PROVIDED.  OTHERWISE,  ERROR
C                   WILL BE SET TO A POSITIVE ERROR CODE (SEE TABLE
C                   GIVEN BELOW).  IF THE REORDERING ALGORITHM HAS BEEN
C                   STOPPED BY LACK OF WORKSPACE, THE SPACE PARAMETER
C                   WILL BE SET TO THE NUMBER OF ADDITIONAL LOCATIONS
C                   REQUIRED TO COMPLETE AT LEAST THE NEXT PHASE OF
C                   THE ALGORITHM.
C
C                   WHENEVER A NON-ZERO VALUE FOR  ERROR  IS GIVEN
C                   PERMUT  WILL RETAIN THE VALUES PROVIDED BY THE USER
C                   AND THE SCALARS  BANDWD  AND  PROFIL  WILL BE SET TO
C                   OUTRAGEOUS VALUES.  IT IS THE USER'S RESPONSIBILITY
C                   TO CHECK THE STATUS OF  ERROR.
C
C         SPACE  -- WILL INDICATE EITHER HOW MUCH SPACE THE REORDERING
C                   ACTUALLY REQUIRED OR HOW MUCH SPACE WILL BE
C                   REQUIRED TO COMPLETE THE NEXT PHASE OF THE
C                   REORDERING ALGORITHM.  THE POSSIBLE OUTCOMES ARE ..
C
C                      ERROR = 0          SPACE IS THE MINIMAL VALUE FOR
C                                         WRKLEN  REQUIRED TO REORDER
C                                         THIS MATRIX AGAIN.
C
C                      ERROR <> 0         SPACE IS THE MINIMUM NUMBER
C                      DUE TO LACK OF     OF EXTRA WORKSPACE REQUIRED
C                      WORKSPACE          TO CONTINUE THE REORDERING
C                                         ALGORITHM ON THIS MATRIX.
C
C                      ERROR <> 0         SPACE = -1
C                      DUE TO ERROR
C                      IN DATA STRUCTURES
C
C
C     ==================================================================
C
C     ----------------------
C     E R R O R    C O D E S
C     ----------------------
C
C         ERROR CODES HAVE THE FORM  0XY  OR  1XY.
C
C         ERRORS OF THE FORM  1XY  RESULT FROM INADEQUATE WORKSPACE.
C
C         ERRORS OF THE FORM  0XY  ARE INTERNAL PROGRAM CHECKS, WHICH
C         MOST LIKELY OCCUR BECAUSE THE CONNECTIVITY STRUCTURE OF THE
C         MATRIX IS REPRESENTED INCORRECTLY (E.G., THE DEGREE OF
C         A NODE IS NOT CORRECT  OR  NODE I IS CONNECTED TO NODE J,
C         BUT NOT CONVERSELY).
C
C         THE LAST DIGIT (Y) IS MAINLY USEFUL FOR DEBUGGING THE
C         THE REORDERING ALGORITHM.  THE MIDDLE DIGIT  (X)  INDICATES
C         HOW MUCH OF THE ALGORITHM HAS BEEN PERFORMED.
C         THE TABLE BELOW GIVES THE CORRESPONDENCE BETWEEN THE
C         VALUES OF  X   AND THE STRUCTURE OF THE ALGORITHM.
C             X = 0     INITIAL PROCESSING
C             X = 1     COMPUTING PSEUDO-DIAMETER  (ALGORITHM I)
C             X = 2     TRANSITION BETWEEN ALGORITHM I AND II
C             X = 3     COMBINING LEVEL STRUCTURES (ALGORITHM II)
C             X = 4     TRANSITION BETWEEN ALGORITHM II AND III
C             X = 5     BANDWIDTH NUMBERING (ALGORITHM IIIA)
C             X = 6     PROFILE NUMBERING (ALGORITHM IIIB)
C             X = 7     FINAL BANDWIDTH/PROFILE COMPUTATION
C
C     ==================================================================
C
C     ---------------------    ---------------
C     A L T E R N A T I V E    V E R S I O N S
C     ---------------------    ---------------
C
C     SHORT INTEGER VERSION
C
C         ON MACHINES WITH TWO OR MORE PRECISIONS FOR INTEGERS,
C         ALL OF THE INPUT ARRAYS EXCEPT 'RSTART' CAN BE CONVERTED
C         TO THE SHORTER PRECISION AS LONG AS THAT SHORTER PRECISION
C         ALLOWS NUMBERS AS LARGE AS 'N'.  A VERSION OF THIS CODE
C         SUITABLE FOR USE ON IBM COMPUTERS (INTEGER * 2) IS EMBEDDED
C         AS COMMENTS IN THIS CODE.  ALL SUCH COMMENTS HAVE THE
C         CHARACTERS 'CIBM' IN THE FIRST FOUR COLUMNS, AND PRECEDE THE
C         EQUIVALENT STANDARD CODE WHICH THEY WOULD REPLACE.
C
C     CONNECTIVITY COMPATIBILITY VERSION
C
C         THE ORIGINAL (1976) TOMS CODE  'REDUCE'  USED A LESS STORAGE
C         EFFICIENT FORMAT FOR THE CONNECTIVITY TABLE  'CONNEC'.
C         THE 1976 CODE USED A RECTANGULAR MATRIX OF DIMENSIONS
C         N  BY  MAXDGR,  WHERE  MAXDGR  IS AT LEAST AS LARGE AS
C         THE MAXIMUM DEGREE OF ANY NODE IN THE GRAPH OF THE MATRIX.
C         THE FORMAT USED IN THE CURRENT CODE IS OFTEN SUBSTANTIALLY
C         MORE EFFICIENT.  HOWEVER, FOR USERS FOR WHOM CONVERSION WILL
C         BE DIFFICULT OR IMPOSSIBLE, TWO ALTERNATIVES ARE ..
C             1.  SIMPLY NOTE THAT CHANGING THE ORDER OF SUBSCRIPTS
C                 IN A RECTANGULAR CONNECTION TABLE WILL ENABLE YOU
C                 TO USE THE NEW VERSION.  THIS SUBROUTINE WILL ACCEPT A
C                 RECTANGULAR CONNECTION TABLE OF DIMENSIONS
C                     MAXDGR BY N,
C                 PROVIDED THAT  RSTART(I)  IS SET TO  (I-1)*MAXDGR + 1.
C             2.  THE AUTHOR WILL MAKE AVAILABLE A VARIANT VERSION
C                 'GPSKRA', WHICH EXPECTS THE ADJACENCY MATRIX OR
C                 CONNECTIVITY TABLE IN THE SAME FORM AS DID  'REDUCE'.
C                 THIS VERSION CAN BE OBTAINED BY WRITING TO ..
C                     JOHN GREGG LEWIS
C                     BOEING COMPUTER SERVICES COMPANY
C                     MAIL STOP 9C-01
C                     P.O. BOX 24346
C                     SEATTLE, WA 98124
C                 PLEASE INCLUDE A DESCRIPTION OF THE COMPUTING
C                 ENVIRONMENT ON WHICH YOU WILL BE USING THE CODE.
C
C     ==================================================================
C
      INTEGER     I, INC1, INC2, AVAIL, NXTNUM, LOWDG, STNODE, NLEFT,
     1            TREE1, TREE2, DEPTH, EMPTY, STOTAL, REQD, CSPACE,
     2            LVLLST, LVLPTR, ACTIVE, RVNODE, WIDTH1, WIDTH2, MXDG
C
      LOGICAL     REVRS1, ONEIS1
C
C     ==================================================================
C
C     << NUMBER ANY DEGREE ZERO NODES >>
C
C     WHILE << SOME NODES YET UNNUMBERED >> DO
C         << FIND A PSEUDO-DIAMETER OF THE MATRIX GRAPH >>
C         << CONVERT FORM OF LEVEL TREES >>
C         << COMBINE LEVEL TREES INTO ONE LEVEL STRUCTURE >>
C         << CONVERT FORM OF LEVEL STRUCTURE >>
C         IF OPTPRO THEN
C             << RENUMBER BY KING ALGORITHM >>
C         ELSE
C             << RENUMBER BY REVERSE CUTHILL-MCKEE ALGORITHM >>
C
C     ==================================================================
C
C     ... INITIALIZE COUNTERS, THEN NUMBER ANY NODES OF DEGREE  0.
C         THE LIST OF NODES, BY NEW NUMBER, WILL BE BUILT IN PLACE AT
C         THE FRONT OF THE WORK AREA.
C
      NXTNUM = 1
      ERROR = 0
      SPACE = 2*N
C
      MXDG = 0
      DO 300 I = 1, N
          IF  (DEGREE(I) .LT. 0) GO TO 6000
          IF  (DEGREE(I) .EQ. 0) GO TO 100
          IF  (DEGREE(I) .GT. 0) GO TO 200
  100         WORK(NXTNUM) = I
              NXTNUM = NXTNUM + 1
              GO TO 300
  200         IF  (DEGREE(I) .GT. MXDG)  MXDG = DEGREE(I)
  300 CONTINUE
C
C
C     ==============================
C     ... WHILE  NXTNUM <= N  DO ...
C     ==============================
C
 1000 IF  ( NXTNUM .GT. N )  GO TO 2000
C
C         ... FIND AN UNNUMBERED NODE OF MINIMAL DEGREE
C
          LOWDG = MXDG + 1
          STNODE = 0
          DO 400 I = 1, N
              IF ( (DEGREE(I) .LE. 0) .OR. (DEGREE(I) .GE. LOWDG) )
     1           GO TO 400
                  LOWDG = DEGREE(I)
                  STNODE = I
  400     CONTINUE
C
          IF ( STNODE .EQ. 0 )  GO TO 6100
C
C         ... SET UP POINTERS FOR THREE LISTS IN WORK AREA, THEN LOOK
C             FOR PSEUDO-DIAMETER, BEGINNING WITH STNODE.
C
          AVAIL = (WRKLEN - NXTNUM + 1) / 3
          NLEFT = N - NXTNUM + 1
          SPACE = MAX0 (SPACE, NXTNUM + 3*N - 1)
          IF ( AVAIL .LT. N )  GO TO 5200
C
          CALL GPSKCB (N, DEGREE, RSTART, CONNEC, AVAIL, NLEFT,
     1                 STNODE, RVNODE, WORK(NXTNUM), TREE1, TREE2,
     2                 ACTIVE, DEPTH, WIDTH1, WIDTH2,
     3                 ERROR, SPACE)
          IF ( ERROR .NE. 0 )  GO TO 5000
          SPACE = MAX0 (SPACE, NXTNUM + 3*(ACTIVE+DEPTH+1) - 1)
C
C         ... DYNAMIC SPACE CHECK FOR MOST OF REMAINDER OF ALGORITHM
C
          REQD = MAX0 (NXTNUM + 2*N + 3*DEPTH - 1, 3*N + 2*DEPTH + 1)
          SPACE = MAX0 (SPACE, REQD)
          IF  ( WRKLEN .LT. REQD )  GO TO 5300
C
C
C         ... OUTPUT FROM GPSKCB IS A PAIR OF LEVEL TREES, IN THE FORM
C             OF LISTS OF NODES BY LEVEL.  CONVERT THIS TO TWO LISTS OF
C             OF LEVEL NUMBER BY NODE.  AT THE SAME TIME PACK
C             STORAGE SO THAT ONE OF THE LEVEL TREE VECTORS IS AT THE
C             BACK END OF THE WORK AREA.
C
          LVLPTR = NXTNUM + AVAIL - DEPTH
          CALL GPSKCE (N, AVAIL, ACTIVE, DEPTH, WRKLEN, WORK(NXTNUM),
     1                 WORK(LVLPTR), WORK(1), NXTNUM, TREE1,
     2                 TREE2, WIDTH1, WIDTH2, ONEIS1, ERROR, SPACE)
          IF ( ERROR .NE. 0 ) GO TO 5000
          IF (( TREE1 .NE. WRKLEN - N + 1 ) .OR. (TREE2 .NE. NXTNUM))
     1      GO TO 6200
C
C         ... COMBINE THE TWO LEVEL TREES INTO A MORE GENERAL
C             LEVEL STRUCTURE.
C
          AVAIL = WRKLEN - NXTNUM + 1 - 2*N - 3*DEPTH
          STOTAL = N + NXTNUM
          EMPTY = STOTAL + DEPTH
          INC1 = TREE1 - DEPTH
          INC2 = INC1 - DEPTH
C
          CALL GPSKCG (N, DEGREE, RSTART, CONNEC, ACTIVE, WIDTH1,
     1                 WIDTH2, WORK(TREE1), WORK(TREE2), WORK(EMPTY),
     2                 AVAIL, DEPTH, WORK(INC1), WORK(INC2),
     3                 WORK(STOTAL), ONEIS1, REVRS1, ERROR, CSPACE)
C
          IF ( ERROR .NE. 0 )  GO TO 5000
          SPACE = MAX0 (SPACE, NXTNUM + CSPACE - 1)
C
C         ... COMBINED LEVEL STRUCTURE IS REPRESENTED BY GPSKCG AS
C             A VECTOR OF LEVEL NUMBERS.  FOR RENUMBERING PHASE,
C             CONVERT THIS ALSO TO THE INVERSE PERMUTATION.
C
          LVLPTR = TREE1 - (DEPTH + 1)
          LVLLST = LVLPTR - ACTIVE
          IF ( STOTAL + DEPTH .GT. LVLPTR )  GO TO 6300
C
          CALL GPSKCI (N, ACTIVE, DEPTH, WORK(TREE1), WORK(LVLLST),
     1                 WORK(LVLPTR), WORK(STOTAL), ERROR, SPACE)
          IF  (ERROR .NE. 0)  GO TO 5000
C
C         ... NOW RENUMBER ALL MEMBERS OF THIS COMPONENT USING
C             EITHER A REVERSE CUTHILL-MCKEE OR A KING STRATEGY,
C             AS PROFILE OR BANDWIDTH REDUCTION IS MORE IMPORTANT.
C
          IF ( OPTPRO )  GO TO 500
              CALL GPSKCJ (N, DEGREE, RSTART, CONNEC, ACTIVE,
     1                     WORK(NXTNUM), STNODE, RVNODE, REVRS1, DEPTH,
     2                     WORK(LVLLST), WORK(LVLPTR), WORK(TREE1),
     3                     ERROR, SPACE)
              IF ( ERROR .NE. 0 )  GO TO 5000
              NXTNUM = NXTNUM + ACTIVE
              GO TO 600
C
  500         CALL GPSKCK (N, DEGREE, RSTART, CONNEC, LVLLST-1, NXTNUM,
     1                     WORK, ACTIVE, DEPTH, WORK(LVLLST),
     2                     WORK(LVLPTR), WORK(TREE1), ERROR, SPACE)
              IF ( ERROR .NE. 0 )  GO TO 5000
C
C         =========================================================
C         ... END OF WHILE LOOP ... REPEAT IF GRAPH IS DISCONNECTED
C         =========================================================
C
  600     GO TO 1000
C
C     ... CHECK WHETHER INITIAL NUMBERING OR FINAL NUMBERING
C         PROVIDES BETTER RESULTS
C
 2000 IF  (WRKLEN .LT. 2*N)  GO TO 5400
C
      IF  (OPTPRO)  GO TO 2100
          CALL GPSKCL (N, DEGREE, RSTART, CONNEC, WORK(1), WORK(N+1),
     *                 PERMUT, BANDWD, PROFIL, ERROR, SPACE)
          GO TO 2200
*
 2100     CALL GPSKCM (N, DEGREE, RSTART, CONNEC, WORK(1), WORK(N+1),
     1                 PERMUT, BANDWD, PROFIL, ERROR, SPACE)
C
 2200 RETURN
C
C
C     . . .  E R R O R   D I A G N O S T I C S
C            ---------------------------------
C
C     ... ERROR DETECTED BY LOWER LEVEL ROUTINE.  MAKE SURE THAT SIGNS
C         OF DEGREE ARE PROPERLY SET
C
 5000 DO 5100 I = 1, N
          IF  (DEGREE(I) .LT. 0)  DEGREE(I) = -DEGREE(I)
 5100 CONTINUE
C
      BANDWD = -1
      PROFIL = -1
      RETURN
C
C     ... STORAGE ALLOCATION ERRORS DETECTED IN THIS ROUTINE
C
 5200 ERROR = 101
      SPACE = -1
      GO TO 5000
C
 5300 ERROR = 102
      SPACE = -1
      GO TO 5000
C
 5400 ERROR =  10
      SPACE = 2*N - WRKLEN
      GO TO 5000
C
C     ... DATA STRUCTURE ERRORS DETECTED IN THIS ROUTINE
C
 6000 ERROR = 1
      GO TO 6900
C
 6100 ERROR = 2
      GO TO 6900
C
 6200 ERROR = 3
      GO TO 6900
C
 6300 ERROR = 4
C
 6900 SPACE = -1
      GO TO 5000
      END
      SUBROUTINE  GPSKCB  (N, DEGREE, RSTART, CONNEC, AVAIL, NLEFT,
     1                     STNODE, RVNODE, WORK, FORWD, BESTBK, NNODES,
     2                     DEPTH, FWIDTH, BWIDTH, ERROR, SPACE)
C
C     ==================================================================
C
C     FIND A PSEUDO-DIAMETER OF THE MATRIX GRAPH ...
C
C         << BUILD A LEVEL TREE FROM STNODE >>
C         REPEAT
C             << BUILD A LEVEL TREE FROM EACH NODE 'BKNODE' IN THE
C                DEEPEST LEVEL OF  STNODE'S TREE >>
C             << REPLACE 'STNODE' WITH 'BKNODE' IF A DEEPER AND
C                NARROWER TREE WAS FOUND. >>
C         UNTIL
C             << NO FURTHER IMPROVEMENT MADE >>
C
C     ... HEURISTIC ABOVE DIFFERS FROM THE ALGORITHM PUBLISHED IN
C         SIAM J. NUMERICAL ANALYSIS, BUT MATCHES THE CODE
C         DISTRIBUTED BY TOMS.
C
C
C     PARAMETERS :
C
C         N, DEGREE, RSTART & CONNEC  DESCRIBE THE MATRIX STRUCTURE
C
C         WORK   -- WORKING SPACE, OF LENGTH  3*AVAIL, USED TO STORE
C         THREE LEVEL TREES.
C
C         STNODE IS INITIALLY THE NUMBER OF A NODE TO BE USED TO
C             START THE PROCESS, TO BE THE ROOT OF THE FIRST TREE.
C             ON OUTPUT, STNODE IS THE END OF THE PSEUDO-DIAMETER WHOSE
C             LEVEL TREE IS NARROWEST.
C
C         RVNODE WILL BE THE OTHER END OF THE PSEUDO-DIAMETER.
C
C         NNODES WILL BE THE NUMBER OF NODES IN THIS CONNECTED
C             COMPONNENT OF THE MATRIX GRAPH, I.E., THE LENGTH OF
C             THE LEVEL TREES.
C
C         DEPTH  -- THE DEPTH OF THE LEVEL TREES BEING RETURNED,
C                   I.E., THE LENGTH OF THE PSEUDO-DIAMETER.
C
C     ==================================================================
C
C     STRUCTURE OF WORKSPACE ...
C
C     ---------------------------------------------------------------
C     : NUMBERED :  TLIST1  PTR1  :  TLIST2  PTR2  :  TLIST3  PTR3  :
C     ---------------------------------------------------------------
C
C     TLISTI IS A LIST OF NODES OF LENGTH  'ACTIVE'
C     PTRI   IS A LIST OF POINTERS INTO TLISTI, OF LENGTH  'DEPTH+1'
C
C     ==================================================================
C
      INTEGER     N, RSTART(N), AVAIL, NLEFT,
     1            STNODE, RVNODE, FORWD, BESTBK, NNODES, DEPTH, FWIDTH,
     4            BWIDTH, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(*), WORK(AVAIL,3)
      INTEGER     DEGREE(N), CONNEC(*), WORK(AVAIL,3)
C
C     ----------------
C
      INTEGER     BACKWD, MXDPTH, WIDTH, FDEPTH, LSTLVL,
     1            NLAST, T, I, BKNODE, LSTLVI
C
      LOGICAL     IMPROV
C
C
C     ... BUILD INITIAL LEVEL TREE FROM 'STNODE'.  FIND OUT HOW MANY
C         NODES LIE IN THE CURRENT CONNECTED COMPONENT.
C
      FORWD = 1
      BACKWD = 2
      BESTBK = 3
C
      CALL GPSKCC (N, DEGREE, RSTART, CONNEC, STNODE, AVAIL, NLEFT,
     1              WORK(1,FORWD), NNODES, DEPTH, WIDTH, ERROR,
     2              SPACE)
      IF ( ERROR .NE. 0 )  GO TO 5000
C
      MXDPTH = AVAIL - NNODES - 1
C
C     ==========================================
C     REPEAT UNTIL NO DEEPER TREES ARE FOUND ...
C     ==========================================
C
 1000     FWIDTH = WIDTH
          FDEPTH = DEPTH
          LSTLVL = AVAIL - DEPTH + 1
          NLAST = WORK (LSTLVL-1, FORWD) - WORK (LSTLVL, FORWD)
          LSTLVL = WORK (LSTLVL, FORWD)
          BWIDTH = N+1
C
C         ... SORT THE DEEPEST LEVEL OF 'FORWD' TREE INTO INCREASING
C             ORDER OF NODE DEGREE.
C
          CALL GPSKCQ (NLAST, WORK(LSTLVL,FORWD), N, DEGREE, ERROR)
          IF  (ERROR .NE. 0)  GO TO 6000
C
C         ... BUILD LEVEL TREE FROM NODES IN 'LSTLVL' UNTIL A DEEPER
C             AND NARROWER TREE IS FOUND OR THE LIST IS EXHAUSTED.
C
          IMPROV = .FALSE.
          DO 1200 I = 1, NLAST
              LSTLVI = LSTLVL + I - 1
              BKNODE = WORK (LSTLVI, FORWD)
              CALL GPSKCD (N, DEGREE, RSTART, CONNEC, BKNODE, AVAIL,
     1                     NNODES, MXDPTH, WORK(1,BACKWD), DEPTH, WIDTH,
     2                     BWIDTH, ERROR, SPACE)
              IF ( ERROR .NE. 0 )  GO TO 5000
C
              IF ( DEPTH .LE. FDEPTH )  GO TO 1100
C
C                 ... NEW DEEPER TREE ... MAKE IT NEW 'FORWD' TREE
C                     AND BREAK OUT OF 'DO' LOOP.
C
                  IMPROV = .TRUE.
                  T = FORWD
                  FORWD = BACKWD
                  BACKWD = T
                  STNODE = BKNODE
                  GO TO 1300
C
C                 ... ELSE CHECK FOR NARROWER TREE.
C
 1100             IF ( WIDTH .GE. BWIDTH )  GO TO 1200
                      T = BESTBK
                      BESTBK = BACKWD
                      BACKWD = T
                      BWIDTH = WIDTH
                      RVNODE = BKNODE
 1200     CONTINUE
C
C         ... END OF REPEAT LOOP
C         ----------------------
C
 1300     IF ( IMPROV )  GO TO 1000
C
      DEPTH = FDEPTH
      RETURN
C
C     ... IN CASE OF ERROR, SIMPLY RETURN ERROR FLAG TO USER.
C
 5000 RETURN
C
 6000 ERROR = 11
      SPACE = -1
      RETURN
C
      END
      SUBROUTINE   GPSKCC   (N, DEGREE, RSTART, CONNEC, STNODE, AVAIL,
     1                       NLEFT, LIST, ACTIVE, DEPTH, WIDTH, ERROR,
     2                       SPACE)
C
C     ==================================================================
C     BUILD THE LEVEL TREE ROOTED AT 'STNODE' IN THE SPACE PROVIDED IN
C     LIST.  CHECK FOR OVERRUN OF SPACE ALLOCATION.
C     ==================================================================
C
      INTEGER     N, RSTART(N), STNODE, AVAIL, NLEFT,
     1            ACTIVE, DEPTH, WIDTH, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(*), LIST(AVAIL)
      INTEGER     DEGREE(N), CONNEC(*), LIST(AVAIL)
C
C     ... PARAMETERS:
C
C         INPUT ...
C
C             N, DEGREE, RSTART, CONNEC -- DESCRIBE THE MATRIX STRUCTURE
C
C             STNODE -- THE ROOT OF THE LEVEL TREE.
C
C             AVAIL  -- THE LENGTH OF THE WORKING SPACE AVAILABLE
C
C             NLEFT  -- THE NUMBER OF NODES YET TO BE NUMBERED
C
C             LIST   -- THE WORKING SPACE.
C
C         OUTPUT ...
C
C             ACTIVE -- THE NUMBER OF NODES IN THE COMPONENT
C
C             DEPTH  -- THE DEPTH OF THE LEVEL TREE ROOTED AT  STNODE.
C
C             WIDTH  -- THE WIDTH OF THE LEVEL TREE ROOTED AT  STNODE.
C
C             ERROR  -- ZERO UNLESS STORAGE WAS INSUFFICIENT.
C
C     ------------------------------------------------------------------
C
      INTEGER         LSTART, NLEVEL, FRONT, J, NEWNOD, PTR, CDGREE,
     1                LFRONT, LISTJ
C
C     ... BUILD THE LEVEL TREE USING  LIST  AS A QUEUE AND LEAVING
C         THE NODES IN PLACE.  THIS GENERATES THE NODES ORDERED BY LEVEL
C         PUT POINTERS TO THE BEGINNING OF EACH LEVEL, BUILDING FROM
C         THE BACK OF THE WORK AREA.
C
      ACTIVE = 1
      DEPTH = 0
      WIDTH = 0
      ERROR = 0
      LSTART = 1
      FRONT = 1
      LIST (ACTIVE) = STNODE
      DEGREE (STNODE) = -DEGREE (STNODE)
      LIST (AVAIL)  = 1
      NLEVEL = AVAIL
C
C     ... REPEAT UNTIL QUEUE BECOMES EMPTY OR WE RUN OUT OF SPACE.
C     ------------------------------------------------------------
C
 1000     IF ( FRONT .LT. LSTART ) GO TO 1100
C
C         ... FIRST NODE OF LEVEL.  UPDATE POINTERS.
C
              LSTART = ACTIVE + 1
              WIDTH = MAX0 (WIDTH, LSTART - LIST(NLEVEL))
              NLEVEL = NLEVEL - 1
              DEPTH = DEPTH + 1
              IF ( NLEVEL .LE. ACTIVE )  GO TO 5000
                  LIST (NLEVEL) = LSTART
C
C         ... FIND ALL NEIGHBORS OF CURRENT NODE, ADD THEM TO QUEUE.
C
 1100     LFRONT = LIST (FRONT)
          PTR = RSTART (LFRONT)
          CDGREE = -DEGREE (LFRONT)
          IF (CDGREE .LE. 0)  GO TO 6000
          DO 1200 J = 1, CDGREE
              NEWNOD = CONNEC (PTR)
              PTR = PTR + 1
C
C             ... ADD TO QUEUE ONLY NODES NOT ALREADY IN QUEUE
C
              IF ( DEGREE(NEWNOD) .LE. 0 )  GO TO 1200
                  DEGREE (NEWNOD) = -DEGREE (NEWNOD)
                  ACTIVE = ACTIVE + 1
                  IF ( NLEVEL .LE. ACTIVE )  GO TO 5000
                  IF ( ACTIVE .GT. NLEFT  )  GO TO 6000
                      LIST (ACTIVE) = NEWNOD
 1200     CONTINUE
          FRONT = FRONT + 1
C
C         ... IS QUEUE EMPTY?
C         -------------------
C
          IF ( FRONT .LE. ACTIVE )  GO TO 1000
C
C     ... YES, THE TREE IS BUILT.  UNDO OUR MARKINGS.
C
      DO 1300 J = 1, ACTIVE
          LISTJ = LIST(J)
          DEGREE (LISTJ) = -DEGREE (LISTJ)
 1300 CONTINUE
C
      RETURN
C
C     ... INSUFFICIENT STORAGE ...
C
 5000 SPACE = 3 * ( (NLEFT+1-ACTIVE)*DEPTH / NLEFT + (NLEFT+1-ACTIVE) )
      ERROR = 110
      RETURN
C
 6000 ERROR = 12
      SPACE = -1
      RETURN
C
      END
      SUBROUTINE   GPSKCD   (N, DEGREE, RSTART, CONNEC, STNODE, AVAIL,
     1                       ACTIVE, MXDPTH, LIST, DEPTH, WIDTH, MAXWID,
     2                       ERROR, SPACE)
C
C     ==================================================================
C     BUILD THE LEVEL TREE ROOTED AT 'STNODE' IN THE SPACE PROVIDED IN
C     LIST.  OVERFLOW CHECK NEEDED ONLY ON DEPTH OF TREE.
C
C     BUILD THE LEVEL TREE TO COMPLETION ONLY IF THE WIDTH OF ALL
C     LEVELS IS SMALLER THAN 'MAXWID'.  IF A WIDER LEVEL IS FOUND
C     TERMINATE THE CONSTRUCTION.
C     ==================================================================
C
      INTEGER     N, RSTART(N), STNODE, AVAIL, ACTIVE, MXDPTH,
     1            DEPTH, WIDTH, MAXWID, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(*), LIST(AVAIL)
      INTEGER     DEGREE(N), CONNEC(*), LIST(AVAIL)
C
C     ... PARAMETERS:
C
C         INPUT ...
C
C             N, DEGREE, RSTART, CONNEC -- DESCRIBE THE MATRIX STRUCTURE
C
C             STNODE -- THE ROOT OF THE LEVEL TREE.
C
C             AVAIL  -- THE LENGTH OF THE WORKING SPACE AVAILABLE
C
C             NLEFT  -- THE NUMBER OF NODES YET TO BE NUMBERED
C
C             ACTIVE -- THE NUMBER OF NODES IN THE COMPONENT
C
C             MXDPTH -- MAXIMUM DEPTH OF LEVEL TREE POSSIBLE IN
C                       ALLOTTED WORKING SPACE
C
C             LIST   -- THE WORKING SPACE.
C
C         OUTPUT ...
C
C             DEPTH  -- THE DEPTH OF THE LEVEL TREE ROOTED AT  STNODE.
C
C             WIDTH  -- THE WIDTH OF THE LEVEL TREE ROOTED AT  STNODE.
C
C             MAXWID -- LIMIT ON WIDTH OF THE TREE.  TREE WILL NOT BE
C                       USED IF WIDTH OF ANY LEVEL IS AS GREAT AS
C                       MAXWID, SO CONSTRUCTION OF TREE NEED NOT
C                       CONTINUE IF ANY LEVEL THAT WIDE IS FOUND.
C             ERROR  -- ZERO UNLESS STORAGE WAS INSUFFICIENT.
C
C     ------------------------------------------------------------------
C
      INTEGER     LSTART, NLEVEL, FRONT, J, NEWNOD, PTR, BACK,
     1            SPTR, FPTR, LFRONT, LISTJ
C
C     ... BUILD THE LEVEL TREE USING  LIST  AS A QUEUE AND LEAVING
C         THE NODES IN PLACE.  THIS GENERATES THE NODES ORDERED BY LEVEL
C         PUT POINTERS TO THE BEGINNING OF EACH LEVEL, BUILDING FROM
C         THE BACK OF THE WORK AREA.
C
      BACK = 1
      DEPTH = 0
      WIDTH = 0
      ERROR = 0
      LSTART = 1
      FRONT = 1
      LIST (BACK) = STNODE
      DEGREE (STNODE) = -DEGREE (STNODE)
      LIST (AVAIL)  = 1
      NLEVEL = AVAIL
C
C     ... REPEAT UNTIL QUEUE BECOMES EMPTY OR WE RUN OUT OF SPACE.
C     ------------------------------------------------------------
C
 1000     IF ( FRONT .LT. LSTART ) GO TO 1100
C
C         ... FIRST NODE OF LEVEL.  UPDATE POINTERS.
C
              LSTART = BACK + 1
              WIDTH = MAX0 (WIDTH, LSTART - LIST(NLEVEL))
              IF  ( WIDTH .GE. MAXWID )  GO TO 2000
              NLEVEL = NLEVEL - 1
              DEPTH = DEPTH + 1
              IF ( DEPTH .GT. MXDPTH )  GO TO 5000
                  LIST (NLEVEL) = LSTART
C
C         ... FIND ALL NEIGHBORS OF CURRENT NODE, ADD THEM TO QUEUE.
C
 1100     LFRONT = LIST (FRONT)
          SPTR = RSTART (LFRONT)
          FPTR = SPTR - DEGREE (LFRONT) - 1
          DO 1200 PTR = SPTR, FPTR
              NEWNOD = CONNEC (PTR)
C
C             ... ADD TO QUEUE ONLY NODES NOT ALREADY IN QUEUE
C
              IF ( DEGREE(NEWNOD) .LE. 0 )  GO TO 1200
                  DEGREE (NEWNOD) = -DEGREE (NEWNOD)
                  BACK = BACK + 1
                  LIST (BACK) = NEWNOD
 1200     CONTINUE
          FRONT = FRONT + 1
C
C         ... IS QUEUE EMPTY?
C         -------------------
C
          IF ( FRONT .LE. BACK )  GO TO 1000
C
C     ... YES, THE TREE IS BUILT.  UNDO OUR MARKINGS.
C
      IF (BACK .NE. ACTIVE)  GO TO 6000
C
 1300 DO 1400 J = 1, BACK
          LISTJ = LIST(J)
          DEGREE (LISTJ) = -DEGREE (LISTJ)
 1400 CONTINUE
C
      RETURN
C
C     ... ABORT GENERATION OF TREE BECAUSE IT IS ALREADY TOO WIDE
C
 2000 WIDTH = N + 1
      DEPTH = 0
      GO TO 1300
C
C     ... INSUFFICIENT STORAGE ...
C
 5000 SPACE = 3 * ( (ACTIVE+1-BACK)*DEPTH / ACTIVE + (ACTIVE+1-BACK) )
      ERROR = 111
      RETURN
C
 6000 ERROR = 13
      SPACE = -1
      RETURN
C
      END
      SUBROUTINE   GPSKCE   (N, AVAIL, ACTIVE, DEPTH, WRKLEN,
     1                       LVLLST, LVLPTR, WORK, NXTNUM, TREE1, TREE2,
     2                       WIDTH1, WIDTH2, ONEIS1, ERROR, SPACE)
C
C     ==================================================================
C
C     TRANSITION BETWEEN ALGORITHM I AND ALGORITHM II OF
C     THE GIBBS-POOLE-STOCKMEYER PAPER.
C
C     IN THIS IMPLEMENTATION ALGORITHM I REPRESENTS LEVEL TREES AS
C     LISTS OF NODES ORDERED BY LEVEL.  ALGORITHM II APPEARS TO REQUIRE
C     LEVEL NUMBERS INDEXED BY NODE -- VECTORS FOR EFFICIENCY.
C     THIS SUBROUTINE CHANGES THE LEVEL TREE REPRESENTATION TO THAT
C     REQUIRED BY ALGORITHM II.  NOTE THAT THE FIRST ALGORITHM CAN BE
C     CARRIED OUT WITH THE LEVEL NUMBER VECTOR FORMAT, PROBABLY REQURING
C     MORE COMPUTATION TIME, BUT PERHAPS LESS STORAGE.
C
C     INPUT:  TWO LEVEL TREES, AS LEVEL LISTS AND LEVEL POINTERS,
C             FOUND IN TWO OF THE THREE COLUMNS OF THE ARRAYS 'LVLLST'
C             AND 'LVLPTR'
C
C     OUTPUT: TWO LEVEL TREES, AS VECTORS OF LEVEL NUMBERS,
C             ONE PACKED TO THE FRONT, ONE TO THE REAR OF THE WORKING
C             AREA 'WORK'.  NOTE THAT 'WORK', 'LVLLST' AND 'LVLPTR'
C             SHARE COMMON LOCATIONS.
C
C     ================================================================
C
C     ... STRUCTURE OF WORKSPACE
C
C         INPUT .. (OUTPUT FROM GPSKCB)
C
C     --------------------------------------------------------------
C     : NUMBERED : TLIST1  PTR1  :  TLIST2  PTR2  :  TLIST3  PTR3  :
C     --------------------------------------------------------------
C
C         OUTPUT .. (GOES TO COMBIN)
C
C     --------------------------------------------------------------
C     : NUMBERED :  TREE2  :           ...               :  TREE1  :
C     --------------------------------------------------------------
C
C     ==================================================================
C
      INTEGER     N, AVAIL, ACTIVE, DEPTH, WRKLEN, NXTNUM,
     1            WIDTH1, WIDTH2, TREE1, TREE2, ERROR, SPACE
C
CIBM  INTEGER *2  LVLLST(AVAIL,3), LVLPTR(AVAIL,3), WORK(WRKLEN)
      INTEGER     LVLLST(AVAIL,3), LVLPTR(AVAIL,3), WORK(WRKLEN)
C
      LOGICAL     ONEIS1
C
C     ------------------------------------------------------------------
C
      INTEGER     I, BTREE, FTREE, FWIDTH, BWIDTH
C
C
C     ... CHECK THAT WE HAVE ENOUGH ROOM TO DO THE NECESSARY UNPACKING
C
      IF (3*AVAIL .GT. WRKLEN)  GO TO 6000
      IF (AVAIL .LT. N)  GO TO 5100
C
C     ... INPUT HAS THREE POSSIBLE CASES:
C             LVLLST(*,1) IS EMPTY
C             LVLLST(*,2) IS EMPTY
C             LVLLST(*,3) IS EMPTY
C
      FTREE = TREE1
      BTREE = TREE2
      FWIDTH = WIDTH1
      BWIDTH = WIDTH2
C
      TREE1 = WRKLEN - N + 1
      TREE2 = NXTNUM
C
      IF ( (FTREE .EQ. 1) .OR. (BTREE .EQ. 1) )  GO TO 300
C
C         ... CASE 1:  1ST SLOT IS EMPTY.  UNPACK 3 INTO 1, 2 INTO 3
C
          IF (FTREE .NE. 2)  GO TO 100
              ONEIS1 = .TRUE.
              WIDTH2 = BWIDTH
              WIDTH1 = FWIDTH
              GO TO 200
C
  100         ONEIS1 = .FALSE.
              WIDTH1 = BWIDTH
              WIDTH2 = FWIDTH
C
  200     CALL GPSKCF (N, ACTIVE, DEPTH, LVLLST(1,3), LVLPTR(1,3),
     1                    WORK(TREE2), ONEIS1)
C
          CALL GPSKCF (N, ACTIVE, DEPTH, LVLLST(1,2), LVLPTR(1,2),
     1                    WORK(TREE1), .NOT. ONEIS1)
C
          GO TO 1000
C
C
  300 IF ( (FTREE .EQ. 2) .OR. (BTREE .EQ. 2) )  GO TO 600
C
C         ... CASE 2:  2ND SLOT IS EMPTY.  TO ENABLE COMPLETE
C              REPACKING, MOVE 3 INTO 2, THEN FALL INTO NEXT CASE
C
          DO 400 I = 1, ACTIVE
              LVLLST(I,2) = LVLLST(I,3)
  400     CONTINUE
C
          DO 500 I = 1, DEPTH
              LVLPTR(I,2) = LVLPTR(I,3)
  500     CONTINUE
C
C         ... CASE 3:  SLOT 3 IS EMPTY.  MOVE 1 INTO 3, THEN 2 INTO 1.
C
  600     IF (FTREE .EQ. 1) GO TO 700
              ONEIS1 = .FALSE.
              WIDTH1 = BWIDTH
              WIDTH2 = FWIDTH
              GO TO 800
C
  700         ONEIS1 = .TRUE.
              WIDTH1 = FWIDTH
              WIDTH2 = BWIDTH
C
  800     CALL GPSKCF (N, ACTIVE, DEPTH, LVLLST(1,1), LVLPTR(1,1),
     1                    WORK(TREE1), .NOT. ONEIS1)
C
          CALL GPSKCF (N, ACTIVE, DEPTH, LVLLST(1,2), LVLPTR(1,2),
     1                    WORK(TREE2), ONEIS1)
 1000 RETURN
C
C     ------------------------------------------------------------------
C
 5100 SPACE = 3 * (N - AVAIL)
      ERROR = 120
      RETURN
C
 6000 ERROR = 20
      SPACE = -1
      RETURN
C
      END
      SUBROUTINE  GPSKCF  (N, ACTIVE, DEPTH, LVLLST, LVLPTR, LVLNUM,
     1                     REVERS)
C
C     ==================================================================
C
C     CONVERT LEVEL STRUCTURE REPRESENTATION FROM A LIST OF NODES
C     GROUPED BY LEVEL TO A VECTOR GIVING LEVEL NUMBER FOR EACH NODE.
C
C     LVLLST, LVLPTR -- LIST OF LISTS
C
C     LVLNUM -- OUTPUT VECTOR OF LEVEL NUMBERS
C
C     REVERS -- IF .TRUE., NUMBER LEVEL STRUCTURE FROM BACK END
C               INSTEAD OF FROM FRONT
C
C     ==================================================================
C
      INTEGER     N, ACTIVE, DEPTH
C
CIBM  INTEGER *2  LVLLST(ACTIVE), LVLPTR(DEPTH), LVLNUM(N)
      INTEGER     LVLLST(ACTIVE), LVLPTR(DEPTH), LVLNUM(N)
      LOGICAL     REVERS
C
C     ------------------------------------------------------------------
C
      INTEGER     I, LEVEL, LSTART, LEND, XLEVEL, PLSTRT, LVLLSI
C
      IF  (ACTIVE .EQ. N)  GO TO 200
C
C         ... IF NOT ALL NODES OF GRAPH ARE ACTIVE, MASK OUT THE
C             NODES WHICH ARE NOT ACTIVE
C
          DO 100 I = 1, N
              LVLNUM(I) = 0
  100     CONTINUE
C
c WFM 7/30/2004 In the last time through the following loop when LEVEL = DEPTH,
c               PLSTRT=1 which means LEND = LVLPTR(0)-1 which is an invalid
c               subscript.  Using LEVEL=1,DEPTH-1 caused the algorithm to fail,
c               so I don't know what the fix is.  Avoid using subscript
c               bounds checking when compiling this file.
  200 DO 400 LEVEL = 1, DEPTH
          XLEVEL = LEVEL
          PLSTRT = DEPTH - LEVEL + 1
          IF (REVERS) XLEVEL = PLSTRT
          LSTART = LVLPTR (PLSTRT)
          LEND = LVLPTR (PLSTRT - 1) - 1
C
          DO 300 I = LSTART, LEND
              LVLLSI = LVLLST(I)
              LVLNUM (LVLLSI) = XLEVEL
  300     CONTINUE
  400 CONTINUE
C
      RETURN
      END
      SUBROUTINE   GPSKCG   (N, DEGREE, RSTART, CONNEC, ACTIVE, WIDTH1,
     1                       WIDTH2, TREE1, TREE2, WORK, WRKLEN, DEPTH,
     2                       INC1, INC2, TOTAL, ONEIS1, REVRS1, ERROR,
     3                       SPACE)
C
C     ==================================================================
C
C     COMBINE THE TWO ROOTED LEVEL TREES INTO A SINGLE LEVEL STRUCTURE
C     WHICH MAY HAVE SMALLER WIDTH THAN EITHER OF THE TREES.  THE NEW
C     STRUCTURE IS NOT NECESSARILY A ROOTED STRUCTURE.
C
C     PARAMETERS:
C
C         N, DEGREE, RSTART, CONNEC -- GIVE THE DIMENSION AND STRUCTURE
C                                      OF THE SPARSE SYMMETRIC MATRIX
C
C         ACTIVE -- THE NUMBER OF NODES IN THIS CONNECTED COMPONENT OF
C                   THE MATRIX GRAPH
C
C         TREE1  -- ON INPUT, ONE OF THE INPUT LEVEL TREES.  ON
C                   OUTPUT, THE COMBINED LEVEL STRUCTURE
C
C         TREE2  -- THE SECOND INPUT LEVEL TREE
C
C         WIDTH1 -- THE MAXIMUM WIDTH OF A LEVEL IN TREE1
C
C         WIDTH2 -- THE MAXIMUM WIDTH OF A LEVEL IN TREE2
C
C         WORK   -- A WORKING AREA OF LENGTH 'WRKLEN'
C
C         INC1,  -- VECTORS OF LENGTH 'DEPTH'
C         INC2,
C         TOTAL
C
C         ONEIS1 -- INDICATES WHETHER TREE1 OR TREE2 REPRESENTS THE
C                   FORWARD TREE OR THE BACKWARDS TREE OF PHASE 1.
C                   USED TO MIMIC ARBITRARY TIE-BREAKING PROCEDURE OF
C                   ORIGINAL GIBBS-POOLE-STOCKMEYER CODE.
C
C         REVRS1 -- OUTPUT PARAMETER INDICATING WHETHER A BACKWARDS
C                   ORDERING WAS USED FOR THE LARGEST COMPONENT OF
C                   THE REDUCED GRAPH
C
C         ERROR  -- NON-ZERO ONLY IF FAILURE OF SPACE ALLOCATION OR
C                   DATA STRUCTURE ERROR FOUND
C
C         SPACE -- MINIMUM SPACE REQUIRED TO RERUN OR COMPLETE PHASE.
C
C     ------------------------------------------------------------------
C
      INTEGER     N, RSTART(N), ACTIVE, WIDTH1, WIDTH2, WRKLEN, DEPTH,
     2            ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(*), TREE1(N), TREE2(N),
      INTEGER     DEGREE(N), CONNEC(*), TREE1(N), TREE2(N),
     1            WORK(WRKLEN), INC1(DEPTH), INC2(DEPTH), TOTAL(DEPTH)
C
      LOGICAL     ONEIS1, REVRS1
C
C     ==================================================================
C
C     << REMOVE ALL NODES OF PSEUDO-DIAMETERS >>
C     << FIND CONNECTED COMPONENTS OF REDUCED GRAPH >>
C     << COMBINE LEVEL TREES, COMPONENT BY COMPONENT >>
C
C     ==================================================================
C
C     STRUCTURE OF WORKSPACE ...
C
C     ------------------------------------------------------------------
C     : NUMBERED : TREE2 : TOTAL : NODES : START : SIZE : INC1 : INC2 :
C     ------------------------------------------------------------------
C
C     --------
C      TREE1 :
C     --------
C
C         NUMBERED  IS THE SET OF  NUMBERED NODES (PROBABLY EMPTY)
C
C         TREE1 AND TREE1 ARE LEVEL TREES (LENGTH N)
C         TOTAL, INC1 AND INC2  ARE VECTORS OF NODE COUNTS PER LEVEL
C             (LENGTH 'DEPTH')
C         NODES IS THE SET OF NODES IN THE REDUCED GRAPH (THE NODES
C             NOT ON ANY SHORTEST PATH FROM ONE END OF THE
C             PSEUDODIAMETER TO THE OTHER)
C         START, SIZE ARE POINTERS INTO 'NODES', ONE OF EACH FOR
C         EACH CONNECTED COMPONENT OF THE REDUCED GRAPH.
C         THE SIZES OF NODES, START AND SIZE ARE NOT KNOWN APRIORI.
C
C     ==================================================================
      INTEGER     I, SIZE, AVAIL, CSTOP, START, COMPON, TREE1I, PCSTRT,
     1            CSTART, MXINC1, MXINC2, COMPNS, MXCOMP, OFFDIA,
     2            CSIZE, PCSIZE, WORKI, TWORKI
C
C     ------------------------------------------------------------------
C
C     ... FIND ALL SHORTEST PATHS FROM START TO FINISH.  REMOVE NODES ON
C         THESE PATHS AND IN OTHER CONNECTED COMPONENTS OF FULL GRAPH
C         FROM FURTHER CONSIDERATION.  SIGN OF ENTRIES IN TREE1 IS USED
C         AS A MASK.
C
      OFFDIA = ACTIVE
C
      DO 100 I = 1, DEPTH
          TOTAL(I) = 0
  100 CONTINUE
C
      DO 200 I = 1, N
          TREE1I = TREE1 (I)
          IF ((TREE1(I) .NE. TREE2(I)) .OR. (TREE1(I) .EQ. 0)) GO TO 200
              TOTAL (TREE1I) = TOTAL (TREE1I) + 1
              TREE1(I) = - TREE1(I)
              OFFDIA = OFFDIA - 1
  200 CONTINUE
C
      IF ( OFFDIA .EQ. 0 )  GO TO 1100
      IF ( OFFDIA .LT. 0 )  GO TO 6000
C
C     ... FIND CONNECTED COMPONENTS OF GRAPH INDUCED BY THE NODES NOT
C         REMOVED.  'MXCOMP' IS THE LARGEST NUMBER OF COMPONENTS
C         REPRESENTABLE IN THE WORKING SPACE AVAILABLE.
C
      AVAIL = WRKLEN - OFFDIA
      MXCOMP = AVAIL/2
      START = OFFDIA + 1
      SIZE = START + MXCOMP
C
      IF  (MXCOMP .LE. 0)  GO TO 5100
C
      CALL GPSKCH (N, DEGREE, RSTART, CONNEC, TREE1, OFFDIA, WORK,
     1             MXCOMP, WORK(START), WORK(SIZE), COMPNS, ERROR,
     2             SPACE)
      IF ( ERROR .NE. 0 )  GO TO 5000
C
C     ... RECORD SPACE ACTUALLY USED  (NOT INCLUDING  NUMBERED )
C
      SPACE = 2*N + 3*(DEPTH) + 2*COMPNS + OFFDIA
C
C     ... SORT THE COMPONENT START POINTERS INTO INCREASING ORDER
C         OF SIZE OF COMPONENT
C
      IF (COMPNS .GT. 1)
     1    CALL GPSKCN (COMPNS, WORK(SIZE), WORK(START), ERROR)
          IF  (ERROR .NE. 0)  GO TO 6200
C
C     ... FOR EACH COMPONENT IN TURN, CHOOSE TO USE THE ORDERING OF THE
C         'FORWARD' TREE1 OR OF THE 'BACKWARD' TREE2 TO NUMBER THE NODES
C         IN THIS COMPONENT.  THE NUMBERING IS CHOSEN TO MINIMIZE THE
C         MAXIMUM INCREMENT TO ANY LEVEL.
C
      DO 1000 COMPON = 1, COMPNS
          PCSTRT = START + COMPON - 1
          CSTART = WORK (PCSTRT)
          PCSIZE = SIZE + COMPON - 1
          CSIZE = WORK (PCSIZE)
          CSTOP  = CSTART + CSIZE - 1
          IF ( ( CSIZE .LT. 0 ) .OR. ( CSIZE .GT. OFFDIA ) )  GO TO 6100
C
          DO 300 I = 1, DEPTH
              INC1(I) = 0
              INC2(I) = 0
  300     CONTINUE
C
          MXINC1 = 0
          MXINC2 = 0
C
          DO 400 I = CSTART, CSTOP
              WORKI = WORK(I)
              TWORKI = -TREE1 (WORKI)
              INC1 (TWORKI) = INC1 (TWORKI) + 1
              TWORKI =  TREE2 (WORKI)
              INC2 (TWORKI) = INC2 (TWORKI) + 1
  400     CONTINUE
C
C         ... BAROQUE TESTS BELOW DUPLICATE THE GIBBS-POOLE-STOCKMEYER-
C             CRANE PROGRAM, *** NOT *** THE PUBLISHED ALGORITHM.
C
          DO 500 I = 1, DEPTH
              IF ((INC1(I) .EQ. 0) .AND. (INC2(I) .EQ. 0))  GO TO 500
                  IF  (MXINC1  .LT.  TOTAL(I) + INC1(I))
     1                 MXINC1 = TOTAL(I) + INC1(I)
                  IF  (MXINC2  .LT.  TOTAL(I) + INC2(I))
     1                 MXINC2 = TOTAL(I) + INC2(I)
  500     CONTINUE
C
C         ... USE ORDERING OF NARROWER TREE UNLESS IT INCREASES
C             WIDTH MORE THAN WIDER TREE.  IN CASE OF TIE, USE TREE 2!
C
          IF ( (MXINC1 .GT. MXINC2)  .OR.
     1         ( (MXINC1 .EQ. MXINC2) .AND. ( (WIDTH1 .GT. WIDTH2) .OR.
     2                                        ( (WIDTH1 .EQ. WIDTH2)
     3                                         .AND. ONEIS1) ) ) )
     4      GO TO 700
C
              IF ( COMPON .EQ. 1 )  REVRS1 = .NOT. ONEIS1
C
              DO 600 I = 1, DEPTH
                  TOTAL(I) = TOTAL(I) + INC1(I)
  600         CONTINUE
              GO TO 1000
C
  700         IF ( COMPON .EQ. 1 )  REVRS1 = ONEIS1
              DO 800 I = CSTART, CSTOP
                  WORKI = WORK(I)
                  TREE1 (WORKI) = - TREE2 (WORKI)
  800         CONTINUE
C
              DO 900 I = 1, DEPTH
                  TOTAL(I) = TOTAL(I) + INC2(I)
  900         CONTINUE
C
 1000 CONTINUE
      GO TO 2000
C
C     ... DEFAULT WHEN THE REDUCED GRAPH IS EMPTY
C
 1100 REVRS1 = .TRUE.
      SPACE = 2*N
C
 2000 RETURN
C
C     ------------------------------------------------------------------
C
C     ERROR FOUND ...
C
 5000 SPACE = -1
      GO TO 2000
C
 5100 SPACE = 2 - AVAIL
      ERROR = 131
      GO TO 2000
C
 6000 ERROR = 30
      GO TO 5000
C
 6100 ERROR = 31
      GO TO 5000
C
 6200 ERROR = 32
      GO TO 5000
C
      END
      SUBROUTINE   GPSKCH   (N, DEGREE, RSTART, CONNEC, STATUS, NREDUC,
     1                       WORK, MXCOMP, START, SIZE, COMPNS, ERROR,
     2                       SPACE)
C
C     ==================================================================
C
C     FIND THE CONNECTED COMPONENTS OF THE GRAPH INDUCED BY THE SET
C     OF NODES WITH POSITIVE 'STATUS'.  WE SHALL BUILD THE LIST OF
C     CONNECTED COMPONENTS IN 'WORK', WITH A LIST OF POINTERS
C     TO THE BEGINNING NODES OF COMPONENTS LOCATED IN 'START'
C
C
      INTEGER     N, RSTART(N), NREDUC, MXCOMP, COMPNS, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(*), STATUS(N), WORK(NREDUC),
      INTEGER     DEGREE(N), CONNEC(*), STATUS(N), WORK(NREDUC),
     1            START(MXCOMP), SIZE(MXCOMP)
C
C
C     PARAMETERS ...
C
C         N      -- DIMENSION OF THE ORIGINAL MATRIX
C         DEGREE, RSTART, CONNEC -- THE STRUCTURE OF THE ORIGINAL MATRIX
C
C         STATUS -- DERIVED FROM A LEVEL TREE. POSITIVE ENTRIES INDICATE
C                   ACTIVE NODES.  NODES WITH STATUS <= 0 ARE IGNORED.
C
C         NREDUC -- THE NUMBER OF ACTIVE NODES
C
C         WORK   -- WORK SPACE, USED AS A QUEUE TO BUILD CONNECTED
C                   COMPONENTS IN PLACE.
C
C         MXCOMP -- MAXIMUM NUMBER OF COMPONENTS ALLOWED BY CURRENT
C                   SPACE ALLOCATION.  MUST NOT BE VIOLATED.
C
C         START  -- POINTER TO BEGINNING OF  I-TH  CONNECTED COMPONENT
C
C         SIZE   -- SIZE OF EACH COMPONENT
C
C         COMPNS -- NUMBER OF COMPONENTS ACTUALLY FOUND
C
C         ERROR  -- SHOULD BE ZERO ON RETURN UNLESS WE HAVE TOO LITTLE
C                   SPACE OR WE ENCOUNTER AN ERROR IN THE DATA STRUCTURE
C
C         SPACE  -- MAXIMUM AMOUNT OF WORKSPACE USED / NEEDED
C
C     ==================================================================
C
      INTEGER     I, J, FREE, JPTR, NODE, JNODE, FRONT, CDGREE, ROOT
C
C     ------------------------------------------------------------------
C
C
C     REPEAT
C         << FIND AN UNASSIGNED NODE AND START A NEW COMPONENT >>
C         REPEAT
C             << ADD ALL NEW NEIGHBORS OF FRONT NODE TO QUEUE, >>
C             << REMOVE FRONT NODE.                            >>
C         UNTIL <<QUEUE EMPTY>>
C     UNTIL << ALL NODES ASSIGNED >>
C
      FREE   = 1
      COMPNS = 0
      ROOT   = 1
C
C     ... START OF OUTER REPEAT LOOP
C
C         ... FIND AN UNASSIGNED NODE
C
  100     DO 200 I = ROOT, N
              IF (STATUS(I) .LE. 0) GO TO 200
                  NODE = I
                  GO TO 300
  200     CONTINUE
          GO TO 6100
C
C         ... START NEW COMPONENT
C
  300     COMPNS = COMPNS + 1
          ROOT   = NODE + 1
          IF (COMPNS .GT. MXCOMP)  GO TO 5000
          START (COMPNS) = FREE
          WORK (FREE) = NODE
          STATUS (NODE) = -STATUS (NODE)
          FRONT = FREE
          FREE = FREE + 1
C
C             ... INNER REPEAT UNTIL QUEUE BECOMES EMPTY
C
  400         NODE = WORK (FRONT)
              FRONT = FRONT + 1
C
              JPTR = RSTART (NODE)
              CDGREE = DEGREE (NODE)
              DO 500 J = 1, CDGREE
                  JNODE = CONNEC (JPTR)
                  JPTR = JPTR + 1
                  IF (STATUS(JNODE) .LT. 0) GO TO 500
                  IF (STATUS(JNODE) .EQ. 0) GO TO 6000
                      STATUS (JNODE) = -STATUS (JNODE)
                      WORK (FREE) = JNODE
                      FREE = FREE + 1
  500         CONTINUE
C
              IF (FRONT .LT. FREE) GO TO 400
C
C         ... END OF INNER REPEAT.  COMPUTE SIZE OF COMPONENT AND
C             SEE IF THERE ARE MORE NODES TO BE ASSIGNED
C
          SIZE (COMPNS) = FREE - START (COMPNS)
          IF (FREE .LE. NREDUC)  GO TO 100
C
      IF (FREE .NE. NREDUC+1)  GO TO 6200
      RETURN
C
C     ------------------------------------------------------------------
C
 5000 SPACE = NREDUC - FREE + 1
      ERROR = 130
      RETURN
C
 6000 ERROR = 33
      SPACE = -1
      RETURN
C
 6100 ERROR = 34
      SPACE = -1
      RETURN
C
 6200 ERROR = 35
      SPACE = -1
      RETURN
      END
      SUBROUTINE   GPSKCI   (N, ACTIVE, DEPTH, LSTRUC, LVLLST, LVLPTR,
     1                       LTOTAL, ERROR, SPACE)
C
C     ==================================================================
C
C     TRANSITIONAL SUBROUTINE, ALGORITHM II TO IIIA OR IIIB.
C
C     CONVERT LEVEL STRUCTURE GIVEN AS VECTOR OF LEVEL NUMBERS FOR NODES
C     TO STRUCTURE AS LIST OF NODES BY LEVEL
C
C     N, ACTIVE, DEPTH -- PROBLEM SIZES
C     LSTRUC -- INPUT LEVEL STRUCTURE
C     LVLLST, LVLPTR -- OUTPUT LEVEL STRUCTURE
C     LTOTAL -- NUMBER OF NODES AT EACH LEVEL (PRECOMPUTED)
C
      INTEGER     N, ACTIVE, DEPTH, ERROR, SPACE
C
CIBM  INTEGER *2  LSTRUC(N), LVLLST(ACTIVE), LVLPTR(*), LTOTAL(DEPTH)
      INTEGER     LSTRUC(N), LVLLST(ACTIVE), LVLPTR(*), LTOTAL(DEPTH)
C
C     ===============================================================
C
C     STRUCTURE OF WORKSPACE ..
C
C         INPUT (FROM COMBIN) ..
C
C     ------------------------------------------------------------------
C     :  NUMBERED  :  ..(N)..  :  TOTAL  :         ...        :  TREE  :
C     ------------------------------------------------------------------
C
C         OUTPUT (TO GPSKCJ OR GPSKCK) ..
C
C     ------------------------------------------------------------------
C     :  NUMBERED  :       ...             :  TLIST  :  TPTR  :  TREE  :
C     ------------------------------------------------------------------
C
C     HERE, NUMBERED IS THE SET OF NODES IN NUMBERED COMPONENTS
C         TOTAL IS A VECTOR OF LENGTH 'DEPTH' GIVING THE NUMBER
C         OF NODES IN EACH LEVEL OF THE 'TREE'.
C         TLIST, TPTR ARE LISTS OF NODES OF THE TREE, ARRANGED
C         BY LEVEL.  TLIST IS OF LENGTH 'ACTIVE', TPTR 'DEPTH+1'.
C
C     =================================================================
C
      INTEGER     I, ACOUNT, START, LEVEL, PLEVEL
C
C     ... ESTABLISH STARTING AND ENDING POINTERS FOR EACH LEVEL
C
      START = 1
      DO 100 I = 1, DEPTH
          LVLPTR(I) = START
          START = START + LTOTAL(I)
          LTOTAL(I) = START
  100 CONTINUE
      LVLPTR(DEPTH+1) = START
C
      ACOUNT = 0
      DO 300 I = 1, N
          IF (LSTRUC(I) .LT. 0) GO TO 200
          IF (LSTRUC(I) .EQ. 0) GO TO 300
          IF (LSTRUC(I) .GT. 0) GO TO 6000
  200         LEVEL = -LSTRUC(I)
              LSTRUC(I) = LEVEL
              PLEVEL = LVLPTR (LEVEL)
              LVLLST (PLEVEL) = I
              LVLPTR (LEVEL) = LVLPTR (LEVEL) + 1
              ACOUNT = ACOUNT + 1
              IF (LVLPTR (LEVEL) .GT. LTOTAL (LEVEL))  GO TO 6100
  300 CONTINUE
C
C     ... RESET STARTING POINTERS
C
      LVLPTR(1) = 1
      DO 400 I = 1, DEPTH
          LVLPTR(I+1) = LTOTAL(I)
  400 CONTINUE
C
      RETURN
C
C     ------------------------------------------------------------------
C
 6000 ERROR = 40
      GO TO 6200
C
 6100 ERROR = 41
C
 6200 SPACE = -1
      RETURN
C
      END
      SUBROUTINE   GPSKCJ   (N, DEGREE, RSTART, CONNEC,
     1                       NCOMPN, INVNUM, SNODE1, SNODE2, REVRS1,
     2                       DEPTH, LVLLST, LVLPTR, LVLNUM, ERROR,
     3                       SPACE)
C
C     ==================================================================
C
C     NUMBER THE NODES IN A GENERALIZED LEVEL STRUCTURE ACCORDING
C     TO A GENERALIZATION OF THE CUTHILL MCKEE STRATEGY.
C
C     N      -- DIMENSION OF ORIGINAL PROBLEM
C     DEGREE, RSTART, CONNEC -- GIVE STRUCTURE OF SPARSE AND
C                               SYMMETRIC MATRIX
C
C     NCOMPN -- NUMBER OF NODES IN THIS COMPONENT OF MATRIX GRAPH
C
C     INVNUM -- WILL BECOME A LIST OF THE ORIGINAL NODES IN THE ORDER
C               WHICH REDUCES THE BANDWIDTH OF THE MATRIX.
C
C     NXTNUM -- THE NEXT INDEX TO BE ASSIGNED (1 FOR FIRST COMPONENT)
C
C     REVRS1 -- IF .TRUE., FIRST COMPONENT OF REDUCED GRAPH WAS NUMBERED
C               BACKWARDS.
C
C     LVLLST -- LIST OF NODES IN LEVEL TREE ORDERED BY LEVEL.
C
C     LVLPTR -- POSITION OF INITIAL NODE IN EACH LEVEL OF LVLLST.
C
C     LVLNUM -- LEVEL NUMBER OF EACH NODE IN COMPONENT
C
C
      INTEGER     N, RSTART(N), NCOMPN, SNODE1, SNODE2, DEPTH,
     1            ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(*), INVNUM(NCOMPN),
      INTEGER     DEGREE(N), CONNEC(*), INVNUM(NCOMPN),
     1            LVLLST(NCOMPN), LVLPTR(DEPTH), LVLNUM(N)
C
      LOGICAL     REVRS1
C
C
C     ==================================================================
C
C     NUMBERING REQUIRES TWO QUEUES, WHICH CAN BE BUILD IN PLACE
C     IN INVNUM.
C
C
C     ==================================================================
C     A L G O R I T H M    S T R U C T U R E
C     ==================================================================
C
C     << SET QUEUE1 TO BE THE SET CONTAINING ONLY THE START NODE. >>
C
C     FOR LEVEL = 1 TO DEPTH DO
C
C         BEGIN
C         LOOP
C
C             REPEAT
C                 BEGIN
C                 << CNODE <- FRONT OF QUEUE1                        >>
C                 << ADD UNNUMBERED NEIGHBORS OF CNODE TO THE BACK   >>
C                 << OF QUEUE1 OR QUEUE2 (USE QUEUE1 IF NEIGHBOR     >>
C                 << AT SAME LEVEL, QUEUE2 IF AT NEXT LEVEL).  SORT  >>
C                 << THE NEWLY QUEUED NODES INTO INCREASING ORDER OF >>
C                 << DEGREE.  NUMBER CNODE, DELETE IT FROM QUEUE1.   >>
C                 END
C             UNTIL
C                 << QUEUE1 IS EMPTY >>
C
C         EXIT IF << ALL NODES AT THIS LEVEL NUMBERED >>
C
C             BEGIN
C             << FIND THE UNNUMBERED NODE OF MINIMAL DEGREE AT THIS >>
C             << LEVEL, RESTART QUEUE1 WITH THIS NODE.              >>
C             END
C
C         END << LOOP LOOP >>
C
C         << PROMOTE QUEUE2 TO BE INITIAL QUEUE1 FOR NEXT ITERATION >>
C         << OF  FOR  LOOP.                                         >>
C
C         END <<FOR LOOP>>
C
C     ==================================================================
C
C     STRUCTURE OF WORKSPACE ..
C
C     --------------------------------------------------------------
C     : NUMBERED :  QUEUE1  :  QUEUE2  : ... : TLIST : TPTR : TREE :
C     --------------------------------------------------------------
C
C     ON COMPLETION, WE HAVE ONLY A NEW, LONGER NUMBERED SET.
C
C     ==================================================================
      INTEGER     I, BQ1, BQ2, FQ1, INC, CPTR, CNODE,
     1            INODE, LEVEL, NLEFT, LSTART, LWIDTH, QUEUE1,
     2            QUEUE2, CDGREE, XLEVEL, STNODE, ILEVEL, SQ1, SQ2,
     3            NSORT, LOWDG, BPTR, LVLLSC, LVLLSB, INVNMI
C
      LOGICAL     FORWRD, RLEVEL
C
C     ------------------------------------------------------------------
C
C     ... GIBBS-POOLE-STOCKMEYER HEURISTIC CHOICE OF ORDER
C
      IF  (DEGREE(SNODE1) .GT. DEGREE(SNODE2))  GO TO 10
          FORWRD = REVRS1
          STNODE = SNODE1
          GO TO 20
C
   10     FORWRD = .NOT. REVRS1
          STNODE = SNODE2
C
C     ... SET UP INITIAL QUEUES AT FRONT OF 'INVNUM' FOR FORWRD ORDER,
C         AT BACK FOR REVERSED ORDER.
C
   20 IF (FORWRD) GO TO 100
          INC = -1
          QUEUE1 = NCOMPN
          GO TO 200
C
  100     INC = +1
          QUEUE1 = 1
C
  200 INVNUM (QUEUE1) = STNODE
      RLEVEL = (LVLNUM(STNODE) .EQ. DEPTH)
      LVLNUM (STNODE) = 0
      FQ1 = QUEUE1
      BQ1 = QUEUE1 + INC
C
C     -------------------------------
C     NUMBER NODES LEVEL BY LEVEL ...
C     -------------------------------
C
      DO 3000 XLEVEL = 1, DEPTH
          LEVEL = XLEVEL
          IF  (RLEVEL)  LEVEL = DEPTH - XLEVEL + 1
C
          LSTART = LVLPTR (LEVEL)
          LWIDTH = LVLPTR (LEVEL+1) - LSTART
          NLEFT = LWIDTH
          QUEUE2 = QUEUE1 + INC*LWIDTH
          BQ2 = QUEUE2
C
C         ==============================================================
C         ... 'LOOP' CONSTRUCT BEGINS AT STATEMENT 1000
C                 THE INNER 'REPEAT' WILL BE DONE AS MANY TIMES AS
C                 IS NECESSARY TO NUMBER ALL THE NODES AT THIS LEVEL.
C         ==============================================================
C
 1000     CONTINUE
C
C             ==========================================================
C             ... REPEAT ... UNTIL QUEUE1 BECOMES EMPTY
C                 TAKE NODE FROM FRONT OF QUEUE1, FIND EACH OF ITS
C                 NEIGHBORS WHICH HAVE NOT YET BEEN NUMBERED, AND
C                 ADD THE NEIGHBORS TO QUEUE1 OR QUEUE2 ACCORDING TO
C                 THEIR LEVELS.
C             ==========================================================
C
 1100             CNODE = INVNUM (FQ1)
                  FQ1 = FQ1 + INC
                  SQ1 = BQ1
                  SQ2 = BQ2
                  NLEFT = NLEFT - 1
C
                  CPTR = RSTART (CNODE)
                  CDGREE = DEGREE (CNODE)
                  DO 1300 I = 1, CDGREE
                      INODE = CONNEC (CPTR)
                      CPTR = CPTR + 1
                      ILEVEL = LVLNUM (INODE)
                      IF (ILEVEL .EQ. 0)  GO TO 1300
                          LVLNUM (INODE) = 0
                          IF ( ILEVEL .EQ. LEVEL ) GO TO 1200
C
                              IF  (IABS(LEVEL-ILEVEL) .NE. 1) GO TO 6400
                                  INVNUM (BQ2) = INODE
                                  BQ2 = BQ2 + INC
                                  GO TO 1300
C
 1200                             INVNUM (BQ1) = INODE
                                  BQ1 = BQ1 + INC
 1300             CONTINUE
C
C                 ==================================================
C                 ... SORT THE NODES JUST ADDED TO QUEUE1 AND QUEUE2
C                     SEPARATELY INTO INCREASING ORDER OF DEGREE.
C                 ==================================================
C
                  IF  (IABS (BQ1 - SQ1) .LE. 1)  GO TO 1500
                      NSORT = IABS (BQ1 - SQ1)
                      IF  (FORWRD)  GO TO 1400
                          CALL GPSKCP (NSORT, INVNUM(BQ1+1), N, DEGREE,
     1                                 ERROR)
                          IF  (ERROR .NE. 0)  GO TO 6600
                          GO TO 1500
C
 1400                     CALL GPSKCQ (NSORT, INVNUM(SQ1), N, DEGREE,
     1                                 ERROR)
                          IF  (ERROR .NE. 0)  GO TO 6600
C
 1500             IF  (IABS (BQ2 - SQ2) .LE. 1)  GO TO 1700
                      NSORT = IABS (BQ2 - SQ2)
                      IF  (FORWRD)  GO TO 1600
                          CALL GPSKCP (NSORT, INVNUM(BQ2+1), N, DEGREE,
     1                                 ERROR)
                          IF  (ERROR .NE. 0)  GO TO 6600
                          GO TO 1700
C
 1600                     CALL GPSKCQ (NSORT, INVNUM(SQ2), N, DEGREE,
     1                                 ERROR)
                          IF  (ERROR .NE. 0)  GO TO 6600
C
C                     ... END OF REPEAT LOOP
C
 1700             IF  (FQ1 .NE. BQ1)  GO TO 1100
C
C         ==============================================================
C         ... QUEUE1 IS NOW EMPTY ...
C             IF THERE ARE ANY UNNUMBERED NODES LEFT AT THIS LEVEL,
C             FIND THE ONE OF MINIMAL DEGREE AND RETURN TO THE
C             REPEAT LOOP ABOVE.
C         ==============================================================
C
 2000     IF  ((BQ1 .EQ. QUEUE2) .AND. (NLEFT .EQ. 0))  GO TO 2900
C
              IF ((NLEFT .LE. 0) .OR. (NLEFT .NE. INC * (QUEUE2 - BQ1)))
     1             GO TO 6200
C
              LOWDG = N + 1
              BPTR  = N + 1
              CPTR  = LSTART - 1
              DO 2800 I = 1, NLEFT
 2600             CPTR   = CPTR + 1
                  LVLLSC = LVLLST (CPTR)
                  IF (LVLNUM (LVLLSC) .EQ. LEVEL)  GO TO 2700
                      IF (LVLNUM (LVLLSC) .NE. 0)  GO TO 6300
                      GO TO 2600
C
 2700             IF  (DEGREE(LVLLSC) .GE. LOWDG)  GO TO 2800
                      LOWDG = DEGREE (LVLLSC)
                      BPTR  = CPTR
C
 2800         CONTINUE
C
C             ... MINIMAL DEGREE UNNUMBERED NODE FOUND ...
C
              IF  (BPTR .GT. N)  GO TO 6500
              LVLLSB = LVLLST (BPTR)
              INVNUM (BQ1) = LVLLSB
              LVLNUM (LVLLSB) = 0
              BQ1 = BQ1 + INC
              GO TO 1000
C
C             =============================================
C             ... ADVANCE QUEUE POINTERS TO MAKE QUEUE2 THE
C                 NEW QUEUE1 FOR THE NEXT ITERATION.
C             =============================================
C
 2900     QUEUE1 = QUEUE2
          FQ1 = QUEUE1
          BQ1 = BQ2
          IF  ((BQ1 .EQ. FQ1) .AND. (XLEVEL .LT. DEPTH))  GO TO 6100
C
 3000 CONTINUE
C
C     ... CHANGE SIGN OF DEGREE TO MARK THESE NODES AS 'NUMBERED'
C
      DO 3100 I = 1, NCOMPN
          INVNMI = INVNUM(I)
          DEGREE (INVNMI) = -DEGREE (INVNMI)
 3100 CONTINUE
C
      RETURN
C
C     ------------------------------------------------------------------
C
 6000 SPACE = -1
      RETURN
C
 6100 ERROR = 51
      GO TO 6000
C
 6200 ERROR = 52
      GO TO 6000
C
 6300 ERROR = 53
      GO TO 6000
C
 6400 ERROR = 54
      GO TO 6000
C
 6500 ERROR = 55
      GO TO 6000
C
 6600 ERROR = 56
      GO TO 6000
C
      END
      SUBROUTINE  GPSKCK  (N, DEGREE, RSTART, CONNEC, WRKLEN, NXTNUM,
     1                     WORK, NCOMPN, DEPTH, LVLLST, LVLPTR, LVLNUM,
     2                     ERROR, SPACE)
C
      INTEGER     N, RSTART(N), WRKLEN, NXTNUM, NCOMPN, DEPTH, ERROR,
     1            SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(*), WORK(WRKLEN), LVLLST(N),
      INTEGER     DEGREE(N), CONNEC(*), WORK(WRKLEN), LVLLST(N),
     1            LVLPTR(DEPTH), LVLNUM(N)
C
C     ==================================================================
C
C     NUMBER NODES IN A GENERALIZED LEVEL STRUCTURE ACCORDING TO
C     A GENERALIZATION OF THE KING ALGORITHM, WHICH REDUCES
C     THE PROFILE OF THE SPARSE SYMMETRIC MATRIX.
C
C     ---------------------
C
C     CODE USES A PRIORITY QUEUE TO CHOOSE THE NEXT NODE TO BE NUMBERED
C     THE PRIORITY QUEUE IS REPRESENTED BY A SIMPLE LINEAR-LINKED LIST
C     TO SAVE SPACE.  THIS WILL REQUIRE MORE SEARCHING THAN A FULLY
C     LINKED REPRESENTATION, BUT THE DATA MANIPULATION IS SIMPLER.
C
C     -------------------
C
C     << ESTABLISH PRIORITY QUEUE 'ACTIVE' FOR LEVEL 1 NODES >>
C
C     FOR I = 1 TO DEPTH DO
C         << SET QUEUE 'QUEUED' TO BE EMPTY, LIST 'NEXT' TO BE >>
C         << SET OF NODES AT NEXT LEVEL.                       >>
C
C         FOR J = 1 TO 'NODES AT THIS LEVEL' DO
C             << FIND FIRST NODE IN ACTIVE WITH MINIMAL CONNECTIONS >>
C             << TO 'NEXT'.  NUMBER THIS NODE AND REMOVE HIM FROM   >>
C             << 'ACTIVE'.  FOR EACH NODE IN 'NEXT' WHICH CONNECTED >>
C             << TO THIS NODE, MOVE IT TO 'QUEUED' AND REMOVE IT    >>
C             << FROM 'NEXT'.                                       >>
C
C         << SET NEW QUEUE 'ACTIVE' TO BE 'QUEUED' FOLLOWED BY ANY >>
C         << NODES STILL IN 'NEXT'.                                >>
C
C     ==================================================================
C
C     DATA STRUCTURE ASSUMPTIONS:
C     THE FIRST 'NXTNUM-1' ELEMENTS OF  WORK  ARE ALREADY IN USE.
C     THE LEVEL STRUCTURE 'LVLLST' IS CONTIGUOUS WITH  WORK, THAT IS,
C     IT RESIDES IN ELEMENTS  WRKLEN+1, ...  OF  WORK.  'LVLPTR' AND
C     'LVLNUM' ARE ALSO EMBEDDED IN WORK, BEHIND 'LVLLST'.  THE
C     THREE VECTORS ARE PASSED SEPARATELY TO CLARIFY THE INDEXING,
C     BUT THE QUEUES DEVELOPED WILL BE ALLOWED TO OVERRUN 'LVLLST'
C     AS NEEDED.
C
C     ... BUILD THE FIRST 'ACTIVE' QUEUE STARTING W1 LOCATIONS FROM
C         THE FRONT OF THE CURRENT WORKING AREA  (W1 IS THE WIDTH OF THE
C         FIRST LEVEL).  BUILD THE FIRST 'QUEUED' QUEUE STARTING FROM
C         THE BACK OF WORK SPACE.  THE LIST 'NEXT' WILL BE REALIZED
C         IMPLICITLY IN 'LVLNUM' AS:
C                  LVLNUM(I) > 0   <== LEVEL NUMBER OF NODE.  'NEXT' IS
C                                      SET WITH LVLNUM(I) = LEVEL+1
C                  LVLNUM(I) = 0   <== I-TH NODE IS IN 'QUEUED' OR IS
C                                      NOT IN THIS COMPONENT OF GRAPH,
C                                      OR HAS JUST BEEN NUMBERED.
C                  LVLNUM(I) < 0   <== I-TH NODE IS IN 'ACTIVE' AND IS
C                                      CONNECTED TO -LVLNUM(I) NODES IN
C                                      'NEXT'.
C
C     ==================================================================
C
C     STRUCTURE OF WORKSPACE ..
C
C     --------------------------------------------------------------
C     : NUMBERED : DONE : ACTIVE : ALEVEL : ... : QUEUED : LVLLST :
C     --------------------------------------------------------------
C
C     -------------------
C       LVLPTR : LVLNUM :
C     -------------------
C
C     IN THE ABOVE,
C         NUMBERED IS THE SET OF NODES ALREADY NUMBERED FROM
C         PREVIOUS COMPONENTS AND EARLIER LEVELS OF THIS COMPONENT.
C         DONE, ACTIVE, ALEVEL  ARE VECTORS OF LENGTH THE WIDTH OF
C         THE CURRENT LEVEL.  ACTIVE IS A SET OF INDICES INTO
C         ALEVEL.  AS THE NODES IN ALEVEL ARE NUMBERED, THEY
C         ARE PLACED INTO 'DONE'.
C         QUEUED IS A QUEUE OF NODES IN THE 'NEXT' LEVEL, WHICH
C         GROWS FROM THE START OF THE 'NEXT' LEVEL IN LVLLST
C         FORWARDS TOWARD 'ALEVEL'.  QUEUED IS OF LENGTH NO MORE
C         THAN THE WIDTH OF THE NEXT LEVEL.
C         LVLLST IS THE LIST OF UNNUMBERED NODES IN THE TREE,
C         ARRANGED BY LEVEL.
C
C     ==================================================================
      INTEGER     I, J, K, PTR, JPTR, KPTR, LPTR, MPTR, PPTR, RPTR,
     1            MPPTR, JNODE, KNODE, CNODE, LEVEL, LOWDG, UNUSED,
     2            MXQUE, NNEXT, ASTART, MINDG, LSTART, LWIDTH, ACTIVE,
     2            QUEUEB, QUEUED, QCOUNT, NCONNC, NACTIV, CDGREE,
     3            LDGREE, NFINAL, JDGREE, STRTIC, ADDED, TWRKLN,
     4            LVLLSL, CONNEJ, CONNER, ASTPTR, ACTPTR, ACTIVI,
     5            ASTRTI, QUEUEI, ACPPTR
C
C     ------------------------------------------------------------------
C
      TWRKLN = WRKLEN + NCOMPN + N + DEPTH + 1
      UNUSED = TWRKLN
C
      ASTART = LVLPTR(1)
      LWIDTH = LVLPTR(2) - ASTART
      ASTART = WRKLEN  + 1
      ACTIVE = NXTNUM + LWIDTH + 1
      NACTIV = LWIDTH
      NFINAL = NXTNUM + NCOMPN
C
      NNEXT = LVLPTR(3) - LVLPTR(2)
      QUEUED = WRKLEN
      QUEUEB = QUEUED
      MXQUE = ACTIVE + LWIDTH
C
C     ... BUILD FIRST PRIORITY QUEUE 'ACTIVE'
C
      LOWDG = - (N + 1)
      LPTR = LVLPTR(1)
      DO 200 I = 1, LWIDTH
          NCONNC = 0
          LVLLSL= LVLLST (LPTR)
          JPTR = RSTART (LVLLSL)
          LDGREE = DEGREE(LVLLSL)
          DO 100 J = 1, LDGREE
              CONNEJ = CONNEC (JPTR)
              IF ( LVLNUM (CONNEJ) .EQ. 2 )  NCONNC = NCONNC - 1
              JPTR = JPTR + 1
  100     CONTINUE
C
          ACTIVI = ACTIVE + I - 1
          WORK (ACTIVI) = I
          LVLNUM (LVLLSL) = NCONNC
          LOWDG = MAX0 (LOWDG, NCONNC)
          LPTR = LPTR + 1
  200 CONTINUE
      WORK (ACTIVE-1) = 0
C
C     -----------------------------------
C     NOW NUMBER NODES LEVEL BY LEVEL ...
C     -----------------------------------
C
      DO 2000 LEVEL = 1, DEPTH
C
C         ... NUMBER ALL NODES IN THIS LEVEL
C
          DO 1100 I = 1, LWIDTH
              PPTR = -1
              PTR = WORK (ACTIVE-1)
              IF (NNEXT .EQ. 0)  GO TO 1000
C
C                 ... IF NODES REMAIN IN NEXT, FIND THE EARLIEST NODE
C                     IN ACTIVE OF MINIMAL DEGREE.
C
                  MINDG = -(N+1)
                  DO 400 J = 1, NACTIV
                      ASTPTR = ASTART + PTR
                      CNODE = WORK (ASTPTR)
                      IF ( LVLNUM (CNODE) .EQ. LOWDG )  GO TO 500
                      IF ( LVLNUM (CNODE) .LE. MINDG )  GO TO 300
                          MPPTR = PPTR
                          MPTR = PTR
                          MINDG = LVLNUM (CNODE)
  300                 PPTR = PTR
                      ACTPTR = ACTIVE + PTR
                      PTR = WORK (ACTPTR)
  400             CONTINUE
C
C                     ... ESTABLISH  PTR  AS FIRST MIN DEGREE NODE
C                         PPTR AS PREDECESSOR IN LIST.
C
                  PTR = MPTR
                  PPTR = MPPTR
C
  500             ASTPTR = ASTART + PTR
                  CNODE = WORK (ASTPTR)
                  LOWDG = LVLNUM (CNODE)
                  LVLNUM (CNODE) = 0
                  JPTR = RSTART (CNODE)
C
C                 ... UPDATE CONNECTION COUNTS FOR ALL NODES WHICH
C                     CONNECT TO  CNODE'S  NEIGHBORS IN  NEXT.
C
                  CDGREE = DEGREE(CNODE)
                  STRTIC = QUEUEB
C
                  DO 700 J = 1, CDGREE
                      JNODE = CONNEC (JPTR)
                      JPTR = JPTR + 1
                      IF (LVLNUM (JNODE) .NE. LEVEL+1 )  GO TO 700
                          IF (QUEUEB .LT. MXQUE)  GO TO 5000
                          WORK (QUEUEB) = JNODE
                          QUEUEB = QUEUEB - 1
                          NNEXT = NNEXT - 1
                          LVLNUM (JNODE) = 0
                          IF  (NACTIV .EQ. 1)  GO TO 700
                            KPTR = RSTART (JNODE)
                            JDGREE = DEGREE (JNODE)
                            DO 600 K = 1, JDGREE
                                KNODE = CONNEC (KPTR)
                                KPTR = KPTR + 1
                                IF (LVLNUM (KNODE) .GE. 0)  GO TO 600
                                    LVLNUM (KNODE) = LVLNUM (KNODE) + 1
                                    IF  (LOWDG .LT. LVLNUM(KNODE))
     1                                   LOWDG = LVLNUM(KNODE)
  600                       CONTINUE
  700             CONTINUE
C
C                 ... TO MIMIC THE ALGORITHM AS IMPLEMENTED BY GIBBS,
C                     SORT THE NODES JUST ADDED TO THE QUEUE INTO
C                     INCREASING ORDER OF ORIGINAL INDEX. (BUT, BECAUSE
C                     THE QUEUE IS STORED BACKWARDS IN MEMORY, THE SORT
C                     ROUTINE IS CALLED FOR DECREASING INDEX.)
C
C                     TREAT  0, 1 OR 2  NODES ADDED AS SPECIAL CASES
C
                  ADDED = STRTIC - QUEUEB
                  IF  (ADDED - 2 .LT. 0)  GO TO 1000
                  IF  (ADDED - 2 .EQ. 0)  GO TO 800
                  IF  (ADDED - 2 .GT. 0)  GO TO 900
C
  800                 IF (WORK(STRTIC-1) .GT. WORK(STRTIC))  GO TO 1000
                          JNODE = WORK(STRTIC)
                          WORK(STRTIC) = WORK(STRTIC-1)
                          WORK(STRTIC-1) = JNODE
                          GO TO 1000
C
  900                 CALL GPSKCO (ADDED, WORK(QUEUEB+1), ERROR)
                      IF  (ERROR .NE. 0)  GO TO 5500
C
C
C                 ... NUMBER THIS NODE AND DELETE IT FROM 'ACTIVE'.
C                     MARK IT UNAVAILABLE BY CHANGING SIGN OF DEGREE
C
 1000         NACTIV = NACTIV - 1
              ASTPTR = ASTART + PTR
              CNODE = WORK (ASTPTR)
              WORK (NXTNUM) = CNODE
              DEGREE (CNODE) = -DEGREE (CNODE)
              NXTNUM = NXTNUM + 1
C
C             ... DELETE LINK TO THIS NODE FROM LIST
C
              ACPPTR = ACTIVE + PPTR
              ACTPTR = ACTIVE + PTR
              WORK (ACPPTR) = WORK (ACTPTR)
 1100     CONTINUE
C
C         ... NOW MOVE THE QUEUE 'QUEUED' FORWARD, AT THE SAME
C             TIME COMPUTING CONNECTION COUNTS FOR ITS ELEMENTS.
C             THEN DO THE SAME FOR THE REMAINING NODES IN 'NEXT'.
C
          UNUSED = MIN0 (UNUSED, QUEUEB - MXQUE)
          IF ( NXTNUM .NE. ACTIVE-1 )  GO TO 5100
          IF ( LEVEL .EQ. DEPTH ) GO TO 2000
              LSTART = LVLPTR (LEVEL+1)
              LWIDTH = LVLPTR (LEVEL+2) - LSTART
              ACTIVE = NXTNUM + LWIDTH + 1
              ASTART = ACTIVE + LWIDTH
              NACTIV = LWIDTH
              MXQUE = ASTART + LWIDTH
              IF ( MXQUE .GT. QUEUEB + 1 )  GO TO 5000
              UNUSED = MIN0 (UNUSED, QUEUEB - MXQUE + 1)
C
              QCOUNT = QUEUED - QUEUEB
              LOWDG = -N-1
              WORK (ACTIVE-1) = 0
C
              PTR = LSTART
              DO 1600 I = 1, LWIDTH
C
C                 ... CHOOSE NEXT NODE FROM EITHER 'QUEUED' OR 'NEXT'
C
                  IF (I .GT. QCOUNT )  GO TO 1200
                      QUEUEI = QUEUED + 1 - I
                      CNODE = WORK (QUEUEI)
                      GO TO 1300
C
 1200                 CNODE = LVLLST (PTR)
                      PTR = PTR + 1
                      IF ( PTR .GT. LVLPTR(LEVEL+2) )  GO TO 5200
                          IF (LVLNUM (CNODE) .GT. 0)  GO TO 1300
                              GO TO 1200
C
 1300             IF ( LEVEL+1 .EQ. DEPTH ) GO TO 1500
C
                      RPTR = RSTART (CNODE)
                      NCONNC = 0
                      JDGREE = DEGREE (CNODE)
                      DO 1400 J = 1, JDGREE
                          CONNER = CONNEC (RPTR)
                          IF ( LVLNUM (CONNER) .EQ. LEVEL+2 )
     1                        NCONNC = NCONNC - 1
                          RPTR = RPTR + 1
 1400                 CONTINUE
                      LVLNUM (CNODE) = NCONNC
                      LOWDG = MAX0 (LOWDG, NCONNC)
C
C             ... ADD CNODE TO NEW 'ACTIVE' QUEUE
C
 1500             ACTIVI = ACTIVE + (I - 1)
                  ASTRTI = ASTART + (I - 1)
                  WORK (ACTIVI) = I
                  WORK (ASTRTI) = CNODE
 1600         CONTINUE
C
              IF (DEPTH .EQ. LEVEL+1 ) GO TO 1700
                  NNEXT = LVLPTR (LEVEL+3) - LVLPTR (LEVEL+2)
                  QUEUED = LSTART - 1 + LWIDTH + WRKLEN
                  QUEUEB = QUEUED
                  GO TO 2000
C
 1700             NNEXT = 0
C
 2000 CONTINUE
C
      IF  (NXTNUM .NE. NFINAL)  GO TO 5300
      SPACE = MAX0 (SPACE, TWRKLN - UNUSED)
      RETURN
C
C
C     ------------------------------------------------------------------
C
 5000 SPACE = NACTIV + NNEXT
      ERROR = 160
      RETURN
C
 5100 ERROR = 61
      GO TO 5400
C
 5200 ERROR = 62
      GO TO 5400
C
 5300 ERROR = 63
C
 5400 RETURN
C
 5500 ERROR = 64
      GO TO 5400
C
      END
      SUBROUTINE   GPSKCL   (N, DEGREE, RSTART, CONNEC, INVNUM, NEWNUM,
     1                       OLDNUM, BANDWD, PROFIL, ERROR, SPACE)
C
C
      INTEGER     N, RSTART(N), BANDWD, PROFIL, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(*), INVNUM(N), NEWNUM(N), OLDNUM(N)
      INTEGER     DEGREE(N), CONNEC(*), INVNUM(N), NEWNUM(N), OLDNUM(N)
C
C     ==================================================================
C
C
C     COMPUTE THE BANDWIDTH AND PROFILE FOR THE RENUMBERING GIVEN
C     BY 'INVNUM' AND ALSO FOR THE RENUMBERING GIVEN BY 'OLDNUM'.
C     'NEWNUM' WILL BE A PERMUTATION VECTOR COPY OF THE NODE
C     LIST 'INVNUM'.
C
C     ==================================================================
C
      INTEGER     I, J, JPTR, IDGREE, OLDBND, OLDPRO, NEWBND, NEWPRO,
     1            OLDRWD, NEWRWD, OLDORG, NEWORG, JNODE, INVNMI
C
C     ------------------------------------------------------------------
C
C     ... CREATE NEWNUM AS A PERMUTATION VECTOR
C
      DO 100 I = 1, N
          INVNMI = INVNUM (I)
          NEWNUM (INVNMI) = I
  100 CONTINUE
C
C     ... COMPUTE PROFILE AND BANDWIDTH FOR BOTH THE OLD AND THE NEW
C         ORDERINGS.
C
      OLDBND = 0
      OLDPRO = 0
      NEWBND = 0
      NEWPRO = 0
C
      DO 300 I = 1, N
          IF (DEGREE(I) .EQ. 0)  GO TO 300
          IF (DEGREE(I) .GT. 0)  GO TO 6000
              IDGREE = -DEGREE(I)
              DEGREE(I) = IDGREE
              NEWORG = NEWNUM(I)
              OLDORG = OLDNUM(I)
              NEWRWD = 0
              OLDRWD = 0
              JPTR = RSTART (I)
C
C             ... FIND NEIGHBOR WHICH IS NUMBERED FARTHEST AHEAD OF THE
C                 CURRENT NODE.
C
              DO 200 J = 1, IDGREE
                  JNODE = CONNEC(JPTR)
                  JPTR = JPTR + 1
                  NEWRWD = MAX0 (NEWRWD, NEWORG - NEWNUM(JNODE))
                  OLDRWD = MAX0 (OLDRWD, OLDORG - OLDNUM(JNODE))
  200         CONTINUE
C
              NEWPRO = NEWPRO + NEWRWD
              NEWBND = MAX0 (NEWBND, NEWRWD)
              OLDPRO = OLDPRO + OLDRWD
              OLDBND = MAX0 (OLDBND, OLDRWD)
  300 CONTINUE
C
C     ... IF NEW ORDERING HAS BETTER BANDWIDTH THAN OLD ORDERING,
C         REPLACE OLD ORDERING BY NEW ORDERING
C
      IF  (NEWBND .GT. OLDBND)  GO TO 500
          BANDWD = NEWBND
          PROFIL = NEWPRO
          DO 400 I = 1, N
              OLDNUM(I) = NEWNUM(I)
  400     CONTINUE
          GO TO 600
C
C     ... RETAIN OLD ORDERING
C
  500     BANDWD = OLDBND
          PROFIL = OLDPRO
C
  600 RETURN
C
C     ------------------------------------------------------------------
C
 6000 SPACE = -1
      ERROR = 70
      RETURN
C
      END
      SUBROUTINE   GPSKCM   (N, DEGREE, RSTART, CONNEC, INVNUM, NEWNUM,
     1                       OLDNUM, BANDWD, PROFIL, ERROR, SPACE)
C
C
      INTEGER     N, RSTART(N), BANDWD, PROFIL, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(N), INVNUM(N), NEWNUM(N), OLDNUM(N)
      INTEGER     DEGREE(N), CONNEC(N), INVNUM(N), NEWNUM(N), OLDNUM(N)
C
C     ==================================================================
C
C
C     COMPUTE THE BANDWIDTH AND PROFILE FOR THE RENUMBERING GIVEN
C     BY 'INVNUM', BY THE REVERSE OF NUMBERING 'INVNUM', AND ALSO
C     BY THE RENUMBERING GIVEN IN 'OLDNUM'.
C     'NEWNUM' WILL BE A PERMUTATION VECTOR COPY OF THE NODE
C     LIST 'INVNUM'.
C
C     ==================================================================
C
      INTEGER     I, J, JPTR, IDGREE, OLDBND, OLDPRO, NEWBND, NEWPRO,
     1            OLDRWD, NEWRWD, OLDORG, NEWORG, JNODE, NRVBND, NRVPRO,
     2            NRVORG, NRVRWD, INVNMI, NMIP1
C
C     ------------------------------------------------------------------
C
C     ... CREATE NEWNUM AS A PERMUTATION VECTOR
C
      DO 100 I = 1, N
          INVNMI = INVNUM (I)
          NEWNUM (INVNMI) = I
  100 CONTINUE
C
C     ... COMPUTE PROFILE AND BANDWIDTH FOR BOTH THE OLD AND THE NEW
C         ORDERINGS.
C
      OLDBND = 0
      OLDPRO = 0
      NEWBND = 0
      NEWPRO = 0
      NRVBND = 0
      NRVPRO = 0
C
      DO 300 I = 1, N
          IF (DEGREE(I) .EQ. 0)  GO TO 300
          IF (DEGREE(I) .GT. 0)  GO TO 6000
              IDGREE = -DEGREE(I)
              DEGREE(I) = IDGREE
              NEWRWD = 0
              OLDRWD = 0
              NRVRWD = 0
              NEWORG = NEWNUM(I)
              OLDORG = OLDNUM(I)
              NRVORG = N - NEWNUM(I) + 1
              JPTR = RSTART (I)
C
C             ... FIND NEIGHBOR WHICH IS NUMBERED FARTHEST AHEAD OF THE
C                 CURRENT NODE.
C
              DO 200 J = 1, IDGREE
                  JNODE = CONNEC(JPTR)
                  JPTR = JPTR + 1
                  NEWRWD = MAX0 (NEWRWD, NEWORG - NEWNUM(JNODE))
                  OLDRWD = MAX0 (OLDRWD, OLDORG - OLDNUM(JNODE))
                  NRVRWD = MAX0 (NRVRWD, NRVORG - N + NEWNUM(JNODE) - 1)
  200         CONTINUE
C
              NEWPRO = NEWPRO + NEWRWD
              NEWBND = MAX0 (NEWBND, NEWRWD)
              NRVPRO = NRVPRO + NRVRWD
              NRVBND = MAX0 (NRVBND, NRVRWD)
              OLDPRO = OLDPRO + OLDRWD
              OLDBND = MAX0 (OLDBND, OLDRWD)
  300 CONTINUE
C
C     ... IF NEW ORDERING HAS BETTER BANDWIDTH THAN OLD ORDERING,
C         REPLACE OLD ORDERING BY NEW ORDERING
C
      IF  ((NEWPRO .GT. OLDPRO)  .OR. (NEWPRO .GT. NRVPRO)) GO TO 500
          BANDWD = NEWBND
          PROFIL = NEWPRO
          DO 400 I = 1, N
              OLDNUM(I) = NEWNUM(I)
  400     CONTINUE
          GO TO 800
C
C     ... CHECK NEW REVERSED ORDERING FOR BEST PROFILE
C
  500 IF  (NRVPRO .GT. OLDPRO)  GO TO 700
          BANDWD = NRVBND
          PROFIL = NRVPRO
          DO 600 I = 1, N
              OLDNUM(I) = N - NEWNUM(I) + 1
              IF  (I .GT. N/2)  GO TO 600
                  J = INVNUM(I)
                  NMIP1 = (N + 1) - I
                  INVNUM(I) = INVNUM (NMIP1)
                  INVNUM (NMIP1) = J
  600     CONTINUE
          GO TO 800
C
C
C     ... RETAIN OLD ORDERING
C
  700     BANDWD = OLDBND
          PROFIL = OLDPRO
C
  800 RETURN
C
C     ------------------------------------------------------------------
C
 6000 ERROR = 71
      SPACE = -1
      RETURN
C
      END
      SUBROUTINE   GPSKCN   (N, KEY, DATA, ERROR)
C
C     ==================================================================
C
C     I N S E R T I O N    S O R T
C
C     INPUT:
C         N    -- NUMBER OF ELEMENTS TO BE SORTED
C         KEY  -- AN ARRAY OF LENGTH  N  CONTAINING THE VALUES
C                 WHICH ARE TO BE SORTED
C         DATA -- A SECOND ARRAY OF LENGTH  N  CONTAINING DATA
C                 ASSOCIATED WITH THE INDIVIDUAL KEYS.
C
C     OUTPUT:
C         KEY  -- WILL BE ARRANGED SO THAT VALUES ARE IN DECREASING
C                 ORDER
C         DATA -- REARRANGED TO CORRESPOND TO REARRANGED KEYS
C         ERROR -- WILL BE ZERO UNLESS THE PROGRAM IS MALFUNCTIONING,
C                  IN WHICH CASE IT WILL BE EQUAL TO 1.
C
C
C     ==================================================================
C
      INTEGER     N, ERROR
C
CIBM  INTEGER *2  KEY(N), DATA(N)
      INTEGER     KEY(N), DATA(N)
C
C     ------------------------------------------------------------------
C
      INTEGER     I, J, D, K, IP1, JM1
C
C     ------------------------------------------------------------------
C
      IF (N .EQ. 1)  RETURN
      IF  (N .LE. 0)  GO TO 6000
C
      ERROR = 0
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2400 I = N - 1
      IP1 = N
C
 2500     IF ( KEY (I) .GE. KEY (IP1) )  GO TO 2800
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              K = KEY (I)
              D = DATA (I)
              J = IP1
              JM1 = I
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR K FOUND'
C
 2600             KEY (JM1) = KEY (J)
                  DATA (JM1) = DATA (J)
                  JM1 = J
                  J = J + 1
                  IF  (J .GT. N)  GO TO 2700
                  IF (KEY (J) .GT. K)  GO TO 2600
C
 2700         KEY (JM1) = K
              DATA (JM1) = D
C
 2800     IP1 = I
          I = I - 1
          IF ( I .GT. 0 )  GO TO 2500
C
 3000 RETURN
C
 6000 ERROR = 1
      GO TO 3000
C
      END
      SUBROUTINE   GPSKCO   (N, KEY, ERROR)
C
C     ==================================================================
C
C     I N S E R T I O N    S O R T
C
C     INPUT:
C         N    -- NUMBER OF ELEMENTS TO BE SORTED
C         KEY  -- AN ARRAY OF LENGTH  N  CONTAINING THE VALUES
C                 WHICH ARE TO BE SORTED
C
C     OUTPUT:
C         KEY  -- WILL BE ARRANGED SO THAT VALUES ARE IN DECREASING
C                 ORDER
C
C     ==================================================================
C
      INTEGER     N, ERROR
C
CIBM  INTEGER *2  KEY(N)
      INTEGER     KEY(N)
C
C     ------------------------------------------------------------------
C
      INTEGER     I, J, K, IP1, JM1
C
C     ------------------------------------------------------------------
C
      IF (N .EQ. 1)  RETURN
      IF  (N .LE. 0)  GO TO 6000
C
      ERROR = 0
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2400 I = N - 1
      IP1 = N
C
 2500     IF ( KEY (I) .GE. KEY (IP1) )  GO TO 2800
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              K = KEY (I)
              J = IP1
              JM1 = I
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR K FOUND'
C
 2600             KEY (JM1) = KEY (J)
                  JM1 = J
                  J = J + 1
                  IF  (J .GT. N)  GO TO 2700
                  IF (KEY (J) .GT. K)  GO TO 2600
C
 2700         KEY (JM1) = K
C
 2800     IP1 = I
          I = I - 1
          IF ( I .GT. 0 )  GO TO 2500
C
 3000 RETURN
C
 6000 ERROR = 1
      GO TO 3000
C
      END
      SUBROUTINE  GPSKCP  (N, INDEX, NVEC, DEGREE, ERROR)
C
C     ==================================================================
C
C     I N S E R T I O N      S O R T
C
C     INPUT:
C         N    -- NUMBER OF ELEMENTS TO BE SORTED
C         INDEX  -- AN ARRAY OF LENGTH  N  CONTAINING THE INDICES
C                 WHOSE DEGREES ARE TO BE SORTED
C         DEGREE -- AN  NVEC  VECTOR, GIVING THE DEGREES OF NODES
C                   WHICH ARE TO BE SORTED.
C
C     OUTPUT:
C         INDEX  -- WILL BE ARRANGED SO THAT VALUES ARE IN DECREASING
C                 ORDER
C         ERROR -- WILL BE ZERO UNLESS THE PROGRAM IS MALFUNCTIONING,
C                  IN WHICH CASE IT WILL BE EQUAL TO 1.
C
C     ==================================================================
C
      INTEGER     N, NVEC, ERROR
C
CIBM  INTEGER *2  INDEX(N), DEGREE(NVEC)
      INTEGER     INDEX(N), DEGREE(NVEC)
C
C     ------------------------------------------------------------------
C
      INTEGER     I, J, V, IP1, JM1, INDEXI, INDXI1, INDEXJ
C
C     ------------------------------------------------------------------
C
      IF (N .EQ. 1)  RETURN
      IF  (N .LE. 0)  GO TO 6000
C
      ERROR = 0
C
C     ------------------------------------------------------------------
C     INSERTION SORT THE ENTIRE FILE
C     ------------------------------------------------------------------
C
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2400 I = N - 1
      IP1 = N
C
 2500     INDEXI = INDEX (I)
          INDXI1 = INDEX (IP1)
          IF ( DEGREE(INDEXI) .GE. DEGREE(INDXI1) )  GO TO 2800
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              V = DEGREE (INDEXI)
              J = IP1
              JM1 = I
              INDEXJ = INDEX (J)
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR V FOUND'
C
 2600             INDEX (JM1) = INDEXJ
                  JM1 = J
                  J = J + 1
                  IF (J .GT. N)  GO TO 2700
                  INDEXJ = INDEX (J)
                  IF (DEGREE(INDEXJ) .GT. V)  GO TO 2600
C
 2700         INDEX (JM1) = INDEXI
C
 2800     IP1 = I
          I = I - 1
          IF ( I .GT. 0 )  GO TO 2500
C
 3000 RETURN
C
 6000 ERROR = 1
      GO TO 3000
C
      END
      SUBROUTINE  GPSKCQ (N, INDEX, NVEC, DEGREE, ERROR)
C
C     ==================================================================
C
C     I N S E R T I O N      S O R T
C
C     INPUT:
C         N    -- NUMBER OF ELEMENTS TO BE SORTED
C         INDEX  -- AN ARRAY OF LENGTH  N  CONTAINING THE INDICES
C                 WHOSE DEGREES ARE TO BE SORTED
C         DEGREE -- AN  NVEC  VECTOR, GIVING THE DEGREES OF NODES
C                   WHICH ARE TO BE SORTED.
C
C     OUTPUT:
C         INDEX  -- WILL BE ARRANGED SO THAT VALUES ARE IN INCREASING
C                   ORDER
C         ERROR -- WILL BE ZERO UNLESS THE PROGRAM IS MALFUNCTIONING,
C                  IN WHICH CASE IT WILL BE EQUAL TO 1.
C
C     ==================================================================
C
      INTEGER     N, NVEC, ERROR
C
CIBM  INTEGER *2  INDEX(N), DEGREE(NVEC)
      INTEGER     INDEX(N), DEGREE(NVEC)
C
C     ------------------------------------------------------------------
C
      INTEGER     I, J, V, INDEXI, INDXI1, INDEXJ, IP1, JM1
C
C     ------------------------------------------------------------------
C
      IF (N .EQ. 1)  RETURN
      IF  (N .LE. 0)  GO TO 6000
C
      ERROR = 0
C
C     ------------------------------------------------------------------
C     INSERTION SORT THE ENTIRE FILE
C     ------------------------------------------------------------------
C
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2400 I = N - 1
      IP1 = N
C
 2500     INDEXI = INDEX (I)
          INDXI1 = INDEX (IP1)
          IF ( DEGREE(INDEXI) .LE. DEGREE(INDXI1) )  GO TO 2800
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              V = DEGREE (INDEXI)
              J = IP1
              JM1 = I
              INDEXJ = INDEX (J)
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR V FOUND'
C
 2600             INDEX (JM1) = INDEXJ
                  JM1 = J
                  J = J + 1
                  IF (J .GT. N)  GO TO 2700
                  INDEXJ = INDEX (J)
                  IF (DEGREE(INDEXJ) .LT. V)  GO TO 2600
C
 2700         INDEX (JM1) = INDEXI
C
 2800     IP1 = I
          I = I - 1
          IF ( I .GT. 0 )  GO TO 2500
C
 3000 RETURN
C
 6000 ERROR = 1
      GO TO 3000
C
      END
