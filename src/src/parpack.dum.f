!---------------------------------------------------------------------!
!                                PHAML                                !
!                                                                     !
! The Parallel Hierarchical Adaptive MultiLevel code for solving      !
! linear elliptic partial differential equations of the form          !
! (PUx)x + (QUy)y + RU = F on 2D polygonal domains with mixed         !
! boundary conditions, and eigenvalue problems where F is lambda*U.   !
!                                                                     !
! PHAML is public domain software.  It was produced as part of work   !
! done by the U.S. Government, and is not subject to copyright in     !
! the United States.                                                  !
!                                                                     !
!     William F. Mitchell                                             !
!     Mathematical and Computational Sciences Division                !
!     National Institute of Standards and Technology                  !
!     william.mitchell@nist.gov                                       !
!     http://math.nist.gov/phaml                                      !
!                                                                     !
!---------------------------------------------------------------------!

! This file contains dummy routines to satisfy PARPACK external references
! when PARPACK is not needed.  The code was lifted from the PARPACK source code.
! Added calls to warning (these routines should never be called),
! use message_passing (to get subroutine warning) and use of the arguments (to
! shut up picky compiler warnings).

      subroutine psnaupd
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &     iparam, ipntr, workd, workl, lworkl, info )
      use message_passing
      integer    comm
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Real 
     &           tol
      integer    iparam(11), ipntr(14)
      Real 
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine psnaupd called.",
     &             "Need to use MPI and link in a PARPACK library.")
      info = comm
      info = ido
      info = ldv
      info = nev
      info = iparam(1)
      info = ipntr(1)
      info = ichar(bmat(1:1))
      info = ichar(which(1:1))
      info = tol
      info = resid(1)
      info = v(1,1)
      info = workd(1)
      info = workl(1)
      info = 0
      return
      end


      subroutine pdnaupd 
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &     iparam, ipntr, workd, workl, lworkl, info )
      use message_passing
      integer    comm
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision 
     &           tol
      integer    iparam(11), ipntr(14)
      Double precision 
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine pdnaupd called.",
     &             "Need to use MPI and link in a PARPACK library.")
      info = comm
      info = ido
      info = ldv
      info = nev
      info = iparam(1)
      info = ipntr(1)
      info = ichar(bmat(1:1))
      info = ichar(which(1:1))
      info = tol
      info = resid(1)
      info = v(1,1)
      info = workd(1)
      info = workl(1)
      info = 0
      return
      end

      subroutine psneupd 
     &         (comm , rvec , howmny, select, dr    , di  ,
     &          z    , ldz  , sigmar, sigmai, workev, bmat,
     &          n    , which, nev   , tol   , resid ,
     &          ncv  , v    , ldv   , iparam, ipntr ,
     &          workd, workl, lworkl, info  )
      use message_passing
      integer   comm
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Real      
     &           sigmar, sigmai, tol
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Real 
     &           dr(nev+1)    , di(nev+1)    , resid(n)  ,
     &           v(ldv,ncv)   , z(ldz,*)     , workd(3*n),
     &           workl(lworkl), workev(3*ncv)
      call warning("Dummy version of ARPACK routine psneupd called.",
     &             "Need to use MPI and link in a PARPACK library.")
      if (rvec .and. select(1)) then
         info = comm
         info = ichar(bmat(1:1))
         info = ichar(howmny(1:1))
         info = ichar(which(1:1))
         info = sigmar
         info = sigmai
         info = tol
         info = iparam(1)
         info = ipntr(1)
         info = dr(1)
         info = di(1)
         info = resid(1)
         info = v(1,1)
         info = z(1,1)
         info = workd(1)
         info = workl(1)
         info = workev(1)
         info = 0
      endif
      return
      end

      subroutine pdneupd  
     &         (comm , rvec , howmny, select, dr    , di  ,
     &          z    , ldz  , sigmar, sigmai, workev, bmat,
     &          n    , which, nev   , tol   , resid ,
     &          ncv  , v    , ldv   , iparam, ipntr ,
     &          workd, workl, lworkl, info  )
      use message_passing
      integer   comm
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision      
     &           sigmar, sigmai, tol
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision 
     &           dr(nev+1)    , di(nev+1)    , resid(n)  ,
     &           v(ldv,ncv)   , z(ldz,*)     , workd(3*n),
     &           workl(lworkl), workev(3*ncv)
      call warning("Dummy version of ARPACK routine pdneupd called.",
     &             "Need to use MPI and link in a PARPACK library.")
      if (rvec .and. select(1)) then
         info = comm
         info = ichar(bmat(1:1))
         info = ichar(howmny(1:1))
         info = ichar(which(1:1))
         info = sigmar
         info = sigmai
         info = tol
         info = iparam(1)
         info = ipntr(1)
         info = dr(1)
         info = di(1)
         info = resid(1)
         info = v(1,1)
         info = z(1,1)
         info = workd(1)
         info = workl(1)
         info = workev(1)
         info = 0
      endif
      return
      end
!---------------------------------------------------------------------!
!                                PHAML                                !
!                                                                     !
! The Parallel Hierarchical Adaptive MultiLevel code for solving      !
! linear elliptic partial differential equations of the form          !
! (PUx)x + (QUy)y + RU = F on 2D polygonal domains with mixed         !
! boundary conditions, and eigenvalue problems where F is lambda*U.   !
!                                                                     !
! PHAML is public domain software.  It was produced as part of work   !
! done by the U.S. Government, and is not subject to copyright in     !
! the United States.                                                  !
!                                                                     !
!     William F. Mitchell                                             !
!     Mathematical and Computational Sciences Division                !
!     National Institute of Standards and Technology                  !
!     william.mitchell@nist.gov                                       !
!     http://math.nist.gov/phaml                                      !
!                                                                     !
!---------------------------------------------------------------------!

      subroutine pdsaupd
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &     iparam, ipntr, workd, workl, lworkl, info )
c
      use message_passing
      character  bmat*1, which*2
      integer    comm, ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
      integer    iparam(11), ipntr(11)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine pdsaupd called.",
     &             "Need to link in an ARPACK library.")
      info = comm
      info = ido
      info = ldv
      info = nev
      info = iparam(1)
      info = ipntr(1)
      info = ichar(bmat(1:1))
      info = ichar(which(1:1))
      info = tol
      info = resid(1)
      info = v(1,1)
      info = workd(1)
      info = workl(1)
      info = 0
      return
      end

      subroutine pdseupd
     &    (comm, rvec, howmny, select, d, z, ldz, sigma, bmat,
     &     n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat, howmny, which*2
      logical    rvec
      integer    comm, info, ldz, ldv, lworkl, n, ncv, nev
      Double precision
     &           sigma, tol
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Double precision
     &           d(nev), resid(n), v(ldv,ncv), z(ldz, nev),
     &           workd(2*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine pdneupd called.",
     &             "Need to link in an ARPACK library.")
      if (rvec .and. select(1)) then
         info = comm
         info = ichar(bmat(1:1))
         info = ichar(howmny(1:1))
         info = ichar(which(1:1))
         info = sigma
         info = tol
         info = iparam(1)
         info = ipntr(1)
         info = d(1)
         info = resid(1)
         info = v(1,1)
         info = z(1,1)
         info = workd(1)
         info = workl(1)
         info = 0
      endif
      return
      end

      subroutine pssaupd
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &     iparam, ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat*1, which*2
      integer    comm, ido, info, ldv, lworkl, n, ncv, nev
      Real
     &           tol
      integer    iparam(11), ipntr(11)
      Real
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine pssaupd called.",
     &             "Need to link in an ARPACK library.")
      info = comm
      info = ido
      info = ldv
      info = nev
      info = iparam(1)
      info = ipntr(1)
      info = ichar(bmat(1:1))
      info = ichar(which(1:1))
      info = tol
      info = resid(1)
      info = v(1,1)
      info = workd(1)
      info = workl(1)
      info = 0
      return
      end

      subroutine psseupd
     &    (comm, rvec, howmny, select, d, z, ldz, sigma, bmat,
     &     n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat, howmny, which*2
      logical    rvec
      integer    comm, info, ldz, ldv, lworkl, n, ncv, nev
      Real
     &           sigma, tol
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Real
     &           d(nev), resid(n), v(ldv,ncv), z(ldz, nev),
     &           workd(2*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine psseupd called.",
     &             "Need to link in an ARPACK library.")
      if (rvec .and. select(1)) then
         info = comm
         info = ichar(bmat(1:1))
         info = ichar(howmny(1:1))
         info = ichar(which(1:1))
         info = sigma
         info = tol
         info = iparam(1)
         info = ipntr(1)
         info = d(1)
         info = resid(1)
         info = v(1,1)
         info = z(1,1)
         info = workd(1)
         info = workl(1)
         info = 0
      endif
      return
      end
