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

! This file contains dummy routines to satisfy ARPACK external references
! when ARPACK is not needed.  The code was lifted from the ARPACK source code.
! Added calls to warning (these routines should never be called),
! use message_passing (to get subroutine warning) and use of the arguments (to
! shut up picky compiler warnings).

      subroutine snaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Real
     &           tol
      integer    iparam(11), ipntr(14)
      Real
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine snaupd called.",
     &             "Need to link in an ARPACK library.")
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

      subroutine dnaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
      integer    iparam(11), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine dnaupd called.",
     &             "Need to link in an ARPACK library.")
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

      subroutine sneupd (rvec, howmny, select, dr, di, z, ldz, sigmar,
     &                   sigmai, workev, bmat, n, which, nev, tol,
     &                   resid, ncv, v, ldv, iparam, ipntr, workd,
     &                   workl, lworkl, info)
      use message_passing
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Real
     &           sigmar, sigmai, tol
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Real
     &           dr(nev+1), di(nev+1), resid(n), v(ldv,ncv), z(ldz,*),
     &           workd(3*n), workl(lworkl), workev(3*ncv)
      call warning("Dummy version of ARPACK routine sneupd called.",
     &             "Need to link in an ARPACK library.")
      if (rvec .and. select(1)) then
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

      subroutine dneupd (rvec, howmny, select, dr, di, z, ldz, sigmar,
     &                   sigmai, workev, bmat, n, which, nev, tol,
     &                   resid, ncv, v, ldv, iparam, ipntr, workd,
     &                   workl, lworkl, info)
      use message_passing
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision
     &           sigmar, sigmai, tol
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision
     &           dr(nev+1), di(nev+1), resid(n), v(ldv,ncv), z(ldz,*),
     &           workd(3*n), workl(lworkl), workev(3*ncv)
      call warning("Dummy version of ARPACK routine dneupd called.",
     &             "Need to link in an ARPACK library.")
      if (rvec .and. select(1)) then
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

      subroutine dsaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
c
      use message_passing
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
      integer    iparam(11), ipntr(11)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine dsaupd called.",
     &             "Need to link in an ARPACK library.")
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

      subroutine dseupd (rvec, howmny, select, d, z, ldz, sigma, bmat,
     &                   n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &                   ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision
     &           sigma, tol
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Double precision
     &           d(nev), resid(n), v(ldv,ncv), z(ldz, nev),
     &           workd(2*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine dneupd called.",
     &             "Need to link in an ARPACK library.")
      if (rvec .and. select(1)) then
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

      subroutine ssaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Real
     &           tol
      integer    iparam(11), ipntr(11)
      Real
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine ssaupd called.",
     &             "Need to link in an ARPACK library.")
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

      subroutine sseupd (rvec, howmny, select, d, z, ldz, sigma, bmat,
     &                   n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &                   ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Real
     &           sigma, tol
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Real
     &           d(nev), resid(n), v(ldv,ncv), z(ldz, nev),
     &           workd(2*n), workl(lworkl)
      call warning("Dummy version of ARPACK routine sseupd called.",
     &             "Need to link in an ARPACK library.")
      if (rvec .and. select(1)) then
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
