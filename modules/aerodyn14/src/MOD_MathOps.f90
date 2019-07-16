MODULE MathOps

  USE NWTC_Library
  USE NWTC_LAPACK

  IMPLICIT NONE

  CONTAINS
      !************************************
      ! Math operations used throughout FVW code
      ! Includes: cross, dot, rms, norm, and pinv
      !************************************
      ! Kelsey Shaler 8/28/14


!=================================================
!> Calculate the inverse of the square matrix A using single value decomposition
!! routines in the LAPACK library
SUBROUTINE Pinv(A, M, Ainv, ErrStat, ErrMsg)

   IMPLICIT NONE

   INTEGER,                intent(in   )  :: M
   REAL(ReKi),             intent(inout)  :: A(M,M)
   REAL(ReKi),             intent(  out)  :: Ainv(M,M)
   INTEGER(IntKi),         intent(  out)  :: ErrStat
   CHARACTER(ErrMsgLen),   intent(  out)  :: ErrMsg

   INTEGER(IntKi)                      :: lwork, lwmax
   INTEGER(IntKi)                      :: r, summation, i

   REAL( ReKi )                        :: tolerance
   REAL( ReKi ),  ALLOCATABLE          :: WORK(:)
   REAL( ReKi )                        :: S(M), U(M,M), VT(M,M), S_mat(M,M)

   INTEGER(IntKi)                      :: ErrStat2
   CHARACTER(ErrMsgLen)                :: ErrMsg2
   CHARACTER(*),  PARAMETER            :: RoutineName='Pinv'

   ErrStat = ErrID_None
   ErrMsg  = ""


!FIXME: To optimize, setup the work array as a miscvar, calculate optimal size at init and set
   !--------------------------
   !  Find size of work array.
   !--------------------------

      ! set the size of the work array to something that we know will be possible to use
   LWMAX = MIN(7*M,1000)
   ALLOCATE(WORK(LWMAX))
   work = 0.0_ReKi

      ! Query the [d,s]gesvd LAPACK routines to find out the optimal size for the work array.
   LWORK = -1
   call LAPACK_gesvd('A', 'A', M, M, A, S, U, Vt, work, lwork, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! If LAPACK (or MKL) suggested a larger work array as optimal, reallocate work array.
   if (LWMAX < work(1) ) then
      LWMAX=work(1)
      deallocate(work)
      allocate(work(LWMAX))
      work = 0.0_ReKi
   endif

   !--------------------------
   !  Compute SVD.
   !--------------------------

   call LAPACK_gesvd('A', 'A', M, M, A, S, U, Vt, work, size(work), ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! To speed up the calculations, find the tolerance and only calculate results
      ! for terms above the tolerance
   tolerance = M*epsilon(maxval(S))

   summation=0
   DO i=1,M
      IF (s(i) .GT. tolerance) THEN
         summation=summation+1;
      END IF
   END DO
   r=summation

      ! Set the diagonal elements of S_mat
   S_mat = 0.0_ReKi
   DO i = 1, M
      IF (i .LE. r)THEN
         S_mat(i,i)=1.0_ReKi/s(i)
      END IF
   END DO

      ! Calculate the inverse of A
   Ainv=transpose(matmul( matmul(U(:,1:r),S_mat(1:r,1:r)), VT(1:r,:)))

   DEALLOCATE(WORK)

END SUBROUTINE Pinv
!=================================================

END MODULE MathOps
