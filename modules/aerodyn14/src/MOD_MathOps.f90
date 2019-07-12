MODULE MathOps

  USE Precision

  IMPLICIT NONE

  CONTAINS
      !************************************
      ! Math operations used throughout FVW code
      ! Includes: cross, dot, rms, norm, and pinv
      !************************************
      ! Kelsey Shaler 8/28/14


!=================================================
SUBROUTINE RMS(r_new,r_old,RMSval)

  USE FVW_Parm, Only : CUTOFF_up, CUTOFF_allocate
  USE MultTurb_Params, Only: NumWakes, NTurb

  IMPLICIT NONE

  INTEGER                                       :: k,n,indx,counter
  REAL( ReKi )                                        :: RMSval, summation
  REAL( ReKi ), DIMENSION( 3, CUTOFF_allocate, NumWakes ) :: r_new, r_old

  summation = 0.0_ReKi
  counter = 0
  NTurb = 1  !! KS -- HACKHACKHACK 7.1.19
  DO n = 1, NumWakes
     DO k = 1, CUTOFF_up( NTurb )
        DO indx = 1, 3
           counter = counter + 1
           summation = summation + ( r_new( indx, k, n ) - r_old( indx, k, n ))**2.0_ReKi
        END DO
     END DO
  END DO

  RMSval = sqrt( summation ) / REAL(counter,ReKi)


END SUBROUTINE RMS
!=================================================


!=================================================
SUBROUTINE Pinv(A,SIZEMAT,Ainv)

  IMPLICIT NONE

  INTEGER                                :: N, M, lwork, info, lda, ldvt, ldu, lwmax
  INTEGER                                :: r, summation, i, j, SIZEMAT

  REAL( ReKi )                                :: tolerance
  REAL( ReKi ), DIMENSION(:), ALLOCATABLE     :: S, WORK
  REAL( ReKi ), DIMENSION(SIZEMAT,SIZEMAT)   :: A, Ainv
  REAL( ReKi ), DIMENSION(:,:), ALLOCATABLE   :: Aold, U, VT, S_mat

  ALLOCATE(Aold(SIZEMAT,SIZEMAT),U(SIZEMAT,SIZEMAT),VT(SIZEMAT,SIZEMAT),S_mat(SIZEMAT,SIZEMAT),&
     &S(SIZEMAT),WORK(1000))

  LWMAX = 1000; M = SIZEMAT; N = SIZEMAT; LDA = M; LDVT = N; LDU = M; Aold = A;

  LWORK = -1
  PRINT*, 'SIZEMAT: ',SIZEMAT
  !!KS -- had to change LWORK b/c newer version of MLK  7.1.19
  !!https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/402473

  CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,WORK, LWORK, INFO )
  PRINT*, 'INFO=',INFO,'LWork asked for by MKL = ', WORK(1)
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
  PRINT*, 'LWORK: ',LWORK
  !LWORK = MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
  !
  !     Compute SVD.
  !

  CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
   !PRINT*, '!!!'
   !PRINT*, 'LWORK = ',LWORK, 'INFO = ',INFO

  tolerance=maxval(shape(A))*epsilon(maxval(S))

  summation=0
  DO i=1,SIZEMAT
     IF(s(i) .GT. tolerance)THEN
        summation=summation+1;
     END IF
  END DO
  r=summation

  DO j = 1, SIZEMAT
     DO i = 1, SIZEMAT
        IF (i .EQ. j .AND. i .LE. r)THEN
           S_mat(i,j)=1.0_ReKi/s(i)
        ELSE
           S_mat(i,j)=0.0_ReKi
        END IF
     END DO

  END DO

  Ainv=transpose(matmul(matmul(U(:,1:r),S_mat(1:r,1:r)),VT(1:r,:)))
  A=Aold

  DEALLOCATE(Aold,U,VT,S_mat,S,WORK)

END SUBROUTINE Pinv
!=================================================

END MODULE MathOps
