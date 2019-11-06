!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
MODULE BeamDyn_Subs

   USE BeamDyn_Types
   USE NWTC_LAPACK

   IMPLICIT NONE

   INTEGER, PARAMETER :: BDKi = R8Ki
   INTEGER, PARAMETER :: FLAG_R1R2   = 0 !<   BD_CrvCompose flag = 0: R(rr) = R    (pp) R    (qq)
   INTEGER, PARAMETER :: FLAG_R1TR2  = 1 !<   BD_CrvCompose flag = 1: R(rr) = R(T) (pp) R    (qq)
   INTEGER, PARAMETER :: FLAG_R1R2T  = 2 !<   BD_CrvCompose flag = 2: R(rr) = R    (pp) R(T) (qq)
   INTEGER, PARAMETER :: FLAG_R1TR2T = 3 !<   BD_CrvCompose flag = 3: R(rr) = R(T) (pp) R(T) (qq)

   INTEGER, PARAMETER :: GAUSS_QUADRATURE = 1 !<   p%quadrature method: gaussian quadrature
   INTEGER, PARAMETER :: TRAP_QUADRATURE  = 2 !<   p%quadrature method: trapeziodal quadrature

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine determines the (N+1) Gauss-Lobatto-Legendre points x and weights w
!!
!! For details, see
!! \@book{Deville-etal:2002,
!!  author =    {M. O. Deville and P. F. Fischer and E. H. Mund},
!!  title =     {High-Order Methods for Incompressible Fluid Flow},
!!  publisher = {Cambridge University Press},
!!  address = {Cambridge},
!!  year =      2002
!!}
SUBROUTINE BD_GenerateGLL(N1, GLL_nodes, ErrStat, ErrMsg)

   ! input variables

   INTEGER(IntKi),          INTENT(IN   ) :: N1             !< 1 + the order of spectral element, or equivalently p%nodes_per_elem
   REAL(BDKi), ALLOCATABLE, INTENT(  OUT) :: GLL_nodes(:)   !< location of GLL nodes
   INTEGER(IntKi),          INTENT(  OUT) :: ErrStat        !< Error status of the operation
   CHARACTER(*),            INTENT(  OUT) :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   REAL(BDKi)      , PARAMETER            :: tol   = 10.0_BDKi*EPSILON(tol) / 2.0_BDKi   ! tolerance for newton-raphson solve (ignores 1 significant digit)
   INTEGER(IntKi)  , PARAMETER            :: maxit = 1000                                ! maximum allowable iterations in newton-raphson solve

   REAL(BDKi)                             :: x_it        ! current NR-iteration value
   REAL(BDKi)                             :: x_old       ! last NR-iteration value
   REAL(BDKi)                             :: dleg(N1)    ! legendre polynomial
   INTEGER(IntKi)                         :: N           ! Order of spectral element
   INTEGER(IntKi)                         :: i           ! do-loop counter
   INTEGER(IntKi)                         :: j           ! do-loop counter
   INTEGER(IntKi)                         :: k           ! do-loop counter

   INTEGER(IntKi)                         :: ErrStat2
   CHARACTER(ErrMsgLen)                   :: ErrMsg2
   character(*), parameter                :: RoutineName = 'BD_GenerateGLL'


   ErrStat = ErrID_None
   ErrMsg  = ""

   N = N1 - 1


   CALL AllocAry(GLL_nodes,N1,'GLL points array',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return


   !bjj: these are getting overwritten later....
   ! enter known endpoints  [-1.0, 1.0]
   GLL_nodes(1) = -1.0_BDKi
   GLL_nodes(N1) = 1.0_BDKi

   DO i = 1, N1

      x_it = -COS(pi_D * REAL(i-1,BDKi) / N) ! initial guess - chebyshev points
      DO j = 1, maxit
         x_old = x_it
         dleg(1) = 1.0_BDKi
         dleg(2) = x_it
         DO k = 2,N
            dleg(k+1) = (  REAL(2*k - 1,BDKi) * dleg(k) * x_it &
                         - REAL(  k - 1,BDKi) * dleg(k-1) ) / REAL(k,BDKi)
         ENDDO
         x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                       (REAL(N1,BDKi) * dleg(N1) )
         IF (ABS(x_it - x_old) .lt. tol) EXIT
      ENDDO

      GLL_nodes(i) = x_it
     !GLL_weights(i) = 2.0_BDKi / (REAL(N * N1, BDKi) * dleg(N1)**2 )

   ENDDO

END SUBROUTINE BD_GenerateGLL
!-----------------------------------------------------------------------------------------------------------------------------------

!> This subroutine computes the rotation tensor (RT=\f$ \phi \f$ )
!! given Wiener-Milenkovic rotation parameters, \f$ \vec{c} \f$,
!! where \f$ \vec{c}=4\tan\left(\frac{\phi}{4}\right)\bar{n} \f$.
SUBROUTINE BD_CrvMatrixR(cc,Rr)
!FIXME: would this routine benefit from matrix mathematics?  The basic formulation is matrix based.

   REAL(BDKi),    INTENT(IN   ):: cc(3)         !< \f$ \vec{c}\f$
   REAL(BDKi),    INTENT(  OUT):: Rr(3,3)       !< returned rotation tensor (transpose of DCM matrix)

   REAL(BDKi)                  :: c0            !< \f$ c_0 = ( 1 - c_1 c_1 - c_2 c_2 - c_3 c_3 ) / 2 \f$ (different from AIAA paper)
   REAL(BDKi)                  :: c1
   REAL(BDKi)                  :: c2
   REAL(BDKi)                  :: c3
   REAL(BDKi)                  :: tr0           !< \f$ 2 / ( 1 - c_0 )^2 \f$ (which is different from AIAA paper)


   !> Compare this routine with equation (14) from AIAA paper, "Geometric Nonlinear Analysis of Composite Beams Using
   !! Wiener-Milenkovic Parameters", Wang, et. al.\n
   !! \f$
   !!  \underline{\underline{R}} (\vec{c}) =    t_{r0}
   !! \begin{bmatrix}
   !! c_0^2 + c_1^2 - c_2^2 - c_3^2    & 2(c_1 c_2 + c_0 c_3)           & 2(c_1 c_3 - c_0 c_2)           \\
   !! 2(c_1 c_2 - c_0 c_3)             & c_0^2 - c_1^2 + c_2^2 - c_3^2  & 2(c_2 c_3 + c_0 c_1)           \\
   !! 2(c_1 c_3 + c_0 c_2)             & 2(c_2 c_3 - c_0 c_1)           & c_0^2 - c_1^2 - c_2^2 + c_3^2  \\
   !! \end{bmatrix}
   !! \f$
   !!
   !! Where \f$ c_0 = 2 - c_1 c_1 / 8 - c_2 c_2 / 8 - c_3 c_3 /8 \f$ and \f$ t_{r0} = 1/(4-c_0)^2 \f$
   !!
   !! _Note:_  The above equation is the transpose of what we have here.  This subroutine returns the transpose of the DCM,
   !!          not the DCM.
   !!
   !! Also, note that the formulation in this routine looks different, but can be shown as algebraically identical.

   c1  = cc(1)/4.0_BDKi
   c2  = cc(2)/4.0_BDKi
   c3  = cc(3)/4.0_BDKi
      ! \f$ c_0 \f$ is equivalent to \f$ 0.5 - (||\vec{c}||)^2 / 32 \f$
   c0  = 0.5_BDKi * (1.0_BDKi-c1*c1-c2*c2-c3*c3)         ! 1/4 the value of the AIAA paper (after plugging in c1, c2, c3 conversions)


   tr0 = 1.0_BDKi - c0                          ! This is 1/4 the value of the AIAA paper, after converting c0.
   tr0 = 2.0_BDKi/(tr0*tr0)                     ! this is 32x the equivalent term from the AIAA paper.   This is well behaved and won't go to zero.

      ! The following terms can be shown to match the transpose of the DCM given in the AIAA paper.
   Rr(1,1) = tr0*(c1*c1 + c0*c0) - 1.0_BDKi
   Rr(2,1) = tr0*(c1*c2 + c0*c3)
   Rr(3,1) = tr0*(c1*c3 - c0*c2)

   Rr(1,2) = tr0*(c1*c2 - c0*c3)
   Rr(2,2) = tr0*(c2*c2 + c0*c0) - 1.0_BDKi
   Rr(3,2) = tr0*(c2*c3 + c0*c1)

   Rr(1,3) = tr0*(c1*c3 + c0*c2)
   Rr(2,3) = tr0*(c2*c3 - c0*c1)
   Rr(3,3) = tr0*(c3*c3 + c0*c0) - 1.0_BDKi

END SUBROUTINE BD_CrvMatrixR
!-----------------------------------------------------------------------------------------------------------------------------------

!> This subroutine calculutates the tangent vector, \f$ \underline{\underline{H}} \f$ defined in equation (35) of Wang, et. al.,
!! "BeamDyn: A hig-fidelity wind turbine blade solver in the FAST modular framework", Wind Energy, 2016.
!!
!! Note that the formulation here is slightly different than in the WE paper, but is algebraically equivalent.
SUBROUTINE BD_CrvMatrixH(cc,Hh)
!FIXME: would this routine benefit from matrix mathematics?  The basic formulation is matrix based.

   REAL(BDKi),INTENT(IN) ::cc(3)
   REAL(BDKi),INTENT(OUT)::Hh(3,3)

   REAL(BDKi):: cf1,cf2,cf3,cq,ocq,aa,cb0,cb1,cb2,cb3

      ! Note that the factor of 4 here is a rescaling that simplifies the equation implementation slightly.
   cf1 = cc(1)/4.0_BDKi
   cf2 = cc(2)/4.0_BDKi
   cf3 = cc(3)/4.0_BDKi
   cq  = cf1 * cf1 + cf2 * cf2 + cf3 * cf3   ! equivalent to \f$ (\vec{c} \dotproduct \vec{c}) / 16 \f$
   ocq = 1.0_BDKi + cq
   aa  = 2.0_BDKi * ocq * ocq
   cb0 = 2.0_BDKi * (1.0_BDKi - cq) / aa
   cb1 = cc(1)/aa
   cb2 = cc(2)/aa
   cb3 = cc(3)/aa

      !> The equation given in the WE paper is given as \n
      !! \f$
      !!    \underline{\underline{H}(\vec{c})} = 2/(4-c_0)^2 \left[ c_0 + \tilde{c} + \frac{1}{4} \underline{\vec{c}\vec{c}^T} \right]
      !! \f$
      !! where \f$ c_0 = 2 - \frac{1}{8} \underline{\vec{c}}^T \underline{\vec{c}} \f$
      !!
      !! This gives the matrix
      !!
      !! \f$
      !!  \underline{\underline{R} (\vec{c}) } =    \frac{2}{(4-c_0)^2}
      !!    \begin{bmatrix}
      !!       c_1 c_1 + c_0        &  c_1 c_2 - c_3     &  c_1 c_3 + c_2  \\
      !!       c_2 c_1 + c_3        &  c_2 c_2 + c_0     &  c_2 c_3 - c_1  \\
      !!       c_3 c_1 - c_2        &  c_3 c_2 + c_1     &  c_3 c_3 + c_0  \\
      !!    \end{bmatrix}
      !!    = 2/(4-c_0)^2 \left[ c_0 + \tilde{c} + \vec{c} \vec{c}^T / 4 \right]
      !! \f$\n
      !! where \f$ c_0 = 2 - \frac{1}{8} \vec{c}^T \vec{c} \f$


   Hh(1,1) = cb1 * cf1 + cb0
   Hh(2,1) = cb2 * cf1 + cb3
   Hh(3,1) = cb3 * cf1 - cb2

   Hh(1,2) = cb1 * cf2 - cb3
   Hh(2,2) = cb2 * cf2 + cb0
   Hh(3,2) = cb3 * cf2 + cb1

   Hh(1,3) = cb1 * cf3 + cb2
   Hh(2,3) = cb2 * cf3 - cb1
   Hh(3,3) = cb3 * cf3 + cb0


END SUBROUTINE BD_CrvMatrixH
!-----------------------------------------------------------------------------------------------------------------------------------
!>   This subroutine composes two Wiener-Milenkovic parameters pp and qq to find the resulting parameter rr
!!   This method is detailed in the paper: Bauchau, O.A., 2008, "Interpolation of finite rotations in flexible
!!   multi-body dynamics simulations", IMechE, Equation (9).
!!   see http://www.dymoresolutions.com/publications/BauchauEppleHeo08.pdf \n
!!   FLAG_R1R2:   flag = 0: R(rr) = R    (pp) R    (qq) \n
!!   FLAG_R1TR2:  flag = 1: R(rr) = R(T) (pp) R    (qq) \n
!!   FLAG_R1R2T:  flag = 2: R(rr) = R    (pp) R(T) (qq) \n
!!   FLAG_R1TR2T: flag = 3: R(rr) = R(T) (pp) R(T) (qq)
!!
SUBROUTINE BD_CrvCompose( rr, pp, qq, flag)
!FIXME: would this routine benefit from matrix mathematics?  The basic formulation is matrix based.

   REAL(BDKi),    INTENT(  OUT):: rr(3)     !< Composed rotation
   REAL(BDKi),    INTENT(IN   ):: pp(3)     !< Input rotation 1
   REAL(BDKi),    INTENT(IN   ):: qq(3)     !< Input rotation 2
   INTEGER       ,INTENT(IN   ):: flag      !< Option flag

   REAL(BDKi)                  :: pp0
   REAL(BDKi)                  :: p(3)
   REAL(BDKi)                  :: qq0
   REAL(BDKi)                  :: q(3)
   REAL(BDKi)                  :: tr1
   REAL(BDKi)                  :: Delta1
   REAL(BDKi)                  :: Delta2
   REAL(BDKi)                  :: dd1
   REAL(BDKi)                  :: dd2


      ! Set the local values pp and qq allowing for the transpose

   IF(flag==FLAG_R1TR2 .OR. flag==FLAG_R1TR2T) THEN ! "transpose" (negative) of first rotation parameter
       p = -pp
   ELSE
       p = pp
   ENDIF

   IF(flag==FLAG_R1R2T .OR. flag==FLAG_R1TR2T) THEN ! "transpose" (negative) of second rotation parameter
       q = -qq
   ELSE
       q = qq
   ENDIF

      !> ## Composing the resulting Wiener-Milenkovic parameter
      !!
      !! From the paper by Bauchau, the resulting parameter \f$ \underline{r} \f$ (denoted as rr in the code)
      !! can be composed from \f$ \underline{p} \f$ and \f$ \underline{q} \f$ (pp and qq, respectively, in the code) as:
      !!
      !!    \f$ \underline{r} = 4 \frac{ q_0 \underline{p} + p_0 \underline{q} + \tilde{p} \underline{q} }
      !!                               { \Delta_1 + \Delta_2 }          \f$
      !!
      !! where
      !!
      !!    \f$ p_0        =  2 - \frac{1}{8} \underline{p}^T \underline{p}
      !!                   = 2 - \frac{1}{8} \left( p_1 p_1 + p_2 p_2 + p_3 p_3 \right)   \f$ \n
      !!    \f$ q_0        =  2 - \frac{1}{8} \underline{q}^T \underline{q}
      !!                   = 2 - \frac{1}{8} \left( q_1 q_1 + q_2 q_2 + q_3 q_3 \right)   \f$ \n
      !!    \f$ \Delta_1   =  \left( 4 - p_0 \right) \left( 4 - q_0 \right)               \f$ \n
      !!    \f$ \Delta_2   =  p_0 q_0  - \underline{p}^T \underline{q}
      !!                   =  p_0 q_0 - \left( p_1 q_1 + p_2 q_2 + p_3 q_3 \right)        \f$ \n
      !!
      !!


      ! Calculate pp0 and qq0. See Bauchau for the mathematics here (equations 8 to 9 and interviening text)

   pp0 = 2.0_BDKi - dot_product(p,p) / 8.0_BDKi   ! p_0
   qq0 = 2.0_BDKi - dot_product(q,q) / 8.0_BDKi   ! q_0

   Delta1 = (4.0_BDKi - pp0) * (4.0_BDKi - qq0)   ! Delta_1 in Bauchau
   Delta2 = pp0 * qq0 - dot_product(p,q)          ! Delta_2 in Bauchau
   dd1 = Delta1 + Delta2                          ! Denomimator term for \Delta_2 >= 0
   dd2 = Delta1 - Delta2                          ! Denomimator term for \Delta_2 <  0

      ! Rescaling to remove singularities at +/- 2 \pi
      ! Note: changed this to test on \Delta_2 (instead of dd1 > dd2) for better consistency with documentation.
   IF ( Delta2 >= 0.0_BDKi ) THEN
       tr1 = 4.0_BDKi / dd1
   ELSE
       tr1 = -4.0_BDKi / dd2
   ENDIF

   rr = tr1 * (qq0*p + pp0*q + cross_product(p,q))


END SUBROUTINE BD_CrvCompose


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the CRV parameters given the rotation matrix
!> The algorithm for this subroutine can be found in Markley, 'Unit Quaternion from Rotation Matrix'
!> https://doi.org/10.2514/1.31730
SUBROUTINE BD_CrvExtractCrv(R, cc, ErrStat, ErrMsg)

   REAL(BDKi),       INTENT(IN   )  :: R(3,3)        !< Rotation Matrix
   REAL(BDKi),       INTENT(  OUT)  :: cc(3)         !< Crv paramters
   INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat       !< Error status of the operation
   CHARACTER(*),     INTENT(  OUT)  :: ErrMsg        !< Error message if ErrStat /= ErrID_None

   REAL(BDKi)                  :: pivot(4) ! Trace of the rotation matrix and diagonal elements
   REAL(BDKi)                  :: sm(0:3)
   REAL(BDKi)                  :: em
   REAL(BDKi)                  :: Rr(3,3)
   INTEGER                     :: i        ! case indicator

   INTEGER(IntKi)              :: ErrStat2 ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2  ! Temporary Error message
   character(*), parameter     :: RoutineName = 'BD_CrvExtractCrv'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! use the local rotation matrix variable to avoid side effects
   Rr = R

   !> Starting with equation (14) from AIAA paper, "Geometric Nonlinear Analysis of Composite Beams Using
   !! Wiener-Milenkovic Parameters", Wang, et. al. \n
   !! \f$
   !!  \underline{\underline{R}} (\vec{c}) =    t_{r0}
   !! \begin{bmatrix}
   !! c_0^2 + c_1^2 - c_2^2 - c_3^2    & 2(c_1 c_2 + c_0 c_3)           & 2(c_1 c_3 - c_0 c_2)           \\
   !! 2(c_1 c_2 - c_0 c_3)             & c_0^2 - c_1^2 + c_2^2 - c_3^2  & 2(c_2 c_3 + c_0 c_1)           \\
   !! 2(c_1 c_3 + c_0 c_2)             & 2(c_2 c_3 - c_0 c_1)           & c_0^2 - c_1^2 - c_2^2 + c_3^2  \\
   !! \end{bmatrix}
   !! \f$
   !!
   !! where \f$ c_0 = 2 - c_1 c_1 / 8 - c_2 c_2 / 8 - c_3 c_3 /8 \f$ and \f$ t_{r0} = 1/(4-c_0)^2 \f$
   !!
   !! _Note:_ The above equation does not match what is in the March 2016 BD manual (it is the transpose).
   !! It does, however, match equation 5.17 in the March 2016 "BeamDyn User's Guide and Theory Manual"

   ! mjs--determine whether R is a valid rotation matrix and correct it if not
   call BD_CheckRotMat(Rr, ErrStat2, ErrMsg2)
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! mjs--find max value of T := Tr(Rr) and diagonal elements of Rr
   ! This tells us which denominator is largest (and less likely to produce numerical noise)
   pivot = (/ Rr(1,1) + Rr(2,2) + Rr(3,3), Rr(1,1), Rr(2,2), Rr(3,3)  /)
   i = maxloc(pivot, 1) - 1 ! our sm array starts at 0, so we need to subtract 1 here to get the correct index

      !> - Condition 1: \f$ \underline{\underline{R(3,3)}} > T \f$
      !! This implies that \f$ c_3^2 > c_0^2 \f$
   select case (i)
   case (3)
      !! We need the condition that c_3 is not zero for this set of equations.
      sm(0)  = Rr(2,1) - Rr(1,2)                           !  4 c_0 c_3 t_{r0}
      sm(1)  = Rr(1,3) + Rr(3,1)                           !  4 c_1 c_3 t_{r0}
      sm(2)  = Rr(2,3) + Rr(3,2)                           !  4 c_2 c_3 t_{r0}
      sm(3)  = 1.0_BDKi - Rr(1,1) - Rr(2,2) + Rr(3,3)      !  4 c_3 c_3 t_{r0}

      !> - Condition 2: \f$ \underline{\underline{R(2,2)}} > T \f$
      !! This implies that \f$ c_2^2 > c_0^2 \f$
   case (2)
      !! We need the condition that c_2 is not zero for this set of equations.
      sm(0)  = Rr(1,3) - Rr(3,1)                           !  4 c_0 c_2 t_{r0}
      sm(1)  = Rr(1,2) + Rr(2,1)                           !  4 c_1 c_2 t_{r0}
      sm(2)  = 1.0_BDKi - Rr(1,1) + Rr(2,2) - Rr(3,3)      !  4 c_2 c_2 t_{r0}
      sm(3)  = Rr(2,3) + Rr(3,2)                           !  4 c_3 c_2 t_{r0}

      !> - Condition 3: \f$ \underline{\underline{R(1,1)}} > T \f$
      !! This implies that \f$ c_1^2 > c_0^2 \f$
   case (1)
      !! We need the condition that c_1 is not zero for this set of equations.
      sm(0)  = Rr(3,2) - Rr(2,3)                           !  4 c_0 c_1 t_{r0}
      sm(1)  = 1.0_BDKi + Rr(1,1) - Rr(2,2) - Rr(3,3)      !  4 c_1 c_1 t_{r0}
      sm(2)  = Rr(1,2) + Rr(2,1)                           !  4 c_2 c_1 t_{r0}
      sm(3)  = Rr(1,3) + Rr(3,1)                           !  4 c_3 c_1 t_{r0}

      !> - Condition 4: all diagonal terms are less than the trace
   case (0)
      !! We need the condition that c_0 is not zero for this set of equations.
      sm(0)  = 1.0_BDKi + Rr(1,1) + Rr(2,2) + Rr(3,3)      !  4 c_0 c_0 t_{r0}
      sm(1)  = Rr(3,2) - Rr(2,3)                           !  4 c_1 c_0 t_{r0}
      sm(2)  = Rr(1,3) - Rr(3,1)                           !  4 c_2 c_0 t_{r0}
      sm(3)  = Rr(2,1) - Rr(1,2)                           !  4 c_3 c_0 t_{r0}
   end select

   em = sm(0) + SIGN( 2.0_BDKi*SQRT(sm(i)), sm(0) )
   em = 4.0_BDKi/em                                        ! 1 / ( 4 t_{r0} c_{i} ), assuming 0 <= c_0 < 4 and c_{i} > 0
   cc = em*sm(1:3)

END SUBROUTINE BD_CrvExtractCrv


SUBROUTINE BD_CheckRotMat(R, ErrStat, ErrMsg)
   !> This subroutine checks for rotation matrix validity.
   !> Returns:
   !>   ErrStat = 0 if valid
   !>   ErrStat = 4 (fatal error) if invalid

   REAL(BDKi),       INTENT(IN   )  :: R(3,3)       !< Rotation Matrix
   INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat      !< Error status of the operation
   CHARACTER(*),     INTENT(  OUT)  :: ErrMsg       !< Error message if ErrStat /= ErrID_None
   REAL(BDKi)                       :: Rr(3,3)      !< Local Rotation Matrix variable
   INTEGER(IntKi),   PARAMETER      :: lwork = 27   !mjs--from LAPACK: dgesvd doc page, lwork >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N))
   REAL(BDKi)                       :: work(lwork)  ! where M x N is dimension of R, and lwork is the dimension of work
   REAL(BDKi)                       :: S(3), U(3,3), VT(3,3) !mjs--these three are the SVD matrices (S is actually a vector)
   INTEGER(IntKi)                   :: ErrStat2     ! Temporary Error status
   CHARACTER(ErrMsgLen)             :: ErrMsg2      ! Temporary Error message
   LOGICAL                          :: ortho        !mjs--logical value indicating whether R is orthogonal
   INTEGER                          :: i
   character(*), parameter          :: RoutineName = 'BD_CheckRotMat'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! use the local rotation matrix variable to avoid side effects
   Rr = R

   ! mjs--Start by determining if R is a valid rotation matrix using the properties:
   ! 1) the eigenvalues of an orthogonal matrix have complex modulus == 1, where
   !    the leading eigenvalue is +1 and the other two are a complex conjugate pair
   ! 2) a valid rotation matrix must have determinant == +1 i.e., the singular values == 1

   call LAPACK_gesvd('A', 'A', 3, 3, Rr, S, U, VT, work, lwork, ErrStat2, ErrMsg2)
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! mjs--If \f$ \underline{\underline{R}} \f$ is not a valid roatation tensor,
   !    and the correction is desired,
   !    compute \f$ \underline{\underline{R_{out}} \f$, the nearest orthogonal tensor
   !    to \f$ \underline{\underline{R}} \f$.
   !    This is done via computing SVD for \f$ \underline{\underline{R}} = USV^T \f$
   !    and setting \f$ \underline{\underline{R_{out}} = UV^T \f$
   !    otherwise, assign \f$ \underline{\underline{R_{out}}}  = \underline{\underline{R}} \f$

   do i = 1, 3
      ortho = equalrealnos(S(i), 1.0_BDKi)
      if (.not. ortho) then
         CALL SetErrStat(ErrID_Fatal, "Passed invalid rotation matrix", ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
      end if
   end do

   ! mjs--after consulting with Mike Sprague, it was decided that instead of fixing the rotation matrix and
   ! notifying the user, the simulation should be stopped if an invalid rotation matrix is passed
   ! To change this and implement the fix, use the following lines
   ! ErrStat2 = ErrID_Info
   ! ErrMsg2 = 'Passed invalid rotation matrix--fixing via SVD'
   ! CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ! if (ErrStat >= AbortErrLev) return
   ! R = matmul(U, VT)

END SUBROUTINE BD_CheckRotMat


!-----------------------------------------------------------------------------------------------------------------------------------
!> calculate Lagrangian interpolant tensor at ns points where basis
!! functions are assumed to be associated with (np+1) GLL points on [-1,1]
SUBROUTINE BD_diffmtc( nodes_per_elem,GLL_nodes,QPtN,nqp,Shp,ShpDer )

   ! See Bauchau equations 17.1 - 17.5

   INTEGER(IntKi),         INTENT(IN   )  :: nodes_per_elem !< Nodes per elemenent
   REAL(BDKi),             INTENT(IN   )  :: GLL_nodes(:)   !< GLL_nodes(p%nodes_per_elem): location of the (p%nodes_per_elem) p%GLL points
   REAL(BDKi),             INTENT(IN   )  :: QPtN(:)        !< Locations of quadrature points ([-1 1])
   INTEGER(IntKi),         INTENT(IN   )  :: nqp            !< number of quadrature points to consider. Should be size of 2nd index of Shp & ShpDer
   REAL(BDKi),             INTENT(INOUT)  :: Shp(:,:)       !< p%Shp    (or another Shp array for when we add outputs at arbitrary locations)
   REAL(BDKi), OPTIONAL,   INTENT(INOUT)  :: ShpDer(:,:)    !< p%ShpDer (or another Shp array for when we add outputs at arbitrary locations)

   REAL(BDKi)                  :: dnum
   REAL(BDKi)                  :: den
   REAL(BDKi),        PARAMETER:: eps = SQRT(EPSILON(eps)) !1.0D-08
   INTEGER(IntKi)              :: l
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: k

   ! See Bauchau equations 17.1 - 17.5

   if (present(ShpDer)) then
      ShpDer(:,:)  = 0.0_BDKi
      do j = 1,nqp
         do l = 1,nodes_per_elem
 
          if (EqualRealNos(QPtN(j),1.0_BDKi).AND.(l.EQ.nodes_per_elem)) then
            ShpDer(l,j) = REAL((nodes_per_elem)*(nodes_per_elem-1), BDKi)/4.0_BDKi
          elseif (EqualRealNos(QPtN(j),1.0_BDKi).AND.(l.EQ.1)) then
            ShpDer(l,j) = -REAL((nodes_per_elem)*(nodes_per_elem-1), BDKi)/4.0_BDKi
          elseif (EqualRealNos(QPtN(j),GLL_nodes(l))) then
            ShpDer(l,j) = 0.0_BDKi
          else
            ShpDer(l,j) = 0.0_BDKi
            den = 1.0_BDKi
            do i = 1,nodes_per_elem
              if (i.NE.l) then
                den = den*(GLL_nodes(l)-GLL_nodes(i))
              endif
              dnum = 1.0_BDKi
              do k = 1,nodes_per_elem
                if ((k.NE.l).AND.(k.NE.i).AND.(i.NE.l)) then
                  dnum = dnum*(QPtN(j)-GLL_nodes(k))
                elseif (i.EQ.l) then
                  dnum = 0.0_BDKi
                endif
              enddo
              ShpDer(l,j) = ShpDer(l,j) + dnum
            enddo
            ShpDer(l,j) = ShpDer(l,j)/den
          endif
        enddo
      enddo
   endif

   Shp(:,:)     = 0.0_BDKi
   do j = 1,nqp
      do l = 1,nodes_per_elem

       if( EqualRealNos(QPtN(j),GLL_nodes(l))) then
         Shp(l,j) = 1.0_BDKi
       else
         dnum = 1.0_BDKi
         den  = 1.0_BDKi
         do k = 1,nodes_per_elem
           if (k.NE.l) then
             den  = den *(GLL_nodes(l) - GLL_nodes(k))
             dnum = dnum*(QPtN(j) - GLL_nodes(k))
           endif
         enddo
         Shp(l,j) = dnum/den
       endif
     enddo
   enddo


 END SUBROUTINE BD_diffmtc


!-----------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!> This subroutine generates n-point gauss-legendre quadrature points and weights
!!
!! Subroutine is based on well-known formulas for generating Legendre polynomials, the
!! roots of which are the Gauss-Legendre quadrature points.  Also used are well-known
!! expressions for the quadrature weights associated with the GL quadrature points.
!! The basic idea of the logic is to use the roots of the Chebyshev polynomial as
!! an initial guess for the roots of the Legendre polynomial, and to then use Newton
!! iteration to find the "exact" roots.
SUBROUTINE BD_GaussPointWeight(n, x, w, QPtWghtDeltaEta, ErrStat, ErrMsg)

   INTEGER(IntKi),INTENT(IN   ):: n       !< Number of Gauss point (p%nqp)
   REAL(BDKi),    INTENT(  OUT):: x(n)    !< Gauss point location (p%QPtN)
   REAL(BDKi),    INTENT(  OUT):: w(n)    !< Gauss point weight (p%QPtWeight)
   REAL(BDKi),    INTENT(  OUT):: QPtWghtDeltaEta(n)    !< Gauss point weight (p%QPtWghtDeltaEta)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat !< Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg  !< Error message if ErrStat /=

   ! local variables

   REAL(BDKi),  PARAMETER:: eps = 1000.0_BDKi*EPSILON(eps) / 2.0_BDKi   ! tolerance for newton-raphson solve (ignores 3 significant digits)


   INTEGER(IntKi) :: n_half

   REAL(BDKi)     :: n_real  ! real version of loop index
   INTEGER(IntKi) :: i       ! loop index
   INTEGER(IntKi) :: j       ! loop index
   REAL(BDKi)     :: j_real  ! real version of loop index

   INTEGER(IntKi) :: newton     ! newton-iteration index
   INTEGER(IntKi), parameter :: newton_max = 10 ! maximum number of newton iterations; should not need many!

   REAL(BDKi)     :: p1  ! legendre polynomial j-2
   REAL(BDKi)     :: p2  ! legendre polynomial j-1
   REAL(BDKi)     :: p3  ! legendre polynomial j

   REAL(BDKi)     :: dp3  ! derivative of legendre polynomial j

   if (n .lt. 2) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'BD_GaussPointWeight: invalid value of n.'
   else
      ErrStat = ErrID_None
      ErrMsg  = ''
   endif

   n_real = REAL(n,BDKi)

   n_half = (n+1) / 2  ! integer for "half" the points: e.g., n=5 => n_half=3, n=4 => n_half = 2

   ! loop through each of the "half" n points
   do i = 1, n_half

      ! Intial guess for ith root; based on ith root of Chebyshev polynomial
      x(i) = - cos( REAL( 2*i - 1, BDKi)  * pi_D / (2.0_BDKi * n_real) )

      ! initialize newton iteration counter
      newton = 0

      p3 = 1. ! some intial big number;  logic will always take at least one newton iteration
      do while (abs(p3) .gt. eps .and. newton .le. newton_max)

         newton = newton + 1  ! keep track of number of newton iterations

         ! build the legendre polynomial and its derivative evaluated at current root guess
         p1 = 1.0_BDKi  ! zeroth order legendre polynomial
         p2 = x(i)      ! first-order legendre polynomial

         do j = 2, n

            j_real = REAL(j, BDKi)

            ! recurrence relation for legendre polynomial
            p3 = ( REAL(2*j - 1, BDKi) * x(i) * p2 - REAL(j - 1, BDKi) * p1 ) / j_real

            ! derivative of legendre polynomial; needed for newton iteration
            dp3 = j_real * ( - x(i) * p3 + p2 ) / ( 1.0_BDKi - x(i)*x(i) )

            ! save values for next iteration
            p1 = p2
            p2 = p3

         enddo

         ! newton iteration: xnew = xold - P_n(xold) / P_n'(xold)
         x(i) = x(i) - p3 / dp3

      enddo

      if (abs(x(i)) .lt. eps) x(i) = 0.0_BDKi

      ! fill in the other "side" of the array
      x(n-i+1) = - x(i)

      ! calculate weight for this Gauss point location
      w(i) = 2.0_BDKi / ( ( 1.0_BDKi - x(i)*x(i) ) * dp3 * dp3 )
      w(n-i+1) = w(i)

   enddo

   ! Find the distance between quadrature points for internal force integration weights
   DO i=1,n-1
      QPtWghtDeltaEta(i)  = (x(i+1)-x(i))/2.0_BDKi      ! right side of QPt
   ENDDO
   QPtWghtDeltaEta(n) = (1.0_BDKi - x(n) ) / 1.0_BDKi

END SUBROUTINE BD_GaussPointWeight
!-----------------------------------------------------------------------------------------------------------------------------------
! This subroutine computes trapezoidal quadrature points and weights, p%QPtN, p%QPtWeight, and p%QPtWghtDeltaEta
!  p%QPtN         -- station location scaled from [-1,1]
!  p%QPtweight    -- weight of QP towards the shape FE nodes
!  p%QPtWghtDeltaEta  -- the delta eta to the right of the node (how much of the
!                    trapezoidal weight is on each side of node when integrating along span).  Used in internal force calculations
SUBROUTINE BD_TrapezoidalPointWeight(p, InputFileData)

   TYPE(BD_ParameterType),INTENT(INOUT):: p              !< BeamDyn parameters
   TYPE(BD_InputFile),    INTENT(IN   ):: InputFileData  !< BeamDyn input-file data

   ! local variables
   REAL(BDKi)                 :: denom ! denominator for quadrature weight computations

   INTEGER(IntKi)             :: indx, temp_id0, temp_id1  ! temporary index into QPtN array
   INTEGER(IntKi), parameter  :: id0 = 1
   INTEGER(IntKi)             :: id1, j

!bjj: this assumes there is only one member

      ! compute the trapezoidal quadrature points, p%QPtN, and scale to range [-1,1]:
      !  If there is refinement, this will add new points between the specified ones. If p%refine == 1, can skip this.
   p%QPtN(1) = InputFileData%InpBl%station_eta(1)
   DO j = 2,p%nqp
      indx =  1+(j-2_IntKi)/p%refine       ! note use of integer math here --> (J-2)/p%refine may not be integer.
      p%QPtN(j) =  InputFileData%InpBl%station_eta(indx) + &
               ((InputFileData%InpBl%station_eta(indx+1) - InputFileData%InpBl%station_eta(indx))/p%refine) * (MOD(j-2,p%refine) + 1)
   ENDDO
   p%QPtN = 2.0_BDKi*p%QPtN - 1.0_BDKi     ! rescale range from [0, 1] to [-1,1]


      ! compute the trapezoidal quadrature weights, p%QPtWeight:
   id1 = InputFileData%InpBl%station_total
   temp_id0 = (id0 - 1)*p%refine + 1            ! Starting index in QPtN --> always going to be 1
   temp_id1 = (id1 - 1)*p%refine + 1            ! ending index in QPtN --> will be  size(p%QPtN)
   denom = p%QPtN(temp_id1) - p%QPtN(temp_id0)  ! This is the range of QPtN --> for single member, is always == 2

   p%QPtWeight(1)    =  (p%QPtN(temp_id0+1) - p%QPtN(temp_id0    ))/denom
   p%QPtWghtDeltaEta(1)  =  (p%QPtN(temp_id0+1) - p%QPtN(temp_id0    ))/denom
   DO j=2,p%nqp-1
      p%QPtWeight(j)  =  (p%QPtN(temp_id0+j) - p%QPtN(temp_id0+j-2))/denom
      p%QPtWghtDeltaEta(j)  = (p%QPtN(j+1)-p%QPtN(j))/denom      ! right side of QPt
   ENDDO
   p%QPtWeight(p%nqp) =  (p%QPtN(temp_id1  ) - p%QPtN(temp_id1-1  ))/denom
   p%QPtWghtDeltaEta(p%nqp) = 0.0_BDKi

END SUBROUTINE BD_TrapezoidalPointWeight


!-----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the weightings required to calculate FE internal forces and moments with the distributed loads from
!! quadrature points.  For both Trap and Gauss quadrature, integration of QP contributions to the next FE inboard uses trapezoidal
!! integration (this works well for both cases).  To expedite the integration since it gets called at each CalcOutput call, the
!! following calculations are done at initialization and stored as parameters. The mathematics in this routine are matched to the
!! algorithm in the BD_InternalForceMomentIGE routine.
SUBROUTINE BD_FEinternalForceQPweights( p, ErrStat, ErrMsg )

   TYPE(BD_ParameterType),       INTENT(INOUT)  :: p              !< Parameters
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat        !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /=

   integer(IntKi)                               :: i              ! generic counter
   integer(IntKi)                               :: idx_FE         ! counter for current FE
   integer(IntKi)                               :: idx_QP         ! counter for current QP
   integer(IntKi)                               :: nelem          ! counter for current element

   real(BDKi)                                   :: weightRange    ! Weighted range of current integration region
   real(BDKi)                                   :: weightQPplus   ! Weighting of QP region outboard
   real(BDKi)                                   :: weightQPminus  ! Weighting of QP region inboard


   integer(intKi)                               :: ErrStat2       ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2        ! temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_FEinternalForceQPweights'

   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Description of arrays
      !  The arrays calculated by this routine are
      !     QPrangeOverlapFE        ( 2, p%nodes_per_elem, p%elem_total )
      !                 This array stores the index to the first and last QP that contribute to a given
      !                 FE point.  This is used as the bounds in an integration loop.
      !                 This is only used for calculating outputs at FE nodes.
      !     QPtWghtIntForceFE       ( p%nqp, 2, p%nodes_per_elem, p%elem_total )
      !                 This array stores the contribution weighting of a given QP.  It is split into
      !                 the inboard contribution of the QP (such as QP-1 to QP) which is stored in
      !                 (:,1,:,:), and the outboard contribution of the QP (such as QP to QP+1) which
      !                 is stored in (:,2,:,:).  This is split this way because the physical range of
      !                 [QP-1:QP] and [QP:QP+1] will be different in the mapping of the force to moment
      !                 at the next inboard location.
      !                 This is only used for calculating outputs at FE nodes.
      !     FEoutboardOfQPt         ( p%nodes_per_elem, p%nqp )
      !                 This array is used in calculating the mapping of point loads at FE points onto
      !                 the QPs. It simpifies some of the logic.
      !                 This is used in calculating outputs at quadrature points.
      !  These arrays are used to speed up the routine that calculates the internal forces.  The data within these
      !  arrays is static, so it is calculated at initialization and stored in the parameters.


   !----------------------------------------------------------------
   ! For calculations of internal forces at FE nodes.

      ! NOTE:  if we are outputting the intenral forces at FE nodes with gaussian quadrature, the integration
      !         is known be inexact.  So we will issue a warning.
   if ((p%quadrature .EQ. GAUSS_QUADRATURE) .and. (p%BldMotionNodeLoc == BD_MESH_FE) .and. p%OutIntForce ) then
      call SetErrStat( ErrID_Warn, 'WARNING: Internal forces are calculated using trapezoidal rule.  This will yield potentially inacurate ' // &
               'results when used with gaussian quadrature. You have been warned.',ErrStat,ErrMsg,RoutineName)
   endif

      ! Initialize arrays to zero
   p%QPrangeOverlapFE   = 0_IntKi
   p%QPtWghtIntForceFE  = 0.0_BDKi

      ! Working from the tip inwards (to match the calculations in the BD_InternalForceMomentIGE routine)
      ! Calculate the limits and weights for each QP to each FE.  Note that many do not contribute.
   do nelem=p%elem_total,1,-1

         !  Last QP.
         !
         !           | <- Region -> |
         !  ---------x--------------#
         !        last QP        last FE
         !
         !  then we need to count the right side of it (if any).  The only time this is
         !  expected to occur is in gaussian quadrature where the first and last QP are
         !  not at the ends. Compare in eta.
      if ( p%QPtN(p%nqp)      < p%GLL_Nodes(p%nodes_per_elem) ) then

            ! weightings -- full range from last QP to last FE is used
         weightRange = 1.0_BDKi
         weightQPplus = p%GLL_Nodes(p%nodes_per_elem) - p%QPtN(p%nqp)

            ! The last QP (p%nqp) has a region to right that contributes the next FE inboard.
            ! NOTE:  This is handled as a special case in the integration, so don't set the qp limits
            !        (it leads to out of array bounds in other calcs when a QP+1 is beyond the end
            !        of the array.  This is due to the algorithm used).
         p%QPtWghtIntForceFE(p%nqp,2,p%nodes_per_elem,nelem) = weightRange * weightQPplus    ! only from the QP to the outboard has this weighting

      endif


         ! For all the remaining QPs, check if it is a range that contributes to the FE.  There
         ! are 4 distinct cases handled by this set of limits.
         !
         !  case 1: region left of QP+1 contributes to FE
         !              | <- Region -> |
         !  ---x--------#--------------x---------#---- ...
         !     QP       FE            QP+1      FE+1
         !
         !  case 2: full region between QP and QP+1 contribute to FE
         !              | <- Region -> |
         !  ---#--------x--------------x---------#---- ...
         !     FE       QP            QP+1      FE+1
         !
         !  case 3: region right of QP to FE+1
         !              | <- Region -> |
         !  ---#--------x--------------#---------x---- ...
         !     FE       QP            FE+1      QP+1
         !
         !  case 4: region between FE and FE+1 contributes to FE
         !              | <- Region -> |
         !  ---x--------#--------------#---------x---- ...
         !     QP       FE            FE+1      QP+1
         !
      do idx_FE=p%nodes_per_elem-1,1,-1
         do idx_qp=p%nqp-1,1,-1
               ! A region between this QP and the next contribute to this FE
            if (     ((p%QPtN(idx_qp)      <  p%GLL_Nodes(idx_FE)) .and. (p%GLL_Nodes(idx_FE)  < p%QPtN(idx_qp+1)     ))   &      ! FE between pair of QP points
               .or.  ((p%GLL_Nodes(idx_FE) <= p%QPtN(idx_qp)     ) .and. (p%QPtN(idx_qp)       < p%GLL_Nodes(idx_FE+1))) ) then   ! QP between pair of FE points

                  ! The following weights are derived by trapezoidal integration over the region of interest
                  ! (weighted average of QP values over region)
                  ! Region of interest divided by range between QP's contributing
               weightRange    = (min(p%GLL_Nodes(idx_FE+1),p%QPtN(idx_qp+1)) - max(p%GLL_Nodes(idx_FE),p%QPtN(idx_qp))) &
                                 / (p%QPtN(idx_qp+1) - p%QPtN(idx_qp))
                  ! contribution of QP+1
               weightQPplus   = (p%QPtN(idx_qp+1) - max(p%GLL_Nodes(idx_FE),p%QPtN(idx_qp)))       &  ! Always have this term
                              + (p%QPtN(idx_qp+1) - min(p%GLL_Nodes(idx_FE+1),p%QPtN(idx_qp+1)))      ! this term is zero if qp+1 < FE+1
                  ! contribution of QP-1
               weightQPminus  = (min(p%GLL_Nodes(idx_FE+1),p%QPtN(idx_qp+1)) - p%QPtN(idx_qp))     &  ! Always have this term
                              + (p%GLL_Nodes(idx_FE) - min(p%GLL_Nodes(idx_FE),p%QPtN(idx_qp)))       ! this term is zero if qp > FE

                  ! Store values
               p%QPrangeOverlapFE(1,idx_FE,nelem) = idx_qp                                            ! min of QP range that contributes
               p%QPrangeOverlapFE(2,idx_FE,nelem) = max(idx_qp,p%QPrangeOverlapFE(2,idx_FE,nelem))    ! max of QP range that contributes
               p%QPtWghtIntForceFE(idx_qp+1,1,idx_FE,nelem) = 0.5_BDKi * weightRange * weightQPminus  ! inboard  of next QP
               p%QPtWghtIntForceFE(idx_qp  ,2,idx_FE,nelem) = 0.5_BDKi * weightRange * weightQPplus   ! outboard of this QP

            endif ! in region contributing
         enddo    ! loop over QP
      enddo       ! loop over FE

         ! First QP:
         !
         !  | <- Region -> |
         !  #--------------x--------
         !  first FE       first QP
         !
         !  then we need to count the left side of it (if any).  The only time this is
         !  expected to occur is in gaussian quadrature where the first and last QP are
         !  not at the ends. Compare in eta.
            ! Need to add the inboard side of idx_qp=1 to idx_FE=1 (expect this only occurs in gaussian)
      if ( p%QPtN(1) >= p%GLL_Nodes(1) ) then

            ! weightings -- full range from first FE to first QP is used
         weightRange = 1.0_BDKi
         weightQPminus = p%QPtN(1) - p%GLL_Nodes(1)

            ! The first QP (1) has a region to left that contributes the first FE
         p%QPrangeOverlapFE(1,1,nelem) = 1
         p%QPtWghtIntForceFE(1,1,1,nelem) = weightRange * weightQPminus
      endif
   enddo       ! loop over elements


   !----------------------------------------------------------------
   ! For calculating outputs at Quadrature points

   !  For each FE, find the next QP inboard.
   !  Points loads frome the driver code are applied at the FE points.  This does not occur when coupled to FAST.
   !  These values are used in the integration of the internal forces at the QP's (trap quadrature).  For the integration
   !  we only consider FE nodes outboard of the QP before the next QP.  To simplify the integration loops, we
   !  set values in the p%FEoutboardOfQPt array to indicate logical status (we may have a distance of 0 and need to apply
   !  the load, but that must be calculated later in 3D).
   p%FEoutboardOfQPt = .FALSE.
   DO idx_FE=1,p%nodes_per_elem
      DO idx_qp=1,p%nqp-1
         if (( p%QPtN(idx_qp) <= p%GLL_Nodes(idx_FE) ) .and. ( p%QPtN(idx_qp+1) > p%GLL_Nodes(idx_FE) )) then
            p%FEoutboardOfQPt(idx_FE,idx_qp) = .TRUE.
         endif
      END DO
      ! last QP -- in the event of a gaussian quadrature, there is no QP at the last FE.
      if ( p%QPtN(p%nqp) <= p%GLL_Nodes(idx_FE) ) then
         p%FEoutboardOfQPt(idx_FE,p%nqp) = .TRUE.
      endif
   END DO
END SUBROUTINE BD_FEinternalForceQPweights

!-----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates y%BldMotion%TranslationDisp, y%BldMotion%Orientation, y%BldMotion%TranslationVel, and
!! y%BldMotion%RotationVel, which depend only on states (and indirectly, u%RootMotion), and parameters.
SUBROUTINE Set_BldMotion_NoAcc(p, x, m, y)

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(BD_MiscVarType),         INTENT(IN   )  :: m           !< misc/optimization variables
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   INTEGER(IntKi)                               :: i
   INTEGER(IntKi)                               :: j
   INTEGER(IntKi)                               :: temp_id
   INTEGER(IntKi)                               :: temp_id2
   REAL(BDKi)                                   :: cc(3)
   REAL(BDKi)                                   :: cc0(3)
   REAL(BDKi)                                   :: temp_R(3,3)
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Set_BldMotion_NoAcc'

   ! The first node on the mesh is just the root location, but since the orientation of the root node may not
   ! RootMotion orientation (due to blade twist), we will calculated it as though it was a regular node

   ! now fill in the other nodes
   SELECT CASE (p%BldMotionNodeLoc)

   CASE (BD_MESH_FE)

      DO i=1,p%elem_total
         DO j=1,p%nodes_per_elem
            temp_id = (i-1)*(p%nodes_per_elem-1)+j      ! The last node of the first element is used as the first node in the second element.
            temp_id2= (i-1)*p%nodes_per_elem+j          ! Index to a node within element i


               ! Calculate the translational displacement of each GLL node in the FAST coordinate system,
               ! referenced against the DCM of the blade root at T=0.
            y%BldMotion%TranslationDisp(1:3,temp_id2) = MATMUL(p%GlbRot,x%q(1:3,temp_id))

!bjj: note differences here compared to BDrot_to_FASTdcm
!adp: in BDrot_to_FASTdcm we are assuming that x%q(4:6,:) is zero because there is no rotatinoal displacement yet
               ! Find the rotation parameter in global coordinates (initial orientation + rotation parameters)
               ! referenced against the DCM of the blade root at T=0.
            CALL BD_CrvCompose( cc, x%q(4:6,temp_id), p%uuN0(4:6,j,i), FLAG_R1R2 )
            CALL BD_CrvCompose( cc0, p%Glb_crv, cc, FLAG_R1R2 )

               ! Create the DCM from the rotation parameters
            CALL BD_CrvMatrixR(cc0,temp_R)  ! returns temp_R (the transpose of the DCM orientation matrix)

               ! Store the DCM for the j'th node of the i'th element (in FAST coordinate system)
            y%BldMotion%Orientation(1:3,1:3,temp_id2) = TRANSPOSE(temp_R)

               ! Calculate the translation velocity and store in FAST coordinate system
               ! referenced against the DCM of the blade root at T=0.
            y%BldMotion%TranslationVel(1:3,temp_id2) = MATMUL(p%GlbRot,x%dqdt(1:3,temp_id))

               ! Calculate the rotational velocity and store in FAST coordinate system
               ! referenced against the DCM of the blade root at T=0.
            y%BldMotion%RotationVel(1:3,temp_id2) = MATMUL(p%GlbRot,x%dqdt(4:6,temp_id))

         ENDDO
      ENDDO

   CASE (BD_MESH_QP)

      DO i=1,p%elem_total
         DO j=1,p%nqp
            temp_id2 = (i-1)*p%nqp + j + p%qp_indx_offset            ! Index to a node within element i

               ! Calculate the translational displacement of each quadrature node in the FAST coordinate system,
               ! referenced against the DCM of the blade root at T=0.
            y%BldMotion%TranslationDisp(1:3,temp_id2) = MATMUL(p%GlbRot,m%qp%uuu(1:3,j,i) )


!bjj: note differences here compared to BDrot_to_FASTdcm
!adp: in BDrot_to_FASTdcm we are assuming that x%q(4:6,:) is zero because there is no rotatinoal displacement yet
               ! Find the rotation parameter in global coordinates (initial orientation + rotation parameters)
               ! referenced against the DCM of the blade root at T=0.
            CALL BD_CrvCompose( cc, m%qp%uuu(4:6,j,i), p%uu0(4:6,j,i), FLAG_R1R2 )
            CALL BD_CrvCompose( cc0, p%Glb_crv, cc, FLAG_R1R2 )

            CALL BD_CrvMatrixR(cc0,temp_R)  ! returns temp_R (the transpose of the DCM orientation matrix)
               ! Store the DCM for the j'th node of the i'th element (in FAST coordinate system)
            y%BldMotion%Orientation(1:3,1:3,temp_id2) = TRANSPOSE(temp_R)

               ! Calculate the translation velocity and store in FAST coordinate system
               ! referenced against the DCM of the blade root at T=0.
            y%BldMotion%TranslationVel(1:3,temp_id2) = MATMUL(p%GlbRot,m%qp%vvv(1:3,j,i))

               ! Calculate the rotational velocity and store in FAST coordinate system
               ! referenced against the DCM of the blade root at T=0.
            y%BldMotion%RotationVel(1:3,temp_id2) = MATMUL(p%GlbRot,m%qp%vvv(4:6,j,i))

         ENDDO
      ENDDO



   CASE (BD_MESH_STATIONS)
   END SELECT


END SUBROUTINE Set_BldMotion_NoAcc
!-----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates values for the y%BldMotion mesh.
SUBROUTINE Set_BldMotion_Mesh(p, u, x, m, y)

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_InputType),           INTENT(IN   )  :: u           !< Inputs at t - in the FAST coordinate system (NOT BD)
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables ! intent(out) so that we can update the accelerations here...
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   INTEGER(IntKi)                               :: i
   INTEGER(IntKi)                               :: j
   INTEGER(IntKi)                               :: j_start !starting node on this element
   INTEGER(IntKi)                               :: idx_node
   INTEGER(IntKi)                               :: temp_id
   INTEGER(IntKi)                               :: temp_id2
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Set_BldMotion_Mesh'


      ! set positions and velocities (not accelerations)
   call Set_BldMotion_NoAcc(p, x, m, y)

   ! Only need this bit for dynamic cases
   IF ( p%analysis_type /= BD_STATIC_ANALYSIS ) THEN

       ! now set the accelerations:

       ! The first node on the mesh is just the root location:
       y%BldMotion%TranslationAcc(:,1)    = u%RootMotion%TranslationAcc(:,1)
       y%BldMotion%RotationAcc(:,1)       = u%RootMotion%RotationAcc(:,1)

       SELECT CASE (p%BldMotionNodeLoc)

       CASE (BD_MESH_FE) ! This is how the NREL-version of BeamDyn works (output nodes at the finite element nodes)

          j_start = 2 ! we'll skip the first node on the first element; otherwise we will start at the first node on the other elements
          DO i=1,p%elem_total
             DO j=j_start,p%nodes_per_elem
                temp_id = (i-1)*(p%nodes_per_elem-1)+j
                temp_id2= (i-1)*p%nodes_per_elem+j

                y%BldMotion%TranslationAcc(1:3,temp_id2) = MATMUL(p%GlbRot, m%RHS(1:3,temp_id) )

                y%BldMotion%RotationAcc(1:3,temp_id2) = MATMUL(p%GlbRot, m%RHS(4:6,temp_id) )
             ENDDO
             j_start = 1
          ENDDO

       CASE (BD_MESH_QP)


          m%qp%aaa(1:3,1,1) = u%RootMotion%TranslationAcc(:,1)
          m%qp%aaa(4:6,1,1) = u%RootMotion%RotationAcc(:,1)

            ! Calculate the and acceleration term at t+dt (OtherState%acc is at t+dt)
            j_start = 2 ! we'll skip the first node on the first element; otherwise we will start at the first node on the other elements
            DO i=1,p%elem_total
                DO j=j_start,p%nqp

                ! Initialize to zero for summation, then recalculate accelerations at quadrature nodes based on accelerations at FE nodes
                m%qp%aaa(:,j,i) = 0.0_BDKi
                DO idx_node=j_start,p%nodes_per_elem
                    m%qp%aaa(:,j,i) = m%qp%aaa(:,j,i) + p%Shp(idx_node,j) * m%RHS(:,p%node_elem_idx(i,1)-1+idx_node)
                ENDDO

                temp_id2 = (i-1)*p%nqp + j + p%qp_indx_offset            ! Index to a node within element i

                    ! Calculate the translational acceleration of each quadrature node in the FAST coordinate system,
                    ! referenced against the DCM of the blade root at T=0.
                y%BldMotion%TranslationAcc(1:3,temp_id2) = MATMUL(p%GlbRot,m%qp%aaa(1:3,j,i) )

                y%BldMotion%RotationAcc(1:3,temp_id2) = MATMUL(p%GlbRot, m%qp%aaa(4:6,j,i) )
                ENDDO
                j_start = 1
            ENDDO

       CASE (BD_MESH_STATIONS)
       END SELECT

   END IF

END SUBROUTINE Set_BldMotion_Mesh
!> This routine calculates values for the y%BldMotion mesh.
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Set_BldMotion_InitAcc(p, u, OtherState, m, y)

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_InputType),           INTENT(IN   )  :: u           !< Inputs at t - in the FAST coordinate system (NOT BD)
   TYPE(BD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables ! intent(out) so that we can update the accelerations here...
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   INTEGER(IntKi)                               :: i
   INTEGER(IntKi)                               :: j
   INTEGER(IntKi)                               :: j_start !starting node on this element
   INTEGER(IntKi)                               :: temp_id
   INTEGER(IntKi)                               :: temp_id2
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Set_BldMotion_InitAcc'


      ! The first node on the mesh is just the root location:
      y%BldMotion%TranslationAcc(:,1)    = u%RootMotion%TranslationAcc(:,1)
      y%BldMotion%RotationAcc(:,1)       = u%RootMotion%RotationAcc(:,1)

      SELECT CASE (p%BldMotionNodeLoc)

      CASE (BD_MESH_FE) ! This is how the NREL-version of BeamDyn works (output nodes at the finite element nodes)

         j_start = 2 ! we'll skip the first node on the first element; otherwise we will start at the first node on the other elements
         DO i=1,p%elem_total
            DO j=j_start,p%nodes_per_elem
               temp_id = (i-1)*(p%nodes_per_elem-1)+j
               temp_id2= (i-1)*p%nodes_per_elem+j

               y%BldMotion%TranslationAcc(1:3,temp_id2) = MATMUL(p%GlbRot, OtherState%Acc(1:3,temp_id) )

               y%BldMotion%RotationAcc(1:3,temp_id2) = MATMUL(p%GlbRot, OtherState%Acc(4:6,temp_id) )
            ENDDO
            j_start = 1
         ENDDO

      CASE (BD_MESH_QP)

         ! Calculate the and acceleration term at t+dt (OtherState%acc is at t+dt)
         j_start = 2 ! we'll skip the first node on the first element; otherwise we will start at the first node on the other elements
         DO i=1,p%elem_total
            DO j=j_start,p%nqp

               temp_id2 = (i-1)*p%nqp + j + p%qp_indx_offset            ! Index to a node within element i

                  ! Calculate the translational acceleration of each quadrature node in the FAST coordinate system,
                  ! referenced against the DCM of the blade root at T=0.
               y%BldMotion%TranslationAcc(1:3,temp_id2) = MATMUL(p%GlbRot,m%qp%aaa(1:3,j,i) )

               y%BldMotion%RotationAcc(1:3,temp_id2) = MATMUL(p%GlbRot, m%qp%aaa(4:6,j,i) )
            ENDDO
            j_start = 1
         ENDDO

      CASE (BD_MESH_STATIONS)
      END SELECT


END SUBROUTINE Set_BldMotion_InitAcc
!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine finds the (initial) nodal position and curvature based on a relative position, eta, along the blade.
!!    If the elementNumber is passed, we consider eta to be on the element.
!!    If the elementNumber is NOT passed, we consider eta to be on the whole blade, so then we must find the appropriate
!!       element before calling the BD_ComputeIniNodalPosition routine. --> This is used by the driver for point loads.
!! The key points are used to obtain the z coordinate of the physical distance along the blade, and POS and CRV vectors are returned
!! from that location.
subroutine Find_IniNode(kp_coordinate, p, member_first_kp, member_last_kp, etaPassed, POS, CRV, ErrStat, ErrMsg, elementNumber)

   REAL(BDKi),                   intent(in   )  :: kp_coordinate(:,:)  !< Key Point coordinates
   type(BD_ParameterType),       intent(in   )  :: p                   !< Parameters
   INTEGER(IntKi),               intent(in   )  :: member_first_kp     !< index of the first key point on a particular member
   INTEGER(IntKi),               intent(in   )  :: member_last_kp      !< index of the last key point on a particular member
   REAL(BDKi),                   intent(in   )  :: etaPassed           !< relative position of desired node, [0,1] on this member (if element number is provided) 
   REAL(BDKi),                   intent(  out)  :: POS(3)              !< position of node (in BD coordinates)
   REAL(BDKi),                   intent(  out)  :: CRV(3)              !< curvature of node (in BD coordinates)
   integer(IntKi),               intent(  out)  :: ErrStat             !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   integer(IntKi), optional,     intent(in   )  :: elementNumber       !< current element

   REAL(BDKi),PARAMETER    :: EPS = 1.0D-10


   ! local variables
   INTEGER(IntKi)          :: kp                ! key point
   INTEGER(IntKi)          :: TmpIdx            ! Temporary index returned by interp routine: we completely ignore this
   INTEGER(IntKi)          :: MaxIdx            ! Maximum index we can try to get info from
   INTEGER(IntKi)          :: elem              ! current element
   REAL(BDKi)              :: Zpos              ! distance (in z coordinate) associated with eta along this member, in meters
   REAL(BDKi)              :: eta               ! relative position of desired node, [0,1] on this member
   REAL(BDKi)              :: temp_twist
   REAL(BDKi)              :: temp_e1(3)

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Find_IniNode'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   TmpIdx  = 1

   if (present(elementNumber)) then
      elem = elementNumber
      eta  = etaPassed
   else
      ! Find which element we are on (assuming eta as the blade eta)
      elem=1
      do while (p%ztoEtaMapping(1,2,elem) > eta)
         elem = elem + 1
      enddo
      ! Now find what the element eta is
      MaxIdx=maxloc(p%zToEtaMapping(:,2,elem),1)
      eta = InterpStp( etaPassed, p%zToEtaMapping(1:MaxIdx,2,elem), p%zToEtaMapping(1:MaxIdx,3,elem), TmpIdx, MaxIdx)
   endif

   ! Find the z value along the beam span from eta
   !     zToEtaMapping -- idx 1: z  value                (there may be zeros at the end)
   !     zToEtaMapping -- idx 2: eta value for blade     (there may be zeros at the end)
   !     zToEtaMapping -- idx 3: eta value for element   (there may be zeros at the end)
   !  InterpStp( XvalOfInterst, Xary, Yary, startSearchIndex, MaxIndex )
   !  NOTE: with zeros at end of arrays, we need to know where to terminate the search
   MaxIdx=maxloc(p%zToEtaMapping(:,1,elem),1)
   Zpos = InterpStp( eta, p%zToEtaMapping(1:MaxIdx,3,elem), p%zToEtaMapping(1:MaxIdx,1,elem), TmpIdx, MaxIdx)


   ! find the first key point that is beyond where this node is on the blade
   ! note that this is the index for p%SP_Coef, so the upper bound is member_last_kp-1 instead of member_last_kp
   ! NOTE: since we are given the starting kp, we know which element we are on here (eta is on the element, not the whole blade)
   ! bjj: to be more efficient, we could probably just start at the kp we found for the previous eta
   kp = member_first_kp
   DO WHILE ( (eta > p%segment_eta(kp) + EPS) .and. kp < (member_last_kp-1) )
      kp = kp + 1
   END DO

   ! using the spline coefficients at this key point, compute the position and orientation of the node
   CALL BD_ComputeIniNodalPosition(p%SP_Coef(kp,:,:),Zpos,POS,temp_e1,temp_twist) ! Compute point physical coordinates (POS) in blade frame
   CALL BD_ComputeIniNodalCrv(temp_e1, temp_twist, CRV, ErrStat2, ErrMsg2)        ! Compute initial rotation parameters (CRV) in blade frame
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return


end subroutine Find_IniNode
!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the initial nodal locations given the coefficients for
!! cubic spline fit. It also computes the unit tangent vector e1 for further use.
SUBROUTINE BD_ComputeIniNodalPosition(SP_Coef,Zpos,PosiVec,e1,Twist_Angle)

   REAL(BDKi),    INTENT(IN   ):: SP_Coef(:,:)  !< Coefficients for cubic spline interpolation
   REAL(BDKi),    INTENT(IN   ):: Zpos           !< z-component of nodal location (that's how SP_Coef was computed), in meters
   REAL(BDKi),    INTENT(  OUT):: PosiVec(:)    !< Physical coordinates of points in blade frame
   REAL(BDKi),    INTENT(  OUT):: e1(:)         !< Tangent vector, normalized
   REAL(BDKi),    INTENT(  OUT):: Twist_Angle   !< Twist angle at PosiVec

   INTEGER(IntKi)              :: i

   DO i=1,3
       PosiVec(i) = SP_Coef(1,i) + SP_Coef(2,i)*Zpos +          SP_Coef(3,i)*Zpos**2 +          SP_Coef(4,i)*Zpos**3 !position
       e1(i)      =                SP_Coef(2,i)      + 2.0_BDKi*SP_Coef(3,i)*Zpos    + 3.0_BDKi*SP_Coef(4,i)*Zpos**2 !tangent (derivative w.r.t. Zpos)
   ENDDO
   e1 = e1/TwoNorm(e1) ! normalize tangent vector

   Twist_Angle = SP_Coef(1,4) + SP_Coef(2,4)*Zpos + SP_Coef(3,4)*Zpos**2 + SP_Coef(4,4)*Zpos**3

END SUBROUTINE BD_ComputeIniNodalPosition
!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes initial CRV parameters
!! given geometry information
SUBROUTINE BD_ComputeIniNodalCrv(e1, phi, cc, ErrStat, ErrMsg)

   REAL(BDKi),     INTENT(IN   )  :: e1(:)         !< Tangent unit vector
   REAL(BDKi),     INTENT(IN   )  :: phi           !< Initial twist angle, in degrees
   REAL(BDKi),     INTENT(  OUT)  :: cc(:)         !< Initial Crv Parameter
   INTEGER(IntKi), INTENT(  OUT)  :: ErrStat       !< Error status of the operation
   CHARACTER(*),   INTENT(  OUT)  :: ErrMsg        !< Error message if ErrStat /= ErrID_None

   REAL(BDKi)                     :: e2(3)         !< Unit normal vector
   REAL(BDKi)                     :: Rr(3,3)       !< Initial rotation matrix
   REAL(BDKi)                     :: PhiRad        !< Phi in radians
   REAL(BDKi)                     :: Delta

   INTEGER(IntKi)                 :: ErrStat2      ! Temporary Error status
   CHARACTER(ErrMsgLen)           :: ErrMsg2       ! Temporary Error message
   CHARACTER(*), PARAMETER        :: RoutineName = 'BD_ComputeIniNodalCrv'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   PhiRad = phi*D2R_D  ! convert to radians

      ! Note that e1 corresponds with the Z-axis direction in our formulation.  For Beam theory, this would be the x-axis.
   Rr(:,3)  =  e1(:)

   e2(3)    = -(e1(1)*COS(PhiRad) + e1(2)*SIN(PhiRad))/e1(3)
   Delta    = SQRT(1.0_BDKi + e2(3)*e2(3))
   e2(1)    = COS(PhiRad)
   e2(2)    = SIN(PhiRad)
   e2       = e2 / Delta
   Rr(:,1)  = e2
   Rr(:,2)  = Cross_Product(e1,e2)

   CALL BD_CrvExtractCrv(Rr, cc, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

END SUBROUTINE BD_ComputeIniNodalCrv
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ExtractRelativeRotation(R, p, rr, ErrStat, ErrMsg)
   real(R8Ki),             INTENT(in   )     :: R(3,3)       !< input rotation matrix (transpose of DCM; in BD coords)
   type(BD_ParameterType), INTENT(in   )     :: p            !< Parameters
   real(BDKi),             INTENT(  OUT)     :: rr(3)        !< W-M parameters of relative rotation
   INTEGER(IntKi),         INTENT(  OUT)     :: ErrStat      !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT)     :: ErrMsg       !< Error message if ErrStat /= ErrID_None

   real(BDKi)                                :: R_WM(3)      ! W-M parameters of R
   real(BDKi)                                :: R_BD(3,3)    ! input rotation matrix in BDKi precision

   INTEGER(IntKi)                            :: ErrStat2     ! Temporary Error status
   CHARACTER(ErrMsgLen)                      :: ErrMsg2      ! Temporary Error message
   CHARACTER(*), PARAMETER                   :: RoutineName = 'ExtractRelativeRotation'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Calculate p(rr) = p(Glb_crv)^- (+) p(R)       where (+) is the curve compose
   !  which is the same as operation as
   !     R(rr) = R(Glb_crv)^T R

   ! note that the u%RootMotion mesh does not contain the initial twist, but p%Glb_crv does not have this twist, either.
   ! The relative rotation will be the same in this case.

   R_BD = R ! possible type conversion (only if BDKi /= R8Ki)

   CALL BD_CrvExtractCrv(R_BD,R_WM, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   CALL BD_CrvCompose(rr,p%Glb_crv,R_WM,FLAG_R1TR2)         ! rr = p%Glb_crv^- composed with R_WM

   ! NOTE: the above calculation is not the inverse of what is in Set_BldMotion_NoAcc.  The reason is that this
   !       routine is only looking at RootMotion.  The Set_BldMotion_NoAcc routine is looking at the blade motion
   !       which at the root differs by the WM values in p%uuN0(4:6,1,1) from the RootMotion mesh.

END SUBROUTINE ExtractRelativeRotation
!-----------------------------------------------------------------------------------------------------------------------------------
FUNCTION BDrot_to_FASTdcm(rr,p) RESULT(dcm)
   real(BDKi),             intent(in) :: rr(3)        !< W-M parameters of relative rotation
   type(BD_ParameterType), intent(in) :: p            !< Parameters
   real(BDKi)                         :: dcm(3,3)     !< input rotation matrix (transpose of DCM; in BD coords)


   REAL(BDKi)                         :: temp_CRV2(3)   ! temp curvature parameters
   real(BDKi)                         :: R(3,3)         ! rotation matrix

! note differences in setting up meshes with Set_BldMotion_NoAcc
!adp: in the case of the meshes in Set_BldMotion_NoAcc, x%q(4:6,:) and m%qp%uuu(4:6,:,:) are not zero.  When this routine is called, they
!     are zero, and the expression in Set_BldMotion_NoAcc simplifies to this expression.

      ! rotate relative W-M rotations to global system?
   CALL BD_CrvCompose(temp_CRV2,p%Glb_crv,rr,FLAG_R1R2) !temp_CRV2 = p%Glb_crv composed with rr

      ! create rotation matrix from W-M parameters:
   CALL BD_CrvMatrixR(temp_CRV2,R) ! returns R (rotation matrix, the transpose of the DCM orientation matrix)

      ! get DCM from rotation matrix:
   dcm = TRANSPOSE(R)


END FUNCTION BDrot_to_FASTdcm
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE BeamDyn_Subs
