MODULE FVW
   USE NWTC_Library
   USE FVW_Types
   USE FVW_Subs
   USE AD14AeroConf_Types


   IMPLICIT NONE

   TYPE(ProgDesc), PARAMETER  :: FVW_Ver = ProgDesc( 'FVW', '', '' )

   PUBLIC   :: FVW_Init             ! Initialization routine
   PUBLIC   :: FVW_End

   PUBLIC   :: FVW_CalcOutput
   PUBLIC   :: FVW_UpdateStates



CONTAINS


!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine FVW_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   type(FVW_InitInputType),         intent(in   )  :: InitInp        !< Input data for initialization routine
   type(FVW_InputType),             intent(  out)  :: u              !< An initial guess for the input; input mesh must be defined
   type(FVW_ParameterType),         intent(  out)  :: p              !< Parameters
   type(FVW_ContinuousStateType),   intent(  out)  :: x              !< Initial continuous states
   type(FVW_DiscreteStateType),     intent(  out)  :: xd             !< Initial discrete states
   type(FVW_ConstraintStateType),   intent(  out)  :: z              !< Initial guess of the constraint states
   type(FVW_OtherStateType),        intent(  out)  :: OtherState     !< Initial other states
   type(FVW_OutputType),            intent(  out)  :: y              !< Initial system outputs (outputs are not calculated;
                                                                     !!   only the output mesh is initialized)
   type(FVW_MiscVarType),           intent(  out)  :: m              !< Initial misc/optimization variables
   real(DbKi),                      intent(inout)  :: interval       !< Coupling interval in seconds: the rate that
                                                                     !!   (1) FVW_UpdateStates() is called in loose coupling &
                                                                     !!   (2) FVW_UpdateDiscState() is called in tight coupling.
                                                                     !!   Input is the suggested time from the glue code;
                                                                     !!   Output is the actual coupling interval that will be used
                                                                     !!   by the glue code.
   type(FVW_InitOutputType),        intent(  out)  :: InitOut        !< Output for initialization routine
   integer(IntKi),                  intent(  out)  :: errStat        !< Error status of the operation
   character(*),                    intent(  out)  :: errMsg         !< Error message if ErrStat /= ErrID_None
   

      ! Local variables
   
   integer(IntKi)                                  :: errStat2       ! temporary error status of the operation
   character(ErrMsgLen)                            :: errMsg2        ! temporary error message 
      
   integer(IntKi)                                  :: UnEcho         ! Unit number for the echo file
   
   character(*), parameter                         :: RoutineName = 'FVW_Init'
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   UnEcho  = -1

      ! Initialize the NWTC Subroutine Library

   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   call DispNVD( FVW_Ver )



      ! Set transformation matrix to go from FVW to AD coordinates
   p%TransformFVWtoAD(1,:)=(/ 0.00_ReKi,  0.00_ReKi,  1.00_ReKi /)
   p%TransformFVWtoAD(2,:)=(/ 0.00_ReKi, -1.00_ReKi,  0.00_ReKi /)
   p%TransformFVWtoAD(3,:)=(/ 1.00_ReKi,  0.00_ReKi,  0.00_ReKi /)

      ! Set transformation matrix to go from FVW to AD coordinates
      ! Note: this is the inverse of (and identical to) the TransformFVWtoAD matrix
   p%TransformADtoFVW(1,:)=(/ 0.00_ReKi,  0.00_ReKi,  1.00_ReKi /)
   p%TransformADtoFVW(2,:)=(/ 0.00_ReKi, -1.00_ReKi,  0.00_ReKi /)
   p%TransformADtoFVW(3,:)=(/ 1.00_ReKi,  0.00_ReKi,  0.00_ReKi /)

end subroutine FVW_Init


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine FVW_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   type(FVW_InputType),             intent(inout)  :: u           !< System inputs
   type(FVW_ParameterType),         intent(inout)  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(inout)  :: x           !< Continuous states
   type(FVW_DiscreteStateType),     intent(inout)  :: xd          !< Discrete states
   type(FVW_ConstraintStateType),   intent(inout)  :: z           !< Constraint states
   type(FVW_OtherStateType),        intent(inout)  :: OtherState  !< Other states
   type(FVW_OutputType),            intent(inout)  :: y           !< System outputs
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:



         ! Destroy the input data:

      call FVW_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      call FVW_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      call FVW_DestroyContState(   x,           ErrStat, ErrMsg )
      call FVW_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      call FVW_DestroyConstrState( z,           ErrStat, ErrMsg )
      call FVW_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
      call FVW_DestroyMisc(        m,           ErrStat, ErrMsg ) 

         ! Destroy the output data:

      call FVW_DestroyOutput( y, ErrStat, ErrMsg )

end subroutine FVW_End



!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine FVW_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
!..................................................................................................................................

   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   integer(IntKi),                  intent(in   )  :: n           !< Current simulation time step n = 0,1,...
   type(FVW_InputType),             intent(inout)  :: u(:)        !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                      intent(in   )  :: utimes(:)   !< Times associated with u(:), in seconds
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(inout)  :: x           !< Input: Continuous states at t;
                                                                  !!   Output: Continuous states at t + Interval
   type(FVW_DiscreteStateType),     intent(inout)  :: xd          !< Input: Discrete states at t;
                                                                  !!   Output: Discrete states at t  + Interval
   type(FVW_ConstraintStateType),   intent(inout)  :: z           !< Input: Constraint states at t;
                                                                  !!   Output: Constraint states at t+dt
   type(FVW_OtherStateType),        intent(inout)  :: OtherState  !< Input: Other states at t;
                                                                  !!   Output: Other states at t+dt
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: errStat     !< Error status of the operation
   character(*),                    intent(  out)  :: errMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(FVW_InputType)                             :: uInterp     ! Interpolated/Extrapolated input
   integer(intKi)                                  :: ErrStat2    ! temporary Error status
   character(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message
   character(*), parameter                         :: RoutineName = 'FVW_UpdateStates'
      
   ErrStat = ErrID_None
   ErrMsg  = ""

end subroutine FVW_UpdateStates 


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine FVW_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
! NOTE: no matter how many channels are selected for output, all of the outputs are calculated
! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(FVW_InputType),             intent(in   )  :: u           !< Inputs at Time t
   type(FVW_ParameterType),         intent(in   )  :: p           !< Parameters
   type(FVW_ContinuousStateType),   intent(in   )  :: x           !< Continuous states at t
   type(FVW_DiscreteStateType),     intent(in   )  :: xd          !< Discrete states at t
   type(FVW_ConstraintStateType),   intent(in   )  :: z           !< Constraint states at t
   type(FVW_OtherStateType),        intent(in   )  :: OtherState  !< Other states at t
   type(FVW_OutputType),            intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                  !!   nectivity information does not have to be recalculated)
   type(FVW_MiscVarType),           intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*), parameter                         :: RoutineName = 'FVW_CalcOutput'

   ErrStat = ErrID_None
   ErrMsg  = ""

end subroutine FVW_CalcOutput

END MODULE FVW
