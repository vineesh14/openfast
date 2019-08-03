MODULE FVW
   USE NWTC_Library
   USE FVW_Types
   USE FVW_Subs

      ! NOTE:  this is a rough format that AD14 stores airfoil info.  This will need
      !        to be replaced by the AirFoilInfo module when we couple FVW into AD15
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


      ! NOTE:  We are copying several things into parameters here that are passed in as values in InitInp.  The
      !        reason for doing this rather than directly putting them into paramaters by the calling code is
      !        that when this code is ported over to AD15, the data structure on the AD15 side will be different.
      !        So, there may be some additional arithmetic involved at that point to convert to the structure
      !        that is present here, or if we upgrade FVW to different data structures, we may need to introduce
      !        some additional data mangling here.

      ! Set some parameters passed in from AeroDyn
   p%DtAero       =  Interval             ! NOTE: depending how FVW gets configured, we might be able to use a different timestep here as dictated by the input file.
   p%TMax         =  InitInp%TMax

      ! Rotor info
   p%HubHt        =  InitInp%HubHt
   p%HubRad       =  InitInp%HubRad
   p%Radius       =  InitInp%Radius

      ! Blade info
   p%NumBl        =  InitInp%NumBl
   p%NElm         =  InitInp%NElm
   IF (.NOT. ALLOCATED( p%C    )) ALLOCATE ( p%C(    p%NElm ))
   IF (.NOT. ALLOCATED( p%RElm )) ALLOCATE ( p%RElm( p%NElm ))
   p%C            =  InitInp%C
   p%RElm         =  InitInp%RElm

      ! Copy over the airfoil parameters.  In AD15, this will need restructuring to use AirFoilInfo 
   CALL AD14AeroConf_CopyParam( InitInp%AirfoilParm, p%AirfoilParm, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN


      ! Read and parse the input file here to get other parameters and info
   ! Read from InitInp%FVWFileName






      ! Set miscvars
   call AllocAry( m%AofA,        p%NElm, p%NumBl, 'Angle of attack',          ErrStat2, ErrMsg2 );    call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
   call AllocAry( m%W2FVW,       p%NElm, p%NumBl, 'Angle of something',       ErrStat2, ErrMsg2 );    call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
   call AllocAry( m%CLFVW,       p%NElm, p%NumBl, 'Coefficient of lift',      ErrStat2, ErrMsg2 );    call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
   call AllocAry( m%VINDFVW,  3, p%NElm, p%NumBl, 'Induced velocity vector',  ErrStat2, ErrMsg2 );    call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN


      ! Return anything in FVW_InitOutput that should be passed back to the calling code here





      ! Multiturbine things
   p%NumTurbs = 1
   call AllocAry( m%loop,     p%NumTurbs, 'Counter of some sort based on number of turbines.', ErrStat2, ErrMsg2 );  call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName )

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


!FIXME: this routine currently does some time propogation and state updating.  This should only be updating outputs, not states.
!        So, anything in here or below that is updating to the next timestep should be relocated to the updatestates routine
!        if this will ever get migrated to AD15.
SUBROUTINE FVW_CalcSomeStuffThatWasInELEMFRC(P, ALPHA, W2, PITNOW, ErrStat, ErrMess, &
                      J, IBlade, VT, VNW, &
                      VNB, Initial, u, m, Time, VINDFW, phi )
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(FVW_ParameterType),       INTENT(INOUT)  :: p           ! Parameters KS--changed from IN to INOUT
   TYPE(FVW_InputType),           INTENT(IN   )  :: u
   TYPE(FVW_MiscVarType),         INTENT(INOUT)  :: m
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   REAL(ReKi),INTENT(  OUT)   :: ALPHA
   REAL(ReKi),INTENT(INOUT)   :: W2
   REAL(ReKi),INTENT(  OUT)   :: PITNOW 


   REAL(ReKi),INTENT(IN)      :: VNB
   REAL(ReKi),INTENT(IN)      :: VNW
   REAL(DbKi),INTENT(IN)      :: Time    !KS
   REAL(ReKi),INTENT(INOUT)   :: VT
   INTEGER, INTENT(IN)        :: J
   INTEGER, INTENT(IN)        :: IBlade
   LOGICAL,   INTENT(IN)      :: Initial
   REAL(ReKi), intent(out)    :: PHI

   REAL(ReKi)                 :: CLFW                   !KS
   REAL(ReKi),INTENT(  OUT)   :: VINDFW(3)              !KS
   REAL(ReKi)                 :: VN_IND                 !KS
   REAL(ReKi)                 :: VT_IND                 !KS
   REAL(ReKi)                 :: CPITCH                 !KS
   REAL(ReKi)                 :: SPITCH                 !KS
   REAL(ReKi)                 :: Pit_tmp                !KS
   REAL(ReKi)                 :: tmpvector(3)           !KS

   ! Local Variables:
   REAL(ReKi)                 :: Vinduced
   REAL(ReKi)                 :: VN

   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.
   CHARACTER(ErrMsgLen)                      :: ErrMessLcl          ! Error message returned by called routines.

   ErrStat = ErrID_None
   ErrMess = ""



     CLFW = 0.0_ReKi
     CALL FVWtest( J, IBlade, Initial, p, u, m, W2, CLFW, VINDFW, Time) !KMK Added FVW call

     Pit_tmp = 0.d0; SPitch = 0.d0; CPitch = 0.d0; tmpVector = 0.d0; VT_IND = 0.d0; VN_Ind = 0.d0
     Pit_tmp    = -1.d0*ATAN2( -1.0_ReKi*DOT_PRODUCT( p%FVWTurbineComponents%Blade(IBlade)%Orientation(1,:),    &
          u%InputMarkers(IBlade)%Orientation(2,:,J)) , &
          DOT_PRODUCT( p%FVWTurbineComponents%Blade(IBlade)%Orientation(1,:),    &
          u%InputMarkers(IBlade)%Orientation(1,:,J)))
     SPitch    = SIN( Pit_tmp  )
     CPitch    = COS( Pit_tmp  )
     tmpVector = -1.d0*SPitch*u%InputMarkers(IBlade)%Orientation(1,:,J) + CPitch*u%InputMarkers(IBlade)%Orientation(2,:,J)
     VT_IND   =     DOT_PRODUCT( tmpVector, VINDFW)
     tmpVector = 0.d0
     tmpVector =     CPitch*u%InputMarkers(IBlade)%Orientation(1,:,J) + SPitch*u%InputMarkers(IBlade)%Orientation(2,:,J)
     VN_IND    =     DOT_PRODUCT( tmpVector, VINDFW )

     VN=VNW +VNB + VN_IND     !KS -- above, it's -Vinduced; why is it (+) here?; 10.13.15 -- I do think the (+) is correct; it's a derivation thing. 
     VT=VT + VT_IND
     PHI   = ATAN2( VN, VT )

     ALPHA = PHI - PITNOW
     CALL MPI2PI ( ALPHA )
     W2 = VN * VN + VT * VT        !KS -- ! This is calculated in FVW code and then reassigned
                                          ! here without ever being used...why?? Same with
                                          ! ALPHA_TMP (which is just never used) and CLFW

END SUBROUTINE FVW_CalcSomeSTuffThatWasInELEMFRC


END MODULE FVW
