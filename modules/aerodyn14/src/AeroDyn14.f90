!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2012-2016  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
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
!
!**********************************************************************************************************************************
!> Module for the old aerodynamic routines. This module is for loose coupling only, without linearization, because it does not
!! fully conform to the FAST framework. This module will eventually be replaced by AeroDyn (i.e., AeroDyn v15 [aerodyn.f90])
MODULE AeroDyn14

   USE AeroDyn14_Types
   USE AeroSubs
   USE FVW
   USE NWTC_Library


   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: AD14_Ver = ProgDesc( 'AeroDyn14', '', '' )

      ! ..... Public Subroutines ............

   PUBLIC :: AD14_Init                           ! Initialization routine
   PUBLIC :: AD14_End                            ! Ending routine (includes clean up)

   PUBLIC :: AD14_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: AD14_CalcOutput                     ! Routine for computing outputs


CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_Init( InitInp, u, p, x, xd, z, O, y, m, Interval, InitOut, ErrStat, ErrMess )
!..................................................................................................................................
   USE               AeroGenSubs,   ONLY: ElemOpen
   USE DWM
   IMPLICIT NONE

   TYPE(AD14_InitInputType),       INTENT(INOUT)  :: InitInp     ! Input data for initialization routine
   TYPE(AD14_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters !KS -- changed from OUT to INOUT
   TYPE(AD14_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(AD14_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(AD14_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD14_OtherStateType),      INTENT(  OUT)  :: O           ! Initial other states
   TYPE(AD14_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                 !   only the output mesh is initialized)
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables  !KS -- changed from OUT to INOUT
   REAL(DbKi),                     INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                 !   (1) AD14_UpdateStates() is called in loose coupling &
                                                                 !   (2) AD14_UpdateDiscState() is called in tight coupling.
                                                                 !   Input is the suggested time from the glue code;
                                                                 !   Output is the actual coupling interval that will be used
                                                                 !   by the glue code.
   TYPE(AD14_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


      ! Internal variables
   REAL(ReKi)                       :: CosPrecone
   REAL(ReKi)                       :: DTip, ElemRad, Dhub, Rhub     ! variables for calculating hub- and tip-loss constants
   REAL(ReKi)                       :: HubRadius
!   REAL(ReKi)                       :: MeanWind
   REAL(ReKi)                       :: TipRadius
   REAL(ReKi)                       :: TmpVar
   REAL(ReKi)                       :: TmpPos(3)
   REAL(ReKi)                                :: TwrNodeHt                     ! The height of the current tower node.

   INTEGER                          :: IB, IE 
   INTEGER                          :: IELM

!   CHARACTER(1024)                  :: Title

   INTEGER                                   :: Elem                          ! Index for mesh element.
   INTEGER                                   :: InterpIndx                 ! Index telling the interpolation routine where to start in the array.
   INTEGER                                   :: Node                          ! Index used to pull points out of the array of values at given node location.
   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.

   CHARACTER(ErrMsgLen)                      :: ErrMessLcl          ! Error message returned by called routines.
   CHARACTER(*), PARAMETER                   :: RoutineName = 'AD14_Init'

         ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMess  = ""
   InterpIndx = 1

   
   CALL NWTC_Init( )
   
   
         ! Display the module information

   CALL DispNVD( AD14_Ver )
   
   InitOut%Ver = AD14_Ver
   m%FirstWarn = .TRUE.
   !-------------------------------------------------------------------------------------------------
   ! Set up AD variables
   !-------------------------------------------------------------------------------------------------

   p%LinearizeFlag     = .FALSE.             ! InitInp%LinearizeFlag
   p%Blade%BladeLength = InitInp%TurbineComponents%BladeLength
   p%DtAero            = Interval            ! set the default DT here; may be overwritten later, when we read the input file in AD14_GetInput()
   p%UseDWM            = InitInp%UseDWM
!FIXME: this should come from the input file, not from things passed to init..
   p%UseFVW            = InitInp%UseFVW   !!KS

         ! Define parameters here:

   p%WrOptFile   = InitInp%WrSumFile

   p%NumBl   = SIZE( InitInp%TurbineComponents%Blade )
   IF ( p%NumBl < 1 ) THEN
      CALL SetErrStat( ErrID_Fatal,'AeroDyn cannot run without blades in the model.',ErrStat,ErrMess,RoutineName)
      RETURN
   END IF
!bjj: what's the difference between p%NumBl, p%Blade%NB, and InitInp%NumBl?
!MLB: Heck if I know!

         ! Define initial system states here:
   !-------------------------------------------------------------------------------------------------
   ! Read the AeroDyn14 input file and open the output file if requested
   ! bjj: these should perhaps be combined
   !-------------------------------------------------------------------------------------------------
   CALL AD14_GetInput(InitInp, P, x, xd, z, m, y, ErrStatLcl, ErrMessLcl )
      CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
      IF (ErrStat >= AbortErrLev ) RETURN 


      ! allocate variables for aerodyn forces
   p%LinearizeFlag     = .FALSE.

   Interval = p%DtAero
      

   IF ( .NOT. ALLOCATED( m%StoredForces  )) THEN
      CALL AllocAry(m%StoredForces, 3,p%Element%NELM,p%NumBl,'m%StoredForces',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF
   IF ( .NOT. ALLOCATED( m%StoredMoments ))  THEN
      CALL AllocAry(m%StoredMoments, 3,p%Element%NELM,p%NumBl,'m%StoredForces',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF
     
   IF (.NOT. ALLOCATED(m%Element%W2) ) THEN
      CALL AllocAry(m%Element%W2, p%Element%NELM, p%NumBl,'m%Element%W2',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF

   IF (.NOT. ALLOCATED(m%Element%Alpha) ) THEN
      CALL AllocAry(m%Element%Alpha, p%Element%NELM, p%NumBl,'m%Element%Alpha',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF
   IF (ErrStat >= AbortErrLev ) RETURN 

   
   P%UnWndOut = -1
   P%UnElem = -1   
   IF ( p%ElemPrn )  THEN
      CALL ElemOpen ( TRIM( InitInp%OutRootName )//'.AD.out', P, m, ErrStat, ErrMess, AD14_Ver )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
         IF (ErrStat >= AbortErrLev ) RETURN 
   END IF
      

   !-------------------------------------------------------------------------------------------------
   ! Calculate the rotor and hub radaii from the input values
   !-------------------------------------------------------------------------------------------------
   HubRadius = DOT_PRODUCT( InitInp%TurbineComponents%Blade(1)%Position(:)        &
                          - InitInp%TurbineComponents%Hub%Position(:),            &
                            InitInp%TurbineComponents%Blade(1)%Orientation(3,:) )

   DO IB = 2,p%NumBl
      TmpVar    = DOT_PRODUCT( InitInp%TurbineComponents%Blade(IB)%Position(:)    &
                             - InitInp%TurbineComponents%Hub%Position(:),         &
                               InitInp%TurbineComponents%Blade(IB)%Orientation(3,:) )
      IF ( ABS( TmpVar - HubRadius ) > 0.001 ) THEN ! within 1 mm
         CALL ProgWarn( ' AeroDyn\AD14_Init() calculated HubRadius is not the same for all '// &
                           'blades. Using value from blade 1.' )
         EXIT
      END IF
   END DO !IB

   TipRadius = InitInp%TurbineComponents%BladeLength + HubRadius

   CosPrecone = ASIN( DOT_PRODUCT( InitInp%TurbineComponents%Blade(1)%Orientation(3,:), &
                                   InitInp%TurbineComponents%Hub%Orientation(1,:) ) )  ! precone angle -- do COS later

   DO IB = 2,p%NumBl
      TmpVar  = ASIN( DOT_PRODUCT( InitInp%TurbineComponents%Blade(IB)%Orientation(3,:), &
                                   InitInp%TurbineComponents%Hub%Orientation(1,:) ) )
      IF ( ABS( TmpVar - CosPrecone ) > 0.009 ) THEN     ! within ~ 1/2 degree
         CALL ProgWarn( ' AeroDyn\AD14_Init() calculated precone angle is not the same for all'// &
                           ' blades. Using value from blade 1.' )
         EXIT
      END IF
   END DO !IBld

   CosPrecone = COS( CosPrecone )

   p%Blade%R = TipRadius * CosPrecone
   RHub = HubRadius * CosPrecone
   p%HubRad = RHub

      ! Check that the AeroDyn input DR and RElm match (use the HubRadius and TipRadius to verify)
      ! before using them to calculate the tip- and hub-loss constants
   CALL CheckRComp( P, x, xd, z, m, y, ErrStat, ErrMess, &
                    InitInp%ADFileName, HubRadius, TipRadius )

   IF ( ErrStat /= ErrID_None ) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Calculate tip-loss constants
   !-------------------------------------------------------------------------------------------------
   DO IElm = 1,p%Element%NElm  ! Loop through all blade elements

      ElemRad = p%Element%RELM(IElm)*CosPrecone

      IF( ElemRad == 0.0 )  THEN  !BJJ: should this be 0.001 (or another small number) instead of exactly 0.0?
         CALL SetErrStat( ErrID_Fatal,'Error calculating tip loss constant for element '//TRIM(Int2LStr(IElm))//&
                          '. Division by zero.',ErrStat,ErrMess,RoutineName)
         
         RETURN
      ELSE
         DTip         = p%Blade%R - ElemRad
         p%Element%TLCNST(IElm) = 0.5 * p%NumBl * DTip / ElemRad
      ENDIF

   ENDDO             ! IElm - all blade elements


   !-------------------------------------------------------------------------------------------------
   ! Calculate hub-loss constants
   !-------------------------------------------------------------------------------------------------
   IF ( RHub > 0.001 )  THEN

      DO Ielm = 1,p%Element%NELM  ! Loop through all blade elements

         ElemRad = p%Element%RELM(Ielm)*CosPrecone  ! Use only the precone angle of blade 1 (assumed very similar to other blades)

         DHub         = ElemRad - RHub
         p%Element%HLCNST(Ielm) = 0.5 * p%NumBl * DHub / RHub

      ENDDO             ! IELM - all blade elements

   ELSE

      p%Element%HLCNST(:) = 0.0

   ENDIF



      !-------------------------------------------------------------------------------------------------
      ! Interpolate the tower diameter at ElastoDyn's tower nodes if we will be computing tower aerodynamics.
      !-------------------------------------------------------------------------------------------------

   IF ( p%TwrProps%CalcTwrAero )  THEN

         !-------------------------------------------------------------------------------------------------
         ! IMPORTANT NOTES:
         !     o  Supposedly, the glue code will not try to do anything with the tower-aero mesh if is is
         !        not created, so the creation is inside the test for CalcTwrAero.
         !     o  The tower properties from AeroDyn's tower file are for heights from the origin (ground or
         !        MSL) to the hub height--not the top of the tower.
         !     o  For now, we are allowing only one set of Cd for the entire tower.
         !     o  InterpIndx is initialize to 1 at compile time.
         !-------------------------------------------------------------------------------------------------


         ! Create the mesh for the tower aerodynamics.

      CALL MeshCreate ( BlankMesh       = u%Twr_InputMarkers      &
                      , IOS             = COMPONENT_INPUT         &
                      , NNodes          = InitInp%NumTwrNodes     &
                      , Orientation     = .TRUE.                  &
                      , TranslationDisp = .TRUE.                  &
                      , TranslationVel  = .TRUE.                  &
                      , ErrStat         = ErrStatLcl              &
                      , ErrMess         = ErrMessLcl              )

      CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
      IF ( ErrStat >= AbortErrLev )  RETURN


         ! Set the positions of the nodes.  MeshCreate() allocated the Position array.

      DO Node = 1,u%Twr_InputMarkers%Nnodes
         CALL MeshPositionNode ( Mesh  = u%Twr_InputMarkers          &
                                ,INode = Node                        &
                                ,Pos   = InitInp%TwrNodeLocs(:,Node) &  
                                ,ErrStat   = ErrStatLcl              &
                                ,ErrMess   = ErrMessLcl              )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN

      END DO         


         ! Construct the tower with Line-2 elements.

      DO Elem=1,u%Twr_InputMarkers%Nnodes-1

         CALL MeshConstructElement ( Mesh     = u%Twr_InputMarkers &
                                   , Xelement = ELEMENT_LINE2      &
                                   , P1       = Elem               &
                                   , P2       = Elem+1             &
                                   , ErrStat  = ErrStatLcl         &
                                   , ErrMess  = ErrMessLcl         )

         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN

      ENDDO


         ! Commit the mesh to the funny farm.

      CALL MeshCommit ( u%Twr_InputMarkers, ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN


         ! Copy the input mesh to create the output mesh.  Does

      CALL MeshCopy ( SrcMesh  = u%Twr_InputMarkers &
                    , DestMesh = y%Twr_OutputLoads  &
                    , CtrlCode = MESH_SIBLING       &
                    , Force    = .TRUE.             &
                    , ErrStat  = ErrStatLcl         &
                    , ErrMess  = ErrMessLcl         )

         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN


         ! Check to ensure that the user did not specify more than one set of Cd(Re) tables.  Temporary restriction.

      IF ( p%TwrProps%NTwrCD /= 1 )  THEN
         CALL SetErrStat(ErrID_Fatal,'You must have one and only one set of drag coefficients for the AeroDyn tower file.',ErrStat,ErrMess,RoutineName )
         RETURN
      END IF


         ! Build the TwrNodeWidth array.

      p%TwrProps%NumTwrNodes = InitInp%NumTwrNodes

      IF (.NOT. ALLOCATED( p%TwrProps%TwrNodeWidth ) ) THEN
         CALL AllocAry( p%TwrProps%TwrNodeWidth, p%TwrProps%NumTwrNodes, "array for tower widths at ED node locations", ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN
      END IF

      DO Node=1,p%TwrProps%NumTwrNodes

         TwrNodeHt = InitInp%TwrNodeLocs(3,Node)/p%Rotor%HH

         p%TwrProps%TwrNodeWidth(Node) = InterpStp( TwrNodeHt, p%TwrProps%TwrHtFr, p%TwrProps%TwrWid, InterpIndx, p%TwrProps%NTwrHT )

      END DO ! Node

   ELSE
      u%Twr_InputMarkers%Nnodes = 0
      y%Twr_OutputLoads%Nnodes  = 0
   END IF ! ( p%TwrProps%CalcTwrAero )


   !-------------------------------------------------------------------------------------------------
   ! Write the summary (opt) file, then close it
   !-------------------------------------------------------------------------------------------------

   IF (p%WrOptFile) THEN

      CALL ADOut(InitInp, P, m, AD14_Ver, TRIM(InitInp%OutRootName)//'.AD.sum', ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN

   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Initialize the inputs from the wind inflow module
   !-------------------------------------------------------------------------------------------------
   CALL AllocAry( u%InflowVelocity, 3, p%Element%NElm*p%NumBl + u%Twr_InputMarkers%Nnodes, 'u%InflowVelocity', ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN
   u%InflowVelocity = 0.0_ReKi
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Calling the DWM
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   IF ( p%UseDWM ) THEN   
      ! InitInp%DWM%IfW%InputFileName is already set in FAST
         
         ! bjj: all this stuff should be put in DWM_Init.....>
      p%DWM%RR              = p%Blade%R
      p%DWM%BNum            = p%NumBl
      !p%DWM%ElementNum      = m%ElOut%NumElOut  !bjj: NumElOut is the number of elements to be printed in an output file. I really think you want the number of blade elements. I guess we should check that NumElOut is the same as p%Element%NElm
      p%DWM%ElementNum      = p%Element%NElm      ! yj: 1/18/2016
      p%DWM%air_density     = p%Wind%Rho
   
      IF (.NOT. ALLOCATED(m%DWM%Nforce    )) ALLOCATE ( m%DWM%Nforce(    p%Element%NElm,p%NumBl),STAT=ErrStatLcl);CALL SetErrStat(ErrStatLcl, 'Error allocating DWM Nforce array', ErrStat,ErrMess,RoutineName )
      IF (.NOT. ALLOCATED(m%DWM%blade_dr  )) ALLOCATE ( m%DWM%blade_dr(  p%Element%NElm),        STAT=ErrStatLcl);CALL SetErrStat(ErrStatLcl, 'Error allocating DWM blade_dr array', ErrStat,ErrMess,RoutineName )
      IF (.NOT. ALLOCATED(p%DWM%ElementRad)) ALLOCATE ( p%DWM%ElementRad(p%Element%NElm),        STAT=ErrStatLcl);CALL SetErrStat(ErrStatLcl, 'Error allocating DWM ElementRad array', ErrStat,ErrMess,RoutineName )
      if (errStat >= AbortErrLev) return
      
      m%DWM%blade_dr = p%Blade%DR(:)
      p%DWM%ElementRad = p%Element%RELM(:)   
   
      CALL DWM_Init( InitInp%DWM, m%DWM_Inputs, p%DWM, x%DWM, xd%DWM, z%DWM, O%DWM, m%DWM_Outputs, m%DWM, Interval, InitOut%DWM, ErrStatLcl, ErrMessLcl)
   

      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      
   END IF !UseDWM
   


   !-------------------------------------------------------------------------------------------------
   ! Turn off dynamic inflow for wind less than 8 m/s (per DJL: 8 m/s is really just an empirical guess)
   ! DJL: Comment out this code when using new proposed GDW check in ELEMFRC
   ! BJJ: FIX THIS!!!!
   !-------------------------------------------------------------------------------------------------

   ! BJJ: can't put this here b/c we need InitInp%MWS from InflowWind
   !IF (p%DynInfl) THEN
   !            
   !   IF ( InitInp%MWS  < 8.0 ) THEN
   !      p%DynInfl = .FALSE.
   !      CALL SetErrStat(ErrID_Info,'Estimated average inflow wind speed is less than 8 m/s. Dynamic Inflow will be turned off.',ErrStat,ErrMess,RoutineName )
   !      IF ( ErrStat >= AbortErrLev )  RETURN      
   !   END IF
   !
   !ENDIF   
   

   !-------------------------------------------------------------------------------------------------
   ! Set initial guesses for inputs:
   !-------------------------------------------------------------------------------------------------
   
   !..........
   ! u%TurbineComponents
   !..........

   CALL AD14AeroConf_CopyInput( InitInp%TurbineComponents, u%TurbineComponents, MESH_NEWCOPY, ErrStatLcl, ErrMessLcl )
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
   
   !..........
   ! u%InputMarkers (blade meshes):
   !..........
   
   ALLOCATE( u%InputMarkers(p%NumBl), STAT=ErrStatLcl )
   IF (ErrStatLcl /= 0 ) THEN
      CALL SetErrStat ( ErrID_Fatal, 'Could not allocate u%InputMarkers (meshes)', ErrStat,ErrMess,RoutineName )
      RETURN
   END IF


   DO IB = 1, p%NumBl
      CALL MeshCreate( BlankMesh      = u%InputMarkers(IB)    &
                     ,IOS            = COMPONENT_INPUT        &
                     ,NNodes         = p%Element%NELM         &
                     ,Orientation    = .TRUE.                 &
                     ,TranslationVel = .TRUE.                 &
                     ,TranslationAcc = .TRUE.                 &  !bjj: added for MHK turbines
                     ,RotationVel    = .TRUE.                 &
                     ,nScalars       = 2                      &  ! scalar 1 is W, scalar 2 is Alpha
                     ,ErrStat        = ErrStatLcl             &
                     ,ErrMess        = ErrMessLcl             )

      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      
      ! create the elements
      DO IE = 1, p%Element%NELM-1 ! construct the blades into Line2 elements
         CALL MeshConstructElement ( Mesh = u%InputMarkers(IB)    &
                                  ,Xelement = ELEMENT_LINE2       &
                                  ,P1       = IE                  &
                                  ,P2       = IE+1                &
                                  ,ErrStat  = ErrStatLcl          &
                                  ,ErrMess  = ErrMessLcl          )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         
      ENDDO
     
      ! position/orient the nodes
      DO IE = 1, p%Element%NELM
         TmpPos(1) = 0.
         TmpPos(2) = 0.
         TmpPos(3) = p%Element%Relm(IE) - HubRadius
         CALL MeshPositionNode ( Mesh = u%InputMarkers(IB)              &
                                 ,INode = IE                            &
                                 ,Pos= TmpPos                           &  ! this info comes from FAST (not yet)
                                 ,ErrStat   = ErrStatLcl                &
                                 ,ErrMess   = ErrMessLcl                )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN

         ! RELATIVE ORIENTATION OF BLADE ELEMENTS
         u%InputMarkers(IB)%Orientation(1,1,IE) = COS( P%Element%TWIST(IE) )
         u%InputMarkers(IB)%Orientation(2,1,IE) = SIN( P%Element%TWIST(IE) )
         u%InputMarkers(IB)%Orientation(3,1,IE) = SIN( P%Element%TWIST(IE) )
         u%InputMarkers(IB)%Orientation(1,2,IE) = -1. * u%InputMarkers(IB)%Orientation(2,1,IE)
         u%InputMarkers(IB)%Orientation(2,2,IE) =       u%InputMarkers(IB)%Orientation(1,1,IE)
         u%InputMarkers(IB)%Orientation(3,2,IE) = 0.0
         u%InputMarkers(IB)%Orientation(1,3,IE) = 0.0
         u%InputMarkers(IB)%Orientation(2,3,IE) = 0.0
         u%InputMarkers(IB)%Orientation(3,3,IE) = 1.0
      ENDDO
     
       CALL MeshCommit ( Mesh = u%InputMarkers(IB)    &
                        ,ErrStat  = ErrStatLcl        &
                        ,ErrMess  = ErrMessLcl        )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
          
   ENDDO
         

   !..........
   ! u%Twr_InputMarkers (tower meshes):
   !..........
     
   !bjj: done above in section for IF Tower Loads is on
   

   !..........
   ! u%MulTabLoc:
   !..........
   
   IF (.NOT. ALLOCATED(u%MulTabLoc)) THEN
      ALLOCATE( u%MulTabLoc(p%Element%NELM, p%NumBl), STAT = ErrStatLcl )
      IF (ErrStatLcl /= 0) THEN
         CALL SetErrStat ( ErrID_Fatal, 'Could not allocate u%MulTabLoc', ErrStat,ErrMess,RoutineName )
         RETURN
      END IF
   END IF

   u%MulTabLoc(:,:) = 0.0
   
   
   !-------------------------------------------------------------------------------------------------
   ! Allocate space for outputs and set up output meshes:
   !-------------------------------------------------------------------------------------------------
   
   !..........
   ! y%OutputLoads (blade meshes):
   !..........
   

   ALLOCATE( y%OutputLoads(p%NumBl), STAT = ErrStatLcl )
      IF (ErrStatLcl /= 0) THEN
         CALL SetErrStat ( ErrID_Fatal, 'Could not allocate y%OutputLoads (meshes)', ErrStat,ErrMess,RoutineName )
         RETURN
      END IF
   
   DO IB = 1, p%NumBl

       CALL MeshCopy ( SrcMesh  = u%InputMarkers(IB)  &
                      ,DestMesh = y%OutputLoads(IB)   &
                      ,CtrlCode = MESH_SIBLING        &
                      ,Force    = .TRUE.              &
                      ,Moment   = .TRUE.              &
                      ,ErrStat  = ErrStatLcl          &
                      ,ErrMess  = ErrMessLcl          )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
   ENDDO
 

   !-------------------------------------------------------------------------------------------------
   ! Initialize FVW module if it is used
   !-------------------------------------------------------------------------------------------------
   if (p%UseFVW ) then

         ! Copy some things to the InitInp.  When FVW is incorporated into a different module, there may
         ! be some additional logic necessary to put it into the correct form for FVW to use (AD15 stores
         ! things differently)
      InitInp%FVW%TMax        = p%TMax
      InitInp%FVW%DtAero      = p%DTAero        ! This could be different than the Interval passed in...
      InitInp%FVW%Radius      = p%Blade%R
      InitInp%FVW%HubHt       = InitInp%HubHt
      InitInp%FVW%HubRad      = p%HubRad
      InitInp%FVW%NumBl       = InitInp%NumBl
      InitInp%FVW%NElm        = p%Element%NElm
      IF (.NOT. ALLOCATED( InitInp%FVW%C    )) ALLOCATE ( InitInp%FVW%C(    p%Element%NElm ))
      IF (.NOT. ALLOCATED( InitInp%FVW%RElm )) ALLOCATE ( InitInp%FVW%RElm( p%Element%NElm ))
      InitInp%FVW%RElm        = p%Element%RElm
      InitInp%FVW%C           = p%Blade%C

         ! The rest of the FVW settings should come from its input file, so pass input file name
      InitInp%FVW%FVWFileName = InitInp%FVWFileName

         ! copy some of the airfoil parameters needed
      CALL AD14AeroConf_CopyParam( p%Airfoil, InitInp%FVW%AirfoilParm, MESH_NEWCOPY, ErrStatLcl, ErrMessLcl )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN

      call FVW_Init( InitInp%FVW, u%FVW, p%FVW, x%FVW, xd%FVW, z%FVW, O%FVW, y%FVW, m%FVW, Interval, InitOut%FVW, ErrStatLcl, ErrMessLcl )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )

         ! If anything is passed back in InitOut%FVW, deal with it here...


         !FIXME   This really probably should be done inside of FVW_Init instead of here.
         !        Not entirely sure how to pass the u%InputMarkers in though.
      if ( p%UseFVW ) then 
         ALLOCATE( u%FVW%InputMarkers(p%NumBl), STAT = ErrStatLcl )
            IF (ErrStatLcl /= 0) THEN
               CALL SetErrStat ( ErrID_Fatal, 'Could not allocate u%FVW%InputMarkers (meshes)', ErrStat,ErrMess,RoutineName )
               RETURN
            END IF

         DO IB = 1, p%NumBl
             CALL MeshCopy ( SrcMesh  = u%InputMarkers(IB)  &
                            ,DestMesh = u%FVW%InputMarkers(IB) &
                            ,CtrlCode = MESH_COUSIN         &
                            ,Orientation    = .TRUE.        &
                            ,TranslationVel = .TRUE.        &
                            ,RotationVel    = .TRUE.        &
                            ,ErrStat  = ErrStatLcl          &
                            ,ErrMess  = ErrMessLcl          )
               CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN
         ENDDO
      endif 
   endif



   !..........
   ! y%Twr_OutputLoads (tower meshes):
   !..........

   !bjj: done above in section for IF Tower Loads is on
   
     
   !-------------------------------------------------------------------------------------------------
   ! Initialize AeroDyn variables not initialized elsewhere (except in module initialization)
   ! and return
   !-------------------------------------------------------------------------------------------------
   m%InducedVel%SumInfl  = 0.0_ReKi
   m%Rotor%AvgInfl       = 0.0_ReKi
   m%OldTime             = 0.0_DbKi
   m%SuperSonic          = .FALSE.   
   m%NoLoadsCalculated   = .TRUE.

   p%TwoPiNB     = TwoPi / REAL( p%NumBl, ReKi )
           
   
   DO ie = 1, maxInfl
      p%DynInflow%xMinv(ie) = PIBY2 / hfunc(MRvector(ie), NJvector(ie))   !bjj: this is really just a Fortran parameter, too.
   END DO !ie   

   
   InitOut%AirDens = p%Wind%Rho


   
   
   RETURN

END SUBROUTINE AD14_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMess )
! This routine is called at the end of the simulation.
!..................................................................................................................................
      USE DWM_Types
      USE DWM

      TYPE(AD14_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(AD14_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(AD14_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(AD14_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other states
      TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


         ! Place any last minute operations or calculations here:
         
      IF (p%UseDWM ) THEN
         !----- Call the DWM ------- 
      
         CALL DWM_End( m%DWM_Inputs, p%DWM, x%DWM, xd%DWM, z%DWM, OtherState%DWM, m%DWM_Outputs, m%DWM, ErrStat, ErrMess )
      END IF ! UseDWM
      
      !--------------------------


      IF (p%UseFVW )    CALL FVW_End( u%FVW, p%FVW, x%FVW, xd%FVW, z%FVW, OtherState%FVW, y%FVW, m%FVW, ErrStat, ErrMess )

      !--------------------------


         ! Close files here:

      ! AD14_IOParams
   IF (P%UnEc > 0)    CLOSE(P%UnEc) ! not currently used

   IF (P%UnWndOut > 0) CLOSE(P%UnWndOut)
   IF (P%UnElem   > 0) CLOSE(P%UnElem)

         ! Destroy the input data:

      CALL AD14_DestroyInput( u, ErrStat, ErrMess )


         ! Destroy the parameter data:

      CALL AD14_DestroyParam( p, ErrStat, ErrMess )


         ! Destroy the state data:

      CALL AD14_DestroyContState(   x,           ErrStat, ErrMess )
      CALL AD14_DestroyDiscState(   xd,          ErrStat, ErrMess )
      CALL AD14_DestroyConstrState( z,           ErrStat, ErrMess )
      CALL AD14_DestroyOtherState(  OtherState,  ErrStat, ErrMess )

      CALL AD14_DestroyMisc( m, ErrStat, ErrMess )

         ! Destroy the output data:

      CALL AD14_DestroyOutput( y, ErrStat, ErrMess )


END SUBROUTINE AD14_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMess )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................

      REAL(DbKi),                           INTENT(IN   ) :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                       INTENT(IN   ) :: n           ! Current simulation time step n = 0,1,...
      TYPE(AD14_InputType),                 INTENT(INOUT) :: u(:)        ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                           INTENT(IN   ) :: utimes(:)   ! Times associated with u(:), in seconds
      TYPE(AD14_ParameterType),             INTENT(IN   ) :: p           ! Parameters
      TYPE(AD14_ContinuousStateType),       INTENT(INOUT) :: x           ! Input: Continuous states at t;
                                                                         !   Output: Continuous states at t + Interval
      TYPE(AD14_DiscreteStateType),         INTENT(INOUT) :: xd          ! Input: Discrete states at t;
                                                                         !   Output: Discrete states at t  + Interval
      TYPE(AD14_ConstraintStateType),       INTENT(INOUT) :: z           ! Input: Constraint states at t;
                                                                         !   Output: Constraint states at t + Interval
      TYPE(AD14_OtherStateType),            INTENT(INOUT) :: OtherState  ! Input: Other states at t;
                                                                         !   Output: Other states at t + Interval
      TYPE(AD14_MiscVarType),               INTENT(INOUT) :: m           ! Misc/optimization variables
      INTEGER(IntKi),                       INTENT(  OUT) :: ErrStat     ! Error status of the operation
      CHARACTER(*),                         INTENT(  OUT) :: ErrMess     ! Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(AD14_ContinuousStateType)                 :: dxdt        ! Continuous state derivatives at Time
      TYPE(AD14_ConstraintStateType)                 :: z_Residual  ! Residual of the constraint state equations (Z)

!      INTEGER(IntKi)                                    :: ErrStat2    ! Error status of the operation (occurs after initial error)
!      CHARACTER(ErrMsgLen)                              :: ErrMess2     ! Error message if ErrStat2 /= ErrID_None

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


   ! AeroDyn v14 DOES actually have states, but they are updated in CalcOutput because no one ever took the time to 
   ! identify which variables are states.
      
END SUBROUTINE AD14_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_CalcOutput( Time, u, p, x, xd, z, O, y, m, ErrStat, ErrMess )
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................
   
      USE               AeroGenSubs,   ONLY: ElemOut
      USE               DWM_Types
      USE               DWM
!FIXME: remove InflowWind from here...
      USE InflowWind !! KS

      REAL(DbKi),                     INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD14_InputType),           INTENT(INOUT)  :: u           ! Inputs at Time !KS changed from IN to INOUT
      TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters !KS changed from IN to INOUT
      TYPE(AD14_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD14_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(AD14_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(AD14_OtherStateType),      INTENT(IN   )  :: O           ! Other states at Time
      TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
      TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


      ! Local variables
   REAL(DbKi), PARAMETER      :: OnePlusEpsilon = 1 + EPSILON(Time)

   REAL(ReKi)                 :: VNElement
   REAL(ReKi)                 :: VelNormalToRotor2
   REAL(ReKi)                 :: VNWind
   REAL(ReKi)                 :: VTTotal
   REAL(ReKi)                 :: DFN
   REAL(ReKi)                 :: DFT
   REAL(ReKi)                 :: PMA
   REAL(ReKi)                 :: SPitch                     ! sine of PitNow
   REAL(ReKi)                 :: CPitch                     ! cosine of PitNow
   REAL(ReKi)                 :: Phi                        ! Local value of Phi

   REAL(ReKi)                 :: AvgVelNacelleRotorFurlYaw
   REAL(ReKi)                 :: AvgVelTowerBaseNacelleYaw
   REAL(ReKi)                 :: AvgVelTowerBaseYaw
   REAL(ReKi)                 :: AzimuthAngle
   REAL(ReKi)                 :: rNacelleHub   (2)
   REAL(ReKi)                 :: rLocal
   REAL(ReKi)                 :: rRotorFurlHub (2)
   REAL(ReKi)                 :: rTowerBaseHub (2)

   REAL(ReKi)                 :: tmpVector     (3)
   REAL(ReKi)                 :: VelocityVec   (3)

   INTEGER                    :: ErrStatLcL        ! Error status returned by called routines.
   INTEGER                    :: IBlade
   INTEGER                    :: IElement
   INTEGER                    :: Node              ! Node index.

   INTEGER                    :: I
   CHARACTER(ErrMsgLen)       :: ErrMessLcl          ! Error message returned by called routines.

   INTEGER                    :: LoopNum     !KS
   LOGICAL                    :: FWAKE             ! Are we doing wake calculations this loop?
   CHARACTER(*), PARAMETER                   :: RoutineName = 'AD14_AeroSubs' !KS Not sure why I added this

    REAL(ReKi)                 :: VIND_FVW(3)   !KS

   ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMess  = ""


   !-------------------------------------------------------------------------------------------------
   ! Determine if loads should be recalculated or just returned
   !-------------------------------------------------------------------------------------------------
      ! NOTE: Time is scaled by OnePlusEps to ensure that loads are calculated at every
   !       time step when DTAero = DT, even in the presence of numerical precision errors.

   p%FVW%Time = Time   ! KS
   IF ( m%NoLoadsCalculated .OR. ( Time*OnePlusEpsilon - m%OldTime ) >= p%DTAERO )  THEN
         ! It's time to update the aero forces

         ! First we reset the DTAERO parameters for next time
      m%DT      = Time - m%OldTime     !bjj: DT = 0 on first step,
                                       !but the subroutines that use DT check for NoLoadsCalculated (or time > 0)
      m%OldTime = Time

   ELSE IF ( .NOT. p%LinearizeFlag ) THEN

         ! Return the previously-calculated loads

!      CurrentOutputs = ADCurrentLoads

      DO IBlade=1,p%NumBl
       DO IElement=1,p%Element%Nelm
         y%OutputLoads(IBlade)%Force(:,IElement)  = m%StoredForces(:,IElement,IBlade)
         y%OutputLoads(IBlade)%Moment(:,IElement) = m%StoredMoments(:,IElement,IBlade)
       ENDDO
      ENDDO

      IF ( m%FirstWarn ) THEN
         CALL SetErrStat ( ErrID_Warn, 'AeroDyn was designed for an explicit-loose coupling scheme. '//&
            'Using last calculated values from AeroDyn on all subsequent calls until time is advanced. '//&
            'Warning will not be displayed again.', ErrStat,ErrMess,'AD14_CalcOutput' )
         m%FirstWarn = .FALSE.       
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF
      END IF
            
      RETURN

   ENDIF

   
   !-------------------------------------------------------------------------------------------------
   ! Calculate the forces and moments for the blade: SUBROUTINE AeroFrcIntrface( FirstLoop, JElemt, DFN, DFT, PMA )
   !-------------------------------------------------------------------------------------------------

      ! calculate rotor speed
      ! note: Subtracting the RotorFurl rotational velocity for REVS is needed to get the
      ! same answers as before v13.00.00. RotorFurl shouldn't be needed.

   m%Rotor%REVS = ABS( DOT_PRODUCT( u%TurbineComponents%Hub%RotationVel(:) - u%TurbineComponents%RotorFurl%RotationVel(:), &
                                    u%TurbineComponents%Hub%Orientation(1,:) ) )


      ! calculate yaw angle
      ! note: YawAng should use the Hub instead of the RotorFurl, but it is calculated this way to
      ! get the same answers as previous version.
   m%Rotor%YawAng = ATAN2( -1.*u%TurbineComponents%RotorFurl%Orientation(1,2), u%TurbineComponents%RotorFurl%Orientation(1,1) )
   m%Rotor%SYaw   = SIN( m%Rotor%YawAng )
   m%Rotor%CYaw   = COS( m%Rotor%YawAng )

      ! tilt angle
      ! note: tilt angle should use the Hub instead of RotorFurl, but it needs hub to get the same
      ! answers as the version before v13.00.00

   m%Rotor%Tilt = ATAN2( u%TurbineComponents%RotorFurl%Orientation(1,3), &
                         SQRT( u%TurbineComponents%RotorFurl%Orientation(1,1)**2 + &
                         u%TurbineComponents%RotorFurl%Orientation(1,2)**2 ) )

   m%Rotor%CTilt     = COS( m%Rotor%Tilt )
   m%Rotor%STilt     = SIN( m%Rotor%Tilt )


      ! HubVDue2Yaw - yaw velocity due solely to yaw

  AvgVelNacelleRotorFurlYaw = u%TurbineComponents%RotorFurl%RotationVel(3) - u%TurbineComponents%Nacelle%RotationVel(3)
  AvgVelTowerBaseNacelleYaw = u%TurbineComponents%Nacelle%RotationVel(3)   - u%TurbineComponents%Tower%RotationVel(3)
  AvgVelTowerBaseYaw        = u%TurbineComponents%Tower%RotationVel(3)

  rRotorFurlHub(1:2)        = u%TurbineComponents%Hub%Position(1:2) - u%TurbineComponents%RotorFurl%Position(1:2)
  rNacelleHub(1:2)          = u%TurbineComponents%Hub%Position(1:2) - u%TurbineComponents%Nacelle%Position(1:2)
  rTowerBaseHub(1:2)        = u%TurbineComponents%Hub%Position(1:2) - u%TurbineComponents%Tower%Position(1:2)

  m%Rotor%YawVel =   ( AvgVelNacelleRotorFurlYaw * rRotorFurlHub(2) + AvgVelTowerBaseNacelleYaw * rNacelleHub(2) &
                         + AvgVelTowerBaseYaw * rTowerBaseHub(2) ) * m%Rotor%SYaw &
                  - ( AvgVelNacelleRotorFurlYaw * rRotorFurlHub(1) + AvgVelTowerBaseNacelleYaw * rNacelleHub(1) &
                         + AvgVelTowerBaseYaw * rTowerBaseHub(1) ) * m%Rotor%CYaw


   !.................................................................................................
   ! start of NewTime routine
   !.................................................................................................

   m%Rotor%AvgInfl = m%InducedVel%SumInfl * 2.0 / (p%Blade%R*p%Blade%R*p%NumBl)  ! Average inflow from the previous time step
   m%InducedVel%SumInfl = 0.0   ! reset to sum for the current time step

   CALL DiskVel(Time, P, m, u%AvgInfVel, ErrStatLcl, ErrMessLcl)  ! Get a sort of "Average velocity" - sets a bunch of stored variables...
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/DiskVel' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
         
   IF ( P%DStall ) CALL BedUpdate( m )

   ! Enter the dynamic inflow routines here

   IF ( p%Wake )  THEN
      CALL Inflow(Time, P, m, ErrStatLcl, ErrMessLcl)
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/Inflow' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF   
   END IF      
   
       !bjj: perhaps we should send NoLoadsCalculated to initialize dynamic inflow [subroutine Infinit()]
       !bjj: instead of the check that time > 0...?

   !.................................................................................................
   ! end of NewTime routine
   !.................................................................................................

print*,'    CalcOutput...'

   LoopNum = 1 !KS
DO LoopNum = 1, 2 !KS   MOVE OVER AFTER DONE ADDING LINES
   IF (LoopNum .EQ. 2 ) THEN! KS
      FWAKE=.TRUE. !KS
   ELSE !KS
      FWAKE=.FALSE. !KS
   ENDIF  !KS
   Node = 0
   DO IBlade = 1,p%NumBl

         ! calculate the azimuth angle ( we add pi because AeroDyn defines 0 as pointing downward)
         ! note: the equation below should use TurbineComponents%Blade markers, but this is used to get the
         ! same answers as the previous version (before v13.00.00)

      AzimuthAngle = ATAN2( -1.*DOT_PRODUCT( u%TurbineComponents%Hub%Orientation(3,:),         &
                                             u%TurbineComponents%RotorFurl%Orientation(2,:) ), &
                                DOT_PRODUCT( u%TurbineComponents%Hub%Orientation(3,:),         &
                                             u%TurbineComponents%RotorFurl%Orientation(3,:) )  ) + pi + (IBlade - 1)*p%TwoPiNB



      DO IElement = 1,p%Element%NElm

            ! calculate element pitch

         m%Element%PitNow    = -1.*ATAN2( -1.*DOT_PRODUCT( u%TurbineComponents%Blade(IBlade)%Orientation(1,:),    &
                                                           u%InputMarkers(IBlade)%Orientation(2,:,IElement) ) , &
                                              DOT_PRODUCT( u%TurbineComponents%Blade(IBlade)%Orientation(1,:),    &
                                                           u%InputMarkers(IBlade)%Orientation(1,:,IElement) )   )

         SPitch    = SIN( m%Element%PitNow )
         CPitch    = COS( m%Element%PitNow )


            ! calculate distance between hub and element

         tmpVector = u%InputMarkers(IBlade)%Position(:,IElement) - u%TurbineComponents%Hub%Position(:)
         rLocal = SQRT(   DOT_PRODUCT( tmpVector, u%TurbineComponents%Hub%Orientation(2,:) )**2  &
                        + DOT_PRODUCT( tmpVector, u%TurbineComponents%Hub%Orientation(3,:) )**2  )

            ! determine if MulTabLoc should be set.  
     
         IF (.not. p%Reynolds) m%AirFoil%MulTabLoc = u%MulTabLoc(IElement,IBlade)
         
         !-------------------------------------------------------------------------------------------
         ! Get wind velocity components; calculate velocity normal to the rotor squared
         ! Save variables for printing in a file later;
         !-------------------------------------------------------------------------------------------
         Node = Node + 1
         VelocityVec(:)    = AD_WindVelocityWithDisturbance( Time, u, p, x, xd, z, m, y, ErrStatLcl, ErrMessLcl, &
                                                             u%InputMarkers(IBlade)%Position(:,IElement), &
                                                             u%InflowVelocity(:,Node) )
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput' )
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF
         
         !-------------------------------------------------------------------------------------------
         ! DWM wind input update phase 1
         !-------------------------------------------------------------------------------------------
         IF (p%UseDWM) THEN
            !bjj: FIX THIS!!!!         
            !bjj: where do p%DWM%RTPD%SimulationOrder_index and p%DWM%RTPD%upwindturbine_number get set?
         
            IF ( p%DWM%RTPD%SimulationOrder_index > 1) THEN
               IF(  p%DWM%RTPD%upwindturbine_number /= 0 ) THEN
                       
                  m%DWM%position_y = u%InputMarkers(IBlade)%Position(2,IElement)
                                         
                  m%DWM%position_z = u%InputMarkers(IBlade)%Position(3,IElement)
              
                  m%DWM%velocity_wake_mean = 1
              
                  DO I = 1,p%DWM%RTPD%upwindturbine_number
                     m%DWM%DWM_tb%Aerodyn_turbine_num = I
                 
                     CALL   DWM_phase1( Time, m%DWM_Inputs, p%DWM, x%DWM, xd%DWM, z%DWM, &
                                           O%DWM, m%DWM_Outputs, m%DWM, ErrStatLcl, ErrMessLcl )
                 
                     CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/DWM_phase1' )
                     IF (ErrStat >= AbortErrLev) THEN
                        CALL CleanUp()
                        RETURN
                     END IF   
                                                           
                     m%DWM%velocity_wake_mean = (1-((1-m%DWM%velocity_wake_mean)**2 + (1-m%DWM%shifted_velocity_aerodyn)**2)**0.5)
                  END DO
              
                  m%DWM%velocity_wake_mean    = m%DWM%velocity_wake_mean * p%DWM%Wind_file_Mean_u
              
                  VelocityVec(1) = (VelocityVec(1) - p%DWM%Wind_file_Mean_u)*(m%DWM_Inputs%Upwind_result%upwind_small_TI(1)/p%DWM%TI_amb) &
                                  + m%DWM%velocity_wake_mean
              
               END IF
            END IF
                                     
           !------------------------DWM PHASE 2-----------------------------------------------
            IF (Time > 50.00 ) THEN
               m%DWM%U_velocity           = VelocityVec(1)
               m%DWM%V_velocity           = VelocityVec(2)
               m%DWM%NacYaw               = m%Rotor%YawAng 
               m%DWM%DWM_tb%Blade_index   = IBlade
               m%DWM%DWM_tb%Element_index = IElement    

               CALL   DWM_phase2( Time, m%DWM_Inputs, p%DWM, x%DWM, xd%DWM, z%DWM, O%DWM, m%DWM_Outputs, m%DWM, ErrStatLcl, ErrMessLcl )
            
                     CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/DWM_phase1' )
                     IF (ErrStat >= AbortErrLev) THEN
                        CALL CleanUp()
                        RETURN
                     END IF   
            
         
               !CALL CalVelScale(VelocityVec(1),VelocityVec(2),m%DWM_Outputs,z%DWM)
         
               !CALL turbine_average_velocity( VelocityVec(1), IBlade, IElement, m%DWM_Outputs,x%DWM,z%DWM)
            END IF
         END IF ! UseDWM
            
         !-----------------------------------------------------------------------------------------------------------------------
         
         
         
         
         VelNormalToRotor2 = ( VelocityVec(3) * m%Rotor%STilt + (VelocityVec(1) * m%Rotor%CYaw               &
                             - VelocityVec(2) * m%Rotor%SYaw) * m%Rotor%CTilt )**2

         !-------------------------------------------------------------------------------------------
         ! reproduce GetVNVT routine:
         !-------------------------------------------------------------------------------------------
         tmpVector =  -1.*SPitch*u%InputMarkers(IBlade)%Orientation(1,:,IElement) &
                        + CPitch*u%InputMarkers(IBlade)%Orientation(2,:,IElement)
         VTTotal   =     DOT_PRODUCT( tmpVector, VelocityVec - u%InputMarkers(IBlade)%TranslationVel(:,IElement)  )

         tmpVector =     CPitch*u%InputMarkers(IBlade)%Orientation(1,:,IElement) &
                       + SPitch*u%InputMarkers(IBlade)%Orientation(2,:,IElement)
         VNWind    =     DOT_PRODUCT( tmpVector, VelocityVec )
         VNElement = -1.*DOT_PRODUCT( tmpVector, u%InputMarkers(IBlade)%TranslationVel(:,IElement ) )

         !-------------------------------------------------------------------------------------------
         ! Get blade element forces and induced velocity
         !-------------------------------------------------------------------------------------------
         IF (p%UseFVW ) THEN! .AND. p%FVW%FVWInit) THEN
            CALL AD14AeroConf_CopyInput( u%TurbineComponents, p%FVW%FVWTurbineComponents, MESH_NEWCOPY, ErrStatLcl, ErrMessLcl )
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
            IF (ErrStat >= AbortErrLev) RETURN
            p%FVW%FVWTurbineComponents = u%TurbineComponents

               ! NOTE:  this isn't really being used as a full mesh, so we only set a few things.
               !        also, if we do a direct copy, we end up with fatal errors at exit since the
               !        sibling/cousin status will be incorrect
            u%FVW%InputMarkers(IBlade)%Position       = u%InputMarkers(IBlade)%Position
            u%FVW%InputMarkers(IBlade)%Orientation    = u%InputMarkers(IBlade)%Orientation
            u%FVW%InputMarkers(IBlade)%TranslationVel = u%InputMarkers(IBlade)%TranslationVel
         END IF

         IF ( .NOT. FWAKE ) THEN
            CALL ELEMFRC( p, m, ErrStatLcl, ErrMessLcl,                             &
                          AzimuthAngle, rLocal, IElement, IBlade, VelNormalToRotor2, VTTotal, VNWind, &
                        VNElement, m%NoLoadsCalculated, phi )
         ENDIF
         IF ( FWAKE ) THEN
            ! The following should moved to FVW_Init
            p%FVW%RotSpeed    = p%RotSpeed
            IF ( p%FVW%FVWInit ) THEN
!FIXME: InflowWind should not be inside AD14 at all...
               ! The following should moved to FVW_Init
                  CALL InflowWind_Init( p%FVW_WindInit, p%FVW%FVW_Wind%InputData, p%FVW%FVW_Wind%ParamData, p%FVW%FVW_Wind%ContData, p%FVW%FVW_Wind%DiscData, p%FVW%FVW_Wind%ConstrData, p%FVW%FVW_Wind%OtherData, &
                            p%FVW%FVW_Wind%OutputData, p%FVW%FVW_Wind%MiscData, p%IfW_DT, p%FVW%FVW_Wind%InitOutputData, ErrStat, ErrMess )
!FIXME:  step 1:  move this into the FVW_CalcSomeStuff routine
!        step 2:  Sort the data needed here to something that can be passed out to the glue code for the InflowWind call... this will require moving some of the FVW components into an FVW_UpdateStates routine instead.
               CALL InflowWind_CalcOutput( Time, p%FVW%FVW_Wind%InputData, p%FVW%FVW_Wind%ParamData, p%FVW%FVW_Wind%ContData, &
                   & p%FVW%FVW_Wind%DiscData, p%FVW%FVW_Wind%ConstrData, p%FVW%FVW_Wind%OtherData, p%FVW%FVW_Wind%OutputData, p%FVW%FVW_Wind%MiscData, &
                   & ErrStat, ErrMess )
               p%FVW%FVWInit = .FALSE.
            END IF !Initial

            ! m%AirFoil%MulTabLoc and m%AirFoil%PMC may have been updated (the rest of the info appears to be static... AD14 is a mess... --ADP)
            p%FVW%AirfoilOut  = m%Airfoil
!FIXME: copy any additional data that is needed in FVW to u%FVW%... and m%FVW%...   NEVER EVER CHANGE PARAMETERS DURING THE SIMULATION!!!!!
               ! Long term, this should be a call to FVW_CalcOutput (or similar like function), and not be doing any state updating...
            CALL FVW_CalcSomeStuffThatWasInELEMFRC( p%FVW, m%Element%ALPHA(IElement,IBlade), m%Element%W2(IElement,IBlade), m%Element%PITNOW, ErrStatLcl, ErrMessLcl,          &
                        IElement, IBlade, VTTotal, VNWind, &
                        VNElement, m%NoLoadsCalculated, u%FVW, m%FVW, Time, VIND_FVW, phi )
               ! Copy over any outputs (y%FVW%) or miscvars (m%FVW%) needed by AD14 and anything else here
         ENDIF
         CALL ELEMFRC2( p, m, ErrStatLcl, ErrMessLcl, IElement, IBlade, &
                     DFN, DFT, PMA, m%NoLoadsCalculated, phi )
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput' )
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF

         IF ( p%UseFVW ) THEN
            VelocityVec = VelocityVec+VIND_FVW
         END IF
            
         !-------------------------------------------------------------------------------------------
         ! Set up dynamic inflow parameters
         !-------------------------------------------------------------------------------------------
         IF ( p%DynInfl .OR. m%DynInit ) THEN
            CALL GetRM (P, m, ErrStatLcl, ErrMessLcl, &
                        rLocal, DFN, DFT, AzimuthAngle, IElement, IBlade)
               CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput' )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL CleanUp()
                  RETURN
               END IF
            
         ENDIF
         
         IF (p%UseDWM) THEN
             m%DWM%Nforce(IElement,IBlade) = DFN              ! 12.4.2014 add by yh
         END IF ! UseDWM

         m%StoredForces(1,IElement,IBlade)  = ( DFN*CPitch + DFT*SPitch ) / p%Blade%DR(IElement)
         m%StoredForces(2,IElement,IBlade)  = ( DFN*SPitch - DFT*CPitch ) / p%Blade%DR(IElement)
         m%StoredForces(3,IElement,IBlade)  = 0.0

         m%StoredMoments(1,IElement,IBlade)  = 0.0
         m%StoredMoments(2,IElement,IBlade)  = 0.0
         m%StoredMoments(3,IElement,IBlade)  = PMA / p%Blade%DR(IElement)

!      DO IBlade=1,p%NumBl
!       DO IElement=1,p%Element%Nelm
!         y%OutputLoads(IBlade)%Force(:,IElement)  = m%StoredForces(:,IElement,IBlade)
!         y%OutputLoads(IBlade)%Moment(:,IElement) = m%StoredMoments(:,IElement,IBlade)
!       ENDDO
!!      ENDDO

            ! save velocities for output, if requested

         IF ( m%ElOut%WndElPrList(IElement) > 0 ) THEN
            m%ElOut%SaveVX( m%ElOut%WndElPrList(IElement), IBlade ) = VelocityVec(1)
            m%ElOut%SaveVY( m%ElOut%WndElPrList(IElement), IBlade ) = VelocityVec(2)
            m%ElOut%SaveVZ( m%ElOut%WndElPrList(IElement), IBlade ) = VelocityVec(3)
         ENDIF


      END DO !IElement

      IF ( IBlade == 1 .AND. p%ElemPrn ) THEN
         m%ElOut%VXSAV  = VelocityVec(1)
         m%ElOut%VYSAV  = VelocityVec(2)
         m%ElOut%VZSAV  = VelocityVec(3)
      ENDIF


   END DO !IBlade
  END DO  !KS -- end of LoopNum loop

   m%NoLoadsCalculated = .FALSE.


   DO IBlade=1,p%NumBl
     DO IElement=1,p%Element%Nelm
       y%OutputLoads(IBlade)%Force(:,IElement)  = m%StoredForces(:,IElement,IBlade)
       y%OutputLoads(IBlade)%Moment(:,IElement) = m%StoredMoments(:,IElement,IBlade)
     ENDDO
   ENDDO
   
   
   !------------------------DWM PHASE 3-----------------------------------------------
   IF (p%UseDWM) THEN
   
      IF (Time > 50.00 ) THEN !BJJ: why is 50 hard-coded here and above???
            
         !m%DWM%Nforce(:,:)    = m%DWM%DFN_DWM(:,:) 
         CALL   DWM_phase3( Time, m%DWM_Inputs, p%DWM, x%DWM, xd%DWM, z%DWM, O%DWM, m%DWM_Outputs, m%DWM, ErrStatLcl, ErrMessLcl )
    
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/DWM_phase3' )
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF   
      
         !CALL filter_average_induction_factor( AD14_ParameterType, DWM_ConstraintStateType, m%DWM_Outputs )
      END IF
   END IF !UseDWM
   
   !-----------------------------------------------------------------------------------


   
       
         

      ! Loop through all the tower nodes to calculate the aerodynamic loads on the tower if aerodynamics were requested.

   IF ( p%TwrProps%CalcTwrAero )  THEN

      DO Node=1,u%Twr_InputMarkers%Nnodes


            ! Calculate the aerodynamic load on this tower node: TwrAeroLoads ( p, Node, NodeDCMGbl, NodeVelGbl, NodeWindVelGbl, NodeFrcGbl )

         CALL TwrAeroLoads ( p, Node, u%Twr_InputMarkers%Orientation(:,:,Node), u%Twr_InputMarkers%TranslationVel(:,Node) &
                           , u%InflowVelocity(:,Node+p%NumBl*p%Element%NElm), y%Twr_OutputLoads%Force(:,Node) )

      END DO ! Node

   END IF ! ( p%TwrProps%CalcTwrAero )

   !................................................................................................
      
   
   CALL ElemOut(time, P, m )
   
   CALL CleanUp (  )

   RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
   SUBROUTINE CleanUp ( )


      ! This subroutine cleans up the parent routine before exiting.


      !   ! Deallocate the IfW_Inputs%Position array if it had been allocated.
      !
      !CALL IfW_DestroyInput( IfW_Inputs, ErrStatLcl, ErrMessLcl )


      RETURN

   END SUBROUTINE CleanUp 
   

END SUBROUTINE AD14_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE AeroDyn14
!**********************************************************************************************************************************

