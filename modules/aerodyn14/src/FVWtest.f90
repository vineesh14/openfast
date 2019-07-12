SUBROUTINE FVWtest(J, IBlade, Initial, p_FVW, W, CLFW, VINDFW, Time )

  USE AeroDyn14_Types, Only: FVW_ParameterType
  USE FVW_Parm, Only: Time_Real!, RotSpeed   !KS -- removed 6.28.19
  USE FVW_vars
  USE NWTC_Library

  IMPLICIT NONE

  TYPE( FVW_ParameterType ), INTENT( INOUT ) :: p_FVW

  REAL(ReKi),                INTENT(   OUT ) :: W, CLFW, VINDFW(3)
  REAL(ReKi),                INTENT( IN    ) :: Time

  LOGICAL,                   INTENT( IN    ) :: Initial

  INTEGER,                   INTENT( IN    ) :: J
  INTEGER,                   INTENT( IN    ) :: IBlade

  INTEGER :: InitVal, a

  Time_Real = Time
  PRINT*, 'RotSpeed in FVWtest: ', p_FVW%RotSpeed
  IF ( Initial .AND. J .EQ. 1 .AND. IBLADE .EQ. 1 ) THEN
     InitVal=1
     CALL FVW_READ_WAKE_PARAM( p_FVW )
     CALL FVW_INITIALIZE_WAKE(  )
    PRINT*, 'a'

        CALL FVW_COMPUTE_WAKE( p_FVW%FVWTurbineComponents, p_FVW%FVWInputMarkers, p_FVW%FVW_Wind )
    PRINT*, 'b'
  ELSE
   PRINT*, 'c'
     InitVal=0

     IF (J .EQ. 1 .AND. IBLADE .EQ. 1) THEN
        AofA=0.0_ReKi
        W2FVW=0.0_ReKi
        CALL FVW_COMPUTE_WAKE( p_FVW%FVWTurbineComponents, p_FVW%FVWInputMarkers, p_FVW%FVW_Wind )
     ENDIF
  ENDIF
  W=W2FVW(J,IBLADE)
  CLFW=CLFVW(J,Iblade)
  VINDFW=VINDFVW(:, J, IBlade)

END SUBROUTINE FVWtest
