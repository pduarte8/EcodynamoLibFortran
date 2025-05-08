MODULE Grazing
      IMPLICIT NONE
      CONTAINS



FUNCTION Hollings(K, Prey, HollingsType) RESULT(GrazingLimitation)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: K, Prey
  INTEGER, INTENT(IN) :: HollingsType
  REAL(8) :: GrazingLimitation, MyK, MyPrey
  INTEGER :: MyHollingsType
  REAL(8), PARAMETER :: TINNY = 1.0D-10

  MyK = MAX(0.0D0, K)
  MyPrey = MAX(0.0D0, Prey)
  MyHollingsType = MAX(0, HollingsType)

  IF (MyPrey > TINNY) THEN
    GrazingLimitation = MAX(0.0D0, (MyPrey**MyHollingsType) / (MyK**MyHollingsType + MyPrey**MyHollingsType))
  ELSE
    GrazingLimitation = 0.0D0
  END IF
END FUNCTION Hollings

FUNCTION IvlevFunction(Lambda, Prey) RESULT(GrazingLimitation)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: Lambda, Prey
  REAL(8) :: GrazingLimitation, MyLambda, MyPrey

  MyLambda = MAX(0.0D0, Lambda)
  MyPrey = MAX(0.0D0, Prey)

  GrazingLimitation = 1.0D0 - EXP(-MyLambda * MyPrey)
END FUNCTION IvlevFunction

END MODULE Grazing
