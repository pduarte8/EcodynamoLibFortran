module NutrientLimitation
  implicit none
  contains

  ! Function to compute Internal Nutrient Limitation
  function InternalNutrientLimitation(CellQuota, MinCellQuota, HalfSaturation) result(LiebigLimitation)
    real(8), intent(in) :: CellQuota, MinCellQuota, HalfSaturation
    real(8) :: LiebigLimitation, Limitation, MyCellQuota, MyMinCellQuota, MyHalfSaturation
    
    MyCellQuota = max(0.0d0, CellQuota)
    MyMinCellQuota = max(0.0d0, MinCellQuota)
    MyHalfSaturation = max(0.0d0, HalfSaturation)
    
    Limitation = max((MyCellQuota - MyMinCellQuota) / (MyCellQuota + MyHalfSaturation - MyMinCellQuota), 0.0d0)
    LiebigLimitation = min(1.0d0, Limitation)
  end function InternalNutrientLimitation

  ! Function to compute Nitrate and Ammonium Limitation
  function NitrateAndAmmoniumLimitation(NH4, KNH4, NO3, KNO3, NO2) result(Limitation)
    real(8), intent(in) :: NH4, KNH4, NO3, KNO3, NO2
    real(8) :: Limitation, L_NH4, L_NO3, cff1, cff2
    real(8) :: MyNH4, MyKNH4, MyNO3, MyKNO3, MyNO2
    real(8), parameter :: TINNY = 1.0d-10
    
    MyNH4 = max(0.0d0, NH4)
    MyKNH4 = max(0.0d0, KNH4)
    MyNO3 = max(0.0d0, NO3)
    MyKNO3 = max(0.0d0, KNO3)
    MyNO2 = max(0.0d0, NO2)
    
    if (MyKNH4 > TINNY) then
      cff1 = MyNH4 / MyKNH4
      L_NH4 = max(cff1 / (1.0d0 + cff1), 0.0d0)
    else
      L_NH4 = 1.0d0
    end if
    
    if (MyKNO3 > TINNY) then
      cff2 = (MyNO3 + MyNO2) / MyKNO3
      L_NO3 = max(cff2 / (1.0d0 + cff2), 0.0d0)
    else
      L_NO3 = 1.0d0
    end if
    
    Limitation = min(L_NO3 + L_NH4, 1.0d0)
  end function NitrateAndAmmoniumLimitation

  ! Function to compute Michaelis-Menten Limitation
  function MichaelisMentenLimitation(N, KN) result(Limitation)
    real(8), intent(in) :: N, KN
    real(8) :: Limitation
    real(8), parameter :: TINNY = 1.0d-10
    
    if (KN > TINNY) then
      Limitation = N / (N + KN)
    else
      Limitation = 1.0d0
    end if
  end function MichaelisMentenLimitation

end module NutrientLimitation