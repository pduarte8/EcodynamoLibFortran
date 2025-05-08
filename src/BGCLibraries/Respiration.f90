module  respiration
  USE TemperatureLimitation      
  implicit none

contains

  function Respiration1(Maintenance, GrossProduction, RespirationCoefficient, &
                         WaterTemperature, TemperatureAugmentationRate, Tmin)
    real(8), intent(in) :: Maintenance, GrossProduction, RespirationCoefficient
    real(8), intent(in) :: WaterTemperature, TemperatureAugmentationRate, Tmin
    real(8) :: Respiration1, MyMaintenance, MyGrossProduction

    MyMaintenance = max(0.0d0, Maintenance)
    MyGrossProduction = max(0.0d0, GrossProduction)

    Respiration1 = MyMaintenance
    Respiration1 = Respiration1 + RespirationCoefficient * MyGrossProduction & 
                  * TemperatureExponentialLimitation(WaterTemperature, TemperatureAugmentationRate, Tmin)
  end function Respiration1

end module respiration
