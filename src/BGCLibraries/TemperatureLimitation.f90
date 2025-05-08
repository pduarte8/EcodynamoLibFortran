module TemperatureLimitation
  implicit none

contains

  function TemperatureExponentialLimitation(Temperature, TemperatureAugmentationRate, Tmin) result(TemperatureLim)
    real(8), intent(in) :: Temperature, TemperatureAugmentationRate, Tmin
    real(8) :: TemperatureLim

    TemperatureLim = exp(TemperatureAugmentationRate * (Temperature - Tmin))

  end function TemperatureExponentialLimitation

end module TemperatureLimitation
