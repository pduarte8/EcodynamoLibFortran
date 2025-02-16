      MODULE BiogeochemicalProcesses
      IMPLICIT NONE
      CONTAINS

      FUNCTION DenitrificationToNH4(NO3, Kdenit, TemperatureLimitation,&
      &		                    OxygenLimitation) RESULT(DenitToNH4)
         REAL(8), INTENT(IN) :: NO3, Kdenit, TemperatureLimitation, &
	                        OxygenLimitation
         REAL(8) :: DenitToNH4, MyNO3, MyKdenit, MyTempLim, MyOxyLim
    
         MyNO3 = MAX(0.0D0, NO3)
         MyKdenit = MAX(0.0D0, Kdenit)
         MyTempLim = MAX(0.0D0, TemperatureLimitation)
         MyOxyLim = MAX(0.0D0, MIN(1.0D0, OxygenLimitation))
    
         DenitToNH4 = MyKdenit * MyNO3 * MyTempLim * (1.0D0 - MyOxyLim)
      END FUNCTION DenitrificationToNH4

      FUNCTION DenitrificationToN2(DenitrificationToNH4, alfa) &
      &                             RESULT(DenitToN2)
         REAL(8), INTENT(IN) :: DenitrificationToNH4, alfa
         REAL(8) :: DenitToN2, MyDenitToNH4, MyAlfa
    
         MyDenitToNH4 = MAX(0.0D0, DenitrificationToNH4)
         MyAlfa = MAX(0.0D0, MIN(1.0D0, alfa))
    
         DenitToN2 = MyDenitToNH4 * MyAlfa
      END FUNCTION DenitrificationToN2

      FUNCTION Mineralization(minR, TemperatureLimitation, &
      &		           OxygenLimitation, Xorganic) RESULT(MinRate)
         REAL(8), INTENT(IN) :: minR, TemperatureLimitation, OxygenLimitation, Xorganic
         REAL(8) :: MinRate, MyminR, MyTempLim, MyOxyLim, MyXorganic
    
         MyminR = MAX(0.0D0, minR)
         MyTempLim = MAX(0.0D0, TemperatureLimitation)
         MyOxyLim = MAX(0.0D0, MIN(1.0D0, OxygenLimitation))
         MyXorganic = MAX(0.0D0, Xorganic)
    
         MinRate = MyminR * MyXorganic * MyTempLim * MyOxyLim
      END FUNCTION Mineralization

      FUNCTION Nitrification(NH4, Knit, TemperatureLimitation, &
      &		             OxygenLimitation, LightLimitation) &
      &                       RESULT(NitRate)
         REAL(8), INTENT(IN) :: NH4, Knit, TemperatureLimitation, &
      &	                        OxygenLimitation, LightLimitation
         REAL(8) :: NitRate, MyNH4, MyKnit, MyTempLim, MyOxyLim,&
      &	            MyLightLim
    
         MyNH4 = MAX(0.0D0, NH4)
         MyKnit = MAX(0.0D0, Knit)
         MyTempLim = MAX(0.0D0, TemperatureLimitation)
         MyOxyLim = MAX(0.0D0, MIN(1.0D0, OxygenLimitation))
         MyLightLim = MAX(0.0D0, MIN(1.0D0, LightLimitation))
    
         NitRate = MyKnit * MyNH4 * MyTempLim * MyOxyLim * MyLightLim
      END FUNCTION Nitrification

      FUNCTION OrganicDissolution(dissR, TemperatureLimitation, &
      &		                  Xorganic) RESULT(DissRate)
         REAL(8), INTENT(IN) :: dissR, TemperatureLimitation, Xorganic
         REAL(8) :: DissRate, MydissR, MyTempLim, MyXorganic
    
         MydissR = MAX(0.0D0, dissR)
         MyTempLim = MAX(0.0D0, TemperatureLimitation)
         MyXorganic = MAX(0.0D0, Xorganic)
    
         DissRate = MydissR * MyXorganic * MyTempLim
      END FUNCTION OrganicDissolution

      FUNCTION PhosphorusAdsorption(PO4InPoreWater, Pads, Pmax, &
      &		      OxygenInPoreWater, OxygenThreshold, Ka1, Ka2) &
      &               RESULT(Adsorption)
         REAL(8), INTENT(IN) :: PO4InPoreWater, Pads, Pmax, &
      &	                  OxygenInPoreWater, OxygenThreshold, Ka1, Ka2
         REAL(8) :: Adsorption, MyPmax, MyPO4, MyPads, MyOxygen, &
      &	            MyOxyThr, MyKa1, MyKa2
         REAL(8), PARAMETER :: TINNY = 1.0D-10
    
         MyPmax = MAX(0.0D0, Pmax)
         MyPads = MIN(MAX(0.0D0, Pads), MyPmax)
         MyPO4 = MAX(0.0D0, PO4InPoreWater)
         MyOxygen = MAX(0.0D0, OxygenInPoreWater)
         MyOxyThr = MAX(0.0D0, OxygenThreshold)
         MyKa1 = MAX(0.0D0, Ka1)
         MyKa2 = MAX(0.0D0, Ka2)
    
         IF (MyPmax > TINNY) THEN
            IF (MyOxygen > MyOxyThr) THEN
               Adsorption = MyKa1 * (1.0D0 - MyPads / MyPmax) * MyPO4
            ELSE
               Adsorption = MyKa2 * (1.0D0 - MyPads / MyPmax) * MyPO4
            END IF
         ELSE
            Adsorption = 0.0D0
         END IF
      END FUNCTION PhosphorusAdsorption

      FUNCTION PhosphorusDesorption(Pads, Kd, Pmax) RESULT(Desorption)
         REAL(8), INTENT(IN) :: Pads, Kd, Pmax
         REAL(8) :: Desorption, MyPads, MyKd, MyPmax
         REAL(8), PARAMETER :: TINNY = 1.0D-10
    
         MyKd = MAX(0.0D0, Kd)
         MyPmax = MAX(0.0D0, Pmax)
         MyPads = MIN(MAX(0.0D0, Pads), MyPmax)
    
         IF (MyPmax > TINNY) THEN
            Desorption = MyKd * MyPads / MyPmax
         ELSE
            Desorption = 0.0D0
         END IF
      END FUNCTION PhosphorusDesorption

      END MODULE BiogeochemicalProcesses
