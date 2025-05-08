module LightFunctions
  implicit none
  contains

  function Platt1(PARtop, KValue, Depth, Pmax, beta, slope, EulerSteps) result(LightLimitation)
    implicit none
    real(8), intent(in) :: PARtop, KValue, Depth, Pmax, beta, slope
    integer, intent(in) :: EulerSteps
    real(8) :: DeltaZ, Soma, TINNY, PAR, MyKValue, MyDepth, MyPmax, Mybeta, Myslope
    integer :: MySteps, Step
    real(8) :: LightLimitation
    
    TINNY = 1.0d-10
    Soma = 0.0d0
    LightLimitation = 0.0d0
    
    PAR = max(0.0d0, PARtop)
    MyKValue = max(0.0d0, KValue)
    MyDepth = max(0.0d0, Depth)
    MyPmax = max(0.0d0, Pmax)
    Mybeta = max(0.0d0, beta)
    Myslope = max(0.0d0, slope)
    MySteps = max(1, EulerSteps)
    
    if ((MyDepth > TINNY) .and. (PAR > TINNY) .and. (MyPmax > TINNY)) then
      DeltaZ = MyDepth / MySteps
      do Step = 1, MySteps
        Soma = Soma + (1.0d0 - exp(-Myslope * PAR / MyPmax)) * exp(-Mybeta * PAR / MyPmax) * DeltaZ
        PAR = PAR * exp(-MyKValue * DeltaZ)
      end do
      LightLimitation = Soma / MyDepth
    end if
  end function Platt1

  function Platt2(PAR, Pmax, beta, slope) result(LightLimitation)
    implicit none
    real(8), intent(in) :: PAR, Pmax, beta, slope
    real(8) :: LightLimitation, TINNY, MyPAR, MyPmax, Mybeta, Myslope
    
    TINNY = 1.0d-10
    MyPAR = max(0.0d0, PAR)
    MyPmax = max(0.0d0, Pmax)
    Mybeta = max(0.0d0, beta)
    Myslope = max(0.0d0, slope)
    
    if ((MyPAR > TINNY) .and. (MyPmax > TINNY)) then
      LightLimitation = (1.0d0 - exp(-Myslope * MyPAR / MyPmax)) * exp(-Mybeta * MyPAR / MyPmax)
    else
      LightLimitation = 0.0d0
    end if
  end function Platt2

  function Steele1(PARtop, KValue, Depth, PARopt) result(LightLimitation)
    implicit none
    real(8), intent(in) :: PARtop, KValue, Depth, PARopt
    real(8) :: LightLimitation, TINNY, MyPARtop, MyKValue, MyDepth, MyPARopt, MyPARbottom
    
    TINNY = 1.0d-10
    MyPARtop = max(0.0d0, PARtop)
    MyKValue = max(0.0d0, KValue)
    MyDepth = max(0.0d0, Depth)
    MyPARopt = max(0.0d0, PARopt)
    MyPARbottom = max(0.0d0, MyPARtop * exp(-MyKValue * MyDepth))
    
    if ((MyKValue * MyDepth > TINNY) .and. (MyPARtop > TINNY)) then
      LightLimitation = exp(1.0d0) / (MyKValue * MyDepth) * (exp(-MyPARbottom / MyPARopt) - exp(-MyPARtop / MyPARopt))
    else
      LightLimitation = 0.0d0
    end if
  end function Steele1

  function Steele2(PAR, PARopt) result(LightLimitation)
    implicit none
    real(8), intent(in) :: PAR, PARopt
    real(8) :: LightLimitation, TINNY, MyPAR, MyPARopt
    
    TINNY = 1.0d-10
    MyPAR = max(0.0d0, PAR)
    MyPARopt = max(0.0d0, PARopt)
    
    if ((MyPAR > TINNY) .and. (MyPARopt > TINNY)) then
      LightLimitation = MyPAR / MyPARopt * exp(1.0d0 - MyPAR / MyPARopt)
    else
      LightLimitation = 0.0d0
    end if
  end function Steele2


  function SteeleSlope(Pmax, Iopt) result (Slope)
    implicit none    
    real(8), intent(in) :: Pmax, Iopt
    real(8) ::  Slope, TINNY, MyPmax, MyIopt
   
    TINNY = 1.0d-10
    MyPmax = max(0.0d0,Pmax)
    MyIopt = max(0.0d0,Iopt)
    if (MyIopt > TINNY) then 
      Slope= MyPmax * exp(1.0) / MyIopt
    else
      Slope = 0.0d0
    end if  
  end function SteeleSlope

!Eilers and Peeters P-I function
!Eilers, P. H. C., Peeters, J. C. H.: A model for the relationship between light intensity and the rate of photosynthesis in phytoplankton,
!Ecological Modelling, 42, 199-215, 1988
!EilersAndPeeters1 returns vertically averaged light inhibition. This is achieved with the analytical integration of the EilersAndPeeters function, which is shown in EilersAndPeeters2.
!EilersAndPeeters2 returns light inhibition at a specified PAR level.

  function EilersAndPeeters1(PARtop, KValue, Depth, a, b, c, Pmax) result (LightLimitation)
    implicit none
    real(8), intent(in) :: PARtop, KValue, Depth, a, b, c, Pmax
    real(8) :: D, B1, B2, LightLimitation, TINNY, P, MyPARtop, MyKValue, MyDepth, MyPARbottom, Mya, Myc, MyPmax

    LightLimitation = 0.0d0
    TINNY = 1.0d-10
    MyPARtop = max(0.0d0,PARtop)
    MyKValue = max(0.0d0,KValue)
    MyDepth =  max(0.0d0,Depth)
    Mya = max(0.0d0,a)
    Myc = max(0.0d0,c)
    MyPmax = max(0.0d0,Pmax)
    MyPARbottom = max(0.0d0, MyPARtop*exp(-MyKValue * MyDepth))
    if (MyPARtop - MyPARbottom > TINNY) then
      if ((MyDepth > TINNY) .and. (MyPARtop > TINNY) .and. (MyPmax > TINNY)) then
         if (Mya /= 0.0) then
            D = b * b - 4.0 * Mya * Myc
            B1 = 2.0 * Mya * MyPARtop + b
            B2 = 2.0 * Mya * MyPARbottom + b
            if (D < 0.0) then
               P = 2.0 / (MyKValue*sqrt(-D)) *(atan(B1/sqrt(-D))- atan(B2/sqrt(-D))) / MyDepth
            else if (D == 0.0) then
               P = 2.0 / MyKValue * (1.0 / B2 - 1.0 / B1) / MyDepth
            else if (D > 0.0) then
               P = 1.0 / (MyKValue*sqrt(D)) * log (((B1 - sqrt(D))*(B2 + sqrt(D)))/ &
                         &  ((B1 + sqrt(D)) * (B2 - sqrt(D)))) / MyDepth
            end if
         else
            P = 1.0 / (MyKValue * b) * log (ABS(b * MyPARtop + Myc)/ &
                       &  ABS(b * MyPARbottom + Myc)) / MyDepth
         end if
         LightLimitation = P/MyPmax
      end if    
    else
      LightLimitation = EilersAndPeeters2(MyPARtop, Mya, b, Myc, MyPmax)!If top and bottom light do not differ significantly, use the standard function, without vertical averaging
    end if    
  end function EilersAndPeeters1

  function EilersAndPeeters2(PAR, a, b, c, Pmax) result (LightLimitation)
    implicit none
    real(8), intent(in) :: PAR, a, b, c, Pmax
    real(8) :: LightLimitation, TINNY, MyPAR, x, P, Mya, Myc, MyPmax

    LightLimitation = 0.0d0
    TINNY = 1.0d-10
    MyPAR = max(0.0d0,PAR)
    Mya = max(0.0d0,a)
    Myc = max(0.0d0,c)
    MyPmax = max(0.0d0,Pmax)
    x = Mya * MyPAR * MyPAR + b * MyPAR + Myc
    if ((x > TINNY) .and. (MyPmax > TINNY)) then
      P = MyPAR / x
      LightLimitation = P / MyPmax
    end if
  end function EilersAndPeeters2

  function EilersAndPeetersSlope(C) result (Slope) 
   implicit none
   real(8), intent(in) :: C
   real(8) :: Slope, TINNY, MyC

   TINNY = 1.0d-10   
   MyC = max(0.0d0,C)
   if (MyC > TINNY) then 
     Slope = 1.0 / MyC
   else
     Slope = 0.0d0
   end if  
  end function EilersAndPeetersSlope


  function LightLimNitr(KI, I0, Light) result(Result)
    implicit none
    real(8), intent(in) :: KI, I0, Light
    real(8) :: MyKI, MyI0, MyLight, Result
    
    MyLight = max(0.0d0, Light)
    MyKI = max(0.0d0, KI)
    MyI0 = max(0.0d0, I0)
    
    if (MyLight > MyI0) then
      Result = 1.0d0 - max(0.0d0, (MyLight - MyI0) / (MyKI + MyLight - MyI0))
    else
      Result = 1.0d0
    end if
  end function LightLimNitr

  function VertLightAvg(LightAtTop, KValue, Depth) result(Result)
    implicit none
    real(8), intent(in) :: LightAtTop, KValue, Depth
    real(8) :: TINNY, MyLight, MyKValue, MyDepth, Result
    
    TINNY = 1.0d-10
    MyLight = max(0.0d0, LightAtTop)
    MyKValue = max(0.0d0, KValue)
    MyDepth = max(0.0d0, Depth)
    
    if (MyKValue * MyDepth > TINNY) then
      Result = MyLight * (1.0d0 - exp(-MyKValue * MyDepth)) / (MyKValue * MyDepth)
    else
      Result = MyLight
    end if
  end function VertLightAvg

end module LightFunctions
