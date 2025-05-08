        program main

        USE BiogeochemicalProcesses
        USE Grazing
        USE LightFunctions
        USE NutrientLimitation
        USE Respiration
        USE TemperatureLimitation

        IMPLICIT NONE
        REAL(8) :: result


        character(30) :: fname
        INTEGER :: EulerSteps,HollingsType
        DOUBLE PRECISION :: PARtop,KValue,Depth
        DOUBLE PRECISION :: PARbottom, PARopt,PlattPARopt
        DOUBLE PRECISION :: Pmax,beta,slope
        DOUBLE PRECISION :: a, b, c
        DOUBLE PRECISION :: LightLim,TempLim
        DOUBLE PRECISION :: PAR
        DOUBLE PRECISION :: Temp, TempAugRate, Tmin
        DOUBLE PRECISION :: NutLim, CellQuota
        DOUBLE PRECISION :: MinCellQuota, HalfSaturation
        DOUBLE PRECISION :: NO3, KNO3,NH4, KNH4,KI,I0
        DOUBLE PRECISION :: NO2, N, KN
        DOUBLE PRECISION :: GrazingLim,K,Prey,Lambda,alfa
        DOUBLE PRECISION :: Kdenit,O2,O2thr,alfa1,Knit
        DOUBLE PRECISION :: PO4, Pads, PmaxSed,OxyPore
        DOUBLE PRECISION :: OxyPoreThr, Ka1, Ka2, Kd,KO2
        DOUBLE PRECISION :: DenitToNH4, DenitToN2
        DOUBLE PRECISION :: Nit, Pdes
        DOUBLE PRECISION :: OxyLim
        DOUBLE PRECISION :: Maintenance,GrossProduction
        DOUBLE PRECISION :: RespirationCoeff,Resp

        PARtop = 5.13509
        PAR = PARtop
        KValue = 0.0703777
        Depth = 2.32157
        Pmax = 0.0476
        beta = 3.64e-05
        slope = 0.000456522
        EulerSteps = 10
        PARopt = 5.13509e-06
        PlattPARopt = (Pmax*log((beta+slope)/beta))/slope
        a = 1.0/(slope*PARopt*PARopt)
        b = 1.0/Pmax-2.0/(slope*PARopt)
        c = 1.0/slope
        write(*,*) 'a = ',a,' b = ',b,' c = ',c
        Temp = 1.0
        TempAugRate = 0.069
        Tmin = 0.0
        CellQuota = 0.18
        MinCellQuota = 0.05
        HalfSaturation = 0.09
        NO3 = 1.0
        NO2 = 0.1
        NH4 = 0.1
        N = NO3
        KNO3 = 2.0
        KNH4 = 0.2
        KN = KNO3
         K = 1.0
        Prey = 2.0
        HollingsType = 0
        Lambda = 1.0
        O2 = 5.0
        O2thr = 1.0
        KO2 = 0.5
        Kdenit = 0.2
        Knit = 0.5
        alfa = 0.2
        PO4 = 1.3
        Pads = 10.0
        PmaxSed = 20.0
        OxyPore = 1.0
        OxyPoreThr = 0.2
        Ka1 = 0.1
        Ka2 = 0.05
        Kd = 0.1
        KI = 0.01
        I0 = 0.0095
        Maintenance = 0.002
        RespirationCoeff = 0.184
        GrossProduction = 1.0

        write(*,*) 'PARtop =',PARtop

        LightLim =  Platt1(PARtop,KValue,Depth, &
        & Pmax,beta,slope,EulerSteps)
     
        write(*,*) 'Platt LightLim1 =',LightLim

        LightLim =  Platt2(PAR, &
        & Pmax,beta,slope)

        write(*,*) 'Platt LightLim2 =',LightLim
      ! Platt done

      ! Steele beginning
        LightLim =  Steele1(PARtop,             &
        &                  KValue,Depth,PARopt)

        write(*,*) 'Steele LightLim1 =',LightLim
         LightLim =  Steele2(PAR,PARopt)

        write(*,*) 'Steele LightLim2 =',LightLim
      ! Steele done

      ! Eilers & Peeters beginning

        LightLim =  EilersAndPeeters1(PARtop,   &
        &                 KValue,Depth,a,b,c,Pmax)

        write(*,*) 'EilersAndPeeters LightLim1 =',LightLim

        LightLim =  EilersAndPeeters2(PAR,   &
        &                 a,b,c,Pmax)

        write(*,*) 'EilersAndPeetersX LightLim2 =',LightLim


      ! Eiler & Peeters done
        
         ! Temperature limitation

        TempLim = TemperatureExponentialLimitation( &
        & Temp,TempAugRate, Tmin)
        write(*,*) 'Temperature limitation =',TempLim

      ! Temperature limitation done

      ! Internal nutrient limitation

        NutLim = InternalNutrientLimitation(CellQuota, &
        &  MinCellQuota, HalfSaturation)

        write(*,*) 'Inter nut limitation =',NutLim

      ! Internal nutrient limitation done

      ! External nutrient limitation
        NutLim = NitrateAndAmmoniumLimitation(NH4,KNH4,&
        & NO3, KNO3, NO2)

        write(*,*) 'NO3+NH4 limitation =',NutLim

      ! External nutrient limitation done

      ! Michaelis-Menten

        NutLim = MichaelisMentenLimitation(N,KN)

        write(*,*) 'Michaelis-Menten lim =',NutLim
      ! Michaelis-Menten done
        
        ! Hollings Type 0, 1, 2....

        write(*,*) 'Prey = ',Prey
        write(*,*) 'K =',K
        write(*,*) 'HollingsType =',HollingsType
        GrazingLim = Hollings(K,Prey,HollingsType)
        write(*,*) 'Grazing Hollings type 0 =',GrazingLim

        HollingsType = 1
        GrazingLim = Hollings(K,Prey,HollingsType)
        write(*,*) 'Grazing Hollings type 1 =',GrazingLim

        HollingsType = 2
        GrazingLim = Hollings(K,Prey,HollingsType)
        write(*,*) 'Grazing Hollings type 2 =',GrazingLim

        GrazingLim = IvlevFunction(Lambda,Prey)
        write(*,*) 'Grazing Ivlev =',GrazingLim

        OxyLim = MichaelisMentenLimitation(O2,KO2)
   
      ! Denitrification
        DenitToNH4 = DenitrificationToNH4 &
        & (NO3,Kdenit,TempLim,OxyLim)
        DenitToN2 = DenitrificationToN2(DenitToNH4,alfa)
        write(*,*) 'DenitToNH4 =',DenitToNH4
        write(*,*) 'DenitToN2 =',DenitToN2
      ! Denitrifcation done
      ! Nitrification
        LightLim = LightLimNitr(KI,I0,PAR);
        Nit = Nitrification(NH4,Knit,TempLim,OxyLim,LightLim)
        write(*,*) 'Nitrification =',Nit
      ! Nitrification done

      ! Respiration
        Resp = Respiration1(Maintenance,GrossProduction,&
                & RespirationCoeff,Temp,&
                & TempAugRate,Tmin)
        write(*,*) 'Respiration =',Resp
      ! Respiration done

        end program Main
