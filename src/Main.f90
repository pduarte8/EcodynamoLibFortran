        program main
        USE BiogeochemicalProcesses
        IMPLICIT NONE
        REAL(8) :: result

        result = DenitrificationToNH4(1.0D0, 0.2D0, 1.07143620729633D0, 0.909090909090909D0)
        PRINT *, 'Denitrification to NH4:', result
        end program Main
