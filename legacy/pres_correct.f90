      SUBROUTINE pres_correct

! Corrects pressures (SIMPLE)

      USE control_parameters
      USE flow
      IMPLICIT NONE

      pres = pres + relax_p*pp

      END SUBROUTINE pres_correct
