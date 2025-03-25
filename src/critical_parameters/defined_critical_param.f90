module defined_critical_parameters
   !! This module defines global arrays and coefficients used for calculating 
   !! the critical properties of a fluid. These properties include the critical
   !! temperature (Tc), critical pressure (Pc), and acentric factor (omega). 
   !! The module provides default values for these properties as well as 
   !! correlation coefficients used in Pedersen correlations and
   !! EOS parameter calculations.
   use constants
   
   implicit none

   ! Allocatable arrays for default critical properties for various components:
   real(pr), allocatable :: tc_def(:)  !! Define critical temperatures [K]
   real(pr), allocatable :: pc_def(:)  !! Define critical pressures [bar]
   real(pr), allocatable :: om_def(:)  !! Define acentric factors

   ! Global coefficients for critical property correlations:
   ! Coefficients for critical temperature (Tc) correlation:
   real(pr) :: c1, c2, c3, c4
   ! Coefficients for critical pressure (Pc) correlation:
   real(pr) :: d1_p, d2, d3, d4, d5
   ! Coefficients for the modification function (m_funtion) used in 
   ! the correlations:
   real(pr) :: e1, e2, e3, e4
   ! Additional coefficients for acentric factor (omega) calculation via a 
   ! quadratic expression:
   real(pr) :: a, b, c0, del_1  
   ! 'del_1' is an additional parameter (currently not used).

contains

   subroutine get_parameteres_for_critical(eos)
      !! This subroutine initializes the define componets critical property 
      !! arrays (tc_def, pc_def,om_def) and sets the global correlation 
      !! coefficients based on the specified equation of state (EOS).
      !!
      !! The default critical arrays are defined element-wise for different 
      !! components.
      !!
      !! Depending on the input string 'eos', the subroutine assigns 
      !! different sets of coefficients:
      !!    - "SRK"  : Soave-Redlich-Kwong EOS coefficients.
      !!    - "PR"   : Peng-Robinson EOS coefficients (original Pedersen values).
      !!    - "RKPR" : Modified Peng-Robinson (RKPR) EOS coefficients.
      !!
      !! Input:
      !!   eos : A character string specifying the equation of state to be used.
      !!         Supported values: "SRK", "PR", "RKPR".
      !!
      !! The coefficients set in this subroutine are used by other modules to
      !! calculate the fluid's critical temperature (Tc), pressure (Pc), and 
      !! acentric factor (omega).
      
      implicit none
      character(len=*) :: eos

      ! Define default critical properties for various components:
      tc_def = [126.2, 304.21, 190.564, 305.32, 369.83, 408.14, 425.12, 460.43, 469.7]
      pc_def = [34.0, 73.83, 45.99, 48.72, 42.48, 36.48, 37.96, 33.81, 33.7]
      om_def = [0.038, 0.224, 0.012, 0.099, 0.152, 0.181, 0.20, 0.228, 0.252]

      ! Set correlation coefficients based on the specified EOS.
      select case(eos)

       case("SRK")
         ! For the Soave-Redlich-Kwong (SRK) EOS:
         c1 = 1.6312d2
         c2 = 8.6052d1
         c3 = 4.3475d-1
         c4 = -1.8774d3
         d1_p = -1.3408d-1
         d2 = 2.5019
         d3 = 2.0846d2
         d4 = -3.9872d3
         d5 = 1.0d0
         e1 = 7.4310d-1
         e2 = 4.8122d-3
         e3 = 9.6707d-3
         e4 = -3.7184d-6
         a = -0.176
         b = 1.574
         c0 = 0.48
         !del_1 = 1.0D0

       case("PR") ! Coeficientes originales Tabla 5.3 Pedersen
         c1 = 7.34043d1
         c2 = 9.73562d1
         c3 = 6.18744d-1
         c4 = -2.05932d3
         d1_p = 7.28462d-2
         d2 = 2.18811d0
         d3 = 1.63910d2
         d4 = -4.04323d3
         d5 = 0.25d0
         e1 = 3.73765d-1
         e2 = 5.49269d-3
         e3 = 1.17934d-2
         e4 = -4.93049d-6
         a = -0.26992
         b = 1.54226
         c0 = 0.37464
         !del_1 = 1.0D0 + sqrt(2.0)

       case("RKPR")
         c1 = 7.34043d1  ! Coeficientes originales Tabla 5.3 Pedersen
         c2 = 9.73562d1
         c3 = 6.18744d-1
         c4 = -2.05932d3
         d1_p = 2.3d0
         d2 = -0.5d0
         d3 = 1.85d2
         d4 = -4.0d3
         d5 = 0.25d0
         e1 = 3.73765d-1
         e2 = 5.49269d-3
         e3 = 1.17934d-2
         e4 = -4.93049d-6
         a = -0.26992
         b = 1.54226
         c0 = 0.37464

      end select

   end subroutine get_parameteres_for_critical

end module defined_critical_parameters