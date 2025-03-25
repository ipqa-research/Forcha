module critical_parameters
   !! This module is responsible for calculating the critical
   !! parameters of a fluid, namely:
   !!   - Critical Temperature (Tc)
   !!   - Critical Pressure (Pc)
   !!   - Acentric Factor (omega)
   !!
   !! The calculations are based on Pedersen correlations
   !! (adapted from the code "C7plusPedersenTcPcOm")
   !! and on the evaluation of EOS parameters taken from CubicParam.

   use constants ! Contains physical and mathematical constants
   use dtypes  !Defines data types, including FluidDataOut
   use defined_critical_parameters
   ! Defines critical parameters and related constants

contains

   subroutine get_critical_constants(characterization, eos)
      !! This subroutine calculates the critical parameters of a fluid, namely:
      !!   - Critical Temperature (Tc)
      !!   - Critical Pressure (Pc)
      !!   - Acentric Factor (omega)
      !!
      !! Additionally, it computes a modification function (m_funtion)
      !! used in the correlations.
      !!
      !! The equations used in the calculations are:
      !!   tc = c1 * rho + c2 * log(mw) + c3 * mw + (c4 / mw)
      !!   pc = exp(d1_p + d2 * (rho ** d5) + (d3 / mw) + (d4 / (mw ** 2)))
      !!   om = (-b + sqrt(b ** 2 - 4 * a * modified_c)) / (2 * a)
      !! where:
      !!   modified_c = c0 - m_funtion
      !!
      !! Depending on the EOS specified (e.g., "SRK", "PR", or "RKPR"),
      !! an appropriate expression for m_funtion is selected.
      !!
      !! Input:
      !!   characterization :
      !! FluidDataOut type (in/out) containing fluid information.
      !! It will be updated with the calculated critical parameters.
      !!
      !!   eos              :
      !! Character string specifying the equation of state to be used.
      !!
      !! Internal variables:
      !!   Allocatable arrays for storing intermediate values for Tc, Pc, omega,
      !!   and m_funtion.

      implicit none
      type(FluidDataOut), intent(inout) :: characterization
      character(len=*), intent(in) :: eos
      ! Allocatable arrays for intermediate calculations
      real(pr), allocatable :: modified_c(:)
      !! Array for the modified constant: modified_c = c0 - m_funtion
      real(pr), allocatable :: tc(:) !! Array for the critical temperature (Tc)
      real(pr), allocatable :: pc(:) !! Array for the critical pressure (Pc)
      real(pr), allocatable :: om(:) !! Array for the acentric factor (omega)
      real(pr), allocatable :: m_funtion(:)
      !! Array for the m function used in the correlations

      associate(&
         rho => characterization%lumped_densities, &
      ! rho: Lumped densities of the fluid
         mw  => characterization%lumped_mw &
      ! mw: Lumped molecular weights of the fluid
         )
         ! Allocate memory for the intermediate arrays
         allocate(modified_c(0))
         allocate(tc(0))
         allocate(pc(0))
         allocate(om(0))
         allocate(m_funtion(0))

         ! Call the subroutine that sets parameters for critical calculations
         ! based on the EOS
         call get_parameteres_for_critical(eos)

         ! Calculate the critical temperature (Tc) using the Pedersen correlation:
         tc = c1 * rho + c2 * log(mw) + c3 * mw + (c4 / mw)

         ! Calculate the critical pressure (Pc) using an exponential correlation:
         pc = exp(d1_p + d2 * (rho ** d5) + (d3 / mw) + (d4 / (mw ** 2)))

         !Determine the modification function (m_funtion) based on the specified EOS
         if (eos == "SRK") then
            ! For the Soave-Redlich-Kwong (SRK) EOS:
            m_funtion = e1 + e2 * mw + e3 * rho + e4 * (mw ** 2)
         else
            if (eos == "PR" .or. eos =="RKPR") then
               ! For the Peng-Robinson (PR) or RKPR EOS:
               m_funtion = e3 * rho + 0.2 + 1.8 * (1 - exp(-mw / 220))
               ! Alternative expression without max with MW - May 2019
            endif
         endif

         ! Calculate the modified constant:
         modified_c = c0 - m_funtion
         ! Calculate the acentric factor (omega) using the quadratic formula:
         om = (-b + sqrt(b ** 2 - 4 * a * modified_c))/(2 * a)
         ! Update the characterization structure with the calculated
         ! critical parameters.
         characterization%critical_temperature = tc
         characterization%critical_pressure = pc
         characterization%acentric_factor = om
         characterization%m_funtion = m_funtion

      end associate

   end subroutine get_critical_constants

end module critical_parameters





