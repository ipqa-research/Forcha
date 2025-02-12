module density
   !! This module is responsible for calculating the density function for a 
   !! fluid. It contains routines to:
   !!   - Compute the experimental volume using density data (density_method).
   !!   - Iteratively adjust the density function parameters so that the 
   !! calculated volume matches the experimental volume (density_funtion).
   !!   - Calculate the volume based on the current density function parameters 
   !! (calculate_volume)
   
   use constants
   use dtypes, only: FluidData, FluidDataOut

contains

   subroutine density_funtion(fluid, rho_method, characterization)
      !!   This subroutine adjust the density function constants (a_d and b_d) 
      !! until the calculated volume matches the experimental volume. This is 
      !! done using an iterative (secant-like) approach.
      !!
      !! Process:
      !!   1. Call density_method to compute the experimental volume 
      !!      (volume_exp) based on the fluid data.
      !!   2. Initialize the density function constants a_d and b_d with an 
      !!      initial guess.
      !!   3. Call calculate_volume to compute the calculated volume 
      !!      (volume_cal) using the current density function parameters.
      !!   4. Compute the difference between volume_cal and volume_exp.
      !!   5. Iteratively adjust a_d (and update b_d accordingly) until the 
      !!      absolute difference is less than the specified tolerance (0.001).

      implicit none
      type(FluidData) :: fluid !! (FluidData) Contains fluid properties.
      type(FluidDataOut) :: characterization
      !! Structure to store density function constants and volumes.
      integer, intent(in) :: rho_method 
      !! Selects the density calculation method
      real(pr) :: a_d_old, aux
      real(pr) :: difference, difference_old


      associate( &
         a_d => characterization%a_d , &
         ! Slope constant of the density function
         b_d => characterization%b_d, &
         ! Intercept constant of the density function
         volume_cal => characterization%volume_cal, &
         ! Calculated volume from the density function.
         volume_exp => characterization%volume_exp &
         ! Experimental volume obtained from the fluid data
         )
         
         ! Compute the experimental volume using the selected density method.
         call density_method (fluid, rho_method,characterization)

         ! Now find constants for the Density function
         ! Initialize the density function constants with initial guesses
         a_d = -0.50_pr ! Initial guess for a_d.
         b_d = 0.685_pr - a_d * exp(-0.60_pr) ! Initial guess for b_d
         
         ! Calculate the volume using the current density function.
         call calculate_volume(fluid, rho_method, characterization)
         difference_old = volume_cal - volume_exp
         a_d_old = a_d
         ! Update the initial guess slightly.
         a_d = -0.49_pr
         b_d = 0.685_pr - a_d * exp(-0.60_pr)
         call calculate_volume(fluid, rho_method, characterization)
         difference = volume_cal - volume_exp
         
         ! Iteratively adjust a_d (and update b_d) until the difference is 
         ! within tolerance.
         do while (abs(difference) > 0.001_pr)
            aux = a_d
            a_d = a_d - difference * (a_d - a_d_old) / (difference - difference_old)
            b_d = 0.685_pr - a_d * exp(-0.60_pr)
            a_d_old = aux
            difference_old = difference
            call calculate_volume(fluid, rho_method, characterization)
            difference = volume_cal - volume_exp
         end do

      end associate
   end subroutine density_funtion

   subroutine density_method (fluid, rho_method,characterization)
      !! This subroutine compute the experimental volume (volume_exp) using the
      !! fluid's density data.
      !!
      !! Details:
      !!   - For rho_method 1 and 2:
      !!       The experimental volume is calculated using SCN densities 
      !!       provided in the fluid data. The volume for the 6+ fraction is
      !!       computed from the product of mole fraction and molecular weight 
      !!       divided by the corresponding density. method 1 use experimental
      !!       densities for scn, while method 2 used reference values of 
      !!       densities for scn cuts.
      !!   - For rho_method 3:
      !!       A single experimental value is used to compute the density 
      !!       distribution.
      
      implicit none
      type(FluidData) :: fluid
      !! Contains fluid properties and density data
      type(FluidDataOut) :: characterization
      !! Structure to store density function constants and volumes.
      integer, intent(in) :: rho_method
      !! Selector for the density computation method
      real(pr) :: volume_exp
      !! Experimental volume obtained from the fluid data

      select case (rho_method)

      case (1)
         ! Case 1:
         !   Use experimental SCN densities reported to calculate the 
         !   experimental volume for the 6+ fraction. The densities are based
         !   on values calculated from the carbon number plus fraction.
         volume_exp = sum((fluid%product_z_mw_scn)/(fluid%scn_density)) + &
         (fluid%product_z_mw_plus/fluid%plus_density)
         characterization%volume_exp = volume_exp

      case (2)
         ! Case 2:
         !   Use reported (reference) SCN densities to calculate the 
         !   experimental volume. The approach is similar to case 1.
         volume_exp = sum((fluid%product_z_mw_scn)/(fluid%scn_density)) + &
         (fluid%product_z_mw_plus/fluid%plus_density)
         characterization%volume_exp = volume_exp

      case (3)
         ! Case 3:
         !   Use a single experimental value (6plus) to compute the density 
         !   distribution.
         volume_exp = (fluid%plus6_z_exp*fluid%plus6_mw_exp)/(fluid%plus6_density_exp)
         characterization%volume_exp = volume_exp
      
      end select

   end subroutine density_method

   subroutine calculate_volume(fluid, rho_method, characterization)
      !! This subroutine calculate the volume (volume_cal) based on the current
      !! density function parameters.
      !!
      !! Process:
      !!   1. Compute estimated densities for the SCN fraction (scn_density_cal)
      !!      and the plus fraction (plus_density_cal) using the density 
      !!      function constants a_d and b_d.
      !!   2. Depending on the selected rho_method, calculate the total volume 
      !!      by summing contributions:
      !!         - The SCN fraction volume is calculated as the sum over 
      !!           adjusted mole fractions (scn_zm) divided by the corresponding
      !!           density.
      !!         - The plus fraction volume is calculated as the sum over the 
      !!           product of mole fraction and molecular weight for the plus 
      !!           fraction divided by the computed plus density.
      !!   3. Update the characterization structure with the calculated volume 
      !!      and the density array for the 6+ fraction.

      implicit none
      type(FluidData), intent(in) :: fluid 
      !! Contains fluid composition and density data.
      type(FluidDataOut), intent(inout) :: characterization
      !!  Structure to update with calculated volume and density arrays
      integer, intent(in) :: rho_method
      !! Density computation method selector.
      real(pr) :: volume_cal 
      ! Calculated volume based on the current density function
      real(pr), allocatable :: plus_density_cal(:) 
      real(pr), allocatable :: plus6_density(:)
      ! Density values for the 6+ fraction 
      real(pr), allocatable :: scn_density_cal(:)

      associate(&
         a_d => characterization%a_d, &  ! Density function slope constant.
         b_d => characterization%b_d, &  ! Density function intercept constant.
         cn_plus => characterization%carbon_number_plus, &  
         ! Carbon number for the plus fraction.
         z_m_plus_i => characterization%product_z_mw_plus_i, &  
         ! Product of mole fraction and MW for the plus fraction.
         scn_z_i => characterization%scn_z, &   ! SCN mole fractions.
         scn_mw_i => characterization%scn_mw, & ! SCN molecular weights.
         scn_zm => characterization%scn_zm & 
         ! Adjusted SCN mole fractions for volume calculation.
         )
         
         ! Compute estimated density for the SCN fraction using the density function:
         scn_density_cal = ((a_d) * (exp(- (real(fluid%scn, pr))/10._pr))) + b_d
         ! Compute estimated density for the plus fraction using the carbon number plus:
         plus_density_cal =((a_d) * (exp(- (real(cn_plus, pr) )/10._pr))) + b_d
         


         select case (rho_method)

            case(1)
               !   Use the fluid's provided SCN densities and the computed 
               !   plus_density_cal.
               volume_cal = sum((scn_zm)/(fluid%scn_density)) + &
                  sum(z_m_plus_i/plus_density_cal)
               plus6_density = [fluid%scn_density, plus_density_cal]

            case(2)
               !   Use the computed SCN densities (scn_density_cal) and 
               !   plus_density_cal.
               volume_cal = sum((scn_zm)/(scn_density_cal)) + &
                  sum(z_m_plus_i/plus_density_cal)
               plus6_density = [scn_density_cal, plus_density_cal]

            case(3)
               !   Similar to case 2: use computed densities for both SCN and
               !   plus fractions
               volume_cal = sum((scn_zm)/(scn_density_cal)) + &
                  sum(z_m_plus_i/plus_density_cal)
               plus6_density = [scn_density_cal, plus_density_cal]

         end select

         ! Update the characterization structure with the calculated volume and
         ! density values
         characterization%volume_cal = volume_cal
         characterization%plus6_density = plus6_density

      end associate

   end subroutine calculate_volume

end module density
