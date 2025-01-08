module density
   use constants
   use dtypes, only: FluidData, FluidDataOut

contains

   subroutine density_funtion(fluid, mw_source, pho_method, characterization)
      !! this subroutine ...
      implicit none
      type(FluidData) :: fluid
      type(FluidDataOut) :: characterization
      character(len=*), intent(in) :: mw_source
      integer, intent(in) :: pho_method
      !real(pr), allocatable :: volume_6plus_cal(:)
      real(pr) :: volume_exp
      real(pr) :: a_d_old, aux
      real(pr) :: difference, difference_old


      !volume_exp = sum((fluid%product_z_mw_scn)/(fluid%scn_density)) + &
      !   (fluid%product_z_mw_plus/fluid%plus_density)




      

      associate( &
         a_d => characterization%a_d , &
         b_d => characterization%b_d, &
         volume_cal => characterization%volume_cal, &
         volume_exp => characterization%volume_exp &
         )

         call density_method (fluid, pho_method,characterization)

         ! Now find constants for the Density function
         a_d = -0.50_pr ! initial guess
         b_d = 0.685_pr - a_d*exp(-0.60_pr)
         call calculate_volume(fluid, mw_source, pho_method, characterization)
         difference_old = volume_cal - volume_exp
         a_d_old = a_d
         a_d = -0.49_pr
         b_d = 0.685_pr - a_d*exp(-0.60_pr)
         call calculate_volume(fluid, mw_source, pho_method, characterization)
         difference = volume_cal - volume_exp

         do while (abs(difference) > 0.001_pr)
            aux = a_d
            a_d = a_d - difference*(a_d - a_d_old)/(difference - difference_old)
            b_d = 0.685_pr - a_d*exp(-0.60_pr)
            a_d_old = aux
            difference_old = difference
            call calculate_volume(fluid, mw_source, pho_method, characterization)
            difference = volume_cal - volume_exp
         end do

      end associate
   end subroutine density_funtion

   subroutine density_method (fluid, pho_method,characterization)
            !! this subroutine ...
      implicit none
      type(FluidData) :: fluid
      type(FluidDataOut) :: characterization
      integer, intent(in) :: pho_method
      real(pr) :: volume_exp

      select case (pho_method) 

      case (1)
         ! case 1 use experimental scn's densities reported to calculate 
         ! experimental volume for 6plus. The densities values are calculated 
         ! since CN plus
         volume_exp = sum((fluid%product_z_mw_scn)/(fluid%scn_density)) + &
         (fluid%product_z_mw_plus/fluid%plus_density)
         characterization%volume_exp = volume_exp

      case (2)
         !This case use scn densities  reported to calculate experimental volume
         ! this densities doesn't experimental, they can be reference values.
         ! the densities values are calculated since 6 plus

         volume_exp = sum((fluid%product_z_mw_scn)/(fluid%scn_density)) + &
         (fluid%product_z_mw_plus/fluid%plus_density)
         characterization%volume_exp = volume_exp

      case (3)
         ! This case use one experimental value (6plus) to calculate the density 
         ! distribution

         volume_exp = (fluid%plus6_z_exp*fluid%plus6_mw_exp)/(fluid%plus6_density_exp)
         characterization%volume_exp = volume_exp
      
      end select








   end subroutine density_method

   subroutine calculate_volume(fluid, mw_source, pho_method, characterization)
      !! this subroutine ...
      implicit none
      type(FluidData), intent(in) :: fluid
      type(FluidDataOut), intent(inout) :: characterization
      character(len=*), intent(in) :: mw_source
      integer, intent(in) :: pho_method
      real(pr) :: volume_cal
      real(pr), allocatable :: plus_density_cal(:), scn_density_cal(:), plus6_density(:)

      !allocate(density_plus_cal(0), density_scn_cal(0), plus6_density(0))


      associate(&
         a_d => characterization%a_d, &
         b_d => characterization%b_d, &
         cn_plus => characterization%carbon_number_plus, &
         z_m_plus_i => characterization%product_z_mw_plus_i, &
         scn_z_i => characterization%scn_z, &
         scn_mw_i => characterization%scn_mw, &
         scn_zm => characterization%scn_zm &
         )

         scn_density_cal = ((a_d)*(exp(- (real(fluid%scn, pr))/10._pr))) + b_d
         plus_density_cal =((a_d)*(exp(- (real(cn_plus, pr) )/10._pr))) + b_d


         select case (pho_method)

            case(1)
               volume_cal = sum((scn_zm)/(fluid%scn_density)) + &
                  sum(z_m_plus_i/plus_density_cal)
               plus6_density = [fluid%scn_density, plus_density_cal]

            case(2)
               volume_cal = sum((scn_zm)/(scn_density_cal)) + &
                  sum(z_m_plus_i/plus_density_cal)
               plus6_density = [scn_density_cal, plus_density_cal]

            case(3)
               volume_cal = sum((scn_zm)/(scn_density_cal)) + &
                  sum(z_m_plus_i/plus_density_cal)
               plus6_density = [scn_density_cal, plus_density_cal]

         end select

         characterization%volume_cal = volume_cal
         characterization%plus6_density = plus6_density

      end associate

   end subroutine calculate_volume


end module density
