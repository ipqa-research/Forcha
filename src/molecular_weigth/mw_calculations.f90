module molecular_weigth
   use constants
   use dtypes, only: FluidData, FluidDataOut
   use bfr_routines
   
contains

   subroutine select_method(fluid, mw_source, method, characterization)
      !! This subroutine defines the calculation method to perform the characterization,&
      !! based on the available experimental data. If the molecular weight data is experimental,
      !! define the mole fractions, molecular weights and densities directly from the input file.
      !! On the other hand, if the molecular weights are assumed, they are recalculated according
      !! to the methodology described by Martin et al and the molar fractions of the fluid are recalculated.

      implicit none
      type(FluidData), intent(inout) :: fluid !! Derived type to save inlet or experimental information of fluid to be characterized.
      type(FluidDataOut), intent(inout) :: characterization !! Derive type variable to save results of characterization methodology.
      character(len=*), intent(in), optional :: method
      !! this option allows choose beetwen obtaining plus_mw or global_mw reported in the experimental information .
      character(len=*), intent(in) :: mw_source !! This variable indicates the source of the input data, that is, whether they are experimental or not.
      real(pr) :: sum_def_comp_z_plus !! compositions sum from (define component +1) to (plus component)
      real(pr), allocatable :: def_comp_moles(:) !!
      real(pr), allocatable :: scn_moles(:) !!
      real(pr) :: plus_moles
      real(pr) :: total_moles

      associate (&
         scn_nc => fluid%scn_nc, def_comp_nc => fluid%def_comp_nc, &
         def_comp_w => fluid%def_comp_w, def_comp_mw => fluid%def_comp_mw, &
         scn_w => fluid%scn_w, plus_w => fluid%plus_w, scn => fluid%scn, &
         plus_mw => characterization%plus_mw, scn_z => characterization%scn_z, &
         plus_z => characterization%plus_z, log_scn_z => characterization%log_scn_z, &
         C => characterization%C, scn_mw => characterization%scn_mw, &
         scn_zm => characterization%scn_zm, &
         plus_zm => characterization%plus_zm &
         )

         allocate(def_comp_moles(def_comp_nc))
         allocate(scn_moles(scn_nc))

         def_comp_moles = (def_comp_w)/(def_comp_mw)

         select case (mw_source)

          case("experimental")
            scn_mw = fluid%scn_mw
            scn_moles = (scn_w)/(scn_mw)
            plus_moles = (plus_w)/(plus_mw)
            total_moles = sum(def_comp_moles)+sum(scn_moles)+(plus_moles)
            scn_z =  scn_moles / total_moles
            plus_z = plus_moles / total_moles

          case("calculated")

            if (method=="global_mw")then
               scn_mw = 84 + C*(scn-6)
               scn_z = fluid%product_z_mw_scn / scn_mw
               sum_def_comp_z_plus = sum(fluid%scn_z) + (fluid%plus_z)
               plus_z = sum_def_comp_z_plus - sum(scn_z) ! new molar fraction for residual cut, based on C
               plus_mw = fluid%product_z_mw_plus / plus_z
            end if

            if (method=="plus_mw")then
               scn_mw = 84 + C*(scn-6)
               scn_moles = (scn_w)/(scn_mw)
               plus_moles = (plus_w)/(plus_mw)
               total_moles = sum(def_comp_moles)+sum(scn_moles)+(plus_moles)
               scn_z =  scn_moles / total_moles
               plus_z = plus_moles / total_moles

            endif
         end select
         scn_zm =  scn_z*scn_mw
         plus_zm = plus_z*plus_mw
         log_scn_z =  log(scn_z)
      end associate

   end subroutine select_method

   subroutine difference_mw_plus(fluid, mw_source, method, start, difference, &
      characterization)

      !! this subroutine ...
      implicit none

      type(FluidData), intent(inout) :: fluid
      type(FluidDataOut), intent(inout) :: characterization
      logical :: start
      character(len=*), intent(in), optional :: method
      character(len=*), intent(in) :: mw_source
      real(pr) :: plus_mw_cal !!  calculated molecular weight of  residual fraction
      real(pr), intent(out) :: difference
      real(pr) :: half
      real(pr) :: r2_old,r2_best, a_old, b_old, a_best, b_best, z_sum, z_aux
      real(pr):: sum_z
      real(pr) :: denom
      !! set to molecular weights calculated by \[M = 84-characterization%C(i-6)\]
      real(pr) :: sum_def_comp_z_plus
      !! compositions sum from (define component +1) to (plus component)
      integer :: j, k , k_old, n_best, x_aux, i_0! internal variables
      integer :: i
      integer, dimension(300) :: carbon_number_plus
      real(pr), dimension(300) :: plus_z_i
      real(pr), dimension(300)  :: product_z_mw_plus_i
      real(pr), dimension(300) :: scn_i

      associate (&
         scn_nc => fluid%scn_nc, scn => fluid%scn,&
         plus_z => characterization%plus_z, plus_mw => characterization%plus_mw, &
         log_scn_z => characterization%log_scn_z, C => characterization%C, &
         a_blr => characterization%a_blr, b_blr => characterization%b_blr, &
         a_60 => characterization%a_60, b_60 => characterization%b_60, &
         a_lim => characterization%a_lim, b_lim => characterization%b_lim, &
         r2 => characterization%r2, n_init => characterization%n_init, &
         c_max_blr => characterization%c_max_lim, a => characterization%a, &
         c_max_lim => characterization%c_max_lim, b => characterization%b, &
         c_max_60 => characterization%c_max_lim, &
         c_max => characterization%c_max &
         )

         1        i_0 = scn(scn_nc)

         call select_method(fluid, mw_source, method, characterization) ! select case
         call Best_Linear_Regression(scn_nc, scn, log_scn_z, &
            plus_z, a_blr, b_blr, r2, n_init, c_max_blr)
         call LimitLine(scn_nc, scn, plus_z, a_blr, b_blr,a_lim,b_lim,c_max_lim,half)
         call Line_C60_max (scn_nc,scn,plus_z,a_blr,b_blr,half,a_60,b_60,c_max_60) ! add line by oscar.

         if(a_blr < a_lim)then
            !Added by oscar 05/12/2023.se elimina por que se agregan las restricciones
            !para characterization%C>12 y characterization%C<12.
            a = a_lim
            b = b_lim
         else  ! this line was removed call Line_C60_max
            if(a_blr > a_60)then
               a = a_60
               b = b_60
            else
               a = a_blr
               b = b_blr
            end if
         end if

         sum_z = 0.0_pr
         i = 0
         do while (sum_z < plus_z .and. i< 300)
            i = i+1
            plus_z_i(i) = exp(a*(i+i_0)+b)
            sum_z = sum_z + plus_z_i(i)
            product_z_mw_plus_i(i) = plus_z_i(i)*(84+C*(i+i_0-6))
            carbon_number_plus(i) = i+i_0
         end do

         if (i == 300.and.sum_z < plus_z)then
            if(start) C = C - 0.07
            if(.not.start) C = C - 0.01
            go to 1
         end if

         c_max = i+i_0
         plus_z_i(i) = plus_z_i(i) - (sum_z-plus_z) !   Adjustment to Zp (z20+)
         product_z_mw_plus_i(i) = plus_z_i(i)*(84+C*(i+i_0-6))

         if (mw_source == "experimental" .and. C>12 .and. C<14 ) then
            ! mayor igual a 12 o menor igual 14.
            denom = sum((plus_z_i(1:i))*((carbon_number_plus(1:i)-6)))
            ! denominator of equation to obtain C value directly
            C = (plus_z*(plus_mw-84))/denom
         endif

         plus_mw_cal = sum(product_z_mw_plus_i(1:i))/plus_z
         difference = plus_mw_cal-plus_mw
         scn_i = [scn, carbon_number_plus(1:i)]
      end associate

      characterization%plus_z_i = plus_z_i(1:i)
      characterization%product_z_mw_plus_i =  product_z_mw_plus_i(1:i)
      characterization%carbon_number_plus =  carbon_number_plus(1:i)
      characterization%nc_plus = i
      characterization%scn_i = scn_i



   end subroutine difference_mw_plus

   subroutine get_c_or_m_plus(fluid, mw_source, method, fix_C, characterization)
      !! This subroutine...
      implicit none
      type(FluidData) :: fluid
      type(FluidDataOut) :: characterization
      logical :: start
      character(len=*), intent(in) :: mw_source
      character(len=*), intent(in), optional :: method
      logical, intent(in) :: fix_C !! If C is < 12 or > 14 then fix to the
      !! closest value.
      !integer :: c_max !! output CN at which characterization%plus_z is reached, as the summation of single z(i) from the best linear distribution (blr)
      real(pr), allocatable :: log_scn_z(:) !! logarithm of corresponding mole fractions of scn cutsn
      real(pr) :: difference !!
      real(pr) :: difference_old, plus_mw_old
      real(pr) :: C_old
      real(pr) :: aux

      associate (C => characterization%C, plus_mw => characterization%plus_mw)

         select case (mw_source)
          case("experimental")
            C = 13.5
            Start = .true.
            plus_mw = fluid%plus_mw
            call difference_mw_plus(fluid, mw_source, method, start, difference, &
               characterization)

          case("calculated")

            if(method == 'global_mw') C=13
            if(method == 'plus_mw')   C=14

            Start = .true.
            plus_mw = fluid%plus_mw
            call difference_mw_plus(fluid, mw_source, method, start, &
               difference_old, characterization)
            C_old = C
            if(method == 'global_mw') C = min(12.8, C-0.07)
            if(method == 'plus_mw') C = 13.5
            Start = .false.
            call difference_mw_plus(fluid, mw_source, method, start, difference, &
               characterization)

            do while (abs(difference) > 0.1)
               aux = C
               C = C - difference*(C-C_old)/(difference-difference_old)
               C_old = aux
               difference_old = difference
               call difference_mw_plus(fluid, mw_source, method, start, &
                  difference, characterization)
            end do
         end select

         if(fix_C .and. C>14 .or. C<12)then


            if(C>14) C=14
            if(C<12) C=12

            Start = .true.
            plus_mw = fluid%plus_mw
            call difference_mw_plus(fluid, mw_source, "plus_mw" , start, &
               difference_old, characterization)
            plus_mw_old = plus_mw
            plus_mw = 0.9*plus_mw
            Start = .false.
            call difference_mw_plus(fluid, mw_source, "plus_mw" , start, &
               difference, characterization)

            do while (abs(difference)>0.00001)
               aux = plus_mw
               plus_mw = plus_mw - difference*(plus_mw-plus_mw_old)/(difference-difference_old)
               plus_mw_old = aux
               difference_old = difference
               call difference_mw_plus(fluid, mw_source, "plus_mw" , start, &
                  difference, characterization)
            end do
         end if
      end associate
   end subroutine get_c_or_m_plus

end module molecular_weigth

