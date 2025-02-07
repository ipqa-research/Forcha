module lumping
   use dtypes

contains

   subroutine lump(fluid, characterization)
      !! This sobroutine ...
      implicit none
      type(FluidData), intent(in) :: fluid
      type(FluidDataOut), intent(inout) :: characterization
      integer :: last_C ! CN of last single cut: typically 19
      integer :: i_last ! Number of elements in the distribution of CN+
      integer :: last   ! Total elements starting after defined components.
      integer ::  scn_nc_input !! inicial number of single cuts being considered in the oil from data input
      integer :: scn_nc_new !! new number of single cuts being considered in the oil defined as from scn_nc_ps variable.
      real(pr) :: plus_w_new
      real(pr), allocatable :: z_m_plus_i(:)
      real(pr), allocatable :: plus_z_i(:)
      real(pr), allocatable :: var_aux_1(:)
      real(pr), allocatable :: var_aux_2(:)
      integer :: i, prev_i
      integer, dimension(15) :: j_ps
      integer :: i_ps       ! auxiliary variables
      real(pr) :: rec_zm, remain_plus_zm, sum_z, sum_zm, sum_volume_ps  ! auxiliary variables
      real(pr), allocatable ::  plus_z_ps(:)
      real(pr), allocatable :: plus_mw_ps(:)
      integer :: numbers_ps
      real(pr), allocatable :: density_ps(:)
      real(pr), allocatable :: w_ps(:)
      real(pr) :: sum_zm_last_ps
      real(pr), allocatable :: moles(:)
      real(pr), allocatable :: lumped_z(:)
      real(pr), allocatable :: lumped_mw(:)
      real(pr), allocatable :: lumped_densities(:)

      associate (&
         def_nc => fluid%def_comp_nc ,  &
         scn_mw => characterization%scn_mw,  &
         plus_zm => characterization%plus_zm,  &
         plus_density => characterization%plus_density,  &
         scn_zm => characterization%scn_zm, &
         plus_z => characterization%plus_z, &
         scn_z => characterization%scn_z, &
         plus_w => fluid%plus_w, &
         w => fluid%w, &
         plus_mw => characterization%plus_mw, &
         carbon_number_plus => characterization%carbon_number_plus, &
         plus6_density => characterization%plus6_density, &
         C => characterization%C, &
         scn_nc_new => characterization%scn_nc_new &
         )

         last_C = fluid%scn(fluid%scn_nc)    !19
         i_last = characterization%nc_plus   !124
         last = fluid%scn_nc + i_last        !138

         allocate (z_m_plus_i(0))
         allocate (plus_z_i(0))
         allocate (var_aux_1(0))
         allocate (var_aux_2(0))
         allocate (plus_z_ps(0))
         allocate (plus_mw_ps(0))
         allocate(density_ps(0))
         allocate (w_ps(0))
         allocate(moles(0))
         allocate(lumped_z(0))
         allocate(lumped_mw(0))
         allocate(lumped_densities(0))


         scn_nc_input =  fluid%scn_nc
         scn_nc_new = fluid%scn_nc_ps - 6
         ! scn_nc_ps : CN from which all SCN fractions will be lumped
         ! into the specified number of pseudos

         characterization%plus_w = plus_w
         z_m_plus_i = characterization%product_z_mw_plus_i
         plus_z_i = characterization%plus_z_i

         if (scn_nc_new < scn_nc_input) then
            ! Plus Fraction needs to be extended to include lower CN's
            plus_zm = plus_zm + sum(scn_zm(scn_nc_new + 1 : scn_nc_input))
            plus_z = plus_z + sum(scn_z(scn_nc_new + 1 : scn_nc_input)) !extended zp
            plus_w_new = plus_w + sum(w(fluid%def_comp_nc+scn_nc_new + 1 &
               : fluid%def_comp_nc+scn_nc_input))
            !! use new variable for plus_w because type fluidata can't be modified.
            plus_mw = plus_zm / plus_z
            z_m_plus_i = scn_zm(scn_nc_new + 1 : scn_nc_input )
            z_m_plus_i = [z_m_plus_i, characterization%product_z_mw_plus_i]
            plus_z_i = scn_z(scn_nc_new + 1 : scn_nc_input )
            plus_z_i = [plus_z_i, characterization%plus_z_i]
            i_last = i_last + scn_nc_input - scn_nc_new
            last_C = last_C - scn_nc_input + scn_nc_new
            characterization%plus_w = plus_w_new
            characterization%product_z_mw_plus_i = z_m_plus_i
            characterization%plus_z_i = plus_z_i
         end if

         !do i= 1, scn_nc_new
         !      print*, scn_z(i), scn_mw(i)
         !end do

         numbers_ps = fluid%numbers_ps
         j_ps = 0
         ! Lumping into Nps pseudos
         rec_zm = plus_zm / numbers_ps
         ! Recommended value for the product z*M (proportional to weight)
         ! for each pseudo according to pedersen
         remain_plus_zm = plus_zm
         i_ps = 1._pr
         sum_z = 0.0_pr
         sum_zm = 0.0_pr
         sum_volume_ps = 0.0_pr

         i = 0.0_pr

         do while (i_ps < numbers_ps)
            i = i + 1
            j_ps(i_ps) = j_ps(i_ps) + 1
            sum_z = sum_z + plus_z_i(i)
            sum_zm = sum_zm + z_m_plus_i(i)
            sum_volume_ps = sum_volume_ps + z_m_plus_i(i) / plus6_density(scn_nc_new+i)
            if (z_m_plus_i(i+1) > 2*(rec_zm - sum_zm)) then
               ! when adding one more would go too far
               plus_z_ps  = [plus_z_ps, sum_z]
               plus_mw_ps = [plus_mw_ps, sum_zm/sum_z]
               w_ps = [w_ps, (characterization%plus_w*(sum_zm/plus_zm))]
               remain_plus_zm = remain_plus_zm - sum_zm
               if (remain_plus_zm < rec_zm) numbers_ps = i_ps + 1
               carbon_number_plus(scn_nc_new + i_ps) = 6+((plus_mw_ps(i_ps)-84)/C)
               density_ps = [density_ps, (sum_zm/sum_volume_ps)]
               i_ps = i_ps + 1
               sum_z = 0._pr
               sum_zm = 0._pr
               sum_volume_ps = 0._pr
            end if
         end do

         plus_z_ps = [plus_z_ps, sum(plus_z_i(i+1:i_last))]
         ! at this point, Nps is the order for the last NON-Asphaltenic pseudo comp.
         ! (e.g. 4 if 5 is for Asphaltenes)
         plus_mw_ps =[plus_mw_ps, sum(z_m_plus_i(i+1:i_last))/plus_z_ps(numbers_ps)]
         w_ps = [w_ps, characterization%plus_w * ((plus_z_ps(numbers_ps))* &
            (plus_mw_ps(numbers_ps))/(plus_zm))]
         carbon_number_plus(scn_nc_new + numbers_ps) = 6 + &
            ((plus_mw_ps(numbers_ps) - 84) / C)
         sum_zm_last_ps = (sum(z_m_plus_i(i+1:i_last)))
         density_ps = [density_ps, (sum_zm_last_ps)/sum((z_m_plus_i(i+1:i_last)) / (plus6_density(scn_nc_new+i+1:last)))]
         j_ps(numbers_ps) = i_last - sum(j_ps(1: numbers_ps-1)) !! revisar manana

         !do i = 1, numbers_ps
         !   print*, 'ps',i, plus_z_ps(i), plus_mw_ps(i), density_ps(i), j_ps(i), carbon_number_plus(scn_nc_new+i)
         !end do

         plus_density = plus_zm / sum(plus_z_ps(1:numbers_ps)* &
            plus_mw_ps(1:numbers_ps)/(density_ps(1:numbers_ps)))

         moles = [(fluid%def_comp_w)/(fluid%def_comp_mw)]
         moles = [moles, fluid%w(def_nc+1:def_nc+scn_nc_new)/ scn_mw(1:scn_nc_new)]
         moles = [moles, w_ps/plus_mw_ps ]

         characterization%mol_fraction = moles/sum(moles) ! mole fractions normalize
         characterization%last_C = last_C
         characterization%i_last = i_last
         characterization%last = last

         ! In this point we create a new array to content final Mis and densities
         ! whit pseudos.

         lumped_z = [scn_z(1:scn_nc_new), plus_z_ps]
         lumped_mw = [scn_mw(1:scn_nc_new), plus_mw_ps]
         lumped_densities = [plus6_density(1:scn_nc_new), density_ps ]


         characterization%lumped_z = lumped_z
         characterization%lumped_mw = lumped_mw
         characterization%lumped_densities = lumped_densities

      end associate

   end subroutine lump

end module lumping

