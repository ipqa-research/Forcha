module characterization
   
   use constants
   use dtypes, only: FluidData, FluidDataOut
   use data_from_input, only: data_from_file, FluidData
   use molecular_weigth
   use density
   use lumping
   use critical_parameters
   
   
contains

   type(FluidDataOut) function characterize(file, mw_source, method, pho_method, fix_C, eos) &
      result(characterization)

      type(FluidData) :: fluid
      character(len=*), intent(in) :: file !! file name
      character(len=*), intent(in), optional :: method !! plus_mw or global_mw
      character(len=*), intent(in) :: mw_source
      logical, intent(in), optional :: fix_C
      character(len=*), intent(in) :: eos
      integer, intent(in) :: pho_method !! selector method of density funtion.



      fluid = data_from_file(file)
      characterization%input_data = fluid
      allocate(characterization%scn_z(fluid%scn_nc))
      allocate(characterization%log_scn_z(fluid%scn_nc))
      allocate(characterization%scn_mw(fluid%scn_nc))
      allocate(characterization%carbon_number_plus(0))
      allocate(characterization%plus_z_i(0))
      allocate(characterization%product_z_mw_plus_i(0))
      allocate(characterization%scn_zm(fluid%scn_nc))
      allocate(characterization%scn_i(0))
      allocate(characterization%plus6_density(0))
      allocate(characterization%mol_fraction(0))
      allocate(characterization%lumped_z(0))
      allocate(characterization%lumped_mw(0))
      allocate(characterization%lumped_densities(0))
      allocate(characterization%critical_temperature(0))
      allocate(characterization%critical_pressure(0))
      allocate(characterization%acentric_factor(0))
      allocate(characterization%m_funtion(0))


      call get_c_or_m_plus(fluid=fluid, mw_source=mw_source, method=method, fix_C=fix_C, characterization= characterization)
      call density_funtion(fluid=fluid, mw_source=mw_source, pho_method= pho_method, characterization=characterization)
      call lump(fluid=fluid, characterization=characterization)
      call get_critical_constants(fluid= fluid, characterization=characterization, eos=eos)


   end function characterize

end module characterization


