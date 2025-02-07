program main

   !use data_from_input, only: data_from_file, FluidData 
   use characterization

   implicit none

   type(FluidDataOut) :: prueba
   integer :: i
   
   prueba = characterize(file='PVT5.nml', mw_source="experimental", method = "plus_mw", pho_method = 3 , fix_C=.true., eos='PR')
   write(*, *) prueba
   print*, prueba%C, prueba%a, prueba%b, prueba%plus_mw, prueba%plus_z
   print*, prueba%c_max
   
   !print*, prueba%a ,  prueba%b , prueba%n_init, prueba%input_data%number_plus_density
   !do i = 1, size(prueba%mol_fraction)
   !   print*, prueba%mol_fraction(i)
   !end do

end program main


