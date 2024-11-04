program main

   use data_from_input, only: data_from_file, FluidData 
   use routines

   implicit none
   type(FluidDataOut) :: prueba
   
   prueba = characterize(file='YPF2.nml', mw_source="calculated", method = "plus_mw", fix_C=.false., eos='PR')
   write(*, *) prueba
   print*, prueba%a ,  prueba%b , prueba%n_init

   !do i = 1, size(prueba%mol_fraction)
   !   print*, prueba%mol_fraction(i)
   !end do
end program main


