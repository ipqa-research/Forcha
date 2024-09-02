program main

   use data_from_input, only: data_from_file, FluidData 
   use routines

   implicit none
   type(FluidDataOut) :: prueba
   type(FluidData) :: fluid
   

   prueba = characterize(file='oil1.nml', mw_source="calculated", method = "plus_mw", fix_C=.true., eos='PR')
   print*, prueba%n_init
   print*, prueba%a , prueba%b
   print*, prueba%C, prueba%plus_mw
   print*, prueba%a_d, prueba%b_d
   print*, "-------------------------------------------------------------------"
   

end program main


