program main

   use data_from_input, only: data_from_file, FluidData 
   use routines

   implicit none
   type(FluidDataOut) :: prueba
   type(FluidData) :: fluid
   

   prueba = characterize(file='oil1.nml', mw_source="calculated", method = "plus_mw", fix_C=.false., eos='PR')

end program main


