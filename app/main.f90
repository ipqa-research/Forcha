program main

   use data_from_input, only: data_from_file, FluidData 
   use routines

   implicit none
   type(FluidData) :: oil

   oil = data_from_file("oil2.nml")
   !call get_c_or_m_plus(fluid=oil, mw_source="calculated", method= "plus_mw")
   call get_c_or_m_plus(fluid=oil, mw_source="experimental")

end program main


