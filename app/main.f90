program main

   use data_from_input, only: data_from_file, FluidData 
   use routines

   implicit none
   type(FluidData) :: oil
   
  
   oil = data_from_file("oil1.nml")
   print *, oil%filename
   call get_c_or_m_plus(fluid=oil, mw_source="experimental")
   
   oil = data_from_file("oil2.nml")
   print *, oil%filename
   call get_c_or_m_plus(fluid=oil, mw_source="experimental")
   

end program main


