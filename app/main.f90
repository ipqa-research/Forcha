program main


   use data_from_input, only: data_from_file, FluidData 
   use pruebas

   implicit none
   logical :: start
   real(pr) :: C
   character(20) :: fluid

   
   
   fluid = "oil2.nml"
   call get_C_or_m_plus(fluid, mw_source="experimental", start=start, C=C)
   print*, C

end program main


