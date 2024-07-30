program main

   use data_from_input, only: data_from_file, FluidData 
   use routines

   implicit none
   integer ::  k, i
   type(FluidDataOut) :: prueba
   type(FluidData) :: fluid
   

   !prueba = characterize(file='oil1.nml', mw_source="experimental")
   !print*, prueba%C, prueba%plus_mw
   !prueba = characterize(file='oil1.nml', mw_source="calculated", method = "plus_mw")
   !print*, prueba%C, prueba%plus_mw
   !prueba = characterize(file='oil1.nml', mw_source="calculated", method = "global_mw")
   !print*, prueba%C, prueba%plus_mw
   !prueba = characterize(file='oil2.nml', mw_source="experimental")
   !print*, prueba%C, prueba%plus_mw
   !prueba = characterize(file='oil2.nml', mw_source="calculated", method = "plus_mw")
   !print*, prueba%C, prueba%plus_mw
   !prueba = characterize(file='oil2.nml', mw_source="calculated", method = "global_mw")
   !print*, prueba%C, prueba%plus_mw


   prueba = characterize(file='YPF1.nml', mw_source="calculated", method = "plus_mw", fix_C=.true.)
   print*, prueba%n_init
   print*, prueba%a , prueba%b
   print*, prueba%C, prueba%plus_mw
   print*, prueba%a_d, prueba%b_d
   print*, "-------------------------------------------------------------------"
   do i = 1, size(prueba%scn_z)
      print*, prueba%scn_z(i)
   end do 
   print*, prueba%plus_z
   print*, "-------------------------------------------------------------------"
   do i = 1, size(prueba%scn_mw)
      print*, prueba%scn_mw(i)
   end do 
   print*, prueba%plus_mw
   print*, "-------------------------------------------------------------------"
   print*, prueba%c_max


   !oil =  data_from_file('YPF1.nml')
   
   !call get_c_or_m_plus(fluid= oil, mw_source="calculated", method="global_mw", result)

   
   !fluid = data_from_file(file='oil1.nml')
   !print*, fluid%scn
   !k = prueba%nc_plus
   !print*, k
   !print*, prueba%nc_plus
   !do i = 1, k
   !   print*, prueba%carbon_number_plus(i)
   !end do

   !print*, prueba%carbon_number_plus
   !print*, prueba%nc_plus
   !print*, prueba%plus_z_i
  

end program main


