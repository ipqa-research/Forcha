program main

   use characterization

   implicit none

   type(FluidDataOut) :: fluid_to_characterize
   
   fluid_to_characterize=characterize(file='YPF2.nml',mw_source="experimental",&
         method = "plus_mw", pho_method = 1 , fix_C=.true., eos='PR')
      
   write(*, *) fluid_to_characterize

   print*, fluid_to_characterize%C, fluid_to_characterize%a, fluid_to_characterize%b, &
           fluid_to_characterize%plus_mw, fluid_to_characterize%plus_z
   print*, fluid_to_characterize%c_max
   
end program main


