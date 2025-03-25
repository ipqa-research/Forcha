program main

   use characterization

   implicit none

   type(FluidDataOut) :: fluid_to_characterize

   fluid_to_characterize=characterize(file='YPF2.nml', mw_source="experimental",&
      method = "plus_mw", rho_method = 1 , fix_C=.TRUE., eos='PR')

   call general_result( fluid_to_characterize)
   write(*, *) fluid_to_characterize

end program main


