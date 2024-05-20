program main


   use data_from_input, only: data_from_file, FluidData 
   use routines
   !use constants, only: pr
   !use routines
   !use data

   implicit none
   logical :: start
   real(pr) :: C
   real(pr) :: difference
   real(pr), allocatable :: scn_z(:)!! calculated composition of residual fraction
   real(pr), allocatable :: log_scn_z(:) !! logarithm of corresponding mole fractions of scn cuts
   character(20) :: fluid
   real(pr) :: plus_mw !!  calculated molecular weight of  residual fraction
   
   
   
   
   type(FluidData) :: asd
   fluid = "oil2.nml"
   !call get_C_or_m_plus(fluid,"global_mw","calculated",start,C)
   !print*, C
   !call get_C_or_m_plus(fluid,"plus_mw","calculated",start,C)
   !print*, C
   call get_C_or_m_plus(fluid,"plus_mw","experimental",start,C)
   print*, C

end program main


