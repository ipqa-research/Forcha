program main


   use data_from_input, only: data_from_file, FluidData 
   use routines
   !use constants, only: pr
   !use routines
   !use data

   implicit none
   logical :: start
   real(pr) :: abe, a_lim, a_60! output real variable. Slope of the best regression line.
   real(pr) :: bbe , b_lim, b_60 !! output real variable. Intercept of the best regression line.
   real(pr):: r2, half!! output real variable. Square correlation coefficient.
   integer :: ninit !! minimum carbon number obtained from the best linear regression
   integer :: c_max_blr, c_max_lim, c_max_60, c_max !
   real(pr), allocatable :: ylog(:)
   real(pr) :: C
   real(pr) :: difference
   real(pr) :: plus_mw !!  calculated molecular weight of  residual fraction
   real(pr) :: plus_z!! calculated composition of residual fraction
   real(pr), allocatable :: scn_z(:)!! calculated composition of residual fraction
   real(pr), allocatable :: log_scn_z(:) !! logarithm of corresponding mole fractions of scn cuts
   character(20) :: fluid
   integer :: i 
   
   
   type(FluidData) :: asd
   asd = data_from_file(file="oil1.nml")

   fluid = "oil2.nml"


   call get_C_or_m_plus(fluid,"global_mw","calculated",start,C)

end program main


