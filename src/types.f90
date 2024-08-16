module dtypes
   use constants, only: pr

   implicit none
   type :: FluidData
      integer :: def_comp_nc !! number of defined components being considered in the oil
      integer :: scn_nc !! number of single cuts being considered in the oil
      integer :: scn_nc_ps !! CN from which all SCN fractions will be lumped into the specified number of pseudos
      integer :: numbers_ps !! number of pseudos in which the scn fractions grouped.
      character(len=:), allocatable :: filename 
      integer, allocatable :: scn (:) !! set of singles cuts being considered in the oil 
      character(len=15), allocatable :: def_components (:) !! set of defined components being considered in the oil  
      character(len=15) :: scn_plus !! name of residual fraction 
      real(pr), allocatable :: def_comp_z(:) !! set of corresponding mole fractions of defined components 
      real(pr), allocatable :: scn_z(:) !!  set of corresponding mole fractions of scn cuts 
      real(pr), allocatable :: def_comp_mw(:) !! set of corresponding molecular weights of defined components 
      real(pr), allocatable :: scn_mw(:) !! set of corresponding molecular weights of scn cuts
      real(pr), allocatable:: product_z_mw_def_comp(:) !! product between composition and molecular weight of defind components
      real(pr), allocatable:: product_z_mw_scn(:) !! product between composition and molecular weight of scn fractions
      real(pr) :: sum_z_mw_i !!  sum of the product between composition and molecular weight of the fluid's compounds
      real(pr), allocatable :: w(:) !! mass fractions of the fluid's compounds
      real(pr) :: plus_z !! composition of residual fraction
      real(pr) :: plus_mw !!  molecular weight of residual fraction 
      real(pr) :: product_z_mw_plus !! product between composition and molecular weight of residual fraction
      real(pr), allocatable :: def_comp_w(:) !! mass fractions of the defined compounds
      real(pr), allocatable :: scn_w(:) !! !! mass fractions of the scn-s compounds
      real(pr) :: plus_w !! mass fractions of the plus fraction
      real(pr), allocatable :: scn_density(:) !! set of corresponding densities of scn cuts
      real(pr) :: plus_density !! experimental density of the plus fraction
   end type FluidData

   type :: FluidDataOut
      ! incluir aqui las variables que quiero que salgan como salida
      real(pr), allocatable :: scn_z(:) !!  set of corresponding mole fractions of scn cuts calculated
      real(pr), allocatable :: log_scn_z(:) !!  set of corresponding mole fractions of scn cuts calculated
      real(pr) :: plus_mw !!  molecular weight of residual fraction
      real(pr) :: plus_z !! composition of residual fraction
      real(pr) :: C !! C constants which is used in equation \[M = 84-characterization%C(i-6)\]
      real(pr), allocatable :: scn_mw(:) !! set of corresponding molecular weights of scn cuts
      real(pr), allocatable :: scn_zm (:)
      real(pr) :: plus_zm
      real(pr) :: a_blr !! A constant for best linear regression line.
      real(pr) :: b_blr !! B constant for best linear regression line.
      real(pr) :: a_60 !! A constant for Cmax60 line.
      real(pr) :: b_60 !! B constant for Cmax60 line.
      real(pr) :: a_lim !! A constant for limit feasible line.
      real(pr) :: b_lim !! B constant for limit feasible line.
      real(pr) :: a !! output real variable. Slope of the best feasible regression line.
      real(pr) :: b !! output real variable. Intercept of the best feasible regression line.
      real(pr) :: r2 !! output real variable. Square correlation coefficient.
      integer :: n_init ! minimum carbon number obtained from the best linear regression
      integer :: c_max_blr, c_max_lim, c_max_60
      integer :: nc_plus !! number of elements in the distribution of C20+
      integer :: c_max  !! output CN at which plus_z is reached, as the summation of single z(i) from the best linear distribution (blr)
      real(pr), allocatable :: carbon_number_plus(:)
      real(pr), allocatable :: plus_z_i(:)
      real(pr), allocatable  :: product_z_mw_plus_i(:)
      real(pr) :: a_d !! ad constant which is used in equation \[rho_i = ad*exp(-i/10) +bd\] for density.
      real(pr) :: b_d !! bd constant which is used in equation \[rho_i = ad*exp(-i/10) +bd\] for density.
      real(pr) :: volume_6plus_cal
      integer :: last_C
      integer :: i_last
      real(pr), allocatable :: scn_i(:)
      real(pr), allocatable :: plus6_density(:)
      real(pr) :: plus_w
      real(pr) :: plus_density
      real(pr), allocatable :: mol_fraction(:)
   end type FluidDataOut

end module dtypes


