module dtypes
   use constants, only: pr

   implicit none
   type :: FluidData
      integer :: def_comp_nc !! number of defined components being considered in the oil
      integer :: scn_nc !! number of single cuts being considered in the oil
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



   end type FluidData


   type :: FluidDataOut
   ! incluir aqui las variables que quiero que salgan como salida
   end type FluidDataOut

end module dtypes


