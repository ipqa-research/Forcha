module dtypes
   use constants, only: pr

   implicit none
   type :: FluidData
      integer :: def_nc
      !! number of defined components before scn from input file
      integer :: scn_nc
      !! number of single carbon number components from input file
      integer, allocatable :: scn_names (:)
      !! names of scn fractions from input file
      character(len=15), allocatable :: def_comp_names(:)
      !! names of defined components before scn from input file
      character(len=15) :: plus_name
      !! name of plus fraction from input file
      real(pr), allocatable :: def_comp_z(:)
      !! compositions of defined components from input file
      real(pr), allocatable :: scn_z(:)
      !! compositions of single carbon numbersfrom input file
      real(pr), allocatable :: def_comp_mw(:)
      !! molecular weights of defined components from input file
      real(pr), allocatable :: scn_mw(:)
      !! molecular weights of single carbon numbers from input file
      real(pr) :: plus_z !! composition of residual fraction from input file
      real(pr) :: plus_mw
      !!  molecular weight of residual fraction from input file
   end type FluidData

end module dtypes


