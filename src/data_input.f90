module data_from_input
   !! This module reads the fluid information from the corresponding input file.
   use dtypes, only: FluidData
   use constants, only: pr
   implicit none
contains

   subroutine read_setup(file, def_comp_nc, scn_nc)
      !! Reads the setup data from input file
      character(len=*), intent(in) :: file !! file name
      integer, intent(out):: def_comp_nc !! number of defined components being considered in the oil
      integer, intent(out):: scn_nc !! number of single cuts being considered in the oil
      integer :: funit
      namelist /nml_setup/ def_comp_nc, scn_nc

      open(newunit=funit, file=file)
      read(funit, nml=nml_setup)
      close(funit)
   end subroutine read_setup

   subroutine read_components(file, def_components, scn, scn_plus)
      !! Reads the component names data from input file
      character(len=*), intent(in) :: file !! file name
      character(len=15), allocatable, intent(out) :: def_components(:) !! names of defined components
      integer, allocatable, intent(out) :: scn(:) !! names of scn fractions
      character(len=15), intent(out) :: scn_plus !! name of residual fraction
      integer :: def_comp_nc !! number of defined components being considered in the oil
      integer :: scn_nc !! number of single cuts being considered in the oil
      integer :: funit
      namelist /nml_components/ def_components, scn, scn_plus

      call read_setup(file, def_comp_nc, scn_nc)
      allocate(def_components(def_comp_nc))
      allocate(scn(scn_nc))

      open(newunit=funit, file=file)
      read(funit, nml=nml_components)
      close(funit)
   end subroutine read_components

   subroutine read_composition(file, def_comp_z, scn_z, plus_z)
      !! Reads the molar compositions of each component from input file
      character(len=*), intent(in) :: file !! file name
      real(pr), allocatable, intent(out) :: def_comp_z(:) !! set of corresponding mole fractions of defined components
      real(pr), allocatable, intent(out) :: scn_z(:) !! set of corresponding mole fractions of scn cuts
      real(pr), intent(out):: plus_z !! composition of residual fraction from input file
      integer :: def_comp_nc !! number of defined components being considered in the oil
      integer :: scn_nc !! number of single cuts being considered in the oil
      integer :: funit
      namelist /nml_composition/ def_comp_z, scn_z, plus_z

      call read_setup(file, def_comp_nc, scn_nc)
      allocate(def_comp_z(def_comp_nc))
      allocate(scn_z(scn_nc))

      open(newunit=funit, file=file)
      read(funit, nml=nml_composition)
      close(funit)
   end subroutine read_composition

   subroutine read_molecular_weight(file, def_comp_mw, scn_mw, plus_mw)
      !! Reads the molecular weights of each component from the input file
      character(len=*), intent(in) :: file !! file name
      real(pr), allocatable, intent(out) :: def_comp_mw(:) !! set of corresponding molecular weights of defined components
      real(pr), allocatable, intent(out) :: scn_mw(:) !! set of corresponding molecular weights of scn cuts
      real(pr), intent(out) :: plus_mw !!  molecular weight of residual fraction
      integer :: def_comp_nc !! number of defined components being considered in the oil
      integer :: scn_nc !! number of single cuts being considered in the oil
      integer :: funit
      namelist /nml_molecular_weight/ def_comp_mw, scn_mw, plus_mw

      call read_setup(file, def_comp_nc, scn_nc)
      allocate(def_comp_mw(def_comp_nc))
      allocate(scn_mw(scn_nc))

      open(newunit=funit, file=file)
      read(funit, nml=nml_molecular_weight)
      close(funit)
   end subroutine read_molecular_weight

   subroutine mass_fractions(file, w, product_z_mw_def_comp, product_z_mw_scn,&
      sum_z_mw_i, product_z_mw_plus, def_comp_w, scn_w, plus_w)
      !! This routine obtains the mass fractions of the fluid's compounds from
      !! the input file, additionally, it calculates the product \[(z*mw)_{i}\]
      character(len=*), intent(in) :: file !! file name
      real(pr), allocatable, intent(out) :: product_z_mw_def_comp(:) !! product between composition and molecular weight of defind components
      real(pr), allocatable, intent(out) :: product_z_mw_scn(:) !! product between composition and molecular weight of scn fractions
      real(pr), intent(out) :: sum_z_mw_i !! sum of the product between composition and molecular weight of the fluid's compounds
      real(pr), allocatable, intent(out) :: w(:) !! mass fractions of the fluid's compounds
      real(pr), intent(out) :: product_z_mw_plus !! product between composition and molecular weight of residual fraction
      integer :: def_comp_nc !! number of defined components being considered in the oil
      integer :: scn_nc !! number of single cuts being considered in the oil
      integer :: funit
      real(pr), allocatable :: def_comp_z(:) !! set of corresponding mole fractions of defined components
      real(pr), allocatable :: scn_z(:) !! set of corresponding mole fractions of scn cuts
      real(pr) :: plus_z !! composition of residual fraction
      real(pr), allocatable :: def_comp_mw(:) !! set of corresponding molecular weights of defined components
      real(pr), allocatable  :: scn_mw(:) !! set of corresponding molecular weights of scn cuts
      real(pr) :: plus_mw !!  molecular weight of residual fraction
      real(pr), allocatable :: product_z_mw_i(:)
      real(pr), allocatable, intent(out) :: def_comp_w(:)
      real(pr), allocatable, intent(out) :: scn_w(:)
      real(pr), intent(out) :: plus_w

      call read_setup(file, def_comp_nc, scn_nc)
      allocate(product_z_mw_def_comp(def_comp_nc))
      allocate(product_z_mw_scn(scn_nc))
      allocate(def_comp_z(def_comp_nc))
      allocate(scn_z(scn_nc))
      allocate(def_comp_mw(def_comp_nc))
      allocate(scn_mw(scn_nc))
      allocate(w(1+def_comp_nc+scn_nc))
      allocate(product_z_mw_i(1+def_comp_nc+scn_nc))
      allocate(def_comp_w(def_comp_nc))
      allocate(scn_w(scn_nc))

      call read_composition(file, def_comp_z, scn_z, plus_z)
      call read_molecular_weight(file, def_comp_mw, scn_mw, plus_mw)

      product_z_mw_def_comp = (def_comp_z)*(def_comp_mw)
      product_z_mw_scn = (scn_z)*(scn_mw)
      product_z_mw_plus = (plus_z)*(plus_mw)
      product_z_mw_i = [product_z_mw_def_comp, product_z_mw_scn, product_z_mw_plus]
      sum_z_mw_i = sum(product_z_mw_i)
      w = (product_z_mw_i) / (sum_z_mw_i)
      def_comp_w = w(1:def_comp_nc)
      scn_w = w(def_comp_nc+1: def_comp_nc+scn_nc)
      plus_w = w(def_comp_nc+scn_nc+1)
   end subroutine mass_fractions

   subroutine read_density()
      !! Reads the density of each component from the input file and calculated molar volume since C6 fraction.




   end subroutine read_density



   type(FluidData) function data_from_file(file) result(data)
      !! This funtion allows to obtain experimental data from data imput
      character(len=*), intent(in) :: file !! file name

      data%filename = file

      call read_setup(file, data%def_comp_nc, data%scn_nc)
      call read_components(file, data%def_components, data%scn, &
         data%scn_plus)
      call read_composition(file, data%def_comp_z, data%scn_z, data%plus_z)
      call read_molecular_weight(file, data%def_comp_mw, data%scn_mw, &
         data%plus_mw)
      call mass_fractions(file, data%w, data%product_z_mw_def_comp, &
         data%product_z_mw_scn, data%sum_z_mw_i, data%product_z_mw_plus, &
         data%def_comp_w, data%scn_w, data%plus_w)
 
         
   end function data_from_file

end module data_from_input









