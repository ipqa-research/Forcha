module data_from_input
   !! This module, data_from_input, is designed to read fluid information from
   !! an input file and organize it into a data structure.
   use dtypes, only: FluidData
   ! The FluidData type is imported, which stores all the fluid information.
   use constants, only: pr 
   ! Pr precision is used to define floating point variables.

   implicit none

contains

   subroutine read_setup(file, def_comp_nc, scn_nc, scn_nc_ps, numbers_ps)
      !! This subroutine reads general fluid parameters from input file

      character(len=*), intent(in) :: file !! file name
      integer, intent(out):: def_comp_nc
      !! Number of defined components being considered in the oil
      integer, intent(out):: scn_nc
      !! Number of single cuts being considered in the oil
      integer, intent(out):: scn_nc_ps
      !! CN from which all SCN fractions will be lumped into the specified
      !! number of pseudos
      integer, intent(out):: numbers_ps
      !! Number of pseudos in which the scn fractions grouped
      integer :: funit ! environment variable

      namelist /nml_setup/ def_comp_nc, scn_nc, scn_nc_ps, numbers_ps

      open(newunit=funit, file=file)
      read(funit, nml=nml_setup)
      close(funit)

   end subroutine read_setup

   subroutine read_components(file, def_components, scn, scn_plus)
      !! This subroutine reads the component names from input file

      character(len=*), intent(in) :: file 
      !! File name of fluid to characterize
      character(len=15), allocatable, intent(out) :: def_components(:)
      !! Names of defined components
      integer, allocatable, intent(out) :: scn(:) !! names of scn fractions
      character(len=15), intent(out) :: scn_plus !! name of residual fraction
      integer :: def_comp_nc
      !! Number of defined components being considered in the oil
      integer :: scn_nc !! number of single cuts being considered in the oil
      integer :: scn_nc_ps
      !! CN from which all SCN fractions will be lumped into the specified
      !! number of pseudos
      integer :: numbers_ps
      !! Number of pseudos in which the scn fractions grouped
      integer :: funit ! Environment variable

      namelist /nml_components/ def_components, scn, scn_plus

      call read_setup(file, def_comp_nc, scn_nc, scn_nc_ps, numbers_ps)
      allocate(def_components(def_comp_nc))
      allocate(scn(scn_nc))

      open(newunit=funit, file=file)
      read(funit, nml=nml_components)
      close(funit)

   end subroutine read_components

   subroutine read_composition(file, def_comp_z, scn_z, plus_z, plus6_z_exp, &
      plus7_z_exp, plus12_z_exp, plus30_z_exp)
      !! This subroutine reads the molar compositions of each component
      !! from input file

      character(len=*), intent(in) :: file 
      !! File name of fluid to characterize
      real(pr), allocatable, intent(out) :: def_comp_z(:)
      !! Set of corresponding mole fractions of defined components
      real(pr), allocatable, intent(out) :: scn_z(:)
      !! Set of corresponding mole fractions of scn cuts
      real(pr), intent(out):: plus_z
      !! Composition of residual fraction from input file, for C20+
      real(pr), intent(out):: plus6_z_exp
      !! Composition of residual fraction from input file, for C6+
      real(pr), intent(out):: plus7_z_exp
      !! Composition of residual fraction from input file, for C7+
      real(pr), intent(out):: plus12_z_exp
      !! Composition of residual fraction from input file, for C12+
      real(pr), intent(out):: plus30_z_exp
      !! Composition of residual fraction from input file, for C30+
      integer :: def_comp_nc
      !! Number of defined components being considered in the oil
      integer :: scn_nc !! number of single cuts being considered in the oil
      integer :: scn_nc_ps
      !! CN from which all SCN fractions will be lumped into the specified
      !! number of pseudos
      integer :: numbers_ps
      !! Number of pseudos in which the scn fractions grouped
      integer :: funit

      namelist /nml_composition/ def_comp_z, scn_z, plus_z, plus6_z_exp, &
         plus7_z_exp, plus12_z_exp, plus30_z_exp

      call read_setup(file, def_comp_nc, scn_nc, scn_nc_ps, numbers_ps)
      allocate(def_comp_z(def_comp_nc))
      allocate(scn_z(scn_nc))

      open(newunit=funit, file=file)
      read(funit, nml=nml_composition)
      close(funit)

   end subroutine read_composition

   subroutine read_molecular_weight(file, def_comp_mw, scn_mw, plus_mw, &
      plus6_mw_exp, plus7_mw_exp, plus12_mw_exp, plus30_mw_exp)
      !! This subroutine reads the molecular weights of each component from
      !! the input file

      character(len=*), intent(in) :: file !! File name
      real(pr), allocatable, intent(out) :: def_comp_mw(:)
      !! Set of corresponding molecular weights of defined components
      real(pr), allocatable, intent(out) :: scn_mw(:)
      !! Set of corresponding molecular weights of scn cuts
      real(pr), intent(out) :: plus_mw
      !!  Molecular weight of residual fraction, normally C20+
      real(pr), intent(out) :: plus6_mw_exp
      !! Molecular weight of the residual fraction, for C6+
      real(pr), intent(out) :: plus7_mw_exp
      !! Molecular weight of the residual fraction, for C7+
      real(pr), intent(out) :: plus12_mw_exp
      !! Molecular weight of the residual fraction, for C12+
      real(pr), intent(out) :: plus30_mw_exp
      !! Molecular weight of the residual fraction, for C30+
      integer :: def_comp_nc
      !! Number of defined components being considered in the oil
      integer :: scn_nc !! number of single cuts being considered in the oil
      integer :: scn_nc_ps
      !! CN from which all SCN fractions will be lumped into the specified
      !! number of pseudos
      integer :: numbers_ps 
      !! Number of pseudos in which the scn fractions grouped
      integer :: funit

      namelist /nml_molecular_weight/ def_comp_mw, scn_mw, plus_mw, &
         plus6_mw_exp, plus7_mw_exp, plus12_mw_exp, plus30_mw_exp

      call read_setup(file, def_comp_nc, scn_nc, scn_nc_ps, numbers_ps)
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
      real(pr), allocatable, intent(out) :: product_z_mw_def_comp(:)
      !! product between composition and molecular weight of defind components
      real(pr), allocatable, intent(out) :: product_z_mw_scn(:)
      !! product between composition and molecular weight of scn fractions
      real(pr), intent(out) :: sum_z_mw_i
      !! sum of the product between composition and molecular weight of the
      !! fluid's compounds
      real(pr), allocatable, intent(out) :: w(:)
      !! mass fractions of the fluid's compounds
      real(pr), intent(out) :: product_z_mw_plus
      !! product between composition and molecular weight of residual fraction
      integer :: def_comp_nc
      !! number of defined components being considered in the oil
      integer :: scn_nc
      !! number of single cuts being considered in the oil
      integer :: scn_nc_ps
      !! CN from which all SCN fractions will be lumped into the specified
      !!number of pseudos
      integer :: numbers_ps
      !! number of pseudos in which the scn fractions grouped
      integer :: funit
      real(pr), allocatable :: def_comp_z(:)
      !! set of corresponding mole fractions of defined components
      real(pr), allocatable :: scn_z(:)
      !! set of corresponding mole fractions of scn cuts
      real(pr) :: plus_z !! composition of residual fraction
      real(pr) :: plus6_z_exp
      !! composition of residual fraction from input file
      !! C6+ is considered the residual fraction
      real(pr) :: plus7_z_exp
      !! composition of residual fraction from input file
      !! C7+ is considered the residual fraction
      real(pr) :: plus12_z_exp
      !! composition of residual fraction from input file
      !! C12+ is considered the residual fraction
      real(pr) :: plus30_z_exp
      !! composition of residual fraction from input file
      !! C30+ is considered the residual fraction
      real(pr), allocatable :: def_comp_mw(:)
      !! set of corresponding molecular weights of defined components
      real(pr), allocatable  :: scn_mw(:)
      !! set of corresponding molecular weights of scn cuts
      real(pr) :: plus_mw !!  molecular weight of residual fraction
      real(pr) :: plus6_mw_exp
      !! molecular weight of residual fraction
      !! C6+ is considered the residual fraction
      real(pr) :: plus7_mw_exp
      !! molecular weight of residual fraction
      !! C7+ is considered the residual fraction
      real(pr) :: plus12_mw_exp
      !! molecular weight of residual fraction
      !! C12+ is considered the residual fraction
      real(pr) :: plus30_mw_exp
      !! molecular weight of residual fraction
      !! C30+ is considered the residual fraction
      real(pr), allocatable :: product_z_mw_i(:)
      ! product of mole fraction and molecular weight of each 
      ! component of the fluid
      real(pr), allocatable, intent(out) :: def_comp_w(:)
      !! mass fractions of defined components
      real(pr), allocatable, intent(out) :: scn_w(:)
      !! mass fractions of scn fractions
      real(pr), intent(out) :: plus_w
      !! mass fraction of residual fraction

      call read_setup(file, def_comp_nc, scn_nc, scn_nc_ps, numbers_ps)
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

      call read_composition(file, def_comp_z, scn_z, plus_z, plus6_z_exp, &
         plus7_z_exp, plus12_z_exp, plus30_z_exp &
         )
      call read_molecular_weight(file, def_comp_mw, scn_mw, plus_mw, &
         plus6_mw_exp, plus7_mw_exp, plus12_mw_exp, plus30_mw_exp &
         )

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

   subroutine read_density(file,scn_density,plus_density, plus6_density_exp, &
      plus7_density_exp, plus12_density_exp, plus30_density_exp &
      )
      !! This subroutine reads the density of each component from the input 
      !! file and calculated molar volume since C6 fraction.

      character(len=*), intent(in) :: file !! file name
      real(pr), allocatable, intent(out) :: scn_density(:) 
      !! Set of corresponding densities of scn cuts
      real(pr), intent(out):: plus_density 
      !! Density of residual fraction from input file
      real(pr), intent(out):: plus6_density_exp !! Experimental density for C6+
      real(pr), intent(out):: plus7_density_exp !! Experimental density for C7+
      real(pr), intent(out):: plus12_density_exp !! Experimental density for C12+
      real(pr), intent(out):: plus30_density_exp
      integer :: def_comp_nc 
      !! Number of defined components being considered in the oil
      integer :: scn_nc 
      !! Number of single cuts being considered in the oil
      integer :: scn_nc_ps 
      !! CN from which all SCN fractions will be lumped into the specified 
      !! number of pseudos
      integer :: numbers_ps 
      !! Number of pseudos in which the scn fractions grouped
      integer :: funit
      
      namelist /nml_density/ scn_density, plus_density, plus6_density_exp, &
         plus7_density_exp, plus12_density_exp, plus30_density_exp
      ! Namelist for reading density values from the input file

      call read_setup(file, def_comp_nc, scn_nc, scn_nc_ps, numbers_ps)
      allocate(scn_density(scn_nc))

      open(newunit=funit, file=file) 
      ! Open the input file and read densities using the namelist
      close(funit)

   end subroutine read_density

   type(FluidData) function data_from_file(file) result(data)
      !! This funtion allows to obtain experimental data from data imput
      character(len=*), intent(in) :: file !! file name

      data%filename = file

      call read_setup(file, data%def_comp_nc, data%scn_nc, data%scn_nc_ps, &
         data%numbers_ps &
         )
      call read_components(file, data%def_components, data%scn, data%scn_plus)
      call read_composition(file, data%def_comp_z, data%scn_z, data%plus_z, &
         data%plus6_z_exp, data%plus7_z_exp, data%plus12_z_exp, data%plus30_z_exp &
         )
      call read_molecular_weight(file, data%def_comp_mw, data%scn_mw, &
         data%plus_mw, data%plus6_mw_exp, data%plus7_mw_exp, data%plus12_mw_exp, &
         data%plus30_mw_exp &
         )
      call mass_fractions(file, data%w, data%product_z_mw_def_comp, &
         data%product_z_mw_scn, data%sum_z_mw_i, data%product_z_mw_plus, &
         data%def_comp_w, data%scn_w, data%plus_w &
         )
      call read_density(file, data%scn_density, data%plus_density, &
         data%plus6_density_exp, data%plus7_density_exp, data%plus12_density_exp, &
         data%plus30_density_exp &
         )

   end function data_from_file


end module data_from_input









