module data_from_input
  !! This module reads the fluid information from the corresponding input file.
  use dtypes, only: FluidData
  use constants, only: pr
  implicit none
contains

  subroutine read_setup(file, def_nc, scn_nc)
     !! Reads the setup data from input file
     character(len=*), intent(in) :: file !! file name
     integer, intent(out):: def_nc !! number of defined components before scn
     integer, intent(out):: scn_nc !! number of single carbon number components
     integer :: funit
     namelist /nml_setup/ def_nc, scn_nc

     open(newunit=funit, file=file)
     read(funit, nml=nml_setup)
     close(funit)
  end subroutine read_setup

  subroutine read_component_name(file, def_comp_names, scn_names, plus_name)
     !! Reads the component names data from input file
     character(len=*), intent(in) :: file !! file name
     character(len=15), allocatable, intent(out) :: def_comp_names(:) !! names of defined components before scn
     integer, allocatable, intent(out) :: scn_names(:) !! names of scn fractions
     character(len=15), intent(out) :: plus_name !! name of plus fraction
     integer :: def_nc !! number of defined components before scn
     integer :: scn_nc !! number of single carbon number components
     integer :: funit
     namelist /nml_components/ def_comp_names, scn_names, plus_name

     call read_setup(file, def_nc, scn_nc)
     allocate(def_comp_names(def_nc))
     allocate(scn_names(scn_nc))

     open(newunit=funit, file=file)
     read(funit, nml=nml_components)
     close(funit)
  end subroutine read_component_name

  subroutine read_composition(file, def_comp_z, scn_z, plus_z)
     !! Reads the molar compositions of each component from input file
     character(len=*), intent(in) :: file !! file name
     real(pr), allocatable, intent(out) :: def_comp_z(:) !! compositions of defined components
     real(pr), allocatable, intent(out) :: scn_z(:) !! compositions of single carbon numbers
     real(pr), intent(out):: plus_z !! composition of residual fraction from input file
     integer :: def_nc !! number of defined components before scn
     integer :: scn_nc !! number of single carbon number components
     integer :: funit
     namelist /nml_composition/ def_comp_z, scn_z, plus_z

     call read_setup(file, def_nc, scn_nc)
     allocate(def_comp_z(def_nc))
     allocate(scn_z(scn_nc))

     open(newunit=funit, file=file)
     read(funit, nml=nml_composition)
     close(funit)
  end subroutine read_composition

  subroutine read_molecular_weight(file, def_comp_mw, scn_mw, plus_mw)
     !! Reads the molecular weights of each component from the input file
     character(len=*), intent(in) :: file !! file name
     real(pr), allocatable, intent(out) :: def_comp_mw(:) !! molecular weights of defined components
     real(pr), allocatable, intent(out) :: scn_mw(:) !! molecular weights of single carbon numbers
     real(pr), intent(out) :: plus_mw !!  molecular weight of residual fraction
     integer :: def_nc !! number of defined components before scn
     integer :: scn_nc !! number of single carbon number components
     integer :: funit
     namelist /nml_molecular_weight/ def_comp_mw, scn_mw, plus_mw
     
     call read_setup(file, def_nc, scn_nc)
     allocate(def_comp_mw(def_nc))
     allocate(scn_mw(scn_nc))

     open(newunit=funit, file=file)
     read(funit, nml=nml_molecular_weight)
     close(funit)
  end subroutine read_molecular_weight

  type(FluidData) function data_from_file(file) result(data)
     character(len=*), intent(in) :: file
     call read_setup(file, data%def_nc, data%scn_nc)
     call read_component_name(file, data%def_comp_names, data%scn_names, &
        data%plus_name)
     call read_composition(file, data%def_comp_z, data%scn_z, data%plus_z)
     call read_molecular_weight(file, data%def_comp_mw, data%scn_mw, &
        data%plus_mw)
  end function data_from_file

end module data_from_input









