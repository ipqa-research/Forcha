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
      type(FluidData) :: input_data
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
      integer :: last
      integer :: scn_nc_new
      real(pr), allocatable :: scn_i(:)
      real(pr), allocatable :: plus6_density(:)
      real(pr) :: plus_w
      real(pr) :: plus_density
      real(pr), allocatable :: mol_fraction(:)
      real(pr), allocatable ::  lumped_z(:)
      real(pr), allocatable ::  lumped_mw(:)
      real(pr), allocatable ::  lumped_densities(:)
      real(pr), allocatable :: critical_temperature(:)
      real(pr), allocatable :: critical_pressure(:)
      real(pr), allocatable :: acentric_factor(:)
      real(pr), allocatable :: m_funtion (:)

   contains
      private
      procedure, pass :: write => write_result
      generic, public :: write (FORMATTED) => write
   end type FluidDataOut

contains

   subroutine write_result(characterization,unit,iotype,v_list,iostat,iomsg)
      use ftools__io, only: str
      use critical_parameters

      implicit none
      class(FluidDataOut), intent(in) :: characterization
      integer :: i, i_prev
      integer, intent(in) :: unit
      integer, intent(out) :: iostat
      character(*), optional, intent(in) :: iotype
      character(*), optional, intent(inout) :: iomsg
      integer, optional, intent(in)  :: v_list(:)
      character(len=*), parameter :: str_fmt_1 = "(A4,7(A15,2x),/)"
      character(len=*), parameter :: num_fmt_1 = "(A3,1x,5(E15.5,2x),/)"
      character(len=*), parameter :: num_fmt_2 = "(A3,1x,7(E15.5,2x),/)"
      character(len=*), parameter :: str_fmt_2 = "(130('-'),/),/)"

      associate(&
         def_nc => characterization%input_data%def_comp_nc, &
         def_name => characterization%input_data%def_components, &
         def_mw => characterization%input_data%def_comp_mw, &
         z => characterization%mol_fraction, &
         mw => characterization%lumped_mw, &
         rho => characterization%lumped_densities, &
         scn => characterization%input_data%scn, &
         scn_nc_new => characterization%scn_nc_new, &
         tc => characterization%critical_temperature, &
         pc => characterization%critical_pressure, &
         omega => characterization%acentric_factor, &
         m_funtion => characterization%m_funtion &
         )



         write(unit, fmt=str_fmt_2, iostat = iostat)                                                                        
         !write(*,"(A64,/)") "----------------------------------------------------------------" 
         !write(*,*) "Best feasible regresion parameters"
         !write(unit, fmt=str_fmt_1, iostat = iostat) ""
         !write(unit, fmt=num_fmt_2, iostat = iostat) "Init_BFR:", characterization%n_init
         !write(*,*) "A:", characterization%a,"      ","B:", characterization%b
         !write(*,*) "C:", characterization%C
         !write(*,*) "MW+:", characterization%plus_mw
         !write(*,*) "Cmax:", characterization%c_max
         !write(*,*) ""
         !write(*,*) "Density funtion parameters parameters"
         !write(*,*) ""
         !write(*,*) "ad:", characterization%a_d,"      ","bd:", characterization%b_d
         !write(*,*) "----------------------------------------------------------------"
         !write(*,*) ""
         !write(*,*) "compositional result"
         !write(*,*) ""







         write(unit, fmt=str_fmt_1, iostat = iostat) "Comp", "Z", "Mw",  "Tc", &
            "Pc", "Omega", "Rho", "M"
         ! defined components
         do i = 1, def_nc
            write(unit, fmt=num_fmt_1, iostat = iostat) def_name(i), z(i), &
               def_mw(i), tc_def(i), pc_def(i), om_def(i)
         end do
         !
         i_prev = scn(1) - 1
         do i = 1, scn_nc_new
            write(unit, fmt=num_fmt_2, iostat = iostat)  adjustl('C' // str(i_prev + i)), &
               z(i+def_nc), mw(i), tc(i), pc(i), omega(i), rho(i), m_funtion(i)
         end do
         !
         do i = 1 + scn_nc_new, size(mw)
            write(unit, fmt=num_fmt_2, iostat = iostat)  'ps' // str(i - scn_nc_new), &
               z(i+def_nc), mw(i), rho(i), tc(i), pc(i), omega(i), m_funtion(i)
         end do

      end associate
   end subroutine write_result


end module dtypes


