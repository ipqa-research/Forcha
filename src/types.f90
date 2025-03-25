module dtypes
   !! This module contains data type definitions and subroutines used for
   !! manipulatingfluid characterization. The main data types are FluidData
   !! and FluidDataOut, which are used to store and handle information about
   !! defined components, cut fractions, and other parameters of fluid.

   use constants, only: pr

   implicit none
   type :: FluidData
      !! Structure to save data from input file
      integer :: def_comp_nc
      !! Number of defined components being considered in the oil
      integer :: scn_nc
      !! Number of single cuts being considered in the oil
      integer :: scn_nc_ps
      !! CN from which all SCN fractions will be lumped into the specified
      !! number of pseudos
      integer :: numbers_ps
      !! number of pseudos in which the scn fractions grouped.
      character(len=:), allocatable :: filename !! Name of input file
      integer, allocatable :: scn (:)
      !! set of singles cuts being considered in the oil
      character(len=15), allocatable :: def_components (:)
      !! set of defined components being considered in the oil
      character(len=15) :: scn_plus !! name of residual fraction
      real(pr), allocatable :: def_comp_z(:)
      !! set of corresponding mole fractions of defined components
      real(pr), allocatable :: scn_z(:)
      !!  set of corresponding mole fractions of scn cuts
      real(pr), allocatable :: def_comp_mw(:)
      !! set of corresponding molecular weights of defined components
      real(pr), allocatable :: scn_mw(:)
      !! set of corresponding molecular weights of scn cuts
      real(pr), allocatable:: product_z_mw_def_comp(:)
      !! product between composition and molecular weight of defind components
      real(pr), allocatable:: product_z_mw_scn(:)
      !! product between composition and molecular weight of scn fractions
      real(pr) :: sum_z_mw_i
      !!  sum of the product between composition and molecular weight of the
      !!  fluid's compounds
      real(pr), allocatable :: w(:) !! mass fractions of the fluid's compounds
      real(pr) :: plus_z !! composition of residual fraction
      real(pr) :: plus6_z_exp !! composition of residual fraction, for C6+
      real(pr) :: plus7_z_exp !! composition of residual fraction, for C7+
      real(pr) :: plus12_z_exp !! composition of residual fraction, for C12+
      real(pr) :: plus30_z_exp !! composition of residual fraction, for C30+
      real(pr) :: plus_mw !!  molecular weight of residual fraction
      real(pr) :: plus6_mw_exp !! molecular weight of residual fraction, for C6+
      real(pr) :: plus7_mw_exp !! molecular weight of residual fraction, for C7+
      real(pr) :: plus12_mw_exp !! molecular weight of residual fraction, for C12+
      real(pr) :: plus30_mw_exp !! molecular weight of residual fraction, for C30+
      real(pr) :: product_z_mw_plus
      !! product between composition and molecular weight of residual fraction
      real(pr), allocatable :: def_comp_w(:)
      !! mass fractions of the defined compounds
      real(pr), allocatable :: scn_w(:)
      !! mass fractions of the scn-s compounds
      real(pr) :: plus_w !! mass fractions of the plus fraction
      real(pr), allocatable :: scn_density(:)
      !! set of corresponding densities of scn cuts
      real(pr) :: plus_density !! experimental density of the plus fraction
      real(pr) :: plus6_density_exp
      !! experimental density of the plus fraction, for C6+
      real(pr) :: plus7_density_exp
      !! experimental density of the plus fraction, for C7+
      real(pr) :: plus12_density_exp
      !! experimental density of the plus fraction, for C12+
      real(pr) :: plus30_density_exp
      !! experimental density of the plus fraction, for C30+

   end type FluidData

   type :: FluidDataOut
      !! Structure to save output data
      type(FluidData) :: input_data !! Structure to save data from input file
      real(pr), allocatable :: scn_z(:)
      !!  set of corresponding mole fractions of scn cuts calculated
      real(pr), allocatable :: log_scn_z(:)
      !!  set of corresponding logarithm of  mole fractions of scn cuts calculated
      real(pr) :: plus_mw !! Molecular weight of residual fraction compute
      real(pr) :: plus_z !! Composition of residual fraction compute
      real(pr) :: C
      !! Compute C constants which is used in equation \[M = 84-characterization%C(i-6)\]
      real(pr), allocatable :: scn_mw(:)
      !! set of corresponding molecular weights of scn cuts compute
      real(pr), allocatable :: scn_zm (:)
      !! Product between composition and molecular weight compute
      real(pr) :: plus_zm
      !! Product between composition and molecular weight for residual fraction compute
      real(pr) :: a_blr !! A constant for best linear regression line.
      real(pr) :: b_blr !! B constant for best linear regression line.
      real(pr) :: a_60 !! A constant for Cmax60 line.
      real(pr) :: b_60 !! B constant for Cmax60 line.
      real(pr) :: a_lim !! A constant for limit feasible line.
      real(pr) :: b_lim !! B constant for limit feasible line.
      real(pr) :: a !!  Slope of the best feasible regression line.
      real(pr) :: b !!  Intercept of the best feasible regression line.
      real(pr) :: r2 !! Square correlation coefficient.
      integer :: n_init !! Minimum carbon number obtained from the best linear regression
      integer :: c_max_blr  !! Maximum carbon number for best linear regression distribution
      integer :: c_max_lim !! Maximum carbon number for limit line distribution
      integer :: c_max_60 !! Maximum carbon number for C60 line distribution
      integer :: nc_plus !! number of elements in the distribution of C20+
      integer :: c_max
      !! output CN at which plus_z is reached, as the summation of single z(i)
      !! from the best linear distribution (blr)
      real(pr), allocatable :: carbon_number_plus(:) !! Carbon numbers of plus fractions
      real(pr), allocatable :: plus_z_i(:) !! Composition of the residual fraction
      real(pr), allocatable  :: product_z_mw_plus_i(:)
      !! Product of composition and molecular weight of the residual fraction
      real(pr) :: a_d
      !! ad constant which is used in equation \[rho_i = ad*exp(-i/10) +bd\] for density.
      real(pr) :: b_d
      !! bd constant which is used in equation \[rho_i = ad*exp(-i/10) +bd\] for density.
      real(pr) :: volume_cal !! Calculated volume based on the current density function
      real(pr) :: volume_exp !! Experimental volume of the fluid based on input data
      integer :: last_C !! Last defined SCN component, typically 19
      integer :: i_last !! Number of plus components, typically 124
      integer :: last !! Total number of elements after lumping
      integer :: scn_nc_new !! Define new SCN count after lumping adjustment
      real(pr), allocatable :: scn_i(:) !! Adjusted number of single cut fractions
      real(pr), allocatable :: plus6_density(:) ! Density values for the 6+ fraction
      real(pr) :: plus_w !!  molecular weight of residual fraction
      real(pr) :: plus_density !!  density of residual fraction
      real(pr), allocatable :: mol_fraction(:) !! Adjusted molar fraction of the fluid
      real(pr), allocatable ::  lumped_z(:)  !! Adjusted molar fractions after lumping
      real(pr), allocatable ::  lumped_mw(:) !! Adjusted molecular weights  after lumping
      real(pr), allocatable ::  lumped_densities(:) !! Adjusted densities  after lumping
      real(pr), allocatable :: critical_temperature(:)
      !! Critical temperature compute
      real(pr), allocatable :: critical_pressure(:)
      !! Critical pressure compute
      real(pr), allocatable :: acentric_factor(:)
      !! Acentric factor compute
      real(pr), allocatable :: m_funtion (:)
      !! M funtion compute to calculate acentric factor

   contains
      private
      procedure, pass :: write => write_result
      generic, public :: write (FORMATTED) => write
   end type FluidDataOut

contains

   subroutine write_result(characterization,unit,iotype,v_list,iostat,iomsg)
      !! This subroutine writes the fluid characterization results to an output file.
      use ftools__io, only: str
      use defined_critical_parameters

      implicit none
      class(FluidDataOut), intent(in) :: characterization !! Structure to save results
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
      character(len=*), parameter :: str_fmt_3 = "(A70,/)"

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

         write(unit, fmt=str_fmt_3, iostat = iostat) "COMPOSITIONAL RESULTS"
         write(unit, fmt=str_fmt_3, iostat = iostat) ""
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

   subroutine general_result(characterization)
      !! This subroutine writes the fluid characterization results to an output file.
      
      implicit none
      class(FluidDataOut), intent(in) :: characterization !! Structure to save results

         write(*,*) "-----------------------------------------------------------&
         --------------------------------------------------------------"
         write(*,*) "Best feasible regresion parameters"
         write(*,*) "Init_BFR:", characterization%n_init
         write(*,*) "A:", characterization%a,"      ","B:", characterization%b
         write(*,*) "C:", characterization%C
         write(*,*) "MW+:", characterization%plus_mw
         write(*,*) "Cmax:", characterization%c_max
         write(*,*) "Plus_z:", characterization%plus_z
         write(*,*) ""
         write(*,*) "Density funtion parameters parameters"
         write(*,*) ""
         write(*,*) "Ad:", characterization%a_d,"      ","Bd:", characterization%b_d
         write(*,*) "-----------------------------------------------------------&
         --------------------------------------------------------------"
         write(*,*) ""

   end subroutine general_result

end module dtypes


