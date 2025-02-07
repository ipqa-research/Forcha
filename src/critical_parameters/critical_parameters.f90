module critical_parameters
   
   use constants
   use dtypes
   use defined_critical_parameters

contains

   subroutine get_critical_constants(fluid, characterization, eos)
      !! this subroutine doing:
      !! - Calculation of Tc, Pc, omega from Pedersen correlations (adapted from code "C7plusPedersenTcPcOm")
      !! - Calculation of EoS parameters (taken from CubicParam)
      !use critical_parameters

      implicit none
      type(FluidData), intent(in) :: fluid
      type(FluidDataOut), intent(inout) :: characterization
      character(len=*), intent(in) :: eos
      real(pr), allocatable :: modified_c(:)
      real(pr), allocatable :: tc(:)
      real(pr), allocatable :: pc(:)
      real(pr), allocatable :: om(:)
      real(pr), allocatable :: m_funtion(:)
      integer :: i

      associate(&
         rho => characterization%lumped_densities, &
         mw  => characterization%lumped_mw &
         )

         allocate(modified_c(0))
         allocate(tc(0))
         allocate(pc(0))
         allocate(om(0))
         allocate(m_funtion(0))

         call get_parameteres_for_critical(eos)

         tc = c1*rho + c2*log(mw) + c3*mw + (c4/mw)
         pc = exp(d1_p + d2*(rho**(d5)) + (d3/mw) + (d4/(mw**2)))

         if (eos == "SRK") then
            m_funtion = e1 + e2*mw + e3*rho + e4*(mw**(2))
         else
            if (eos == "PR" .or. eos =="RKPR") then
               m_funtion = e3*rho + 0.2 + 1.8*(1-exp(-mw/220))
               ! Alternative expression without max with MW - May 2019
            endif
         endif

         modified_c = c0 - m_funtion
         om = (-b + sqrt(b**2 - 4*a*modified_c))/(2*a)

         !do i = 1, size(mw)
         !   print*, mw(i), rho(i), tc(i), pc(i), om(i), m_funtion(i)
         !end do

         characterization%critical_temperature = tc
         characterization%critical_pressure = pc
         characterization%acentric_factor = om
         characterization%m_funtion = m_funtion

      end associate

   end subroutine get_critical_constants

end module critical_parameters





