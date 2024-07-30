module routines
   use constants
   use dtypes, only: FluidData, FluidDataOut
   use data_from_input, only: data_from_file, FluidData

contains

   subroutine Linear_Regression(x, y, a, b, r2)
      !! This subroutine computes the regression line for a data set of x, y variables.
      implicit none

      integer ::  i
      integer ::  n !! total number of data set.
      real(pr),  intent(in) :: x(:) !! x: input array of length n which contains the set of independent variable
      real(pr),  intent(in) :: y(:) !! y: input array of length n which contains the set of dependent variable
      real(pr) :: t1,t2,t3,t4, aux1,aux2,aux3,aux4,aux5 ! internal variables
      real(pr), intent(out) :: a !! a: output real variable. Slope of the regression line
      real(pr), intent(out) :: b !! b: output real variable. Intercept of the regression line
      real(pr), intent(out) :: r2 !! r2: output real variable. Square correlation coefficient

      n = size(x)
      !Calculation of a y b --> y=a*x+b
      t1=0.;t2=0.;t3=0.;t4=0.;

      do i=1,n
         t1=t1+x(i)*y(i)
         t2=t2+x(i)
         t3=t3+y(i)
         t4=t4+x(i)**2
      end do

      a=(n*t1-t2*t3)/(n*t4-t2**2)
      b=(t3-a*t2)/n
      !coefficient calculation of correlations r2
      aux1=0.;aux2=0.;aux3=0.;aux4=0.;aux5=0.;

      do i=1, n
         aux1= aux1 + x(i)*y(i)
         aux2= aux2 + x(i)
         aux3= aux3 + y(i)
         aux4= aux4 + x(i)**2
         aux5= aux5 + y(i)**2
      end do

      r2=(aux1-aux2*aux3/n)**2 /((aux4-aux2**2/n)*(aux5-aux3**2/n))

   end subroutine Linear_Regression

   subroutine Best_Linear_Regression(scn_nc, scn, scn_z, plus_z, a, b, r2, &
         n_init, c_max_blr)
      !! This subroutine calculates the best regression line for an fluid.
      implicit none

      integer, intent(in) :: scn_nc !! integer input variable set to the total number of single cuts being considered in the fluid
      integer, intent(in) :: scn(:) !! set of singles cuts being considered in the fluid
      integer, intent(out) :: c_max_blr !! output CN at which plus_z is reached, as the summation of single z(i) from the best linear distribution (blr)
      real(pr), intent(in) :: scn_z(:) !! set of corresponding mole fractions of scn cuts
      real(pr), intent(in) :: plus_z !! composition of residual fraction from input file
      real(pr), intent(out) :: a !! output real variable. Slope of the best regression line.
      real(pr), intent(out) :: b !! output real variable. Intercept of the best regression line.
      real(pr), intent(out) :: r2 !! output real variable. Square correlation coefficient.
      integer, intent(out) :: n_init !! minimum carbon number obtained from the best linear regression
      integer :: i, j, k , k_old, n_best, x_aux ! internal variables
      real(pr), dimension(scn_nc) ::  x_blr, y_blr ! x_blr: carbon number vector, y_blr: logarithm of mole fraction vector; for each step
      real(pr) :: r2_old,r2_best, a_old, b_old, a_best, b_best, z_sum, z_aux ! auxiliary variables

      k=5
      r2=0.0001_pr
      r2_old=0.00001_pr
      r2_best=0.0001_pr

      do while (r2.gt.r2_old.or.r2_old.lt.0.9)
         k_old=k
         r2_old=r2
         a_old=a
         b_old=b

         if (r2.gt.r2_best) then
            r2_best=r2
            a_best=a
            b_best=b
            n_best=scn(scn_nc-k+2)
         end if

         if (k.gt.scn_nc) then
            if (r2.gt.r2_best)then
               n_init=scn(scn_nc-k+2)
               go to 22
            else
               r2=r2_best
               a=a_best
               b=b_best
               n_init=n_best
               go to 22
            end if
         end if

         j=1
         x_blr=0.
         y_blr=0.

         do i=scn_nc-k+1, scn_nc
            x_blr(j)=scn(i)
            y_blr(j)=scn_z(i)
            j=j+1
         end do

         call Linear_Regression(x_blr(:k), y_blr(:k), a, b, r2)
         k=k+1
      end do

      r2=r2_old
      a=a_old
      b=b_old
      n_init=scn(scn_nc-k_old+2)

      22    continue

      z_sum = 0d0
      x_aux=scn(scn_nc)

      do while (z_sum.lt.plus_z.and.x_aux<300)
         x_aux=x_aux+1
         z_aux= exp(a*x_aux+b)
         z_sum=z_sum+z_aux
      end do

      c_max_blr=x_aux

      !print*, a, b, c_max_blr, r2, n_init

   end subroutine Best_Linear_Regression

   subroutine LimitLine(scn_nc, scn, plus_z, a_blr, b_blr, a_lim, b_lim, & 
         c_max_lim, half)
      !! This subroutine obtains the limit line constants for an fluid.
      implicit none

      integer, intent(in) :: scn_nc !! integer input variable set to the total number of single cuts being considered in the fluid.
      integer, intent(in) :: scn(:) !! set of singles cuts being considered in the fluid.
      integer, intent(out) :: c_max_lim !! maximum carbon number obtained for limit line.
      !! output CN at which plus_z is reached, as the summation of single z(i) from the limit distribution (blr).
      real(pr), intent(in) :: a_blr  !! input constant from the Best linear regression.
      real(pr), intent(in) :: b_blr  !! input constant from the Best linear regression.
      real(pr), intent(in) :: plus_z !! composition of residual fraction from input file.
      real(pr), intent(out) :: a_lim !! output real variable. Slope of the limit line.
      real(pr), intent(out) :: b_lim !! output real variable. Intercept of the limit line.
      real(pr), intent(out) :: half ! variable defined in the numerical method to converge the area under the curve of the function.
      real(pr) :: z_lim, cross_cn, z_cross, z_aux 
      integer :: x_aux ! internal variable

      !A and B limit calculation.
      z_lim = 0.0_pr
      half = 0.506_pr ! modified by oscar, the variable half was removed from the argument
      
      do while (z_lim.lt.plus_z)
         half = half + 0.002d0
         cross_cn = scn(scn_nc)+half   ! typically 19.508 in first try
         z_cross = exp((a_blr*cross_cn) + b_blr)
         a_lim = -(z_cross)/plus_z
         b_lim = log(z_cross)-(a_lim*cross_cn)
         x_aux=scn(scn_nc)
         z_lim = 0._pr
         do while (z_lim.lt.plus_z.and.x_aux<300)
            x_aux=x_aux+1
            z_aux= exp(a_lim*x_aux+b_lim)
            z_lim=z_lim+z_aux
         end do
         continue
      end do

      c_max_lim = x_aux
      
   end subroutine LimitLine

   subroutine Line_C60_max (scn_nc, scn, plus_z, a_blr, b_blr, half, a_60, &
         b_60, c_max_60)
      !! This subroutine obtains the Cmax60 line constants for an fluid.

      implicit none
      integer, intent(in) :: scn_nc !! integer input variable set to the total number of single cuts being considered in the fluid
      integer, intent(in) :: scn(:) !! set of singles cuts being considered in the fluid
      integer, intent(out) :: c_max_60 !!output CN at which Zp is reached, as the summation of single z(i) from the Cmax60 distribution.
      real(pr), intent(in) :: a_blr  !! input constant from the Best linear regression
      real(pr), intent(in) :: b_blr  !! input constant from the Best linear regression
      real(pr), intent(in) :: plus_z !! composition of residual fraction from input file
      real(pr), intent(out) :: a_60 !! output real variable. Slope of the limit line
      real(pr), intent(out) :: b_60 !! output real variable. Intercept of the C60max line.
      real(pr), intent(in) :: half !variable defined in the numerical method to converge the area under the curve of the function
      integer ::  x_aux ! auxiliary variable
      real(pr) :: z_sum, cross_cn, z_cross, z_aux, F_tol, a_tol, var_range ! internal variables
      real(pr) :: a_old, Zp_60, F, dF_dA ! internal variables

      cross_cn = scn(scn_nc)+half   ! typically 19.51
      z_cross= exp((a_blr*cross_cn)+b_blr)
      F_tol=1_pr
      a_tol=1_pr
      a_60=a_blr
      var_range = 60.5_pr - cross_cn

      do while (a_tol.gt.1e-6_pr.or.F_tol.gt.1e-6_pr)
         a_old=a_60
         Zp_60 = (exp(a_60*var_range+log(z_cross))-z_cross)/a_60
         F = plus_z - Zp_60
         dF_dA = - (var_range*(a_60*Zp_60+z_cross)-Zp_60) / a_60 ! modified by oscar, error in the derivative in previous version is corrected
         a_60=a_old-(F/dF_dA)
         a_tol=abs(a_old-a_60)
         F_tol=abs(F)
      end do

      b_60=log(z_cross)-a_60*cross_cn
      z_sum=0.0_pr
      x_aux=scn(scn_nc)

      do while (z_sum.lt.plus_z)
         x_aux=x_aux+1
         z_aux= exp(a_60*x_aux+b_60)
         z_sum=z_sum+z_aux
      end do

      c_max_60=x_aux

      !print*, a_60,b_60,c_max_60
   end subroutine Line_C60_max

   subroutine select_method(fluid, mw_source, method, characterization)
      !! This subroutine defines the calculation method to perform the characterization,&
      !! based on the available experimental data. If the molecular weight data is experimental,
      !! define the mole fractions, molecular weights and densities directly from the input file. 
      !! On the other hand, if the molecular weights are assumed, they are recalculated according 
      !! to the methodology described by Martin et al and the molar fractions of the fluid are recalculated.

      implicit none
      type(FluidData), intent(inout) :: fluid !! Derived type to save inlet or experimental information of fluid to be characterized.
      type(FluidDataOut), intent(inout) :: characterization !! Derive type variable to save results of characterization methodology.
      character(len=*), intent(in), optional :: method 
      !! this option allows choose beetwen obtaining plus_mw or global_mw reported in the experimental information .
      character(len=*), intent(in) :: mw_source !! This variable indicates the source of the input data, that is, whether they are experimental or not.
      real(pr) :: sum_def_comp_z_plus !! compositions sum from (define component +1) to (plus component)
      real(pr), allocatable :: def_comp_moles(:) !!
      real(pr), allocatable :: scn_moles(:) !!
      real(pr) :: plus_moles
      real(pr) :: total_moles
      
      associate (&
         scn_nc => fluid%scn_nc, def_comp_nc => fluid%def_comp_nc, &
         def_comp_w => fluid%def_comp_w, def_comp_mw => fluid%def_comp_mw, &
         scn_w => fluid%scn_w, plus_w => fluid%plus_w, scn => fluid%scn, &  
         plus_mw => characterization%plus_mw, scn_z => characterization%scn_z, &
         plus_z => characterization%plus_z, log_scn_z => characterization%log_scn_z, & 
         C => characterization%C, scn_mw => characterization%scn_mw &
         )

      allocate(def_comp_moles(def_comp_nc))
      allocate(scn_moles(scn_nc))
      
      def_comp_moles = (def_comp_w)/(def_comp_mw)

      select case (mw_source)

       case("experimental")
         scn_mw = fluid%scn_mw 
         scn_moles = (scn_w)/(scn_mw)
         plus_moles = (plus_w)/(plus_mw)
         total_moles = sum(def_comp_moles)+sum(scn_moles)+(plus_moles)
         scn_z =  scn_moles / total_moles
         plus_z = plus_moles / total_moles
         
       case("calculated")
         
         if (method=="global_mw")then
            scn_mw = 84 + C*(scn-6)
            scn_z = fluid%product_z_mw_scn / scn_mw
            sum_def_comp_z_plus = sum(fluid%scn_z) + (fluid%plus_z)
            plus_z = sum_def_comp_z_plus - sum(scn_z) ! new molar fraction for residual cut, based on C
            plus_mw = fluid%product_z_mw_plus / plus_z
         end if

         if (method=="plus_mw")then
            scn_mw = 84 + C*(scn-6)
            scn_moles = (scn_w)/(scn_mw)
            plus_moles = (plus_w)/(plus_mw)
            total_moles = sum(def_comp_moles)+sum(scn_moles)+(plus_moles)
            scn_z =  scn_moles / total_moles
            plus_z = plus_moles / total_moles
            
         endif
      end select
      log_scn_z =  log(scn_z)
      end associate

   end subroutine select_method

   subroutine difference_mw_plus(fluid, mw_source, method, start, difference, &
         characterization)

      !! this subroutine ...
      implicit none

      type(FluidData), intent(inout) :: fluid
      type(FluidDataOut), intent(inout) :: characterization
      logical :: start
      character(len=*), intent(in), optional :: method
      character(len=*), intent(in) :: mw_source
      real(pr) :: plus_mw_cal !!  calculated molecular weight of  residual fraction
      real(pr), intent(out) :: difference 
      real(pr) :: half
      real(pr) :: r2_old,r2_best, a_old, b_old, a_best, b_best, z_sum, z_aux
      real(pr):: sum_z
      real(pr) :: denom 
      !! set to molecular weights calculated by \[M = 84-characterization%C(i-6)\]
      real(pr) :: sum_def_comp_z_plus 
      !! compositions sum from (define component +1) to (plus component)
      integer :: j, k , k_old, n_best, x_aux, i_0! internal variables
      integer :: i
      integer, dimension(300) :: carbon_number_plus
      real(pr), dimension(300) :: plus_z_i
      real(pr), dimension(300)  :: product_z_mw_plus_i

      associate (&
         scn_nc => fluid%scn_nc, scn => fluid%scn,&
         plus_z => characterization%plus_z, plus_mw => characterization%plus_mw, &
         log_scn_z => characterization%log_scn_z, C => characterization%C, &
         a_blr => characterization%a_blr, b_blr => characterization%b_blr, &
         a_60 => characterization%a_60, b_60 => characterization%b_60, &
         a_lim => characterization%a_lim, b_lim => characterization%b_lim, &
         r2 => characterization%r2, n_init => characterization%n_init, &
         c_max_blr => characterization%c_max_lim, a => characterization%a, &
         c_max_lim => characterization%c_max_lim, b => characterization%b, &
         c_max_60 => characterization%c_max_lim, &
         c_max => characterization%c_max &
         )
         
      1     i_0 = scn(scn_nc)

      call select_method(fluid, mw_source, method, characterization) ! select case
      call Best_Linear_Regression(scn_nc, scn, log_scn_z, &
            plus_z, a_blr, b_blr, r2, n_init, c_max_blr)
      call LimitLine(scn_nc, scn, plus_z, a_blr, b_blr,a_lim,b_lim,c_max_lim,half)
      call Line_C60_max (scn_nc,scn,plus_z,a_blr,b_blr,half,a_60,b_60,c_max_60) ! add line by oscar.
      
      if(a_blr < a_lim)then 
         !Added by oscar 05/12/2023.se elimina por que se agregan las restricciones 
         !para characterization%C>12 y characterization%C<12.
         a = a_lim
         b = b_lim
      else  ! this line was removed call Line_C60_max
         if(a_blr > a_60)then
            a = a_60
            b = b_60
         else
            a = a_blr
            b = b_blr
         end if
      end if

      sum_z = 0.0_pr
      i = 0
      do while (sum_z < plus_z .and. i< 300)
         i = i+1
         plus_z_i(i) = exp(a*(i+i_0)+b)
         sum_z = sum_z + plus_z_i(i)
         product_z_mw_plus_i(i) = plus_z_i(i)*(84+C*(i+i_0-6))
         carbon_number_plus(i) = i+i_0  
      end do

      if (i == 300.and.sum_z < plus_z)then
         if(start) C = C - 0.07
         if(.not.start) C = C - 0.01
         go to 1
      end if

      c_max = i+i_0
      plus_z_i(i) = plus_z_i(i) - (sum_z-plus_z) !   Adjustment to Zp (z20+)
      product_z_mw_plus_i(i) = plus_z_i(i)*(84+C*(i+i_0-6))

      if (mw_source == "experimental" .and. C>12 .and. C<14 ) then 
         ! mayor igual a 12 o menor igual 14.
         denom = sum((plus_z_i(1:i))*((carbon_number_plus(1:i)-6))) 
         ! denominator of equation to obtain C value directly
         C = (plus_z*(plus_mw-84))/denom
      endif

      plus_mw_cal = sum(product_z_mw_plus_i(1:i))/plus_z
      difference = plus_mw_cal-plus_mw

      end associate
      
      characterization%plus_z_i = plus_z_i(1:i)
      characterization%product_z_mw_plus_i =  product_z_mw_plus_i(1:i)
      characterization%carbon_number_plus =  carbon_number_plus(1:i)
      characterization%nc_plus = i


   end subroutine difference_mw_plus

   subroutine get_c_or_m_plus(fluid, mw_source, method, fix_C, characterization)
      !! This subroutine...
      implicit none
      type(FluidData) :: fluid
      type(FluidDataOut) :: characterization
      logical :: start
      character(len=*), intent(in) :: mw_source
      character(len=*), intent(in), optional :: method
      logical, intent(in) :: fix_C !! If C is < 12 or > 14 then fix to the
                                   !! closest value.
      !integer :: c_max !! output CN at which characterization%plus_z is reached, as the summation of single z(i) from the best linear distribution (blr)
      real(pr), allocatable :: log_scn_z(:) !! logarithm of corresponding mole fractions of scn cutsn
      real(pr) :: difference !!
      real(pr) :: difference_old, plus_mw_old
      real(pr) :: C_old
      real(pr) :: aux

      associate (C => characterization%C, plus_mw => characterization%plus_mw)

      select case (mw_source)
       case("experimental")
         C = 13.5
         Start = .true.
         plus_mw = fluid%plus_mw
         call difference_mw_plus(fluid, mw_source, method, start, difference, &
               characterization)
         
       case("calculated")

         if(method == 'global_mw') C=13
         if(method == 'plus_mw')   C=14

         Start = .true.
         plus_mw = fluid%plus_mw
         call difference_mw_plus(fluid, mw_source, method, start, & 
               difference_old, characterization)
         C_old = C
         if(method == 'global_mw') C = min(12.8, C-0.07)
         if(method == 'plus_mw') C = 13.5
         Start = .false.
         call difference_mw_plus(fluid, mw_source, method, start, difference, &
               characterization)

         do while (abs(difference) > 0.1)
            aux = C
            C = C - difference*(C-C_old)/(difference-difference_old)
            C_old = aux
            difference_old = difference
            call difference_mw_plus(fluid, mw_source, method, start, &
                  difference, characterization)
         end do
      end select

      if(fix_C .and. C>14 .or. C<12)then 
         
         
         if(C>14) C=14
         if(C<12) C=12

         Start = .true.
         plus_mw = fluid%plus_mw
         call difference_mw_plus(fluid, mw_source, "plus_mw" , start, &
               difference_old, characterization)
         plus_mw_old = plus_mw
         plus_mw = 0.9*plus_mw
         Start = .false.
         call difference_mw_plus(fluid, mw_source, "plus_mw" , start, & 
                  difference, characterization)

         do while (abs(difference)>0.00001)
            aux = plus_mw
            plus_mw = plus_mw - difference*(plus_mw-plus_mw_old)/(difference-difference_old)
            plus_mw_old = aux
            difference_old = difference
            call difference_mw_plus(fluid, mw_source, "plus_mw" , start, &
                  difference, characterization)
         end do  
      end if
      end associate
   end subroutine get_c_or_m_plus

   subroutine density_funtion(fluid, mw_source, characterization)
      !! this subroutine ...
      implicit none
      type(FluidData) :: fluid
      type(FluidDataOut) :: characterization
      character(len=*), intent(in) :: mw_source
      !real(pr), allocatable :: volume_6plus_cal(:)
      real(pr) :: volume_6plus_exp
      real(pr) :: a_d_old, aux
      real(pr) :: difference, difference_old


      volume_6plus_exp = sum((fluid%product_z_mw_scn)/(fluid%scn_density)) + &
            (fluid%product_z_mw_plus/fluid%plus_density)

      associate( &
         a_d => characterization%a_d , &
         b_d => characterization%b_d, &
         volume_6plus_cal => characterization%volume_6plus_cal &
         )
      ! Now find constants for the Density function
      a_d = -0.50 ! initial guess
      b_d = 0.685 - a_d*exp(-0.60)
      call calculate_volume_6plus(fluid, mw_source, characterization)
      difference_old = volume_6plus_cal - volume_6plus_exp
      a_d_old = a_d
      a_d = -0.49
      b_d = 0.685 - a_d*exp(-0.60)
      call calculate_volume_6plus(fluid, mw_source, characterization)
      difference = volume_6plus_cal - volume_6plus_exp
      
      do while (abs(difference) > 0.001)
         aux = a_d
         a_d = a_d - difference*(a_d - a_d_old)/(difference - difference_old)
         b_d = 0.685 - a_d*exp(-0.60)
         a_d_old = aux
         difference_old = difference
         call calculate_volume_6plus(fluid, mw_source, characterization)
         difference = volume_6plus_cal - volume_6plus_exp
     end do
      
     
      end associate
   end subroutine density_funtion

   subroutine calculate_volume_6plus(fluid, mw_source, characterization)
            !! this subroutine ...
      implicit none
      type(FluidData), intent(in) :: fluid
      type(FluidDataOut), intent(inout) :: characterization
      character(len=*), intent(in) :: mw_source
      real(pr) :: volume_6plus_cal
      real(pr), allocatable :: plus_density_cal(:), scn_density_cal(:), plus6_density(:)

      !allocate(density_plus_cal(0), density_scn_cal(0), plus6_density(0))


      associate(&
         a_d => characterization%a_d, &
         b_d => characterization%b_d, &
         cn_plus => characterization%carbon_number_plus, &
         z_m_plus_i => characterization%product_z_mw_plus_i, &
         scn_z_i => characterization%scn_z, &
         scn_mw_i => characterization%scn_mw &
         )

      scn_density_cal  = a_d*(exp(- real(fluid%scn, pr) )/10) + b_d
      plus_density_cal = a_d*(exp(- real(cn_plus, pr) )/10) + b_d
      plus6_density = [scn_density_cal, plus_density_cal]
   
      select case (mw_source)
         case("experimental")
            volume_6plus_cal = sum((scn_z_i*scn_mw_i)/(fluid%scn_density)) + &
            sum(z_m_plus_i/plus_density_cal)    
         case("calculated")
            volume_6plus_cal = sum((scn_z_i*scn_mw_i)/(scn_density_cal)) + &
            sum(z_m_plus_i/plus_density_cal)
      end select

      characterization%volume_6plus_cal = volume_6plus_cal

      end associate
   
   end subroutine calculate_volume_6plus


   type(FluidDataOut) function characterize(file, mw_source, method, fix_C) &
         result(characterization)

      type(FluidData) :: fluid
      character(len=*), intent(in) :: file !! file name
      character(len=*), intent(in), optional :: method !! plus_mw or global_mw
      character(len=*), intent(in) :: mw_source
      logical, intent(in), optional :: fix_C


      fluid = data_from_file(file)
      allocate(characterization%scn_z(fluid%scn_nc))
      allocate(characterization%log_scn_z(fluid%scn_nc))
      allocate(characterization%scn_mw(fluid%scn_nc))
      allocate(characterization%carbon_number_plus(0))
      allocate(characterization%plus_z_i(0))
      allocate(characterization%product_z_mw_plus_i(0))

      call get_c_or_m_plus(fluid=fluid, mw_source=mw_source, method=method, fix_C=fix_C, characterization= characterization)
      call density_funtion(fluid=fluid, mw_source=mw_source, characterization=characterization)

   end function characterize

end module routines


