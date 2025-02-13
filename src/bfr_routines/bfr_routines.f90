module bfr_routines
   !! This module calculates the best possible linear regression for the
   !! distribution of mole fractions as a function of carbon number.
   use constants

contains

   subroutine Linear_Regression(x, y, a, b, r2)
      !! This subroutine computes the linear regression parameters
      !! (slope, intercept, and R²) for a given dataset of independent (x)
      !! and dependent (y) variables.
      implicit none
      ! Variable declarations
      integer :: i                ! Loop index
      integer :: n                ! Total number of data points
      real(pr), intent(in) :: x(:)
      !! Input array containing independent variable values (size n)
      real(pr), intent(in) :: y(:)
      !! Input array containing dependent variable values (size n)
      real(pr), intent(out) :: a
      !! Output variable: Slope of the regression line
      real(pr), intent(out) :: b
      !! Output variable: Intercept of the regression line
      real(pr), intent(out) :: r2
      !! Output variable: Coefficient of determination (R²)

      ! Internal variables for summations
      real(pr) :: t1, t2, t3, t4
      ! Summation terms for slope and intercept calculation
      real(pr) :: aux1, aux2, aux3, aux4, aux5
      ! Auxiliary variables for R² calculation

      ! Get the number of data points
      n = size(x)

      ! Initialize summation variables
      t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0

      ! Compute summations needed for slope (a) and intercept (b)
      do i = 1, n
         t1 = t1 + x(i) * y(i)
         t2 = t2 + x(i)
         t3 = t3 + y(i)
         t4 = t4 + x(i) ** 2
      end do

      ! Compute slope (a) and intercept (b) using least squares formula
      a = (n * t1 - t2 * t3) / (n * t4 - t2 ** 2)
      b = (t3 - a * t2) / n

      ! Initialize auxiliary variables for correlation coefficient (R²)
      aux1 = 0.0; aux2 = 0.0; aux3 = 0.0; aux4 = 0.0; aux5 = 0.0

      ! Compute summations needed for R² calculation
      do i = 1, n
         aux1 = aux1 + x(i) * y(i)
         aux2 = aux2 + x(i)
         aux3 = aux3 + y(i)
         aux4 = aux4 + x(i) ** 2
         aux5 = aux5 + y(i) ** 2
      end do

      ! Compute coefficient of determination (R²)
      r2 = (aux1 - aux2 * aux3 / n) ** 2 / &
         ((aux4 - aux2 ** 2 / n) * (aux5 - aux3 ** 2 / n))

   end subroutine Linear_Regression

   subroutine Best_Linear_Regression(scn_nc, scn, scn_z, plus_z, a, b, r2, &
      n_init, c_max_blr)
      !! This subroutine calculates the best linear regression for a
      !! fluid dataset. It determines the optimal regression line parameters
      !! and the carbon number at which plus_z is reached.

      implicit none

      ! Input variables
      integer, intent(in) :: scn_nc
      !! Total number of single cuts considered in the fluid
      integer, intent(in) :: scn(:)
      !! Array of single cuts considered in the fluid
      real(pr), intent(in) :: scn_z(:)
      !! Corresponding mole fractions of SCN cuts
      real(pr), intent(in) :: plus_z
      !! Composition of the residual fraction from input file

      ! Output variables
      real(pr), intent(out) :: a !! Slope of the best regression line
      real(pr), intent(out) :: b !! Intercept of the best regression line
      real(pr), intent(out) :: r2 !! Coefficient of determination (R²)
      integer, intent(out) :: n_init
      !! Minimum carbon number obtained from the best linear regression
      integer, intent(out) :: c_max_blr
      !! CN at which plus_z is reached, computed from the best linear
      !! distribution (BLR)

      ! Internal variables
      integer :: i, j, k, k_old, n_best, x_aux
      real(pr), dimension(scn_nc) :: x_blr, y_blr
      !! x_blr: Carbon number vector, y_blr: Logarithm of mole fraction vector
      real(pr) :: r2_old, r2_best, a_old, b_old, a_best, b_best, z_sum, z_aux
      !! Auxiliary variables

      ! Initialization
      k = 5 ! This value sets the carbon number from which the regressions start
      r2 = 0.0001_pr
      r2_old = 0.00001_pr
      r2_best = 0.0001_pr

      ! Iterative regression refinement
      do while (r2 > r2_old .or. r2_old < 0.9)
         k_old = k
         r2_old = r2
         a_old = a
         b_old = b

         if (r2 > r2_best) then
            r2_best = r2
            a_best = a
            b_best = b
            n_best = scn(scn_nc - k + 2)
         end if

         if (k > scn_nc) then
            if (r2 > r2_best) then
               n_init = scn(scn_nc - k + 2)
               go to 22
            else
               r2 = r2_best
               a = a_best
               b = b_best
               n_init = n_best
               go to 22
            end if
         end if

         j = 1
         x_blr = 0.0
         y_blr = 0.0

         do i = scn_nc - k + 1, scn_nc
            x_blr(j) = scn(i)
            y_blr(j) = scn_z(i)
            j = j + 1
         end do

         ! Perform linear regression on selected data
         call Linear_Regression(x_blr(:k), y_blr(:k), a, b, r2)
         k = k + 1
      end do

      ! Store best values
      r2 = r2_old
      a = a_old
      b = b_old
      n_init = scn(scn_nc - k_old + 2)

      22    continue

      ! Compute CN at which plus_z is reached
      z_sum = 0.0d0
      x_aux = scn(scn_nc)

      do while (z_sum < plus_z .and. x_aux < 300)
         x_aux = x_aux + 1
         z_aux = exp(a * x_aux + b)
         z_sum = z_sum + z_aux
      end do

      c_max_blr = x_aux

   end subroutine Best_Linear_Regression

   subroutine LimitLine(scn_nc, scn, plus_z, a_blr, b_blr, a_lim, b_lim, &
      c_max_lim, half)
      !! This subroutine obtains the limit line constants for a fluid dataset.
      !! It determines the parameters of the limit line and the carbon number
      !! at which plus_z is reached.
      implicit none

      ! Input variables
      integer, intent(in) :: scn_nc
      !! Total number of single cuts considered in the fluid
      integer, intent(in) :: scn(:)
      !! Array of single cuts considered in the fluid
      real(pr), intent(in) :: a_blr
      !! Slope from the best linear regression
      real(pr), intent(in) :: b_blr
      !! Intercept from the best linear regression
      real(pr), intent(in) :: plus_z
      !! Composition of the residual fraction from input file

      ! Output variables
      real(pr), intent(out) :: a_lim !! Slope of the limit line
      real(pr), intent(out) :: b_lim !! Intercept of the limit line
      integer, intent(out) :: c_max_lim
      !! Maximum carbon number obtained for the limit line
      real(pr), intent(out) :: half
      !! Parameter used in the numerical method to converge the area under the curve

      ! Internal variables
      real(pr) :: z_lim, cross_cn, z_cross, z_aux
      integer :: x_aux

      ! A and B limit calculation
      z_lim = 0.0_pr
      half = 0.506_pr  !! Initial half-value

      do while (z_lim < plus_z)
         half = half + 0.002d0
         cross_cn = scn(scn_nc) + half  !! Typically 19.508 in first iteration
         z_cross = exp((a_blr * cross_cn) + b_blr)
         a_lim = -(z_cross) / plus_z
         b_lim = log(z_cross) - (a_lim * cross_cn)
         x_aux = scn(scn_nc)
         z_lim = 0._pr

         do while (z_lim < plus_z .and. x_aux < 300)
            x_aux = x_aux + 1
            z_aux = exp(a_lim * x_aux + b_lim)
            z_lim = z_lim + z_aux
         end do
      end do

      c_max_lim = x_aux

   end subroutine LimitLine

   subroutine Line_C60_max(scn_nc, scn, plus_z, a_blr, b_blr, half, a_60, &
      b_60, c_max_60)
      !! This subroutine obtains the Cmax60 line constants for a fluid.
      implicit none

      ! Input variables
      integer, intent(in) :: scn_nc
      !! Total number of single cuts considered in the fluid
      integer, intent(in) :: scn(:)
      !! Array of single cuts considered in the fluid
      real(pr), intent(in) :: plus_z
      !! Composition of residual fraction from input file
      real(pr), intent(in) :: a_blr
      !! Slope from the best linear regression
      real(pr), intent(in) :: b_blr
      !! Intercept from the best linear regression
      real(pr), intent(in) :: half
      !! Numerical method variable to adjust the integration range

      ! Output variables
      real(pr), intent(out) :: a_60 !! Slope of the Cmax60 line
      real(pr), intent(out) :: b_60 !! Intercept of the Cmax60 line
      integer, intent(out) :: c_max_60
      !! Carbon number at which plus_z is reached from the Cmax60 distribution

      ! Internal variables
      integer :: x_aux
      real(pr) :: z_sum, cross_cn, z_cross, z_aux, F_tol, a_tol, var_range
      real(pr) :: a_old, Zp_60, F, dF_dA

      ! Compute initial values
      cross_cn = scn(scn_nc) + half
      z_cross = exp((a_blr * cross_cn) + b_blr)
      F_tol = 1_pr
      a_tol = 1_pr
      a_60 = a_blr
      var_range = 60.5_pr - cross_cn

      ! Newton-Raphson iteration for a_60
      do while (a_tol > 1e-6_pr .or. F_tol > 1e-6_pr)
         a_old = a_60
         Zp_60 = (exp(a_60 * var_range + log(z_cross)) - z_cross) / a_60
         F = plus_z - Zp_60
         dF_dA = - (var_range * (a_60 * Zp_60 + z_cross) - Zp_60) / a_60
         a_60 = a_old - (F / dF_dA)
         a_tol = abs(a_old - a_60)
         F_tol = abs(F)
      end do

      b_60 = log(z_cross) - a_60 * cross_cn
      z_sum = 0.0_pr
      x_aux = scn(scn_nc)

      ! Compute c_max_60
      do while (z_sum < plus_z)
         x_aux = x_aux + 1
         z_aux = exp(a_60 * x_aux + b_60)
         z_sum = z_sum + z_aux
      end do

      c_max_60 = x_aux

   end subroutine Line_C60_max


end module bfr_routines

