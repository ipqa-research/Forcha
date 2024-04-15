module routines
   use constants
contains

   subroutine Linear_Regression(x,y,a,b,r2)
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

   subroutine Best_Linear_Regression(scn_nc,scn,scn_z,plus_z,a,b,r2,n_init,c_max_blr)
      !! This subroutine calculates the best regression line for an oil.
      implicit none

      integer, intent(in) :: scn_nc !! integer input variable set to the total number of single cuts being considered in the oil
      integer, allocatable, intent(in) :: scn(:) !! set of singles cuts being considered in the oil
      integer, intent(out) :: c_max_blr
      !! output CN at which plus_z is reached, as the summation of single z(i) from the best linear distribution (blr)
      real(pr), allocatable, intent(in) :: scn_z(:) !! set of corresponding mole fractions of scn cuts
      real(pr), intent(in) :: plus_z !! composition of residual fraction from input file
      real(pr), intent(out) :: a !! output real variable. Slope of the best regression line.
      real(pr), intent(out) :: b !! output real variable. Intercept of the best regression line.
      real(pr), intent(out) :: r2 !! output real variable. Square correlation coefficient.
      integer, intent(out) :: n_init !! minimum carbon number obtained from the best linear regression
      integer :: i, j, k , k_old, n_best, x_aux ! internal variables
      real(pr), dimension(scn_nc) ::  x_blr, y_blr
      real(pr) :: r2_old,r2_best, a_old, b_old, a_best, b_best, z_sum, z_aux

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

      print*, a, b, c_max_blr, r2, n_init

   end subroutine Best_Linear_Regression

   subroutine LimitLine(scn_nc,scn,plus_z,a_blr,b_blr,a_lim,b_lim,c_max_lim,half)
      !! This subroutine obtains the limit line constants for an oil.
      implicit none

      integer, intent(in) :: scn_nc !! integer input variable set to the total number of single cuts being considered in the oil
      integer, allocatable, intent(in) :: scn(:) !! set of singles cuts being considered in the oil
      integer, intent(out) :: c_max_lim
      !! output CN at which plus_z is reached, as the summation of single z(i) from the limit distribution (blr)
      real(pr), intent(in) :: a_blr  !! input constant from the Best linear regression
      real(pr), intent(in) :: b_blr  !! input constant from the Best linear regression
      real(pr), intent(in) :: plus_z !! composition of residual fraction from input file
      real(pr), intent(out) :: a_lim !! output real variable. Slope of the limit line
      real(pr), intent(out) :: b_lim !! output real variable. Intercept of the limit line.
      real(pr), intent(out) :: half
      real(pr) :: z_lim, cross_cn, z_cross, z_aux ! cn: carbon number
      integer :: x_aux

      !A and B limit calculation.
      z_lim = 0.0_pr
      half = 0.506_pr ! modified by oscar, the variable half was removed fron the argument
      print*, scn
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
      
      print*, a_lim,b_lim,c_max_lim
   end subroutine LimitLine


   subroutine Line_C60_max (scn_nc,scn,plus_z,a_blr,b_blr,half,a_60,b_60,c_max_60)
      !! This subroutine obtains the Cmax60 line constants for an oil

      implicit none
      integer, intent(in) :: scn_nc !! integer input variable set to the total number of single cuts being considered in the oil
      integer, allocatable, intent(in) :: scn(:) !! set of singles cuts being considered in the oil
      integer, intent(out) :: c_max_60 !!output CN at which Zp is reached, as the summation of single z(i) from the Cmax60 distribution.
      real(pr), intent(in) :: a_blr  !! input constant from the Best linear regression
      real(pr), intent(in) :: b_blr  !! input constant from the Best linear regression
      real(pr), intent(in) :: plus_z !! composition of residual fraction from input file
      real(pr), intent(out) :: a_60 !! output real variable. Slope of the limit line
      real(pr), intent(out) :: b_60 !! output real variable. Intercept of the C60max line.
      real(pr), intent(in) :: half
      integer ::  x_aux
      real(pr) :: z_sum, cross_cn, z_cross, z_aux, F_tol, a_tol, var_range
      real(pr) :: a_old, Zp_60, F, dF_dA

      
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

      print*, a_60,b_60,c_max_60

   end subroutine Line_C60_max

end module routines

