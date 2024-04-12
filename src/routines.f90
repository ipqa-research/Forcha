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

   subroutine Best_Linear_Regression(scn_nc,scn,scn_z,plus_z,a,b,r2,ninit)
      !! This subroutine calculates the best regression line for an oil.
      implicit none

      integer, intent(in) :: scn_nc !! integer input variable set to the total number of single cuts being considered in the oil
      integer, allocatable, intent(in) :: scn(:) !! set of singles cuts being considered in the oil
      integer :: i, j, k , kold, nbest, xaux,cbmax ! internal variables
      real(pr), allocatable, intent(in) :: scn_z(:) !! set of corresponding mole fractions of scn cuts
      real(pr), intent(in) :: plus_z
      real(pr), intent(out) :: a !! output real variable. Slope of the best regression line.
      real(pr), intent(out) :: b !! output real variable. Intercept of the best regression line.
      real(pr), intent(out) :: r2 !! output real variable. Square correlation coefficient.
      integer, intent(out) :: ninit !! minimum carbon number obtained from the best linear regression
      real(pr), dimension(scn_nc) ::  xBR, yBR
      real(pr) :: r2old,r2best, aold, bold, abest, bbest, zsum,zaux

      k=5
      r2=0.0001d0
      r2old=0.00001d0
      r2best=0.0001d0

      do while (r2.gt.r2old.or.r2old.lt.0.9)
         kold=k
         r2old=r2
         aold=a
         bold=b

         if (r2.gt.r2best) then
            r2best=r2
            abest=a
            bbest=b
            nbest=scn(scn_nc-k+2)
         end if

         if (k.gt.scn_nc) then
            if (r2.gt.r2best)then
               ninit=scn(scn_nc-k+2)
               go to 22
            else
               r2=r2best
               a=abest
               b=bbest
               ninit=nbest
               go to 22
            end if
         end if

         j=1
         xBR=0.
         yBR=0.

         do i=scn_nc-k+1, scn_nc
            xBR(j)=scn(i)
            yBR(j)=scn_z(i)
            j=j+1
         end do

         call Linear_Regression(xBR(:k), yBR(:k), a, b, r2)
         k=k+1
      end do

      r2=r2old
      a=aold
      b=bold
      ninit=scn(scn_nc-kold+2)

      22 continue

      zsum = 0d0
      xaux=scn(scn_nc)

      do while (zsum.lt.plus_z.and.xaux<300)
         xaux=xaux+1
         zaux= exp(a*xaux+b)
         zsum=zsum+zaux
      end do

      cbmax=xaux

      print*, a, b, cbmax, r2, Ninit

   end subroutine Best_Linear_Regression

end module routines


