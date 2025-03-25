module my_objective1
   !! This module used ForTimize package to optimize the parameters
   !! for pedersen distribution funtion
   use ForTimize, only: pr
   type :: ExpData
      real(pr), allocatable :: x(:) !! Vector of parameters to optimize
      real(pr), allocatable :: zplus(:) !! Experimental mole fraction of residual fraction
      real(pr), allocatable :: mplus(:) !! Experimental molecular weight of residual fraction
   end type ExpData

contains

   subroutine foo(x, F, dF, data)
      real(pr), intent(in) :: x(:)
      !! Vector of parameters to optimize
      real(pr), intent(out) :: F
      !! Value of the objetive function after optimization.
      real(pr), optional, intent(out) :: dF(:) !!
      !! Optional gradient.
      class(*), optional, intent(in out) :: data
      !! Special data that the function could use.

      real(pr) :: a, b !! parameters of pedersen distribution funtion 
      integer :: i !! Iteration variable
      integer, parameter :: n = 74 !! Parameter for defined maximum carbon number
      real(pr) :: cn(n) !! Vector of carbon numbers
      real(pr) :: z_i(n) !! Vector of compositions of residual fraction
      real(pr) :: m_i(n) !! Vector of molecular weights of residual fraction
      real(pr) :: zplus_cal !! Mole fraction of plus fraction compute
      real(pr) :: mplus_cal !! Molecular weight of plus fraction compute

      ! compute the carbon number vector
      do i = 1, n
         cn(i) = 7.0d0 + i - 1
      end do

      ! define parameters to optimizar into x vector 
      a = x(1)
      b = x(2)
      
      ! compute mole fraction using pedersen linear funtion
      z_i = exp(a*cn+b)
      ! compute molecular weighth using pedersen linear funtion
      m_i = 14*cn-4
      ! Sum the values to obtain the residual value
      zplus_cal = sum(z_i)
      mplus_cal = sum(z_i*m_i)/zplus_cal

      F = 0

      select type(data)
       type is (ExpData)
         F = sum((data%zplus - (zplus_cal))**2 + (data%mplus - (mplus_cal))**2)
      end select
   end subroutine foo

end module my_objective1


program pedersen_characterize
   !! This program optimizes the parameter of pedersen distribution funtion
   use ForTimize, only: pr, minimize
   use my_objective1, only: foo, ExpData

   real(pr) :: x(2), F
   type(ExpData) :: exp

   exp%zplus = [0.0545]
   exp%mplus = [158]

   ! Initial guess
   x = [   -9.3214408744470914E-002,  -2.7557414112453325  ]

   ! Minimize uses the Nelder-Mead algorithm as a default
   call minimize(foo, x, F, data=exp)

   ! Print results
   print *, x
   print *, F
   
end program pedersen_characterize