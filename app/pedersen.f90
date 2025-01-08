module my_objective1
   use ForTimize, only: pr
   

   type :: ExpData
      real(pr), allocatable :: x(:)
      real(pr), allocatable :: zplus(:)
      real(pr), allocatable :: mplus(:)
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

      real(pr) :: a, b
      integer :: i
      integer, parameter :: n = 181
      real(pr) :: cn(n)
      real(pr) :: z_i(n)
      real(pr) :: m_i(n)
      real(pr) :: zplus_cal
      real(pr) :: mplus_cal


      do i = 1, n
         cn(i) = 20.0d0 + i - 1
      end do

      a = x(1)
      b = x(2)
      
      z_i = exp(a*cn+b)
      m_i = 14*cn-4
      zplus_cal = sum(z_i)
      mplus_cal = sum(z_i*m_i)/zplus_cal

      F = 0

      select type(data)
       type is (ExpData)
         F = sum((data%zplus - (zplus_cal))**2 + (data%mplus - (mplus_cal))**2)
      end select
   end subroutine foo

end module my_objective1


program oscar
   use ForTimize, only: pr, minimize
   use my_objective1, only: foo, ExpData

   real(pr) :: x(2), F
   type(ExpData) :: exp

   exp%zplus = [0.1107]
   exp%mplus = [419.3]

   ! Initial guess
   x = [   -9.3214408744470914E-002,  -2.7557414112453325  ]

   ! Minimize uses the Nelder-Mead algorithm as a default
   call minimize(foo, x, F, data=exp)

   ! Print results
   print *, x
   print *, F
   !print*, cn
end program oscar
