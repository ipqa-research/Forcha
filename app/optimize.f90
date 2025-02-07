module my_objective
   use ForTimize, only: pr

   type :: ExpData
      real(pr), allocatable :: x(:)
      real(pr), allocatable :: y(:)
   end type

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

       real(pr) :: a, b, c

       a = x(1)
       b = x(2)
       c = x(3)

       F = 0

       select type(data)
       type is (ExpData)
         F = sum( (data%y - (quad(data%x, a, b, c)))**2 )
         print *, F, data%x
       end select
   end subroutine

   elemental function quad(x, a, b, c)
      real(pr), intent(in) :: x, a, b, c
      real(pr) :: quad

      quad = a*x**2 + b*x + c
   end function
end module


program oscar
   use ForTimize, only: pr, minimize
   use my_objective, only: foo, ExpData

   real(pr) :: x(3), F
   type(ExpData) :: exp

   exp%x = [2, 5, 7]
   exp%y = [4, 22, 50]

   ! Initial guess
   x = [1, 4, 3]

   ! Minimize uses the Nelder-Mead algorithm as a default    
   call minimize(foo, x, F, data=exp)
   ! call foo(x, F, data=exp)
   
   ! Print results
   print *, x
   print *, F
end program oscar
