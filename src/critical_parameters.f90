module critical_parameters
   use constants
   use dtypes, only: FluidData

   implicit none
   
   real(pr), allocatable :: tc_def(:), pc_def(:), om_def(:)
   real(pr):: c1, c2, c3, c4, d1_p, d2, d3, d4, d5, e1, e2, e3, e4
   real(pr):: a, b, c0, del_1

contains

   subroutine get_parameteres_for_critical(eos)
      implicit none
      character(len=*) :: eos

      tc_def = [126.2, 304.21, 190.564, 305.32, 369.83, 408.14, 425.12, 460.43, 469.7]
      pc_def = [34.0, 73.83, 45.99, 48.72, 42.48, 36.48, 37.96, 33.81, 33.7]
      om_def = [0.038, 0.224, 0.012, 0.099, 0.152, 0.181, 0.20, 0.228, 0.252]

      select case(eos) 

         case("SRK")
            c1 = 1.6312d2
            c2 = 8.6052d1
            c3 = 4.3475d-1
            c4 = -1.8774d3
            d1_p = -1.3408d-1
            d2 = 2.5019
            d3 = 2.0846d2
            d4 = -3.9872d3
            d5 = 1.0d0
            e1 = 7.4310d-1
            e2 = 4.8122d-3
            e3 = 9.6707d-3
            e4 = -3.7184d-6
            a = -0.176
            b = 1.574
            c0 = 0.48
            !del_1 = 1.0D0 

         case("PR") ! Coeficientes originales Tabla 5.3 Pedersen
            c1 = 7.34043d1   
            c2 = 9.73562d1    
            c3 = 6.18744d-1
            c4 = -2.05932d3
            d1_p = 7.28462d-2
            d2 = 2.18811d0
            d3 = 1.63910d2
            d4 = -4.04323d3
            d5 = 0.25d0
            e1 = 3.73765d-1
            e2 = 5.49269d-3
            e3 = 1.17934d-2
            e4 = -4.93049d-6
            a = -0.26992
            b = 1.54226
            c0 = 0.37464
            !del_1 = 1.0D0 + sqrt(2.0)

         case("RKPR")
            c1 = 7.34043d1  ! Coeficientes originales Tabla 5.3 Pedersen
            c2 = 9.73562d1    
            c3 = 6.18744d-1
            c4 = -2.05932d3
            d1_p = 2.3d0
            d2 = -0.5d0
            d3 = 1.85d2
            d4 = -4.0d3
            d5 = 0.25d0
            e1 = 3.73765d-1
            e2 = 5.49269d-3
            e3 = 1.17934d-2
            e4 = -4.93049d-6
            a = -0.26992
            b = 1.54226
            c0 = 0.37464
            
      end select
   
   end subroutine get_parameteres_for_critical
   
   end module





