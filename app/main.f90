program main


   use data_from_input, only: data_from_file, FluidData 
   use routines
   !use constants, only: pr
   !use routines
   !use data

   implicit none
   logical :: start
   real(pr) :: abe, a_lim, a_60! output real variable. Slope of the best regression line.
   real(pr) :: bbe , b_lim, b_60 !! output real variable. Intercept of the best regression line.
   real(pr):: r2, half!! output real variable. Square correlation coefficient.
   integer :: ninit !! minimum carbon number obtained from the best linear regression
   integer :: c_max_blr, c_max_lim, c_max_60, c_max !
   real(pr), allocatable :: ylog(:)
   real(pr) :: C
   real(pr) :: difference
   real(pr) :: plus_mw !!  calculated molecular weight of  residual fraction
   real(pr) :: plus_z!! calculated composition of residual fraction
   real(pr), allocatable :: scn_z(:)!! calculated composition of residual fraction
   real(pr), allocatable :: log_scn_z(:) !! logarithm of corresponding mole fractions of scn cuts
  
   

   !integer :: i,j, file_unit_1, file_unit_2, file_unit_3, Nfluid, Ncut, Ndef, Nps, maxC
   !integer, parameter :: maxD=15, imax=48
   !integer, dimension(330) :: rCN
   !real(pr) :: rMWplus, DenPlus , watPlus, zMp, Z6p, sumV , C, rMp, a, b, Mglobal
   !real(pr) :: aux, dold, dif, Cold,zp
   !real(pr), dimension(300):: zpi, zMpi
   !real(pr), dimension(maxD) :: zdef, rMdef
   !real(pr), dimension(imax) :: zcomp, rMW, Den, wat, zM,z, w, rn,rM
   !character(len=30) :: infile, outfile
   !character(len=7) :: Oil, Plus
   !character(len=4) :: DefComp(maxD)
   !logical :: start


   !integer, parameter :: n=5
   !real(pr), dimension(n):: array
   character(20) :: fluid
   integer :: i 
   type(FluidData) :: asd
   asd = data_from_file(file="oil1.nml")
   !allocate(ylog(asd%scn_nc))
   !print*, asd%def_comp_names
   !print*, asd%scn_z
   !ylog = log(asd%scn_z)
   
   !call Best_Linear_Regression(asd%scn_nc,asd%scn,ylog,asd%plus_z,abe,bbe,r2,ninit,c_max_blr)
   !call LimitLine(asd%scn_nc,asd%scn,asd%plus_z,abe,bbe,a_lim,b_lim,c_max_lim,half)
   !call  Line_C60_max (asd%scn_nc,asd%scn,asd%plus_z,abe,bbe,half,a_60,b_60,c_max_60)
   !print*, half
   
   !print*, "-------------------------------"
   !do i = 1, asd%def_comp_nc+asd%scn_nc+1
   !   print*, asd%w(i)
   !end do
   print*, "-------------------------------"
   !do i = 1, asd%scn_nc
   !   print*, asd%scn_z(i)
   !end do
   print*, "-------------------------------"
   !do i = 1, asd%scn_nc
   !   print*, asd%scn(i)
   !end do
   !print*, asd%plus_z

   print*, "-------------------------------"
   !allocate(scn_z(asd%scn_nc))
   Start = .true.
   C= 14
   fluid = "oil1.nml"
   call difference_M_from_C(fluid,"global_mw","calculate",start,plus_z,plus_mw,scn_z,C,c_max,difference)
   print*, difference

   !call data_in(Nfluid,Oil,Ncut,Ndef,DefComp,zdef,rMdef,rCN,zcomp,rMW,&
   !   Den,Plus, rMWplus, DenPlus,Nps,wat,watPlus,zM,zMp,Z6p,sumV,w,rn)

   !call get_C(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef)
   !print*, z
   !print*, zp
   !print*, C, rMp,zp
   ! add new conditions for C constant if this value is between 12 and 14
   !if (C > 14)then
   !   C = 14
   !   call GetNewMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef,rMW,rMdef, Mglobal)
   !else

  !    if (C < 12)then
  !       C = 12
  !       call GetNewMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef,rMW,rMdef, Mglobal)
  !    endif
  ! endif
   !print*, z(1:14)
  ! print*, C, rMp, zp

  ! print*,'----------------------------------------------------'


   !array = [12.0_pr , 12.5_pr , 13._pr  , 13.5_pr , 14._pr]

   !do i = 1, n
   !  C = array(i)
   !  call GetNewMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef,rMW,rMdef, Mglobal)

   !print*, C, rMp, Mglobal, zp, maxC
   !print*, '-----------------------------------------'

   !end do
   !allocate(log_scn_z(oil1%scn_nc))
   call select_method(fluid,"global_mw","calculate",C,log_scn_z,plus_z,plus_mw )
      
   print*, plus_z, C
end program main


