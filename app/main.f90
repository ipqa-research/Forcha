program main
    
  use precision_mod, only: pr
  use routines
  use data
  
  implicit none 

  integer :: i,j, file_unit_1, file_unit_2, file_unit_3, Nfluid, Ncut, Ndef, Nps, maxC
  integer, parameter :: maxD=15, imax=48 
  integer, dimension(330) :: rCN
  real(pr) :: rMWplus, DenPlus , watPlus, zMp, Z6p, sumV , C, rMp, a, b, Mglobal
  real(pr) :: aux, dold, dif, Cold,zp
  real(pr), dimension(300):: zpi, zMpi
  real(pr), dimension(maxD) :: zdef, rMdef
  real(pr), dimension(imax) :: zcomp, rMW, Den, wat, zM,z, w, rn,rM
  character(len=30) :: infile, outfile
  character(len=7) :: Oil, Plus
  character(len=4) :: DefComp(maxD)
  logical :: start
  
  
  
  integer, parameter :: n=5
  real(pr), dimension(n):: array
    
  call data_in(Nfluid,Oil,Ncut,Ndef,DefComp,zdef,rMdef,rCN,zcomp,rMW,&
  Den,Plus, rMWplus, DenPlus,Nps,wat,watPlus,zM,zMp,Z6p,sumV,w,rn)
  
  call get_C(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef)
  !print*, z
  !print*, zp
  !print*, C, rMp,zp
  ! add new conditions for C constant if this value is between 12 and 14
  if (C > 14)then
    C = 14
   call GetNewMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef,rMW,rMdef, Mglobal)
  else

    if (C < 12)then
        C = 12
        call GetNewMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef,rMW,rMdef, Mglobal)
    endif
  endif
  !print*, z(1:14)
  print*, C, rMp, zp
  
  print*,'----------------------------------------------------'

  
  array = [12.0_pr , 12.5_pr , 13._pr  , 13.5_pr , 14._pr]

  do i = 1, n
    C = array(i)
    call GetNewMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef,rMW,rMdef, Mglobal)
    
    print*, C, rMp, Mglobal, zp, maxC
    !print*, '-----------------------------------------'
    
  end do


end program main


