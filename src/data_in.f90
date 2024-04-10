module BE
    implicit none
    integer :: Ninit
end module BE

module MassFracAndMoles
    use constants, only: pr
    implicit none
    integer, parameter :: imax = 48
    real(pr), dimension(imax) :: w, rn
end module MassFracAndMoles


module data

    use constants, only: pr
    !use types, only: FluidData

    implicit none
    integer :: nv, nv1, nNcut, nNdef
    real(pr) :: Zplus
    

contains

    subroutine data_in (Nfluid,Oil,Ncut,Ndef,DefComp,zdef,rMdef,rCN,zcomp,rMW,&
        Den,Plus, rMWplus, DenPlus,Nps,wat,watPlus,zM,zMp,Z6p,sumV,w, rn )

        implicit none

        integer, parameter :: maxD=15, imax=48
        integer :: i,j!, nv1
        integer :: file_unit_1, file_unit_2, file_unit_3, Nfluid, Ncut, Ndef, Nps
        integer, dimension(330) :: rCN
        real(pr), dimension(maxD) :: zdef, rMdef
        real(pr), dimension(imax) :: zcomp, rMW, Den, wat, zM, zMdef, w, rn,MWt
        real(pr) :: rMWplus, DenPlus , watPlus, zMp, Z6p, sumV, zMtotal
        character(len=30) :: infile, outfile
        character(len=7) :: Oil, Plus
        character(len=4) :: DefComp(maxD)

        write(*,*) 'ARE MOLECULAR WEIGTHS FROM DATAIN EXPERIMENTAL?'
        write(*,*) '1: YES   2: NO'
        read(*,*) nv

        write(*,*) 'ENTER THE VERSION NUMBER'
        write(*,*) '1: C.M   2: CMplus'
        read(*,*) nv1

        !nv2 = nv1

        ! open input file
        write(*,*) 'ENTER INFILE'
        read(*,*) infile
        write(*, *) "READING INFILE: " // infile
        open(newunit=file_unit_1, file =infile)
        write(*,*) 'ENTER OUTFILE'
        read(*,*) outfile
        open(newunit=file_unit_2, file = outfile)
        open(newunit=file_unit_3, file = 'TfPsatIn.txt')
        ! read input file:read number fluid, single carbon number for each fluid.
        read(file_unit_1,*) Nfluid
        do j = 1, Nfluid
            read (file_unit_1,*) Oil
            read (file_unit_1,*) Ncut, Ndef
        end do

        nNcut = Ncut
        nNdef = Ndef


        ! read defined compounds with  molar fraction and M by fluid j and compound i.
        do i = 1, Ndef
            read(file_unit_1,*) DefComp(i), zdef(i), rMdef(i)
        end do
        ! read carbon number and molar fraction of fluid j and comp. i
        do i = 1, Ncut
            read(file_unit_1,*) rCN(i), zcomp(i), rMW(i),Den(i)
        end do

        read(file_unit_1,*) Plus, zPlus, rMWplus, DenPlus ! plus fraction properties
        read (file_unit_1,*) Nps !number of pseudos for the C20+ fraction
        read (file_unit_1,*) wat(1:Ndef+Ncut),watPlus ! atmospheric oil weight fractions

        !zpl  =  zPlus ! added by oscar
        ! There get Zi*Mi and volume from cut to plus fraction

        zMdef(1:Ndef) = zdef(1:Ndef)*rMdef(1:Ndef)
        zM(1:Ncut) = zcomp(1:Ncut) * rMW(1:Ncut)
        zMp = zPlus * rMWplus
        zMtotal = sum(zMdef(1:Ndef)) + sum(zM(1:Ncut)) + zMp
        w(1:Ndef) = zMdef(1:Ndef) / zMtotal
        w(Ndef+1:Ndef+Ncut) = zM(1:Ncut) / zMtotal
        w(Ndef+Ncut+1) = zMp / zMtotal
        rn(1:Ndef) = zdef(1:Ndef)/zMtotal
        Z6p = sum(zcomp(1:Ncut)) + zPlus
        sumV = sum(zM(1:Ncut)/Den(1:Ncut)) + zMp/DenPlus
        
        !agregado para calcular el peso molecular global
        !MWt(1:Ndef) = rMdef(1:Ndef) 
        !MWt(Ndef+1:Ndef+Ncut) = rMW(1:Ncut)
        
        !do i= 1, Ndef+Ncut
        !    print*, MWt(i)
        !end do

 


        !print*, Zplus
    end subroutine data_in

end module data




module routines

    use constants, only: pr
    use data
contains


    subroutine get_C(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef)

        implicit none

        integer, parameter :: imax=48, maxD=15
        logical :: start
        real(pr) :: C, Z6p, zMp, a, b, rMp, rMWplus
        real(pr) :: aux, dold, dif, Cold,zp
        real(pr), dimension(imax) :: zM  ! inputs coming from main program
        integer, dimension(imax) :: rCN ! inputs coming from main program
        real(pr), dimension(imax) :: z, rM,w,rn, zcomp
        integer :: Ncut, maxC, Ndef
        real(pr), dimension(300):: zpi, zMpi
        real(pr), dimension(maxD) :: zdef

        !write(*,*) 'ENTER THE VERSION NUMBER'
        !write(*,*) '1: CMOD   2: CMM'
        !read(*,*) nv


        if(nv1==1)then

            ! Find the C value for molecular weights line
            C = 13.0  ! initial guess
            Start = .true.
            rMp = rMWplus

            call difMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,dold,a,b,rMp,maxC,z,zpi,zMpi,zcomp,zdef,zp)
            Cold = C
            C = min(12.8,C-0.07)
            Start = .false.

            

            call difMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,dif,a,b,rMp,maxC,z,zpi,zMpi,zcomp,zdef,zp)

            do while (abs(dif) > 0.1)
                aux = C
                C = C - dif*(C-Cold)/(dif-dold)
                Cold = aux
                dold = dif
                call difMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,dif,a,b,rMp,maxC,z,zpi,zMpi,zcomp,zdef,zp)
            end do

            !print*, C,maxC, rMp
            

        else

            ! Find the C value for molecular weights line
            C = 14.0  ! initial guess  and maximum C allowed
            Start = .true.
            rMp = rMWplus

            call difMpfromCfull (Start,C,Ndef,Ncut,rCN,dold,a,b,rM,rMp,maxC,z,zp,zpi,zMpi,w,rn,zcomp,zdef)
            Cold = C
            C =  13.5 !min(12.8,C-0.07)
            Start = .false.

            call difMpfromCfull (Start,C,Ndef,Ncut,rCN,dif,a,b,rM,rMp,maxC,z,zp,zpi,zMpi,w,rn, zcomp,zdef)

            do while (abs(dif) > 0.1)
                aux = C
                C = C - dif*(C-Cold)/(dif-dold)
                Cold = aux
                dold = dif
                call difMpfromCfull (Start,C,Ndef,Ncut,rCN,dif,a,b,rM,rMp,maxC,z,zp,zpi,zMpi,w,rn,zcomp,zdef)
            end do

            !print*, C, maxC, rMp, zp

        end if

    end subroutine get_C

    subroutine difMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,dif,a,b,rMp,maxC,z,zpi,zMpi,zcomp,zdef,zp)

        use data

        implicit none

        integer, parameter :: imax=48,  maxD=15
        logical :: start
        real(pr) :: C, Z6p, zMp, dif, a, b, rMp, zp,aBE,bBE,r2, alim, blim, half
        real(pr) :: a60,b60, sumz, rMpcalc
        real(pr), dimension(imax) :: zM   ! inputs coming from main program
        integer, dimension(imax) :: rCN,x   ! inputs coming from main program
        real(pr), dimension(imax) :: ylog, rM, z, zcomp
        integer :: i, Ncut, maxC, i0, CmaxL, C60max
        real(pr), dimension(300):: zpi, zMpi
        real(pr), dimension(maxD) :: zdef

        1       i0 = rCN(Ncut)

        if(nv == 1)then !introduce 'if' to choise where experimental data from
            x = rCN
            ylog = log(zcomp)
            zp = 1 - sum(zdef(1:nNdef)) - sum(zcomp(1:nNcut))
            !zp = zpl
        end if


        if(nv == 2)then
            rM(1:Ncut) = 84 + C*(rCN(1:Ncut)-6)
            z(1:Ncut) = zM(1:Ncut) / rM(1:Ncut)
            Zp = Z6p - sum(z(1:Ncut)) ! new molar fraction for residual cut, based on C
            rMp = zMp / Zp
            x = rCN
            ylog = log(z)
        end if

        print*, zp

        call BestLinearRegression(Ncut,x,ylog,Zp,aBE,bBE,r2)
        call LimitLine(Ncut,rCN,Zp,aBE,bBE,half,alim,blim,CmaxL)
        call LineC60max (Ncut,Zp,rCN,aBE,bBE,half,a60,b60,C60max)

        if(aBE < alim)then ! added by oscar 05/12/2023.se elimina por que se agregan las restriccones para C>12 y C<12.
            a = alim
            b = blim
        else
            !call LineC60max (Ncut,Zp,rCN,aBE,bBE,half,a60,b60,C60max)
            if(aBE > a60)then
                a = a60
                b = b60
            else
                a = aBE
                b = bBE
            end if
        end if

        sumz = 0.0d0
        i = 0
        do while (sumz<Zp.and.i<300)
            i = i+1
            zpi(i) = exp(a*(i+i0)+b)
            sumz = sumz + zpi(i)
            zMpi(i) = zpi(i)*(84+C*(i+i0-6))
        end do

        if(i==300.and.sumz<Zp)then
            if(start) C = C - 0.07
            if(.not.start) C = C - 0.01
            go to 1
        end if
        !   Adjustment to Zp (z20+)
        zpi(i) = zpi(i) - (sumz-Zp)
        zMpi(i) = zpi(i)*(84+C*(i+i0-6))
        ! ===============================
        rMpcalc = sum(zMpi(1:i))/Zp
        dif = rMpcalc-rMp
        maxC = i+i0

    end subroutine difMpfromC
  

    subroutine difMpfromCfull(Start,C,Ndef,Ncut,rCN,dif,a,b,rM,rMp,maxC,z,zp,zpi,zMpi,w,rn,zcomp,zdef)

        !use MassFracAndMoles, only: w,rn
        use data
        implicit none

        integer, parameter :: imax=48, maxD=15
        logical :: start
        real(pr) :: C, dif, a, b, rMp, zp ,aBE ,bBE ,r2 , alim, blim, half
        real(pr) :: a60,b60, sumz, rMpcalc, rntot, aold
        integer, dimension(imax) :: rCN,x   ! inputs coming from main program
        real(pr), dimension(imax) :: ylog, rM, z, w, rn, zcomp, ztotal
        integer :: i,Ndef, Ncut, maxC, i0, CmaxL, C60max
        real(pr), dimension(300):: zpi, zMpi
        real(pr), dimension(maxD) :: zdef

        1       i0 = rCN(Ncut)

        if(nv == 1)then !introduce 'if' to choise where experimental data from
            x = rCN
            ylog = log(zcomp)
            zp = 1 - sum(zdef(1:nNdef)) - sum(zcomp(1:nNcut))
            !zp = zpl
             
        end if

        if(nv == 2)then
            rM(1:Ncut) = 84 + C*(rCN(1:Ncut)-6)
            rn(Ndef+1:Ndef+Ncut) = w(Ndef+1:Ndef+Ncut)/rM(1:Ncut)
            rn(Ndef+Ncut+1) = w(Ndef+Ncut+1)/rMp
            rntot = sum(rn(1:Ndef+Ncut+1))
            ztotal(1:Ndef+Ncut) = rn(1:Ndef+Ncut) / rntot
            z(1:Ncut) = rn(Ndef+1:Ndef+Ncut) / rntot
            Zp = rn(Ndef+Ncut+1) / rntot
            x = rCN
            ylog = log(z)
        end if
        !print*, zp
        !print*, z

        call BestLinearRegression(Ncut,x,ylog,Zp,aBE,bBE,r2)
        call LimitLine(Ncut,rCN,Zp,aBE,bBE,half,alim,blim,CmaxL)
        call LineC60max (Ncut,Zp,rCN,aBE,bBE,half,a60,b60,C60max)

        if(aBE < alim)then !, para C se excede el valor de 14 en la mayoria de los casos.
            a = alim
            b = blim
        else
            !call LineC60max (Ncut,Zp,rCN,aBE,bBE,half,a60,b60,C60max)
            if(aBE > a60)then
                a = a60
                b = b60
            else
                a = aBE
                b = bBE
            end if
        end if


        sumz = 0.0d0
        i = 0
        do while (sumz<Zp.and.i<300)
            i = i+1
            zpi(i) = exp(a*(i+i0)+b)
            sumz = sumz + zpi(i)
            zMpi(i) = zpi(i)*(84+C*(i+i0-6))
        end do
        !agregar condicion cuando C = 14 o C = 12 no olvidar consultar con martin
        if(i==300.and.sumz<Zp)then
            if(start) C = C - 0.07
            if(.not.start) C = C - 0.01
            go to 1
        end if

        !    do while (i<40.and.a>alim) eliminado por oscar 20/03/2024
        !        aold = a
        !        a = max(1.0001*a, alim)
        !        b = b - (a-aold)*19.51
        !        sumz = 0.0d0
        !        i = 0
        !        do while (sumz<Zp.and.i<300)
        !            i = i+1
        !            zpi(i) = exp(a*(i+i0)+b)
        !            sumz = sumz + zpi(i)
        !            zMpi(i) = zpi(i)*(84+C*(i+i0-6))
        !        end do
        !    end do

        !   Adjustment to Zp (z20+)
        zpi(i) = zpi(i) - (sumz-Zp)
        zMpi(i) = zpi(i)*(84+C*(i+i0-6))
        ! ===============================
        rMpcalc = sum(zMpi(1:i))/Zp
        dif = rMpcalc-rMp
        maxC = i+i0

    end subroutine difMpfromCfull

   
    subroutine BestLinearRegression(Ncut,rCN,zcomp,Zp,a,b,r2)
        !================================================================================================================
        !This subroutine calculates the best regression line for an oil.
        !Ncut: integer input variable set to the total number of single cuts being considered in the oil.
        !rCN: input array of length Ncuts which contains the set of carbon numbers.
        !zcomp: input array of length Ncuts which contains the set of corresponding mole fractions.
        !a: output real variable. Slope of the best regression line.
        !b: output real variable. Intercept of the best regression line.
        !r2: output real variable. Square correlation coefficient.
        !================================================================================================================
        use be, only: ninit

        implicit none

        integer:: k,kold, Nbest, j,i, Ncut, xaux,CBmax
        integer, dimension(Ncut) :: rCN
        real(pr), dimension(Ncut) :: zcomp, xBR, yBR
        real(pr) :: r2,r2old,r2best, aold, bold, abest, bbest,a,b, zsum,zaux, Zp

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
                Nbest=rCN(Ncut-k+2)
            end if

            if (k.gt.Ncut) then
                if (r2.gt.r2best)then
                    Ninit=rCN(Ncut-k+2)
                    go to 22
                else
                    r2=r2best
                    a=abest
                    b=bbest
                    Ninit=Nbest
                    go to 22
                end if
            end if

            j=1
            xBR=0.
            yBR=0.

            do i=Ncut-k+1, Ncut
                xBR(j)=rCN(i)
                yBR(j)=zcomp(i)
                j=j+1
            end do

            call LinearRegression(xBR(:k), yBR(:k), a, b, r2)
            k=k+1

        end do

        r2=r2old
        a=aold
        b=bold
        Ninit=rCN(Ncut-kold+2)

        22      continue

        zsum = 0d0
        xaux=rCN(Ncut)

        do while (zsum.lt.Zp.and.xaux<300)
            xaux=xaux+1
            zaux= exp(a*xaux+b)
            zsum=zsum+zaux
        end do

        CBmax=xaux

        !print*, a, b, CBmax, r2, Ninit

    end subroutine BestLinearRegression


    subroutine LinearRegression(x,y,a,b,r2)
        !================================================================================================================
        !This subroutine computes the regression line for a data set of x, y variables.
        !n: integer input variable set to the total number of data set.
        !x: input array of length n which contains the set of independent variable.
        !y: input array of length n which contains the set of dependent variable.
        !a: output real variable. Slope of the regression line.
        !b: output real variable. Intercept of the regression line.
        !r2: output real variable. Square correlation coefficient.
        !================================================================================================================
        implicit none

        integer ::  i, n
        real(pr), dimension(:) :: x,y
        real(pr) :: t1,t2,t3,t4, aux1,aux2,aux3,aux4,aux5,a,b,r2

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

    end subroutine LinearRegression


    subroutine LimitLine(Ncut,rCN,Zp,aBE,bBE,half,alim,blim,Cmax)
        !================================================================================================================
        !This subroutine obtains the limit line constants for an oil.
        !Ncut: integer input number of single cuts being considered in the oil.
        !rCN: input array which contains the set of carbon numbers.
        !Zp: input mole fractions of the residual fraction.
        !aBE & bBE: input constants from the Best Extrapolation
        !alim: output real variable. Slope of the limit line.
        !blim: output real variable. Intercept of the limit line.
        !Cmax: output CN at which Zp is reached, as the summation of single z(i) from the limit distribution.
        !================================================================================================================

        implicit none

        integer, parameter :: imax=48
        integer::  Ncut, xaux,Cmax
        integer, dimension(imax) :: rCN
        real(pr) :: Zp,aBE,bBE,half,alim,blim, zlim, crossCN, zcross, zaux

        !A and B limit calculation.
        zlim = 0.d0
        half = 0.506d0

        do while (zlim.lt.Zp)
            half = half + 0.002d0
            crossCN = rCN(Ncut)+half   ! typically 19.508 in first try
            zcross = exp(aBE*crossCN + bBE)
            alim = -zcross/Zp
            blim = log(zcross)-alim*crossCN
            xaux=rCN(Ncut)
            zlim = 0.d0
            do while (zlim.lt.Zp.and.xaux<300)
                xaux=xaux+1
                zaux= exp(alim*xaux+blim)
                zlim=zlim+zaux
            end do
            continue
        end do

        Cmax = xaux
        !print*, alim,blim,Cmax
    end subroutine LimitLine


    subroutine LineC60max (Ncut,Zp,rCN,aBE,bBE,half,a60,b60,C60max)
        !================================================================================================================
        !This subroutine obtains the Cmax60 line constants for an oil.
        !Ncut: integer input number of single cuts being considered in the oil.
        !rCN: input array which contains the set of carbon numbers.
        !Zp: input mole fractions of the residual fraction.
        !aBE & bBE: input constants from the Best Extrapolation
        !a60: output real variable. Slope of the limit line.
        !b60: output real variable. Intercept of the limit line.
        !C60max: output CN at which Zp is reached, as the summation of single z(i) from the Cmax60 distribution.
        !================================================================================================================
        implicit none
        integer, parameter :: imax=48
        integer::  Ncut,xaux,C60max
        integer, dimension(imax) :: rCN
        real(pr) :: Zp,aBE,bBE,half,a60,b60, zsum, crossCN, zcross, zaux, Ftol, atol, range
        real(pr) :: a_old, Zp60, F, dFdA

        crossCN = rCN(Ncut)+half   ! typically 19.51
        zcross= exp(aBE*crossCN+bBE)
        Ftol=1d0
        atol=1d0
        a60=aBE
        range = 60.5d0 - crossCN

        do while (atol.gt.1e-6.or.Ftol.gt.1e-6)
            a_old=a60
            Zp60 = (exp(a60*range+log(zcross))-zcross)/a60
            F = Zp - Zp60
            !dFdA = - (60.5*(a60*Zp60+zcross)-Zp60) / a60
            dFdA = - (range*(a60*Zp60+zcross)-Zp60) / a60
            a60=a_old-F/dFdA
            atol=abs(a_old-a60)
            Ftol=abs(F)
        end do

        b60=log(zcross)-a60*crossCN
        zsum=0d0
        xaux=rCN(Ncut)

        do while (zsum.lt.Zp)
            xaux=xaux+1
            zaux= exp(a60*xaux+b60)
            zsum=zsum+zaux
        end do

        C60max=xaux
        !print*, a60,b60,C60max

    end subroutine LineC60max


    subroutine GetNewMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,a,b,rMp,maxC,z,zpi,zMpi,Ndef,zp,w,rn,rMWplus,zcomp,zdef,rMW,rMdef,Mglobal)
      !! this subroutines for a C value returns the correponding  new M20+ value
      
      implicit none

      integer, parameter :: imax=48, maxD=15
      logical :: start
      real(pr) :: C, Z6p, zMp, a, b, rMp, rMWplus,Mglobal
      real(pr) :: aux, dold, dif, rMpold,zp
      real(pr), dimension(imax) :: zM  ! inputs coming from main program
      integer, dimension(imax) :: rCN ! inputs coming from main program
      real(pr), dimension(imax) :: z, rM,w,rn, zcomp,rMW
      integer :: Ncut, maxC, Ndef
      real(pr), dimension(300):: zpi, zMpi
      real(pr), dimension(maxD) :: zdef,rMdef

      
      if(nv1==1)then

          ! Find the new rMp for a fixed C constant 
          !C = 14  ! fixed value
          call difMpfromC(start,C,Ncut,rCN,zM,zMp,Z6p,dif,a,b,rMp,maxC,z,zpi,zMpi,zcomp,zdef,zp)

          !print*, C, maxC, rMp


      else

        ! Find the new rMp for a fixed C constant 
        !C = 14  ! fixed value
        Start = .true.
        rMp = rMWplus !inicial guess

          call difnewMp(Start,C,Ndef,Ncut,rCN,dold,a,b,rM,rMp,maxC,z,zp,zpi,zMpi,w,rn,zcomp,zdef,rMW,rMdef,Mglobal)
          rMpold = rMp
          rMp = 0.9*rMp
          Start = .false.

          call difnewMp(Start,C,Ndef,Ncut,rCN,dif,a,b,rM,rMp,maxC,z,zp,zpi,zMpi,w,rn,zcomp,zdef,rMW,rMdef, Mglobal)

          do while (abs(dif) > 0.00001)
            aux = rMp
            rMp = rMp - dif*(rMp-rMpold)/(dif-dold)
            rMpold = aux
            dold = dif
              call difnewMp(Start,C,Ndef,Ncut,rCN,dif,a,b,rM,rMp,maxC,z,zp,zpi,zMpi,w,rn,zcomp,zdef,rMW,rMdef, Mglobal)
          end do


          

          !print*, C, maxC, rMp
          !print*, zp , z6p
          !print*, z, zp
      end if


    end subroutine GetNewMpfromC


    subroutine difnewMp(Start,C,Ndef,Ncut,rCN,dif,a,b,rM,rMp,maxC,z,zp,zpi,zMpi,w,rn,zcomp,zdef,rMW,rMdef, Mglobal)
        !! this subroutines for a C value returns the correponding Mp value
        !use MassFracAndMoles, only: w,rn
        use data
        implicit none

        integer, parameter :: imax=48, maxD=15
        logical :: start
        real(pr) :: C, dif, a, b, rMp, zp ,aBE ,bBE ,r2 , alim, blim, half
        real(pr) :: a60,b60, sumz, rMpcalc, rntot, aold, Mglobal
        integer, dimension(imax) :: rCN,x   ! inputs coming from main program
        real(pr), dimension(imax) :: ylog, rM, z, w, rn, zcomp, ztotal,rMW, MWt
        integer :: i,Ndef, Ncut, maxC, i0, CmaxL, C60max
        real(pr), dimension(300):: zpi, zMpi,ZM_1
        real(pr), dimension(maxD) :: zdef,rMdef

        1       i0 = rCN(Ncut)

        if(nv == 1)then !introduce 'if' to choise where experimental data from
            x = rCN
            rn(Ndef+1:Ndef+Ncut) = w(Ndef+1:Ndef+Ncut)/rMW(1:Ncut)
            rn(Ndef+Ncut+1) = w(Ndef+Ncut+1)/rMp
            rntot = sum(rn(1:Ndef+Ncut+1))
            ztotal(1:Ndef+Ncut)  = rn(1:Ndef+Ncut) / rntot
            z(1:Ncut) = rn(Ndef+1:Ndef+Ncut) / rntot
            Zp = rn(Ndef+Ncut+1) / rntot
            !ylog = log(zcomp)
            ylog = log(z)
            
            !calculo m global
            MWt(1:Ndef) = rMdef(1:Ndef) 
            MWt(Ndef+1:Ndef+Ncut) = rMW(1:Ncut)
            ZM_1(1:Ndef+Ncut) =  ztotal(1:Ndef+Ncut)*MWt(1:Ndef+Ncut)
            
            !do i= 1, Ndef+Ncut
            !    print*, MWt(i)
            !end do

        end if

        if(nv == 2)then
            rM(1:Ncut) = 84 + C*(rCN(1:Ncut)-6)
            rn(Ndef+1:Ndef+Ncut) = w(Ndef+1:Ndef+Ncut)/rM(1:Ncut)
            rn(Ndef+Ncut+1) = w(Ndef+Ncut+1)/rMp
            rntot = sum(rn(1:Ndef+Ncut+1))
            ztotal(1:Ndef+Ncut)  = rn(1:Ndef+Ncut) / rntot
            z(1:Ncut) = rn(Ndef+1:Ndef+Ncut) / rntot
            Zp = rn(Ndef+Ncut+1) / rntot
            x = rCN
            ylog = log(z)
            !calculo m global
            MWt(1:Ndef) = rMdef(1:Ndef) 
            MWt(Ndef+1:Ndef+Ncut) = rM(1:Ncut)
            ZM_1(1:Ndef+Ncut) =  ztotal(1:Ndef+Ncut)*MWt(1:Ndef+Ncut)
            
            !do i= 1, Ndef+Ncut
            !    print*, MWt(i)
            !end do
        
        end if
        !print*, zp
        !print*, ztotal

        call BestLinearRegression(Ncut,x,ylog,Zp,aBE,bBE,r2)
        call LimitLine(Ncut,rCN,Zp,aBE,bBE,half,alim,blim,CmaxL)
        call LineC60max (Ncut,Zp,rCN,aBE,bBE,half,a60,b60,C60max)

        if(aBE < alim)then !, para C se excede el valor de 14 en la mayoria de los casos.
            a = alim
            b = blim
        else
            !call LineC60max (Ncut,Zp,rCN,aBE,bBE,half,a60,b60,C60max)
            if(aBE > a60)then
                a = a60
                b = b60
            else
                a = aBE
                b = bBE
            end if
        end if


        sumz = 0.0d0
        i = 0
        do while (sumz<Zp.and.i<300)
            i = i+1
            zpi(i) = exp(a*(i+i0)+b)
            sumz = sumz + zpi(i)
            zMpi(i) = zpi(i)*(84+C*(i+i0-6))
        end do
        !agregar condicion cuando C = 14 o C = 12 no olvidar consultar con martin
        if(i==300.and.sumz<Zp)then
            if(start) C = C - 0.07
            if(.not.start) C = C - 0.01
            go to 1
        end if

        do while (i<40.and.a>alim)
            aold = a
            a = max(1.0001*a, alim)
            b = b - (a-aold)*19.51
            sumz = 0.0d0
            i = 0
            do while (sumz<Zp.and.i<300)
                i = i+1
                zpi(i) = exp(a*(i+i0)+b)
                sumz = sumz + zpi(i)
                zMpi(i) = zpi(i)*(84+C*(i+i0-6))
            end do
        end do

        !   Adjustment to Zp (z20+)
        zpi(i) = zpi(i) - (sumz-Zp)
        zMpi(i) = zpi(i)*(84+C*(i+i0-6))
        ! ===============================
        rMpcalc = sum(zMpi(1:i))/Zp
        dif = rMpcalc-rMp
        maxC = i+i0

        ! calculo del peso molecular global
        Mglobal =  sum(ZM_1(1:Ndef+Ncut))+ sum(zMpi(1:i))
        !rint*, Mglobal
        !print*, maxC

        !do i = 1, maxC-i0
        !    print*, zpi(i)
        !end do
        !print*, '----------'
        !do i = 1, Ndef+Ncut
        !    print*, ztotal(i)
        !end do

    end subroutine difnewMp



    
      


end module routines
