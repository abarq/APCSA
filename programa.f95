program programa
use randomNum
use datosEntrada
use funcionCosto
implicit none

integer :: m, i, nMovAcep, nMovRechz, ii
integer, parameter :: N = 7500,  ciclosTempMax = 250, iteracionesMax = ciclosTempMax*N
real*8 :: alfa, Tfin, Tini
real*8 :: T, deltaE, E
real*8 :: L(nParam), U(nParam), deltaP(nParam), stepP(nParam), x(nParam), xOld(nParam)
real*8 :: xn, EOld, probAcep, bestCostFnTemp(ciclosTempMax), desEstan, desviacion
real*8 :: start, finish, EMem, xMem(nParam), TAPCSA, energiaMovAcep(N), aceptanciaHistorica(N)
real*8 :: acceptanceRateTemp(ciclosTempMax)
real*8 :: energy(nSamples), epsUnoE(nSamples), epsDosE(nSamples)
real*8 :: dif, numAlea



semilla = 312312166
alfa = 0.9D0
Tini = 50
Tfin = 1.0D-9


call cpu_time(start)
call readData(energy, epsUnoE, epsDosE)
call getBounds(U,L)

deltaP = U-L
stepP = deltaP/dble(1000)
x = 0.0D0

do i=1,nParam
    !call rng%generate_single(xn)
    call randNum(xn)
    xOld(i) = L(i) + xn*deltaP(i)
end do
call fn(EOld, xold, epsUnoE, epsDosE, energy)

EMem = EOld
xMem = xOld
T = Tini
m = 1
bestCostFnTemp(m) = EOld
TAPCSA = Tini


do while (m.lt.ciclosTempMax)

    nMovAcep = 0
    nMovRechz = 0
    do i=1,N
        
        !call crearMovimiento(xOld, x, L, U, stepP)
        
        do ii=1, nParam
            call randNum(xn)
            !call rng%generate_single(xn)
            numAlea = -1.0D0 + 2.0D0*xn
            x(ii) = xOld(ii) +  numAlea*stepP(ii)
            
            if (x(ii).gt.U(ii)) then
                dif = x(ii) - U(ii)
                x(ii) = U(ii)-dif
            end if
            
            if (x(ii).lt.L(ii)) then
                dif = abs(x(ii)-L(ii))
                x(ii) = L(ii)+dif
            end if
        
        end do
        
        
        call fn(E, x, epsUnoE, epsDosE, energy)
        
        deltaE = E-EOld
        
        if (deltaE.lt.0.0D0) then
        
            EOld = E
            xOld = x
            call actualizarRegistro(xMem, EMem, xOld, EOld)
            nMovAcep = nMovAcep + 1
            energiaMovAcep(nMovAcep) = EOld
            aceptanciaHistorica(nMovAcep) = 1.0D0
        else 
            probAcep = exp(-deltaE/T)
            call randNum(xn)
            
            if (xn.lt.probAcep) then
                EOld = E
                xOld = x
                call actualizarRegistro(xMem, EMem, xOld, EOld)
                nMovAcep = nMovAcep + 1
                energiaMovAcep(nMovAcep) = EOld
                aceptanciaHistorica(nMovAcep) = probAcep
            else
                nMovRechz = nMovRechz + 1 
            end if
        end if
    
    
    
    end do
    
    acceptanceRateTemp(m) = float(nMovAcep)/(float(nMovAcep) + float(nMovRechz))
    bestCostFnTemp(m) = EOld
    call calcularDesviacion(energiaMovAcep,nMovAcep, desviacion)
    call modificarPaso(x,stepP)
    call calcularDesviacionEstandar(aceptanciaHistorica, nMovAcep, desEstan)
    TAPCSA = -desviacion/(log(0.9D0) - ((dble(m)**2)/(2*desEstan**2)))
    T = alfa*T
    m = m+1
end do



call cpu_time(finish)
write(*,*) EMem
write(*,*) finish-start

end program


subroutine getBounds(U,L)
use funcionCosto
implicit none

! integer, parameter :: nParam = 23
real*8 :: L(nParam), U(nParam)


L(1) = 1.0D0
U(1) = 10D0

L(2) = 0.0D0
U(2) = 1.0D0

L(3) = 0.0D0
U(3) = 0.1D0

L(4) = 0.0D0
U(4) = 1.0D0

L(5) = 0.5D0
U(5) = 1.5D0

L(6) = 0.5D0
U(6) = 2.5D0

L(7) = 2.0D0
U(7) = 3.0D0

L(8) = 0.0D0
U(8) = 5.0D0

L(9) = 0.0D0
U(9) = 5.0D0

L(10) = 0.40D0
U(10) = 0.42D0

L(11) = 0.0D0
U(11) = 5.0D0

L(12) = 0.0D0
U(12) = 5.0D0

L(13) = 1.214D0
U(13) = 1.226D0

L(14) = 0.0D0
U(14) = 5.0D0

L(15) = 0.0D0
U(15) = 5.0D0

L(16) = 2.7D0
U(16) = 2.9D0

L(17) = 0.0D0
U(17) = 5.0D0

L(18) = 0.0D0
U(18) = 5.0D0

L(19) = 4.9D0
U(19) = 5.0D0

L(20) = 0.0D0
U(20) = 5.0D0

L(21) = 0.0D0
U(21) = 5.0D0

L(22) = 11.0D0
U(22) = 11.1D0

L(23) = 0.1D0
U(23) = 0.4D0



end

subroutine solWilliam(xold)
implicit none
integer, parameter :: nParam = 23
real*8 :: xold(nParam)

xold(1) = 5.145D0
xold(2) = 0.552D0
xold(3) = 0.009D0
xold(4) = 0.244D0
xold(5) = 0.965D0
xold(6) = 1.287D0
xold(7) = 2.265D0
xold(8) = 0.344D0
xold(9) = 0.666D0
xold(10) = 0.411D0
xold(11) = 0.710D0
xold(12) = 1.969D0
xold(13) = 1.216D0
xold(14) = 0.372D0
xold(15) = 3.423D0
xold(16) = 2.777D0
xold(17) = 0.358D0
xold(18) = 3.093D0
xold(19) = 4.974D0
xold(20) = 0.999D0
xold(21) = 4.352D0
xold(22) = 11.06D0
xold(23) = 0.37D0

end subroutine


subroutine crearMovimiento(xOld, x, L, U, stepP)
use funcionCosto
use randomNum
implicit none
real*8 :: numAlea, xn, dif
real*8 :: L(nParam), U(nParam), stepP(nParam), x(nParam), xOld(nParam)
integer i 

do i=1, nParam
    call randNum(xn)
    numAlea = -1.0D0 + 2.0D0*xn
    x(i) = xOld(i) +  numAlea*stepP(i)
    
    if (x(i).gt.U(i)) then
        dif = x(i) - U(i)
        x(i) = U(i)-dif
    end if
    
    if (x(i).lt.L(i)) then
        dif = abs(x(i)-L(i))
        x(i) = L(i)+dif
    end if
    
end do




end subroutine

subroutine modificarPaso(x,stepP)
use funcionCosto
implicit none
integer i
real*8 :: dummy1,dummy2
real*8 :: stepP(nParam), x(nParam)

do i=1,nParam
    
    dummy2 = 0.005D0
    dummy1 = stepP(i)/x(i)
    
    if (dummy1.gt.dummy2) then
        stepP(i) = stepP(i)/1.0D0
    end if
    
end do

end subroutine

subroutine actualizarRegistro(xMem, EMem, xOld, EOld)
use funcionCosto
implicit none
real*8 :: EOld, EMem
real*8 :: xMem(nParam), xOld(nParam)

    if (EOld.lt.EMem) then
        EMem = EOld
        xMem = xold
    end if

end subroutine


subroutine calcularDesviacion(energiaMovAcep,nMovAcep, desviacion)
implicit none

integer, parameter :: N = 7500
real*8 :: energiaMovAcep(N), suma, prom, desviacion
integer :: nMovAcep, i

suma = 0.0D0
do i=1,nMovAcep
    suma = suma + energiaMovAcep(i)
end do

prom = suma/float(nMovAcep)

suma = 0.0D0

do i=1,nMovAcep
    suma = suma + abs(energiaMovAcep(i)-prom)
end do

desviacion = suma/float(nMovAcep)



end subroutine


subroutine calcularDesviacionEstandar(aceptanciaHistorica, nMovAcep, desEstan)
implicit none
integer, parameter :: N = 7500
real*8 :: aceptanciaHistorica(N), suma, prom, desEstan
integer :: nMovAcep, i

suma = 0.0D0
do i=1,nMovAcep
    suma = suma + aceptanciaHistorica(i)
end do

prom = suma/float(nMovAcep)

suma = 0.0D0

do i=1,nMovAcep
    suma = suma + (aceptanciaHistorica(i)-prom)**2
end do

desEstan = sqrt(suma/float(nMovAcep))



end subroutine




