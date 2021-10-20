module funcionCosto
use randomNum
use datosEntrada
implicit none
real*8 :: fnCosto
integer :: nOs = 5
integer, parameter :: nParam = 23

contains

subroutine fn(epsilonSumaTotal, x, epsUnoE, epsDosE, energy)
use datosEntrada
implicit none
real*8 :: omegaPe, f0, gammaOE, beta, eta,epshf, omegaPe2,OPh2, gammaE, gammaH 
real*8 :: x(nParam)
real*8 :: omega
real*8 :: dummy1, dummy2, dummy3, dummy4,suma1, suma2, epsilonSumaTotal 
real*8 :: suma1a, suma1b, suma1c, suma1d, suma1e
real*8 :: suma2a, suma2b, suma2c, suma2d, suma2e
real*8 :: energy(nSamples), epsUnoE(nSamples), epsDosE(nSamples)
real*8 :: drude1(nSamples), drude2(nSamples), lorentz1(nSamples), lorentz2(nSamples), epsilon1(nSamples)
real*8 :: epsilonTotal(nSamples), epsilon2(nSamples),difEpsUno(nSamples), difEpsDos(nSamples)
integer :: i, k, kk, nDrudeParameters 
complex*16 :: II

nDrudeParameters = 7
omegaPe = x(1)
f0 = x(2)
gammaOE = x(3)
beta = x(4)
eta = x(5)
epshf = 1/f0
omegaPe2 = omegaPe**2
OPh2 = omegaPe2*beta
II = dcmplx(0.0D0,1.D0)


dummy1 = 1/x(6) 
dummy2 = 1/x(7)
dummy3 = 1/x(nParam) 

do i=1,nSamples
    omega = energy(i)
    gammaE = gammaOE*(1+(omega*dummy1)**2)
    gammaH = eta*gammaOE*(1+(omega*dummy2)**2)
    drude1(i) = epshf - omegaPe2/(omega**2 + gammaE**2) - OPh2/(omega**2 + gammaH**2)
    drude2(i) = (gammaE*omegaPe2)/(omega*(omega**2 + gammaE**2)) + (OPh2*gammaH)/(omega*(omega**2 + gammaH**2))
	
    suma1 = 0.0D0
    suma2 = 0.0D0
	
    ! do j=1,nOs
        ! k = j*3
        ! kk = nDrudeParameters + k
        ! suma1 = suma1 +  ((x(kk-2)*omegaPe2)*(x(kk)**2-omega**2))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2); 
        ! suma2 = suma2 +  ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2); 		
    ! end do	

        k = 1*3
        kk = nDrudeParameters + k
        suma1a =  ((x(kk-2)*omegaPe2)*(x(kk)**2-omega**2))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2) 
        suma2a =  ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2)

        k = 2*3
        kk = nDrudeParameters + k
        suma1b = ((x(kk-2)*omegaPe2)*(x(kk)**2-omega**2))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2) 
        suma2b = ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2)

        k = 3*3
        kk = nDrudeParameters + k
        suma1c = ((x(kk-2)*omegaPe2)*(x(kk)**2-omega**2))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2) 
        suma2c = ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2)


        k = 4*3
        kk = nDrudeParameters + k
        suma1d = ((x(kk-2)*omegaPe2)*(x(kk)**2-omega**2))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2) 
        suma2d = ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2)
		
        k = 5*3
        kk = nDrudeParameters + k
        suma1e =  ((x(kk-2)*omegaPe2)*(x(kk)**2-omega**2))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2) 
        suma2e =  ((x(kk-2)*omegaPe2)*(omega*x(kk-1)))/((x(kk)**2 - omega**2)**2 + (omega*x(kk-1))**2)

    suma1 = suma1a+suma1b+suma1c+suma1d+suma1e
    suma2 = suma2a+suma2b+suma2c+suma2d+suma2e
	
    lorentz1(i) = dummy3*suma1
    lorentz2(i) = dummy3*suma2
end do


epsilon1 = drude1 + lorentz1
epsilon2 = drude2 + lorentz2


! do i=1,nSamples
	! omega = energy(i)
	! gammaE = gammaOE*(1+(omega*dummy1)**2)
	! gammaH = eta*gammaOE*(1+(omega*dummy2)**2)
	! eps = epshf - 1/(omega)*(omegaPe2/(omega + II*gammaE) + OPh2/(omega + II*gammaH))
	
	! do j=1,nOs
		! k = j*3
		! kk = nDrudeParameters + k
		! eps = eps +  dummy3*(((x(kk-2)*omegaPe2))/(x(kk)**2 - omega**2 - II*(omega*x(kk-1)))); 
		
	! end do
	
	! epsilon1(i) = dreal(eps)
	! epsilon2(i) = dimag(eps)
! end do





dummy4 = 1/dble(2*nSamples-nParam-1);

difEpsUno = ((epsilon1 - epsUnoE)/epsUnoE)**2
difEpsDos = ((epsilon2 - epsDosE)/epsDosE)**2

epsilonTotal = difEpsUno + difEpsDos

epsilonSumaTotal = sum(epsilonTotal)*dummy4

end subroutine





end module funcionCosto