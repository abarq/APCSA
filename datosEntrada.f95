module datosEntrada
implicit none
integer, parameter :: nSamples = 71
character (len=8) :: expData = 'DATCRU'
contains

subroutine readData(energy, epsUnoE, epsDosE)
implicit none
real*8 :: energy(nSamples), epsUnoE(nSamples), epsDosE(nSamples)
real*8 :: coefN(nSamples), coefK(nSamples)
integer :: j


    open(3,file=expData)

    do j=1,nSamples
        read(3,*) energy(j), coefN(j), coefK(j)
    end do
    close(3)

    epsUnoE = coefN**2 - coefK**2
    epsDosE = 2*coefN*coefK

end subroutine
end module datosEntrada