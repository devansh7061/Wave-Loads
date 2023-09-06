! Code to generate wavelength iteratively using Dispersion relation

module wavelength
    implicit none
    real, parameter, public:: pi = 4*atan(1.0), g = 9.81
    contains
    function findWavelength(T,d) result(L)
        real, intent(in) :: T, d
        real :: L, k, w, sum1, sum2, sum3, temp
        k = 1/d
        w = 2*pi/T
        temp = 1
        do while (temp .ne. 0)
            sum1 = tanh(k*d)
            sum2 = w**2/(g*k)
            sum3 = sum1 - sum2
            if (abs(sum3) .le. 0.000001) then
                temp = 0
            else if (sum3 .ge. 0) then
                k = k - abs(sum3)*0.001
            else
                k = k + abs(sum3)*0.001
            end if
        end do
        L = 2*pi/k
    end function findWavelength
end module wavelength

! Code to generate velocity using Airy Wave theory and wavelength which we got from wavelength module
module velocity
    use wavelength
    implicit none
    real:: H, d, x, T, k, w
    contains
    function findVelocity (time, z) result(u)
        real, intent(in):: time, z
        real :: u
        ! Values given in question
        H = 18
        d = 85
        x = 5
        T = 14
        w = 2*pi/T
        k = 2*pi / findWavelength(T, d)
        u = (H/2)*g*k*cosh(k*(d+z))*cos(k*x - w*time)/(w*cosh(k*d))
    end function findVelocity
end module velocity

! Code to generate acceleration using Airy Wave theory
module acceleration
    use velocity
    implicit none
contains
    function findAcceleration(time, z) result(acc)
        real, intent(in) :: time, z
        real acc
        acc = (H/2)*g*k* cosh(k*(d+z)) * sin(k*x-w*time) / cosh(k*d)
    end function findAcceleration
end module acceleration

! Code to generate force using Wheeler Stretching Method
program forces
    use acceleration
    implicit none
    real:: div, ux, ax, zbar, force, cd, cm, rho, dia, area, eta, ti, div_height,height, sum, sumex, extra, maxtime
    div = 100
    rho = 1025
    cd = 0.7
    cm = 2
    rho = 1025
    dia = 0.5
    maxtime = 42
    ti = 0
    open(3, file='forces_wheeler.dat')

    do while (ti <= maxtime)
        sum = 0
        extra = 0
        div_height = 0.1
        area = (pi*dia**2)/4
        eta = (H/2)*cos(k*x - w*ti)
        height = eta
        if (eta >= 0) then
            height = 0
            extra = 1
        end if
        do while (height > -(d))
            ax = findAcceleration(ti, height)
            ux = findVelocity(ti, height)
            force = 0.5*rho*cd*dia*ux*abs(ux) + rho*cm*area*ax
            sum = sum + force
            height = height - div_height
        end do

        if (extra == 1) then
            sumex = sum + (0.5*rho*cd*dia*findVelocity(ti,0.0)*abs(findVelocity(ti, 0.0)) &
            & + rho*cm*findAcceleration(ti, 0.0)*area)*eta/0.1
        else
            sumex = sum
        end if
        sum = 0
        height = eta
        do while (height > -(d))
            ax = findAcceleration(ti, height)
            ux = findVelocity(ti, height)
            zbar = (height - eta) * (d/(d+eta))
            ax = findAcceleration(ti, zbar)
            ux = findVelocity(ti, zbar)
            force = 0.5*rho*cd*dia*ux*abs(ux) + rho*cm*area*ax
            sum = sum + force
            height = height - div_height
        end do
        write(3, *) ti, sum
        ti = ti + 0.1
    end do
    call execute_command_line('gnuplot -p plot_forces_wheeler.plt')
end program forces