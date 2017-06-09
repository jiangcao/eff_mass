program main

use types, only : dp 
use constants, only : hbar, m0
use eff_mass
use input

implicit none 

real(dp) :: x(100), y(100), me, k1(3), k2(3)
real(dp), allocatable :: dat(:,:)
integer  :: i 





call input_2c('energy.dat', dat)
! print *, dat(:,2)

k1 = (/0.0_dp,0.0_dp,0.0_dp  /)
k2 = (/3.98403e+6_dp, 7.96806e+6_dp, 3.98403e+6_dp  /)

me = calc_eff_mass(energies = dat(:,2), point1 = k1, point2 = k2, position = 0.0_dp)

print *, me







end program main
