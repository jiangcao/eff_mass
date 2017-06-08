program main

use types, only : dp 
use constants, only : hbar, m0
use eff_mass

implicit none 

real(dp) :: x(100), y(100), me
integer  :: i 

do i = 1, 100
	x(i) = dble(i) / dble(100) * 1e6
	y(i) = hbar**2 * x(i)**2 /2.0_dp/m0/0.5
enddo 
write(10,*) x 
write(20,*) y
me = calc_eff_mass(energies=y, point1=(/0.0_dp,0.0_dp,0.0_dp/), &
						point2=(/0.0_dp,0.0_dp,1.0_dp/)*1e6, position=0.5_dp)
print *, me

end program main