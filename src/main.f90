program main

use types, only : dp 
use constants, only : hbar, m0
use eff_mass

implicit none 

real(dp) :: x(200), y(200)
integer  :: i 

do i = 1, 200
	x(i) = dble(i) / dble(100) * 1e6 - 1e6 
	y(i) = hbar**2 * x(i)**2 /2.0_dp/m0/0.2
enddo 



end program main