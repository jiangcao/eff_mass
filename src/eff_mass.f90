module eff_mass

use types, only : dp
use constants, only : hbar, m0

private

public calc_eff_mass

contains

! energy in [eV]
! effective mass without unit
function calc_eff_mass(energies, point1, point2, position) result(me)
implicit none    
real(dp) , intent(in) :: energies(:)
real(dp) , intent(in) :: point1(3), point2(3), position
real(dp) :: x(size(energies)), vk(3), x0
integer  :: i, ind(6), indv(6)
	vk = point2 - point1
	x  = dble((/(i, i=1,size(x),1)/)) / dble(size(x)) * sqrt(dot_product(vk,vk))
	x0 = position * sqrt(dot_product(vk,vk))
	ind= find_neighbor(x0,3,x)   ! find the neighbor points around the point we want the derivative
	if (count(ind /= 0) < 3) then 
		print *, "problem: not enough points"
		call abort
	else
		indv(1:count(ind /= 0)) = pack(ind, ind /= 0)
		! then use x(indv(1:3)) and energy(indv(1:3)) to calculate the 2nd derivative 
        me = 1.0_dp / second_derivative(x(indv(1:3)), energies(indv(1:3))) * hbar**2 / 2.0_dp / m0
    endif
end function calc_eff_mass





function second_derivative(x,y) result(d2)
implicit none    
real(dp) , intent(in) :: x(3), y(3)
real(dp) :: d2 , h 
    h = (x(3)-x(1)) / 2.0_dp
    d2 = (y(3) - 2.0_dp*y(2) + y(1)) / h**2
end function second_derivative





        ! function returns the indices of the neighbor points around `x` in the data list `lst`
        function find_neighbor(x,n,lst) result(ind)
            implicit none
            real(dp) , intent(in) :: x, lst(:) 
            integer  , intent(in) :: n
            integer :: ind(2*n) , nx 
            real(dp) :: xnp(2*n) ! neighbor points value
            integer  :: npi(2*n) ! neighbor points indices
            integer  :: i, m
            nx = size(lst)
            npi = 0
            xnp(1:n) = -HUGE(1.0_dp)
            xnp(n+1:2*n) =  HUGE(1.0_dp)
            do i = 1, nx
                if ( (lst(i) <= x ) .and. (lst(i) > xnp(1)) ) then 
                    m = find_insert(lst(i), lst=xnp(1:n)) 
                    xnp(1:m-1) = xnp(2:m)
                    xnp(m) = lst(i)
                    npi(1:m-1) = npi(2:m)
                    npi(m) = i
                endif
                if ( (lst(i) > x) .and. (lst(i) < xnp(2*n))) then 
                    m = find_insert(lst(i), lst=xnp(n+1:2*n)) 
                    m = m+n
                    xnp(m+2:2*n) = xnp(m+1:n*2-1)
                    xnp(m+1) = lst(i)
                    npi(m+2:2*n) = npi(m+1:n*2-1)
                    npi(m+1) = i
                endif
            enddo
            ind = npi
        end function find_neighbor
 

        function find_insert(x, lst) result(ind)
            implicit none
            real(dp), intent(in) :: x , lst(:)
            integer :: ind
            integer :: i, nx
            nx = size(lst)
            i = 1
            do while ( ( lst(i) < x ) .and. (i<= nx) )
                i = i+1
            enddo
            ind = i-1
        end function find_insert



end module eff_mass