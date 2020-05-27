program expect_max
  !
  implicit none
  !
  integer :: n, m, i
  real(8),allocatable :: a(:), x(:)
  !
  read(*,*) n
  allocate(a(n), x(n))
  read(*,*) a(1:n)
  !
  do m = 1, n
     x(m) = a(m)
     do i = 1, n - m
        x(m) = a(m + i) + dble(i) / dble(m + i - 1) * x(m)
     end do
     x(m) = dble(m) / dble(n) * x(m)
  end do
  !
  do m = 1, n
     write(*,*) x(m)
  end do
  !
end program expect_max
