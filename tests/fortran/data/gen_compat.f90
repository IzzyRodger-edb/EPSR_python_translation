! Generates compat_expected.txt: reference strings and bit patterns from gfortran for the
! Python compatibility layer tests (tests/fortran/test_compat.py).
! Build: gfortran -O2 -fno-tree-vectorize gen_compat.f90 -o gen_compat && ./gen_compat > compat_expected.txt
program gen_compat
  implicit none
  real :: v(24), x, y
  real*8 :: d(8)
  integer :: i, j, k, iv(6)
  character*11 :: c11
  v = (/ 1.0, -0.123456789, 0.0, 1.82, -2.44, 0.1, 0.0999999, 1.0e8, 1.0e9, -0.000123, &
         -40.6027222, 99999999.0, -0.0005, 12.3456, 3.3333333e-5, 1.0e-30, 7.5e35, -1.5, &
         0.5, 2.5, -2.5, 0.49999997, 123456.789, -9.87654e-12 /)
  d = (/ 1.0d0, 0.1d0, -3.3333333333333d-5, 1.0d20, -2.5d0, 123456789.123d0, 1.0d-300, 0.0d0 /)
  iv = (/ 0, 7, -7, 12345, -2147483647, 42 /)
  do i = 1, 24
    write(*,'(a,i0)') '#v ', i
    write(*,'(e14.7,1x,e13.6,1x,e12.5,1x,e9.3,1x,e10.3e3)') v(i), v(i), v(i), v(i), v(i)
    write(*,'(f10.5,1x,f8.3,1x,f6.3,1x,f14.8,1x,f10.2,1x,f7.0)') v(i), v(i), v(i), v(i), v(i), v(i)
    write(*,*) v(i)
    write(*,*) i, v(i), 'ab', 'cd', v(i)
  end do
  do i = 1, 8
    write(*,'(a,i0)') '#d ', i
    write(*,*) d(i)
    write(*,'(e25.16,1x,f30.10)') d(i), d(i)
  end do
  write(*,'(a)') '#i'
  write(*,'(i0,1x,i5,1x,i2.2,1x,i3.3,1x,i4,1x,i11)') iv(1), iv(2), iv(3), iv(4), iv(4), iv(5)
  write(*,*) iv
  write(*,'(1x,i4,(5(1x,i4,1x,e9.3)))') 7, (i, real(i)+0.5, i=1,7)
  c11 = 'abc'
  write(*,'(a11,1x,a,a5,a2)') c11, 'x', 'yz', 'pqr'
  write(*,'(3x,a,1h-,a,09x,a,1h-,a,05x)') 'S', 'C', 'C', 'S'
  write(*,'(l1,1x,l3)') .true., .false.
  write(*,*) .true., 3, 'text'
  ! intrinsics: print bit patterns of results
  write(*,'(a)') '#math'
  do i = 1, 24
    x = v(i)
    y = abs(x) + 0.25
    write(*,'(10(1x,i0))') transfer(exp(x/100.0),k), transfer(log(y),k), transfer(sin(x),k), &
      transfer(cos(x),k), transfer(atan2(x,y),k), transfer(sqrt(y),k), transfer(y**1.37,k), &
      transfer(y**(i-12),k), nint(x), int(x)
  end do
  write(*,'(a)') '#imod'
  write(*,'(8(1x,i0))') mod(7,3), mod(-7,3), mod(7,-3), mod(-7,-3), 7/3, -7/3, 7/(-3), (-7)/(-3)
end program gen_compat
