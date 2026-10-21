! Driver used to generate ran1_reference.txt (tests/library/test_ran1.py).
!
! Build and run (from the repository root), with the flags of
! tools/build_fortran_reference.sh:
!   gfortran -O2 -fno-tree-vectorize -fno-backslash -ffixed-line-length-none -fopenmp -w \
!     -Ifortran_src fortran_src/lngamma_mod.f90 fortran_src/parse.f \
!     fortran_src/auxiliary_routines.f tests/library/data/ran1_driver.f90 -o ran1_driver
!   ./ran1_driver > tests/library/data/ran1_reference.txt
!
! For each starting seed (idum, with iy = 0 so that ran1 initialises), 1000 deviates
! are printed as the integer bit pattern of the real*4 value, followed by the final
! state written by saveran1 (list-directed).  The last block continues the generator
! from a state restored with initran1.
program ran1_driver
   implicit none
   integer, parameter :: ntab = 32
   integer :: mcall, idum, iy, iv(ntab)
   common /ran1com/ mcall, idum, iy, iv
   real, external :: ran1
   integer, parameter :: nseed = 6
   integer :: seeds(nseed) = (/ -1, -12345, 0, -2147483647, 123, -987654321 /)
   integer :: is, i
   real :: r

   do is = 1, nseed
      mcall = 0
      idum = seeds(is)
      iy = 0
      iv = 0
      write(6, '(a,i12)') 'seed ', seeds(is)
      do i = 1, 1000
         r = ran1()
         write(6, '(i12)') transfer(r, i)
      end do
      call saveran1(6)
   end do
   ! save the state to a file, read it back with initran1 and continue
   open(10, file='ran1_state.tmp', status='unknown')
   call saveran1(10)
   close(10)
   mcall = 0
   idum = 0
   iy = 0
   iv = 0
   open(10, file='ran1_state.tmp', status='old')
   call initran1(10)
   close(10, status='delete')
   write(6, '(a)') 'continued'
   do i = 1, 1000
      r = ran1()
      write(6, '(i12)') transfer(r, i)
   end do
   call saveran1(6)
end program ran1_driver
