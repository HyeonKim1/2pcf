program main
  use two_pcf
  use read_snapshot
  use hdf5
  implicit none
  include 'mpif.h'

  character(len=100) :: filename
  integer :: bins, ierr
  real, allocatable :: my_data1(:,:), my_data2(:,:)
  real :: range_min, range_max
  integer :: N
  integer, parameter :: sp = kind(1.0)
  real(sp), allocatable, intent(out) :: coords(:,:)

  call MPI_Init(ierr)

  if (ThisTask == 0) then
    print *, "Enter the input file name:"
  endif

  read(*,*) filename

  call read_CDM_pos(filename,coords)
  call cal_corr(coords,coords)

  

  call MPI_Finalize(ierr)
  
end program main