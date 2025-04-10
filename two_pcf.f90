module two_pcf
    implicit none
    include 'mpif.h'
    contains
    

    subroutine cal_corr(my_data1,my_data2,bins,range_min,range_max,hist_all)

        real, allocatable :: distance(:) 
        real, dimension(:, :) :: my_data1(:,:),my_data2(:,:)
        real :: range_min, range_max
        integer :: ierr, ThisTask, NTaskWithN
        integer :: data_Number, each_Number, extra
        integer :: bins
        integer :: hist(bins), hist_sum(bins)
        integer :: start_idx, end_idx, ii
        real :: distance_x2(size(my_data2,1)), distance_y2(size(my_data2,1)), distance_z2(size(my_data2,1))
        integer ,intent(out) :: hist_all(bins)


        call MPI_Comm_rank(MPI_COMM_WORLD, ThisTask, ierr)

        call MPI_Comm_size(MPI_COMM_WORLD, NTaskWithN, ierr)

        data_Number = size(my_data1, dim=1)

        each_Number = data_Number / NTaskWithN
        extra       = mod(data_Number, NTaskWithN)

        if (ThisTask < extra) then
            each_Number = each_Number + 1
            start_idx = ThisTask * each_Number + 1
        else
            start_idx = ThisTask * each_Number + extra + 1
        end if

        end_idx = start_idx + each_Number - 1

        hist_sum=0

        do ii=start_idx, end_idx
            distance_x2= (my_data2(:,1)-my_data1(ii,1))**2
            distance_y2= (my_data2(:,2)-my_data1(ii,2))**2
            distance_z2= (my_data2(:,3)-my_data1(ii,3))**2

            distance=sqrt(distance_x2+distance_y2+distance_z2)

            call compute_histogram(distance,bins,range_min,range_max,hist)
            hist_sum = hist_sum+hist

        enddo

        call MPI_Allreduce(hist_sum, hist_all, bins, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)


    end subroutine cal_corr

    subroutine compute_histogram(data, nbins, range_min, range_max,hist)
        implicit none
        real, intent(in) :: data(:)
        integer, intent(in) :: nbins
        real, intent(in) :: range_min, range_max
        integer, intent(out) :: hist(nbins)
        
        integer :: i, bin
        real :: bin_width, value
        
        hist = 0
        bin_width = (range_max - range_min) / real(nbins)
        
        do i = 1, size(data)
            value = data(i)
        
            ! 범위 내에 있을 때만 처리
            if (value >= range_min .and. value < range_max) then
            bin = int((value - range_min) / bin_width) + 1
            hist(bin) = hist(bin) + 1
            else if (value == range_max) then
            ! 가장 마지막 bin 포함 (edge case)
            hist(nbins) = hist(nbins) + 1
            end if
        end do
    end subroutine compute_histogram
          
    
end module two_pcf