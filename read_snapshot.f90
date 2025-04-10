module read_snapshot
    use hdf5
    implicit none
    contains
    

    subroutine read_CDM_pos(filename,coords)
        character(len=*), intent(in) :: filename
        integer(HID_T) :: file_id, group_id, dset_id, space_id
        integer :: hdferr  
        integer(HSIZE_T), dimension(2) :: dims2d

        integer, parameter :: sp = kind(1.0)

        real(sp), allocatable, intent(out) :: coords(:,:)

        call h5open_f(hdferr)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferr)

        call h5gopen_f(file_id, "PartType1", group_id, hdferr)
        call h5dopen_f(group_id, "Coordinates", dset_id, hdferr)
        call h5dget_space_f(dset_id, space_id, hdferr)
        call h5sget_simple_extent_dims_f(space_id, dims2d, rank=2, hdferr)
      
        allocate(coords(dims2d(1), dims2d(2)))

        call h5dread_f(dset_id, H5T_NATIVE_REAL, coords, dims2d, hdferr)
      
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
        call h5gclose_f(group_id, hdferr)
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)

    end subroutine read_CDM_pos
          
    
end module two_pcf