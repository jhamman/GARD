module output_mod
    use io_routines, only : io_write
    use data_structures
    use string, only : str

    implicit none

    interface shift_z_dim
        module procedure shift_z_dim_3d, shift_z_dim_4d
    end interface shift_z_dim
contains
    subroutine write_output(output, options)
        implicit none
        type(config),  intent(in)   :: options
        type(results), intent(in)   :: output

        character(len=MAXFILELENGTH)            :: filename
        character(len=MAXFILELENGTH)            :: dimnames(3)
        integer                                 :: dimids(2)
        real, dimension(:,:,:),     allocatable :: output_data
        real, dimension(:,:,:,:),   allocatable :: output_data_4d
        integer :: nvars, i, nx, ny, nt, nv
        integer :: Mem_Error
        character(len=19) :: todays_date_time
        CHARACTER(32) :: username, hostname
        character(len=MAXFILELENGTH) :: err

        character (len = *), parameter :: units = "units"
        character (len = *), parameter :: long_name = "long_name"
        character (len = *), parameter :: standard_name = "standard_name"
        character (len = *), parameter :: comment = "comment"
        character (len = *), parameter :: note = "note"
        character (len = *), parameter :: coordinates = "coordinates"

        character (len = *), parameter :: lat_longname = "latitude"
        character (len = *), parameter :: lat_units = "degrees_north"
        character (len = *), parameter :: lat_standard_name = "latitude"
        character (len = *), parameter :: lon_longname = "longitude"
        character (len = *), parameter :: lon_units = "degrees_east"
        character (len = *), parameter :: lon_standard_name = "longitude"
        character (len = *), parameter :: mask_longname = "domain mask"

        nvars = size(output%variables)

        nt = size(output%variables(1)%data, 1)
        nx = size(output%variables(1)%data, 2)
        ny = size(output%variables(1)%data, 3)
        nv = size(output%variables(1)%predictors, 4)

        allocate(output_data(nx,ny,nt), stat=Mem_Error )
        if (Mem_Error /= 0) call memory_error(Mem_Error, "output_data", [nx,ny,nt])

        ! Metadata for file
        ! now
        call date_and_time(values=date_time,zone=UTCoffset)
        date_format='(I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)'
        write(todays_date_time,date_format) date_time(1:3),date_time(5:7)
        ! username
        call getlog(username)
        ! hostname
        status = hostnm(hostname)

        ! Create the file.
        call check(nf90_create(filename, ior(NF90_HDF5, NF90_CLOBBER), ncid), trim(err))

        ! Define the dimensions. The time dimension is defined to have
        ! unlimited length - it can grow as needed.
        call check(nf90_def_dim(ncid, y_dim_name, ny, y_dimid), trim(err))
        call check(nf90_def_dim(ncid, x_dim_name, nx, x_dimid), trim(err))
        call check(nf90_def_dim(ncid, "time", nf90_unlimited, time_dimid), trim(err))
        ! call check(nf90_def_dim(ncid, "vars", nv, vars_dimid), trim(err))

        ! Define the coordinate variables.
        if (options%coords_2d) then
            dimids = (/ y_dimid, x_dimid /)
            call check(nf90_def_var(ncid, lat_name, nf90_real, dimids, lat_varid), trim(err))
            call check(nf90_def_var(ncid, lon_nate, nf90_real, dimids, lon_varid), trim(err))
            coords = (/ lat_name, lon_name /)
        else
            call check(nf90_def_var(ncid, lat_name, nf90_real, lat_dimid, lat_varid), trim(err))
            call check(nf90_def_var(ncid, lon_name, nf90_real, lon_dimid, lon_varid), trim(err))
        endif

        call check(nf90_def_var(ncid, "time", NF90_REAL, time_dimid, time_varid), trim(err))
        ! TODO: add vars coordinate
        ! call check(nf90_def_var(ncid, "vars", type?, vars_dimid, vars_varid), trim(err))

        ! Assign units attributes to coordinate variables.
        call check(nf90_put_att(ncid, lat_varid, longname, lat_longname), trim(err))
        call check(nf90_put_att(ncid, lat_varid, units, lat_units), trim(err))
        call check(nf90_put_att(ncid, lat_varid, standard_name, lat_standard_name), trim(err))

        call check(nf90_put_att(ncid, lon_varid, longname, lon_longname), trim(err))
        call check(nf90_put_att(ncid, lon_varid, units, lon_units), trim(err))
        call check(nf90_put_att(ncid, lon_varid, standard_name, lon_standard_name), trim(err))

        call check(nf90_put_att(ncid, time_varid, longname, time_longname), trim(err))
        call check(nf90_put_att(ncid, time_varid, units, time_units), trim(err))
        call check(nf90_put_att(ncid, time_varid, calendar, time_calendar), trim(err))

        ! mask variable
        dimids = (/ y_dimid, x_dimid /)
        call check(nf90_def_var(ncid, "mask", nf90_int, dimids, mask_varid), trim(err))
        call check(nf90_put_att(ncid, mask_varid, longname, mask_longname), trim(err))
        call check(nf90_put_att(ncid, mask_varid, note, "unitless"), trim(err))
        call check(nf90_put_att(ncid, mask_varid, comment, "0 value indicates cell is not active"), trim(err))

        ! Add coordinates attribute
        if (options%coords_2d) then
            call check(nf90_put_att(ncid, mask_varid, coordinates, coords), trim(err))
        endif

        ! loop over each prediction variable and
        do i=1,nvars
            dimids = (/ time_dimid, y_dimid, x_dimid /)

            call check(nf90_def_var(ncid, trim(output%variables(i)%name),
                                    NF90_REAL, dimids, varid), trim(err))
            call check(nf90_def_var_deflate(ncid, varid,
                       shuffle=options%comp_level, deflate=options%comp_level,
                       deflate_level=options%comp_level), trim(err))

            natts = size(output%variables(i)%attributes_names)
            do j=1,natts
                if (output%variables(i)%attributes_names(j) /= kNULL_CHAR) then
                    ! only put character attributes for now
                    call check(nf90_put_att(ncid, varid,
                                            output%variables(i)%attributes_names(j),
                                            output%variables(i)%attributes_values(j)), trim(err))
                endif
            enddo

            call check(nf90_put_att(ncid, varid, "_FillValue", options%fill_val), trim(err))

            ! Add coordinates attribute
            if (options%coords_2d) then
                call check(nf90_put_att(ncid, varid, coordinates, coords), trim(err))
            endif

            ! errors
            call check(nf90_def_var(ncid, trim(output%variables(i)%name)//"_error",
                                    NF90_REAL, dimids, varid), trim(err))
            call check(nf90_def_var_deflate(ncid, varid,
                       shuffle=options%comp_level, deflate=options%comp_level,
                       deflate_level=options%comp_level), trim(err))
            call check(nf90_put_att(ncid, varid, long_name, trim(output%variables(i)%name)//"_error"), trim(err))
            call check(nf90_put_att(ncid, varid, comment, "error term"), trim(err))
            call check(nf90_put_att(ncid, varid, "transform", options%input_transformations(i)), trim(err))

            ! logistic
            if (options%logistic_threshold /= kFILL_VALUE) then
                call check(nf90_def_var(ncid, trim(output%variables(i)%name)//"_exceedence_probability",
                                        NF90_REAL, dimids, varid), trim(err))
                call check(nf90_def_var_deflate(ncid, varid,
                           shuffle=options%comp_level, deflate=options%comp_level,
                           deflate_level=options%comp_level), trim(err))
                call check(nf90_put_att(ncid, varid, long_name, trim(output%variables(i)%name)//"_exceedence_probability"), trim(err))
                call check(nf90_put_att(ncid, varid, units, "1"), trim(err))
                call check(nf90_put_att(ncid, varid, comment, "logistic probability term"), trim(err))
            endif
        enddo

        ! Add global attributes
        call check(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "1.6"), trim(err))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "title", "GARD (Generalized Analog Regression Downscaling)"), trim(err))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "history", "Created:"//todays_date_time//UTCoffset)//" by "//username//" on "//trim(hostname), trim(err))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "institution", "National Center for Atmospheric Research"), trim(err))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "source", "GARD Version "//trim(kVERSION_STRING)), trim(err))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "Git Version", kGIT_VERSION), trim(err))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "references", "E. Gutmann, J Hamman, & A Wood. (2017). NCAR/GARD: v0.4 [Data set]. Zenodo. http://doi.org/10.5281/zenodo.376874"), trim(err))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "contact", "Ethan Gutmann : gutmann@ucar.edu"), trim(err))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "comment", "GARD is still in development, these results are not final."), trim(err))

        call check(nf90_enddef(ncid), trim(err))

        ! Put coordinate data
        call check(nf90_put_var(ncid, lat_varid, obs_data%lat), trim(err))
        call check(nf90_put_var(ncid, lon_varid, obs_data%lon), trim(err))
        ! call check(nf90_put_var(ncid, time_varid, time_data), trim(err))
        ! TODO: put time/lat/lon vars

        ! mask variable
        call check(nf90_put_var(ncid, mask_varid, obs_data%mask), trim(err))

        do i=1,nvars
            call check(nf90_inq_varid(ncid, trim(output%variables(i)%name), varid), trim(err))
            call check(nf90_put_var(ncid, varid, output_data), trim(err))

            ! errors
            call check(nf90_inq_varid(ncid, trim(output%variables(i)%name)//"_error", varid), trim(err))
            call check(nf90_put_var(ncid, varid, output_data), trim(err))

            ! logistic
            if (options%logistic_threshold /= kFILL_VALUE) then
                call check(nf90_inq_varid(ncid, trim(output%variables(i)%name)//"_exceedence_probability", varid), trim(err))
                call check(nf90_put_var(ncid, varid, output_data), trim(err))
            endif
        enddo

    end subroutine write_output

    ! for whatever reason reshape(data, [nx, ny, nt], order=[2,3,1]) doesn't work!
    subroutine shift_z_dim_3d(input, output)
        implicit none
        real, intent(in), dimension(:,:,:) :: input
        real, intent(inout), dimension(:,:,:), allocatable :: output
        integer :: i,j,k
        integer :: nx,ny,nt
        integer :: Mem_Error

        nt = size(input, 1)
        nx = size(input, 2)
        ny = size(input, 3)

        if ((size(output,1) /= nx).or.(size(output,2) /= ny).or.(size(output,3) /= nt)) then
            deallocate(output)
            allocate(output(nx,ny,nt), stat=Mem_Error)
            if (Mem_Error /= 0) call memory_error(Mem_Error, "output (shift_z)", [nx,ny,nt])
        endif

        do j=1,ny
            do i=1,nx
                output(i,j,:) = input(:,i,j)
            enddo
        enddo

    end subroutine shift_z_dim_3d

    subroutine shift_z_dim_4d(input, output)
        implicit none
        real, intent(in), dimension(:,:,:,:) :: input
        real, intent(inout), dimension(:,:,:,:), allocatable :: output
        integer :: i,j,k
        integer :: nx,ny,nt, nv
        integer :: Mem_Error

        nt = size(input, 1)
        nx = size(input, 2)
        ny = size(input, 3)
        nv = size(input, 4)

        if ((size(output,1) /= nx).or.(size(output,2) /= ny).or.(size(output,3) /= nt).or.(size(output,4) /= nv)) then
            deallocate(output)
            allocate(output(nx,ny,nt,nv), stat=Mem_Error)
            if (Mem_Error /= 0) call memory_error(Mem_Error, "output_4d (shift_z)", [nx,ny,nt,nv])
        endif

        do k=1,nv
            do j=1,ny
                do i=1,nx
                    output(i,j,1:nt,k) = input(:,i,j,k)
                enddo
            enddo
        enddo

    end subroutine shift_z_dim_4d

    subroutine memory_error(error, variable_name, dims)
        implicit none
        integer,          intent(in)                :: error
        character(len=*), intent(in)                :: variable_name
        integer,          intent(in), dimension(:)  :: dims

        write(*,*) "Error allocating memory for variable: ", trim(variable_name)
        write(*,*) "  ERROR        = ", error
        write(*,*) "  Dimensions   = ", dims

        stop "MEMORY ALLOCATION ERROR"

    end subroutine memory_error


end module output_mod
