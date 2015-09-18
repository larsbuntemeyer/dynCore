
      subroutine nc_file_init( )
c-----------------------------------------------------------------------
c     ... initialize netcdf output file
c-----------------------------------------------------------------------

c---------------------------------------------------------------------
c	... local variables
c---------------------------------------------------------------------
      integer :: m
      integer :: ncid, nc_stat, varid
      integer :: lon_id, lat_id, zm_id, zi_id, time_id
      integer :: dims(4), strt_ndx(4), length(4)
      integer :: cvals(max(moie,moje,moke))

      character(len=132) :: mess

      include 'parorg.h'
      include 'netcdf.inc'

c---------------------------------------------------------------------
c	... only root process
c---------------------------------------------------------------------
      if( myid == 0 ) then
c---------------------------------------------------------------------
c	... open file
c---------------------------------------------------------------------
      nc_stat = nf_create( 'remo_out.nc', nf_clobber, ncid )
      if( nc_stat /= nf_noerr ) then
         write(mess,*) 
     1  myid,': nc_file_init: failed to create remo_out.nc'
         call handle_ncerr
      endif

c---------------------------------------------------------------------
c     	... define the dimensions
c---------------------------------------------------------------------
      nc_stat = nf_def_dim( ncid, 'lon', moie, lon_id )
      if( nc_stat /= nf_noerr ) then
         write(mess,*) 
     1   myid,': nc_file_init: failed to define longitude dimension'
         call handle_ncerr
      endif
      nc_stat = nf_def_dim( ncid, 'lat', moje, lat_id )
      if( nc_stat /= nf_noerr ) then
         write(mess,*) 
     1  myid,': nc_file_init: failed to define latitude dimension'
         call handle_ncerr
      endif
      nc_stat = nf_def_dim( ncid, 'lev', moke, zm_id)
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define zm dimension'
         call handle_ncerr
      endif
      nc_stat = nf_def_dim( ncid, 'levi', moke+1, zi_id)
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define zi dimension'
         call handle_ncerr
      endif
      nc_stat = nf_def_dim( ncid, 'time', nf_unlimited, time_id )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define time dimension'
         call handle_ncerr
      endif

c---------------------------------------------------------------------
c     	... define the coordinate variables
c---------------------------------------------------------------------
      dims(1) = lon_id
      nc_stat = 
     1   nf_def_var( ncid, 'lons', nf_int, 1, dims(1), varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable lons'
         call handle_ncerr
      endif
      dims(1) = lat_id
      nc_stat = 
     1  nf_def_var( ncid, 'lats', nf_int, 1, dims(1), varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable lats'
         call handle_ncerr
      endif
      dims(1) = zm_id
      nc_stat = 
     1  nf_def_var( ncid, 'levs', nf_int, 1, dims(1), varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable levs'
         call handle_ncerr
      endif
      dims(1) = zi_id
      nc_stat = nf_def_var( ncid, 'levis', nf_int, 1,
     1          dims(1), varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable levis'
         call handle_ncerr
      endif
      dims(1) = time_id
      nc_stat = nf_def_var( ncid, 'times', nf_int, 
     1                      1, dims(1), varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable times'
         call handle_ncerr
      endif

c---------------------------------------------------------------------
c     	... define the data variables
c---------------------------------------------------------------------
      dims(1:3) = (/ lon_id, lat_id, time_id /)
      nc_stat = nf_def_var( ncid, 'PS', nf_double, 3, dims(1:3), varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable PS'
         call handle_ncerr
      endif
      nc_stat = 
     1   nf_def_var( ncid, 'APRL', nf_double, 3, dims(1:3), varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable APRL'
         call handle_ncerr
      endif
      nc_stat = 
     1   nf_def_var( ncid, 'TSRF', nf_double, 3, dims(1:3), varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable TSRF'
         call handle_ncerr
      endif
      nc_stat = 
     1   nf_def_var( ncid, 'PINT_SRF', nf_double, 3, dims(1:3), varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable PINT_SRF'
         call handle_ncerr
      endif

      dims(:) = (/ lon_id, lat_id, zm_id, time_id /)
      nc_stat = nf_def_var( ncid, 'T', nf_double, 4, dims, varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable T'
         call handle_ncerr
      endif
      nc_stat = nf_def_var( ncid, 'DWDT', nf_double, 4, dims, varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable DWDT'
         call handle_ncerr
      endif

      dims(:) = (/ lon_id, lat_id, zi_id, time_id /)
      nc_stat = nf_def_var( ncid, 'PINT', nf_double, 4, dims, varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable PINT'
         call handle_ncerr
      endif
      nc_stat = nf_def_var( ncid, 'W', nf_double, 4, dims, varid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to define variable W'
         call handle_ncerr
      endif
c---------------------------------------------------------------------
c     	... leave define mode
c---------------------------------------------------------------------
      nc_stat = nf_enddef( ncid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to leave define mode'
         call handle_ncerr
      endif
c---------------------------------------------------------------------
c     	... set and write the coordinate variables
c---------------------------------------------------------------------
      nc_stat = nf_inq_varid( ncid, 'lons', varid )
      cvals(:moie) = (/ (m,m=1,moie) /)
      nc_stat = nf_put_var_int( ncid, varid, cvals(:moie) )
      nc_stat = nf_inq_varid( ncid, 'lats', varid )
      cvals(:moje) = (/ (m,m=1,moje) /)
      nc_stat = nf_put_var_int( ncid, varid, cvals(:moje) )
      nc_stat = nf_inq_varid( ncid, 'levs', varid )
      cvals(:moke) = (/ (m,m=1,moke) /)
      nc_stat = nf_put_var_int( ncid, varid, cvals(:moke) )
      nc_stat = nf_inq_varid( ncid, 'levis', varid )
      cvals(:moke) = (/ (m,m=1,moke+1) /)
      nc_stat = nf_put_var_int( ncid, varid, cvals(:moke+1) )
c---------------------------------------------------------------------
c     	... close the file
c---------------------------------------------------------------------
      nc_stat = nf_close( ncid )
      if( nc_stat /= nf_noerr ) then
         mess = 'nc_file_init: failed to close remo_out.nc'
         call handle_ncerr
      endif

      endif

      CONTAINS

      subroutine handle_ncerr
c---------------------------------------------------------------------
c	... netcdf error handling routine
c---------------------------------------------------------------------

      write(*,*) 'handle_ncerr: ' // trim(mess)
      write(*,*) nf_strerror( ret )
      call mpi_finalize( ret )
      stop 'netcdf error'

      end subroutine handle_ncerr

      end subroutine nc_file_init
