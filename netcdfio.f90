      MODULE netcdfio

      use netcdf
      implicit none
      include 'netcdf.inc'

      real(8), parameter :: doubleval = -99999.d0



      CONTAINS



      SUBROUTINE read_netcdf_hotfile(the_meshdata, the_hotdata,        &
     &    global_fort14_dir)

      use adcirpolate, only: allocate_hotdata, hotdata,  meshdata

      implicit none

      type(meshdata), intent(in)   :: the_meshdata
      type(hotdata),  intent(out)  :: the_hotdata
      character(len=*), intent(in) :: global_fort14_dir

      integer :: ncdim_nele
      integer :: ncdim_node
      integer :: ncid
      integer :: ncvar_iestp
      integer :: ncvar_igep
      integer :: ncvar_igpp
      integer :: ncvar_igvp
      integer :: ncvar_igwp
      integer :: ncvar_imhs
      integer :: ncvar_ipstp
      integer :: ncvar_iths
      integer :: ncvar_ivstp
      integer :: ncvar_iwstp
      integer :: ncvar_nodecode
      integer :: ncvar_noff
      integer :: ncvar_nscoue
      integer :: ncvar_nscouge
      integer :: ncvar_nscougv
      integer :: ncvar_nscougw
      integer :: ncvar_nscoum
      integer :: ncvar_nscouv
      integer :: ncvar_time
      integer :: ncvar_uvel
      integer :: ncvar_vvel
      integer :: ncvar_zeta1
      integer :: ncvar_zeta2
      integer :: ncvar_zetad

      integer,allocatable :: getIntEl(:)
      integer,allocatable :: getIntNd(:)
      integer :: ie
      integer :: iret
      integer :: iv
      integer :: numGlobalEl
      integer :: numGlobalNd

      real(8),allocatable :: getRealNd(:)
      real(8) :: getArray(1)

      !
      ! ALLOCATE ARRAYS
      !

      call allocate_hotdata(the_hotdata, the_meshdata)

      ! INITIALIZE

      ! Open the file.
      iret = nf90_open( trim(global_fort14_dir)//'/fort.67.nc',        &
                        nf90_nowrite, ncid )

      !
      ! READ DIMENSIONS
      !

      ! node

      iret = nf90_inq_dimid( ncid, "node", ncdim_node )
      iret = nf90_inquire_dimension( ncid, ncdim_node, len=numGlobalNd ) 
      allocate(getIntNd(1:numGlobalNd))
      allocate(getRealNd(1:numGlobalNd))

      ! nele

      iret = nf90_inq_dimid( ncid, "nele", ncdim_nele )
      iret = nf90_inquire_dimension( ncid, ncdim_nele, len=numGlobalEl ) 
      allocate(getIntEl(1:numGlobalEl))

      !
      ! READ VARIABLES
      !

      ! time

      iret = nf90_inq_varid( ncid, "time", ncvar_time )
      iret = nf90_get_var( ncid, ncvar_time, getArray )
      the_hotdata%TimeLoc = getArray(1)

      ! zeta1

      iret = nf90_inq_varid( ncid, "zeta1", ncvar_zeta1 )
      iret = nf90_get_var( ncid, ncvar_zeta1, getRealNd,               &
     &                     start=(/1,1/),count=(/numGlobalNd,1/) )
      do iv=1,the_meshdata%NumNd
        the_hotdata%ETA1(iv) = getRealNd(the_meshdata%NdIDs(iv))
      enddo

      ! zeta2

      iret = nf90_inq_varid( ncid, "zeta2", ncvar_zeta2 )
      iret = nf90_get_var( ncid, ncvar_zeta2, getRealNd,               &
     &                     start=(/1,1/),count=(/numGlobalNd,1/) )
      do iv=1,the_meshdata%NumNd
        the_hotdata%ETA2(iv) = getRealNd(the_meshdata%NdIDs(iv))
      enddo

      ! zeta2

      iret = nf90_inq_varid( ncid, "zetad", ncvar_zetad )
      iret = nf90_get_var( ncid, ncvar_zetad, getRealNd,               &
     &                     start=(/1,1/),count=(/numGlobalNd,1/) )
      do iv=1,the_meshdata%NumNd
        the_hotdata%ETADisc(iv) = getRealNd(the_meshdata%NdIDs(iv))
      enddo

      ! u-vel

      iret = nf90_inq_varid( ncid, "u-vel", ncvar_uvel )
      iret = nf90_get_var( ncid, ncvar_uvel, getRealNd,                &
     &                     start=(/1,1/),count=(/numGlobalNd,1/) )
      do iv=1,the_meshdata%NumNd
        the_hotdata%UU2(iv) = getRealNd(the_meshdata%NdIDs(iv))
      enddo

      ! v-vel

      iret = nf90_inq_varid( ncid, "v-vel", ncvar_vvel )
      iret = nf90_get_var( ncid, ncvar_vvel, getRealNd,                &
     &                     start=(/1,1/),count=(/numGlobalNd,1/) )
      do iv=1,the_meshdata%NumNd
        the_hotdata%VV2(iv) = getRealNd(the_meshdata%NdIDs(iv))
      enddo

      ! nodecode

      iret = nf90_inq_varid( ncid, "nodecode", ncvar_nodecode )
      iret = nf90_get_var( ncid, ncvar_nodecode, getIntNd,             &
     &                     start=(/1,1/),count=(/numGlobalNd,1/) )
      do iv=1,the_meshdata%NumNd
        the_hotdata%NNODECODE(iv) = getIntNd(the_meshdata%NdIDs(iv))
        the_hotdata%realNODECODE(iv) = dble(the_hotdata%NNODECODE(iv))
      enddo

      ! noff

      iret = nf90_inq_varid( ncid, "noff", ncvar_noff )
      iret = nf90_get_var( ncid, ncvar_noff, getIntEl,                 &
     &                     start=(/1,1/),count=(/numGlobalEl,1/) )
      do ie=1,the_meshdata%NumEl
        the_hotdata%NOFF(ie) = getIntEl(the_meshdata%ElIDs(ie))
      enddo

      ! imhs

      iret = nf90_inq_varid( ncid, "imhs", ncvar_imhs )
      iret = nf90_get_var( ncid, ncvar_imhs, the_hotdata%IMHS )

      ! iths

      iret = nf90_inq_varid( ncid, "iths", ncvar_iths )
      iret = nf90_get_var( ncid, ncvar_iths, the_hotdata%ITHS )

      ! iestp

      iret = nf90_inq_varid( ncid, "iestp", ncvar_iestp )
      iret = nf90_get_var( ncid, ncvar_iestp, the_hotdata%IESTP )

      ! nscoue

      iret = nf90_inq_varid( ncid, "nscoue", ncvar_nscoue )
      iret = nf90_get_var( ncid, ncvar_nscoue, the_hotdata%NSCOUE )

      ! ivstp

      iret = nf90_inq_varid( ncid, "ivstp", ncvar_ivstp )
      iret = nf90_get_var( ncid, ncvar_ivstp, the_hotdata%IVSTP )

      ! nscouv

      iret = nf90_inq_varid( ncid, "nscouv", ncvar_nscouv )
      iret = nf90_get_var( ncid, ncvar_nscouv, the_hotdata%NSCOUV )

      ! ipstp

      iret = nf90_inq_varid( ncid, "ipstp", ncvar_ipstp )
      iret = nf90_get_var( ncid, ncvar_ipstp, the_hotdata%IPSTP )

      ! iwstp

      iret = nf90_inq_varid( ncid, "iwstp", ncvar_iwstp )
      iret = nf90_get_var( ncid, ncvar_iwstp, the_hotdata%IWSTP )

      ! nscoum

      iret = nf90_inq_varid( ncid, "nscoum", ncvar_nscoum )
      iret = nf90_get_var( ncid, ncvar_nscoum, the_hotdata%NSCOUM )

      ! igep

      iret = nf90_inq_varid( ncid, "igep", ncvar_igep )
      iret = nf90_get_var( ncid, ncvar_igep, the_hotdata%IGEP )

      ! nscouge

      iret = nf90_inq_varid( ncid, "nscouge", ncvar_nscouge )
      iret = nf90_get_var( ncid, ncvar_nscouge, the_hotdata%NSCOUGE )

      ! igvp

      iret = nf90_inq_varid( ncid, "igvp", ncvar_igvp )
      iret = nf90_get_var( ncid, ncvar_igvp, the_hotdata%IGVP )

      ! nscougv

      iret = nf90_inq_varid( ncid, "nscougv", ncvar_nscougv )
      iret = nf90_get_var( ncid, ncvar_nscougv, the_hotdata%NSCOUGV )

      ! igpp

      iret = nf90_inq_varid( ncid, "igpp", ncvar_igpp )
      iret = nf90_get_var( ncid, ncvar_igpp, the_hotdata%IGPP )

      ! igwp

      iret = nf90_inq_varid( ncid, "igwp", ncvar_igwp )
      iret = nf90_get_var( ncid, ncvar_igwp, the_hotdata%IGWP )

      ! nscougw

      iret = nf90_inq_varid( ncid, "nscougw", ncvar_nscougw )
      iret = nf90_get_var( ncid, ncvar_nscougw, the_hotdata%NSCOUGW )

      !
      ! FINALIZE
      !

      !
      ! Close the file.
      !
      iret = nf90_close( ncid )

      ENDSUBROUTINE



      SUBROUTINE write_netcdf_hotfile(the_meshdata, the_hotdata)

      use adcirpolate, only: meshdata, hotdata

      implicit none

      type(hotdata),  intent(in) :: the_hotdata
      type(meshdata), intent(in) :: the_meshdata

      character(len=1000) :: variable_name
      character(len=1000) :: ncatts(10,2)

      integer :: ncdim_nele
      integer :: ncdim_node
      integer :: ncdim_time
      integer :: ncid
      integer :: ncvar_iestp
      integer :: ncvar_igep
      integer :: ncvar_igpp
      integer :: ncvar_igvp
      integer :: ncvar_igwp
      integer :: ncvar_imhs
      integer :: ncvar_ipstp
      integer :: ncvar_iths
      integer :: ncvar_ivstp
      integer :: ncvar_iwstp
      integer :: ncvar_nodecode
      integer :: ncvar_noff
      integer :: ncvar_nscoue
      integer :: ncvar_nscouge
      integer :: ncvar_nscougv
      integer :: ncvar_nscougw
      integer :: ncvar_nscoum
      integer :: ncvar_nscouv
      integer :: ncvar_time
      integer :: ncvar_uvel
      integer :: ncvar_vvel
      integer :: ncvar_zeta1
      integer :: ncvar_zeta2
      integer :: ncvar_zetad

      integer :: iret

      real(8) :: putArray(1)

      ! INITIALIZE

      ! Open the file.
      iret = nf90_create( 'fine/fort.67.nc',                           &
                          ior(nf90_classic_model,nf90_hdf5),           &
                          ncid )

      ! DIMENSIONS

      !
      ! time -- read in netcdfio/getDimensions
      !
      iret = nf90_def_dim( ncid, 'time', nf_unlimited, ncdim_time )
      !
      ! node -- read in netcdfio/getDimensions
      !
      iret = nf90_def_dim( ncid, 'node', the_meshdata%NumNd,           &
                           ncdim_node )
      !
      ! nele -- read in netcdfio/getDimensions
      !
      iret = nf90_def_dim( ncid, 'nele', the_meshdata%NumEl,           &
                           ncdim_nele )

      ! VARIABLES & ATTRIBUTES

      ncatts = " "

      !
      ! time -- read in netcdfio/readNetCDFHotstart
      !
      write(variable_name,'(a)') "time"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "model time"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "time"
      write(ncatts(3,1),'(a)') "units"
      write(ncatts(3,2),'(a)') "seconds since"                         &
                               //" 2000-01-01 00:00:00 UTC"
      write(ncatts(4,1),'(a)') "base_date"
      write(ncatts(4,2),'(a)') "2000-01-01 00:00:00 UTC"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_time, (/ ncdim_time, 0 /),                 &
                .FALSE., .FALSE., ncatts )
      !
      ! zeta1
      !
      write(variable_name,'(a)') "zeta1"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "water surface elevation at"            &
                               //" previous time step"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "water surface elevation at"            &
                               //" previous time step"
      write(ncatts(3,1),'(a)') "units"
      write(ncatts(3,2),'(a)') "m"
      write(ncatts(4,1),'(a)') "_FillValue"
      write(ncatts(4,2),'(a)') " "
      write(ncatts(5,1),'(a)') "positive"
      write(ncatts(5,2),'(a)') "up"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_zeta1, (/ ncdim_node, ncdim_time /),       &
                .FALSE., .TRUE., ncatts )
      !
      ! zeta2
      !
      write(variable_name,'(a)') "zeta2"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "water surface elevation at"            &
                               //" current time step"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "water surface elevation at"            &
                               //" current time step"
      write(ncatts(3,1),'(a)') "units"
      write(ncatts(3,2),'(a)') "m"
      write(ncatts(4,1),'(a)') "_FillValue"
      write(ncatts(4,2),'(a)') " "
      write(ncatts(5,1),'(a)') "positive"
      write(ncatts(5,2),'(a)') "up"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_zeta2, (/ ncdim_node, ncdim_time /),       &
                .FALSE., .TRUE., ncatts )
      !
      ! zetad
      !
      write(variable_name,'(a)') "zetad"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "water elevation at"                    &
                               //" flux specified boundary"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "water elevation at"                    &
                               //" flux specified boundary"
      write(ncatts(3,1),'(a)') "units"
      write(ncatts(3,2),'(a)') "m"
      write(ncatts(4,1),'(a)') "_FillValue"
      write(ncatts(4,2),'(a)') " "
      write(ncatts(5,1),'(a)') "positive"
      write(ncatts(5,2),'(a)') "up"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_zetad, (/ ncdim_node, ncdim_time /),       &
                .FALSE., .TRUE., ncatts )
      !
      ! u-vel
      !
      write(variable_name,'(a)') "u-vel"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "vertically averaged e/w velocity"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "u_velocity"
      write(ncatts(3,1),'(a)') "positive"
      write(ncatts(3,2),'(a)') "east"
      write(ncatts(4,1),'(a)') "units"
      write(ncatts(4,2),'(a)') "m s-1"
      write(ncatts(5,1),'(a)') "_FillValue"
      write(ncatts(5,2),'(a)') " "
      write(ncatts(6,1),'(a)') "dry_Value"
      write(ncatts(6,2),'(a)') " "
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_uvel, (/ ncdim_node, ncdim_time /),        &
                .FALSE., .TRUE., ncatts )
      !
      ! v-vel
      !
      write(variable_name,'(a)') "v-vel"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "vertically averaged n/s velocity"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "v_velocity"
      write(ncatts(3,1),'(a)') "positive"
      write(ncatts(3,2),'(a)') "north"
      write(ncatts(4,1),'(a)') "units"
      write(ncatts(4,2),'(a)') "m s-1"
      write(ncatts(5,1),'(a)') "_FillValue"
      write(ncatts(5,2),'(a)') " "
      write(ncatts(6,1),'(a)') "dry_Value"
      write(ncatts(6,2),'(a)') " "
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_vvel, (/ ncdim_node, ncdim_time /),        &
                .FALSE., .TRUE., ncatts )
      !
      ! nodecode
      !
      write(variable_name,'(a)') "nodecode"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "wet or dry state of node where 1"      &
                               //" indicates that the node is wet and" &
                               //" 0 indicates that the node is dry"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "node_wet_or_dry"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_nodecode, (/ ncdim_node, ncdim_time /),    &
                .TRUE., .TRUE., ncatts )
      !
      ! noff
      !
      write(variable_name,'(a)') "noff"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "wet or dry state of element where 1"   &
                               //" indicates that the element is wet"  &
                               //" and 0 indicates that it is dry"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "element_wet_or_dry"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_noff, (/ ncdim_nele, ncdim_time /),        &
                .TRUE., .TRUE., ncatts )
      !
      ! imhs
      !
      write(variable_name,'(a)') "imhs"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "model_type"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "model_type"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_imhs, (/ 0, 0 /),                          &
                .TRUE., .FALSE., ncatts )
      !
      ! iths
      !
      write(variable_name,'(a)') "iths"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "model time step number since the"      &
                               //" beginning of the model run"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "model_time_step"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_iths, (/ 0, 0 /),                          &
                .TRUE., .FALSE., ncatts )
      !
      ! iestp
      !
      write(variable_name,'(a)') "iestp"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "line number (for ASCII output) or"     &
                               //" record number (for binary output)"  &
                               //" of the most recent entry in the"    &
                               //" elevation time series at specified" &
                               //" elevation recording stations output"&
                               //" file"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "line/record_number_of_last_entry_in"   &
                               //"_elev_rec_stations"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_iestp, (/ 0, 0 /),                         &
                .TRUE., .FALSE., ncatts )
      !
      ! nscoue
      !
      write(variable_name,'(a)') "nscoue"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "time step counter to determine when"   &
                               //" the next entry will be written to"  &
                               //" the elevation time series at"       &
                               //" specified elevation recording"      &
                               //" stations output file"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "time_step_counter_for_next_entry_elev" &
                               //"_rec_stations"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_nscoue, (/ 0, 0 /),                        &
                .TRUE., .FALSE., ncatts )
      !
      ! ivstp
      !
      write(variable_name,'(a)') "ivstp"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "line number (for ASCII output) or"     &
                               //" record number (for binary output)"  &
                               //" of the most recent entry in the"    &
                               //" depth-averaged velocity time series"&
                               //" at specified velocity recording"    &
                               //" stations output file"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "record_number_of_last_entry_in_vel_rec"&
                               //"_stations"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_ivstp, (/ 0, 0 /),                         &
                .TRUE., .FALSE., ncatts )
      !
      ! nscouv
      !
      write(variable_name,'(a)') "nscouv"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "time step counter to determine when"   &
                               //" the next entry will be written to"  &
                               //" the depth-averaged velocity time"   &
                               //" series at specified velocity"       &
                               //" recording stations output file"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "time_step_counter_for_next_entry_vel"  &
                               //"_rec_stations"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_nscouv, (/ 0, 0 /),                        &
                .TRUE., .FALSE., ncatts )
      !
      ! ipstp
      !
      write(variable_name,'(a)') "ipstp"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "line number (for ASCII output) or"     &
                               //" record number (for binary output)"  &
                               //" of the most recent entry in the"    &
                               //" atmospheric pressure time series at"&
                               //" specified meteorological recording" &
                               //" stations"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "record_number_of_last_entry_of_atm"    &
                               //"_press_at_rec_stations"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_ipstp, (/ 0, 0 /),                         &
                .TRUE., .FALSE., ncatts )
      !
      ! iwstp
      !
      write(variable_name,'(a)') "iwstp"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "line number (for ASCII output) or"     &
                               //" record number (for binary output)"  &
                               //" of the most recent entry in the"    &
                               //" wind velocity time series at"       &
                               //" specified meteorological recording" &
                               //" stations"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "record_number_of_last_entry_of_wind"   &
                               //"_vel_at_rec_stations"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_iwstp, (/ 0, 0 /),                         &
                .TRUE., .FALSE., ncatts )
      !
      ! nscoum
      !
      write(variable_name,'(a)') "nscoum"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "time step counter to determine when"   &
                               //" the next entry will be written to"  &
                               //" the atmospheric pressure time"      &
                               //" series at specified meteorological" &
                               //" recording stations and wind"        &
                               //" velocity time series at specified"  &
                               //" meteorological recording stations"  &
                               //" output files"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "time_step_counter_of_atm_press_and"    &
                               //"_wind_vel_at_rec_stations"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_nscoum, (/ 0, 0 /),                        &
                .TRUE., .FALSE., ncatts )
      !
      ! igep
      !
      write(variable_name,'(a)') "igep"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "line number (for ASCII output) or"     &
                               //" record number (for binary output)"  &
                               //" of the most recent entry in the"    &
                               //" elevation time series at all nodes" &
                               //" in the model grid output file"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "record_number_of_last_entry_of_elev_at"&
                               //"_model_nodes"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_igep, (/ 0, 0 /),                          &
                .TRUE., .FALSE., ncatts )
      !
      ! nscouge
      !
      write(variable_name,'(a)') "nscouge"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "time step counter to determine when"   &
                               //" the next entry will be written to"  &
                               //" the elevation time series at all"   &
                               //" nodes in the model grid output file"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "time_step_counter_of_elev_at_model"    &
                               //"_nodes"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_nscouge, (/ 0, 0 /),                       &
                .TRUE., .FALSE., ncatts )
      !
      ! igvp
      !
      write(variable_name,'(a)') "igvp"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "line number (for ASCII output) or"     &
                               //" record number (for binary output)"  &
                               //" of the most recent entry in the"    &
                               //" depth-averaged velocity time series"&
                               //" at all nodes in the model grid"     &
                               //" output file"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "record_number_of_last_entry_of_vel_at" &
                               //"_model_nodes"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_igvp, (/ 0, 0 /),                          &
                .TRUE., .FALSE., ncatts )
      !
      ! nscougv
      !
      write(variable_name,'(a)') "nscougv"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "time step counter to determine when"   &
                               //" the next entry will be written to"  &
                               //" the depth-averaged velocity time"   &
                               //" series at all nodes in the model"   &
                               //" grd output file"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "time_step_counter_of_vel_at_model_nodes"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_nscougv, (/ 0, 0 /),                       &
                .TRUE., .FALSE., ncatts )
      !
      ! igpp
      !
      write(variable_name,'(a)') "igpp"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "line number (for ASCII output) or"     &
                               //" record number (for binary output)"  &
                               //" of the most recent entry in the"    &
                               //" atmospheric pressure time series at"&
                               //" all nodes in the model grid output" &
                               //" file" 
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "record_number_of_last_entry_of_atm"    &
                               //"_press_at_model_nodes"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_igpp, (/ 0, 0 /),                          &
                .TRUE., .FALSE., ncatts )
      !
      ! igwp
      !
      write(variable_name,'(a)') "igwp"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "line number (for ASCII output) or"     &
                               //" record number (for binary output)"  &
                               //" of the most recent entry in the"    &
                               //" wind stress or velocity time series"&
                               //" at all nodes in the model grid"     &
                               //" output file" 
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "record_number_of_last_entry_of_wind"   &
                               //"_vel_at_model_nodes"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_igwp, (/ 0, 0 /),                          &
                .TRUE., .FALSE., ncatts )
      !
      ! nscougw
      !
      write(variable_name,'(a)') "nscougw"
      write(ncatts(1,1),'(a)') "long_name"
      write(ncatts(1,2),'(a)') "time step counter to determine when"   &
                               //" the next entry will be written to"  &
                               //" the atmospheric pressure time"      &
                               //" series at all nodes in the model"   &
                               //" grid and wind stress or velocity"   &
                               //" time series at all nodes in the"    &
                               //" model grid output files"
      write(ncatts(2,1),'(a)') "standard_name"
      write(ncatts(2,2),'(a)') "time_step_counter_of_atm_press_and"    &
                               //"_wind_vel_at_model_nodes"
      call nf90_define_variable( variable_name,                        &
                ncid, ncvar_nscougw, (/ 0, 0 /),                       &
                .TRUE., .FALSE., ncatts )

      !
      ! End the variable definitions & attributes.
      !
      iret = nf90_enddef( ncid )

      !
      ! PUT VARIABLES
      !

      !
      ! time -- read in netcdfio/readNetCDFHotstart
      !
      putArray(1) = the_hotdata%TimeLoc
      iret = nf90_put_var( ncid, ncvar_time, putArray,                 &
             start=(/1/), count=(/1/) )
      !
      ! zeta1 -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_zeta1, the_hotdata%ETA1,        &
             start=(/1,1/), count=(/the_meshdata%NumNd,1/) )
      !
      ! zeta2 -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_zeta2, the_hotdata%ETA2,        &
             start=(/1,1/), count=(/the_meshdata%NumNd,1/) )
      !
      ! zetad -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_zetad, the_hotdata%ETADisc,     &
             start=(/1,1/), count=(/the_meshdata%NumNd,1/) )
      !
      ! u-vel -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_uvel, the_hotdata%UU2,          &
             start=(/1,1/), count=(/the_meshdata%NumNd,1/) )
      !
      ! v-vel -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_vvel, the_hotdata%VV2,          &
             start=(/1,1/), count=(/the_meshdata%NumNd,1/) )
      !
      ! nodecode -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_nodecode, the_hotdata%NNODECODE,&
             start=(/1,1/), count=(/the_meshdata%NumNd,1/) )
      !
      ! noff -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_noff, the_hotdata%NOFF,         &
             start=(/1,1/), count=(/the_meshdata%NumEl,1/) )
      !
      ! imhs -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_imhs, the_hotdata%IMHS ) 
      !
      ! iths -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_iths, the_hotdata%ITHS ) 
      !
      ! iestp -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_iestp, the_hotdata%IESTP ) 
      !
      ! nscoue -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_nscoue, the_hotdata%NSCOUE ) 
      !
      ! ivstp -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_ivstp, the_hotdata%IVSTP ) 
      !
      ! nscouv -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_nscouv, the_hotdata%NSCOUV ) 
      !
      ! ipstp -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_ipstp, the_hotdata%IPSTP ) 
      !
      ! iwstp -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_iwstp, the_hotdata%IWSTP ) 
      !
      ! nscoum -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_nscoum, the_hotdata%NSCOUM ) 
      !
      ! igep -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_igep, the_hotdata%IGEP ) 
      !
      ! nscouge -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_nscouge, the_hotdata%NSCOUGE ) 
      !
      ! igvp -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_igvp, the_hotdata%IGVP ) 
      !
      ! nscougv -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_nscougv, the_hotdata%NSCOUGV ) 
      !
      ! igpp -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_igpp, the_hotdata%IGPP ) 
      !
      ! igwp -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_igwp, the_hotdata%IGWP ) 
      !
      ! nscougw -- read in netcdfio/readNetCDFHotstart
      !
      iret = nf90_put_var( ncid, ncvar_nscougw, the_hotdata%NSCOUGW ) 

      !
      ! FINALIZE
      !

      !
      ! Close the file.
      !
      iret = nf90_close( ncid )

      ENDSUBROUTINE write_netcdf_hotfile



      SUBROUTINE nf90_define_variable( variable_name,                  &
                                       ncid,                           &
                                       ncvar,                          &
                                       ncdims,                         &
                                       whether_int,                    &
                                       deflate,                        &
                                       ncatts )

      IMPLICIT NONE

      character(len=1000), intent(in)    :: variable_name
      integer,             intent(in)    :: ncid
      integer,             intent(inout) :: ncvar
      integer,             intent(in)    :: ncdims(2)
      logical,             intent(in)    :: whether_int
      logical,             intent(in)    :: deflate
      character(len=1000), intent(inout) :: ncatts(10,2)

      integer :: ia
      integer :: iatt
      integer :: iret
      integer,allocatable :: ncdim(:)

      if(ncdims(1).eq.0)then
        allocate(ncdim(0))
      elseif(ncdims(2).eq.0)then
        allocate(ncdim(1))
        ncdim(1) = ncdims(1)
      else
        allocate(ncdim(2))
        ncdim = ncdims
      endif

      if(whether_int)then
        iret = nf90_def_var( ncid, trim(adjustl(variable_name)),       &
                             nf90_int, ncdim , varid=ncvar )
      else
        iret = nf90_def_var( ncid, trim(adjustl(variable_name)),       &
                             nf90_double, ncdim , varid=ncvar )
      endif
      if(deflate)then
        iret = nf90_def_var_deflate( ncid, ncvar, 1, 1, 2 )
      endif
      do ia=1,10
        if(len_trim(ncatts(ia,1)).gt.0)then
          if( (index(ncatts(ia,1),'start_index').gt.0) .or.            &
              (index(ncatts(ia,1),'topology_dimension').gt.0) )then
            read(ncatts(ia,2),'(i)') iatt
            iret = nf90_put_att( ncid, ncvar,                          &
                        trim(adjustl(ncatts(ia,1))),                   &
                        iatt )
          elseif( (index(ncatts(ia,1),'_FillValue').gt.0).or.          &
                  (index(ncatts(ia,1),'dry_Value').gt.0) )then
            iret = nf90_put_att( ncid, ncvar,                          &
                        trim(adjustl(ncatts(ia,1))),                   &
                        doubleval )
          else
            iret = nf90_put_att( ncid, ncvar,                          &
                        trim(adjustl(ncatts(ia,1))),                   &
                        trim(adjustl(ncatts(ia,2))) )
          endif
        endif
      enddo
      ncatts = " "

      deallocate(ncdim)

      ENDSUBROUTINE


      ENDMODULE
