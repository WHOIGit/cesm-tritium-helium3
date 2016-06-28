!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module tr3he_mod

!BOP
! !MODULE: tr3he_mod
!
!  Module for tritium (3H) and helium 3 (3He) 
!
!  Created by Ivan Lima <ivan@whoi.edu> on Tue Nov 24 2015 10:39:33 -0500
!  Last modified on Wed May 11 2016 17:20:01 -0400
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod

   use kinds_mod
   use blocks, only: nx_block, ny_block, block
   use domain_size, only: max_blocks_clinic, km
   use domain, only: nblocks_clinic, distrb_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use constants
   use io_types, only: stdout
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field
   use passive_tracer_tools, only: forcing_monthly_every_ts, &
       ind_name_pair, tracer_read, read_field
   use broadcast
   use netcdf

   implicit none
   save

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   private

! !PUBLIC MEMBER FUNCTIONS:

   public ::               &
       tr3he_tracer_cnt,   &
       tr3he_init,         &
       tr3he_set_interior, &
       tr3he_set_sflux,    &
       tr3he_tavg_forcing

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracers
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
       tr3he_tracer_cnt = 2

!-----------------------------------------------------------------------
!  flags controlling which portion of code are executed
!-----------------------------------------------------------------------

  logical (log_kind) :: &
     lflux_gas_3he      ! gas flux for helium 3

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      tritium_ind =  1,  & ! tritium
      helium3_ind =  2     ! helium 3

!-----------------------------------------------------------------------
!  gas flux parameters
!-----------------------------------------------------------------------

   real (r8), parameter ::  &
      R = 8.314425_r8     ! ideal gas constant in (m^3 Pa)/(K mol)

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(tr3he_tracer_cnt) :: &
      ind_name_table = (/ &
      ind_name_pair(tritium_ind, 'TRITIUM'), &
      ind_name_pair(helium3_ind, 'HELIUM3') /)

!-----------------------------------------------------------------------
!  mask that eases avoidance of computation over land
!-----------------------------------------------------------------------

   logical (log_kind), dimension(:,:,:), allocatable :: &
      LAND_MASK

!-----------------------------------------------------------------------
!  forcing related variables
!-----------------------------------------------------------------------

   character(char_len) :: &
      tr3he_formulation,     & ! how to calculate flux (file or model)
      tritium_file             ! filename for atm tritium

   integer (int_kind) ::  &
      model_year,          & ! arbitrary model year
      data_year              ! year in data that corresponds to model_year

   real (r8), dimension(:,:,:,:), allocatable :: &
      INTERP_WORK            ! temp array for interpolate_forcing output

   type(forcing_monthly_every_ts) :: &
      fice_file,           & ! ice fraction, if read from file
      xkw_file,            & ! a * wind-speed ** 2, if read from file
      ap_file                ! atmoshperic pressure, if read from file

!-----------------------------------------------------------------------
!  define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable ::   &
      tr3he_SFLUX_TAVG

   integer (int_kind) :: &
      tavg_3He_IFRAC,      bufind_3He_IFRAC,      & ! tavg id for ice fraction
      tavg_3He_XKW,        bufind_3He_XKW,        & ! tavg id for xkw
      tavg_3He_ATM_PRESS,  bufind_3He_ATM_PRESS,  & ! tavg id for atmospheric pressure
      tavg_3He_SCHMIDT,    bufind_3He_SCHMIDT,    & ! tavg id for 3He Schmidt number
      tavg_3He_PV,         bufind_3He_PV,         & ! tavg id for 3He piston velocity
      tavg_3He_SURF_SAT,   bufind_3He_SURF_SAT,   & ! tavg id for 3He surface saturation
      tavg_3He_GAS_FLUX,   buf_ind_3He_GAS_FLUX     ! 3He gas flux

!-----------------------------------------------------------------------
!  define tavg id for nonstandard 3d fields
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_TRITIUM_DECAY        ! tavg id for tritium decay

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: tr3he_sflux_timer

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: tr3he_init
! !INTERFACE:

 subroutine tr3he_init(init_ts_file_fmt, read_restart_filename, &
                     tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize tr3he module. This involves setting metadata, reading
!  the modules namelist and setting initial conditions.

! !REVISION HISTORY:
!  same as module

! !USES:

   use prognostic, only: curtime, oldtime
   use grid, only: KMT, n_topo_smooth, fill_points
   use grid, only: REGION_MASK
   use io_types, only: nml_in, nml_filename
   use prognostic, only: tracer_field
   use timers, only: get_timer
   use passive_tracer_tools, only: init_forcing_monthly_every_ts, &
       rest_read_tracer_block, file_read_tracer_block

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(tr3he_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real (r8), dimension(nx_block,ny_block,km,tr3he_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: sub_name = 'tr3he_mod:tr3he_init'

   character(char_len) ::       &
      init_tr3he_option,        & ! option for initialization of bgc
      init_tr3he_init_file,     & ! filename for option 'file'
      init_tr3he_init_file_fmt    ! file format for option 'file'

   integer (int_kind) ::      &
      n,                      & ! index for looping over tracers
      k,                      & ! index for looping over depth levels
      iblock,                 & ! index for looping over blocks
      nml_error                 ! namelist i/o error flag

   type(tracer_read), dimension(tr3he_tracer_cnt) :: &
      tracer_init_ext           ! namelist variable for initializing tracers

   type(tracer_read) ::       &
      gas_flux_fice,          & ! ice fraction for gas fluxes
      gas_flux_ws,            & ! wind speed for gas fluxes
      gas_flux_ap               ! atmospheric pressure for gas fluxes

   namelist /tr3he_nml/ &
      init_tr3he_option, init_tr3he_init_file, init_tr3he_init_file_fmt, &
      tracer_init_ext, tritium_file, model_year, data_year, &
      tr3he_formulation, gas_flux_fice, gas_flux_ws, gas_flux_ap

   character (char_len) ::  &
      tr3he_restart_filename      ! modified file name for restart file

!-----------------------------------------------------------------------
!  initialize forcing_monthly_every_ts variables
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call init_forcing_monthly_every_ts(fice_file)
   call init_forcing_monthly_every_ts(xkw_file)
   call init_forcing_monthly_every_ts(ap_file)

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   tracer_d_module(tritium_ind)%long_name  = 'tritium'
   tracer_d_module(helium3_ind)%long_name  = 'helium 3'
   do n = 1, tr3he_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%units      = 'pmol/m^3'
      tracer_d_module(n)%tend_units = 'pmol/m^3'
      tracer_d_module(n)%flux_units = 'pmol/m^2/s'
   end do

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_tr3he_option        = 'unknown'
   init_tr3he_init_file     = 'unknown'
   init_tr3he_init_file_fmt = 'bin'

   do n = 1, tr3he_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   tritium_file      = 'unknown'
   model_year        = 1
   data_year         = 1931
   tr3he_formulation = 'model'

   gas_flux_fice%filename     = 'unknown'
   gas_flux_fice%file_varname = 'FICE'
   gas_flux_fice%scale_factor = c1
   gas_flux_fice%default_val  = c0
   gas_flux_fice%file_fmt     = 'bin'

   gas_flux_ws%filename     = 'unknown'
   gas_flux_ws%file_varname = 'XKW'
   gas_flux_ws%scale_factor = c1
   gas_flux_ws%default_val  = c0
   gas_flux_ws%file_fmt     = 'bin'

   gas_flux_ap%filename     = 'unknown'
   gas_flux_ap%file_varname = 'P'
   gas_flux_ap%scale_factor = c1
   gas_flux_ap%default_val  = c0
   gas_flux_ap%file_fmt     = 'bin'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=tr3he_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(sub_name, 'tr3he_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ sub_name)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_tr3he_option, master_task)
   call broadcast_scalar(init_tr3he_init_file, master_task)
   call broadcast_scalar(init_tr3he_init_file_fmt, master_task)

   do n = 1, tr3he_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

   call broadcast_scalar(tritium_file, master_task)
   call broadcast_scalar(model_year, master_task)
   call broadcast_scalar(data_year, master_task)
   call broadcast_scalar(tr3he_formulation, master_task)

   call broadcast_scalar(gas_flux_fice%filename, master_task)
   call broadcast_scalar(gas_flux_fice%file_varname, master_task)
   call broadcast_scalar(gas_flux_fice%scale_factor, master_task)
   call broadcast_scalar(gas_flux_fice%default_val, master_task)
   call broadcast_scalar(gas_flux_fice%file_fmt, master_task)

   fice_file%input = gas_flux_fice

   call broadcast_scalar(gas_flux_ws%filename, master_task)
   call broadcast_scalar(gas_flux_ws%file_varname, master_task)
   call broadcast_scalar(gas_flux_ws%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ws%default_val, master_task)
   call broadcast_scalar(gas_flux_ws%file_fmt, master_task)

   xkw_file%input = gas_flux_ws

   call broadcast_scalar(gas_flux_ap%filename, master_task)
   call broadcast_scalar(gas_flux_ap%file_varname, master_task)
   call broadcast_scalar(gas_flux_ap%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ap%default_val, master_task)
   call broadcast_scalar(gas_flux_ap%file_fmt, master_task)

   ap_file%input = gas_flux_ap

!-----------------------------------------------------------------------
!   initialize tracers
!-----------------------------------------------------------------------

   select case (init_tr3he_option)

   case ('zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d tritium and 3He set to all zeros'
          write(stdout,delim_fmt)
      endif

   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      tr3he_restart_filename = char_blank

      if (init_tr3he_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(sub_name, 'no restart file to read tritium and 3He from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ sub_name)
         endif
         tr3he_restart_filename = read_restart_filename
         init_tr3he_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         tr3he_restart_filename = trim(init_tr3he_init_file)

      endif

      call rest_read_tracer_block(init_tr3he_init_file_fmt, &
                                  tr3he_restart_filename,   &
                                  tracer_d_module,        &
                                  TRACER_MODULE)

   case ('file', 'ccsm_startup')

      call document(sub_name, 'tritium and 3He being read from separate file')

      call file_read_tracer_block(init_tr3he_init_file_fmt, &
                                  init_tr3he_init_file,     &
                                  tracer_d_module,        &
                                  ind_name_table,         &
                                  tracer_init_ext,        &
                                  TRACER_MODULE)

      if (n_topo_smooth > 0) then
         do n = 1, tr3he_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'tr3he_init: error in fill_points')
                  return
               endif
            end do
         end do
      endif

   case default
      call document(sub_name, 'init_tr3he_option', init_tr3he_option)
      call exit_POP(sigAbort, 'unknown init_tr3he_option')

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock = 1, nblocks_clinic
   do n = 1, tr3he_tracer_cnt
      do k = 1, km
         where (k > KMT(:,:,iblock))
            TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
            TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
         end where
      end do
   end do
   end do

!-----------------------------------------------------------------------
!  allocate and initialize LAND_MASK (true for ocean points)
!-----------------------------------------------------------------------

   allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
   LAND_MASK = (KMT.gt.0)

   call get_timer(tr3he_sflux_timer, 'tr3he_SFLUX', 1, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  call other initialization subroutines
!-----------------------------------------------------------------------

   call tr3he_init_tavg
   call tr3he_init_sflux

!-----------------------------------------------------------------------
!EOC

 end subroutine tr3he_init

!***********************************************************************
!BOP
! !IROUTINE: tr3he_init_tavg
! !INTERFACE:

 subroutine tr3he_init_tavg

! !DESCRIPTION:
!  Define tavg fields not automatically handled by the base model.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      var_cnt             ! how many tavg variables are defined

!-----------------------------------------------------------------------
!  2D fields related to surface fluxes
!-----------------------------------------------------------------------

   var_cnt = 0

   call define_tavg_field(tavg_3He_IFRAC,'HELIUM3_IFRAC',2,           &
                          long_name='Ice fraction for 3He fluxes',&
                          units='fraction', grid_loc='2110')
   var_cnt = var_cnt+1
   bufind_3He_IFRAC = var_cnt

   call define_tavg_field(tavg_3He_XKW,'HELIUM3_XKW',2,               &
                          long_name='XKW for 3He fluxes',         &
                          units='m/s', grid_loc='2110')
   var_cnt = var_cnt+1
   bufind_3He_XKW = var_cnt

   call define_tavg_field(tavg_3He_ATM_PRESS,'HELIUM3_ATM_PRESS',2,           &
                          long_name='Atmospheric pressure for 3He fluxes',&
                          units='Pascals', grid_loc='2110')
   var_cnt = var_cnt+1
   bufind_3He_ATM_PRESS = var_cnt

   call define_tavg_field(tavg_3He_SCHMIDT,'HELIUM3_SCHMIDT',2,     &
                          long_name='3He Schmidt number',       &
                          units='none', grid_loc='2110')
   var_cnt = var_cnt+1
   bufind_3He_SCHMIDT = var_cnt

   call define_tavg_field(tavg_3He_PV,'HELIUM3_PV',2,               &
                          long_name='3He piston velocity',      &
                          units='m/s', grid_loc='2110')
   var_cnt = var_cnt+1
   bufind_3He_PV = var_cnt

   call define_tavg_field(tavg_3He_SURF_SAT,'HELIUM3_SURF_SAT',2,   &
                          long_name='3He saturation',           &
                          units='pmol/m^3', grid_loc='2110')
   var_cnt = var_cnt+1
   bufind_3He_SURF_SAT = var_cnt

   call define_tavg_field(tavg_3He_GAS_FLUX,'HELIUM3_FG',2,   &
                          long_name='3He surface gas flux',           &
                          units='pmol/m^3/s', grid_loc='2110')
   var_cnt = var_cnt+1
   buf_ind_3He_GAS_FLUX = var_cnt

!-----------------------------------------------------------------------

   allocate(tr3he_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   tr3he_SFLUX_TAVG = c0

!-----------------------------------------------------------------------
!  nonstandard 3D fields
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_TRITIUM_DECAY,'TRITIUM_DECAY',3,        &
                          long_name='tritium decay',                   &
                          units='pmol/m^3/s', grid_loc='3111',               &
                          coordinates='TLONG TLAT z_t time')

!-----------------------------------------------------------------------
!EOC

 end subroutine tr3he_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: tr3he_init_sflux
! !INTERFACE:

 subroutine tr3he_init_sflux

! !USES:

   use forcing_tools, only: find_forcing_times

! !DESCRIPTION:
!  Initialize surface flux computations for tr3he tracer module.
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: sub_name = 'tr3he_mod:tr3he_init_sflux'

   integer (int_kind) :: &
      n,                 & ! index for looping over tracers
      iblock               ! index for looping over blocks

   real (r8), dimension (nx_block,ny_block) :: WORK

   real (r8), dimension (nx_block,ny_block,12,max_blocks_clinic), target :: &
      WORK_READ            ! temporary space to read in fields

!-----------------------------------------------------------------------

   select case (tr3he_formulation)

   case ('file')

!-----------------------------------------------------------------------
!  allocate space for interpolate_forcing
!-----------------------------------------------------------------------

      allocate(INTERP_WORK(nx_block,ny_block,max_blocks_clinic,1))

!-----------------------------------------------------------------------
!  first, read ice file
!-----------------------------------------------------------------------

      allocate(fice_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(fice_file%input%file_fmt, &
                      fice_file%input%filename, &
                      fice_file%input%file_varname, &
                      WORK_READ)
      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         fice_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            fice_file%DATA(:,:,iblock,1,n) = c0
         fice_file%DATA(:,:,iblock,1,n) = &
            fice_file%DATA(:,:,iblock,1,n) * fice_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(fice_file%data_time, &
                              fice_file%data_inc, fice_file%interp_type, &
                              fice_file%data_next, fice_file%data_time_min_loc, &
                              fice_file%data_update, fice_file%data_type)

!-----------------------------------------------------------------------
!  next, read piston velocity file
!-----------------------------------------------------------------------

      allocate(xkw_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(xkw_file%input%file_fmt, &
                      xkw_file%input%filename, &
                      xkw_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         xkw_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            xkw_file%DATA(:,:,iblock,1,n) = c0
         xkw_file%DATA(:,:,iblock,1,n) = &
            xkw_file%DATA(:,:,iblock,1,n) * xkw_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(xkw_file%data_time, &
                              xkw_file%data_inc, xkw_file%interp_type, &
                              xkw_file%data_next, xkw_file%data_time_min_loc, &
                              xkw_file%data_update, xkw_file%data_type)

!-----------------------------------------------------------------------
!  last, read atmospheric pressure file
!-----------------------------------------------------------------------

      allocate(ap_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(ap_file%input%file_fmt, &
                      ap_file%input%filename, &
                      ap_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         ap_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            ap_file%DATA(:,:,iblock,1,n) = c0
         ap_file%DATA(:,:,iblock,1,n) = &
            ap_file%DATA(:,:,iblock,1,n) * ap_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(ap_file%data_time, &
                              ap_file%data_inc, ap_file%interp_type, &
                              ap_file%data_next, ap_file%data_time_min_loc, &
                              ap_file%data_update, ap_file%data_type)

   case ('model')

      if (my_task == master_task) then
         write(stdout,*)  &
            'Using fields from model forcing for calculating Tr 3He fluxes'
      endif

   case default
      call document(sub_name, 'tr3he_formulation', tr3he_formulation)

      call exit_POP(sigAbort, &
                    'tr3he_init_sflux: Unknown value for tr3he_formulation')

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine tr3he_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: tr3he_set_sflux
! !INTERFACE:

 subroutine tr3he_set_sflux(U10_SQR,IFRAC,PRESS,SST,SSS,RHO, &
                          SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  Compute tritium and helium 3 surface fluxes and store related tavg 
!  fields for subsequent accumulating.

! !REVISION HISTORY:
!  same as module

! !USES:

   use time_management, only: thour00
   use forcing_tools, only: update_forcing_data, interpolate_forcing
   use timers, only: timer_start, timer_stop

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      U10_SQR,   & ! 10m wind speed squared (cm/s)**2
      IFRAC,     & ! sea ice fraction (non-dimensional)
      PRESS,     & ! sea level atmospheric pressure (dyne/cm**2)
      SST,       & ! sea surface temperature (C)
      SSS,       & ! sea surface salinity (psu)
      RHO          ! sea surface density (kg/m^3)

   real (r8), dimension(nx_block,ny_block,tr3he_tracer_cnt,max_blocks_clinic), &
         intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,tr3he_tracer_cnt,max_blocks_clinic), &
         intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock             ! block index

   integer (int_kind) :: i, j

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      IFRAC_USED,      & ! used ice fraction (non-dimensional)
      XKW_USED,        & ! part of piston velocity (cm/s)
      AP_USED            ! used atm pressure (converted from dyne/cm**2 to Pa)

   real (r8), dimension(nx_block,ny_block) :: &
      SURF_VALS,       & ! filtered surface tracer values
      He3_SCHMIDT,     & ! helium 3 Schmidt number
      He4_SOL_0,       & ! solubility of helium 4 (mol/m^3/Pa)
      He3_SURF_SAT,    & ! helium 3 surface saturation (mol/m^3)
      He_ALPHA_SOL,    & ! temperature-dependent solubility fractionation
      XKW_ICE,         & ! common portion of piston vel., (1-fice)*xkw (m/s)
      PV,              & ! piston velocity (m/s)
      FLUX               ! tracer flux (pmol/m^2/s)

   character (char_len) :: &
      tracer_data_label          ! label for what is being updated

   character (char_len), dimension(1) :: &
      tracer_data_names          ! short names for input data fields

   integer (int_kind), dimension(1) :: &
      tracer_bndy_loc,         & ! location and field type for ghost
      tracer_bndy_type           !    cell updates

   logical (log_kind), save :: &
      first = .true.

!-----------------------------------------------------------------------
!  local parameters
!-----------------------------------------------------------------------

   real (r8), parameter ::         &
      !xkw_coeff  = 8.6e-7_r8,      & ! xkw_coeff = 0.31 cm/hr s^2/m^2 in (s/m)
      xkw_coeff  = 6.02e-7_r8,     & ! 0.7 * 8.6e-7 from Stanley et. al 2015  (s/m)
      Xhe        = 5.24e-6_r8,     & ! atmospheric helium 4 mole fraction (mol/mol)
      Ir         = 1.384e-6_r8,    & ! 3He/4He isotopic ratio
      alpha_diff = 1.0496_r8,      & ! 3He/4He diffusivity ratio (Bourg & Sposito 2008)
      MH2O       = 18.01528e-3_r8, & ! molecular weight of water in kg
      g          = 9.806_r8,       & ! gravitational acceleration in m/s^2
      mol2pmol   = 1.e+12_r8,      & ! mol -> pmol conversion factor
      pmol2mol   = 1.e-12_r8         ! pmol -> mol conversion factor

!-----------------------------------------------------------------------

   call timer_start(tr3he_sflux_timer)

   do iblock = 1, nblocks_clinic
      IFRAC_USED(:,:,iblock) = c0
      XKW_USED(:,:,iblock) = c0
      AP_USED(:,:,iblock) = c0
   end do

!-----------------------------------------------------------------------
!  Interpolate gas flux forcing data if necessary
!-----------------------------------------------------------------------

   if (tr3he_formulation == 'file') then
       if (thour00 >= fice_file%data_update) then
          tracer_data_names = fice_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Ice Fraction'
          call update_forcing_data(          fice_file%data_time,   &
               fice_file%data_time_min_loc,  fice_file%interp_type, &
               fice_file%data_next,          fice_file%data_update, &
               fice_file%data_type,          fice_file%data_inc,    &
               fice_file%DATA(:,:,:,:,1:12), fice_file%data_renorm, &
               tracer_data_label,            tracer_data_names,     &
               tracer_bndy_loc,              tracer_bndy_type,      &
               fice_file%filename,           fice_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            fice_file%DATA(:,:,:,:,1:12), &
            fice_file%data_time,         fice_file%interp_type, &
            fice_file%data_time_min_loc, fice_file%interp_freq, &
            fice_file%interp_inc,        fice_file%interp_next, &
            fice_file%interp_last,       0)
       IFRAC_USED = INTERP_WORK(:,:,:,1)

       if (thour00 >= xkw_file%data_update) then
          tracer_data_names = xkw_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Piston Velocity'
          call update_forcing_data(         xkw_file%data_time,   &
               xkw_file%data_time_min_loc,  xkw_file%interp_type, &
               xkw_file%data_next,          xkw_file%data_update, &
               xkw_file%data_type,          xkw_file%data_inc,    &
               xkw_file%DATA(:,:,:,:,1:12), xkw_file%data_renorm, &
               tracer_data_label,           tracer_data_names,    &
               tracer_bndy_loc,             tracer_bndy_type,     &
               xkw_file%filename,           xkw_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            xkw_file%DATA(:,:,:,:,1:12), &
            xkw_file%data_time,         xkw_file%interp_type, &
            xkw_file%data_time_min_loc, xkw_file%interp_freq, &
            xkw_file%interp_inc,        xkw_file%interp_next, &
            xkw_file%interp_last,       0)
       XKW_USED = INTERP_WORK(:,:,:,1)

       if (thour00 >= ap_file%data_update) then
          tracer_data_names = ap_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Atmospheric Pressure'
          call update_forcing_data(        ap_file%data_time,   &
               ap_file%data_time_min_loc,  ap_file%interp_type, &
               ap_file%data_next,          ap_file%data_update, &
               ap_file%data_type,          ap_file%data_inc,    &
               ap_file%DATA(:,:,:,:,1:12), ap_file%data_renorm, &
               tracer_data_label,          tracer_data_names,   &
               tracer_bndy_loc,            tracer_bndy_type,    &
               ap_file%filename,           ap_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            ap_file%DATA(:,:,:,:,1:12), &
            ap_file%data_time,         ap_file%interp_type, &
            ap_file%data_time_min_loc, ap_file%interp_freq, &
            ap_file%interp_inc,        ap_file%interp_next, &
            ap_file%interp_last,       0)
       AP_USED = INTERP_WORK(:,:,:,1)
   endif

   do iblock = 1, nblocks_clinic

      if (tr3he_formulation == 'file') then
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < 0.2000_r8) &
            IFRAC_USED(:,:,iblock) = 0.2000_r8
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > 0.9999_r8) &
            IFRAC_USED(:,:,iblock) = 0.9999_r8
      endif

      if (tr3he_formulation == 'model') then
         where (LAND_MASK(:,:,iblock))
            IFRAC_USED(:,:,iblock) = IFRAC(:,:,iblock)
            XKW_USED(:,:,iblock) = xkw_coeff * U10_SQR(:,:,iblock) * 1.e-4_r8 !(cm/s)**2 -> (m/s)**2
            AP_USED(:,:,iblock) = PRESS(:,:,iblock)
         endwhere
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < c0) &
            IFRAC_USED(:,:,iblock) = c0
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > c1) &
            IFRAC_USED(:,:,iblock) = c1
      endif

!-----------------------------------------------------------------------
!  convert pressure from dyne/cm**2 to Pascals
!  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
!  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
!-----------------------------------------------------------------------

      !AP_USED(:,:,iblock) = AP_USED(:,:,iblock) * (c1 / 1013.25e+3_r8)
      AP_USED(:,:,iblock) = AP_USED(:,:,iblock) / c10

!-----------------------------------------------------------------------
!  Compute difusive gas exchange for helium 3 (pmol/m^3/s)
!  Stanley et al. Biogeosciences 12, 5199-5210 (2015) 
!-----------------------------------------------------------------------

      He4_SOL_0 = he4_henry_sol_0(LAND_MASK(:,:,iblock), SST(:,:,iblock), &
          SSS(:,:,iblock))

      He_ALPHA_SOL = alpha_sol_he(LAND_MASK(:,:,iblock), SST(:,:,iblock))

      He3_SCHMIDT = schmidt_he4(LAND_MASK(:,:,iblock), SST(:,:,iblock)) &
          / alpha_diff

      where (LAND_MASK(:,:,iblock))
         tr3he_SFLUX_TAVG(:,:,bufind_3He_IFRAC,iblock)     = IFRAC_USED(:,:,iblock)
         tr3he_SFLUX_TAVG(:,:,bufind_3He_XKW,iblock)       = XKW_USED(:,:,iblock)
         tr3he_SFLUX_TAVG(:,:,bufind_3He_ATM_PRESS,iblock) = AP_USED(:,:,iblock)
         tr3he_SFLUX_TAVG(:,:,bufind_3He_SCHMIDT,iblock)   = He3_SCHMIDT

         XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_USED(:,:,iblock)

         He3_SURF_SAT = He4_SOL_0 * Xhe * AP_USED(:,:,iblock) * Ir &
             * He_ALPHA_SOL
         tr3he_SFLUX_TAVG(:,:,bufind_3He_SURF_SAT,iblock) = He3_SURF_SAT

         PV = XKW_ICE * sqrt(660.0_r8 / He3_SCHMIDT)
         tr3he_SFLUX_TAVG(:,:,bufind_3He_PV,iblock) = PV

         SURF_VALS = p5*(SURF_VALS_OLD(:,:,helium3_ind,iblock) + &
                         SURF_VALS_CUR(:,:,helium3_ind,iblock))
         FLUX = PV * (He3_SURF_SAT - SURF_VALS)
         STF_MODULE(:,:,helium3_ind,iblock) = FLUX
         tr3he_SFLUX_TAVG(:,:,buf_ind_3He_GAS_FLUX,iblock) = FLUX
      elsewhere
         STF_MODULE(:,:,helium3_ind,iblock) = c0
      endwhere

   end do
   !$OMP END PARALLEL DO

   call timer_stop(tr3he_sflux_timer)

!-----------------------------------------------------------------------
!EOC

 end subroutine tr3he_set_sflux

!***********************************************************************
!BOP
! !IROUTINE: tr3he_set_interior
! !INTERFACE:

 subroutine tr3he_set_interior(k, TRACER_MODULE_OLD, TRACER_MODULE_CUR, &
         DTRACER_MODULE, this_block)

! !DESCRIPTION:
!  Compute time derivatives for tritium and helium 3
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use grid, only: KMT

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: &
      k                    ! vertical level index

   real (r8), dimension(nx_block,ny_block,km,tr3he_tracer_cnt), intent(in) :: &
      TRACER_MODULE_OLD, & ! old tracer values
      TRACER_MODULE_CUR    ! current tracer values

   type (block), intent(in) :: &
      this_block          ! block info for the current block

! !OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,tr3he_tracer_cnt), intent(out) :: &
      DTRACER_MODULE       ! computed source/sink term

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: &
      subname = 'tr3he_mod:tr3he_set_interior'

   real (r8), parameter ::          &
      spd    = 86400.0_r8,          & ! number of seconds in a day
      spy    = 365.0_r8*spd,        & ! number of seconds in a year
      hlife  = 12.31_r8,            & ! tritium half-life in years
      lambda = -log(p5)/(hlife*spy)   ! tritium decay constant (1/sec)

   real (r8), dimension(nx_block,ny_block) :: &
      tritium_loc,  & ! local copy of model tritium
      helium3_loc,  & ! local copy of model helium 3
      TRITIUM_DECAY   ! tritium decay term

   integer (int_kind) :: &
      bid             ! local_block id

!-----------------------------------------------------------------------

   bid = this_block%local_id

   DTRACER_MODULE = c0

!-----------------------------------------------------------------------
!  create local copies of model tracers
!  treat negative values as zero
!-----------------------------------------------------------------------

   tritium_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,tritium_ind) + &
                              TRACER_MODULE_CUR(:,:,k,tritium_ind)))
!   helium3_loc      = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,helium3_ind) + &
!                              TRACER_MODULE_CUR(:,:,k,helium3_ind)))

!-----------------------------------------------------------------------
!  apply LAND_MASK to local copies 
!-----------------------------------------------------------------------

   where (.not. LAND_MASK(:,:,bid) .or. k > KMT(:,:,bid))
       tritium_loc = c0
       !helium3_loc = c0
   end where

!-----------------------------------------------------------------------
!  compute time derivatives
!-----------------------------------------------------------------------

    TRITIUM_DECAY = lambda * tritium_loc
    DTRACER_MODULE(:,:,tritium_ind) = -TRITIUM_DECAY
    DTRACER_MODULE(:,:,helium3_ind) = TRITIUM_DECAY

    call accumulate_tavg_field(TRITIUM_DECAY, tavg_TRITIUM_DECAY,bid,k)

!-----------------------------------------------------------------------
!EOC

 end subroutine tr3he_set_interior

!***********************************************************************
!BOP
! !IROUTINE: he4_bunsen_sol_0
! !INTERFACE:

 function he4_bunsen_sol_0(LAND_MASK, SST, SSS)

! !DESCRIPTION:
!  Compute the Bunsen solubility coefficient for helium 4 
!  Ref : Weiss (1971), Journal of Chemical & Engineering Data
!        Vol 16, No. 2 (Table II and Equation 1)
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK          ! land mask for this block

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SST,             & ! sea surface temperature (C)
      SSS                ! sea surface salinity (psu)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: &
      he4_bunsen_sol_0           ! Bunsen solubility of helium 4 

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter ::     &
       A1  = -34.6261_r8,      &
       A2  =  43.0285_r8,      &
       A3  =  14.1391_r8,      &
       B1  = -0.042340_r8,     &
       B2  =  0.022624_r8,     &
       B3  = -0.003312_r8

   real (r8), dimension(nx_block,ny_block) :: &
       Tk     ! sea surface temperature in Kelvin

!-----------------------------------------------------------------------
!  compute Bunsen solubility as in Weiss 1971
!-----------------------------------------------------------------------

   where (LAND_MASK)
       Tk = SST + T0_Kelvin
       he4_bunsen_sol_0 = exp(A1 + A2 * 100.0_r8/Tk + A3 * log(Tk/100.0_r8) &
           + SSS * (B1 + B2 * Tk/100.0_r8 + B3 * (Tk/100.0_r8)**2))
   elsewhere
       he4_bunsen_sol_0 = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end function he4_bunsen_sol_0

!***********************************************************************
!BOP
! !IROUTINE: he4_henry_sol_0
! !INTERFACE:

 function he4_henry_sol_0(LAND_MASK, SST, SSS)

! !DESCRIPTION:
!  Compute the Henry's Law solubility coefficient for helium 4 in pmol/(m^3 Pa)
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK          ! land mask for this block

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SST,             & ! sea surface temperature (C)
      SSS                ! sea surface salinity (psu)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: &
      he4_henry_sol_0    ! Henry's Law solubility coefficient of helium 4 

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block) :: &
       Bsol     ! Bunsen solubility

!-----------------------------------------------------------------------
!  compute Bunsen solubility as in Weiss 1971
!-----------------------------------------------------------------------

    Bsol = he4_bunsen_sol_0(LAND_MASK, SST, SSS)

!-----------------------------------------------------------------------
!  convert Bunsen solubility to Henry's law coefficient (pmol/(m^3 atm))
!-----------------------------------------------------------------------

   where (LAND_MASK)
       he4_henry_sol_0 = Bsol / (R * T0_Kelvin) * 1.e+12_r8 ! mol -> pmol
   elsewhere
       he4_henry_sol_0 = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end function he4_henry_sol_0

!***********************************************************************
!BOP
! !IROUTINE: alpha_sol_he
! !INTERFACE:

 function alpha_sol_he(LAND_MASK, SST)

! !DESCRIPTION:
!  Temperature-dependent solubility fractionation factor for He
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK          ! land mask for this block

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SST                ! sea surface temperature (C)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: &
      alpha_sol_he       ! solubility fractionation factor

!-----------------------------------------------------------------------

   where (LAND_MASK)
       alpha_sol_he = 0.98144_r8 + 8.666e-5_r8 * SST
   elsewhere
       alpha_sol_he = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end function alpha_sol_he

!***********************************************************************
!BOP
! !IROUTINE: schmidt_he4
! !INTERFACE:

 function schmidt_he4(LAND_MASK, SST)

! !DESCRIPTION:
!  Compute Schmidt number for helium 4 (Wanninkhof 1992) 
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK          ! land mask for this block

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SST                ! sea surface temperature (C)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: &
      schmidt_he4       ! Schmidt number for helium 4

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
       A = 410.14_r8,      &
       B = 20.503_r8,      &
       C = 0.53175_r8,     &
       D = 0.0060111_r8

!-----------------------------------------------------------------------

   where (LAND_MASK)
       schmidt_he4 = A - B*SST + C*SST**2 - D*SST**3
   elsewhere
       schmidt_he4 = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end function schmidt_he4

!***********************************************************************
!BOP
! !IROUTINE: diff_he4
! !INTERFACE:

 function diff_he4(LAND_MASK, SST)

! !DESCRIPTION:
!  Compute diffusivity coefficient for helium 4 in m^2/s 
!  Jahne et al. JGR 92 C10 (1987)
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK          ! land mask for this block

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SST                ! sea surface temperature (C)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: &
      diff_he4           ! diffusivity coefficient

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
       A = 818.e-9_r8,     &  ! m^2/s
       E = 11.70e+3_r8        ! kJ/mol -> J/mol = (Pa*m^3)/mol

   real (r8), dimension(nx_block,ny_block) :: &
       Tk                     ! sea surface temperature in Kelvin

!-----------------------------------------------------------------------

   where (LAND_MASK)
       Tk = SST + T0_Kelvin
       diff_he4 = A * exp(-E/(R*Tk))
   elsewhere
       diff_he4 = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end function diff_he4

!***********************************************************************
!BOP
! !IROUTINE: tr3he_tavg_forcing
! !INTERFACE:

 subroutine tr3he_tavg_forcing

! !DESCRIPTION:
!  Make accumulation calls for forcing related tavg fields. This is
!  necessary because the forcing routines are called before tavg flags
!  are set.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! block loop index

!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1, nblocks_clinic
         call accumulate_tavg_field(tr3he_SFLUX_TAVG(:,:,bufind_3He_IFRAC,iblock),     &
             tavg_3He_IFRAC,iblock,1)
         call accumulate_tavg_field(tr3he_SFLUX_TAVG(:,:,bufind_3He_XKW,iblock),       &
             tavg_3He_XKW,iblock,1)
         call accumulate_tavg_field(tr3he_SFLUX_TAVG(:,:,bufind_3He_ATM_PRESS,iblock), &
             tavg_3He_ATM_PRESS,iblock,1)
         call accumulate_tavg_field(tr3he_SFLUX_TAVG(:,:,bufind_3He_SCHMIDT,iblock),   &
             tavg_3He_SCHMIDT,iblock,1)
         call accumulate_tavg_field(tr3he_SFLUX_TAVG(:,:,bufind_3He_PV,iblock),        &
             tavg_3He_PV,iblock,1)
         call accumulate_tavg_field(tr3he_SFLUX_TAVG(:,:,bufind_3He_SURF_SAT,iblock),  &
             tavg_3He_SURF_SAT,iblock,1)
         call accumulate_tavg_field(tr3he_SFLUX_TAVG(:,:,buf_ind_3He_GAS_FLUX,iblock),  &
             tavg_3He_GAS_FLUX,iblock,1)
   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine tr3he_tavg_forcing

!***********************************************************************

end module tr3he_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
