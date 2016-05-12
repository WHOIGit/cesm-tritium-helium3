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
   !use thermodyn_mod

   implicit none
   save

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: &
       tr3he_tracer_cnt, &
       tr3he_init, &
       tr3he_set_interior
!       tr3he_set_sflux,  &
!       tr3he_tavg_forcing

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracers
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
       tr3he_tracer_cnt = 2

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      tritium_ind =  1,  & ! tritium
      helium3_ind =  2     ! helium 3

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
      tr3he_formulation,     & ! how to calculate flux (ocmip or model)
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
      tavg_HELIUM3_IFRAC,      & ! tavg id for ice fraction
      tavg_HELIUM3_XKW,        & ! tavg id for xkw
      tavg_HELIUM3_ATM_PRESS,  & ! tavg id for atmospheric pressure
      tavg_HELIUM3_SCHMIDT,    & ! tavg id for 3He Schmidt number
      tavg_HELIUM3_PV,         & ! tavg id for 3He piston velocity
      tavg_HELIUM3_surf_sat      ! tavg id for 3He surface saturation

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

   use constants, only: char_blank, delim_fmt
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
      !tracer_d_module(n)%long_name  = ind_name_table(n)%name
      tracer_d_module(n)%units      = 'TU'
      tracer_d_module(n)%tend_units = 'TU'
      tracer_d_module(n)%flux_units = 'mmol/cm^2/s'
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

   tritium_file    = 'unknown'
   model_year      = 1
   data_year       = 1931
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

   case ('one', 'ccsm_startup_spunup')
      TRACER_MODULE = c1
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d tritium and 3He set to all ones'
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

   call define_tavg_field(tavg_HELIUM3_IFRAC,'HELIUM3_IFRAC',2,           &
                          long_name='Ice Fraction for 3He fluxes',&
                          units='fraction', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_HELIUM3_XKW,'HELIUM3_XKW',2,               &
                          long_name='XKW for 3He fluxes',         &
                          units='cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_HELIUM3_ATM_PRESS,'HELIUM3_ATM_PRESS',2,           &
                          long_name='Atmospheric Pressure for 3He fluxes',&
                          units='atmospheres', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_HELIUM3_SCHMIDT,'HELIUM3_SCHMIDT',2,     &
                          long_name='3He Schmidt Number',       &
                          units='none', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_HELIUM3_PV,'HELIUM3_PV',2,               &
                          long_name='3He piston velocity',      &
                          units='cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_HELIUM3_surf_sat,'HELIUM3_surf_sat',2,   &
                          long_name='3He Saturation',           &
                          units='mmol/m^3', grid_loc='2110')
   var_cnt = var_cnt+1

!-----------------------------------------------------------------------

   allocate(tr3he_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   tr3he_SFLUX_TAVG = c0

!-----------------------------------------------------------------------
!  nonstandard 3D fields
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_TRITIUM_DECAY,'TRITIUM_DECAY',3,        &
                          long_name='tritium decay',                   &
                          units='TU/s', grid_loc='3111',               &
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
! !IROUTINE: tr3he_set_interior
! !INTERFACE:

 subroutine tr3he_set_interior(k, TRACER_MODULE_OLD, TRACER_MODULE_CUR, &
         DTRACER_MODULE, this_block)

! !DESCRIPTION:
!  Compute time derivatives for tritium and helium 3
!
! !REVISION HISTORY:
!  same as module

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
!  apply mask to local copies ???
!  TODO: define LAND_MASK as in ecosys module and define block_id stuff
!  NOTE: LAND_MASK is TRUE at ocean points
!-----------------------------------------------------------------------

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

end module tr3he_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
