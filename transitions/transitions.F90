module land_transitions_mod
#include <fms_platform.h>

#include "../shared/debug.inc"

use constants_mod, only : PI

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif
use mpp_io_mod, only : mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY

use fms_mod, only : string, error_mesg, FATAL, WARNING, NOTE, &
     mpp_pe, lowercase, file_exist, close_file, &
     check_nml_error, stdlog, mpp_root_pe, fms_error_handler

use time_manager_mod, only : time_type, set_date, get_date, set_time, &
     operator(+), operator(-), operator(>), operator(<), operator(<=), operator(/), &
     operator(//), operator(==), days_in_year, print_date, increment_date, get_time, &
     valid_calendar_types, get_calendar_type
use get_cal_time_mod, only : get_cal_time
use horiz_interp_mod, only : horiz_interp_type, horiz_interp_init, &
     horiz_interp_new, horiz_interp_del
use time_interp_mod, only : time_interp
use diag_manager_mod, only : register_diag_field, send_data, diag_field_add_attribute

use nfu_mod, only : nfu_validtype, nfu_inq_var, nfu_get_dim_bounds, nfu_get_rec, &
     nfu_get_dim, nfu_get_var, nfu_get_valid_range, nfu_is_valid

use vegn_data_mod, only : &
     N_LU_TYPES, LU_PAST, LU_CROP, LU_NTRL, LU_SCND, LU_RANGE, LU_URBN, &
     landuse_name, landuse_longname, &
     M_LU_TYPES, LU_IRRIG

use cana_tile_mod, only : cana_tile_heat
use snow_tile_mod, only : snow_tile_heat
use vegn_tile_mod, only : vegn_tile_heat, vegn_tile_type, vegn_tile_bwood
use soil_tile_mod, only : soil_tile_heat

use land_tile_mod, only : land_tile_map, &
     land_tile_type, land_tile_list_type, land_tile_enum_type, new_land_tile, delete_land_tile, &
     first_elmt, tail_elmt, loop_over_tiles, operator(==), current_tile, &
     land_tile_list_init, land_tile_list_end, nitems, elmt_at_index, &
     erase, remove, insert, merge_land_tile_into_list, &
     get_tile_water, land_tile_carbon, land_tile_heat, &
     empty
use land_tile_io_mod, only : print_netcdf_error
use land_tile_diag_mod, only : cmor_name

use land_data_mod, only : lnd, log_version, horiz_interp_ug
use vegn_harvesting_mod, only : vegn_cut_forest

use land_debug_mod, only : set_current_point, is_watch_cell, &
     get_current_point, check_var_range, log_date
use land_numerics_mod, only : rank_descending

implicit none
private

! ==== public interface =====================================================
public :: land_transitions_init
public :: land_transitions_end
public :: save_land_transitions_restart

public :: land_transitions
public :: land_irrigatedareas_init
! ==== end of public interface ==============================================

! ==== module constants =====================================================
character(len=*), parameter :: module_name = 'land_transitions_mod'
character(len=*), parameter :: diag_mod_name = 'landuse'
#include "../shared/version_variable.inc"

! selectors for overshoot handling options, for efficiency
integer, parameter :: &
     OPT_IGNORE = 0, &
     OPT_STOP   = 1, &
     OPT_REPORT = 2
integer, parameter :: &
     DISTR_LM3 = 0, &
     DISTR_MIN = 1
! order of transitions (resulting land use types, hight to low priority) for the
! min-n-tiles transition distribution option. ALL land use types MUST be present in this
! array, otherwise some transitions may be missed -- except perhaps LU_NTRL, since we
! assume there are no transitions to LU_NTRL
integer, parameter :: tran_order(N_LU_TYPES) = (/LU_URBN, LU_CROP, LU_PAST, LU_RANGE, LU_SCND, LU_NTRL/)

! TODO: describe differences between data sets

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

! ==== data types ===========================================================
! set of variables that are summed up on input
type :: var_set_type
   character(64) :: name  = '' ! internal lm3 name of the field
   integer       :: nvars = 0  ! number of variable ids
   integer, allocatable :: id(:)    ! ids of the input fields
end type

! a description of single transition
type :: tran_type
   integer :: donor    = 0  ! kind of donor tile
   integer :: acceptor = 0  ! kind of acceptor tile
   real    :: frac     = 0  ! area of transition
end type tran_type

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.

integer :: tran_ncid  = -1 ! netcdf id of the input file
integer :: state_ncid = -1 ! netcdf id of the input file, if any
integer :: nlon_in, nlat_in

type(var_set_type) :: input_tran  (N_LU_TYPES,N_LU_TYPES) ! input transition rate fields
type(var_set_type) :: input_state (N_LU_TYPES,N_LU_TYPES) ! input state field (for initial transition only)

integer :: diag_ids  (N_LU_TYPES,N_LU_TYPES)
real, allocatable :: norm_in  (:,:) ! normalizing factor to convert input data to
        ! units of [fractions of vegetated area per year]
type(time_type), allocatable :: time_in(:) ! time axis in input data
type(time_type), allocatable :: state_time_in(:) ! time axis in input data
type(horiz_interp_type), save :: interp ! interpolator for the input data
type(time_type) :: time0 ! time of previous transition calculations

integer :: tran_distr_opt = -1 ! selector for transition distribution option, for efficiency
integer :: overshoot_opt = -1 ! selector for overshoot handling options, for efficiency
integer :: conservation_opt = -1 ! selector for non-conservation handling options, for efficiency

! translation of luh2 names and LM3 landuse types
character(5) :: luh2name(12)
integer      :: luh2type(12)
integer :: idata
data (luh2name(idata), luh2type(idata), idata = 1, 12) / &
   'primf', LU_NTRL, &
   'primn', LU_NTRL, &
   'secdf', LU_SCND, &
   'secdn', LU_SCND, &
   'urban', LU_CROP, &
   'c3ann', LU_CROP, &
   'c4ann', LU_CROP, &
   'c3per', LU_CROP, &
   'c4per', LU_CROP, &
   'c3nfx', LU_CROP, &
   'pastr', LU_PAST, &
   'range', LU_RANGE /

! variables for LUMIP diagnostics
integer, parameter :: N_LUMIP_TYPES = 4, &
   LUMIP_PSL = 1, LUMIP_PST = 2, LUMIP_CRP = 3, LUMIP_URB = 4
character(4), parameter :: lumip_name(N_LUMIP_TYPES) = ['psl ','pst ','crop','urbn']
integer :: &
   id_frac_in (N_LUMIP_TYPES) = -1, &
   id_frac_out(N_LUMIP_TYPES) = -1
! translation table: model land use types -> LUMIP types: for each of the model
! LU types it lists the corresponding LUMIP type.
integer, parameter :: lu2lumip(N_LU_TYPES) = [LUMIP_PST, LUMIP_CRP, LUMIP_PSL, LUMIP_PSL, LUMIP_URB, LUMIP_PST]

! variables for irrigation
integer :: nlon_in_manag, nlat_in_manag
type(var_set_type) :: input_manag  (1), input_flood(1) ! input management fields
integer :: diag_area, diag_floodarea
integer :: manag_ncid = -1
type(time_type), allocatable :: manag_time_in(:) ! time axis in input data
type(horiz_interp_type), save :: interp_manag ! interpolator for the input data
type(time_type) :: time0_manag ! time of previous transition calculations
character(len=5), public, parameter  :: &
     crop_name (1) = (/'c3ann'/) !,'c4per', 'c3nfx' /)

! ---- namelist variables ---------------------------------------------------
logical, protected, public :: do_landuse_change = .FALSE. ! if true, then the landuse changes with time
character(len=1024) :: input_file  = '' ! input data set of transition dates
character(len=1024) :: state_file  = '' ! input data set of LU states (for initial transition only)
character(len=1024) :: static_file = '' ! static data file, for input land fraction
character(len=16)  :: data_type  = 'luh1' ! or 'luh2'
! distribute_transitions sets how the land use transitions are distributed among
! tiles within grid cells. 'lm3' is traditional (transitions applied to every
! tile in equal measure, except secondary-to-secondary); 'min-tiles' applies
! transitions to tiles in the order of priority, thereby minimizing the number
! of resulting tiles
logical :: rangeland_is_pasture = .FALSE. ! if true, rangeland is combined with pastures.
! This only applies to luh2 transitions, since there is no rangeland in luh1 anyway.
character(len=16)  :: distribute_transitions  = 'lm3' ! or 'min-n-tiles'
! sets how to handle transition overshoot: that is, the situation when transition
! is larger than available area of the given land use type.
character(len=16) :: overshoot_handling = 'report' ! or 'stop', or 'ignore'
real :: overshoot_tolerance = 1e-4 ! tolerance interval for overshoots
! specifies how to handle non-conservation
character(len=16) :: conservation_handling = 'stop' ! or 'report', or 'ignore'

! for irrigation
character(len=1024) :: manag_file = '' ! management data file, for input land fraction irrigated area
real :: cost(M_LU_TYPES, M_LU_TYPES)=reshape((/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, &
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, &
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, &
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, &
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, &
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, &
                                               2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0 /), (/M_LU_TYPES,M_LU_TYPES/), order =(/ 2, 1 /)) 
logical :: irrigation_on = .TRUE.

namelist/landuse_nml/do_landuse_change, input_file, state_file, static_file, data_type, &
     rangeland_is_pasture, distribute_transitions, &
     overshoot_handling, overshoot_tolerance, &
     conservation_handling, &
     manag_file, irrigation_on


contains ! ###################################################################

! ============================================================================
subroutine land_transitions_init(id_ug, id_cellarea)
  integer, intent(in) :: id_ug !<Unstructured axis id.
  integer, intent(in) :: id_cellarea !<id of cell area diagnostic fields

  ! ---- local vars
  integer        :: unit, ierr, io, ncid1
  integer        :: year,month,day,hour,min,sec
  integer        :: k1,k2,k3, id, n1,n2

  real, allocatable :: lon_in(:,:),lat_in(:,:) ! horizontal grid of input data
  real, allocatable :: buffer_in(:,:) ! buffers for input data reading
  real, allocatable :: mask_in  (:,:) ! valid data mask on the input data grid

  integer :: dimids(NF_MAX_VAR_DIMS), dimlens(NF_MAX_VAR_DIMS)
  type(nfu_validtype) :: v ! valid values range
  character(len=12) :: fieldname

  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce

  if(module_is_initialized) return
  module_is_initialized = .TRUE.
  call log_version(version, module_name, &
  __FILE__)

  call horiz_interp_init

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=landuse_nml, iostat=io)
  ierr = check_nml_error(io, 'landuse_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=landuse_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'landuse_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=landuse_nml)
  endif

  ! read restart file, if any
  if (file_exist('INPUT/landuse.res')) then
     call error_mesg('land_transitions_init','reading restart "INPUT/landuse.res"',&
          NOTE)
     call mpp_open(unit,'INPUT/landuse.res', action=MPP_RDONLY, form=MPP_ASCII)
     read(unit,*) year,month,day,hour,min,sec
     time0 = set_date(year,month,day,hour,min,sec)
     call mpp_close(unit)
  else
     call error_mesg('land_transitions_init','cold-starting land transitions',&
          NOTE)
     time0 = set_date(0001,01,01);
  endif

  ! parse the transition distribution option
  select case(trim(lowercase(distribute_transitions)))
  case ('lm3')
     tran_distr_opt = DISTR_LM3
  case ('min-n-tiles')
     tran_distr_opt = DISTR_MIN
  case default
     call error_mesg('land_transitions_init','distribute_transitions value "'//&
          trim(distribute_transitions)//'" is incorrect, use "lm3" or "min-n-tiles"',&
          FATAL)
  end select

  ! parse the overshoot handling option
  if (trim(overshoot_handling)=='stop') then
     overshoot_opt = OPT_STOP
  else if (trim(overshoot_handling)=='ignore') then
     overshoot_opt = OPT_IGNORE
  else if (trim(overshoot_handling)=='report') then
     overshoot_opt = OPT_REPORT
  else
     call error_mesg('land_transitions_init','overshoot_handling value "'//&
          trim(overshoot_handling)//'" is illegal, use "stop", "report", or "ignore"',&
          FATAL)
  endif

  ! parse the non-conservation handling option
  if (trim(conservation_handling)=='stop') then
     conservation_opt = OPT_STOP
  else if (trim(conservation_handling)=='ignore') then
     conservation_opt = OPT_IGNORE
  else if (trim(conservation_handling)=='report') then
     conservation_opt = OPT_REPORT
  else
     call error_mesg('land_transitions_init','conservation_handling value "'//&
          trim(conservation_handling)//'" is illegal, use "stop", "report", or "ignore"',&
          FATAL)
  endif

  ! initialize diagnostics
  diag_ids(:,:) = 0

  do k1 = 1,size(diag_ids,1)
  do k2 = 1,size(diag_ids,2)
     ! skip unnamed tiles
     if(landuse_name(k1)=='')cycle
     if(landuse_name(k2)=='')cycle
     ! construct a name of input field and register the field
     fieldname = trim(landuse_name(k1))//'2'//trim(landuse_name(k2))
     diag_ids(k1,k2) = register_diag_field(diag_mod_name,fieldname,(/id_ug/), lnd%time, &
          'rate of transition from '//trim(landuse_longname(k1))//' to '//trim(landuse_longname(k2)),&
          units='1/year', missing_value=-1.0)
  enddo
  enddo
  ! register CMIP/LUMIP transition fields
  do k1 = 1,N_LUMIP_TYPES
     id_frac_in(k1) = register_diag_field(cmor_name, &
         'fracInLut_'//trim(lumip_name(k1)), (/id_ug/), lnd%time, &
         'Gross Fraction That Was Transferred into This Tile From Other Land Use Tiles', &
         units='%', standard_name='area_fraction', area = id_cellarea)
     call diag_field_add_attribute(id_frac_in(k1),'ocean_fillvalue',0.0)
     id_frac_out(k1) = register_diag_field(cmor_name, &
         'fracOutLut_'//trim(lumip_name(k1)), (/id_ug/), lnd%time, &
         'Gross Fraction of Land Use Tile That Was Transferred into Other Land Use Tiles', &
         units='%', standard_name='area_fraction', area = id_cellarea)
     call diag_field_add_attribute(id_frac_out(k1),'ocean_fillvalue',0.0)
  enddo

  ! change rangeland to pasture.
  if (rangeland_is_pasture) then
     ! change the type of transitions
     do k1 = 1,size(luh2type)
        if (luh2type(k1)==LU_RANGE) luh2type(k1)=LU_PAST
     enddo
     ! change land use type in existing tiles
     ce = first_elmt(land_tile_map, ls=lnd%ls )
     do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%vegn)) cycle
        if (tile%vegn%landuse == LU_RANGE) tile%vegn%landuse = LU_PAST
     enddo
  endif

  if (.not.do_landuse_change) return ! do nothing more if no land use requested

  ! stop if landuse.res looks inconsistent
  if (time0>lnd%time) then
     call error_mesg('land_transitions_init',&
          'current time is earlier than the time of last land use transition application',&
          FATAL)
  endif

  ! check that we are starting from potential vegetation
  if (time0==set_date(0001,01,01)) then
     ce = first_elmt(land_tile_map, ls=lnd%ls )
     do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%vegn)) cycle
        if (tile%vegn%landuse /= LU_NTRL) then
            call error_mesg('land_transitions_init', &
                'starting land use transitions, but the land use tiles already exist in the initial conditions', &
                FATAL)
        endif
     enddo
  endif

  if (trim(input_file)=='') call error_mesg('land_transitions_init', &
       'do_landuse_change is requested, but landuse transition file is not specified', &
       FATAL)

  ierr=nf_open(input_file,NF_NOWRITE,tran_ncid)
  if(ierr/=NF_NOERR) call error_mesg('land_transitions_init', &
       'do_landuse_change is requested, but landuse transition file "'// &
       trim(input_file)//'" could not be opened because '//nf_strerror(ierr), FATAL)
  call get_time_axis(tran_ncid,time_in)

  ! initialize arrays of input fields
  select case (trim(lowercase(data_type)))
  case('luh1')
     do k1 = 1,size(input_tran,1)
     do k2 = 1,size(input_tran,2)
        ! construct a name of input field and register the field
        fieldname = trim(landuse_name(k1))//'2'//trim(landuse_name(k2))
        if(trim(fieldname)=='2') cycle ! skip unspecified tiles
        input_tran(k1,k2)%name=fieldname
        call add_var_to_varset(input_tran(k1,k2),tran_ncid,input_file,fieldname)
     enddo
     enddo

  case('luh2')
     ! LUH2 data set has more land use types and transitions than LM3,
     ! therefore several transitions need to be aggregated on input to get
     ! the transitions among LM3 land use types
     do n1 = 1,size(luh2type)
     do n2 = 1,size(luh2type)
        k1 = luh2type(n1)
        k2 = luh2type(n2)
        input_tran(k1,k2)%name=trim(landuse_name(k1))//'2'//trim(landuse_name(k2))
        if (k1==k2.and.k1/=LU_SCND) cycle ! skip transitions to the same LM3 LU type, except scnd2scnd
        call add_var_to_varset(input_tran(k1,k2),tran_ncid,input_file,trim(luh2name(n1))//'_to_'//trim(luh2name(n2)))
     enddo
     enddo
     ! add transitions that are not part "state1_to_state2" variable set
     call add_var_to_varset(input_tran(LU_NTRL,LU_SCND),tran_ncid,input_file,'primf_harv')
     call add_var_to_varset(input_tran(LU_NTRL,LU_SCND),tran_ncid,input_file,'primn_harv')
     call add_var_to_varset(input_tran(LU_SCND,LU_SCND),tran_ncid,input_file,'secmf_harv')
     call add_var_to_varset(input_tran(LU_SCND,LU_SCND),tran_ncid,input_file,'secyf_harv')
     call add_var_to_varset(input_tran(LU_SCND,LU_SCND),tran_ncid,input_file,'secnf_harv')

     if (time0==set_date(0001,01,01)) then
        call error_mesg('land_transitions_init','setting up initial land use transitions', NOTE)
        ! initialize state input for initial transition from all-natural state.
        if (trim(state_file)=='') call error_mesg('land_transitions_init',&
            'starting land use transitions, but land use state file is not specified',FATAL)

        ! open state file
        ierr=nf_open(state_file,NF_NOWRITE,state_ncid)
        if(ierr/=NF_NOERR) call error_mesg('land_transitions_init', 'landuse state file "'// &
             trim(state_file)//'" could not be opened because '//nf_strerror(ierr), FATAL)
        call get_time_axis(state_ncid, state_time_in)
        ! initialize state variable array
        do n2 = 1,size(luh2type)
           k2 = luh2type(n2)
           if (k2==LU_NTRL) cycle
           input_state(LU_NTRL,k2)%name='initial '//trim(landuse_name(LU_NTRL))//'2'//trim(landuse_name(k2))
           call add_var_to_varset(input_state(LU_NTRL,k2),state_ncid,state_file,luh2name(n2))
        enddo
     endif
  case default
     call error_mesg('land_transitions_init','unknown data_type "'&
                    //trim(data_type)//'", use "luh1" or "luh2"', FATAL)
  end select
  if (mpp_pe()==mpp_root_pe()) then
     write(*,*)'land_transitions_init: summary of land use transitions'
     do k1 = 1,size(input_tran,1)
     do k2 = 1,size(input_tran,2)
        if(input_tran(k1,k2)%name/='') &
             write(*,'(a)') varset_descr(tran_ncid,input_tran(k1,k2))
     enddo
     enddo
  endif
  ! initialize the input data grid and horizontal interpolator
  ! find any field that is defined in input data
  id = -1
l1:do k1 = 1,size(input_tran,1)
  do k2 = 1,size(input_tran,2)
     if (.not.allocated(input_tran(k1,k2)%id)) cycle
     do k3 = 1,size(input_tran(k1,k2)%id(:))
        if (input_tran(k1,k2)%id(k3)>0) then
           id = input_tran(k1,k2)%id(k3)
           exit l1 ! from all loops
        endif
     enddo
  enddo
  enddo l1

  if (id<=0) call error_mesg('land_transitions_init',&
         'could not find any land transition fields in the input file', FATAL)

  ! we assume that all transition rate fields are specified on the same grid,
  ! in both horizontal and time "directions". Therefore there is a single grid
  ! for all fields, initialized only once.

  __NF_ASRT__(nfu_inq_var(tran_ncid,id,dimids=dimids,dimlens=dimlens))
  nlon_in = dimlens(1); nlat_in=dimlens(2)
  ! allocate temporary variables
  allocate(buffer_in(nlon_in,nlat_in), &
           mask_in(nlon_in,nlat_in),   &
           lon_in(nlon_in+1,1), lat_in(1,nlat_in+1) )
  ! allocate module data
  allocate(norm_in(nlon_in,nlat_in))

  ! get the boundaries of the horizontal axes and initialize horizontal
  ! interpolator
  __NF_ASRT__(nfu_get_dim_bounds(tran_ncid, dimids(1), lon_in(:,1)))
  __NF_ASRT__(nfu_get_dim_bounds(tran_ncid, dimids(2), lat_in(1,:)))

  ! get the first record from variable and obtain the mask of valid data
  ! assume that valid mask does not change with time
  __NF_ASRT__(nfu_get_rec(tran_ncid,id,1,buffer_in))
  ! get the valid range for the variable
  __NF_ASRT__(nfu_get_valid_range(tran_ncid,id,v))
  ! get the mask
  where (nfu_is_valid(buffer_in,v))
     mask_in = 1
  elsewhere
     mask_in = 0
  end where

  ! calculate the normalizing factor to convert input data to units of
  ! [fraction of vegetated area per year]
  select case (trim(lowercase(data_type)))
  case ('luh1')
     ! LUH1 (CMIP5) data were converted on pre-processing
     norm_in = 1.0
  case ('luh2')
     ! read static file and calculate normalizing factor
     ! LUH2 data are in [fraction of cell area per year]
     if (trim(static_file)=='') call error_mesg('land_transitions_init', &
          'using LUH2 data set, but static data file is not specified', FATAL)
     ierr=nf_open(static_file,NF_NOWRITE,ncid1)
     if(ierr/=NF_NOERR) call error_mesg('land_transitions_init', &
          'using LUH2 data set, but static data file "'// &
          trim(static_file)//'" could not be opened because '//nf_strerror(ierr), FATAL)
     __NF_ASRT__(nfu_get_var(ncid1,'landfrac',buffer_in))
     where (buffer_in > 0.0)
        norm_in = 1.0/buffer_in
     elsewhere
        norm_in = 0.0
        mask_in = 0
     end where
     ierr = nf_close(ncid1)
  case default
     call error_mesg('land_transitions_init','unknown data_type "'&
                    //trim(data_type)//'", use "luh1" or "luh2"', FATAL)
  end select

  ! initialize horizontal interpolator
  call horiz_interp_new(interp, lon_in*PI/180,lat_in*PI/180, &
       lnd%sg_lonb, lnd%sg_latb, &
       interp_method='conservative',&
       mask_in=mask_in, is_latlon_in=.TRUE. )

  ! get rid of temporary allocated data
  deallocate(buffer_in, mask_in,lon_in,lat_in)

end subroutine land_transitions_init

!===========================================================================
subroutine land_irrigatedareas_init(id_ug)
  integer, intent(in) :: id_ug ! the IDs of land diagnostic axes
  ! ---- local vars
  integer        :: unit, ierr, io, ncid1
  integer        :: year,month,day,hour,min,sec
  integer        :: k1,k2,k3, id, n1,n2
 ! type(horiz_interp_type) :: interp
  real, allocatable :: lon_in(:,:),lat_in(:,:) ! horizontal grid of input data
  real, allocatable :: buffer_in(:,:) ! buffers for input data reading
  real, allocatable :: mask_in  (:,:) ! valid data mask on the input data grid

  integer :: dimids(NF_MAX_VAR_DIMS), dimlens(NF_MAX_VAR_DIMS)
  type(nfu_validtype) :: v ! valid values range
  character(len=12) :: fieldname

! if(module_is_initialized) return
!  module_is_initialized = .TRUE.
!  call log_version(version, module_name, __FILE__, tagname)

!  call horiz_interp_init

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=landuse_nml, iostat=io)
  ierr = check_nml_error(io, 'landuse_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=landuse_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'landuse_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=landuse_nml)
  endif

!call error_mesg('land_irrarea_init', 'could not find any land transition fields in the input file', FATAL)

     time0_manag = set_date(0001,01,01);

   diag_area = register_diag_field(diag_mod_name,'irr_area',(/id_ug/), lnd%time, &
         'fraction irrigated area',units='', missing_value=-1.0)

   !diag_floodarea = register_diag_field(diag_mod_name,'flood_area',(/id_ug/), lnd%time, &
   !      'flooded fraction of C3 annual crop area',units='', missing_value=-1.0)

  if (.not.do_landuse_change) return ! do nothing more if no land use requested

  ierr=nf_open(manag_file,NF_NOWRITE,manag_ncid)
  call get_time_axis(manag_ncid,manag_time_in)

  
   do n2 = 1,size(crop_name)
        call add_var_to_varset(input_manag(n2),manag_ncid,manag_file,'irrig'//'_'//crop_name(n2))
        input_manag(n2)%name='irrig'//'_'//crop_name(n2)
!        print*, input_manag(n2)%id, input_manag(n2)%name, 'test', crop_name(n2), 'irrig'//'_'//crop_name(n2), manag_ncid, size(input_manag) 
   enddo
  !call add_var_to_varset(input_flood(1),manag_ncid,manag_file,'flood')
  !input_flood(1)%name='flood'        

  ! initialize the input data grid and horizontal interpolator
  ! find any field that is defined in input data
  id = -1
  l1:do k1 = 1,size(input_manag)
      if (.not.allocated(input_manag(k1)%id)) cycle
        do k3 = 1,size(input_manag(k1)%id(:))
          if (input_manag(k1)%id(k3)>0) then
            id = input_manag(k1)%id(k3)
            exit l1 ! from all loops
        endif
     enddo
   enddo l1
! print*, 'id', id
  ! we assume that all transition rate fields are specified on the same grid,
  ! in both horizontal and time "directions". Therefore there is a single grid
  ! for all fields, initialized only once.

  __NF_ASRT__(nfu_inq_var(manag_ncid,id,dimids=dimids,dimlens=dimlens))
  nlon_in_manag = dimlens(1); nlat_in_manag=dimlens(2)
  ! allocate temporary variables
  allocate(buffer_in(nlon_in_manag,nlat_in_manag), &
           mask_in(nlon_in_manag,nlat_in_manag),   &
           lon_in(nlon_in_manag+1,1), lat_in(1,nlat_in_manag+1) )
  ! allocate module data
!  allocate(norm_in(nlon_in,nlat_in))

  ! get the boundaries of the horizontal axes and initialize horizontal
  ! interpolator
  __NF_ASRT__(nfu_get_dim_bounds(manag_ncid, dimids(1), lon_in(:,1)))
  __NF_ASRT__(nfu_get_dim_bounds(manag_ncid, dimids(2), lat_in(1,:)))
  if(lat_in(1,1).gt.90.) lat_in(1,1)=90.
  if(lat_in(1,1).lt.-90.) lat_in(1,1)=-90. 
  if(lat_in(1,nlat_in_manag+1).lt.-90.) lat_in(1,nlat_in_manag+1)=-90.
  if(lat_in(1,nlat_in_manag+1).gt.90.) lat_in(1,nlat_in_manag+1)=90.   
  ! get the first record from variable and obtain the mask of valid data
  ! assume that valid mask does not change with time
  __NF_ASRT__(nfu_get_rec(manag_ncid,id,1,buffer_in))
  ! get the valid range for the variable
  __NF_ASRT__(nfu_get_valid_range(manag_ncid,id,v))
  ! get the mask
  where (nfu_is_valid(buffer_in,v))
     mask_in = 1
  elsewhere
     mask_in = 0
  end where


!  ! initialize horizontal interpolator
  call horiz_interp_new(interp_manag, lon_in*PI/180,lat_in*PI/180, &
       lnd%sg_lonb, lnd%sg_latb, &
       interp_method='conservative',&
       mask_in=mask_in, is_latlon_in=.TRUE. )

  ! get rid of temporary allocated data
  deallocate(buffer_in, mask_in,lon_in,lat_in)
!  deallocate(in_lonb, in_latb, in_cover, x)

end subroutine land_irrigatedareas_init

! ============================================================================
subroutine get_time_axis(ncid, time_in)
  integer, intent(in) :: ncid
  type(time_type), allocatable :: time_in(:)

  integer :: timedim ! id of the record (time) dimension
  integer :: timevar ! id of the time variable
  character(len=NF_MAX_NAME) :: timename  ! name of the time variable
  character(len=256)         :: timeunits ! units ot time in the file
  character(len=24) :: calendar ! model calendar
  real, allocatable :: time(:)  ! real values of time coordinate
  integer :: i, nrec

  ! get the time axis
  __NF_ASRT__(nf_inq_unlimdim(ncid, timedim))
  __NF_ASRT__(nf_inq_dimlen(ncid, timedim, nrec))
  allocate(time(nrec), time_in(nrec))
  __NF_ASRT__(nfu_get_dim(ncid, timedim, time))
  ! get units of time
  __NF_ASRT__(nf_inq_dimname(ncid, timedim, timename))
  __NF_ASRT__(nf_inq_varid(ncid, timename, timevar))
  timeunits = ' '
  __NF_ASRT__(nf_get_att_text(ncid,timevar,'units',timeunits))
  ! get model calendar
  calendar=valid_calendar_types(get_calendar_type())

  ! loop through the time axis and get time_type values in time_in
  if (index(lowercase(timeunits),'calendar_year')>0) then
     do i = 1,size(time)
        time_in(i) = set_date(nint(time(i)),1,1,0,0,0) ! uses model calendar
     end do
  else
     do i = 1,size(time)
        time_in(i) = get_cal_time(time(i),timeunits,calendar)
     end do
  endif
  deallocate(time)
end subroutine get_time_axis

! ============================================================================
subroutine land_transitions_end()

  module_is_initialized=.FALSE.
  if (do_landuse_change) call horiz_interp_del(interp)
  if(allocated(time_in)) deallocate(time_in)

end subroutine land_transitions_end

! ============================================================================
subroutine add_var_to_varset(varset,ncid,filename,varname)
   type(var_set_type), intent(inout) :: varset
   integer     , intent(in) :: ncid     ! id of netcdf file
   character(*), intent(in) :: filename ! name of the file (for reporting problems only)
   character(*), intent(in) :: varname  ! name of the variable

   integer, allocatable :: id(:)
   integer :: varid, ierr

   if (.not.allocated(varset%id)) then
      allocate(varset%id(10))
      varset%id(:) = -1
   endif
   if (varset%nvars >= size(varset%id)) then
      ! make space for new variables
      allocate(id(size(varset%id)+10))
      id(:) = -1
      id(1:varset%nvars) = varset%id(1:varset%nvars)
      call move_alloc(id,varset%id)
   endif

   ierr = nfu_inq_var(ncid, trim(varname), id=varid)
   select case(ierr)
   case (NF_NOERR)
      call error_mesg('land_transitions_init',&
           'adding field "'//trim(varname)//'" from file "'//trim(filename)//'"'//&
           ' to transition "'//trim(varset%name)//'"',&
           NOTE)
      varset%nvars = varset%nvars+1
      varset%id(varset%nvars) = varid
   case (NF_ENOTVAR)
!       call error_mesg('land_transitions_init',&
!            'field "'//trim(varname)//'" not found in file "'//trim(filename)//'"',&
!            NOTE)
   case default
      call error_mesg('land_transitions_init',&
           'error initializing field "'//varname//&
           '" from file "'//trim(filename)//'" : '//nf_strerror(ierr), FATAL)
   end select
end subroutine add_var_to_varset

! ============================================================================
! read, aggregate, and interpolate set of transitions
subroutine get_varset_data(ncid,varset,rec,frac)
   integer, intent(in) :: ncid
   type(var_set_type), intent(in) :: varset
   integer, intent(in) :: rec
   real, intent(out) :: frac(:)

   real :: buff0(nlon_in,nlat_in)
   real :: buff1(nlon_in,nlat_in)
   integer :: i

   frac = 0.0
   buff1 = 0.0
   do i = 1,varset%nvars
     if (varset%id(i)>0) then
        __NF_ASRT__(nfu_get_rec(ncid,varset%id(i),rec,buff0))
        buff1 = buff1 + buff0
     endif
   enddo
   call horiz_interp_ug(interp,buff1*norm_in,frac)
end subroutine get_varset_data

! ============================================================================
! returns a string representing the parts of the transition
function varset_descr(ncid,varset) result(str)
  character(:), allocatable :: str
  integer, intent(in) :: ncid
  type(var_set_type), intent(in) :: varset

  character(NF_MAX_NAME) :: varname
  integer :: i

  str = trim(varset%name)//' = '
  if (varset%nvars == 0) then
     str = str//'0'
  else
     do i = 1, varset%nvars
        __NF_ASRT__(nf_inq_varname(ncid,varset%id(i),varname))
        if (i==1) then
           str = str//trim(varname)
        else
           str = str//' + '//trim(varname)
        endif
     enddo
  endif
end function varset_descr

! ============================================================================
subroutine save_land_transitions_restart(timestamp)
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  integer :: unit,year,month,day,hour,min,sec

  call mpp_open( unit, 'RESTART/'//trim(timestamp)//'landuse.res', nohdrs=.TRUE. )
  if (mpp_pe() == mpp_root_pe()) then
     call get_date(time0, year,month,day,hour,min,sec)
     write(unit,'(6i6,8x,a)') year,month,day,hour,min,sec, &
          'Time of previous landuse transition calculation'
  endif
  call mpp_close(unit)

end subroutine save_land_transitions_restart


! =============================================================================
subroutine land_transitions (time)
  type(time_type), intent(in) :: time

  ! ---- local vars.
  integer :: i,k,k1,k2,k3,i1,i2,l,m, m1,m2, n
  real    :: frac(lnd%ls:lnd%le)
  type(tran_type), pointer :: transitions(:,:)
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  real    :: w
  real    :: diag(lnd%ls:lnd%le)
  logical :: used
  real    :: irr_area(lnd%ls:lnd%le)
  real    :: irr_frac(lnd%ls:lnd%le,size(crop_name))  
  type(tran_type), pointer :: transitions2(:,:)  
  real :: tran1 (M_LU_TYPES, M_LU_TYPES) ! output array of transitions
  real :: tran0 (N_LU_TYPES, N_LU_TYPES) ! output array of transitions
  real :: fi1(lnd%ls:lnd%le)    
  real :: area0 (lnd%ls:lnd%le, M_LU_TYPES) ! fraction of each land use type before transitions
  real :: temp_area (lnd%ls:lnd%le, M_LU_TYPES) ! fraction of each landuse type before transitions
  real :: atot (lnd%ls:lnd%le)
  real :: check_area  
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile  
  integer :: land_flag = 1
  integer :: ntime


  if (.not.do_landuse_change) &
       return ! do nothing if landuse change not requested
  ! NB: in this case file/interp/data are not initialized, so it is
  ! not even possible to use the code below

  call get_date(time,             year0,month0,day0,hour,minute,second)
  call get_date(time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
  if(year0 == year1) &
!!$  if(day0 == day1) &
       return ! do nothing during a year

  if (mpp_pe()==mpp_root_pe()) &
       call log_date('land_transitions: applying land use transitions on ', time)

  ! get transition rates for current time: read map of transitions, and accumulate
  ! as many time steps in array of transitions as necessary. Note that "transitions"
  ! array gets reallocated inside add_to_transitions as necessary, it has only as many
  ! layers as the max number of transitions occurring at a point at the time.
  transitions => NULL()
  do k1 = 1,N_LU_TYPES
  do k2 = 1,N_LU_TYPES
     ! get transition rate for this specific transition
     frac(:) = 0.0
     if (time0==set_date(0001,01,01).and.state_ncid>0) then
        ! read initial transition from state file
        call time_interp(time, state_time_in, w, i1,i2)
        call get_varset_data(state_ncid,input_state(k1,k2),i1,frac)
     else
        if (any(input_tran(k1,k2)%id(:)>0)) then
           call integral_transition(time0,time,input_tran(k1,k2),frac)
        endif
     endif
     call add_to_transitions(frac,time0,time,k1,k2,transitions)
  enddo
  enddo

  ! save the "in" and "out" diagnostics for the transitions
  do k1 = 1, N_LUMIP_TYPES
     if (id_frac_out(k1) > 0) then
        diag(:) = 0.0
        do k2 = 1, size(transitions,2)
        do i = lnd%ls,lnd%le
           if (transitions(i,k2)%donor>0) then
              if (lu2lumip(transitions(i,k2)%donor) == k1) &
                    diag(i) = diag(i) + transitions(i,k2)%frac
           endif
        enddo
        enddo
        used=send_data(id_frac_out(k1), diag*lnd%ug_landfrac*100.0, time)
     endif
  enddo
  do k1 = 1, N_LUMIP_TYPES
     if (id_frac_in(k1) > 0) then
        diag(:) = 0.0
        do k2 = 1, size(transitions,2)
        do i = lnd%ls,lnd%le
           if (transitions(i,k2)%acceptor>0) then
              if (lu2lumip(transitions(i,k2)%acceptor) == k1) &
                    diag(i) = diag(i) + transitions(i,k2)%frac
           endif
        enddo
        enddo
        used=send_data(id_frac_in(k1), diag*lnd%ug_landfrac*100.0, time)
     endif
  enddo


  ! perform the transitions
if (irrigation_on) then 
!-----irr code begin-----
   irr_frac(:,:) = 0.0
   do k2 = 1,size(crop_name)
      ! get fraction irrigation
      irr_area(:) = 0.0
      ! read initial transition from state file
      call fractions_irr_area(time, input_manag(k2), irr_area) !, input_area(k2), crop_area)
      do l = lnd%ls,lnd%le
         irr_frac(l,k2) = irr_area(l)
      enddo   
   enddo

   fi1(:)=0.0
   do l = lnd%ls,lnd%le
     fi1(l) = irr_frac(l,1)
   enddo

   used = send_data(diag_area,fi1,time)  
   
   atot(:) =0.
   area0(:,:) = 0.0
   do l = lnd%ls,lnd%le
      ce = first_elmt(land_tile_map(l))
      do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%vegn)) cycle ! skip non-vegetated tiles
         atot(l) = atot(l)+tile%frac
         do n = 1,M_LU_TYPES
           if(tile%vegn%landuse == n) & 
            area0(l,n) = area0(l,n)+tile%frac
         enddo
      enddo
   enddo

   allocate(transitions2(lnd%ls:lnd%le,(M_LU_TYPES*M_LU_TYPES-M_LU_TYPES-N_LU_TYPES)))
   temp_area = area0
   tran0 = 0.
   do l = lnd%ls,lnd%le !i,j

       tran0 = 0.
       select case (tran_distr_opt)
       case (DISTR_LM3)
         do k1 = 1,size(transitions,2)!1,N_LU_TYPES
           if (transitions(l,k1)%frac .gt. 0..and.atot(l).gt.0.) then
            check_area = min(transitions(l,k1)%frac, temp_area(l,transitions(l,k1)%donor)/atot(l))! area0(i,j,transitions(i,j,k1)%donor))
            tran0(transitions(l,k1)%donor,transitions(l,k1)%acceptor) = check_area
            temp_area(l,transitions(l,k1)%donor) = temp_area(l,transitions(l,k1)%donor) - check_area*atot(l)
           endif
         enddo
       case (DISTR_MIN)
            ! d_kinds and a_kinds are the arrays of initial and final LU types for each of
             ! the transitions. The arrays are of equal size. For each initial and final
               ! LU types src and dst, there is only one element src->dst in these arrays.
              ! We go in order (U,C,P,S) through the final LU types, and apply all transitions
              ! that convert land to this type. Since initial type for each o transitions
              ! are different, there should not be dependence on the order of operations.
            ! An alternative algorithm would be to arrange d_kinds, a_kinds, and areain
            ! the above order (x->U, x->C, x->P, x->S for any x), and go through the arranged array.
         do k = 1,size(tran_order)
          do k1 = 1,size(transitions,2)!1,N_LU_TYPES
           if (transitions(l,k1)%acceptor == tran_order(k)) then
             if (transitions(l,k1)%frac .gt. 0..and.atot(l).gt.0.) then
               check_area = min(transitions(l,k1)%frac, temp_area(l,transitions(l,k1)%donor)/atot(l))!  area0(i,j,transitions(i,j,k1)%donor))
               tran0(transitions(l,k1)%donor,transitions(l,k1)%acceptor) = check_area
               temp_area(l,transitions(l,k1)%donor) = temp_area(l,transitions(l,k1)%donor) - check_area*atot(l)
             endif
           endif
         enddo
        enddo
       end select
!       do k1 = 1,size(transitions,3)!1,N_LU_TYPES
 !         if (transitions(i,j,k1)%frac .gt. 0.) then
!             !write(*,*) 'test', transitions(i,j,k1)%frac
!             check_area = min(transitions(i,j,k1)%frac, temp_area(i,j,transitions(i,j,k1)%donor))!  area0(i,j,transitions(i,j,k1)%donor))
!             tran0(transitions(i,j,k1)%donor,transitions(i,j,k1)%acceptor) = check_area
!             temp_area(i,j,transitions(i,j,k1)%donor) = temp_area(i,j,transitions(i,j,k1)%donor) - check_area*atot(i,j)
!          !   tran0(transitions(i,j,k1)%donor,transitions(i,j,k1)%acceptor) = transitions(i,j,k1)%frac
!         endif
!       enddo
       call set_current_point(l,1)
       call add_irrigation_transitions(area0(l,:),tran0,cost,fi1(l), atot(l), tran1,verbose=.FALSE.)
!       allocate(transitions2(lnd%is:lnd%ie,lnd%js:lnd%je,count(tran1 .gt. 0.)))
       k3=0
       do m1 = 1,M_LU_TYPES
        do m2 = 1,M_LU_TYPES
          if (tran1(m1,m2) .gt. 0.) then
            if (tran1(m1,m2) .lt. 1e-10) tran1(m1,m2) = 0.
            k3=k3+1
            transitions2(l,k3)%donor = m1
            transitions2(l,k3)%acceptor = m2
            transitions2(l,k3)%frac = tran1(m1,m2)
          endif
        enddo
       enddo

   enddo
   ! perform the transitions
   do l = lnd%ls,lnd%le
       if(empty(land_tile_map(l))) cycle ! skip cells where there is no land
       ! set current point for debugging
       call set_current_point(l,1)
       ! transition land area between different tile types
       call land_transitions_0d(land_tile_map(l), &
          transitions2(l,:)%donor, &
          transitions2(l,:)%acceptor,&
          transitions2(l,:)%frac )
   enddo

   ! deallocate array of transitions
   if (associated(transitions)) deallocate(transitions)
   if (associated(transitions2)) deallocate(transitions2)
 else
  do l = lnd%ls,lnd%le
     ! set current point for debugging
     call set_current_point(l,1)
     ! transition land area between different tile types
     call land_transitions_0d(land_tile_map(l), &
          transitions(l,:)%donor, &
          transitions(l,:)%acceptor,&
          transitions(l,:)%frac )
  enddo
 endif 

  ! deallocate array of transitions
  if (associated(transitions)) deallocate(transitions)

  ! store current time for future reference
  time0=time

end subroutine land_transitions


! =============================================================================
! performs tile transitions in a given grid cell
subroutine land_transitions_0d(d_list,d_kinds,a_kinds,area)
  type(land_tile_list_type), intent(inout) :: d_list ! list of tiles
  integer, intent(in) :: d_kinds(:) ! array of donor tile kinds
  integer, intent(in) :: a_kinds(:) ! array of acceptor tile kinds
  real   , intent(in) :: area(:)    ! array of areas changing from donor tiles to acceptor tiles

  ! ---- local vars
  integer :: i, k
  type(land_tile_type), pointer :: ptr
  type(land_tile_list_type) :: a_list
  type(land_tile_enum_type) :: ts, te
  real :: atot ! total fraction of tiles that can be involved in transitions
  real :: htot ! total fraction heat, for debugging only
  ! variable used for conservation check:
  real :: lmass0, fmass0, cmass0, heat0, &
       soil_heat0, vegn_heat0, cana_heat0, snow_heat0 ! pre-transition values
  real :: lmass1, fmass1, cmass1, heat1, &
       soil_heat1, vegn_heat1, cana_heat1, snow_heat1 ! post-transition values
  real :: lm, fm ! buffers for transition calculations

  ! conservation check code, part 1: calculate the pre-transition grid
  ! cell totals
  lmass0 = 0 ; fmass0 = 0 ; cmass0 = 0 ; heat0 = 0
  soil_heat0 = 0 ;  vegn_heat0 = 0 ; cana_heat0 = 0 ; snow_heat0 = 0
  ts = first_elmt(d_list)
  do while (loop_over_tiles(ts, ptr))
     call get_tile_water(ptr,lm,fm)
     lmass0 = lmass0 + lm*ptr%frac ; fmass0 = fmass0 + fm*ptr%frac

     heat0  = heat0  + land_tile_heat  (ptr)*ptr%frac
     cmass0 = cmass0 + land_tile_carbon(ptr)*ptr%frac

     if(associated(ptr%soil)) soil_heat0 = soil_heat0 + soil_tile_heat(ptr%soil)*ptr%frac
     if(associated(ptr%vegn)) vegn_heat0 = vegn_heat0 + vegn_tile_heat(ptr%vegn)*ptr%frac
     if(associated(ptr%cana)) cana_heat0 = cana_heat0 + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow)) snow_heat0 = snow_heat0 + snow_tile_heat(ptr%snow)*ptr%frac
  enddo

  ! calculate the area that can participate in land transitions
  atot = 0 ; ts = first_elmt(d_list)
  do while (loop_over_tiles(ts,ptr))
     if (associated(ptr%vegn)) atot = atot + ptr%frac
  enddo

  if (is_watch_cell()) then
     write(*,*)'### land_transitions_0d: input parameters ###'
     do i = 1, size(d_kinds)
        __DEBUG4__(i,d_kinds(i),a_kinds(i),area(i))
     enddo

     write(*,*)'### land_transitions_0d: land fractions before transitions (initial state) ###'
     ts = first_elmt(d_list); htot = 0.0; k = 1
     do while (loop_over_tiles(ts,ptr))
        if (associated(ptr%vegn)) then
            write(*,'(i2.2,2x)', advance='no') k; k = k+1
            call dpri('landuse',ptr%vegn%landuse)
            call dpri('area',ptr%frac)
            call dpri('heat',vegn_tile_heat(ptr%vegn))
            call dpri('heat*frac',vegn_tile_heat(ptr%vegn)*ptr%frac)
            write(*,*)
            htot = htot+vegn_tile_heat(ptr%vegn)*ptr%frac
        endif
     enddo
     call dpri('total area=',atot)
     call dpri('total heat=',htot)
     write(*,*)
  endif

  ! split each donor tile and gather the parts that undergo a
  ! transition into a separate list. Note that the kind of the landuse is
  ! changed during this transition, including forest harvesting if necessary.
  ! This has to occur at some time before the tiles are merged, and it seems
  ! to be the most convenient place as both original and final landuse kind
  ! is known for each part.
  call land_tile_list_init(a_list)
  select case (tran_distr_opt)
  case (DISTR_LM3)
     do i = 1,size(d_kinds)
        call split_changing_tile_parts(d_list,d_kinds(i),a_kinds(i),area(i)*atot,a_list)
        ! the factor atot normalizes the transitions to the total area in the grid cell
        ! available for the land use, that is, the area of land excluding lakes and glaciers
     enddo
  case (DISTR_MIN)
     ! d_kinds and a_kinds are the arrays of initial and final LU types for each of
     ! the transitions. The arrays are of equal size. For each initial and final
     ! LU types src and dst, there is only one element src->dst in these arrays.
     !
     ! We go in order (U,C,P,R,S) through the final LU types, and apply all transitions
     ! that convert land to this type. Since initial type for each of these transitions
     ! are different, there should not be dependence on the order of operations.
     !
     ! An alternative algorithm would be to arrange d_kinds, a_kinds, and area in
     ! the above order (x->U, x->C, x->P, x->S for any x), and go through the
     ! arranged array.
     do k = 1,size(tran_order)
        do i = 1,size(a_kinds)
           if (a_kinds(i)==tran_order(k)) then
              call split_changing_tile_parts_by_priority( &
                         d_list,d_kinds(i),a_kinds(i),area(i)*atot,a_list)
           endif
        enddo
     enddo
  end select
  if (is_watch_cell()) then
     write(*,*)'### land_transitions_0d: land fractions after splitting changing parts ###'
     atot = 0 ; ts = first_elmt(d_list)
     do while (loop_over_tiles(ts,ptr))
        if (.not.associated(ptr%vegn)) cycle
        write(*,'(2(a,g23.16,2x))')'   donor: landuse=',ptr%vegn%landuse,' area=',ptr%frac
        atot = atot + ptr%frac
     enddo
     ts = first_elmt(a_list)
     do while (loop_over_tiles(ts, ptr))
        if (.not.associated(ptr%vegn)) cycle
        write(*,'(2(a,g23.16,2x))')'acceptor: landuse=',ptr%vegn%landuse,' area=',ptr%frac
        atot = atot + ptr%frac
     enddo
     write(*,'(a,g23.16)')'total area=',atot
  endif

  ! move all tiles from the donor list to the acceptor list -- this will ensure
  ! that all the tiles that can be merged at this time will be
  te = tail_elmt(d_list)
  do
     ts=first_elmt(d_list)
     if(ts==te) exit ! reached the end of the list
     ptr=>current_tile(ts)
     if(ptr%frac <= 0.0) then
        call erase(ts) ! if area of the tile is zero, free it
     else
        ! otherwise, move it to a_list
        call remove(ts)
        call insert(ptr,a_list)
     endif
  enddo
  ! d_list is empty at this point

  ! merge all generated tiles into the source (donor) list
  te = tail_elmt(a_list)
  do
     ts=first_elmt(a_list)
     if(ts==te) exit ! break out of loop
     ptr=>current_tile(ts)
     call remove(ts)
     call merge_land_tile_into_list(ptr,d_list)
  enddo
  ! a_list is empty at this point
  call land_tile_list_end(a_list)

  if (is_watch_cell()) then
     write(*,*)'### land_transitions_0d: land fractions final state ###'
     ts = first_elmt(d_list); htot = 0.0; k=0
     do while (loop_over_tiles(ts,ptr))
        if (associated(ptr%vegn)) then
            write(*,'(i2.2,2x)', advance='no') k; k = k+1
            call dpri('landuse',ptr%vegn%landuse)
            call dpri('area',ptr%frac)
            call dpri('heat',vegn_tile_heat(ptr%vegn))
            call dpri('heat*frac',vegn_tile_heat(ptr%vegn)*ptr%frac)
            write(*,*)
            htot = htot+vegn_tile_heat(ptr%vegn)*ptr%frac
        endif
     enddo
     call dpri('total area=',atot)
     call dpri('total heat=',htot)
     write(*,*)
  endif

  ! conservation check part 2: calculate grid cell totals in final state, and
  ! compare them with pre-transition totals
  lmass1 = 0 ; fmass1 = 0 ; cmass1 = 0 ; heat1 = 0
  soil_heat1 = 0 ;  vegn_heat1 = 0 ; cana_heat1 = 0 ; snow_heat1 = 0
  ts = first_elmt(d_list)
  do while (loop_over_tiles(ts,ptr))
     call get_tile_water(ptr,lm,fm)
     lmass1 = lmass1 + lm*ptr%frac ; fmass1 = fmass1 + fm*ptr%frac

     heat1  = heat1  + land_tile_heat  (ptr)*ptr%frac
     cmass1 = cmass1 + land_tile_carbon(ptr)*ptr%frac

     if(associated(ptr%soil)) soil_heat1 = soil_heat1 + soil_tile_heat(ptr%soil)*ptr%frac
     if(associated(ptr%vegn)) vegn_heat1 = vegn_heat1 + vegn_tile_heat(ptr%vegn)*ptr%frac
     if(associated(ptr%cana)) cana_heat1 = cana_heat1 + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow)) snow_heat1 = snow_heat1 + snow_tile_heat(ptr%snow)*ptr%frac
  enddo
  call check_conservation ('liquid water', lmass0, lmass1, 1e-6)
  call check_conservation ('frozen water', fmass0, fmass1, 1e-6)
  call check_conservation ('carbon'      , cmass0, cmass1, 1e-6)
  call check_conservation ('canopy air heat content', cana_heat0 , cana_heat1 , 1e-6)
! heat content of vegetation may not conserve because of the cohort merging issues
!  call check_conservation ('vegetation heat content', vegn_heat0 , vegn_heat1 , 1e-6)
  call check_conservation ('snow heat content',       snow_heat0 , snow_heat1 , 1e-6)
  call check_conservation ('soil heat content',       soil_heat0 , soil_heat1 , 1e-4)
  call check_conservation ('heat content', heat0 , heat1 , 1e-4)

end subroutine land_transitions_0d

! =============================================================================
! check that the requested area of transitions is not larger than available area
! in tiles
subroutine check_area_overshoot(area, d_kind, a_kind, dfrac)
  real,    intent(in) :: area   ! total area of donor tiles
  integer, intent(in) :: d_kind ! LU type of donor tiles
  integer, intent(in) :: a_kind ! LU type of acceptor tiles
  real,    intent(in) :: dfrac  ! fraction of land area that changes LU type

  integer :: severity ! severity of overshoot errors
  integer :: i,j,k,face ! coordinates of current point, for overshoot diagnostics

  ! check for overshoot situation: that is, a case where the transition area is
  ! larger than the available area
  if(overshoot_opt /= OPT_IGNORE.and.dfrac>area+overshoot_tolerance) then
     severity = WARNING
     if (overshoot_opt==OPT_STOP) severity = FATAL
     call get_current_point(i,j,k,face)
     call error_mesg('landuse',&
          'transition at ('//trim(string(i))//','//trim(string(j))//&
          ',face='//trim(string(face))//&
          ') from "'//trim(landuse_name(d_kind))// &
          '" to "'  //trim(landuse_name(a_kind))//&
          '" ('//trim(string(dfrac))//') is larger than area of "'&
          //trim(landuse_name(d_kind))//'" ('//trim(string(area))//')', &
          severity)
  endif
end subroutine check_area_overshoot

! =============================================================================
! splits changing parts of donor tiles into a separate tile list, performing
! land use changes in the process
subroutine split_changing_tile_parts_by_priority(d_list,d_kind,a_kind,dfrac,a_list)
  type(land_tile_list_type), intent(in) :: d_list ! list of donor tiles
  integer, intent(in) :: d_kind ! LU type of donor tiles
  integer, intent(in) :: a_kind ! LU type of acceptor tiles
  real,    intent(in) :: dfrac  ! fraction of land area that changes LU type
  type(land_tile_list_type), intent(inout) :: a_list ! list of acceptors

  ! ---- local vars
  type(land_tile_enum_type) :: ct
  type(land_tile_type), pointer :: tile, temp
  real :: area, darea, tfrac
  real,    allocatable :: priority(:) ! priority of the land use transition fro each tile
  integer, allocatable :: idx(:)      ! array of tile indices in the descending priority order
  integer :: k
  integer :: ntiles ! number of tiles in d_list

  ! calculate total area of the tiles that should be transitioned to another kind
  ct = first_elmt(d_list); area = 0.0
  do while (loop_over_tiles(ct, tile))
     if (.not.associated(tile%vegn)) cycle
     if (tile%vegn%landuse == d_kind) area = area + tile%frac
  enddo

  call check_area_overshoot(area,d_kind,a_kind,dfrac)

  ! calculate transition priorities
  ntiles = nitems(d_list)
  allocate(priority(ntiles), idx(ntiles))
  priority(:) = -HUGE(1.0)
  k = 0; ct = first_elmt(d_list)
  do while (loop_over_tiles(ct,tile))
     k = k+1
     if(.not.associated(tile%vegn))  cycle ! skip non-vegetated tiles
     if(tile%vegn%landuse /= d_kind) cycle ! skip tiles that do not match donor LU type
     priority(k) = landuse_priority(tile, a_kind)
  enddo

  ! sort landuse transition priorities in descending order
  call rank_descending(priority, idx)

  ! transition cannot be more than current total area of specified kind
  tfrac = min(dfrac,area)
  do k = 1, ntiles
     if (tfrac==0) exit ! from loop, no more area to transition
     tile=>elmt_at_index(d_list, idx(k))
     if (.not.associated(tile%vegn)) cycle ! landuse cannot be applied to non-vegetated tiles
     if(tile%vegn%landuse /= d_kind) cycle ! skip tiles that do not match donor LU type
     darea = min(tile%frac, tfrac)
     if (darea>0) then
        ! make a copy of current tile
        temp => new_land_tile(tile)
        temp%frac = darea
        tile%frac = tile%frac-darea
        ! convert land use type of the tile: cut the forest, if necessary
        if(temp%vegn%landuse==LU_NTRL.or.temp%vegn%landuse==LU_SCND.or.temp%vegn%landuse==LU_RANGE) &
                call vegn_cut_forest(temp, a_kind)
        ! change landuse type of the tile
        temp%vegn%landuse = a_kind
        ! reset time elapsed since last disturbance and time elapsed since last land use
        ! event in the new tile
        temp%vegn%age_since_disturbance = 0.0
        temp%vegn%age_since_landuse     = 0.0
        ! add the new tile to the resulting list
        call insert(temp, a_list) ! insert tile into output list
        ! calculate remaining area of transition
        tfrac = tfrac-darea
     endif
  enddo

end subroutine split_changing_tile_parts_by_priority

! ============================================================================
! returns priority of the land use tile: tiles with highest number will be
! consumed first by the land use transition
function landuse_priority(tile, dst) result(P); real P
  type(land_tile_type), intent(in) :: tile
  integer, intent(in) :: dst ! land use types we are transitioning to

  integer :: src ! land use type of the tile

  P = -HUGE(1.0) ! very low priority
  if (.not.associated(tile%vegn)) return
  src = tile%vegn%landuse

  if ((src==LU_SCND.or.src==LU_NTRL).and.dst==LU_SCND) then
     ! for wood harvesting (NTRL->SCND or SCND->SCND), first
     ! consume tiles with highest wood biomass
     P = vegn_tile_bwood(tile%vegn)
  else if (dst==LU_SCND) then
     ! for abandonment, we first consume top-of-the-hill tiles
     ! hidx_j is the index of the hillslope tile; the higher the index the
     ! higher the tile in the hillslope
     P = tile%soil%hidx_j
  else if (src==LU_CROP.and.dst==LU_PAST) then
     ! for CROP->PAST conversion, start from the top of the hill
     P = tile%soil%hidx_j
  else
     ! for everything else, start from the bottom
     P = -tile%soil%hidx_j
  endif
end function landuse_priority

! =============================================================================
! splits changing parts of donor tiles into a separate tile list, performing
! land use changes in the process
subroutine split_changing_tile_parts(d_list,d_kind,a_kind,dfrac,a_list)
  type(land_tile_list_type), intent(in) :: d_list ! list of donor tiles
  integer, intent(in) :: d_kind ! LU type of donor tiles
  integer, intent(in) :: a_kind ! LU type of acceptor tiles
  real,    intent(in) :: dfrac  ! fraction of land area that changes LU type
  type(land_tile_list_type), intent(inout) :: a_list ! list of acceptors

  ! ---- local vars
  type(land_tile_enum_type) :: ct
  type(land_tile_type), pointer :: tile, temp
  real :: area, darea, area0, area1
  real :: x0,x1,x2 ! values of transition intensity
  real, parameter :: eps = 1e-6 ! area calculation precision
  real, parameter :: factor = 1.6 ! multiplier for solution bracketing
  integer :: iter

  ! calculate total area of the tiles that should be transitioned to another kind
  ct = first_elmt(d_list); area = 0.0
  do while (loop_over_tiles(ct, tile))
     if (.not.associated(tile%vegn)) cycle
     if (tile%vegn%landuse == d_kind)  &
          area = area + tile%frac
  enddo

  call check_area_overshoot(area,d_kind,a_kind,dfrac)

  ! if area of the tiles of requested kind is zero we cannot transition
  ! anything, so just return
  if (area==0) return

  ! transition cannot be more than current total area of specified kind
  darea = min(dfrac, area)

  ! solve equation to get transition intensity
  ! (1) bracket transition intensity interval so that requested area is within it
  x0=0.0; area0 = total_transition_area(d_list, d_kind, a_kind, x0)
  x1=1.0; area1 = total_transition_area(d_list, d_kind, a_kind, x1)
  iter = 0
  do
     if ((area0<=darea).and.(area1>=darea)) exit
     if (area0>darea) then
        x0 = x0-(x1-x0)*factor
        area0 = total_transition_area(d_list, d_kind, a_kind, x0)
     else
        x1 = x1+(x1-x0)*factor
        area1 = total_transition_area(d_list, d_kind, a_kind, x1)
     endif
     iter = iter+1
     if (iter>50) then
        call error_mesg('veg_tile_transitions',&
             'cannot braket transition intensity interval after 50 iterations',&
             FATAL)
     endif
  enddo

  ! find solution for transition intensity by binary search
  do iter = 1,50
     x2 = (x0+x1)/2
     area = total_transition_area(d_list, d_kind, a_kind, x2)
     if (abs(x1-x2)<eps) exit
     if (area>darea) then
        x1=x2
     else
        x0=x2
     endif
  enddo

  ! do tile transitions to destination list
  ct = first_elmt(d_list)
  do while (loop_over_tiles(ct, tile))
     if(.not.associated(tile%vegn))  cycle ! skip all non-vegetation tiles
     if(tile%vegn%landuse /= d_kind) cycle ! skip all tiles that doe not match "donor" LU kind
     darea = vegn_tran_priority(tile%vegn, a_kind, x2)
     if(tile%frac*darea > 0) then
        ! make a copy of current tile
        temp => new_land_tile(tile)
        temp%frac = tile%frac*darea
        tile%frac = tile%frac*(1.0-darea)
        ! convert land use type of the tile:
        ! cut the forest, if necessary
        if(temp%vegn%landuse==LU_NTRL.or.temp%vegn%landuse==LU_SCND.or.temp%vegn%landuse==LU_RANGE) &
             call vegn_cut_forest(temp, a_kind)
        ! change landuse type of the tile
        temp%vegn%landuse = a_kind
        ! reset time elapsed since last disturbance and time elapsed since last land use
        ! event in the new tile
        temp%vegn%age_since_disturbance = 0.0
        temp%vegn%age_since_landuse     = 0.0
        ! add the new tile to the resulting list
        call insert(temp, a_list) ! insert tile into output list
     endif
  enddo

end subroutine split_changing_tile_parts


! ============================================================================
! calculates total area (fraction of grid cell area) participating in
! vegetation transition from src_kind to dst_kind for given transition
! intensity tau
function total_transition_area(list,src_kind,dst_kind,tau) result (total_area)
  real :: total_area
  type(land_tile_list_type), intent(in) :: list ! list of tiles
  integer , intent(in) :: src_kind, dst_kind ! source and destination kinds
  real    , intent(in) :: tau                ! transition intensity

  ! ---- local vars
  type(land_tile_enum_type) :: ct
  type(land_tile_type), pointer :: tile

  total_area = 0
  ct = first_elmt(list)
  do while (loop_over_tiles(ct, tile))
     if (.not.associated(tile%vegn)) cycle ! skip non-vegetated tiles
     if(tile%vegn%landuse == src_kind) &
          total_area = total_area + tile%frac*vegn_tran_priority(tile%vegn,dst_kind,tau)
  enddo

end function total_transition_area


! ============================================================================
! given a vegetation patch, destination kind of transition, and "transition
! intensity" value, this function returns a fraction of tile that will parti-
! cipate in transition.
!
! this function must be contiguous, monotonic, its value must be within
! interval [0,1]
!
! this function is used to determine what part of each tile is to be converted
! to another land use kind; the equation is solved to get "transition intensity"
! tau for which total area is equal to requested. Tau is, therefore, a dummy
! parameter, and only relative values of the priority functions for tiles
! participating in transition have any meaning. For most transitions the priority
! function is just equal to tau: therefore there is no preference, and all tiles
! contribute equally to converted area. For secondary vegetation harvesting,
! however, priority also depends on wood biomass, and therefore tiles
! with high wood biomass are harvested first.
function vegn_tran_priority(vegn, dst_kind, tau) result(P); real :: P
  type(vegn_tile_type), intent(in) :: vegn
  integer             , intent(in) :: dst_kind
  real                , intent(in) :: tau

  real :: vegn_bwood

  if (vegn%landuse==LU_SCND.and.dst_kind==LU_SCND) then ! secondary biomass harvesting
     vegn_bwood = vegn_tile_bwood(vegn)
     P = max(min(tau+vegn_bwood,1.0),0.0)
  else
     P = max(min(tau,1.0),0.0)
  endif
end function vegn_tran_priority


! ============================================================================

subroutine add_to_transitions(frac, time0,time1,k1,k2,tran)
  real, intent(in) :: frac(lnd%ls:lnd%le)
  type(time_type), intent(in) :: time0       ! time of previous calculation of
    ! transitions (the integral transitions will be calculated between time0
    ! and time)
  type(time_type), intent(in) :: time1       ! current time
  integer, intent(in) :: k1,k2               ! kinds of tiles
  type(tran_type), pointer :: tran(:,:)    ! transition info

  ! ---- local vars
  integer :: k,sec,days,l
  type(tran_type), pointer :: ptr(:,:) => NULL()
  real    :: part_of_year
  logical :: used

  ! allocate array of transitions, if necessary
  if (.not.associated(tran)) allocate(tran(lnd%ls:lnd%le,1))

  do l = lnd%ls, lnd%le
     if(frac(l) == 0) cycle ! skip points where transition rate is zero
     ! find the first empty transition element for the current indices
     k = 1
     do while ( k <= size(tran,2) )
        if(tran(l,k)%donor == 0) exit
        k = k+1
     enddo

     if (k>size(tran,2)) then
        ! if there is no room, make the array of transitions larger
        allocate(ptr(lnd%ls:lnd%le,size(tran,2)*2))
        ptr(:,1:size(tran,2)) = tran
        deallocate(tran)
        tran => ptr
        nullify(ptr)
     end if

     ! store the transition element
     tran(l,k) = tran_type(k1,k2,frac(l))
  enddo

  ! send transition data to diagnostics
  if(diag_ids(k1,k2)>0) then
     call get_time(time1-time0, sec,days)
     part_of_year = (days+sec/86400.0)/days_in_year(time0)
     used = send_data(diag_ids(k1,k2), frac/part_of_year, time1)
  endif

end subroutine add_to_transitions


! ==============================================================================
! given boundaries of time interval [t1,t2], calculates total transition (time
! integral of transition rates) over the specified interval
subroutine integral_transition(t1, t2, tran, frac, err_msg)
  type(time_type), intent(in)  :: t1,t2 ! time boundaries
  type(var_set_type), intent(in)  :: tran ! id of the field
  real           , intent(out) :: frac(:)
  character(len=*),intent(out), optional :: err_msg

  ! ---- local vars
  integer :: n ! size of time axis
  type(time_type) :: ts,te
  integer         :: i1,i2
  real :: w  ! time interpolation weight
  real :: dt ! current time interval, in years
  real :: sum(size(frac(:)))
  integer :: l
  character(len=256) :: msg

  msg = ''
  ! adjust the integration limits, in case they are out of range
  n = size(time_in)
  ts = t1;
  if (ts<time_in(1)) ts = time_in(1)
  if (ts>time_in(n)) ts = time_in(n)
  te = t2
  if (te<time_in(1)) te = time_in(1)
  if (te>time_in(n)) te = time_in(n)

  call time_interp(ts, time_in, w, i1,i2, err_msg=msg)
  if(msg /= '') then
    if(fms_error_handler('integral_transition','Message from time_interp: '//trim(msg),err_msg)) return
  endif
  call get_varset_data(tran_ncid,tran,i1,frac)

  dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
  sum = -frac*w*dt
  do while(time_in(i2)<=te)
     call get_varset_data(tran_ncid,tran,i1,frac)
     dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
     sum = sum+frac*dt
     i2 = i2+1
     i1 = i2-1
     if(i2>size(time_in)) exit ! from loop
  enddo

  call time_interp(te,time_in,w,i1,i2, err_msg=msg)
  if(msg /= '') then
    if(fms_error_handler('integral_transition','Message from time_interp: '//trim(msg),err_msg)) return
  endif
  call get_varset_data(tran_ncid,tran,i1,frac)
  dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
  frac = sum+frac*w*dt
  ! check the transition rate validity
  do l = 1,size(frac(:))
     call set_current_point(l+lnd%ls-1,1)
     call check_var_range(frac(l),0.0,HUGE(1.0),'integral_transition',tran%name, FATAL)
  enddo
end subroutine integral_transition

!====for irrigation=============================================================
! read, aggregate, and interpolate set of transitions
subroutine get_varset_data_manag(ncid,varset,rec,frac)
   integer, intent(in) :: ncid
   type(var_set_type), intent(in) :: varset
   integer, intent(in) :: rec
   real, intent(out) :: frac(:)

   real :: buff0(nlon_in_manag,nlat_in_manag)
   real :: buff1(nlon_in_manag,nlat_in_manag)
   integer :: i

   frac = 0.0
   buff1 = 0.0
   do i = 1,varset%nvars
     if (varset%id(i)>0) then
        __NF_ASRT__(nfu_get_rec(ncid,varset%id(i),rec,buff0))
        buff1 = buff1 + buff0
     endif
   enddo
   call horiz_interp_ug(interp_manag,buff1,frac)
end subroutine get_varset_data_manag

! ============================================================================
subroutine fractions_irr_area(t2, manag, irr_area,  err_msg)
  type(time_type), intent(in)  :: t2 ! time boundaries
  type(var_set_type), intent(in)  :: manag ! id of the field
  real           , intent(out) :: irr_area(:)
  real :: irr_area1(lnd%ls:lnd%le)
  real :: irr_area2(lnd%ls:lnd%le)
  character(len=*),intent(out), optional :: err_msg
  integer :: ntime

  ! ---- local vars
  integer         :: i1,i2
  real :: w  ! time interpolation weight
  character(len=256) :: msg

  msg = ''
  ntime=size(manag_time_in)
  if(t2>manag_time_in(ntime))then
    call time_interp(manag_time_in(ntime), manag_time_in, w, i1,i2, err_msg=msg)    
  else
    call time_interp(t2, manag_time_in, w, i1,i2, err_msg=msg)
  endif
  if(msg /= '') then
    if(fms_error_handler('fractions_irr_area','Message from time_interp:'//trim(msg),err_msg)) return
  endif
  call get_varset_data_manag(manag_ncid,manag,i1,irr_area1)
  call get_varset_data_manag(manag_ncid,manag,i2,irr_area2)
 
  irr_area = irr_area1*(1-w)+irr_area2*w

end subroutine
!=============================== 
!=================================================================
subroutine add_irrigation_transitions(area0,tran0,cost,fi1,atot,tran1,verbose)
  real, intent(in)  :: area0(M_LU_TYPES)
  real, intent(inout)  :: tran0(N_LU_TYPES,   N_LU_TYPES)   ! initial transition matrix
  real, intent(in)  :: cost (M_LU_TYPES, M_LU_TYPES) ! cost of transitions
  real, intent(in)  :: fi1 ! fraction of irrigated area after transition
  real, intent(in)  :: atot
  real, intent(out) :: tran1(M_LU_TYPES, M_LU_TYPES) ! resulting transition matrix
  logical, intent(in), optional :: verbose

  ! local vars
  integer :: map2 (M_LU_TYPES, M_LU_TYPES) ! mapping transition to var number
  integer, allocatable :: map1i(:), map1j(:) ! mapping var number to transitions
  integer :: n  ! number of variables
  integer :: m1 ! Number of <= inequalities
  integer :: m2 ! Number of >= inequalities
  integer :: m3 ! Number of == equalities
  integer :: m  ! total number of constraints
  integer :: eq ! equation number
  real,    allocatable :: a(:,:) ! input matrix for simplex method
  integer, allocatable :: iposv(:), izrov(:)
  real    :: area00(N_LU_TYPES)
  real    :: area1 (M_LU_TYPES) ! fraction of each land use type after transitions
  integer :: i,j,k,icase
  logical :: verbose_
  real    :: fi0 ! fraction of irrigated area before transition, for diagnostics only
  real    :: tran_temp(N_LU_TYPES, N_LU_TYPES)

  verbose_ = .FALSE.
  if (present(verbose)) verbose_ = verbose

  do i = 1,N_LU_TYPES
    do j = 1,N_LU_TYPES
     if (tran0(i,j) < 1e-16) tran0(i,j)= 0. 
    enddo
  enddo

  tran_temp = tran0*atot
  area00(:) = area0(1:N_LU_TYPES)
  area00(LU_CROP) = area00(LU_CROP)+area0(LU_IRRIG)

!  if (verbose_) then
 if (is_watch_cell()) then
     write(*,*)'INPUT DATA'
     write(*,*)'initial land use fractions'
     do i = 1,M_LU_TYPES
        !write(*,'(a," : ",g)') landuse_name(i),area00(i)
     enddo
     write(*,*)
     write(*,*)'initial transition matrix:'
     do i = 1,N_LU_TYPES
        write(*,'(a5)',advance='NO') landuse_name(i)
        do j = 1,N_LU_TYPES
           write(*,'(x,g10.3)',advance='NO') tran_temp(i,j)
        enddo
        write(*,*)
     enddo
     write(*,'(2x,99(x,a10))') (landuse_name(i),i=1,N_LU_TYPES)
     if ((area00(LU_CROP)) > 0) then
        fi0 = area0(LU_IRRIG)/(area00(LU_CROP))
     else 
        fi0 = 0.0
     endif
     !write(*,'(x,a,99(/2x,a,g:))') 'irrigated cropland fraction', 'before transition:',&
               !fi0,'after transition:',fi1
     write(*,*)
     write(*,*)'relative cost of transitions:'
     do i = 1,M_LU_TYPES
        write(*,'(a5)',advance='NO') landuse_name(i)
        do j = 1,M_LU_TYPES
           write(*,'(x,g10.3)',advance='NO') cost(i,j)
        enddo
        write(*,*)
     enddo
     write(*,'(2x,99(x,a10))') (landuse_name(i),i=1,M_LU_TYPES)
     write(*,*)'END OF INPUT DATA'
  endif

  ! build the translation between indices in output transition matrix and variable number
  map2(:,:) = -1
  n = 0
  do i = 1,M_LU_TYPES
  do j = 1,M_LU_TYPES
     if (i==LU_CROP.or.j==LU_CROP .or.i==LU_IRRIG.or.j==LU_IRRIG) then
        if (i==j) cycle       ! skip c->c and i->i transition
        if (j==LU_NTRL) cycle ! skip x->n transitions
        n = n+1
        map2(i,j) = n
     endif
  enddo
  enddo

  allocate(map1i(n),map1j(n))
  do i = 1,M_LU_TYPES
  do j = 1,M_LU_TYPES
     if (map2(i,j)>0) then
        map1i(map2(i,j)) = i; map1j(map2(i,j)) = j
     endif
  enddo
  enddo

!  if (verbose_) then
 if (is_watch_cell()) then
     write(*,*)
     write(*,*)
     write(*,'(x,a)')'transition index encoding:'
     do i = 1,M_LU_TYPES
        write(*,'(a5)',advance='NO') landuse_name(i)
        do j = 1,M_LU_TYPES
           if (map2(i,j)>0) then
              write(*,'(x,i5)',advance='NO') map2(i,j)
           else
              write(*,'(6x)',advance='NO')
           endif
        enddo
        write(*,*)
     enddo
     write(*,'(8x,99(x,a5))') (landuse_name(i),i=1,M_LU_TYPES)

     do i = 1,size(map1i)
        write(*,'("x",i2.2," : ",a5," -> ",a5)') i, landuse_name(map1i(i)),landuse_name(map1j(i))
     enddo
  endif

  m1 = 0
  m2 = 0
  m3 = 2 + 2*(N_LU_TYPES-1) ! eq 1,11, and N_LU_TYPES of eq 8,9
  ! calculate # of equations (17) and 18
  do i = 1,N_LU_TYPES
     if (i==LU_CROP) cycle
     do j = i+1,N_LU_TYPES
        if (j==LU_CROP) cycle
        if ((cost(i,LU_IRRIG)==cost(j,LU_IRRIG)).and.(cost(i,LU_CROP)==cost(j,LU_CROP))) m3 = m3+1
        if (j==LU_NTRL) cycle
        if ((cost(LU_IRRIG,i)==cost(LU_IRRIG,j)).and.(cost(LU_CROP,i)==cost(LU_CROP,j))) m3 = m3+1
     enddo
  enddo
  
  m = m1 + m2 + m3  ! Total number of constraints

!  if (verbose_) then
 if (is_watch_cell()) then
     write(*,*)
     write(*,'(99(a,I2))'),'Number of variables       (n) :',N
     write(*,'(99(a,I2))'),'Number of constraints     (m) :',M
     write(*,'(99(a,I2))'),'Number of <= inequalities (m1):',M1
     write(*,'(99(a,I2))'),'Number of >= inequalities (m2):',M2
     write(*,'(99(a,I2))'),'Number of == equalities   (m3):',M3
  endif

  ! calculate areas after transition
  area1(1:N_LU_TYPES) = area00(:)
 ! area1(:) = area0(:)
 ! area1(LU_CROP)=area0(LU_CROP)+area0(LU_IRRIG)
  do i = 1,N_LU_TYPES
  do j = 1,N_LU_TYPES
     area1(i) = area1(i) + tran_temp(j,i) - tran_temp(i,j)
  enddo
  enddo

!  if (verbose_) then
  if (is_watch_cell()) then
     write(*,*)
     write(*,*)'total cropland area after transitions:',area1(LU_CROP)
     write(*,*)'irrigated cropland area after transitions:',area1(LU_CROP)*fi1
     write(*,*)'non-irrigated cropland area after transitions:',area1(LU_CROP)*(1-fi1)
  endif
  ! fill up matrix for the simplex method
  allocate (a(m+2,n+1))
  a(:,:) = 0.0

  eq = 1  ! equation (14) in the notes: objective function (i.e. negative cost)
  do k = 1,n
     a(eq,k+1) = -cost(map1i(k),map1j(k))
  enddo
  a(eq,1) = 0.0

  eq = eq+1 ! equation (10) in the notes
  a(eq,1) = area1(LU_CROP)*fi1 - area0(LU_IRRIG)
  do k = 1, N_LU_TYPES
     if (map2(k,LU_IRRIG)>0) a(eq,map2(k,LU_IRRIG)+1) = -1
     if (map2(LU_IRRIG,k)>0) a(eq,map2(LU_IRRIG,k)+1) = +1
  enddo
  if (a(eq,1)<0) a(eq,:) = -a(eq,:)

  eq = eq+1 ! equation (11) in the notes
  a(eq,1) = area1(LU_CROP)*(1-fi1) - area0(LU_CROP)
  do k = 1, N_LU_TYPES
     if (map2(k,LU_CROP)>0) a(eq,map2(k,LU_CROP)+1) = -1
     if (map2(LU_CROP,k)>0) a(eq,map2(LU_CROP,k)+1) = +1
  enddo
  if (a(eq,1)<0) a(eq,:) = -a(eq,:)

  ! equations (8): sum of transitions to irrigated and un-irrigated cropland is equal
  ! to the total transition to cropland
  do k = 1,N_LU_TYPES
     if (k==LU_CROP) cycle
     eq = eq+1
     a(eq,1) = tran_temp(k,LU_CROP)
     a(eq,map2(k,LU_CROP)+1)  = -1.0
     a(eq,map2(k,LU_IRRIG)+1) = -1.0
  enddo
  ! eq (9): sum of transitions from irrigated and un-irrigated cropland is equal
  ! to the total transition to cropland
  do i = 1,N_LU_TYPES
     if (i==LU_CROP) cycle
     eq = eq+1
     a(eq,1) = tran_temp(LU_CROP,i)
     if (map2(LU_CROP, i)>0) a(eq,map2(LU_CROP, i)+1) = -1.0
     if (map2(LU_IRRIG,i)>0) a(eq,map2(LU_IRRIG,i)+1) = -1.0
  enddo
  ! eq (17): proportionality assumption to resolve ambiguities in transitions
  ! to irrigated and non-irrigated cropland
  do i = 1,N_LU_TYPES
     if (i==LU_CROP) cycle
     do j = i+1,N_LU_TYPES
        if (j==LU_CROP) cycle
        if ((cost(i,LU_IRRIG)==cost(j,LU_IRRIG)).and.(cost(i,LU_CROP)==cost(j,LU_CROP))) then
           !write(*,*) i,j,map2(i,LU_IRRIG),map2(j,LU_IRRIG)
           eq = eq+1
           a(eq,1) = 0
           a(eq,map2(i,LU_IRRIG)+1) =  tran_temp(j,LU_CROP)
           a(eq,map2(j,LU_IRRIG)+1) = -tran_temp(i,LU_CROP)
        endif
     enddo
  enddo
  ! eq (18): proportionality assumption to resolve ambiguities in transitions
  ! from irrigated and non-irrigated cropland
  do i = 1,N_LU_TYPES
     if (i==LU_CROP.or.i==LU_NTRL) cycle
     do j = i+1,N_LU_TYPES
        if (j==LU_CROP.or.j==LU_NTRL) cycle
        if ((cost(LU_IRRIG,i)==cost(LU_IRRIG,j)).and.(cost(LU_CROP,i)==cost(LU_CROP,j))) then
          ! write(*,*) i,j,map2(i,LU_IRRIG),map2(j,LU_IRRIG)
           eq = eq+1
           a(eq,1) = 0
           a(eq,map2(LU_IRRIG,i)+1) =  tran_temp(LU_CROP,j)
           a(eq,map2(LU_IRRIG,j)+1) = -tran_temp(LU_CROP,i)
        endif
     enddo
  enddo

!  if (verbose_) then
  if (is_watch_cell()) then
     print *,' Input Table for simplx:'
     write(*,'(8x)',advance='NO')
     do i = 1,size(map1i)
        write(*,'(4x,a1,"->",a1)',advance='NO')landuse_name(map1i(i)),landuse_name(map1j(i))
     enddo
     write(*,*)
     do i = 1, m+1
        write(*,'(99f8.2)') (a(i,j),j=1,n+1)
     enddo
  endif

  allocate(izrov(n),iposv(m))
  call simplx(a,m,n,m1,m2,m3,icase,izrov,iposv)
  if (icase == 0) then
     if (verbose_) write(*,*) 'simplx finished successfully'
  else
     ! TODO: make it a FATAL error
     if (verbose_) write(*,*) 'simplx finished un-successfully, with code',icase
  endif
! write(*,'(x,a,99i3)') 'izrov=',izrov
! write(*,'(x,a,99i3)') 'iposv=',iposv

! print *,' Output Table:'
! write(*,'(8x)',advance='NO')
! do i = 1,size(map1i)
!    write(*,'(4x,a1,"->",a1)',advance='NO')landuse_name(map1i(i)),landuse_name(map1j(i))
! enddo
! write(*,*)
! do i = 1, m+1
!    write(*,'(99f8.2)') (a(i,j),j=1,n+1)
! enddo

!  if (verbose_) then
  if (is_watch_cell()) then
     print *,' '
     print *,' Maximum of objective function = ', A(1,1)
  endif

  do I=1, N
    do J=1, M
      if (IPOSV(J).eq.I) then
        if (is_watch_cell()) then
          write(*,'("  x",i2.2," = ",g10.3,x,a1,"->",a1)') I, A(J+1, 1),landuse_name(map1i(i)),landuse_name(map1j(i))
        endif 
        goto 3
      end if
    end do
    if (is_watch_cell()) then
      write(*,'("  y",i2.2," = ",g10.3,x,a1,"->",a1)') I, 0.0, landuse_name(map1i(i)),landuse_name(map1j(i))
    endif 
3 end do
!  print *,' '

  tran1(:,:) = 0.0
  tran1(1:N_LU_TYPES,1:N_LU_TYPES) = tran_temp(:,:)
  do i=1,n
     do j=1,m
        if (iposv(j)==i) then
           tran1(map1i(i),map1j(i)) = a(j+1, 1)
           goto 4
        end if
     end do
     tran1(map1i(i),map1j(i)) = 0.0
4 end do

  tran1 = tran1/atot

!  if (verbose_) then
  if (is_watch_cell()) then
    write(*,*)
    write(*,*)'final transition matrix:'
    do i = 1,M_LU_TYPES
       write(*,'(a5)',advance='NO') landuse_name(i)
       do j = 1,M_LU_TYPES
          write(*,'(x,g10.3)',advance='NO') tran1(i,j)
       enddo
       write(*,*)
    enddo
    write(*,'(2x,99(x,a10))') (landuse_name(i),i=1,M_LU_TYPES)
  endif

  deallocate (map1i, map1j)
  deallocate (a, izrov, iposv)
end subroutine add_irrigation_transitions

!----------------------------------------------------------------------------------------
! USES simp1,simp2,simp3
!Simplex method for linear programming. Input parameters a, m, n, mp, np, m1, m2, and m3,
!and output parameters a, icase, izrov, and iposv are described above.
!number of variables expected; EPS is the absolute precision, which should be adjusted to
!the scale of your variables.
subroutine simplx(a,m,n,m1,m2,m3,icase,izrov,iposv)
  real,    intent(inout) :: a(:,:)
  integer, intent(in) :: n,m,m1,m2,m3
  integer, intent(out) :: icase
  integer, intent(out) :: izrov(:)
  integer, intent(out) :: iposv(:)

  real, parameter :: EPS = 1e-6

  integer :: m12, nl2, ir
  INTEGER :: i,ip,is,k,kh,kp,nl1,l1(n),l2(m),l3(m)
  REAL bmax,q1
  character(len=512) :: mesg ! for error message


  ! TODO: check size of input tableau a
  if (m.ne.m1+m2+m3) call error_mesg('simplx','Bad input constraint counts in simplx.', FATAL)
  if (size(izrov)<n) call error_mesg('simplx','ERROR: size of izrov is less then number of variables (N)', FATAL)
  if (size(iposv)<m) call error_mesg('simplx','ERROR: size of ipozv is less then number of constraints (M)', FATAL)

  nl1=n
  do k=1,n
    l1(k)=k    !Initialize index list of columns admissible for exchange.
    izrov(k)=k !Initially make all variables right-hand.
  end do
  nl2=m
  do i=1,m
    if(a(i+1,1).lt.0.) then
        call check_var_range(a(i+1,1),0.0,HUGE(1.0),'Bad input tableau in simplx, Constants bi must be nonnegative.','a(i+1,1)', FATAL)
        write(mesg,*) 'Bad input tableau in simplx, Constants bi must be nonnegative.', a(i+1,1)
        call error_mesg('simplx',mesg, FATAL)!pause ' Bad input tableau in simplx, Constants bi must be nonnegative.'
    endif
    l2(i)=i
    iposv(i)=n+i
  !-------------------------------------------------------------------------------------------------
  !Initial left-hand variables. m1 type constraints are represented by having their slackv ariable
  !initially left-hand, with no artificial variable. m2 type constraints have their slack
  !variable initially left-hand, with a minus sign, and their artificial variable handled implicitly
  !during their first exchange. m3 type constraints have their artificial variable initially
  !left-hand.
  !-------------------------------------------------------------------------------------------------
  end do
  do i=1,m2
    l3(i)=1
  end do
  ir=0
  if(m2+m3.eq.0) goto 30 !The origin is a feasible starting solution. Go to phase two.
  ir=1
  do k=1,n+1             !Compute the auxiliary objective function.
    q1=0.
    do i=m1+1,m
      q1=q1+a(i+1,k)
    end do
    a(m+2,k)=-q1
  end do
  10 call simp1(a,m+1,l1,nl1,0,kp,bmax) !Find max. coeff. of auxiliary objective fn
  if(bmax.le.EPS.and.a(m+2,1).lt.-EPS)then
    !write(*,*) 'bmax', bmax, 'a', a(m+2,1)
    icase=-1        !Auxiliary objective function is still negative and canât be improved,
    return          !hence no feasible solution exists.
  else if(bmax.le.EPS.and.a(m+2,1).le.EPS)then
  !Auxiliary objective function is zero and canât be improved; we have a feasible starting vector.
  !Clean out the artificial variables corresponding to any remaining equality constraints by
  !goto 1âs and then move on to phase two by goto 30.
    m12=m1+m2+1
    if (m12.le.m) then
      do ip=m12,m
        if(iposv(ip).eq.ip+n)then !Found an artificial variable for an equalityconstraint.
          call simp1(a,ip,l1,nl1,1,kp,bmax)
          if(bmax.gt.EPS) goto 1  !Exchange with column corresponding to maximum
        end if                    !pivot element in row.
      end do
    end if
    ir=0
    m12=m12-1
    if (m1+1.gt.m12) goto 30
    do i=m1+1,m1+m2               !Change sign of row for any m2 constraints
                                  !still present from the initial basis.
      if(l3(i-m1).eq.1)then
        do k=1,n+1
          a(i+1,k)=-a(i+1,k)
        end do
      end if
    end do
    goto 30                        !Go to phase two.
  end if
  call simp2(a,m,n,l2,nl2,ip,kp,q1) !Locate a pivot element (phase one).
  if(ip.eq.0)then                  !Maximum of auxiliary objective function is
                                   !unbounded, so no feasible solution exists.
    icase=-2
    return
  end if
  1 call simp3(a,m+1,n,ip,kp)
  !Exchange a left- and a right-hand variable (phase one), then update lists.
  if(iposv(ip).ge.n+m1+m2+1)then   !Exchanged out an artificial variable for an
                                   !equality constraint. Make sure it stays
                                   !out by removing it from the l1 list.
    do k=1,nl1
      if(l1(k).eq.kp) goto 2
    end do
  2 nl1=nl1-1
    do is=k,nl1
      l1(is)=l1(is+1)
    end do
  else
    if(iposv(ip).lt.n+m1+1) goto 20
    kh=iposv(ip)-m1-n
    if(l3(kh).eq.0) goto 20      !Exchanged out an m2 type constraint.
    l3(kh)=0                     !If itâs the first time, correct the pivot column
                                 !or the minus sign and the implicit
                                 !artificial variable.
  end if
  a(m+2,kp+1)=a(m+2,kp+1)+1.
  do i=1,m+2
    a(i,kp+1)=-a(i,kp+1)
  end do
  20 is=izrov(kp)                !Update lists of left- and right-hand variables.
  izrov(kp)=iposv(ip)
  iposv(ip)=is
  if (ir.ne.0) goto 10           !if still in phase one, go back to 10.
  !End of phase one code for finding an initial feasible solution. Now, in phase two, optimize it.
  30 call simp1(a,0,l1,nl1,0,kp,bmax) !Test the z-row for doneness.
  if(bmax.le.EPS)then            !Done. Solution found. Return with the good news.
    icase=0
    return
  end if
  call simp2(a,m,n,l2,nl2,ip,kp,q1)  !Locate a pivot element (phase two).
  if(ip.eq.0)then                !Objective function is unbounded. Report and return.
    icase=1
    return
  end if
  call simp3(a,m,n,ip,kp)  !Exchange a left- and a right-hand variable (phase two),
  goto 20                        !update lists of left- and right-hand variables and
                                 !return for another iteration.
end subroutine simplx

!The preceding routine makes use of the following utility subroutines:

! ==============================================================================
! Determines the maximum of those elements whose index is contained in the
! supplied list ll, either with or without taking the absolute value, as flagged
! by iabf.
SUBROUTINE simp1(a,mm,ll,nll,iabf,kp,bmax)
  real,    intent(in)  :: a(:,:) ! tableau
  integer, intent(in)  :: mm,ll(:),nll,iabf
  integer, intent(out) :: kp
  real,    intent(out) :: bmax

  integer :: k
  real    :: test

  kp=ll(1)
  bmax=a(mm+1,kp+1)
  if(nll.lt.2) return
  do k=2,nll
    if(iabf.eq.0)then
      test=a(mm+1,ll(k)+1)-bmax
    else
      test=abs(a(mm+1,ll(k)+1))-abs(bmax)
    endif
    if(test.gt.0.)then
      bmax=a(mm+1,ll(k)+1)
      kp=ll(k)
    endif
  end do
end subroutine simp1

! ==============================================================================
! Locate a pivot element, taking degeneracy into account.
subroutine simp2(a,m,n,l2,nl2,ip,kp,q1)
  real,    intent(in)  :: a(:,:)
  integer, intent(in)  :: m,n,l2(:),nl2, kp
  integer, intent(out) :: ip
  real,    intent(out) :: q1

  real, parameter :: EPS=1.e-6
  integer :: i,k,ii
  real    :: q,q0,qp
  ip=0
  if(nl2.lt.1) return
  do i=1,nl2
    if(a(i+1,kp+1).lt.-EPS) goto 2
  end do
  return  ! No possible pivots. Return with message.
2 q1=-a(l2(i)+1,1)/a(l2(i)+1,kp+1)
  ip=l2(i)
  if(i+1.gt.nl2) return
  do i=i+1, nl2
    ii=l2(i)
    if(a(ii+1,kp+1).lt.-EPS)then
      q=-a(ii+1,1)/a(ii+1,kp+1)
      if(q.lt.q1)then
        ip=ii
        q1=q
        else if (q.eq.q1) then !We have a degeneracy.
        do k=1,n
          qp=-a(ip+1,k+1)/a(ip+1,kp+1)
          q0=-a(ii+1,k+1)/a(ii+1,kp+1)
          if(q0.ne.qp)goto 6
        end do
6       if(q0.lt.qp) ip=ii
      end if
    end if
  end do
end subroutine simp2

! ==============================================================================
! Matrix operations to exchange a left-hand and right-hand variable (see text).
subroutine simp3(a,i1,k1,ip,kp)
  real, intent(inout) :: a(:,:)
  integer, intent(in) :: i1,k1,ip,kp

  integer :: ii,kk
  real    :: piv

  piv=1./a(ip+1,kp+1)
  if (i1.ge.0) then
    do ii=1,i1+1
      if(ii-1.ne.ip)then
        a(ii,kp+1)=a(ii,kp+1)*piv
        do kk=1,k1+1
          if(kk-1.ne.kp)then
            a(ii,kk)=a(ii,kk)-a(ip+1,kk)*a(ii,kp+1)
          end if
        end do
      end if
    end do
  end if
  do kk=1,k1+1
    if(kk-1.ne.kp) a(ip+1,kk)=-a(ip+1,kk)*piv
  end do
  a(ip+1,kp+1)=piv
end subroutine

! ==============================================================================
! checks conservation and aborts with fatal error if tolerance is exceeded
subroutine check_conservation(name, d1, d2, tolerance)
  character(*), intent(in) :: name ! name of the component
  real, intent(in) :: d1,d2 ! values to check
  real, intent(in) :: tolerance ! tolerance of the test

  integer :: curr_i, curr_j, face
  integer :: severity ! severity of the generated message
  character(256) :: message

  if (conservation_opt == OPT_IGNORE) return ! do nothing

  severity = WARNING
  if (conservation_opt==OPT_STOP) severity = FATAL

  if (abs(d1-d2)>tolerance) then
     call get_current_point(i=curr_i,j=curr_j,face=face)
     write(message,'(a,3(x,a,i4), 3(x,a,g23.16))')&
          'conservation of '//trim(name)//' is violated', &
          'at i=',curr_i,'j=',curr_j,'face=',face, &
          'value before=', d1, 'after=', d2, 'diff=',d2-d1
     call error_mesg('land_transitions',message,severity)
  endif
end subroutine check_conservation

end module
