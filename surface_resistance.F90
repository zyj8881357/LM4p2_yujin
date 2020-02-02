module surface_resistance_mod
! calculations of soil resistance to water vapor and heat
#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use constants_mod, only : dens_h2o, pi, grav
use fms_mod, only: error_mesg, NOTE, WARNING, FATAL, file_exist, &
     close_file, check_nml_error, stdlog, string, lowercase
use mpp_mod, only: mpp_pe, mpp_root_pe
use land_constants_mod, only : diffusivity_h2o, kin_visc_air, heat_cond_air
use land_debug_mod, only : is_watch_point
use land_data_mod, only : lnd, log_version
use land_numerics_mod, only : gamma
use land_tile_diag_mod, only : register_tiled_diag_field, &
     send_tile_data, diag_buff_type, set_default_diag_filter
use vegn_tile_mod, only: vegn_tile_type, vegn_tile_bwood, vegn_tile_LAI, vegn_tile_SAI
use soil_tile_mod, only: soil_tile_type, soil_data_hydraulic_properties, dz, &
        get_soil_litter_C

implicit none
private

! ==== public interfaces =====================================================
public :: surface_resistance_init
public :: surface_resistance_end

public :: surface_resistances
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'surface_resistance_mod'
character(len=*), parameter :: diag_mod_name = 'sfcres'
#include "../shared/version_variable.inc"

integer, parameter :: &
   RESIST_NONE   = 0, & ! no resistance
   RESIST_HO2013 = 1    ! soil resistance based on Haghighi and Or (2013) and related papers

! ==== module variables ======================================================
logical :: module_is_initialized = .FALSE.

character(32) :: soil_resistance_to_use = 'none' ! or 'HO2013'
  ! resistances in soil upper layer and viscous sublayer
real :: rav_lit_0         = 0.0 ! constant litter resistance to vapor
real :: rav_lit_vi        = 0.0 ! litter resistance to vapor per LAI+SAI
real :: rav_lit_fsc       = 0.0 ! litter resistance to vapor per fsc
real :: rav_lit_ssc       = 0.0 ! litter resistance to vapor per ssc
real :: rav_lit_deadmic   = 0.0 ! litter resistance to vapor per dead microbe C
real :: rav_lit_bwood     = 0.0 ! litter resistance to vapor per bwood

namelist /surface_resistance_nml/ &
    soil_resistance_to_use, &
    rav_lit_0, rav_lit_vi, rav_lit_fsc, rav_lit_ssc, rav_lit_deadmic, rav_lit_bwood

integer :: soil_resistance_option = -1

! ---- diag field IDs
integer :: id_r_litt_evap, id_r_bl_sens, id_r_bl_evap, id_r_sv_evap, id_d_visc, &
           id_ustar_sfc, id_u_sfc

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ---------------------------------------------------------------------------------------
subroutine surface_resistance_init(id_ug)
  integer, intent(in) :: id_ug   !<Unstructured axis id.

  integer :: unit, io, iostat, ierr
  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=surface_resistance_nml, iostat=io)
    ierr = check_nml_error(io, 'surface_resistance_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=surface_resistance_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'surface_resistance_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  unit=stdlog()

  if (mpp_pe() == mpp_root_pe()) then
     write(unit, nml=surface_resistance_nml)
  endif

  ! convert symbolic names of surface resistance options into numeric IDs to
  ! speed up selection during run-time
  if (trim(lowercase(soil_resistance_to_use))=='none') then
     soil_resistance_option = RESIST_NONE
  else if (trim(lowercase(soil_resistance_to_use))=='ho2013') then
     soil_resistance_option = RESIST_HO2013
  else
     call error_mesg('surface_resistance_init',&
          'soil resistance option soil_resistance_to_use="'//&
          trim(soil_resistance_to_use)//'" is invalid, use "none" or "HO"',&
          FATAL)
  endif

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')
  id_r_litt_evap = register_tiled_diag_field( diag_mod_name, 'r_litt_evap', &
       (/id_ug/), lnd%time, 'resistance of litter layer for water vapor', 's/m', missing_value=-1.0 )
  id_r_bl_sens = register_tiled_diag_field( diag_mod_name, 'r_bl_sens', &
       (/id_ug/), lnd%time, 'resistance of viscous boundary layer for heat', 's/m', missing_value=-1.0 )
  id_r_bl_evap = register_tiled_diag_field( diag_mod_name, 'r_bl_evap', &
       (/id_ug/), lnd%time, 'resistance of viscous boundary layer for water vapor', 's/m', missing_value=-1.0 )
  id_r_sv_evap = register_tiled_diag_field( diag_mod_name, 'r_sv_evap', &
       (/id_ug/), lnd%time, 'resistance of soil surface layer for water vapor', 's/m', missing_value=-1.0 )
  id_d_visc = register_tiled_diag_field( diag_mod_name, 'd_visc', &
       (/id_ug/), lnd%time, 'thickness of viscous sublayer', 'm', missing_value=-1.0 )
  id_u_sfc = register_tiled_diag_field( diag_mod_name, 'u_sfc', &
       (/id_ug/), lnd%time, 'near-surface wind velocity', 'm/s', missing_value=-1.0 )
  id_ustar_sfc = register_tiled_diag_field( diag_mod_name, 'ustar_sfc', &
       (/id_ug/), lnd%time, 'friction velocity at the surface', 'm/s', missing_value=-1.0 )

  module_is_initialized = .TRUE.
end subroutine surface_resistance_init

! ---------------------------------------------------------------------------------------
subroutine surface_resistance_end()
  ! do nothing
end subroutine surface_resistance_end

! ---------------------------------------------------------------------------------------
subroutine surface_resistances(soil, vegn, diag, T_sfc, u_sfc, ustar_sfc, p, snow_active, &
       r_evap, r_sens)
  type(soil_tile_type), intent(in) :: soil
  type(vegn_tile_type), intent(in) :: vegn
  type(diag_buff_type), intent(inout) :: diag ! diagnostic buffer
  real, intent(in) :: T_sfc     ! surface temperature, K
  real, intent(in) :: u_sfc     ! near-surface wind velocity, m/s
  real, intent(in) :: ustar_sfc ! friction velocity at the soil surface, m/s
  real, intent(in) :: p         ! surface pressure, N/m2
  logical, intent(in) :: snow_active ! if TRUE, ground is covered by snow
  ! output
  real, intent(out) :: r_evap ! surface resistance for evaporation
  real, intent(out) :: r_sens ! surface resistance for sensible heat

  real :: r_litt_evap ! litter resistance, s/m
  real :: r_sv_evap   ! soil surface resistance to evaporation, s/m
  real :: r_bl_evap   ! viscous sublayer resistance to evaporation, s/m
  real :: r_bl_sens   ! viscous sublayer resistance to heat flux, s/m
  real :: d_visc      ! thickness of viscous sublayer, m

  if (snow_active) then
     r_litt_evap = 0
  else
     r_litt_evap = evap_resistance_litter(soil,vegn)
  endif

  select case(soil_resistance_option)
  case(RESIST_NONE)
      r_sv_evap = 0
      r_bl_evap = 0
      r_bl_sens = 0
      d_visc    = 0
  case(RESIST_HO2013)
      r_sv_evap = soil_evap_sv_resistance(soil)
      d_visc    = sfc_visc_bl_depth(u_sfc, ustar_sfc, T_sfc, p)
      r_bl_evap = soil_evap_bl_resistance(soil, T_sfc, p, d_visc)
      r_bl_sens = d_visc/heat_cond_air(T_sfc)
  case default
     call error_mesg(module_name, 'invalid surface resistance option', FATAL)
  end select

  r_evap = r_bl_evap + r_sv_evap + r_litt_evap
  r_sens = r_bl_sens
  if (is_watch_point()) then
     __DEBUG4__(r_bl_evap, r_bl_sens, r_sv_evap, d_visc)
  endif

  ! diagnostic section
  call send_tile_data(id_r_litt_evap, r_litt_evap, diag)
  call send_tile_data(id_r_bl_sens,   r_bl_sens,   diag)
  call send_tile_data(id_r_bl_evap,   r_bl_evap,   diag)
  call send_tile_data(id_r_sv_evap,   r_sv_evap,   diag)
  call send_tile_data(id_d_visc,      d_visc,      diag)
  call send_tile_data(id_u_sfc,       u_sfc,       diag)
  call send_tile_data(id_ustar_sfc,   ustar_sfc,   diag)

end subroutine surface_resistances

! ---------------------------------------------------------------------------------------
! additional resistance of litter to the water vapor flux.
! not a good parameterization, but just using for sensitivity analyses now.
! ignores differing biomass and litter turnover rates.
real function evap_resistance_litter(soil,vegn) result(rav_lit)
  type(soil_tile_type), intent(in) :: soil
  type(vegn_tile_type), intent(in) :: vegn

  real :: litter_fast_C, litter_slow_C, litter_deadmic_C ! litter carbon pools

  call get_soil_litter_C(soil, litter_fast_C, litter_slow_C, litter_deadmic_C)

  associate(cc=>vegn%cohorts)
  rav_lit = rav_lit_0 + rav_lit_vi * (vegn_tile_LAI(vegn)+vegn_tile_SAI(vegn)) &
                      + rav_lit_fsc * litter_fast_C &
                      + rav_lit_ssc * litter_slow_C &
                      + rav_lit_deadmic * litter_deadmic_C &
                      + rav_lit_bwood * vegn_tile_bwood(vegn)
  end associate
end function evap_resistance_litter

! ---------------------------------------------------------------------------------------
! internal soil viscous resistance to evaporation.
! Haghighi et al. (2013): Evaporation rates across a convective air boundary layer
!     are dominated by diffusion. Water Resources Research, 49, No.3, 1602–1610,
!     doi:10.1002/wrcr.20166.
real function soil_evap_sv_resistance(soil) result(r_sv)
   type(soil_tile_type), intent(in) :: soil ! soil properties and parameters

   real, parameter :: gam = 1.73e-5 ! unit conversion constant (Haghighi et al. 2013) eq(13)

   real :: vlc(1), vsc(1) ! volumetric soil and ice content
   real :: psi(1)   ! soil matric potential
   real :: DThDP(1), DKDP(1), DPsi_min, DPsi_max ! unused
   real :: K_z(1)   ! hydraulic conductivity in vertical, kg/(m2 s)
   real :: K_x(1)   ! hydraulic conductivity in horizontal, unused
   ! perhaps we can make soil_data_hydraulic_properties return only requested parameters?

   vlc(1) = max(0.0, soil%wl(1) / (dens_h2o * dz(1)))
   vsc(1) = max(0.0, soil%ws(1) / (dens_h2o * dz(1)))
   call soil_data_hydraulic_properties (soil, vlc, vsc, &
                    psi, DThDP, K_z, K_x, DKDP, DPsi_min, DPsi_max )
   ! calculate resistance to transport within soil upper layer
   r_sv = gam*dens_h2o/(4*K_z(1))
end function soil_evap_sv_resistance

! ---------------------------------------------------------------------------------------
! soil resistance to evaporation in the viscous sublayer
! Haghighi and Or (2013): Evaporation from porous surfaces into turbulent airflows: Coupling
!      eddy characteristics with pore scale vapor diffusion. Water Resources Research, 49,
!      8432–8442, doi:10.1002/2012WR013324View.
real function soil_evap_bl_resistance(soil, T_sfc, p, d_bl) result(r_bl)
   type(soil_tile_type), intent(in) :: soil
   real, intent(in) :: T_sfc     ! surface temperature, K
   real, intent(in) :: p         ! pressure, N/m2
   real, intent(in) :: d_bl      ! thickness of viscous sublayer, m

   real, parameter :: sfc_tension_h2o = 0.071 ! surface tension of liquid water, J/m2

   real :: f_diff   ! surface wetness dependency factor, unitless
   real :: psi_sat_sfc ! saturated matric water potential at the surface, m
   real :: theta_sfc ! relative soil wetness at the surface, unitless
   real :: r_pores  ! surface pore radius, m
   real :: diff_h2o ! diffusivity of water vapor

   ! relative wetness of the surface:
   theta_sfc = max(0.0, soil%wl(1) / (dens_h2o * dz(1)))/soil%pars%vwc_sat
   if (theta_sfc > 0) then
      ! surface wetness dependency model, Schlunder (1988):
      f_diff = (sqrt(pi/(4*theta_sfc))-1)/(pi*sqrt(theta_sfc))
      ! to set reasonable value for theta_sfc > pi/4: f_diff = 0 means that water diffuses
      ! like from a completely wet surface
      f_diff = max(f_diff,0.0)
      ! pore radius (should really be moved into initialization or soil properties update):
      psi_sat_sfc = abs(soil%pars%psi_sat_ref/soil%alpha(1)) ! saturated matric water potential at the surface, m
      r_pores = 2*sfc_tension_h2o/(dens_h2o*grav*psi_sat_sfc)
      ! diffusivity of water vapor for current conditions:
      diff_h2o = diffusivity_h2o(T_sfc,p)
      ! finally, calculate resistance
      r_bl = (d_bl + r_pores*sqrt(pi)*f_diff)/diff_h2o
      if (is_watch_point()) then
         __DEBUG5__(theta_sfc,f_diff,r_pores,diff_h2o,r_bl)
      endif
   else ! theta_sfc <= 0
      r_bl = HUGE(r_bl)
      if (is_watch_point()) then
         __DEBUG2__(theta_sfc,r_bl)
      endif
   endif
end function soil_evap_bl_resistance

! ---------------------------------------------------------------------------------------
! average depth of viscous sublayer
! Haghighi and Or (2013): Evaporation from porous surfaces into turbulent airflows: Coupling
!      eddy characteristics with pore scale vapor diffusion. Water Resources Research, 49,
!      8432–8442, doi:10.1002/2012WR013324View.
real function sfc_visc_bl_depth(u_sfc, ustar_sfc, T, p) result(d_bl)
   real, intent(in) :: u_sfc     ! near-surface wind velocity, m/s
   real, intent(in) :: ustar_sfc ! friction velocity at the soil surface, m/s
   real, intent(in) :: T         ! temperature, degK
   real, intent(in) :: p         ! pressure, N/m2

   real, parameter :: c2 = 2.2   ! parameter of viscous sublayer depth (Haghighi & Or 2013),
                                 ! eq (11), (Haghighi & Or 2015) eq (10a)
   real, parameter :: c3 = 112.0 ! parameter of average eddy exposure time (Haghighi & Or,
                                 ! 2013) eq (15)

   real :: alpha    ! parameter of eddy exposure time probability distribution, unitless
   real :: visc

   ! parameter of eddy exposure time probability distribution, (Haghighi & Or 2013) eq (A4)
   ! for ustar_sfc > 0.3 u_sfc we set alpha to zero, effectively turning gamma-distribution
   ! to exponential distribution.
   alpha = max(0.3 * u_sfc/ustar_sfc-1.0,0.0)
   ! kinematic viscosity of air
   visc = kin_visc_air(T,p)
   d_bl = visc/ustar_sfc * c2*sqrt(c3)/sqrt(alpha+1) * gamma(alpha+1.5)/gamma(alpha+1)
   if (is_watch_point()) then
   __DEBUG5__(d_bl, alpha, visc, u_sfc, ustar_sfc)
   endif
end function sfc_visc_bl_depth

end module surface_resistance_mod