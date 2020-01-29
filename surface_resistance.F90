module surface_resistance_mod
! calculations of soil resistance to water vapor and heat

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use constants_mod, only : dens_h2o, pi, grav
use fms_mod, only: error_mesg, NOTE, WARNING, FATAL, file_exist, &
     close_file, check_nml_error, stdlog, string, lowercase
use mpp_mod, only: mpp_pe, mpp_root_pe
use land_constants_mod, only : diffusivity_h2o, kin_visc_air
use land_data_mod, only : log_version
use land_numerics_mod, only : gamma
use vegn_tile_mod, only: vegn_tile_type, vegn_tile_bwood, vegn_tile_LAI, vegn_tile_SAI
use soil_tile_mod, only: soil_tile_type, soil_data_hydraulic_properties, dz, &
        get_soil_litter_C

implicit none
private

! ==== public interfaces =====================================================
public :: surface_resistance_init
public :: surface_resistance_end

public :: evap_resistance_litter
public :: evap_resistance_soil

public :: soil_evap_sv_resistance
public :: soil_evap_bl_resistance
public :: sfc_visc_bl_depth
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'surface_resistance_mod'
#include "../shared/version_variable.inc"

! ==== module variables ======================================================
logical :: module_is_initialized = .FALSE.

logical :: do_evap_resistance_soil = .FALSE. ! if TRUE, evaporation is slowed down by
  ! resistances in soil upper layer and viscous sublayer
real :: rav_lit_0         = 0.0 ! constant litter resistance to vapor
real :: rav_lit_vi        = 0.0 ! litter resistance to vapor per LAI+SAI
real :: rav_lit_fsc       = 0.0 ! litter resistance to vapor per fsc
real :: rav_lit_ssc       = 0.0 ! litter resistance to vapor per ssc
real :: rav_lit_deadmic   = 0.0 ! litter resistance to vapor per dead microbe C
real :: rav_lit_bwood     = 0.0 ! litter resistance to vapor per bwood

namelist /surface_resistance_nml/ &
    do_evap_resistance_soil, &
    rav_lit_0, rav_lit_vi, rav_lit_fsc, rav_lit_ssc, rav_lit_deadmic, rav_lit_bwood

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ---------------------------------------------------------------------------------------
subroutine surface_resistance_init()
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

  module_is_initialized = .TRUE.
end subroutine surface_resistance_init

! ---------------------------------------------------------------------------------------
subroutine surface_resistance_end()
  ! do nothing
end subroutine surface_resistance_end

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
real function evap_resistance_soil(soil,u_sfc,ustar_sfc,p) result(res)
  type(soil_tile_type), intent(in) :: soil
  real, intent(in) :: u_sfc     ! near-surface wind velocity, m/s
  real, intent(in) :: ustar_sfc ! friction velocity at the soil surface, m/s
  real, intent(in) :: p         ! surface pressure, N/m2

  if (do_evap_resistance_soil) then
     res = soil_evap_sv_resistance(soil) &
         + soil_evap_bl_resistance(soil, u_sfc, ustar_sfc, p)
  else
     res = 0.0
  endif
end function evap_resistance_soil

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
real function soil_evap_bl_resistance(soil, u_sfc, ustar_sfc, p) result(r_bl)
   type(soil_tile_type), intent(in) :: soil
   real, intent(in) :: u_sfc     ! near-surface wind velocity, m/s
   real, intent(in) :: ustar_sfc ! friction velocity at the soil surface, m/s
   real, intent(in) :: p         ! pressure, N/m2

   real, parameter :: sfc_tension_h2o = 0.071 ! surface tension of liquid water, J/m2

   real :: f_diff   ! surface wetness dependency factor, unitless
   real :: alpha    ! parameter of eddy exposure time probability distribution, unitless
   real :: d_bl     ! depth of viscous sublayer, m
   real :: psi_sat_sfc ! saturated matric water potential at the surface, m
   real :: theta_sfc ! relative soil wetness at the surface, unitless
   real :: r_pores  ! surface pore radius, m
   real :: T_sfc    ! temperature, degK

   theta_sfc = max(0.0, soil%wl(1) / (dens_h2o * dz(1)))/soil%pars%vwc_sat
   T_sfc = soil%T(1)
   f_diff = (sqrt(pi/theta_sfc)-1)/(pi*sqrt(theta_sfc)) ! surface wetness dependency model, Schlunder (1988)
   d_bl = sfc_visc_bl_depth(u_sfc, ustar_sfc, T_sfc, p) ! depth of viscous sublayer
   psi_sat_sfc = abs(soil%pars%psi_sat_ref/soil%alpha(1)) ! saturated matric water potential at the surface, m
   r_pores = 2*sfc_tension_h2o/(dens_h2o*grav*psi_sat_sfc)
   r_bl = (d_bl + r_pores*f_diff)/diffusivity_h2o(T_sfc,p)
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

   alpha = 0.3 * u_sfc/ustar_sfc ! parameter of eddy exposure time probablity distribution,
                                 ! (Haghighi & Or 2013) eq (A4)
   d_bl = c2*sqrt(c3)*kin_visc_air(T,p)/(sqrt(alpha+1)*ustar_sfc)*gamma(alpha+1.5)/gamma(alpha+0.5)
end function sfc_visc_bl_depth

end module surface_resistance_mod