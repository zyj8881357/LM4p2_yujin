module tiling_input_types_mod

! public type -- used to temporarily store all the input data for a cell
type :: tile_parameters_type
 !miscellanous
 integer :: ntile,nc_grpid,nband
 real,allocatable,dimension(:) :: frac
 !soil
 real,allocatable,dimension(:) :: dat_w_sat
 real,allocatable,dimension(:) :: dat_awc_lm2
 real,allocatable,dimension(:) :: dat_k_sat_ref
 real,allocatable,dimension(:) :: dat_psi_sat_ref
 real,allocatable,dimension(:) :: dat_chb
 real,allocatable,dimension(:) :: dat_heat_capacity_dry
 real,allocatable,dimension(:) :: dat_thermal_cond_dry
 real,allocatable,dimension(:) :: dat_thermal_cond_sat
 real,allocatable,dimension(:) :: dat_thermal_cond_exp
 real,allocatable,dimension(:) :: dat_thermal_cond_scale
 real,allocatable,dimension(:) :: dat_thermal_cond_weight
 real,allocatable,dimension(:,:) :: dat_refl_dry_dir 
 real,allocatable,dimension(:,:) :: dat_refl_dry_dif
 real,allocatable,dimension(:,:) :: dat_refl_sat_dir 
 real,allocatable,dimension(:,:) :: dat_refl_sat_dif 
 real,allocatable,dimension(:) :: dat_emis_dry
 real,allocatable,dimension(:) :: dat_emis_sat
 real,allocatable,dimension(:) :: dat_z0_momentum
 real,allocatable,dimension(:) :: dat_tf_depr
 real,allocatable,dimension(:) :: rsa_exp_global
 real,allocatable,dimension(:) :: gw_res_time
 real,allocatable,dimension(:) :: gw_hillslope_length
 real,allocatable,dimension(:) :: gw_scale_length
 real,allocatable,dimension(:) :: gw_hillslope_zeta_bar
 real,allocatable,dimension(:) :: gw_hillslope_relief
 real,allocatable,dimension(:) :: gw_scale_relief
 real,allocatable,dimension(:) :: gw_soil_e_depth
 real,allocatable,dimension(:) :: gw_scale_soil_depth
 real,allocatable,dimension(:) :: gw_hillslope_a
 real,allocatable,dimension(:) :: gw_hillslope_n
 real,allocatable,dimension(:) :: gw_perm
 real,allocatable,dimension(:) :: gw_scale_perm
 real,allocatable,dimension(:) :: microtopo
 real,allocatable,dimension(:) :: tile_hlsp_length
 real,allocatable,dimension(:) :: tile_hlsp_slope
 real,allocatable,dimension(:) :: tile_hlsp_elev
 real,allocatable,dimension(:) :: tile_hlsp_hpos
 real,allocatable,dimension(:) :: tile_hlsp_width
 !hillslope tiling
 integer,allocatable,dimension(:) :: hidx_k
 integer,allocatable,dimension(:) :: hidx_j
 !vegetation
 integer,allocatable,dimension(:) :: vegn

end type

end module
