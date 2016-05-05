module predefined_tiles_mod

 use constants_mod     , only : pi
 use land_data_mod, only: land_state_type
 use vegn_cohort_mod, only : vegn_cohort_type
 use land_tile_mod, only : first_elmt, insert, new_land_tile_predefined
 use land_tile_mod, only : land_tile_list_type,land_tile_type,&
                           land_tile_enum_type
 use tiling_input_types_mod, only : tile_parameters_type,lake_predefined_type
 use time_manager_mod, only : time_type
 use time_interp_mod, only : time_interp

 implicit none

 interface get_parameter_data
   module procedure get_parameter_data_1d
   module procedure get_parameter_data_2d
 end interface


contains

subroutine land_cover_cold_start_0d_predefined_tiles(tiles,lnd,i,j)
  
  use netcdf
  type(land_tile_list_type),intent(inout) :: tiles
  type(land_state_type),intent(in) :: lnd
  integer,intent(in) :: i,j
  type(land_tile_type), pointer :: tile
  integer :: itile
  integer :: parent_id = 0
  integer :: ncid,status,varid,grpid,dimid,cell_grpid,cellid
  character(100) :: cellid_string
  real :: lat,lon
  real,allocatable,dimension(:) :: tmp
  type(tile_parameters_type) :: tile_parameters
  !ncid = lnd%ncid

  !Determine the lat/lon of the grid cell (degrees)
  lon = 180.0*lnd%lon(i,j)/pi
  lat = 180.0*lnd%lat(i,j)/pi

  !Open access to the model input database
  status = nf90_open('INPUT/land_model_input_database.nc', NF90_NOWRITE, ncid)

  !Determine the cell id
  status = nf90_inq_grp_ncid(ncid,"metadata",grpid)
  status = nf90_inq_varid(grpid,"mapping",varid)
  status = nf90_get_var(grpid,varid,cellid,start=(/j,i/))
  print*,i,j,cellid
  print*,lon,lat

  !Open access to the cell's group
  status = nf90_inq_grp_ncid(ncid,"grid_data",grpid)
  write(cellid_string,'(I10)') cellid
  cellid_string = trim('g' // trim(adjustl(cellid_string)))
  status = nf90_inq_grp_ncid(grpid,cellid_string,grpid)

  !Retrieve the number of tiles
  status = nf90_inq_dimid(grpid,"tile",dimid)
  status = nf90_inquire_dimension(grpid,dimid,len=tile_parameters%ntile)

  !Retrieve the number of bands
  status = nf90_inq_dimid(grpid,"band",dimid)
  status = nf90_inquire_dimension(grpid,dimid,len=tile_parameters%nband)
  tile_parameters%lake%nband = tile_parameters%nband
 
  !Allocate memory for the land container
  status = nf90_inq_grp_ncid(grpid,"parameters",tile_parameters%nc_grpid)
  status = nf90_inq_grp_ncid(grpid,"lake",tile_parameters%lake%nc_grpid)

  !Miscellanous parameters
  allocate(tile_parameters%frac(tile_parameters%ntile))
  tile_parameters%frac(:) = get_parameter_data(tile_parameters%nc_grpid,&
              "grid_cell_fraction",tile_parameters%ntile)
  !Retrieve the number of lakes
  status = nf90_inq_dimid(tile_parameters%lake%nc_grpid,"lake",dimid)
  status = nf90_inquire_dimension(tile_parameters%lake%nc_grpid,dimid,&
           len=tile_parameters%lake%nlake)
  allocate(tile_parameters%lake%frac(tile_parameters%lake%nlake))
  tile_parameters%lake%frac(:) = get_parameter_data(tile_parameters%lake%nc_grpid,&
              "frac",tile_parameters%lake%nlake)

  !Normalize the fractions (This should not be here)
  !tile_parameters%lake%frac(:) = 0.0
  tile_parameters%frac(:) = (1.0-sum(tile_parameters%lake%frac(:)))*tile_parameters%frac(:)

  !Soil and hillslope parameters
  call retrieve_soil_parameters(tile_parameters)

  !Lake parameters
  call retrieve_lake_parameters(tile_parameters%lake)

  !Vegetation parameters

  !Define the lake tiles
  do itile = 1,tile_parameters%lake%nlake
   if (tile_parameters%lake%frac(itile) .eq. 0.0)cycle
   tile => new_land_tile_predefined(frac=tile_parameters%lake%frac(itile),lake=itile,&
           lake_predefined=tile_parameters%lake,itile=itile)
   call insert(tile,tiles)
  enddo

  !Define the soil tiles
  do itile = 1,tile_parameters%ntile
   print*,tile_parameters%frac(itile)
   tile => new_land_tile_predefined(frac=tile_parameters%frac(itile),&
           soil=1,vegn=tile_parameters%vegn(itile),htag_j=tile_parameters%hidx_j(itile),&
           htag_k=tile_parameters%hidx_k(itile),&
           tile_parameters=tile_parameters,itile=itile)
   tile%parent_id = itile
   tile%ncid = ncid
   call insert(tile,tiles)
  enddo
  
  !Free up all the memory from the input database

  !Memorize some info from the land derived data type
  !lnd%ntile = tile_parameters%ntile

end subroutine

subroutine retrieve_lake_parameters(lake)

  type(lake_predefined_type),intent(inout) :: lake
  integer :: nlake,grpid,nband
  nlake = lake%nlake
  nband = lake%nband
  grpid = lake%nc_grpid

  allocate(lake%connected_to_next(nlake))
  lake%connected_to_next(:) = get_parameter_data(grpid,&
              "connected_to_next",nlake)
  allocate(lake%whole_area(nlake))
  lake%whole_area(:) = get_parameter_data(grpid,&
              "whole_lake_area",nlake)
  allocate(lake%depth_sill(nlake))
  lake%depth_sill(:) = get_parameter_data(grpid,&
              "lake_depth_sill",nlake)
  allocate(lake%width_sill(nlake))
  lake%width_sill(:) = get_parameter_data(grpid,&
              "lake_width_sill",nlake)
  allocate(lake%backwater(nlake))
  lake%backwater(:) = get_parameter_data(grpid,&
              "lake_backwater",nlake)
  allocate(lake%backwater_1(nlake))
  lake%backwater_1(:) = get_parameter_data(grpid,&
              "lake_backwater_1",nlake)
  allocate(lake%refl_dry_dir(nlake,nband))
  lake%refl_dry_dir(:,:) = get_parameter_data(grpid,&
              "refl_dry_dir",nlake,nband)
  allocate(lake%awc_lm2(nlake))
  lake%awc_lm2(:) = get_parameter_data(grpid,&
              "awc_lm2",nlake)
  allocate(lake%w_sat(nlake))
  lake%w_sat(:) = get_parameter_data(grpid,&
              "w_sat",nlake)
  allocate(lake%refl_dry_dif(nlake,nband))
  lake%refl_dry_dif(:,:) = get_parameter_data(grpid,&
              "refl_dry_dif",nlake,nband)
  allocate(lake%chb(nlake))
  lake%chb(:) = get_parameter_data(grpid,&
              "chb",nlake)
  allocate(lake%z0_momentum(nlake))
  lake%z0_momentum(:) = get_parameter_data(grpid,&
              "z0_momentum",nlake)
  allocate(lake%psi_sat_ref(nlake))
  lake%psi_sat_ref(:) = get_parameter_data(grpid,&
              "psi_sat_ref",nlake)
  allocate(lake%refl_sat_dir(nlake,nband))
  lake%refl_sat_dir(:,:) = get_parameter_data(grpid,&
              "refl_sat_dir",nlake,nband)
  allocate(lake%emis_dry(nlake))
  lake%emis_dry(:) = get_parameter_data(grpid,&
              "emis_dry",nlake)
  allocate(lake%z0_momentum_ice(nlake))
  lake%z0_momentum_ice(:) = get_parameter_data(grpid,&
              "z0_momentum_ice",nlake)
  allocate(lake%heat_capacity_ref(nlake))
  lake%heat_capacity_ref(:) = get_parameter_data(grpid,&
              "heat_capacity_ref",nlake)
  allocate(lake%refl_sat_dif(nlake,nband))
  lake%refl_sat_dif(:,:) = get_parameter_data(grpid,&
              "refl_sat_dif",nlake,nband)
  allocate(lake%alpha(nlake))
  lake%alpha(:) = get_parameter_data(grpid,&
              "alpha",nlake)
  allocate(lake%thermal_cond_ref(nlake))
  lake%thermal_cond_ref(:) = get_parameter_data(grpid,&
              "thermal_cond_ref",nlake)
  allocate(lake%emis_sat(nlake))
  lake%emis_sat(:) = get_parameter_data(grpid,&
              "emis_sat",nlake)
  allocate(lake%k_sat_ref(nlake))
  lake%k_sat_ref(:) = get_parameter_data(grpid,&
              "k_sat_ref",nlake)

end subroutine retrieve_lake_parameters

subroutine retrieve_soil_parameters(tile_parameters)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer :: ntile,grpid,nband
  ntile = tile_parameters%ntile
  grpid = tile_parameters%nc_grpid
  nband = tile_parameters%nband

  allocate(tile_parameters%dat_w_sat(ntile))
  tile_parameters%dat_w_sat(:) = get_parameter_data(grpid,&
              "dat_w_sat",ntile)
  allocate(tile_parameters%dat_awc_lm2(ntile))
  tile_parameters%dat_awc_lm2(:) = get_parameter_data(grpid,&
              "dat_awc_lm2",ntile)
  allocate(tile_parameters%dat_k_sat_ref(ntile))
  tile_parameters%dat_k_sat_ref(:) = get_parameter_data(grpid,&
              "dat_k_sat_ref",ntile)
  allocate(tile_parameters%dat_psi_sat_ref(ntile))
  tile_parameters%dat_psi_sat_ref(:) = get_parameter_data(grpid,&
              "dat_psi_sat_ref",ntile)
  allocate(tile_parameters%dat_chb(ntile))
  tile_parameters%dat_chb(:) = get_parameter_data(grpid,&
              "dat_chb",ntile)
  allocate(tile_parameters%dat_heat_capacity_dry(ntile))
  tile_parameters%dat_heat_capacity_dry(:) = get_parameter_data(grpid,&
              "dat_heat_capacity_dry",ntile)
  allocate(tile_parameters%dat_thermal_cond_dry(ntile))
  tile_parameters%dat_thermal_cond_dry(:) = get_parameter_data(grpid,&
              "dat_thermal_cond_dry",ntile)
  allocate(tile_parameters%dat_thermal_cond_sat(ntile))
  tile_parameters%dat_thermal_cond_sat(:) = get_parameter_data(grpid,&
              "dat_thermal_cond_sat",ntile)
  allocate(tile_parameters%dat_thermal_cond_exp(ntile))
  tile_parameters%dat_thermal_cond_exp(:) = get_parameter_data(grpid,&
              "dat_thermal_cond_exp",ntile)
  allocate(tile_parameters%dat_thermal_cond_scale(ntile))
  tile_parameters%dat_thermal_cond_scale(:) = get_parameter_data(grpid,&
              "dat_thermal_cond_scale",ntile)
  allocate(tile_parameters%dat_thermal_cond_weight(ntile))
  tile_parameters%dat_thermal_cond_weight(:) = get_parameter_data(grpid,&
              "dat_thermal_cond_weight",ntile)
  allocate(tile_parameters%dat_refl_dry_dir(ntile,nband))
  tile_parameters%dat_refl_dry_dir(:,:) = get_parameter_data(grpid,&
              "dat_refl_dry_dir",ntile,nband)
  allocate(tile_parameters%dat_refl_dry_dif(ntile,nband))
  tile_parameters%dat_refl_dry_dif(:,:) = get_parameter_data(grpid,&
              "dat_refl_dry_dif",ntile,nband)
  allocate(tile_parameters%dat_refl_sat_dir(ntile,nband))
  tile_parameters%dat_refl_sat_dir(:,:) = get_parameter_data(grpid,&
              "dat_refl_sat_dir",ntile,nband)
  allocate(tile_parameters%dat_refl_sat_dif(ntile,nband))
  tile_parameters%dat_refl_sat_dif(:,:) = get_parameter_data(grpid,&
              "dat_refl_sat_dif",ntile,nband)
  allocate(tile_parameters%dat_emis_dry(ntile))
  tile_parameters%dat_emis_dry(:) = get_parameter_data(grpid,&
              "dat_emis_dry",ntile)
  allocate(tile_parameters%dat_emis_sat(ntile))
  tile_parameters%dat_emis_sat(:) = get_parameter_data(grpid,&
              "dat_emis_sat",ntile)
  allocate(tile_parameters%dat_z0_momentum(ntile))
  tile_parameters%dat_z0_momentum(:) = get_parameter_data(grpid,&
              "dat_z0_momentum",ntile)
  allocate(tile_parameters%dat_tf_depr(ntile))
  tile_parameters%dat_tf_depr(:) = get_parameter_data(grpid,&
              "dat_tf_depr",ntile)
  allocate(tile_parameters%rsa_exp_global(ntile))
  tile_parameters%rsa_exp_global(:) = get_parameter_data(grpid,&
              "rsa_exp_global",ntile)
  allocate(tile_parameters%gw_res_time(ntile))
  tile_parameters%gw_res_time(:) = get_parameter_data(grpid,&
              "gw_res_time",ntile)
  allocate(tile_parameters%gw_hillslope_length(ntile))
  tile_parameters%gw_hillslope_length(:) = get_parameter_data(grpid,&
              "gw_hillslope_length",ntile)
  allocate(tile_parameters%gw_scale_length(ntile))
  tile_parameters%gw_scale_length(:) = get_parameter_data(grpid,&
              "gw_scale_length",ntile)
  allocate(tile_parameters%gw_hillslope_zeta_bar(ntile))
  tile_parameters%gw_hillslope_zeta_bar(:) = get_parameter_data(grpid,&
              "gw_hillslope_zeta_bar",ntile)
  allocate(tile_parameters%gw_hillslope_relief(ntile))
  tile_parameters%gw_hillslope_relief(:) = get_parameter_data(grpid,&
              "gw_hillslope_relief",ntile)
  allocate(tile_parameters%gw_scale_relief(ntile))
  tile_parameters%gw_scale_relief(:) = get_parameter_data(grpid,&
              "gw_scale_relief",ntile)
  allocate(tile_parameters%gw_soil_e_depth(ntile))
  tile_parameters%gw_soil_e_depth(:) = get_parameter_data(grpid,&
              "gw_soil_e_depth",ntile)
  allocate(tile_parameters%gw_scale_soil_depth(ntile))
  tile_parameters%gw_scale_soil_depth(:) = get_parameter_data(grpid,&
              "gw_scale_soil_depth",ntile)
  allocate(tile_parameters%gw_hillslope_a(ntile))
  tile_parameters%gw_hillslope_a(:) = get_parameter_data(grpid,&
              "gw_hillslope_a",ntile)
  allocate(tile_parameters%gw_hillslope_n(ntile))
  tile_parameters%gw_hillslope_n(:) = get_parameter_data(grpid,&
              "gw_hillslope_n",ntile)
  allocate(tile_parameters%gw_perm(ntile))
  tile_parameters%gw_perm(:) = get_parameter_data(grpid,&
              "gw_perm",ntile)
  allocate(tile_parameters%gw_scale_perm(ntile))
  tile_parameters%gw_scale_perm(:) = get_parameter_data(grpid,&
              "gw_scale_perm",ntile)
  allocate(tile_parameters%microtopo(ntile))
  tile_parameters%microtopo(:) = get_parameter_data(grpid,&
              "microtopo",ntile)
  allocate(tile_parameters%tile_hlsp_length(ntile))
  tile_parameters%tile_hlsp_length(:) = get_parameter_data(grpid,&
              "tile_hlsp_length",ntile)
  allocate(tile_parameters%tile_hlsp_slope(ntile))
  tile_parameters%tile_hlsp_slope(:) = get_parameter_data(grpid,&
              "tile_hlsp_slope",ntile)
  allocate(tile_parameters%tile_hlsp_elev(ntile))
  tile_parameters%tile_hlsp_elev(:) = get_parameter_data(grpid,&
              "tile_hlsp_elev",ntile)
  allocate(tile_parameters%tile_hlsp_hpos(ntile))
  tile_parameters%tile_hlsp_hpos(:) = get_parameter_data(grpid,&
              "tile_hlsp_hpos",ntile)
  allocate(tile_parameters%tile_hlsp_width(ntile))
  tile_parameters%tile_hlsp_width(:) = get_parameter_data(grpid,&
              "tile_hlsp_width",ntile)
  allocate(tile_parameters%hidx_k(ntile))
  tile_parameters%hidx_k(:) = get_parameter_data(grpid,&
              "hidx_k",ntile)
  allocate(tile_parameters%hidx_j(ntile))
  tile_parameters%hidx_j(:) = get_parameter_data(grpid,&
              "hidx_j",ntile)
  allocate(tile_parameters%vegn(ntile))
  tile_parameters%vegn(:) = get_parameter_data(grpid,&
              "vegn",ntile)

end subroutine retrieve_soil_parameters

function get_parameter_data_1d(grpid,var,nx) result(tmp)

 use netcdf
 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx
 real :: tmp(nx)
 integer :: itile,varid,status
 
 status = nf90_inq_varid(grpid,var,varid)
 status = nf90_get_var(grpid,varid,tmp)

end function

function get_parameter_data_2d(grpid,var,nx,ny) result(tmp)

 use netcdf
 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx,ny
 real :: tmp(nx,ny)
 integer :: itile,varid,status

 status = nf90_inq_varid(grpid,var,varid)
 status = nf90_get_var(grpid,varid,tmp)

end function

subroutine set_tile_meteorology(tile,itime)

  use netcdf
  type(land_tile_type), pointer :: tile
  integer :: status,itime,varid,grpid
  real :: tmp

  !Read in and set the meteorology for each tile
  !status = nf90_inq_grp_ncid(tile%ncid,"meteorology",grpid)
  !status = nf90_inq_varid(grpid,"prec",varid)
  !status = nf90_get_var(grpid,varid,tmp,start=(/tile%parent_id,itime/))
  !tile%cana%prec = tmp/3600.0

end subroutine 

subroutine read_static_vegn_nwc(tile,itime)
 
 use netcdf
 type(land_tile_type), pointer :: tile
 type(vegn_cohort_type), pointer :: cohort
 integer :: status,itime,varid,grpid

 !Define the cohort
 cohort => tile%vegn%cohorts(1)

 !Read in the vegetation data for each tile
 status = nf90_inq_grp_ncid(tile%ncid,"vegetation",grpid)
 status = nf90_inq_varid(grpid,"bl",varid)
 status = nf90_get_var(grpid,varid,cohort%bl,start=(/tile%parent_id,itime/))
 status = nf90_inq_varid(grpid,"blv",varid)
 status = nf90_get_var(grpid,varid,cohort%blv,start=(/tile%parent_id,itime/))
 status = nf90_inq_varid(grpid,"br",varid)
 status = nf90_get_var(grpid,varid,cohort%br,start=(/tile%parent_id,itime/))
 status = nf90_inq_varid(grpid,"bsw",varid)
 status = nf90_get_var(grpid,varid,cohort%bsw,start=(/tile%parent_id,itime/))
 status = nf90_inq_varid(grpid,"bwood",varid)
 status = nf90_get_var(grpid,varid,cohort%bwood,start=(/tile%parent_id,itime/))
 status = nf90_inq_varid(grpid,"bliving",varid)
 status = nf90_get_var(grpid,varid,cohort%bliving,start=(/tile%parent_id,itime/))
 status = nf90_inq_varid(grpid,"status",varid)
 status = nf90_get_var(grpid,varid,cohort%status,start=(/tile%parent_id,itime/))
 status = nf90_inq_varid(grpid,"species",varid)
 status = nf90_get_var(grpid,varid,cohort%species,start=(/tile%parent_id,itime/))

end subroutine

end module predefined_tiles_mod
