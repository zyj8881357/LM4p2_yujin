module predefined_tiles_mod

 use netcdf
 use constants_mod     , only : pi
 use land_data_mod, only: land_state_type
 use vegn_cohort_mod, only : vegn_cohort_type
 use land_tile_mod, only : first_elmt, insert, new_land_tile_predefined
 use land_tile_mod, only : land_tile_list_type,land_tile_type,&
                           land_tile_enum_type
 use tiling_input_types_mod, only : tile_parameters_type,lake_predefined_type
 use tiling_input_types_mod, only : glacier_predefined_type,soil_predefined_type
 use tiling_input_types_mod, only : metadata_predefined_type
 use time_manager_mod, only : time_type
 use time_interp_mod, only : time_interp

 implicit none

 interface get_parameter_data
   module procedure get_parameter_data_1d_integer
   module procedure get_parameter_data_1d_real
   module procedure get_parameter_data_2d
 end interface

contains

subroutine land_cover_cold_start_0d_predefined_tiles(tiles,lnd,i,j)
  
  !use netcdf
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
  !status = nf90_inq_dimid(grpid,"tile",dimid)
  !status = nf90_inquire_dimension(grpid,dimid,len=tile_parameters%soil%nsoil)

  !Retrieve the number of bands
  !status = nf90_inq_dimid(grpid,"band",dimid)
  !status = nf90_inquire_dimension(grpid,dimid,len=tile_parameters%soil%nband)
  !tile_parameters%lake%nband = tile_parameters%soil%nband
  !tile_parameters%glacier%nband = tile_parameters%soil%nband
 
  !Allocate memory for the land container
  !status = nf90_inq_grp_ncid(grpid,"soil",tile_parameters%soil%nc_grpid)
  !status = nf90_inq_grp_ncid(grpid,"lake",tile_parameters%lake%nc_grpid)
  !status = nf90_inq_grp_ncid(grpid,"glacier",tile_parameters%glacier%nc_grpid)

  !Miscellanous parameters
  !allocate(tile_parameters%soil%frac(tile_parameters%soil%nsoil))
  !tile_parameters%soil%frac(:) = get_parameter_data(tile_parameters%soil%nc_grpid,&
  !            "grid_cell_fraction",tile_parameters%soil%nsoil)
  !Retrieve the number of lakes
  !status = nf90_inq_dimid(tile_parameters%lake%nc_grpid,"lake",dimid)
  !status = nf90_inquire_dimension(tile_parameters%lake%nc_grpid,dimid,&
  !         len=tile_parameters%lake%nlake)
  !allocate(tile_parameters%lake%frac(tile_parameters%lake%nlake))
  !tile_parameters%lake%frac(:) = get_parameter_data(tile_parameters%lake%nc_grpid,&
  !            "frac",tile_parameters%lake%nlake)
  !Retrieve the number of glaciers
  !status = nf90_inq_dimid(tile_parameters%glacier%nc_grpid,"glacier",dimid)
  !status = nf90_inquire_dimension(tile_parameters%glacier%nc_grpid,dimid,&
  !         len=tile_parameters%glacier%nglacier)
  !allocate(tile_parameters%glacier%frac(tile_parameters%glacier%nglacier))
  !tile_parameters%glacier%frac(:) = get_parameter_data(tile_parameters%glacier%nc_grpid,&
  !            "frac",tile_parameters%glacier%nglacier)

  !Normalize the fractions (This should not be here)
  !tile_parameters%lake%frac(:) = 0.0
  !tile_parameters%soil%frac(:) = (1.0-(sum(tile_parameters%lake%frac(:))+sum(tile_parameters%glacier%frac(:))))*tile_parameters%soil%frac(:)

  !Metadata
  call retrieve_metadata(tile_parameters,grpid)

  !Soil and hillslope parameters
  !call retrieve_soil_parameters(tile_parameters%soil)

  !Lake parameters
  call retrieve_lake_parameters(tile_parameters,grpid)

  !Glacier parameters
  call retrieve_glacier_parameters(tile_parameters,grpid)
  stop

  !Vegetation parameters

  !Define the glacier tiles
  do itile = 1,tile_parameters%glacier%nglacier
   if (tile_parameters%glacier%frac(itile) .eq. 0.0)cycle
   print*,'glac',tile_parameters%glacier%frac(itile)
   tile => new_land_tile_predefined(frac=tile_parameters%glacier%frac(itile),glac=itile,&
           glacier_predefined=tile_parameters%glacier,itile=itile)
   call insert(tile,tiles)
  enddo

  !Define the lake tiles
  do itile = 1,tile_parameters%lake%nlake
   if (tile_parameters%lake%frac(itile) .eq. 0.0)cycle
   print*,'lake',tile_parameters%lake%frac(itile)
   tile => new_land_tile_predefined(frac=tile_parameters%lake%frac(itile),lake=itile,&
           lake_predefined=tile_parameters%lake,itile=itile)
   call insert(tile,tiles)
  enddo

  !Define the soil tiles
  do itile = 1,tile_parameters%soil%nsoil
   print*,tile_parameters%soil%frac(itile)
   print*,'soil',tile_parameters%soil%frac(itile)
   tile => new_land_tile_predefined(frac=tile_parameters%soil%frac(itile),&
           soil=1,vegn=tile_parameters%soil%vegn(itile),&
           htag_j=tile_parameters%soil%hidx_j(itile),&
           htag_k=tile_parameters%soil%hidx_k(itile),&
           soil_predefined=tile_parameters%soil,itile=itile)
   call insert(tile,tiles)
  enddo
  
end subroutine

subroutine retrieve_metadata(tile_parameters,cid)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer,intent(inout) :: cid
  type(metadata_predefined_type),pointer :: metadata
  integer :: dimid,grpid,ntile,nband,status
  allocate(tile_parameters%metadata)
  metadata => tile_parameters%metadata

  !Retrieve the number of tiles
  status = nf90_inq_dimid(cid,"tile",dimid)
  status = nf90_inquire_dimension(cid,dimid,len=ntile)
  metadata%ntile = ntile

  !Retrieve the number of bands
  status = nf90_inq_dimid(cid,"band",dimid)
  status = nf90_inquire_dimension(cid,dimid,len=nband)
  metadata%nband = nband

  !Retrieve the group id
  status = nf90_inq_grp_ncid(cid,"metadata",grpid)

  !Retrieve the info
  call get_parameter_data(grpid,'frac',ntile,metadata%frac)
  call get_parameter_data(grpid,'tile',ntile,metadata%tile)
  call get_parameter_data(grpid,'type',ntile,metadata%ttype)

end subroutine retrieve_metadata

subroutine retrieve_glacier_parameters(tile_parameters,cid)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer,intent(inout) :: cid
  type(glacier_predefined_type),pointer :: glacier
  integer :: dimid,grpid,nglacier,nband,status
  allocate(tile_parameters%glacier)
  glacier => tile_parameters%glacier

  !Retrieve the group id
  status = nf90_inq_grp_ncid(cid,"glacier",grpid)

  !Retrieve the number of glaciers
  status = nf90_inq_dimid(grpid,"glacier",dimid)
  status = nf90_inquire_dimension(grpid,dimid,len=nglacier)
  glacier%nglacier = nglacier

  !Retrieve the parameters
  call get_parameter_data(grpid,"wsat",nglacier,glacier%w_sat)
  call get_parameter_data(grpid,"awc_lm2",nglacier,glacier%awc_lm2)
  call get_parameter_data(grpid,"k_sat_ref",nglacier,glacier%k_sat_ref)
  call get_parameter_data(grpid,"psi_sat_ref",nglacier,glacier%psi_sat_ref)
  call get_parameter_data(grpid,"chb",nglacier,glacier%chb)
  call get_parameter_data(grpid,"alpha",nglacier,glacier%alpha)
  call get_parameter_data(grpid,"heat_capacity_ref",nglacier,glacier%heat_capacity_ref)
  call get_parameter_data(grpid,"thermal_cond_ref",nglacier,glacier%thermal_cond_ref)
  call get_parameter_data(grpid,"emis_dry",nglacier,glacier%emis_dry)
  call get_parameter_data(grpid,"emis_sat",nglacier,glacier%emis_sat)
  call get_parameter_data(grpid,"z0_momentum",nglacier,glacier%z0_momentum)
  call get_parameter_data(grpid,"tfreeze",nglacier,glacier%tfreeze)
  call get_parameter_data(grpid,"refl_max_dir",nglacier,nband,glacier%refl_max_dir)
  call get_parameter_data(grpid,"refl_max_dif",nglacier,nband,glacier%refl_max_dif)
  call get_parameter_data(grpid,"refl_min_dir",nglacier,nband,glacier%refl_min_dir)
  call get_parameter_data(grpid,"refl_min_dif",nglacier,nband,glacier%refl_min_dif)

end subroutine retrieve_glacier_parameters

subroutine retrieve_lake_parameters(tile_parameters,cid)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer,intent(inout) :: cid
  type(lake_predefined_type),pointer :: lake
  integer :: dimid,grpid,nlake,nband,status
  allocate(tile_parameters%lake)
  lake => tile_parameters%lake

  !Retrieve the group id
  status = nf90_inq_grp_ncid(cid,"lake",grpid)

  !Retrieve the number of glaciers
  status = nf90_inq_dimid(grpid,"lake",dimid)
  status = nf90_inquire_dimension(grpid,dimid,len=nlake)
  lake%nlake = nlake

  !Retrieve the parameters
  call get_parameter_data(grpid,"connected_to_next",nlake,lake%connected_to_next)
  call get_parameter_data(grpid,"whole_lake_area",nlake,lake%whole_area)
  call get_parameter_data(grpid,"lake_depth_sill",nlake,lake%depth_sill)
  call get_parameter_data(grpid,"lake_width_sill",nlake,lake%width_sill)
  call get_parameter_data(grpid,"lake_backwater",nlake,lake%backwater)
  call get_parameter_data(grpid,"lake_backwater_1",nlake,lake%backwater_1)
  call get_parameter_data(grpid,"awc_lm2",nlake,lake%awc_lm2)
  call get_parameter_data(grpid,"w_sat",nlake,lake%w_sat)
  call get_parameter_data(grpid,"chb",nlake,lake%chb)
  call get_parameter_data(grpid,"refl_dry_dif",nlake,nband,lake%refl_dry_dif)
  call get_parameter_data(grpid,"refl_dry_dir",nlake,nband,lake%refl_dry_dir)
  call get_parameter_data(grpid,"refl_sat_dif",nlake,nband,lake%refl_sat_dif)
  call get_parameter_data(grpid,"refl_sat_dir",nlake,nband,lake%refl_sat_dir)
  call get_parameter_data(grpid,"z0_momentum",nlake,lake%z0_momentum)
  call get_parameter_data(grpid,"psi_sat_ref",nlake,lake%psi_sat_ref)
  call get_parameter_data(grpid,"emis_dry",nlake,lake%emis_dry)
  call get_parameter_data(grpid,"z0_momentum_ice",nlake,lake%z0_momentum_ice)
  call get_parameter_data(grpid,"heat_capacity_ref",nlake,lake%heat_capacity_ref)
  call get_parameter_data(grpid,"alpha",nlake,lake%alpha)
  call get_parameter_data(grpid,"thermal_cond_ref",nlake,lake%thermal_cond_ref)
  call get_parameter_data(grpid,"emis_sat",nlake,lake%emis_sat)
  call get_parameter_data(grpid,"k_sat_ref",nlake,lake%k_sat_ref)
!
end subroutine retrieve_lake_parameters
!
!subroutine retrieve_soil_parameters(tile_parameters)
!
!  type(soil_predefined_type),intent(inout) :: tile_parameters
!  integer :: ntile,grpid,nband
!  ntile = tile_parameters%nsoil
!  grpid = tile_parameters%nc_grpid
!  nband = tile_parameters%nband
!
!  allocate(tile_parameters%dat_w_sat(ntile))
!  tile_parameters%dat_w_sat(:) = get_parameter_data(grpid,&
!              "dat_w_sat",ntile)
!  allocate(tile_parameters%dat_awc_lm2(ntile))
!  tile_parameters%dat_awc_lm2(:) = get_parameter_data(grpid,&
!              "dat_awc_lm2",ntile)
!  allocate(tile_parameters%dat_k_sat_ref(ntile))
!  tile_parameters%dat_k_sat_ref(:) = get_parameter_data(grpid,&
!              "dat_k_sat_ref",ntile)
!  allocate(tile_parameters%dat_psi_sat_ref(ntile))
!  tile_parameters%dat_psi_sat_ref(:) = get_parameter_data(grpid,&
!              "dat_psi_sat_ref",ntile)
!  allocate(tile_parameters%dat_chb(ntile))
!  tile_parameters%dat_chb(:) = get_parameter_data(grpid,&
!              "dat_chb",ntile)
!  allocate(tile_parameters%dat_heat_capacity_dry(ntile))
!  tile_parameters%dat_heat_capacity_dry(:) = get_parameter_data(grpid,&
!              "dat_heat_capacity_dry",ntile)
!  allocate(tile_parameters%dat_thermal_cond_dry(ntile))
!  tile_parameters%dat_thermal_cond_dry(:) = get_parameter_data(grpid,&
!              "dat_thermal_cond_dry",ntile)
!  allocate(tile_parameters%dat_thermal_cond_sat(ntile))
!  tile_parameters%dat_thermal_cond_sat(:) = get_parameter_data(grpid,&
!              "dat_thermal_cond_sat",ntile)
!  allocate(tile_parameters%dat_thermal_cond_exp(ntile))
!  tile_parameters%dat_thermal_cond_exp(:) = get_parameter_data(grpid,&
!              "dat_thermal_cond_exp",ntile)
!  allocate(tile_parameters%dat_thermal_cond_scale(ntile))
!  tile_parameters%dat_thermal_cond_scale(:) = get_parameter_data(grpid,&
!              "dat_thermal_cond_scale",ntile)
!  allocate(tile_parameters%dat_thermal_cond_weight(ntile))
!  tile_parameters%dat_thermal_cond_weight(:) = get_parameter_data(grpid,&
!              "dat_thermal_cond_weight",ntile)
!  allocate(tile_parameters%dat_refl_dry_dir(ntile,nband))
!  tile_parameters%dat_refl_dry_dir(:,:) = get_parameter_data(grpid,&
!              "dat_refl_dry_dir",ntile,nband)
!  allocate(tile_parameters%dat_refl_dry_dif(ntile,nband))
!  tile_parameters%dat_refl_dry_dif(:,:) = get_parameter_data(grpid,&
!              "dat_refl_dry_dif",ntile,nband)
!  allocate(tile_parameters%dat_refl_sat_dir(ntile,nband))
!  tile_parameters%dat_refl_sat_dir(:,:) = get_parameter_data(grpid,&
!              "dat_refl_sat_dir",ntile,nband)
!  allocate(tile_parameters%dat_refl_sat_dif(ntile,nband))
!  tile_parameters%dat_refl_sat_dif(:,:) = get_parameter_data(grpid,&
!              "dat_refl_sat_dif",ntile,nband)
!  allocate(tile_parameters%dat_emis_dry(ntile))
!  tile_parameters%dat_emis_dry(:) = get_parameter_data(grpid,&
!              "dat_emis_dry",ntile)
!  allocate(tile_parameters%dat_emis_sat(ntile))
!  tile_parameters%dat_emis_sat(:) = get_parameter_data(grpid,&
!              "dat_emis_sat",ntile)
!  allocate(tile_parameters%dat_z0_momentum(ntile))
!  tile_parameters%dat_z0_momentum(:) = get_parameter_data(grpid,&
!              "dat_z0_momentum",ntile)
!  allocate(tile_parameters%dat_tf_depr(ntile))
!  tile_parameters%dat_tf_depr(:) = get_parameter_data(grpid,&
!              "dat_tf_depr",ntile)
!  allocate(tile_parameters%rsa_exp_global(ntile))
!  tile_parameters%rsa_exp_global(:) = get_parameter_data(grpid,&
!              "rsa_exp_global",ntile)
!  allocate(tile_parameters%gw_res_time(ntile))
!  tile_parameters%gw_res_time(:) = get_parameter_data(grpid,&
!              "gw_res_time",ntile)
!  allocate(tile_parameters%gw_hillslope_length(ntile))
!  tile_parameters%gw_hillslope_length(:) = get_parameter_data(grpid,&
!              "gw_hillslope_length",ntile)
!  allocate(tile_parameters%gw_scale_length(ntile))
!  tile_parameters%gw_scale_length(:) = get_parameter_data(grpid,&
!              "gw_scale_length",ntile)
!  allocate(tile_parameters%gw_hillslope_zeta_bar(ntile))
!  tile_parameters%gw_hillslope_zeta_bar(:) = get_parameter_data(grpid,&
!              "gw_hillslope_zeta_bar",ntile)
!  allocate(tile_parameters%gw_hillslope_relief(ntile))
!  tile_parameters%gw_hillslope_relief(:) = get_parameter_data(grpid,&
!              "gw_hillslope_relief",ntile)
!  allocate(tile_parameters%gw_scale_relief(ntile))
!  tile_parameters%gw_scale_relief(:) = get_parameter_data(grpid,&
!              "gw_scale_relief",ntile)
!  allocate(tile_parameters%gw_soil_e_depth(ntile))
!  tile_parameters%gw_soil_e_depth(:) = get_parameter_data(grpid,&
!              "gw_soil_e_depth",ntile)
!  allocate(tile_parameters%gw_scale_soil_depth(ntile))
!  tile_parameters%gw_scale_soil_depth(:) = get_parameter_data(grpid,&
!              "gw_scale_soil_depth",ntile)
!  allocate(tile_parameters%gw_hillslope_a(ntile))
!  tile_parameters%gw_hillslope_a(:) = get_parameter_data(grpid,&
!              "gw_hillslope_a",ntile)
!  allocate(tile_parameters%gw_hillslope_n(ntile))
!  tile_parameters%gw_hillslope_n(:) = get_parameter_data(grpid,&
!              "gw_hillslope_n",ntile)
!  allocate(tile_parameters%gw_perm(ntile))
!  tile_parameters%gw_perm(:) = get_parameter_data(grpid,&
!              "gw_perm",ntile)
!  allocate(tile_parameters%gw_scale_perm(ntile))
!  tile_parameters%gw_scale_perm(:) = get_parameter_data(grpid,&
!              "gw_scale_perm",ntile)
!  allocate(tile_parameters%microtopo(ntile))
!  tile_parameters%microtopo(:) = get_parameter_data(grpid,&
!              "microtopo",ntile)
!  allocate(tile_parameters%tile_hlsp_length(ntile))
!  tile_parameters%tile_hlsp_length(:) = get_parameter_data(grpid,&
!              "tile_hlsp_length",ntile)
!  allocate(tile_parameters%tile_hlsp_slope(ntile))
!  tile_parameters%tile_hlsp_slope(:) = get_parameter_data(grpid,&
!              "tile_hlsp_slope",ntile)
!  allocate(tile_parameters%tile_hlsp_elev(ntile))
!  tile_parameters%tile_hlsp_elev(:) = get_parameter_data(grpid,&
!              "tile_hlsp_elev",ntile)
!  allocate(tile_parameters%tile_hlsp_hpos(ntile))
!  tile_parameters%tile_hlsp_hpos(:) = get_parameter_data(grpid,&
!              "tile_hlsp_hpos",ntile)
!  allocate(tile_parameters%tile_hlsp_width(ntile))
!  tile_parameters%tile_hlsp_width(:) = get_parameter_data(grpid,&
!              "tile_hlsp_width",ntile)
!  allocate(tile_parameters%hidx_k(ntile))
!  tile_parameters%hidx_k(:) = get_parameter_data(grpid,&
!              "hidx_k",ntile)
!  allocate(tile_parameters%hidx_j(ntile))
!  tile_parameters%hidx_j(:) = get_parameter_data(grpid,&
!              "hidx_j",ntile)
!  allocate(tile_parameters%vegn(ntile))
!  tile_parameters%vegn(:) = get_parameter_data(grpid,&
!              "vegn",ntile)
!
!end subroutine retrieve_soil_parameters

subroutine get_parameter_data_1d_integer(grpid,var,nx,tmp)

 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx
 integer,dimension(:),pointer :: tmp
 integer :: itile,varid,status
 allocate(tmp(nx))

 status = nf90_inq_varid(grpid,var,varid)
 status = nf90_get_var(grpid,varid,tmp)

end subroutine

subroutine get_parameter_data_1d_real(grpid,var,nx,tmp)

 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx
 real,dimension(:),pointer :: tmp
 integer :: itile,varid,status
 allocate(tmp(nx))
 
 status = nf90_inq_varid(grpid,var,varid)
 status = nf90_get_var(grpid,varid,tmp)

end subroutine 

subroutine get_parameter_data_2d(grpid,var,nx,ny,tmp)

 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx,ny
 real :: tmp(nx,ny)
 integer :: itile,varid,status

 status = nf90_inq_varid(grpid,var,varid)
 status = nf90_get_var(grpid,varid,tmp)

end subroutine

end module predefined_tiles_mod
