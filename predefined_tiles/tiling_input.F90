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
   module procedure get_parameter_data_2d_integer
   module procedure get_parameter_data_2d_real
 end interface

contains

subroutine open_database_predefined_tiles(ncid)

  integer :: status,ncid
  !Open access to the model input database
  status = nf90_open('INPUT/land_model_input_database.nc', NF90_NOWRITE, ncid)

end subroutine open_database_predefined_tiles

subroutine close_database_predefined_tiles(ncid)

  integer :: status,ncid
  !Close access to the model input database
  status = nf90_close(ncid)

end subroutine close_database_predefined_tiles

subroutine land_cover_cold_start_0d_predefined_tiles(tiles,lnd,i,j,ncid)
  
  type(land_tile_list_type),intent(inout) :: tiles
  type(land_state_type),intent(in) :: lnd
  integer,intent(in) :: i,j,ncid
  type(land_tile_type), pointer :: tile
  integer :: itile,tid,is,js
  integer :: parent_id = 0
  integer :: status,varid,grpid,dimid,cell_grpid,cellid
  character(100) :: cellid_string
  real :: lat,lon,t0,t1
  real,allocatable,dimension(:) :: tmp
  type(tile_parameters_type) :: tile_parameters

  !Determine the lat/lon of the grid cell (degrees)
  is = i+lnd%is-1
  js = j+lnd%js-1
  lon = 180.0*lnd%lon(is,js)/pi
  lat = 180.0*lnd%lat(is,js)/pi

  !call cpu_time(t0)
  !Open access to the model input database
  !status = nf90_open('INPUT/land_model_input_database.nc', NF90_NOWRITE, ncid)
  !call cpu_time(t1)
  !print*,'opening the file',t1-t0

  !call cpu_time(t0)
  !Determine the cell id
  status = nf90_inq_grp_ncid(ncid,"metadata",grpid)
  status = nf90_inq_varid(grpid,"mapping",varid)
  !status = nf90_get_var(grpid,varid,cellid,start=(/j,i/))
  status = nf90_get_var(grpid,varid,cellid,start=(/js,is/))
  !print*,i,j,cellid
  !print*,is,js,cellid
  print*,lon,lat

  !Open access to the cell's group
  status = nf90_inq_grp_ncid(ncid,"grid_data",grpid)
  write(cellid_string,'(I10)') cellid
  cellid_string = trim('g' // trim(adjustl(cellid_string)))
  status = nf90_inq_grp_ncid(grpid,cellid_string,grpid)
  !call cpu_time(t1)
  !print*,'retrieving metadata',t1-t0 

  !call cpu_time(t0)
  !Metadata
  call retrieve_metadata(tile_parameters,grpid)

  !Soil and hillslope parameters
  call retrieve_soil_parameters(tile_parameters,grpid)

  !Lake parameters
  call retrieve_lake_parameters(tile_parameters,grpid)

  !Glacier parameters
  call retrieve_glacier_parameters(tile_parameters,grpid)
  !call cpu_time(t1)
  !print*,'retrieving parameters',t1-t0

  !Vegetation parameters

  !Create the tiles
  do itile = 1,tile_parameters%metadata%ntile
   if (tile_parameters%metadata%frac(itile) .eq. 0.0)cycle
   tid = tile_parameters%metadata%tid(itile)
   select case (tile_parameters%metadata%ttype(itile))
    case(1)
     tile => new_land_tile_predefined(frac=tile_parameters%metadata%frac(itile),&
            glac=tid,glacier_predefined=tile_parameters%glacier,&
            itile=tid)
    case(2)
     tile => new_land_tile_predefined(frac=tile_parameters%metadata%frac(itile),&
            lake=tid,lake_predefined=tile_parameters%lake,&
            itile=tid)
    case(3)
     tile => new_land_tile_predefined(frac=tile_parameters%metadata%frac(itile),&
           soil=1,vegn=tile_parameters%soil%vegn(tid),&
           htag_j=tile_parameters%soil%hidx_j(tid),&
           htag_k=tile_parameters%soil%hidx_k(tid),&
           soil_predefined=tile_parameters%soil,itile=tid)
   end select
   call insert(tile,tiles)
  enddo

  !Close the netcdf file
  !status = nf90_close(ncid)
  
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
  call get_parameter_data(grpid,'tid',ntile,metadata%tid)

  !Clean up (This should go in the database creation)
  where ((metadata%frac .lt. 1.e-8) .and. (metadata%ttype .eq. 2))metadata%frac = 0.0
  metadata%frac = metadata%frac/sum(metadata%frac)

end subroutine retrieve_metadata

subroutine retrieve_glacier_parameters(tile_parameters,cid)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer,intent(inout) :: cid
  type(glacier_predefined_type),pointer :: glacier
  integer :: dimid,grpid,nglacier,nband,status
  allocate(tile_parameters%glacier)
  glacier => tile_parameters%glacier
  nband = tile_parameters%metadata%nband

  !Retrieve the group id
  status = nf90_inq_grp_ncid(cid,"glacier",grpid)

  !Retrieve the number of glaciers
  status = nf90_inq_dimid(grpid,"glacier",dimid)
  status = nf90_inquire_dimension(grpid,dimid,len=nglacier)
  glacier%nglacier = nglacier

  !Retrieve the parameters
  call get_parameter_data(grpid,"w_sat",nglacier,glacier%w_sat)
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
  nband = tile_parameters%metadata%nband

  !Retrieve the group id
  status = nf90_inq_grp_ncid(cid,"lake",grpid)

  !Retrieve the number of lakes
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

end subroutine retrieve_lake_parameters

subroutine retrieve_soil_parameters(tile_parameters,cid)

  type(tile_parameters_type),intent(inout) :: tile_parameters
  integer,intent(inout) :: cid
  type(soil_predefined_type),pointer :: soil
  integer :: dimid,grpid,nsoil,nband,status
  allocate(tile_parameters%soil)
  soil => tile_parameters%soil
  nband = tile_parameters%metadata%nband

  !Retrieve the group id
  status = nf90_inq_grp_ncid(cid,"soil",grpid)

  !Retrieve the number of soil tiles
  status = nf90_inq_dimid(grpid,"soil",dimid)
  status = nf90_inquire_dimension(grpid,dimid,len=nsoil)
  soil%nsoil = nsoil

  !Retrieve the parameters
  call get_parameter_data(grpid,"dat_w_sat",nsoil,soil%dat_w_sat)
  call get_parameter_data(grpid,"dat_awc_lm2",nsoil,soil%dat_awc_lm2)
  call get_parameter_data(grpid,"dat_k_sat_ref",nsoil,soil%dat_k_sat_ref)
  call get_parameter_data(grpid,"dat_psi_sat_ref",nsoil,soil%dat_psi_sat_ref)
  call get_parameter_data(grpid,"dat_chb",nsoil,soil%dat_chb)
  call get_parameter_data(grpid,"dat_heat_capacity_dry",nsoil,soil%dat_heat_capacity_dry)
  call get_parameter_data(grpid,"dat_thermal_cond_dry",nsoil,soil%dat_thermal_cond_dry)
  call get_parameter_data(grpid,"dat_thermal_cond_sat",nsoil,soil%dat_thermal_cond_sat)
  call get_parameter_data(grpid,"dat_thermal_cond_exp",nsoil,soil%dat_thermal_cond_exp)
  call get_parameter_data(grpid,"dat_thermal_cond_scale",nsoil,soil%dat_thermal_cond_scale)
  call get_parameter_data(grpid,"dat_thermal_cond_weight",nsoil,soil%dat_thermal_cond_weight)
  call get_parameter_data(grpid,"dat_refl_dry_dir",nsoil,nband,soil%dat_refl_dry_dir)
  call get_parameter_data(grpid,"dat_refl_dry_dif",nsoil,nband,soil%dat_refl_dry_dif)
  call get_parameter_data(grpid,"dat_refl_sat_dir",nsoil,nband,soil%dat_refl_sat_dir)
  call get_parameter_data(grpid,"dat_refl_sat_dif",nsoil,nband,soil%dat_refl_sat_dif)
  call get_parameter_data(grpid,"dat_emis_dry",nsoil,soil%dat_emis_dry)
  call get_parameter_data(grpid,"dat_emis_sat",nsoil,soil%dat_emis_sat)
  call get_parameter_data(grpid,"dat_z0_momentum",nsoil,soil%dat_z0_momentum)
  call get_parameter_data(grpid,"dat_tf_depr",nsoil,soil%dat_tf_depr)
  call get_parameter_data(grpid,"rsa_exp_global",nsoil,soil%rsa_exp_global)
  call get_parameter_data(grpid,"gw_res_time",nsoil,soil%gw_res_time)
  call get_parameter_data(grpid,"gw_hillslope_length",nsoil,soil%gw_hillslope_length)
  call get_parameter_data(grpid,"gw_scale_length",nsoil,soil%gw_scale_length)
  call get_parameter_data(grpid,"gw_hillslope_zeta_bar",nsoil,soil%gw_hillslope_zeta_bar)
  call get_parameter_data(grpid,"gw_hillslope_relief",nsoil,soil%gw_hillslope_relief)
  call get_parameter_data(grpid,"gw_scale_relief",nsoil,soil%gw_scale_relief)
  call get_parameter_data(grpid,"gw_soil_e_depth",nsoil,soil%gw_soil_e_depth)
  call get_parameter_data(grpid,"gw_scale_soil_depth",nsoil,soil%gw_scale_soil_depth)
  !call get_parameter_data(grpid,"gw_hillslope_a",nsoil,soil%gw_hillslope_a)
  !call get_parameter_data(grpid,"gw_hillslope_n",nsoil,soil%gw_hillslope_n)
  call get_parameter_data(grpid,"gw_perm",nsoil,soil%gw_perm)
  call get_parameter_data(grpid,"gw_scale_perm",nsoil,soil%gw_scale_perm)
  call get_parameter_data(grpid,"gw_res_time",nsoil,soil%gw_res_time)
  call get_parameter_data(grpid,"microtopo",nsoil,soil%microtopo)
  call get_parameter_data(grpid,"tile_hlsp_length",nsoil,soil%tile_hlsp_length)
  call get_parameter_data(grpid,"tile_hlsp_slope",nsoil,soil%tile_hlsp_slope)
  call get_parameter_data(grpid,"tile_hlsp_elev",nsoil,soil%tile_hlsp_elev)
  call get_parameter_data(grpid,"tile_hlsp_hpos",nsoil,soil%tile_hlsp_hpos)
  call get_parameter_data(grpid,"tile_hlsp_width",nsoil,soil%tile_hlsp_width)
  call get_parameter_data(grpid,"hidx_k",nsoil,soil%hidx_k)
  call get_parameter_data(grpid,"hidx_j",nsoil,soil%hidx_j)
  call get_parameter_data(grpid,"vegn",nsoil,soil%vegn)

end subroutine retrieve_soil_parameters

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

subroutine get_parameter_data_2d_integer(grpid,var,nx,ny,tmp)

 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx,ny
 integer,dimension(:,:),pointer :: tmp
 integer :: itile,varid,status
 allocate(tmp(nx,ny))

 status = nf90_inq_varid(grpid,var,varid)
 status = nf90_get_var(grpid,varid,tmp)

end subroutine

subroutine get_parameter_data_2d_real(grpid,var,nx,ny,tmp)

 character(len=*),intent(in) :: var
 integer,intent(in) :: grpid,nx,ny
 real,dimension(:,:),pointer :: tmp
 integer :: itile,varid,status
 allocate(tmp(nx,ny))

 status = nf90_inq_varid(grpid,var,varid)
 status = nf90_get_var(grpid,varid,tmp)

end subroutine

end module predefined_tiles_mod
