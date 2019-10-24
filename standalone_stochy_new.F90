program  standalone_stochy_new

use standalone_stochy_module
use stochastic_physics_new, only : init_stochastic_physics_atm,run_stochastic_physics_atm
use stochy_namelist_def_new, only : do_sppt,do_shum,do_skeb,do_sfcperts

use fv_mp_mod,           only: mp_start, domain_decomp
use atmosphere_stub_mod, only: Atm,atmosphere_init_stub
use fv_arrays_mod,       only: fv_atmos_type
use fv_control_mod,      only: setup_pointers
!use mpp_domains
use mpp_mod,             only: mpp_set_current_pelist,mpp_get_current_pelist,mpp_init,mpp_pe,mpp_npes ,mpp_declare_pelist
use mpp_domains_mod,     only: mpp_broadcast_domain,MPP_DOMAIN_TIME,mpp_domains_init ,mpp_domains_set_stack_size
use fms_mod,             only:  fms_init
!use time_manager_mod,    only: time_type
use xgrid_mod,           only: grid_box_type
use netcdf

implicit none
type(GFS_control_type)  :: Model
type(GFS_init_type)     :: Init_parm
integer, parameter      :: nlevs=64
integer                 :: ntasks,fid
integer                 :: nthreads,omp_get_num_threads
integer                 :: ncid,xt_dim_id,yt_dim_id,time_dim_id,xt_var_id,yt_var_id,time_var_id,sppt_id
character*1             :: strid
type(GFS_grid_type),allocatable     :: Grid(:)
type(GFS_coupling_type),allocatable :: Coupling(:)
include 'mpif.h'
include 'netcdf.inc'
real :: ak(nlevs+1),bk(nlevs+1)
real(kind=4) :: ts,undef

data ak(:) /0.000, 0.000, 0.575, 5.741, 21.516, 55.712, 116.899, 214.015, 356.223, 552.720, 812.489, &
   1143.988, 1554.789, 2051.150, 2637.553, 3316.217, 4086.614, 4945.029, 5884.206, 6893.117,    &
   7956.908, 9057.051, 10171.712, 11276.348, 12344.490, 13348.671, 14261.435, 15056.342,        &
   15708.893, 16197.315, 16503.145, 16611.604, 16511.736, 16197.967, 15683.489, 14993.074,      &
   14154.316, 13197.065, 12152.937, 11054.853, 9936.614, 8832.537, 7777.150, 6804.874, 5937.050,&
   5167.146, 4485.493, 3883.052, 3351.460, 2883.038, 2470.788, 2108.366, 1790.051, 1510.711,    &
   1265.752, 1051.080, 863.058, 698.457, 554.424, 428.434, 318.266, 221.958, 137.790, 64.247,0.0 /
data bk(:) /1.00000000, 0.99467117, 0.98862660, 0.98174226, 0.97386760, 0.96482760, 0.95443410, 0.94249105, &
  0.92879730, 0.91315103, 0.89535499, 0.87522358, 0.85259068, 0.82731885, 0.79930973, 0.76851469, &
  0.73494524, 0.69868290, 0.65988702, 0.61879963, 0.57574666, 0.53113484, 0.48544332, 0.43921080, &
  0.39301825, 0.34746850, 0.30316412, 0.26068544, 0.22057019, 0.18329623, 0.14926878, 0.11881219, &
  0.09216691, 0.06947458, 0.05064684, 0.03544162, 0.02355588, 0.01463712, 0.00829402, 0.00410671, &
  0.00163591, 0.00043106, 0.00003697, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, &
  0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, &
  0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, &
  0.00000000 /
integer     :: cres,blksz,nblks,ierr,my_id,i,j,nx2,ny2,nx,ny,id
integer,target :: npx,npy
integer     :: ng,layout(2),io_layout(2),commID,grid_type,ntiles
integer :: halo_update_type = 1
real        :: dx,dy,pi
logical,target :: nested
integer  :: pe,npes,stackmax=4000000

real(kind=4),allocatable,dimension(:,:) :: workg
real(kind=8),pointer    ,dimension(:,:) :: area
type(grid_box_type)           :: grid_box
!type(time_type)               :: Time               ! current time
!type(time_type)               :: Time_step          ! atmospheric time step.
!type(time_type)               :: Time_init          ! reference time.

! define stuff
ng=3
pi=3.14159265359
undef=9.99e+20

call fms_init()
call mpp_init()
call fms_init
my_id=mpp_pe()

call atmosphere_init_stub (grid_box, area)
isd=Atm(1)%bd%isd
ied=Atm(1)%bd%ied
jsd=Atm(1)%bd%jsd
jed=Atm(1)%bd%jed
isc=Atm(1)%bd%isc
iec=Atm(1)%bd%iec
jsc=Atm(1)%bd%jsc
jec=Atm(1)%bd%jec
nx=Atm(1)%npx-1
ny=Atm(1)%npy-1
allocate(workg(nx,ny))
nblks=ny
blksz=nx
nthreads = omp_get_num_threads()
Model%me=my_id
Model%phour=0
Model%kdt=1
Model%dtp=900
Model%fn_nml='input.nml'
Model%levs=nlevs
allocate(Init_parm%blksz(nblks))
Init_parm%blksz(:)=blksz
! setup GFS_init parameters
allocate(Init_parm%ak(nlevs+1))
allocate(Init_parm%bk(nlevs+1))
Init_parm%ak=ak
Init_parm%bk=bk
Init_parm%nlunit=21

!define model grid
Model%nx=nx
Model%ny=ny
dx=360.0/Model%nx
dy=180.0/Model%ny
allocate(Init_parm%xlon(Model%nx,Model%ny))
allocate(Init_parm%xlat(Model%nx,Model%ny))
Init_parm%xlon(:,:)=Atm(1)%gridstruct%agrid(:,:,1)
Init_parm%xlat(:,:)=Atm(1)%gridstruct%agrid(:,:,2)

!setup GFS_coupling
allocate(Coupling(nblks))
write(strid,'(I1.1)') my_id+1
fid=30+my_id
ierr=nf90_create('workg.tile'//strid//'.nc',cmode=NF90_CLOBBER,ncid=ncid)
ierr=NF90_DEF_DIM(ncid,"grid_xt",nx,xt_dim_id)
ierr=NF90_DEF_DIM(ncid,"grid_yt",ny,yt_dim_id)
ierr=NF90_DEF_DIM(ncid,"time",NF90_UNLIMITED,time_dim_id)
  !> - Define the dimension variables.
ierr=NF90_DEF_VAR(ncid,"grid_xt",NF90_FLOAT,(/ xt_dim_id /), xt_var_id)
ierr=NF90_PUT_ATT(ncid,xt_var_id,"long_name","T-cell longitude")
ierr=NF90_PUT_ATT(ncid,xt_var_id,"cartesian_axis","X")
ierr=NF90_PUT_ATT(ncid,xt_var_id,"units","degrees_E")
ierr=NF90_DEF_VAR(ncid,"grid_yt",NF90_FLOAT,(/ yt_dim_id /), yt_var_id)
ierr=NF90_PUT_ATT(ncid,yt_var_id,"long_name","T-cell latitude")
ierr=NF90_PUT_ATT(ncid,yt_var_id,"cartesian_axis","Y")
ierr=NF90_PUT_ATT(ncid,yt_var_id,"units","degrees_N")
ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
ierr=NF90_PUT_ATT(ncid,time_var_id,"long_name","time")
ierr=NF90_PUT_ATT(ncid,time_var_id,"units","hours since 2014-08-01 00:00:00")
ierr=NF90_PUT_ATT(ncid,time_var_id,"cartesian_axis","T")
ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar_type","JULIAN")
ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar","JULIAN")
ierr=NF90_DEF_VAR(ncid,"sppt_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), sppt_id)
ierr=NF90_PUT_ATT(ncid,sppt_id,"long_name","random pattern")
ierr=NF90_PUT_ATT(ncid,sppt_id,"units","None")
ierr=NF90_PUT_ATT(ncid,sppt_id,"missing_value",undef)
ierr=NF90_PUT_ATT(ncid,sppt_id,"_FillValue",undef)
ierr=NF90_PUT_ATT(ncid,sppt_id,"cell_methods","time: point")
ierr=NF90_ENDDEF(ncid)

call init_stochastic_physics_atm(Model, Init_parm, ntasks, nthreads)
do i=1,nblks
   if (do_sppt) allocate(Coupling(i)%sppt_wts(blksz,nlevs))
   if (do_shum) allocate(Coupling(i)%shum_wts(blksz,nlevs))
   if (do_skeb) allocate(Coupling(i)%skebu_wts(blksz,nlevs))
   if (do_skeb) allocate(Coupling(i)%skebv_wts(blksz,nlevs))
enddo
do i=1,10
   ts=i/4.0
   call run_stochastic_physics_atm(Model, Grid, Coupling, nthreads)
   do j=1,ny
      workg(:,j)=Coupling(j)%sppt_wts(:,30)   
   enddo
   !write(fid) workg
   ierr=NF90_PUT_VAR(ncid,sppt_id,workg,(/1,1,i/))
   ierr=NF90_PUT_VAR(ncid,time_var_id,ts,(/i/))
   if (my_id.EQ.0) write(6,fmt='(a,i5,4f6.3)') 'sppt_wts=',i,Coupling(1)%sppt_wts(1:4,30)
enddo
!close(fid)
ierr=NF90_CLOSE(ncid)
end
