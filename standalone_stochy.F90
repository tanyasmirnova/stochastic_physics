program  standalone_stochy

use standalone_stochy_module
use stochastic_physics, only : init_stochastic_physics,run_stochastic_physics
use stochy_namelist_def, only : do_sppt,do_shum,do_skeb,do_sfcperts
implicit none
type(GFS_control_type)  :: Model
type(GFS_init_type)     :: Init_parm
integer, parameter      :: nlevs=64
integer                 :: ntasks
integer                 :: nthreads,omp_get_num_threads
type(GFS_grid_type),allocatable     :: Grid(:)
type(GFS_coupling_type),allocatable :: Coupling(:)
include 'mpif.h'
real :: ak(nlevs+1),bk(nlevs+1)
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
integer                  :: cres,blksz,nblks,ierr,my_id,i

! define stuff
print*,'1'
call MPI_INIT(ierr)
print*,'2',ierr
call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
print*,'3',ierr
call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
print*,'4',ierr
!cres=24
!blksz=32
!nblks=cres*cres*6/ntasks
nblks=1
blksz=4 
!nblks=cres*cres*6/ntasks
nthreads = omp_get_num_threads()
print*,'nthreads=',nthreads
! setup GFS_control parameters
Model%me=my_id
Model%phour=0
Model%kdt=1
Model%dtp=900
Model%fn_nml='input.nml'
Model%levs=nlevs
allocate(Model%blksz(nblks))
Model%blksz(:)=blksz
! setup GFS_init parameters
allocate(Init_parm%ak(nlevs+1))
allocate(Init_parm%bk(nlevs+1))
Init_parm%ak=ak
Init_parm%bk=bk
Init_parm%nlunit=21

!setup GFS_grid parameters
allocate(Grid(nblks))
do i=1,nblks
   allocate(Grid(i)%xlon(blksz))
   allocate(Grid(i)%xlat(blksz))
enddo
!define model grid
Grid(1)%xlon(1)=0.0
Grid(1)%xlon(2)=90.0
Grid(1)%xlon(3)=180.0
Grid(1)%xlon(4)=270.0
Grid(1)%xlat(1)=-40.0
Grid(1)%xlat(2)=0.0
Grid(1)%xlat(3)=40.0
Grid(1)%xlat(4)=80.0

Grid(1)%xlon=Grid(1)%xlon*3.1415928/180.0
Grid(1)%xlat=Grid(1)%xlat*3.1415928/180.0

!setup GFS_coupling
allocate(Coupling(nblks))
print*,'calling init'
print*,'my_id=',my_id
if (my_id.EQ.0) open(30,file='workg.out',form='unformatted')
call init_stochastic_physics(Model, Init_parm, ntasks, nthreads)
do i=1,nblks
   if (do_sppt) allocate(Coupling(i)%sppt_wts(blksz,nlevs))
   if (do_shum) allocate(Coupling(i)%shum_wts(blksz,nlevs))
   if (do_skeb) allocate(Coupling(i)%skebu_wts(blksz,nlevs))
   if (do_skeb) allocate(Coupling(i)%skebv_wts(blksz,nlevs))
enddo
do i=1,1000
   call run_stochastic_physics(Model, Grid, Coupling, nthreads)
   if (my_id.EQ.0) write(6,fmt='(a,4f6.3)') 'sppt_wts=',Coupling(1)%sppt_wts(:,40)
   if (my_id.EQ.0) print*,'time-step=',i
   Model%kdt=Model%kdt+1
enddo
call mpi_barrier(mpi_comm_world,ierr)
if(my_id.EQ.0)close(30)
end
