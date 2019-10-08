module  standalone_stochy_module

use machine
implicit none
public
type GFS_control_type
   integer  :: levs,me
   integer,allocatable  :: blksz(:)        !< for explicit data blocking: block sizes of all blocks
   real(kind=kind_phys) :: dtp             !< physics timestep in seconds
   real(kind=kind_phys) :: phour           !< previous forecast hour
   integer              :: kdt             !< current forecast iteration
   logical  :: do_sppt,do_shum,do_skeb,do_sfcperts,use_zmtnblck
   integer  ::  skeb_npass,nsfcpert
   character(len=65) :: fn_nml                   !< namelist filename
   character(len=256),allocatable :: input_nml_file(:) !< character string containing full namelist
   real(kind=kind_phys) ::  pertz0(5),pertzt(5),pertshc(5),pertlai(5),pertvegf(5),pertalb(5)
end type GFS_control_type

type GFS_init_type
   integer :: nlunit
   real(kind=kind_phys),allocatable :: ak(:),bk(:)
end type GFS_init_type

type GFS_grid_type
    real (kind=kind_phys),allocatable :: xlat   (:)    !< grid latitude in radians, default to pi/2 ->
    real (kind=kind_phys),allocatable :: xlon   (:)    !< grid longitude in radians, default to pi/2 ->
end type GFS_grid_type

type GFS_coupling_type
    real (kind=kind_phys),allocatable :: shum_wts  (:,:)
    real (kind=kind_phys),allocatable :: sppt_wts  (:,:)
    real (kind=kind_phys),allocatable :: skebu_wts (:,:)
    real (kind=kind_phys),allocatable :: skebv_wts (:,:)
    real (kind=kind_phys),allocatable :: sfc_wts   (:,:)
    integer              :: nsfcpert=6                             !< number of sfc perturbations
end type GFS_coupling_type
end module standalone_stochy_module

