module stochastic_physics_new
#ifdef STOCHY_UNIT_TEST
use atmosphere_stub_mod, only: atmosphere_smooth_noise,atmosphere_return_winds
#else
use fv_mp_mod, only : isc,iec,jsc,jec,isd,ied,jsd,jed
use atmosphere_mod, only: atmosphere_smooth_noise,atmosphere_return_winds
#endif
use stochy_patterngenerator_new
use stochy_namelist_def_new
use mpp_mod
use MPI
#ifdef STOCHY_UNIT_TEST
 use standalone_stochy_module,   only: GFS_control_type, GFS_init_type,GFS_coupling_type,GFS_grid_type,isc,iec,jsc,jec,isd,ied,jsd,jed
#else
use GFS_typedefs,       only: GFS_control_type, GFS_init_type,GFS_control_type,GFS_grid_type
#endif
use mersenne_twister, only: random_gauss
implicit none
private

public :: init_stochastic_physics_atm, run_stochastic_physics_atm

contains

subroutine init_stochastic_physics_atm(Model,Init_parm, ntasks, nthreads)

type(GFS_control_type),   intent(inout) :: Model
type(GFS_init_type),      intent(in) :: Init_parm
integer,                  intent(in)    :: ntasks
integer,                  intent(in)    :: nthreads
real :: PRSI(Model%levs),PRSL(Model%levs),dx,dtp
real, allocatable :: skeb_vloc(:),rnda_vloc(:)
real(kind=kind_phys),allocatable :: tnoise(:)
integer :: i,j,ix,k,kflip,latghf,nodes,blk,k2,me
character*2::proc

! replace
dtp=Model%dtp
me=Model%me
nodes=ntasks
gis_stochy_atm%me=me
gis_stochy_atm%nodes=nodes
gis_stochy_atm%is=isc
gis_stochy_atm%ie=iec
gis_stochy_atm%js=jsc
gis_stochy_atm%je=jec
gis_stochy_atm%isd=isd
gis_stochy_atm%ied=ied
gis_stochy_atm%jsd=jsd
gis_stochy_atm%jed=jed

!call init_stochdata(dtp,Model%input_nml_file,Model%fn_nml,Init_parm%nlunit,Init_parm%xlon,Init_parm%xlat,gis_stochy_atm)
call init_stochdata(dtp,Init_parm%xlon,Init_parm%xlat,gis_stochy_atm)
gis_stochy_atm%nblks = size(Init_parm%blksz)
allocate (gis_stochy_atm%blksz(gis_stochy_atm%nblks))
gis_stochy_atm%blksz(:) = Init_parm%blksz(:)
! check to see decomposition
Model%do_sppt=do_sppt
Model%use_zmtnblck=use_zmtnblck
Model%do_shum=do_shum
Model%do_skeb=do_skeb
Model%do_rnda=do_rnda
Model%do_sfcperts=do_sfcperts             ! mg, sfc-perts
Model%nsfcpert=nsfcpert         ! mg, sfc-perts
Model%pertz0=pertz0         ! mg, sfc-perts
Model%pertzt=pertzt         ! mg, sfc-perts
Model%pertshc=pertshc         ! mg, sfc-perts
Model%pertlai=pertlai         ! mg, sfc-perts
Model%pertalb=pertalb         ! mg, sfc-perts
Model%pertvegf=pertvegf         ! mg, sfc-perts
if ( (.NOT. do_sppt) .AND. (.NOT. do_skeb) .AND. (.NOT. do_shum) .AND. (.NOT. do_rnda)  .AND. (.NOT. do_sfcperts) ) return
allocate(sl(Model%levs))
do k=1,Model%levs
   sl(k)= 0.5*(Init_parm%ak(k)/101300.+Init_parm%bk(k)+Init_parm%ak(k+1)/101300.0+Init_parm%bk(k+1)) ! si are now sigmas
!   if(me==0)print*,'sl(k)',k,sl(k),Init_parm%ak(k),Init_parm%bk(k)
enddo
if (do_sppt) then
   allocate(vfact_sppt(Model%levs))
   do k=1,Model%levs
      if (sl(k) .lt. sppt_sigtop1 .and. sl(k) .gt. sppt_sigtop2) then
         vfact_sppt(k) = (sl(k)-sppt_sigtop2)/(sppt_sigtop1-sppt_sigtop2)
      else if (sl(k) .lt. sppt_sigtop2) then
          vfact_sppt(k) = 0.0
      else
          vfact_sppt(k) = 1.0
      endif
   enddo
   if (sppt_sfclimit) then
       vfact_sppt(2)=vfact_sppt(3)*0.5
       vfact_sppt(1)=0.0
   endif
   if (me==0) then
      do k=1,Model%levs
         print *,'sppt vert profile',k,sl(k),vfact_sppt(k)
      enddo
   endif
   !call advance_pattern(rpattern_sppt(1),gis_stochy_atm)
endif
if (do_skeb) then
   allocate(vfact_skeb(Model%levs))
   allocate(skeb_vloc(gis_stochy_atm%skeblevs)) ! local
   allocate(skeb_vwts(Model%levs,2)) ! save for later
   allocate(skeb_vpts(Model%levs,2)) ! save for later
   do k=1,Model%levs
      if (sl(k) .lt. skeb_sigtop1 .and. sl(k) .gt. skeb_sigtop2) then
         vfact_skeb(k) = (sl(k)-skeb_sigtop2)/(skeb_sigtop1-skeb_sigtop2)
      else if (sl(k) .lt. skeb_sigtop2) then
          vfact_skeb(k) = 0.0
      else
          vfact_skeb(k) = 1.0
      endif
      if (me==0)  print *,'skeb vert profile',k,sl(k),vfact_skeb(k) 
   enddo
! calculate vertical interpolation weights
   do k=1,gis_stochy_atm%skeblevs
      skeb_vloc(k)=sl(1)-real(k-1)/real(gis_stochy_atm%skeblevs-1.0)*(sl(1)-sl(Model%levs))
   enddo
! surface
skeb_vwts(1,2)=0
skeb_vpts(1,1)=1
! top
skeb_vwts(Model%levs,2)=1
skeb_vpts(Model%levs,1)=gis_stochy_atm%skeblevs-1
! internal
do k=2,Model%levs-1
   do k2=1,gis_stochy_atm%skeblevs-1
      IF (sl(k) .LE. skeb_vloc(k2) .AND. sl(k) .GT. skeb_vloc(k2+1)) THEN
        skeb_vpts(k,1)=k2
        skeb_vwts(k,2)=(skeb_vloc(k2)-sl(k))/(skeb_vloc(k2)-skeb_vloc(k2+1))
      ENDIF
   enddo 
enddo  
deallocate(skeb_vloc)
if (me==0) then
do k=1,Model%levs
   print*,'skeb vpts ',skeb_vpts(k,1),skeb_vwts(k,2)
enddo
endif
skeb_vwts(:,1)=1.0-skeb_vwts(:,2)
skeb_vpts(:,2)=skeb_vpts(:,1)+1.0
endif

if (do_rnda) then
   allocate(vfact_rnda(Model%levs))
   do k=1,Model%levs
      if (sl(k) .lt. rnda_sigtop1 .and. sl(k) .gt. rnda_sigtop2) then
         vfact_rnda(k) = (sl(k)-rnda_sigtop2)/(rnda_sigtop1-rnda_sigtop2)
      else if (sl(k) .lt. rnda_sigtop2) then
          vfact_rnda(k) = 0.0
      else
          vfact_rnda(k) = 1.0
      endif
      if (me==0)  print *,'rnda vert profile',k,sl(k),vfact_rnda(k) 
   enddo
   allocate(rndau_save(isc:iec+1,jsc:jec))
   allocate(rndav_save(isc:iec  ,jsc:jec+1))
   allocate(tnoise(isc:iec))
   do j=jsc,jec
      call random_gauss(tnoise,rpattern_rnda(1)%rstate2)
      rpattern_rnda(1)%wnoise(isc:ied,j,1)=tnoise
   enddo
   deallocate(tnoise)
   call atmosphere_smooth_noise (rpattern_rnda(1)%wnoise,rpattern_rnda(1)%nsmth,ns_type,0) !  do smoothing for RNDA
   rpattern_rnda(1)%rnoise(isc:iec,jsc:jec)= rpattern_rnda(1)%wnoise(isc:iec,jsc:jec,1)*rpattern_rnda(1)%stdev(1)
   call atmosphere_return_winds(rpattern_rnda(1)%rnoise(:,:),rndau_save(:,:),rndav_save(:,:),1,km=Model%levs,vwts=vfact_rnda)
endif

if (do_shum) then
   allocate(vfact_shum(Model%levs))
   do k=1,Model%levs
      vfact_shum(k) = exp((sl(k)-1.)/shum_sigefold)
      if (sl(k).LT. 2*shum_sigefold) then
         vfact_shum(k)=0.0
      endif
      if (me==0)  print *,'shum vert profile',k,sl(k),vfact_shum(k)
   enddo
endif

!print *,'done with init_stochastic_physics'

end subroutine init_stochastic_physics_atm

subroutine run_stochastic_physics_atm(Model,Grid,Coupling, nthreads)

type(GFS_control_type),   intent(in) :: Model
type(GFS_coupling_type),  intent(inout) :: Coupling(gis_stochy_atm%nblks)
type(GFS_grid_type),  intent(inout) :: Grid(gis_stochy_atm%nblks)
integer,                  intent(in)    :: nthreads
real, allocatable :: tmp_wts(:,:)
real(kind=kind_phys),allocatable :: tnoise(:)
         
!D-grid
integer :: k
integer j,ierr,i
integer :: blk,nb,ix
character*120 :: sfile
character*6   :: STRFH
if ( (.NOT. do_sppt) .AND. (.NOT. do_skeb) .AND. (.NOT. do_shum) ) return
! check to see if it is time to write out random patterns
!if (Model%phour .EQ. fhstoch) then
!   write(STRFH,FMT='(I6.6)') nint(Model%fhour)
!   sfile='stoch_out.F'//trim(STRFH)
!   call dump_patterns(sfile)
!endif

allocate(tmp_wts(gis_stochy_atm%nlons,gis_stochy_atm%nlats))
allocate(tnoise(isc:iec))
if (do_sppt) then
   if (stoch_ptype.NE.1) then
      call make_gridded(rpattern_sppt(1),tmp_wts,gis_stochy_atm)
      print*,'after make gridded',maxval(tmp_wts),minval(tmp_wts)
      call advance_pattern(rpattern_sppt(1),gis_stochy_atm)
   endif
   if (stoch_ptype.NE.0) then
      if (first_call .OR. rpattern_sppt(1)%tau(1) .GE. 0) then
         do j=jsc,jec
            call random_gauss(tnoise,rpattern_sppt(1)%rstate2)
            rpattern_sppt(1)%wnoise(isc:iec,j,1)=tnoise
         enddo
         !print*,'1st',maxval(rpattern_sppt(1)%wnoise),minval(rpattern_sppt(1)%wnoise)
      else
          rpattern_sppt(1)%wnoise(isc:iec,jsc:jec,1)=rpattern_sppt(1)%rnoise
      endif
! remove for real runs
!      if(.not.first_call) call atmosphere_smooth_noise (rpattern_sppt(1)%wnoise,rpattern_sppt(1)%nsmth,ns_type,1) !  do smoothing for sppt
      call atmosphere_smooth_noise (rpattern_sppt(1)%wnoise,rpattern_sppt(1)%nsmth,ns_type,1) !  do smoothing for sppt
   endif
   nb = 1
   ix = 0
   do j = 1,Model%ny
      do i = 1,Model%nx
         ix=ix+1
         if (ix > gis_stochy_atm%blksz(nb)) then
            nb = nb + 1
            ix = 1
         endif
         if (stoch_ptype.NE.0) then
            if (first_call.OR. rpattern_sppt(1)%tau(1) .LT. 0) then
               rpattern_sppt(1)%rnoise(i+isc-1,j+jsc-1)= rpattern_sppt(1)%wnoise(i+isc-1,j+jsc-1,1)*rpattern_sppt(1)%stdev(1)
            else
               rpattern_sppt(1)%rnoise(i+isc-1,j+jsc-1)=rpattern_sppt(1)%phi(1)*rpattern_sppt(1)%rnoise(i+isc-1,j+jsc-1)+  & ! AR(1)
               sqrt(1-rpattern_sppt(1)%phi(1)**2)*rpattern_sppt(1)%stdev(1)*rpattern_sppt(1)%wnoise(i+isc-1,j+jsc-1,1)
            endif
         endif
         do k=1,Model%levs
            if (stoch_ptype.EQ.0) then
               Coupling(nb)%sppt_wts(ix,k)=tmp_wts(i,j)*vfact_sppt(k)
            endif
            if (stoch_ptype.EQ.1) then
               Coupling(nb)%sppt_wts(ix,k)= rpattern_sppt(1)%rnoise(i+isc-1,j+jsc-1)*vfact_sppt(k)
            endif
            if (stoch_ptype.EQ.2) then
               Coupling(nb)%sppt_wts(ix,k)= (rpattern_sppt(1)%rnoise(i+isc-1,j+jsc-1) + tmp_wts(i,j))*vfact_sppt(k)
            endif
         enddo
      enddo
   enddo
   do i=1,gis_stochy_atm%nblks
      Coupling(i)%sppt_pattern(:)= Coupling(i)%sppt_wts(:,30)
   enddo
   if (sppt_logit) then
      do i=1,gis_stochy_atm%nblks
         Coupling(i)%sppt_wts(:,:) = (2./(1.+exp(Coupling(i)%sppt_wts(:,:))))
      enddo
   else
      do i=1,gis_stochy_atm%nblks
          Coupling(i)%sppt_wts(:,:)= Coupling(i)%sppt_wts(:,:)+1.0
      enddo
   endif
endif
if (do_shum) then
   if (stoch_ptype.NE.1) then
      call make_gridded(rpattern_shum(1),tmp_wts,gis_stochy_atm)
      call advance_pattern(rpattern_shum(1),gis_stochy_atm)
   endif
   if (stoch_ptype.NE.0) then
      if (first_call .OR. rpattern_shum(1)%tau(1) .GE. 0) then
         do j=jsc,jec
            call random_gauss(tnoise,rpattern_shum(1)%rstate2)
            rpattern_shum(1)%wnoise(isc:iec,j,1)=tnoise
         enddo
      else
          rpattern_shum(1)%wnoise(:,:,1)=rpattern_shum(1)%rnoise
      endif
      call atmosphere_smooth_noise (rpattern_shum(1)%wnoise,rpattern_shum(1)%nsmth,ns_type,1) !  do smoothing for SHUM
   endif
   nb = 1
   ix = 0
   if (first_call.OR. rpattern_shum(1)%tau(1) .LT. 0) then
      do j = 1,Model%ny
         do i = 1,Model%nx
            ix=ix+1
            if (ix > gis_stochy_atm%blksz(nb)) then
               nb = nb + 1
               ix = 1
            endif
            do k=1,Model%levs
               if (stoch_ptype.EQ.0) Coupling(nb)%shum_wts(ix,k)=tmp_wts(i,j)*vfact_shum(k)
               if (stoch_ptype.EQ.1) Coupling(nb)%shum_wts(ix,k)=rpattern_shum(1)%stdev(1)*rpattern_shum(1)%wnoise(i+isc-1,j+jsc-1,1)*vfact_shum(k)
            enddo
         enddo
      enddo
   else
      do j = 1,Model%ny
         do i = 1,Model%nx
            ix=ix+1
            if (ix > gis_stochy_atm%blksz(nb)) then
               nb = nb + 1
               ix = 1
            endif
            do k=1,Model%levs
               if (stoch_ptype.EQ.0) Coupling(nb)%shum_wts(ix,k)=tmp_wts(i,j)*vfact_shum(k)
               if (stoch_ptype.EQ.1) Coupling(nb)%shum_wts(ix,k)=rpattern_shum(1)%phi(1)*Coupling(nb)%shum_wts(ix,k)+  &
                                      sqrt(1-rpattern_shum(1)%phi(1)**2)*rpattern_shum(1)%stdev(1)*rpattern_shum(1)%wnoise(i+isc-1,j+jsc-1,1)*vfact_shum(k)
            enddo
         enddo
      enddo
   endif
endif
if (do_skeb) then
   if (mod(Model%kdt,nsskeb).EQ.0 .OR. first_call) then
      if (first_call) then
         allocate(skebu_save(isc:iec,jsc:jec,gis_stochy_atm%skeblevs))
         allocate(skebv_save(isc:iec,jsc:jec,gis_stochy_atm%skeblevs))
         !skebu_save(isc:iec,jsc:jec,gis_stochy_atm%skeblevs)=0.0
         !skebv_save(isc:iec,jsc:jec,gis_stochy_atm%skeblevs)=0.0
         do k=1,gis_stochy_atm%skeblevs-1
            if (stoch_ptype.NE.1) then
               call make_gridded_vecta(rpattern_skeb(1),skebu_save(:,:,k),skebv_save(:,:,k),gis_stochy_atm)
               call advance_pattern(rpattern_skeb(1),gis_stochy_atm)
            endif
            if (stoch_ptype.NE.0) then
               if (first_call .OR. rpattern_skeb(1)%tau(1) .GE. 0) then
                  do j=jsc,jec
                     call random_gauss(tnoise,rpattern_skeb(1)%rstate2)
                     rpattern_skeb(1)%wnoise(isc:iec,j,1)=tnoise
                  enddo
               else
                  rpattern_skeb(1)%wnoise(:,:,1)=rpattern_skeb(1)%rnoise
               endif
               call atmosphere_smooth_noise (rpattern_skeb(1)%wnoise,rpattern_skeb(1)%nsmth,ns_type,0) !  do smoothing for SKEB
               if (k.EQ.1.OR. rpattern_skeb(1)%tau(1) .LT. 0) then
                  rpattern_skeb(1)%rnoise(isc:iec,jsc:jec)= rpattern_skeb(1)%wnoise(isc:iec,jsc:jec,1)*rpattern_skeb(1)%stdev(1)
               else
                  rpattern_skeb(1)%rnoise(:,:)=rpattern_skeb(1)%phi(1)*rpattern_skeb(1)%rnoise(:,:)+  &
                  sqrt(1-rpattern_skeb(1)%phi(1)**2)*rpattern_skeb(1)%stdev(1)*rpattern_skeb(1)%wnoise(:,:,1)
               endif
               call atmosphere_return_winds(rpattern_skeb(1)%rnoise(:,:),skebu_save(:,:,k),skebv_save(:,:,k),0)
            endif
         enddo
      else
      ! shift patterns
         do k=1,gis_stochy_atm%skeblevs-1
            skebu_save(:,:,k)=skebu_save(:,:,k+1)
            skebv_save(:,:,k)=skebv_save(:,:,k+1)
         enddo
      endif
      if (stoch_ptype.NE.1) then
         call make_gridded_vecta(rpattern_skeb(1),skebu_save(:,:,gis_stochy_atm%skeblevs),skebv_save(:,:,gis_stochy_atm%skeblevs),gis_stochy_atm)
         call advance_pattern(rpattern_skeb(1),gis_stochy_atm)
      endif
      if (stoch_ptype.NE.0) then
         !print*,'calling random_gauss'
         if (rpattern_skeb(1)%tau(1) .GE. 0) then
            do j=jsc,jec
               call random_gauss(tnoise,rpattern_skeb(1)%rstate2)
               rpattern_skeb(1)%wnoise(isc:iec,j,1)=tnoise
            enddo
         else
             rpattern_skeb(1)%wnoise(:,:,1)=rpattern_skeb(1)%rnoise
         endif
         call atmosphere_smooth_noise (rpattern_skeb(1)%wnoise,rpattern_skeb(1)%nsmth,ns_type,0) !  do smoothing for SKEB
         if (rpattern_skeb(1)%tau(1) .LT. 0) then
             rpattern_skeb(1)%rnoise(isc:iec,jsc:jec)= rpattern_skeb(1)%wnoise(isc:iec,jsc:jec,1)*rpattern_skeb(1)%stdev(1)
         else
            rpattern_skeb(1)%rnoise(isc:iec,jsc:jec)=rpattern_skeb(1)%phi(1)*rpattern_skeb(1)%rnoise(isc:iec,jsc:jec)+  &
            sqrt(1-rpattern_skeb(1)%phi(1)**2)*rpattern_skeb(1)%stdev(1)*rpattern_skeb(1)%wnoise(isc:iec,jsc:jec,1)
         endif
         call atmosphere_return_winds(rpattern_skeb(1)%rnoise(:,:),skebu_save(:,:,gis_stochy_atm%skeblevs),skebv_save(:,:,gis_stochy_atm%skeblevs),0)
         
      endif
 
      !print*,'vertical interpolation',Model%levs
      do k=1,Model%levs
         nb = 1
         ix = 0
         do j = 1,Model%ny
            do i = 1,Model%nx
               ix=ix+1
               if (ix > gis_stochy_atm%blksz(nb)) then
                  nb = nb + 1
                  ix = 1
               endif
               Coupling(nb)%skebu_wts(ix,k)=(skeb_vwts(k,1)*skebu_save(i+isc-1,j+jsc-1,skeb_vpts(k,1))+skeb_vwts(k,2)*skebu_save(i+isc-1,j+jsc-1,skeb_vpts(k,2)))*vfact_skeb(k)
               Coupling(nb)%skebv_wts(ix,k)=(skeb_vwts(k,1)*skebv_save(i+isc-1,j+jsc-1,skeb_vpts(k,1))+skeb_vwts(k,2)*skebv_save(i+isc-1,j+jsc-1,skeb_vpts(k,2)))*vfact_skeb(k)
               ! debugif (k.EQ.1) Coupling(nb)%shum_wts(ix,k)=rpattern_skeb(1)%rnoise(i+isc-1,j+jsc-1)
            enddo
         enddo
      enddo
   endif
endif

if (do_rnda) then
   if (rpattern_rnda(1)%tau(1) .GE. 0) then
      do j=jsc,jec
           call random_gauss(tnoise,rpattern_rnda(1)%rstate2)
           rpattern_rnda(1)%wnoise(isc:iec,j,1)=tnoise
      enddo
   else
      rpattern_rnda(1)%wnoise(:,:,1)=rpattern_rnda(1)%rnoise
   endif
   call atmosphere_smooth_noise (rpattern_rnda(1)%wnoise,rpattern_rnda(1)%nsmth,ns_type,0) !  do smoothing for RNDA
   if (rpattern_rnda(1)%tau(1) .LT. 0) then
      rpattern_rnda(1)%rnoise(isc:iec,jsc:jec)= rpattern_rnda(1)%wnoise(isc:iec,jsc:jec,1)*rpattern_rnda(1)%stdev(1)
   else
      rpattern_rnda(1)%rnoise(isc:iec,jsc:jec)=rpattern_rnda(1)%phi(1)*rpattern_rnda(1)%rnoise(isc:iec,jsc:jec)+  &
      sqrt(1-rpattern_rnda(1)%phi(1)**2)*rpattern_rnda(1)%stdev(1)*rpattern_rnda(1)%wnoise(isc:iec,jsc:jec,1)
   endif
   call atmosphere_return_winds(rpattern_rnda(1)%rnoise(:,:),rndau_save(:,:),rndav_save(:,:),1,km=Model%levs,vwts=vfact_rnda)
endif
if (first_call) then
if (do_sfcperts) then
   do k=1,Model%nsfcpert
      if (stoch_ptype.NE.1) then
         call make_gridded(rpattern_sfc(k),tmp_wts,gis_stochy_atm)
         call advance_pattern(rpattern_sfc(k),gis_stochy_atm)
      endif
      if (stoch_ptype.NE.0) then
         do j=jsc,jec
            call random_gauss(tnoise,rpattern_sfc(k)%rstate2)
            rpattern_sfc(k)%wnoise(isc:iec,j,1)=tnoise
         enddo
         call atmosphere_smooth_noise (rpattern_sfc(k)%wnoise,rpattern_sfc(k)%nsmth,ns_type,1) !  do smoothing for SFC 
      endif
      nb = 1
      ix = 0
      if (stoch_ptype.NE.0) then
        do j = 1,Model%ny
           do i = 1,Model%nx
               ix=ix+1
               if (ix > gis_stochy_atm%blksz(nb)) then
                  nb = nb + 1
                  ix = 1
               endif
               if (stoch_ptype.EQ.0) then
                  Coupling(nb)%sfc_wts(ix,k)=tmp_wts(i,j)
               endif
               if (stoch_ptype.EQ.1) then
                  Coupling(nb)%sfc_wts(ix,k)= rpattern_sfc(k)%wnoise(i+isc-1,j+jsc-1,1)
               endif
            enddo
         enddo
      endif
   enddo
endif
endif
deallocate(tmp_wts)
deallocate(tnoise)
first_call=.false.

end subroutine run_stochastic_physics_atm
end module stochastic_physics_new    
