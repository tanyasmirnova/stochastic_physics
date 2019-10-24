module stochy_patterngenerator_new

! set up and initialize stochastic random patterns.

 use stochy_namelist_def_new
 use constants_mod, only : radius
#ifndef STOCHY_UNIT_TEST
 use mpp_mod
#endif
 use machine
 use mersenne_twister, only: random_setseed,random_gauss,random_stat

 implicit none 
 private
 public :: init_stochdata
!derived types
 type stochy_internal_state
    real, public, allocatable, dimension(:,:,:) :: legP    
    real, public, allocatable, dimension(:,:,:) :: legPtop
    real, public, allocatable, dimension(:) :: kenorm  
    real, public, allocatable, dimension(:) :: eon     
    real, public, allocatable, dimension(:) :: elonn1  
    real, public, allocatable, dimension(:) :: epstop  
    real, public, allocatable, dimension(:) :: eontop  
    real, public, allocatable, dimension(:,:) :: lons    
    real, public, allocatable, dimension(:,:) :: lats    
    integer            , public, allocatable, dimension(:) :: blksz   
    integer            , public, allocatable, dimension(:) :: indxm,indxn
    integer, public :: ntrunc,nlats,nlons,skeblevs,nblks,is,ie,js,je,me,nodes
    integer, public :: isd,ied,jsd,jed
 end type stochy_internal_state

 type random_pattern
    real, allocatable, dimension(:) :: stdev
    real, allocatable, dimension(:) :: tau
    real, allocatable, dimension(:) :: phi
    real, allocatable, dimension(:) :: var

    real, allocatable, dimension(:,:) :: varspectrum
    real, public, allocatable, dimension(:,:,:) :: wnoise   ! white noise holder
    real, public, allocatable, dimension(:,:) :: rnoise   ! red noise holder

    integer, public :: npatterns,nsmth
    real, allocatable, dimension(:,:,:), public :: spec_coeffu !  un-scaled coefficiends
    real, allocatable, dimension(:,:), public :: spec_coeff
    type(random_stat), public :: rstate,rstate2
 end type random_pattern
! available methods
 public :: make_gridded,make_gridded_vecta,make_gridded_vectc,advance_pattern

 integer, public :: nsppt=0
 integer, public :: nshum=0
 integer, public :: nskeb=0
 integer, public :: nrnda=0
 integer,allocatable:: sgn(:)
 real*8, public,allocatable :: sl(:)

 real,public, allocatable :: vfact_sppt(:),vfact_shum(:),vfact_skeb(:),vfact_rnda(:)
 real,public, allocatable :: skeb_vwts(:,:),skeb_vpts(:,:)
 real,public, allocatable :: skebu_save(:,:,:),skebv_save(:,:,:),rndau_save(:,:),rndav_save(:,:)
 logical,public :: first_call=.true.
 
 public :: random_pattern,stochy_internal_state


 type(random_pattern),       public,save :: rpattern_sppt(1),rpattern_shum(1),rpattern_skeb(1),rpattern_rnda(1),rpattern_sfc(6)
 type(stochy_internal_state),public,save :: gis_stochy_atm

 contains
 !subroutine init_stochdata(delt,input_nml_file,fn_nml,nlunit,lons,lats,gis_stochy)
 subroutine init_stochdata(delt,lons,lats,gis_stochy)

! initialize random patterns.  A spinup period of spinup_efolds times the
! temporal time scale is run for each pattern.  
!   integer, intent(in) :: nlunit
!   character(len=*),  intent(in) :: input_nml_file(:)
!   character(len=64), intent(in) :: fn_nml
#ifdef STOCHY_UNIT_TEST
      include 'mpif.h'
#endif
   real, intent(in) :: delt
   real(kind=kind_phys), intent(in) :: lons(:,:)
   real(kind=kind_phys), intent(in) :: lats(:,:)
   type(stochy_internal_state), intent(inout) :: gis_stochy

   real :: rn,f1,f2,amptmp(5)
   integer :: nn,k,nm,stochlun,ierr,iret,nx,ny,me,m,n,i,nlm,nsmthtmp,nperts
   me=gis_stochy%me

   if(me==0) then
       print*,'in init stochdata'
   endif
   !call stochy_namelist(me,size(input_nml_file,1),input_nml_file(:),fn_nml,nlunit,delt,iret)
   call stochy_namelist(me,delt,iret)
   if (do_sppt.EQ. .false. .AND. do_shum.EQ. .false..AND.do_skeb.EQ..false.) return
   stochlun=99
   nx=size(lons,1)
   ny=size(lons,2)
   gis_stochy%nlons=nx
   gis_stochy%nlats=ny
   allocate(gis_stochy%lons(nx,ny))
   allocate(gis_stochy%lats(nx,ny))
   gis_stochy%lons=lons
   gis_stochy%lats=lats
   print*,'in stochy_patterngenerator',gis_stochy%lons(2,2),gis_stochy%lats(2,2)
   if (stoch_ptype.NE.1) then
      nlm=(ntrunc+1)*(ntrunc+2)/2
      gis_stochy%ntrunc=ntrunc
      allocate(gis_stochy%indxm(nlm))
      allocate(gis_stochy%indxn(nlm))
      allocate(sgn(0:gis_stochy%ntrunc))
      i=1
      do m=0,ntrunc
         do n=m,ntrunc
            gis_stochy%indxm(i)=m
            gis_stochy%indxn(i)=n
            i=i+1
         enddo
      enddo
   endif
         
   iret=0
! determine number of random patterns to be used for each scheme.
   do n=1,size(sppt)
     if (sppt(n) > 0) then
        nsppt=nsppt+1
     else
        exit
     endif
   enddo
   if (me==0) print *,'nsppt = ',nsppt
   do n=1,size(shum)
     if (shum(n) > 0) then
        nshum=nshum+1
     else
        exit
     endif
   enddo
   if (me==0) print *,'nshum = ',nshum
   do n=1,size(skeb)
     if (skeb(n) > 0) then
        nskeb=nskeb+1
     else
        exit
     endif
   enddo
   if (me==0) print *,'nskeb = ',nskeb
   do n=1,size(rnda)
     if (rnda(n) > 0) then
        nrnda=nrnda+1
     else
        exit
     endif
   enddo
   if (me==0) print *,'nrnda = ',nrnda


!  if stochini is true, then read in pattern from a file
   if (me==0) then
      if (stochini) then
         print*,'opening stoch_ini'
         OPEN(stochlun,file='stoch_ini',form='unformatted',iostat=ierr,status='old')
#ifdef STOCHY_UNIT_TEST
         if (ierr .NE. 0) call mpi_abort(MPI_COMM_WORLD,12)
#else  
         if (ierr .NE. 0) call mpp_error(FATAL,'error opening stoch_ini')
#endif
      endif
   endif
   ! no spinup needed if initial patterns are defined correctly.
   if (nsppt > 0) then
      if (me==0) print *, 'Initialize random pattern for SPPT',nsppt
      call patterngenerator_init(sppt_lscale(1:nsppt),delt,sppt_tau(1:nsppt),sppt(1:nsppt),iseed_sppt,rpattern_sppt(1),nsppt,sppt_nsmth,gis_stochy)
      if (me==0) print *, 'returned'
      if (stochini) then
         !call read_pattern(rpattern_sppt(n),1,stochlun)
         print*,'ERROR not supported yet'
         STOP
      endif
   endif
   if (nshum > 0) then
      if (me==0) print *, 'Initialize random pattern for SHUM',nshum
      call patterngenerator_init(shum_lscale(1:nshum),delt,shum_tau(1:nshum),shum(1:nshum),iseed_shum,rpattern_shum(1),nshum,shum_nsmth,gis_stochy)
      if (stochini) then
         !call read_pattern(rpattern_shum(n),1,stochlun)
         print*,'ERROR not supported yet'
         STOP
      endif
   endif
   if (nrnda > 0) then
      if (me==0) print *, 'Initialize random pattern for RNDA',nrnda
      call patterngenerator_init(rnda_lscale(1:nrnda),delt,rnda_tau(1:nrnda),rnda(1:nrnda),iseed_rnda,rpattern_rnda(1),nrnda,rnda_nsmth,gis_stochy)
      if (stochini) then
         !call read_pattern(rpattern_rnda(n),1,stochlun)
         print*,'ERROR not supported yet'
         STOP
      endif
   endif
   if (nsfcpert > 0) then
      do k=1,6
        select case(k)
        case(1)
        amptmp=pertz0
        nsmthtmp=z0_nsmth

        case(2)
        amptmp=pertzt
        nsmthtmp=zt_nsmth

        case(3)
        amptmp=pertshc
        nsmthtmp=shc_nsmth

        case(4)
        amptmp=pertlai
        nsmthtmp=lai_nsmth

        case(5)
        amptmp=pertalb
        nsmthtmp=alb_nsmth

        case(6)
        amptmp=pertvegf
        nsmthtmp=vegf_nsmth

        end select

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
         nperts=0
         do n=1,size(amptmp)
            if (amptmp(n) > 0) then
               nperts=nperts+1
            else
               exit
            endif
         enddo
         if (me==0) print *, 'Initialize random pattern for SFCPERTS'
         call patterngenerator_init(sfc_lscale(1:nperts),delt,sfc_tau(1:nperts),amptmp(1:nperts),iseed_sfc,rpattern_sfc(k),nperts,nsmthtmp,gis_stochy)
         if (stochini) then
            !call read_pattern(rpattern_sfc(k),1,stochlun)
            print*,'ERROR not supported yet'
            STOP
         endif
      enddo 
   endif

   if (nskeb > 0) then
!  determine number of skeb levels to deal with temperoal/vertical correlations
   gis_stochy%skeblevs=nint(skeb_tau(1)/skebint*skeb_vdof)
! backscatter noise.
   if (me==0) print *, 'Initialize random pattern for SKEB',gis_stochy%skeblevs
   call patterngenerator_init(skeb_lscale(1:nskeb),delt,skeb_tau(1:nskeb),skeb(1:nskeb),iseed_skeb,rpattern_skeb(1),nskeb,skeb_nsmth,gis_stochy)
   if (stochini) then
      do k=1,gis_stochy%skeblevs
          !call read_pattern(rpattern_skeb(n),k,stochlun)
          print*,'ERROR not supported yet'
          STOP
       enddo
   endif
   if (stoch_ptype.NE.1) then
      allocate(gis_stochy%kenorm(0:ntrunc))
      gis_stochy%kenorm=1.
      if (skebnorm==0) then
         do n=0,gis_stochy%ntrunc
            rn = n*(n+1.)
            !gis_stochy%kenorm(n) = rn/radius**2
            gis_stochy%kenorm(n) = rn/radius**2
         enddo
         if (me==0) print*,'using streamfunction ',maxval(gis_stochy%kenorm(:)),minval(gis_stochy%kenorm(:))
      endif
      if (skebnorm==1) then
        do n=0,gis_stochy%ntrunc
           rn = n*(n+1.)
           gis_stochy%kenorm(n) = sqrt(rn)/radius
        enddo
        if (me==0) print*,'using kenorm ',maxval(gis_stochy%kenorm(:)),minval(gis_stochy%kenorm(:))
      endif
   endif
endif ! skeb > 0
if (me==0 .and. stochini) CLOSE(stochlun)

 if (stoch_ptype.NE.1) call generate_legendre_polynomials(gis_stochy)

end subroutine init_stochdata


  
 subroutine patterngenerator_init(lscale, delt, tscale, stdev, iseed, rpattern, npatterns,nsmth,gis_stochy)
#ifdef STOCHY_UNIT_TEST
      include 'mpif.h'
#endif
   real, intent(in) :: lscale(npatterns),tscale(npatterns),stdev(npatterns)
   real, intent(in) :: delt
   integer, intent(in) :: npatterns,nsmth
   integer(8), intent(inout) :: iseed
   type(random_pattern), intent(out) :: rpattern
   type(stochy_internal_state), intent(inout) :: gis_stochy
! locals
   integer i,j,m,n,nlm
   integer(8) count, count_rate, count_max, count_trunc
   integer(8) :: iscale = 10000000000
   integer :: npes
   integer,                allocatable :: pelist(:)
   integer count4,seed,me,ierr
    
!  propagate seed supplied from namelist to all patterns... (can remove PJP)
   me=gis_stochy%me
   if (me==0) print*,'in patterngenerator_init',nsmth,stoch_ptype
   ! seed computed on root, then bcast to all tasks and set.
   if (me==0) then
      if (iseed == 0) then
        ! generate a random seed from system clock and ens member number
        call system_clock(count, count_rate, count_max)
        ! iseed is elapsed time since unix epoch began (secs)
        ! truncate to 4 byte integer
        count_trunc = iscale*(count/iscale)
        count4 = count - count_trunc !+ member_id
        print *,'using seed',count4
      else
        count4 = mod(iseed + 2147483648, 4294967296) - 2147483648
        print *,'using seed',count4,iseed!,member_id
      endif
   endif
   ! broadcast seed to all tasks.
#ifdef STOCHY_UNIT_TEST
      call mpi_bcast( count4,1, mpi_integer,0,mpi_comm_world,ierr )
#else
      npes = mpp_npes()
      allocate(pelist(0:npes-1))
      call mpp_get_current_pelist(pelist)
      call mpp_broadcast( count4, mpp_root_pe(),pelist )
      deallocate(pelist)
#endif
   seed = count4
   rpattern%npatterns=npatterns 
   !print*,'allocated spec_coeffu',nlm,2,npatterns
   allocate(rpattern%tau(npatterns))
   allocate(rpattern%stdev(npatterns))
   allocate(rpattern%phi(npatterns))
   rpattern%tau(1:npatterns) = tscale(1:npatterns)
   rpattern%phi(1:npatterns) = exp(-delt/tscale(1:npatterns))
   rpattern%stdev(1:npatterns) = stdev(1:npatterns)
   if (stoch_ptype.NE.1) then
      nlm=(ntrunc+1)*(ntrunc+2)/2.0
      allocate(rpattern%spec_coeffu(nlm,2,npatterns))
      allocate(rpattern%spec_coeff(nlm,2))
      allocate(rpattern%var(npatterns))
      allocate(rpattern%varspectrum(0:gis_stochy%ntrunc,npatterns))

      rpattern%spec_coeffu(:,:,:)=0.
      rpattern%spec_coeff(:,:)=0.
   
      ! set seed (to be the same) on all tasks. Save random state.
      call random_setseed(seed,rpattern%rstate)
   
      do i=1,npatterns
         do n=0,gis_stochy%ntrunc
            if (n.GT. 0)rpattern%varspectrum(n,i) = gis_stochy%ntrunc*exp((lscale(i)*0.25)**2*(-n*(n+1))/(radius**2))
            !rpattern%varspectrum(n,i) = gis_stochy%ntrunc*exp((lscale(i)*0.25)**2*(-n*(n+1))/(radius**2))
            !rpattern%varspectrum(n,i) = sqrt(gis_stochy%ntrunc*exp((lscale(i)**2)*(-n*(n+1))/(4.*radius**2)))
            !rpattern%varspectrum(n,i) = sqrt(gis_stochy%ntrunc*exp((lscale(i)**2)*(-n*(n+1))/(4.*radius**2)))
         enddo
      enddo

      call computevarspec(rpattern,gis_stochy%ntrunc)
      call populate_pattern(rpattern,gis_stochy)
      print*,'varspectrun=',rpattern%varspectrum(:,1),rpattern%var(1)
   endif
   if (stoch_ptype.GT.0) then 
      seed=seed+nint(100*(gis_stochy%lons(2,2)+gis_stochy%lats(2,2)))
      print*,'setting seeds',seed,me,gis_stochy%lons(2,2),gis_stochy%lats(2,2)
      call random_setseed(seed,rpattern%rstate2)
       rpattern%nsmth=nsmth
       print*,'allocating noise',gis_stochy%isd,gis_stochy%ied,gis_stochy%jsd,gis_stochy%jed
       allocate(rpattern%wnoise(gis_stochy%isd:gis_stochy%ied,gis_stochy%jsd:gis_stochy%jed,1))
       allocate(rpattern%rnoise(gis_stochy%is:gis_stochy%ie,gis_stochy%js:gis_stochy%je))
   endif


 end subroutine patterngenerator_init


subroutine generate_legendre_polynomials(gis_stochy)
type(stochy_internal_state), intent(inout) :: gis_stochy
! locals
integer                         :: m,n,ntrunc,nlm,i,im1
real, allocatable :: theta(:,:),mu(:,:),eps(:),epstop(:)
real :: rerth

rerth=6.3712e6

ntrunc=gis_stochy%ntrunc
nlm=(ntrunc+1)*(ntrunc+2)/2.0
allocate(gis_stochy%legP(nlm,gis_stochy%nlons,gis_stochy%nlats))
allocate(gis_stochy%legPtop(0:ntrunc,gis_stochy%nlons,gis_stochy%nlats))
allocate(gis_stochy%eon(nlm))
allocate(gis_stochy%elonn1(nlm))
allocate(gis_stochy%eontop(0:ntrunc))
allocate(gis_stochy%epstop(0:ntrunc))
allocate(theta(gis_stochy%nlons,gis_stochy%nlats))
allocate(mu(gis_stochy%nlons,gis_stochy%nlats))
allocate(eps(nlm))
theta=sin(gis_stochy%lats)
mu=cos(gis_stochy%lats)

eps(:)=0.0
gis_stochy%eon(:)=0.0
gis_stochy%elonn1(:)=0.0


DO m=0,ntrunc
   i=m*(2*ntrunc-(m-1))/2+m+1
   gis_stochy%elonn1(i)=rerth/real(m+1.0)
ENDDO

DO m=0,ntrunc
   DO n=m+1,ntrunc
      i=m*(2*ntrunc-(m-1))/2+n+1
      gis_stochy%elonn1(i)=rerth*m/(real(n)*(n+1.0))
      gis_stochy%eon(i)=rerth/real(n)*sqrt(real(n**2-m**2)/real(4.0*n**2-1.0))
   ENDDO
ENDDO
      
DO m=0,ntrunc
   n=ntrunc+1
   gis_stochy%epstop(m)=sqrt((n**2-m**2)/real(4.0*n**2.0-1.0))
   gis_stochy%eontop(m)=rerth/(ntrunc+1.0)*gis_stochy%epstop(m)
ENDDO

DO m=0,ntrunc-1
   DO n=m+1,ntrunc
      i=m*(2*ntrunc-(m-1))/2+n+1
      eps(i)=sqrt((n**2-m**2)/real(4.0*n**2-1.0))
   ENDDO
ENDDO


gis_stochy%legP(:,:,:)=0
!  ITERATIVELY COMPUTE PLN(L,L) (BOTTOM HYPOTENUSE OF DOMAIN)
gis_stochy%legP(1,:,:)=sqrt(0.5)
i=1
do m=1,ntrunc
   im1=i
   i=m*(2*ntrunc-(m-1))/2+m+1
   gis_stochy%legP(i,:,:)=mu(:,:)*sqrt((2.0*m+1.0)/real(2.0*m))*gis_stochy%legP(im1,:,:)
enddo
! COMPUTE PLN(L,L+1) (DIAGONAL NEXT TO BOTTOM HYPOTENUSE OF DOMAIN)
do m=0,ntrunc-1
   i=m*(2*ntrunc-(m-1))/2+m+2
   gis_stochy%legP(i,:,:)=theta*gis_stochy%legP(i-1,:,:)/eps(i)
enddo

!  COMPUTE REMAINING pln IN SPECTRAL DOMAIN
do n=2,ntrunc
   do m=0,ntrunc-n
      i=m*(2*ntrunc-(m-1))/2+m+n+1
      gis_stochy%legP(i,:,:)=(theta*gis_stochy%legP(i-1,:,:)-eps(i-1)*gis_stochy%legP(i-2,:,:))/eps(i)
   enddo
enddo

!  COMPUTE POLYNOMIALS OVER TOP OF SPECTRAL DOMAIN
do m=0,ntrunc
   n=ntrunc+1-m
   i=m*(2*ntrunc-(m-1))/2+m+n+1
   gis_stochy%legPtop(m,:,:)=(theta*gis_stochy%legP(i-1,:,:)-eps(i-1)*gis_stochy%legP(i-2,:,:))/gis_stochy%epstop(m)
enddo
   
deallocate(theta)
deallocate(mu)
deallocate(eps)
!   sign handling
if (mod(ntrunc,2)==0) gis_stochy%legPtop(:,:,:)=-1*gis_stochy%legPtop(:,:,:)
if (mod(ntrunc,2)==1) then
   DO i=0,ntrunc-1,2
      sgn(i)=1
      sgn(i+1)=-1
   ENDDO
else
   DO i=0,ntrunc-2,2
      sgn(i)=1
      sgn(i+1)=-1
   ENDDO
   sgn(ntrunc)=1
endif
do i=1,nlm
   gis_stochy%legP(i,:,:)=gis_stochy%legP(i,:,:)*sgn(gis_stochy%indxn(i))
enddo
return
end subroutine generate_legendre_polynomials

! populate spectral pattern with unit variance for all wave numbers
subroutine populate_pattern(rpattern,gis_stochy)
   type(random_pattern), intent(inout) :: rpattern
   type(stochy_internal_state), intent(inout) :: gis_stochy
!locals
   real(kind_dbl_prec)  :: noise((gis_stochy%ntrunc+1)*(gis_stochy%ntrunc+2)),spec_coeff_tmp((gis_stochy%ntrunc+1)*(gis_stochy%ntrunc+2)/2,2,rpattern%npatterns)
   integer:: i,nlm,np

   nlm=(gis_stochy%ntrunc+1)*(gis_stochy%ntrunc+2)/2.0
   DO np=1,rpattern%npatterns
      call random_gauss(noise,rpattern%rstate)
      rpattern%spec_coeffu(:,1,np)=noise(1:nlm)
      rpattern%spec_coeffu(:,2,np)=noise(nlm+1:2*nlm)
      spec_coeff_tmp(:,:,np)=rpattern%spec_coeffu(:,:,np)/sqrt(1.0/float(gis_stochy%ntrunc))
   ENDDO
   spec_coeff_tmp(1,:,:)=0
   rpattern%spec_coeff(:,:)=0
   DO np=1,rpattern%npatterns
      do i=2,nlm
         spec_coeff_tmp(i,1,np) = spec_coeff_tmp(i,1,np)/sqrt(2.*gis_stochy%indxn(i)+1)*rpattern%varspectrum(gis_stochy%indxn(i),np)
         spec_coeff_tmp(i,2,np) = spec_coeff_tmp(i,2,np)/sqrt(2.*gis_stochy%indxn(i)+1)*rpattern%varspectrum(gis_stochy%indxn(i),np)
      enddo
      spec_coeff_tmp(2:gis_stochy%ntrunc+1,:,np)=spec_coeff_tmp(2:gis_stochy%ntrunc+1,:,np)*sqrt(2.0)
      spec_coeff_tmp(gis_stochy%ntrunc+2:nlm,:,np)=spec_coeff_tmp(gis_stochy%ntrunc+2:nlm,:,np)*2
      rpattern%spec_coeff(:,:)=rpattern%spec_coeff(:,:)+spec_coeff_tmp(:,:,np)/sqrt(rpattern%var(np))*rpattern%stdev(np)
   ENDDO
return
end subroutine populate_pattern

! advance spectral pattern with unit variance for all wave numbers
subroutine advance_pattern(rpattern,gis_stochy)
   type(random_pattern), intent(inout) :: rpattern
   type(stochy_internal_state), intent(inout) :: gis_stochy
!locals
   real(kind_dbl_prec)           :: noise((gis_stochy%ntrunc+1)*(gis_stochy%ntrunc+2)),spec_coeff_tmp((gis_stochy%ntrunc+1)*(gis_stochy%ntrunc+2)/2.0,2,rpattern%npatterns)
   integer:: i,nlm,np

   nlm=(gis_stochy%ntrunc+1)*(gis_stochy%ntrunc+2)/2.0
   DO np=1,rpattern%npatterns
      call random_gauss(noise,rpattern%rstate)
      ! adviance ar(1) with autocorrelation varying by wavenumber
      do i=2,nlm
         rpattern%spec_coeffu(i,1,np)=rpattern%phi(np)*rpattern%spec_coeffu(i,1,np)+sqrt(1-rpattern%phi(np)**2)*noise(i)
         rpattern%spec_coeffu(i,2,np)=rpattern%phi(np)*rpattern%spec_coeffu(i,2,np)+sqrt(1-rpattern%phi(np)**2)*noise(i+nlm)
      enddo
   ENDDO
   spec_coeff_tmp=rpattern%spec_coeffu/sqrt(1.0/float(gis_stochy%ntrunc))
   spec_coeff_tmp(1,:,:)=0
   rpattern%spec_coeff(:,:)=0
   DO np=1,rpattern%npatterns
      do i=2,nlm
         spec_coeff_tmp(i,1,np) = spec_coeff_tmp(i,1,np)/sqrt(2.*gis_stochy%indxn(i)+1)*rpattern%varspectrum(gis_stochy%indxn(i),np)
         spec_coeff_tmp(i,2,np) = spec_coeff_tmp(i,2,np)/sqrt(2.*gis_stochy%indxn(i)+1)*rpattern%varspectrum(gis_stochy%indxn(i),np)
      enddo
      spec_coeff_tmp(2:gis_stochy%ntrunc+1,:,np)=spec_coeff_tmp(2:gis_stochy%ntrunc+1,:,np)*sqrt(2.0)
      spec_coeff_tmp(gis_stochy%ntrunc+2:nlm,:,np)=spec_coeff_tmp(gis_stochy%ntrunc+2:nlm,:,np)*2
      rpattern%spec_coeff(:,:)=rpattern%spec_coeff(:,:)+spec_coeff_tmp(:,:,np)/sqrt(rpattern%var(np))*rpattern%stdev(np)
   ENDDO
   if (gis_stochy%me==0) print*,'in advance',SUM(rpattern%var(1:rpattern%npatterns)),SUM(rpattern%stdev(1:rpattern%npatterns)),SUM(rpattern%spec_coeff(:,:)**2)
return
end subroutine advance_pattern


! do spectral transform directly to physics grid
subroutine make_gridded(rpattern,grid_out,gis_stochy)
   type(random_pattern), intent(inout) :: rpattern
   real,intent(inout) :: grid_out(:,:)
   type(stochy_internal_state), intent(in) :: gis_stochy
!locals
   real        :: CoeffA,CoeffB
   integer:: i,j,m,n,ntrunc,x,y,i1,i2
   
   ntrunc=gis_stochy%ntrunc
   grid_out(:,:)=0.0
   do x=1,gis_stochy%nlons
      do y=1,gis_stochy%nlats
         !m=0
         i1=1
         i2=ntrunc+1
         grid_out(x,y)=SUM(gis_stochy%legP(2:ntrunc+1,x,y)*rpattern%spec_coeff(2:ntrunc+1,1))
         do m=1,ntrunc
            i1=i2+1
            i2=i1+ntrunc-m
            CoeffA=sgn(m)*SUM(gis_stochy%legP(i1:i2,x,y)*rpattern%spec_coeff(i1:i2,1))
            CoeffB=sgn(m)*SUM(gis_stochy%legP(i1:i2,x,y)*rpattern%spec_coeff(i1:i2,2))
            grid_out(x,y)=grid_out(x,y)+CoeffA*cos(m*gis_stochy%lons(x,y))-CoeffB*sin(m*gis_stochy%lons(x,y))   
         enddo
      enddo
   enddo
   if (gis_stochy%me==0) print*,'in make_gridded',SUM(rpattern%var(1:rpattern%npatterns)),SUM(rpattern%stdev(1:rpattern%npatterns)),SUM(grid_out(:,:)**2)
return
end subroutine make_gridded

subroutine make_gridded_vecta(rpattern,gridu_out,gridv_out,gis_stochy)
! do spectral transform of winds directly to physics grid (need to modify for
! randtran, which uses winds on staggered grid, and cube relative
   type(random_pattern), intent(inout) :: rpattern
   real,intent(inout) :: gridu_out(:,:),gridv_out(:,:)
   type(stochy_internal_state), intent(in) :: gis_stochy
!locals
   real     :: specu((gis_stochy%ntrunc+1)*(gis_stochy%ntrunc+2)/2,2)
   real     :: specv((gis_stochy%ntrunc+1)*(gis_stochy%ntrunc+2)/2,2)
   real     :: specutop(0:gis_stochy%ntrunc,2)
   real     :: CoeffA,CoeffB
   integer:: i,j,m,n,ntrunc,i1,i2,nlm

   ntrunc=gis_stochy%ntrunc
   nlm=(ntrunc+1)*(ntrunc+2)/2

! get zonal winds from vorticy pattern
   do i=2,nlm
      rpattern%spec_coeff(i,1)=rpattern%spec_coeff(i,1)*gis_stochy%kenorm(gis_stochy%indxn(i))
      rpattern%spec_coeff(i,2)=rpattern%spec_coeff(i,2)*gis_stochy%kenorm(gis_stochy%indxn(i))
   enddo
!  m=0, n=0
   specu(1,1)=gis_stochy%eon(2)*rpattern%spec_coeff(2,1)
   specu(1,2)=gis_stochy%eon(2)*rpattern%spec_coeff(2,2)
   i=1
   DO m=0,ntrunc
      DO  n=m,ntrunc-1
         IF(n>0) THEN
           specu(i,1)=gis_stochy%eon(i+1)*rpattern%spec_coeff(i+1,1)-gis_stochy%eon(i)*rpattern%spec_coeff(i-1,1)
           specu(i,2)=gis_stochy%eon(i+1)*rpattern%spec_coeff(i+1,2)-gis_stochy%eon(i)*rpattern%spec_coeff(i-1,2)
         ENDIF
         i=i+1 ! increment i
      ENDDO
      specu(i,1)=-1*gis_stochy%eon(i)*rpattern%spec_coeff(i-1,1)
      specu(i,2)=-1*gis_stochy%eon(i)*rpattern%spec_coeff(i-1,2)
      specutop(m,1)=-1*gis_stochy%eontop(m)*rpattern%spec_coeff(i,1)
      specutop(m,2)=-1*gis_stochy%eontop(m)*rpattern%spec_coeff(i,2)
      i=i+1 ! increment i for ntrunc row
   ENDDO
   DO i=ntrunc+2,nlm
      specv(i,1)=gis_stochy%elonn1(i)*rpattern%spec_coeff(i,2)
      specv(i,2)=-1*gis_stochy%elonn1(i)*rpattern%spec_coeff(i,1)
   ENDDO

   gridu_out(:,:)=0.0
   gridv_out(:,:)=0.0
   do j=1,gis_stochy%nlats
      do i=1,gis_stochy%nlons
         m=0
         i1=1
         i2=ntrunc+1
         gridu_out(i,j)=SUM(gis_stochy%legP(i1:i2,i,j)*specu(i1:i2,1))+gis_stochy%legPtop(m,i,j)*specutop(m,1)

         do m=1,ntrunc
            i1=i2+1
            i2=i1+ntrunc-m
            CoeffA=sgn(m)*SUM(gis_stochy%legP(i1:i2,i,j)*specv(i1:i2,1))
            CoeffB=sgn(m)*SUM(gis_stochy%legP(i1:i2,i,j)*specv(i1:i2,2))
            gridv_out(i,j)=gridv_out(i,j)+CoeffA*cos(m*gis_stochy%lons(i,j))-CoeffB*sin(m*gis_stochy%lons(i,j))   

            CoeffA=sgn(m)*(SUM(gis_stochy%legP(i1:i2,i,j)*specu(i1:i2,1))+gis_stochy%legPtop(m,i,j)*specutop(m,1))
            CoeffB=sgn(m)*(SUM(gis_stochy%legP(i1:i2,i,j)*specu(i1:i2,2))+gis_stochy%legPtop(m,i,j)*specutop(m,2))
            gridu_out(i,j)=gridu_out(i,j)+CoeffA*cos(m*gis_stochy%lons(i,j))-CoeffB*sin(m*gis_stochy%lons(i,j))   
         enddo
         if (cos(gis_stochy%lats(i,j))>0.0) then
            gridu_out(i,j)=gridu_out(i,j)/cos(gis_stochy%lats(i,j))
            gridv_out(i,j)=-1*gridv_out(i,j)/cos(gis_stochy%lats(i,j))
         else 
            print*,'Have cos_lat=0',i,j
         endif
      enddo
   enddo
return
end subroutine make_gridded_vecta

subroutine make_gridded_vectc(rpattern,gridu_out,gridv_out)
! do spectral transform of winds directly to physics grid (need to modify for
! randtran, which uses winds on staggered grid, and cube relative
   type(random_pattern), intent(inout) :: rpattern
   real,intent(inout) :: gridu_out(:,:),gridv_out(:,:)
!locals
! stub
end subroutine make_gridded_vectc

subroutine computevarspec(rpattern,ntrunc)
   ! compute globally integrated variance from spectral coefficients
   type(random_pattern), intent(inout) :: rpattern
   integer, intent(in)  :: ntrunc
   real :: sqrt2x2
   integer m,n,np
   sqrt2x2=2.0*sqrt(2.0)
   rpattern%var = 0.
   do np=1,rpattern%npatterns
      do m=0,ntrunc
         do n=m,ntrunc
            if (m == 0) then
                rpattern%var(np) = rpattern%var(np) + sqrt2x2*rpattern%varspectrum(n,np)**2/((2.*n+1)/real(ntrunc))
            else
                rpattern%var(np) = rpattern%var(np) + 2.0*rpattern%varspectrum(n,np)**2/((2.*n+1)/real(ntrunc))
            endif
         enddo
      enddo
   enddo
end subroutine computevarspec

end module stochy_patterngenerator_new
