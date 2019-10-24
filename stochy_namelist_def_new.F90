      module stochy_namelist_def_new
!
      use machine
      implicit none
      
      public
      integer :: nsrnda,nsskeb,ntrunc

! pjp stochastic phyics
      integer :: skeb_npass,sppt_nsmth,shum_nsmth,rnda_nsmth,skeb_nsmth,ns_type
      integer :: stoch_ptype
      logical :: sppt_sfclimit

      real :: skeb_sigtop1,skeb_sigtop2,          &
              rnda_sigtop1,rnda_sigtop2,          &
              sppt_sigtop1,sppt_sigtop2,shum_sigefold, &
              skeb_vdof,rnda_vdof
      real fhstoch,skeb_diss_smooth,skebint,skebnorm,rndaint,rndanorm
      real, dimension(5) :: skeb,skeb_lscale,skeb_tau
      real, dimension(5) :: rnda,rnda_lscale,rnda_tau
      real, dimension(5) :: sppt,sppt_lscale,sppt_tau
      real, dimension(5) :: ocn_rp,ocn_rp_lscale,ocn_rp_tau
      real, dimension(5) :: shum,shum_lscale,shum_tau
      real, dimension(5) :: sfc_lscale,sfc_tau
      real(kind=kind_dbl_prec), dimension(5) :: pertz0,pertshc,pertzt
      real(kind=kind_dbl_prec), dimension(5) :: pertlai,pertvegf,pertalb
      integer :: z0_nsmth,shc_nsmth,zt_nsmth
      integer :: lai_nsmth,vegf_nsmth,alb_nsmth
      integer::  nsfcpert
      integer :: lon_s,lat_s ! legacy
      integer(8) ::iseed_sfc
      logical :: sppt_land
      logical :: do_sfcperts
      integer,dimension(5) ::skeb_vfilt,rnda_vfilt
      integer(8) ::iseed_sppt,iseed_shum,iseed_skeb,iseed_rnda,iseed_ocn_rp
      logical :: stochini,sppt_logit,sppt_blend
      logical :: do_shum,do_sppt,do_skeb,do_rnda,use_zmtnblck

      contains

      !subroutine stochy_namelist(me,sz_nml,input_nml_file,fn_nml,nlunit,deltim,iret)
      subroutine stochy_namelist(me,deltim,iret)
      
      implicit none

 
      integer,              intent(out)   :: iret
      integer,              intent(in)    ::  me
      !integer,              intent(in)    :: nlunit,me,sz_nml
      !character(len=*),     intent(in)    :: input_nml_file(sz_nml)
      !character(len=64),    intent(in)    :: fn_nml
      real,                 intent(in)    :: deltim
      real tol
      integer k,ios

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      namelist /nam_stochy/ntrunc,sppt,sppt_tau,sppt_lscale,sppt_logit, &
      iseed_shum,iseed_sppt,shum,shum_tau,& 
      shum_lscale,fhstoch,stochini,sppt_sfclimit, &
      skeb,skeb_tau,skeb_vdof,skeb_lscale,iseed_skeb,skeb_vfilt,skeb_diss_smooth, &
      rnda,rnda_tau,rnda_vdof,rnda_lscale,iseed_rnda,rnda_vfilt, &
      ocn_rp,ocn_rp_lscale,ocn_rp_tau,iseed_ocn_rp,              &
      skeb_sigtop1,skeb_sigtop2,skebnorm,sppt_sigtop1,sppt_sigtop2,&
      rnda_sigtop1,rnda_sigtop2,rndanorm,                          &
      shum_sigefold,skebint,skeb_npass,use_zmtnblck,sppt_blend,sppt_nsmth,shum_nsmth,skeb_nsmth,ns_type,&
      rndaint,rnda_nsmth,stoch_ptype,nsfcpert,pertz0,pertshc,pertzt,pertlai, & ! mg, sfcperts
      z0_nsmth,shc_nsmth,zt_nsmth,lai_nsmth,vegf_nsmth,alb_nsmth,&
      pertvegf,pertalb,iseed_sfc,sfc_tau,sfc_lscale,sppt_land,lon_s,lat_s
       tol=0.01  ! tolerance for calculations
!     spectral resolution defintion
      ntrunc=-999
      ! can specify up to 5 values for the stochastic physics parameters 
      ! (each is an array of length 5)
      sppt             = -999.  ! stochastic physics tendency amplitude
      shum             = -999.  ! stochastic boundary layer spf hum amp   
      skeb             = -999.  ! stochastic KE backscatter amplitude
      rnda             = -999.  ! stochastic KE backscatter amplitude
! logicals
      do_sppt = .false.
      sppt_blend = .false.
      use_zmtnblck = .false.
      do_shum = .false.
      do_skeb = .false.
      do_rnda = .false.
      ! mg, sfcperts
      do_sfcperts = .false.
      sppt_land = .false.
      nsfcpert = 0
! for sfcperts random patterns
      sfc_lscale  = -999.       ! length scales
      sfc_tau     = -999.       ! time scales
      iseed_sfc   = 0           ! random seeds (if 0 use system clock)
! for SKEB random patterns.
      stoch_ptype       = 0  ! pattern type 0 = spectral, 1 = smooth noise, 2 = both
      skeb_vfilt       = 0
      rnda_vfilt       = 0
      skebint          = 0
      rndaint          = 0
      skeb_npass       = 11  ! number of passes of smoother for dissipation estiamte
      sppt_nsmth       = -1  ! number of passes of smoother for dissipation estiamte
      shum_nsmth       = -1 ! number of passes of smoother for dissipation estiamte
      skeb_nsmth       = -1  ! number of passes of smoother for dissipation estiamte
      rnda_nsmth       = -1  ! number of passes of smoother for dissipation estiamte
      z0_nsmth         = -1  ! number of passes of smoother for dissipation estiamte
      zt_nsmth         = -1  ! number of passes of smoother for dissipation estiamte
      shc_nsmth        = -1  ! number of passes of smoother for dissipation estiamte
      lai_nsmth        = -1  ! number of passes of smoother for dissipation estiamte
      vegf_nsmth       = -1  ! number of passes of smoother for dissipation estiamte
      alb_nsmth        = -1  ! number of passes of smoother for dissipation estiamte
      sppt_tau         = -999.  ! time scales
      shum_tau         = -999.
      skeb_tau         = -999.
      rnda_tau         = -999.
      skeb_vdof        = 5 ! proxy for vertical correlation, 5 is close to 40 passes of the 1-2-1 filter in the GFS
      rnda_vdof        = 0 ! proxy for vertical correlation, 5 is close to 40 passes of the 1-2-1 filter in the GFS
      skebnorm         = 0  ! 0 - random pattern is stream function, 1- pattern is kenorm, 2- pattern is vorticity
      rndanorm         = 0  ! 0 - random pattern is stream function, 1- pattern is kenorm, 2- pattern is vorticity
      sppt_lscale      = -999.  ! length scales
      shum_lscale      = -999.
      skeb_lscale      = -999.
      rnda_lscale      = -999.
      iseed_sppt       = 0      ! random seeds (if 0 use system clock)
      iseed_shum       = 0
      iseed_skeb       = 0
      iseed_rnda       = 0

      ns_type          = 0  ! 1 = 2nd order lapalacian,  anything else = 9-point box mean
! parameters to control vertical tapering of stochastic physics with
! height
      sppt_sigtop1 = 0.1
      sppt_sigtop2 = 0.025
      skeb_sigtop1 = 0.1
      skeb_sigtop2 = 0.025
      rnda_sigtop1 = 0.1
      rnda_sigtop2 = 0.025
      shum_sigefold = 0.2
! reduce amplitude of sppt near surface (lowest 2 levels)
      sppt_sfclimit = .false.
! gaussian or power law variance spectrum for skeb (0: gaussian, 1:
! power law). If power law, skeb_lscale interpreted as a power not a
! length scale.
      sppt_logit        = .false. ! logit transform for sppt to bounded interval [-1,+1]
      fhstoch           = -999.0  ! forecast hour to dump random patterns
      stochini          = .false. ! true= read in pattern, false=initialize from seed

!!#ifdef INTERNAL_FILE_NML
!!      read(input_nml_file, nml=nam_stochy)
!!#else
!!      rewind (nlunit)
!!      open (unit=nlunit, file=fn_nml, READONLY, status='OLD', iostat=ios)
!!      read(nlunit,nam_stochy)
!!#endif
      open (unit=565, file='input.nml', READONLY, status='OLD', iostat=ios)
      read(565,nam_stochy)

      if (me == 0) then
      print *,' in stochy_namelist'
      print*,'skeb=',skeb
      endif

! PJP stochastic physics additions
      IF (sppt(1) > 0 ) THEN
        do_sppt=.true.
      ENDIF
      IF (shum(1) > 0 ) THEN
        do_shum=.true.
!     shum parameter has units of 1/hour, to remove time step
!     dependence.
!     change shum parameter units from per hour to per timestep
         DO k=1,5
            IF (shum(k) .gt. 0.0) shum(k)=shum(k)*deltim/3600.0
         ENDDO
      ENDIF
      IF (rnda(1) > 0 ) THEN
         do_rnda=.true.
      ENDIF
      IF (skeb(1) > 0 ) THEN
         do_skeb=.true.
     !    if (skebnorm==0) then ! stream function norm
     !       skeb=skeb*1.111e3*sqrt(deltim)
     !    endif
     !    if (skebnorm==1) then ! ke norm
     !       skeb=skeb*0.00222*sqrt(deltim)
     !    endif
     !    if (skebnorm==2) then ! vorticty function norm
     !       skeb=skeb*1.111e-9*sqrt(deltim)
     !    endif
      ENDIF
!    compute frequencty to estimate dissipation timescale
      IF (skebint == 0.) skebint=deltim
      nsskeb=nint(skebint/deltim)                              ! skebint in seconds
      IF(nsskeb<=0 .or. abs(nsskeb-skebint/deltim)>tol) THEN
         WRITE(0,*) "SKEB interval is invalid",skebint
        iret=9
        return
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (rndaint == 0.) rndaint=deltim
      nsrnda=nint(rndaint/deltim)                              ! rndaint in seconds
      IF(nsrnda<=0 .or. abs(nsrnda-rndaint/deltim)>tol) THEN
         WRITE(0,*) "RNDA interval is invalid",rndaint
        iret=9
        return
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (skebint == 0.) skebint=deltim
      nsskeb=nint(skebint/deltim)                              ! skebint in seconds
      IF(nsskeb<=0 .or. abs(nsskeb-skebint/deltim)>tol) THEN
         WRITE(0,*) "SKEB interval is invalid",skebint
        iret=9
        return
      ENDIF
     ! mg, sfcperts
      IF (pertz0(1) > 0 .OR. pertshc(1) > 0 .OR. pertzt(1) > 0 .OR. &
          pertlai(1) > 0 .OR. pertvegf(1) > 0 .OR. pertalb(1) > 0) THEN
        do_sfcperts=.true.
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  All checks are successful.
!
      if (me == 0) then
         print *, 'stochastic physics'
         print *, ' do_sppt : ', do_sppt
         print *, ' do_shum : ', do_shum
         print *, ' do_skeb : ', do_skeb
         print *, ' do_rnda : ', do_rnda
         print *, ' do_sfcperts : ', do_sfcperts
      endif
      iret = 0
!
      return
      end subroutine stochy_namelist 
end module stochy_namelist_def_new
