program lw_driver

 use netcdf
 use rrtmg_lw_rad, only: rrtmg_lw  !  RRTMG Code
 use rrtmg_lw_init, only: rrtmg_lw_ini

 implicit none

 !Data dimensions
 integer :: np, im, jm, lm, nb_rrtmg
 integer :: is, js, ie, je

 !Configuration for the rrtmg scheme
 integer :: corespernode, partition_size
 integer :: icld, idrv, inflglw, iceflglw, liqflglw
 integer :: doy
 integer :: i, j, k, ij, lv

 !Constants
 real, parameter :: MAPL_GRAV  = 9.80665
 real, parameter :: MAPL_RUNIV = 8314.47
 real, parameter :: MAPL_AIRMW = 28.965
 real, parameter :: MAPL_H2OMW = 18.015
 real, parameter :: MAPL_O3MW  = 47.9982
 real, parameter :: MAPL_RDRY  = MAPL_RUNIV/MAPL_AIRMW
 real, parameter :: MAPL_RGAS  = MAPL_RDRY
 real, parameter :: MAPL_CPDRY = 3.5*MAPL_RDRY
 real, parameter :: MAPL_KAPPA = MAPL_RDRY/MAPL_CPDRY

 !Inputs from the NetCDF file
 character(len=2048) :: filename
 real, allocatable, dimension(:,:,:)   :: dp, t, q, qi, ql, ri, rl, o3, fcld
 real, allocatable, dimension(:,:)     :: ts, emis, lats

 integer :: ncid, eid

 !Intermediate variable transforms
 integer, parameter :: kice     = 1
 integer, parameter :: kliquid  = 2
 real, allocatable, dimension(:,:,:,:) :: reff, cwc
 real, allocatable, dimension(:,:,:)   :: ple, pl, tlev

 !Main inputs for the scheme
 real :: xx
 real, allocatable, dimension(:,:)     :: pl_r, ple_r, t_r, tlev_r
 real, allocatable, dimension(:)       :: tsfc
 real, allocatable, dimension(:,:)     :: q_r, o3_r, emiss, fcld_r, cicewp, cliqwp
 real, allocatable, dimension(:,:)     :: reice, reliq, zl_r
 real, allocatable, dimension(:)       :: alat

 !Trace gases
 real, allocatable, dimension(:,:)     :: co2_r, o2_r, ch4_r, n2o_r, cfc11_r, cfc12_r, cfc22_r, ccl4_r
 real, allocatable, dimension(:,:,:)   :: taucld, tauaer

 !Model levels separating low and middle clouds and high and middle clouds
 integer :: lcldmh, lcldlm

 !Outputs
 real, allocatable, dimension(:,:)     :: uflx, dflx, uflxc, dflxc, duflx_dt, duflxc_dt
 real, allocatable, dimension(:,:)     :: hr, hrc
 integer, allocatable, dimension(:,:)  :: cloudflag

 ! User input
 ! ----------
 is = 1                   !Starting point i direction
 js = 1                   !Starting point j direction
 ie = 1                   !End point i direction
 je = 1                   !End point j direction

 filename = '/gpfsm/dnb31/drholdaw/Victor/IRRADTrainingData/f522_dh.trainingdata_in.lcv.20190401_0000z.nc4'

 !Open file
 eid = nf90_open(trim(filename), NF90_NOWRITE, ncid)

 ! Data dimensions
 ! ---------------
 im = (ie-is+1)           !Number of profiles in the i direction
 jm = (je-js+1)           !Number of profiles in the j direction
 lm = 72                  !Number of vertical levels

 np = im*jm               !Number of profiles
 nb_rrtmg = 16            !Number of bands used in rrtmg

 ! Configuration for the rrtmg scheme
 ! ----------------------------------
 corespernode = 1
 partition_size = 4
 icld = 4
 idrv = 1
 inflglw = 2
 iceflglw = 3
 liqflglw = 1

 ! Read inputs from the NetCDF file
 ! --------------------------------

 !dp, t, q, qi, ql, ri, rl, o3, fcld, ts, emis, lats

 allocate(  dp(im,jm,lm))
 allocate(   t(im,jm,lm))
 allocate(   q(im,jm,lm))
 allocate(  qi(im,jm,lm))
 allocate(  ql(im,jm,lm))
 allocate(  ri(im,jm,lm))
 allocate(  rl(im,jm,lm))
 allocate(  o3(im,jm,lm))
 allocate(fcld(im,jm,lm))

 allocate(  ts(im,jm))
 allocate(emis(im,jm))
 allocate(lats(im,jm))


 ! Intermediate variable transforms
 ! --------------------------------
 
 !Cloud water content and effective radii
 allocate( cwc(im,jm,lm,2))
 allocate(reff(im,jm,lm,2))

 cwc (:,:,:,kice   ) = qi
 cwc (:,:,:,kliquid) = ql
 reff(:,:,:,kice   ) = ri * 1.0e6
 reff(:,:,:,kliquid) = rl * 1.0e6

 !Pressures
 allocate(ple(im,jm,0:lm))
 allocate(pl (im,jm,lm))

 ple(:,:,0) = 1
 do k=1,lm
   do j=js,je
     do i=is,ie
        ple(i,j,k) = ple(i,j,k-1) + dp(i,j,k) 
     enddo
   enddo
 enddo
 pl = 0.5*(ple(:,:,1:lm)+ple(:,:,0:lm-1))

 !Temperature at levels
 allocate(tlev(im,jm,lm+1))

 do k = 2, lm
    tlev(:,:,k)=  ( t(:,:,k-1)*dp(:,:,k) + t(:,:,k)*dp(:,:,k-1) ) / (dp(:,:,k-1) + dp(:,:,k))
 enddo

 tlev(:,:,lm+1) = t(:,:,lm)*(0.5*(1.0 + ple(:,:,lm-1)/ple(:,:,lm)))**(-mapl_kappa) !t2m adiabatic
 tlev(:,:,   1) = tlev(:,:,2)


 ! Main inputs for the scheme
 ! --------------------------
 allocate(pl_r(np,lm))
 allocate(ple_r(np,0:lm))
 allocate(t_r(np,lm))
 allocate(tlev_r(np,0:lm))
 allocate(tsfc(np))
 allocate(q_r(np,lm))
 allocate(o3_r(np,lm))
 allocate(emiss(np,nb_rrtmg))
 allocate(fcld_r(np,lm))
 allocate(cicewp(np,lm))
 allocate(cliqwp(np,lm))
 allocate(reice(np,lm))
 allocate(reliq(np,lm))
 allocate(zl_r(np,lm))
 allocate(alat(np))

 ! Loop to fill main inputs (note that rrtmg wants inputs with index 1 at the surface)
 ij = 0
 do j=js,je
   do i=is,ie

     !Horizontal counter
     ij = ij + 1

     !2D fields
     tsfc(ij)    = ts(i,j)
     emiss(ij,:) = emis(i,j)
     alat(ij)    = lats(i,j)

     do k = 1,lm

       lv = lm-k+1 !Flip levels

       PL_R(IJ,k) = PL(I,J,LV)/100.
       PLE_R(IJ,k) = PLE(I,J,LV)/100.

       T_R (IJ,k) = T(I,J,LV)

       xx = 1.02*100*DP(i,j,LV)
       cicewp(ij,k) = xx*cwc(i,j,lv,kice)
       cliqwp(ij,k) = xx*cwc(i,j,lv,kliquid)
       reice (ij,k) =   reff(i,j,lv,kice)
       reliq (ij,k) =   reff(i,j,lv,kliquid)

       if (liqflglw.eq.0) then
          reliq(ij,k) = min(max(reliq(ij,k),5.0),10.0)
       elseif (liqflglw.eq.1) then
          reliq(ij,k) = min(max(reliq(ij,k),2.5),60.0)
       endif

       if (iceflglw.eq.0) then
          reice(ij,k) = min(max(reice(ij,k),10.0),30.0)
       elseif (iceflglw.eq.1) then
          reice(ij,k) = min(max(reice(ij,k),13.0),130.0)
       elseif (iceflglw.eq.2) then
          reice(ij,k) = min(max(reice(ij,k), 5.0),131.0)
       elseif (iceflglw.eq.3) then
          reice(ij,k) = min(max(reice(ij,k), 5.0),140.0)
       elseif (iceflglw.eq.4) then
          reice(ij,k) = min(max(reice(ij,k)*2., 1.0),200.0)
       endif

       ple_r  (ij,k-1) = ple(i,j,lv)/100.
       tlev_r (ij,k-1) = tlev(i,j,lv+1)
       pl_r   (ij,k) = pl(i,j,lv)/100.
       t_r    (ij,k) = t(i,j,lv)
       q_r    (ij,k) = q(i,j,lv) / (1.-q(i,j,lv)) * (mapl_airmw/mapl_h2omw)
       o3_r   (ij,k) = o3(i,j,lv) * (mapl_airmw/mapl_o3mw)
       fcld_r (ij,k) = fcld(i,j,lv)

     enddo

     ple_r (ij,lm) = ple(i,j,0)/100.
     tlev_r(ij,lm) = tlev(i,j,1)

     zl_r(ij,1) = 0.
     do k=2,lm
        zl_r(ij,k) = zl_r(ij,k-1)+mapl_rgas*tlev_r(ij,k)/mapl_grav*(pl_r(ij,k-1)-pl_r(ij,k))/ple_r(ij,k)
     enddo

   enddo
 enddo

 ! Trace gases
 ! -----------
 allocate(co2_r(np,lm))
 allocate(ch4_r(np,lm))
 allocate(n2o_r(np,lm))
 allocate(o2_r(np,lm))
 allocate(cfc11_r(np,lm))
 allocate(cfc12_r(np,lm))
 allocate(cfc22_r(np,lm))
 allocate(ccl4_r(np,lm))
 allocate(taucld(np,nb_rrtmg,lm))
 allocate(tauaer(np,lm,nb_rrtmg))

 co2_r   = 0.0
 ch4_r   = 0.0
 n2o_r   = 0.0
 o2_r    = 0.0
 cfc11_r = 0.0
 cfc12_r = 0.0
 cfc22_r = 0.0
 ccl4_r  = 0.0
 taucld  = 0.0
 tauaer  = 0.0

 ! Model levels separating low and middle clouds and high and middle clouds
 ! ------------------------------------------------------------------------
 lcldlm = 50
 lcldmh = 50

 ! Outputs
 ! -------
 allocate(uflx(np,lm+1))
 allocate(dflx(np,lm+1))
 allocate(uflxc(np,lm+1))
 allocate(dflxc(np,lm+1))
 allocate(duflx_dt(np,lm+1))
 allocate(duflxc_dt(np,lm+1))
 allocate(hr(np,lm+1))
 allocate(hrc(np,lm+1))
 allocate(cloudflag(np,4))

 !Initialize the RRTMG scheme
 ! --------------------------
 call RRTMG_LW_INI(1.004e3)

 ! Main call to the RRTMG scheme
 ! -----------------------------
 call RRTMG_LW ( np,            &  ! Number of atmospheric profiles
                 lm,            &  ! Number of vertical levels
                 icld,          &  ! Cloud overlap method (hard wired to 4)
                 idrv,          &  ! Flag for calculation of dFdT, the change 
                 PL_R,          &  !
                 PLE_R,         &  !
                 T_R,           &  !
                 TLEV_R,        &  !
                 TSFC,          &  !
                 Q_R,           &  !
                 O3_R,          &  !
                 CO2_R,         &  ! Carbon dioxide (0.0 for data assimilation)
                 CH4_R,         &  ! Methane (0.0 for data assimilation)
                 N2O_R,         &  ! Nitrous oxide (0.0 for data assimilation)
                 O2_R,          &  ! Oxygen (0.0 for data assimilation)
                 CFC11_R,       &  ! Trichlorofluoromethane (0.0 for data assimilation)
                 CFC12_R,       &  ! Dichlorodifluoromethane (0.0 for data assimilation)
                 CFC22_R,       &  ! Chlorodifluoromethane (0.0 for data assimilation)
                 CCL4_R,        &  ! Carbon tetrachloride (0.0 for data assimilation)
                 EMISS,         &  !
                 INFLGLW,       &  ! Flag for cloud optical properties
                 ICEFLGLW,      &  ! Flag for ice particle specification
                 LIQFLGLW,      &  ! Flag for liquid droplet specification
                 FCLD_R,        &  ! 
                 TAUCLD,        &  ! In-cloud optical depth
                 CICEWP,        &  ! 
                 CLIQWP,        &  ! 
                 REICE,         &  ! 
                 RELIQ,         &  ! 
                 TAUAER,        &  ! Optical depth
                 ZL_R,          &  ! 
                 LCLDLM,        &  ! 
                 LCLDMH,        &  ! 
                 UFLX,          &  ! Total sky longwave upward flux (W/m2)
                 DFLX,          &  ! Total sky longwave downward flux (W/m2)
                 HR,            &  ! Total sky longwave radiative heating rate (K/d)
                 UFLXC,         &  ! Clear sky longwave upward flux (W/m2)
                 DFLXC,         &  ! Clear sky longwave downward flux (W/m2)
                 HRC,           &  ! Clear sky longwave radiative heating rate (K/d)
                 DUFLX_DT,      &  ! Change in upward longwave flux (w/m2/K)
                 DUFLXC_DT,     &  ! Change in clear sky upward longwave flux (w/m2/K)
                 CLOUDFLAG,     &  ! Optional output of cloudflag
                 DOY,           &  ! 
                 ALAT,          &  ! 
                 corespernode,  &  ! Number of procs per node
                 partition_size )  ! Partition size (hardwired to 4)

 ! Deallocate memory
 ! -----------------
 deallocate(dp, t, q, qi, ql, ri, rl, o3, fcld)
 deallocate(ts, emis, lats)

 deallocate(reff, cwc)
 deallocate(ple, pl, tlev)

 deallocate(pl_r, ple_r, t_r, tlev_r)
 deallocate(tsfc)
 deallocate(q_r, o3_r, emiss, fcld_r, cicewp, cliqwp)
 deallocate(reice, reliq, zl_r)
 deallocate(alat)

 deallocate(co2_r, ch4_r, n2o_r, cfc11_r, cfc12_r, cfc22_r, ccl4_r)
 deallocate(taucld, tauaer)

 deallocate(uflx, dflx, uflxc, dflxc, duflx_dt, duflxc_dt)
 deallocate(hr, hrc)
 deallocate(cloudflag)

 contains

! ------------------------------------------------------------------------------

! subroutine nccheck(status,iam)
! 
!  implicit none
!  integer, intent ( in) :: status
!  character(len=*), optional :: iam 
!  
!  character(len=1024) :: error_descr
! 
!  if(status /= nf90_noerr) then
! 
!    error_descr = "NetCDF error, aborting ... "
! 
!    if (present(iam)) then
!      error_descr = trim(error_descr)//", "//trim(iam)
!    endif
! 
!    error_descr = trim(error_descr)//". Error code: "//trim(nf90_strerror(status))
! 
!    !call abor1_ftn(trim(error_descr))
! 
!  end if
! 
! end subroutine nccheck

! ------------------------------------------------------------------------------

end program
