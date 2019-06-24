module rtmg_lw_driver_mod

use rtmg_lw_tools_mod
use rrtmg_lw_rad, only: rrtmg_lw
use rrtmg_lw_init, only: rrtmg_lw_ini
use mapl_constantsmod

implicit none

private

public rrtmg_lw_driver

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine rrtmg_lw_driver(config,fields,fluxes)

 implicit none
 type(configuration), intent(in)    :: config !Configuration
 type(datafields),    intent(in)    :: fields !Input
 type(fluxfields),    intent(inout) :: fluxes !Outputs

 !Configuration for the rrtmg scheme
 integer :: corespernode, partition_size
 integer :: icld, idrv, inflglw, iceflglw, liqflglw
 integer :: i, j, k, ij, lv, np, nb_rrtmg

 !Intermediate variable transforms
 integer, parameter :: kice     = 1
 integer, parameter :: kliquid  = 2
 real, allocatable, dimension(:,:,:,:) :: reff, cwc
 real, allocatable, dimension(:,:,:)   :: ple

 real, allocatable, dimension(:)       :: TLEV
 real, allocatable, dimension(:)       :: DP

 !Main inputs for the scheme
 real :: xx
 real, allocatable, dimension(:,:)     :: pl_r, ple_r, t_r, tlev_r, t2m
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
 real, allocatable, dimension(:,:,:) :: flxu_int
 real, allocatable, dimension(:,:,:) :: flxd_int
 real, allocatable, dimension(:,:,:) :: flcu_int
 real, allocatable, dimension(:,:,:) :: flcd_int
 real, allocatable, dimension(:,:,:) :: dfdts
 real, allocatable, dimension(:,:,:) :: dfdtsc

 !Debugging write input/output
 logical :: write_inout = .false.


 np = config%im*config%jm               !Number of profiles
 nb_rrtmg = 16            !Number of bands used in rrtmg

 ! Configuration for the rrtmg scheme
 ! ----------------------------------
 corespernode = 24
 partition_size = 4
 icld = 4
 idrv = 1
 inflglw = 2
 iceflglw = 3
 liqflglw = 1
 lcldlm = 18
 lcldmh = 26

  ! Intermediate variable transforms
 ! --------------------------------

 !Cloud water content and effective radii
 allocate( cwc(config%im,config%jm,config%lm,2))
 allocate(reff(config%im,config%jm,config%lm,2))

 cwc (:,:,:,kice   ) = fields%qi
 cwc (:,:,:,kliquid) = fields%ql
 reff(:,:,:,kice   ) = fields%ri * 1.0e6
 reff(:,:,:,kliquid) = fields%rl * 1.0e6

 !Pressures
 allocate(ple(config%im,config%jm,0:config%lm))

 ple(:,:,0) = 1.0
 do k=1,config%lm
   ple(:,:,k) = 2.0*fields%pl(:,:,k) - ple(:,:,k-1)
 enddo

 allocate(t2m(config%im,config%jm))
 t2m = fields%t(:,:,config%lm)*(0.5*(1.0 + ple(:,:,config%lm-1)/ple(:,:,config%lm)))**(-mapl_kappa)

 ! Main inputs for the scheme
 ! --------------------------

 allocate(pl_r(np,config%lm))
 allocate(ple_r(np,0:config%lm))
 allocate(t_r(np,config%lm))
 allocate(tlev_r(np,0:config%lm))
 allocate(tsfc(np))
 allocate(q_r(np,config%lm))
 allocate(o3_r(np,config%lm))
 allocate(emiss(np,nb_rrtmg))
 allocate(fcld_r(np,config%lm))
 allocate(cicewp(np,config%lm))
 allocate(cliqwp(np,config%lm))
 allocate(reice(np,config%lm))
 allocate(reliq(np,config%lm))
 allocate(zl_r(np,config%lm))
 allocate(alat(np))

 allocate(tlev(config%lm+1))
 allocate(dp(config%lm))

 ! Loop to fill main inputs (note that rrtmg wants inputs with index 1 at the surface)
 ij = 0
 do j=config%js,config%je
   do i=config%is,config%ie

     !Horizontal counter
     ij = ij + 1

     !2D fields
     tsfc(ij)    = fields%ts(i,j)
     emiss(ij,:) = fields%emis(i,j)
     alat(ij)    = MAPL_DEGREES_TO_RADIANS*fields%lats(i,j)

     dp(1) = (ple(i,j,1)-ple(i,j,0))
     do k = 2, config%lm
        dp(k) = (ple(i,j,k)-ple(i,j,k-1) )
        tlev(k)=  (fields%t(i,j,k-1)* dp(k) + fields%t(i,j,k) * dp(k-1)) &
                 /            (dp(k-1) + dp(k))
     enddo

     tlev(config%lm+1) = t2m(i,j)
     tlev(   1) = tlev(2)


     do k = 1,config%lm

       lv = config%lm-k+1 !Flip levels

       xx = 1.02*100*DP(LV)
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
       tlev_r (ij,k-1) = tlev(lv+1)
       pl_r   (ij,k) = fields%pl(i,j,lv)/100.
       t_r    (ij,k) = fields%t(i,j,lv)
       q_r    (ij,k) = fields%q(i,j,lv) / (1.-fields%q(i,j,lv)) * (mapl_airmw/mapl_h2omw)
       o3_r   (ij,k) = fields%o3(i,j,lv) * (mapl_airmw/mapl_o3mw)
       fcld_r (ij,k) = fields%fcld(i,j,lv)

     enddo

     ple_r (ij,config%lm) = ple(i,j,0)/100.
     tlev_r(ij,config%lm) = tlev(1)

     zl_r(ij,1) = 0.
     do k=2,config%lm
        zl_r(ij,k) = zl_r(ij,k-1)+mapl_rgas*tlev_r(ij,k)/mapl_grav*(pl_r(ij,k-1)-pl_r(ij,k))/ple_r(ij,k)
     enddo

   enddo
 enddo

 ! Trace gases
 ! -----------

 allocate(co2_r(np,config%lm))
 allocate(ch4_r(np,config%lm))
 allocate(n2o_r(np,config%lm))
 allocate(o2_r(np,config%lm))
 allocate(cfc11_r(np,config%lm))
 allocate(cfc12_r(np,config%lm))
 allocate(cfc22_r(np,config%lm))
 allocate(ccl4_r(np,config%lm))
 allocate(taucld(np,nb_rrtmg,config%lm))
 allocate(tauaer(np,config%lm,nb_rrtmg))

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

 ! Outputs
 ! -------

 allocate(uflx(np,config%lm+1))
 allocate(dflx(np,config%lm+1))
 allocate(uflxc(np,config%lm+1))
 allocate(dflxc(np,config%lm+1))
 allocate(duflx_dt(np,config%lm+1))
 allocate(duflxc_dt(np,config%lm+1))
 allocate(hr(np,config%lm+1))
 allocate(hrc(np,config%lm+1))
 allocate(cloudflag(np,4))

 uflx = 0.0
 dflx = 0.0
 uflxc = 0.0
 dflxc = 0.0
 duflx_dt = 0.0
 duflxc_dt = 0.0
 hr = 0.0
 hrc = 0.0
 cloudflag = 0.0

 ! Main call to the RRTMG scheme
 ! -----------------------------

 if (write_inout) call write_input()

 call RRTMG_LW_INI(1.004e3)

 call RRTMG_LW ( np,            &  ! Number of atmospheric profiles
                 config%lm,            &  ! Number of vertical levels
                 icld,          &  ! Cloud overlap method (hard wired to 4)
                 idrv,          &  ! Flag for calculation of dFdT, the change
                 PL_R,          &  ! Pressure at level mid points (hPa)
                 PLE_R,         &  ! Pressure at level edge points (hPa)
                 T_R,           &  ! Temperature at level mid points (K)
                 TLEV_R,        &  ! Temperature at level edge points (K)
                 TSFC,          &  ! Temperature at surface (K)
                 Q_R,           &  ! Moisture volume mixing ratio
                 O3_R,          &  ! Ozone volume mixing ratio
                 CO2_R,         &  ! Carbon dioxide (0.0 for data assimilation)
                 CH4_R,         &  ! Methane (0.0 for data assimilation)
                 N2O_R,         &  ! Nitrous oxide (0.0 for data assimilation)
                 O2_R,          &  ! Oxygen (0.0 for data assimilation)
                 CFC11_R,       &  ! Trichlorofluoromethane (0.0 for data assimilation)
                 CFC12_R,       &  ! Dichlorodifluoromethane (0.0 for data assimilation)
                 CFC22_R,       &  ! Chlorodifluoromethane (0.0 for data assimilation)
                 CCL4_R,        &  ! Carbon tetrachloride (0.0 for data assimilation)
                 EMISS,         &  ! Surface emissivity
                 INFLGLW,       &  ! Flag for cloud optical properties
                 ICEFLGLW,      &  ! Flag for ice particle specification
                 LIQFLGLW,      &  ! Flag for liquid droplet specification
                 FCLD_R,        &  ! Cloud fraction
                 TAUCLD,        &  ! In-cloud optical depth
                 CICEWP,        &  ! In-cloud ice water path (g/m2)
                 CLIQWP,        &  ! In-cloud liquid water path (g/m2)
                 REICE,         &  ! Cloud ice particle effective size (microns)
                 RELIQ,         &  ! Cloud water drop effective radius (microns)
                 TAUAER,        &  ! Optical depth
                 ZL_R,          &  ! Heights of level midpoints (m)
                 LCLDLM,        &  ! cloud layer heights for cloudFlag
                 LCLDMH,        &  ! cloud layer heights for cloudFlag
                 UFLX,          &  ! Total sky longwave upward flux (W/m2)
                 DFLX,          &  ! Total sky longwave downward flux (W/m2)
                 HR,            &  ! Total sky longwave radiative heating rate (K/d)
                 UFLXC,         &  ! Clear sky longwave upward flux (W/m2)
                 DFLXC,         &  ! Clear sky longwave downward flux (W/m2)
                 HRC,           &  ! Clear sky longwave radiative heating rate (K/d)
                 DUFLX_DT,      &  ! Change in upward longwave flux (w/m2/K)
                 DUFLXC_DT,     &  ! Change in clear sky upward longwave flux (w/m2/K)
                 CLOUDFLAG,     &  ! Optional output of cloudflag
                 config%doy,    &  ! Day of the year
                 ALAT,          &  ! Latitude of the profile in radians
                 corespernode,  &  ! Number of procs per node
                 partition_size )  ! Partition size (hardwired to 4)

 if (write_inout) call write_output()

! Process outputs and write to file
! ---------------------------------

 allocate(flxu_int(config%im,config%jm,0:config%lm))
 allocate(flxd_int(config%im,config%jm,0:config%lm))

 ij = 0
 do j=config%js,config%je
   do i=config%is,config%ie
     ij = ij + 1
     do k=0,config%lm
        lv = config%lm-k+1
        flxu_int(i,j,k)     =-uflx     (ij,lv)
        flxd_int(i,j,k)     = dflx     (ij,lv)
        fluxes%dfdts(i,j,k) =-duflx_dt (ij,lv)
     enddo
   enddo
 enddo

 fluxes%flx = flxd_int + flxu_int

 ! Deallocate memory
 ! -----------------

 deallocate(reff, cwc)
 deallocate(ple, dp, tlev)

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

 deallocate(flxu_int)
 deallocate(flxd_int)

 contains


 ! Write input
 ! -----------
 subroutine write_input

 !Write the result with what GEOS produced
 open (unit = 101, file = "input_offline.txt")
 do ij=1,np
   write(101,*) 'icld', icld
   write(101,*) 'idrv', idrv
   write(101,*) 'TSFC', TSFC(ij)
   do k=1,nb_rrtmg
     write(101,*) 'EMISS', EMISS(ij,k)
   enddo
   write(101,*) 'INFLGLW', INFLGLW
   write(101,*) 'ICEFLGLW', ICEFLGLW
   write(101,*) 'LIQFLGLW', LIQFLGLW
   write(101,*) 'TAUCLD', maxval(abs(TAUCLD))
   write(101,*) 'TAUAER', maxval(abs(TAUAER))
   write(101,*) 'LCLDLM', LCLDLM
   write(101,*) 'LCLDMH', LCLDMH
   write(101,*) 'DOY', config%doy
   write(101,*) 'ALAT', ALAT
   write(101,*) 'corespernode', corespernode
   write(101,*) 'partition_size', partition_size
   do k=0,config%lm
     write(101,*) 'PLE_R', PLE_R(ij,k)
   enddo
   do k=0,config%lm
     write(101,*) 'TLEV_R', TLEV_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'PL_R', PL_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'T_R', T_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'Q_R', Q_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'O3_R', O3_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'CO2_R', CO2_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'CH4_R', CH4_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'N2O_R', N2O_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'O2_R', O2_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'CFC11_R', CFC11_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'CFC12_R', CFC12_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'CFC22_R', CFC22_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'CCL4_R', CCL4_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'FCLD_R', FCLD_R(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'CICEWP', CICEWP(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'CLIQWP', CLIQWP(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'REICE', REICE(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'RELIQ', RELIQ(ij,k)
   enddo
   do k=1,config%lm
       write(101,*) 'ZL_R', ZL_R(ij,k)
   enddo
 enddo
 close(101)

 end subroutine write_input

 ! Write output
 ! ------------
 subroutine write_output

 !Write the result with what GEOS produced
 open (unit = 101, file = "output_offline.txt")
 do ij=1,np
   do k=1,config%lm+1
     write(101,*) 'UFLX', UFLX(ij,k)
   enddo
   do k=1,config%lm+1
     write(101,*) 'DFLX', DFLX(ij,k)
   enddo
   do k=1,config%lm+1
     write(101,*) 'HR', HR(ij,k)
   enddo
   do k=1,config%lm+1
     write(101,*) 'UFLXC', UFLXC(ij,k)
   enddo
   do k=1,config%lm+1
     write(101,*) 'DFLXC', DFLXC(ij,k)
   enddo
   do k=1,config%lm+1
     write(101,*) 'HRC', HRC(ij,k)
   enddo
   do k=1,config%lm+1
     write(101,*) 'DUFLX_DT', DUFLX_DT(ij,k)
   enddo
   do k=1,config%lm+1
     write(101,*) 'DUFLXC_DT', DUFLXC_DT(ij,k)
   enddo
   do k=1,4
     write(101,*) 'CLOUDFLAG', CLOUDFLAG(ij,k)
   enddo
 enddo
 close(101)

 end subroutine write_output

end subroutine rrtmg_lw_driver

! ------------------------------------------------------------------------------

end module rtmg_lw_driver_mod

