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
 logical :: write_inout

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
 character(len=2048) :: filename_in
 character(len=2048) :: filename_out
 character(len=20) :: field
 real, allocatable, dimension(:,:,:)   :: pl, t, q, qi, ql, ri, rl, o3, fcld
 real, allocatable, dimension(:,:)     :: ts, emis, lats, lons

 real, allocatable, dimension(:,:,:)   :: flx_geos
 real, allocatable, dimension(:,:,:)   :: flxu_geos
 real, allocatable, dimension(:,:,:)   :: flxd_geos
 real, allocatable, dimension(:,:,:)   :: flcu_geos
 real, allocatable, dimension(:,:,:)   :: flcd_geos
 real, allocatable, dimension(:,:,:)   :: dfdts_geos
 real, allocatable, dimension(:,:,:)   :: dfdtsc_geos

 integer :: ncid, varid
 integer, allocatable, dimension(:)    :: istart3, icount3
 integer, allocatable, dimension(:)    :: istart2, icount2

 !Intermediate variable transforms
 integer, parameter :: kice     = 1
 integer, parameter :: kliquid  = 2
 real, allocatable, dimension(:,:,:,:) :: reff, cwc
 real, allocatable, dimension(:,:,:)   :: ple, dp, tlev

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
 real, allocatable, dimension(:,:,:)   :: flxu_int,flxd_int,dfdts
 real, allocatable, dimension(:,:,:)   :: flcu_int,flcd_int,dfdtsc
 real, allocatable, dimension(:,:,:)   :: flx_int
 real, allocatable, dimension(:,:,:)   :: flc_int

 print*, ' '
 print*, 'Standalone RRTMG scheme'
 print*, ' '

 ! User input
 ! ----------
 is = 1                   !Starting point i direction
 js = 1                   !Starting point j direction
 ie = 1                   !End point i direction
 je = 1                   !End point j direction

 filename_in  = '/gpfsm/dnb31/drholdaw/wrk.f522dh/f522_dh.trainingdata_in.lcv.20190401_0000z.nc4'  !Training data in
 filename_out = '/gpfsm/dnb31/drholdaw/wrk.f522dh/f522_dh.trainingdata_out.lcv.20190401_0000z.nc4' !Training data out

 doy = 90                 !Which day of the year is it?

 write_inout = .false.     !Option to write out the input and output of the scheme (for debugging)

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
 lcldlm = 18
 lcldmh = 26

 ! Read inputs/outputs from the NetCDF files
 ! -----------------------------------------

 print*, 'Read inputs/outputs from the NetCDF files'
 print*, ' '

 allocate(  pl(im,jm,lm))
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
 allocate(lons(im,jm))

 allocate(flx_geos(im,jm,0:lm))
 allocate(flxu_geos(im,jm,0:lm))
 allocate(flxd_geos(im,jm,0:lm))
 allocate(flcu_geos(im,jm,0:lm))
 allocate(flcd_geos(im,jm,0:lm))
 allocate(dfdts_geos(im,jm,0:lm))
 allocate(dfdtsc_geos(im,jm,0:lm))

 !Read ranges
 allocate(istart3(4), icount3(4))
 allocate(istart2(3), icount2(3))

 istart3(1) = is; icount3(1) = ie - is + 1
 istart3(2) = js; icount3(2) = je - js + 1
 istart3(3) = 1 ; icount3(3) = lm
 istart3(4) = 1 ; icount3(4) = 1

 istart2(1) = is; icount2(1) = ie - is + 1
 istart2(2) = js; icount2(2) = je - js + 1
 istart2(3) = 1 ; icount2(3) = 1

 print*, '   input file: ', trim(filename_in)
 print*, ' '

 !Open input file
 call nccheck ( nf90_open(trim(filename_in), NF90_NOWRITE, ncid), "nf90_open"//trim(filename_in) )

 field = "pl"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, pl, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "t"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, t, istart3, icount3 ),    "nf90_get_var "//trim(field) )

 field = "q"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, q, istart3, icount3 ),    "nf90_get_var "//trim(field) )

 field = "qi"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, qi, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "ql"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, ql, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "ri"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, ri, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "rl"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, rl, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "o3"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, o3, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "fcld"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fcld, istart3, icount3 ), "nf90_get_var "//trim(field) )

 field = "ts"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, ts, istart2, icount2 ),   "nf90_get_var "//trim(field) )

 field = "emis"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, emis, istart2, icount2 ), "nf90_get_var "//trim(field) )

 field = "lats"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, lats, istart2, icount2 ), "nf90_get_var "//trim(field) )

 field = "lons"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, lons, istart2, icount2 ), "nf90_get_var "//trim(field) )

 call nccheck ( nf90_close(ncid), "nf90_close" )


 !Open output file
 print*, '   output file: ', trim(filename_out)
 print*, ' '

 call nccheck ( nf90_open(trim(filename_out), NF90_NOWRITE, ncid), "nf90_open"//trim(filename_out) )

 icount3(3) = lm+1

 field = "flx"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, flx_geos, istart3, icount3 ),    "nf90_get_var "//trim(field) )

 field = "flxu"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, flxu_geos, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "flxd"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, flxd_geos, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "flcu"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, flcu_geos, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "flcd"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, flcd_geos, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "dfdts"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, dfdts_geos, istart3, icount3 ),  "nf90_get_var "//trim(field) )

 field = "dfdtsc"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, dfdtsc_geos, istart3, icount3 ), "nf90_get_var "//trim(field) )

 call nccheck ( nf90_close(ncid), "nf90_close" )

 deallocate(istart3, istart2)
 deallocate(icount3, icount2)


 ! Intermediate variable transforms
 ! --------------------------------

 print*, 'Intermediate variable transforms'
 print*, ' '

 !Cloud water content and effective radii
 allocate( cwc(im,jm,lm,2))
 allocate(reff(im,jm,lm,2))

 cwc (:,:,:,kice   ) = qi
 cwc (:,:,:,kliquid) = ql
 reff(:,:,:,kice   ) = ri * 1.0e6
 reff(:,:,:,kliquid) = rl * 1.0e6

 !Pressures
 allocate(ple(im,jm,0:lm))
 allocate(dp (im,jm,lm))

 ple(:,:,0) = 1
 do k=1,lm
   ple(:,:,k) = 2*pl(:,:,k) - ple(:,:,k-1)
   dp(:,:,k) = ple(:,:,k)-ple(:,:,k-1)
 enddo

 !Temperature at levels
 allocate(tlev(im,jm,lm+1))

 do k = 2, lm
    tlev(:,:,k)=  ( t(:,:,k-1)*dp(:,:,k) + t(:,:,k)*dp(:,:,k-1) ) / (dp(:,:,k-1) + dp(:,:,k))
 enddo

 tlev(:,:,lm+1) = t(:,:,lm)*(0.5*(1.0 + ple(:,:,lm-1)/ple(:,:,lm)))**(-mapl_kappa) !t2m adiabatic
 tlev(:,:,   1) = tlev(:,:,2)


 ! Main inputs for the scheme
 ! --------------------------

 print*, 'Main inputs for the scheme'
 print*, ' '

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
     alat(ij)    = 0.0174533*lats(i,j)

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

 print*, 'Trace gases'
 print*, ' '

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

 ! Outputs
 ! -------

 print*, 'Prepare outputs'
 print*, ' '

 allocate(uflx(np,lm+1))
 allocate(dflx(np,lm+1))
 allocate(uflxc(np,lm+1))
 allocate(dflxc(np,lm+1))
 allocate(duflx_dt(np,lm+1))
 allocate(duflxc_dt(np,lm+1))
 allocate(hr(np,lm+1))
 allocate(hrc(np,lm+1))
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

 print*, 'Main call to the RRTMG scheme'
 print*, ' '

 if (write_inout) call write_input()

 call RRTMG_LW_INI(1.004e3)

 call RRTMG_LW ( np,            &  ! Number of atmospheric profiles
                 lm,            &  ! Number of vertical levels
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
                 DOY,           &  ! Day of the year
                 ALAT,          &  ! Latitude of the profile in radians
                 corespernode,  &  ! Number of procs per node
                 partition_size )  ! Partition size (hardwired to 4)

 if (write_inout) call write_output()

! Process outputs and write to file
! ---------------------------------

 print*, 'Process outputs and write to file'
 print*, ' '

 allocate(flxu_int(im,jm,0:lm))
 allocate(flxd_int(im,jm,0:lm))
 allocate(flcu_int(im,jm,0:lm))
 allocate(flcd_int(im,jm,0:lm))
 allocate(dfdts   (im,jm,0:lm))
 allocate(dfdtsc  (im,jm,0:lm))

 ij = 0
 do j=js,je
   do i=is,ie
     ij = ij + 1
     do k=0,lm
        lv = lm-k+1
        flxu_int(i,j,k) =-uflx     (ij,lv)
        flxd_int(i,j,k) = dflx     (ij,lv)
        flcu_int(i,j,k) =-uflx     (ij,lv)
        flcd_int(i,j,k) = dflx     (ij,lv)
        dfdts   (i,j,k) =-duflx_dt (ij,lv)
        dfdtsc  (i,j,k) =-duflxc_dt(ij,lv)
     enddo
   enddo
 enddo

 allocate(flx_int(im,jm,0:lm))
 allocate(flc_int(im,jm,0:lm))

 flx_int = flxu_int + flxd_int
 flc_int = flcu_int + flcd_int

 !Write the result with what GEOS produced
 open (unit = 101, file = "output.txt")
 do j=js,je
   do i=is,ie
     write(101,*) 'Profile print versus GEOS'
     write(101,*) 'Latitude', lats(i,j), ',  longitude', lons(i,j)
     write(101,*) ' '
     write(101,*) ' flxu '
     do k=0,lm
       write(101,*) flxu_geos(i,j,k), flxu_int(i,j,k)
     enddo
     write(101,*) ' flxd '
     do k=0,lm
       write(101,*) flxd_geos(i,j,k), flxd_int(i,j,k)
     enddo
     write(101,*) ' dfdts '
     do k=0,lm
       write(101,*) dfdts_geos(i,j,k), dfdts(i,j,k)
     enddo
   enddo
 enddo
 close(101)

 ! Deallocate memory
 ! -----------------
 deallocate(pl, t, q, qi, ql, ri, rl, o3, fcld)
 deallocate(ts, emis, lats, lons)

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

 deallocate(flxu_int,flxd_int,dfdts)
 deallocate(flcu_int,flcd_int,dfdtsc)

 deallocate(flx_int)
 deallocate(flc_int)

 contains

! ------------------------------------------------------------------------------

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
   write(101,*) 'DOY', DOY
   write(101,*) 'ALAT', ALAT
   write(101,*) 'corespernode', corespernode
   write(101,*) 'partition_size', partition_size
   do k=0,lm
     write(101,*) 'PLE_R', PLE_R(ij,k)
   enddo
   do k=0,lm
     write(101,*) 'TLEV_R', TLEV_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'PL_R', PL_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'T_R', T_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'Q_R', Q_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'O3_R', O3_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'CO2_R', CO2_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'CH4_R', CH4_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'N2O_R', N2O_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'O2_R', O2_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'CFC11_R', CFC11_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'CFC12_R', CFC12_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'CFC22_R', CFC22_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'CCL4_R', CCL4_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'FCLD_R', FCLD_R(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'CICEWP', CICEWP(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'CLIQWP', CLIQWP(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'REICE', REICE(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'RELIQ', RELIQ(ij,k)
   enddo
   do k=1,lm
       write(101,*) 'ZL_R', ZL_R(ij,k)
   enddo
 enddo
 close(101)

end subroutine write_input

! ------------------------------------------------------------------------------

subroutine write_output

 !Write the result with what GEOS produced
 open (unit = 101, file = "output_offline.txt")
 do ij=1,np
   do k=1,lm+1
     write(101,*) 'UFLX', UFLX(ij,k)
   enddo
   do k=1,lm+1
     write(101,*) 'DFLX', DFLX(ij,k)
   enddo
   do k=1,lm+1
     write(101,*) 'HR', HR(ij,k)
   enddo
   do k=1,lm+1
     write(101,*) 'UFLXC', UFLXC(ij,k)
   enddo
   do k=1,lm+1
     write(101,*) 'DFLXC', DFLXC(ij,k)
   enddo
   do k=1,lm+1
     write(101,*) 'HRC', HRC(ij,k)
   enddo
   do k=1,lm+1
     write(101,*) 'DUFLX_DT', DUFLX_DT(ij,k)
   enddo
   do k=1,lm+1
     write(101,*) 'DUFLXC_DT', DUFLXC_DT(ij,k)
   enddo
   do k=1,4
     write(101,*) 'CLOUDFLAG', CLOUDFLAG(ij,k)
   enddo
 enddo
 close(101)

end subroutine write_output

! ------------------------------------------------------------------------------

subroutine nccheck(status,iam)

  implicit none
  integer, intent ( in) :: status
  character(len=*), optional :: iam

  character(len=1024) :: error_descr

  if(status /= nf90_noerr) then

    error_descr = "NetCDF error, aborting ... "

    if (present(iam)) then
      error_descr = trim(error_descr)//", "//trim(iam)
    endif

    error_descr = trim(error_descr)//". Error code: "//trim(nf90_strerror(status))

    print*, "Aborting: ", trim(error_descr)

    call abort()

  end if

 end subroutine nccheck

! ------------------------------------------------------------------------------

end program
