program lw_driver

 use rtmg_lw_tools_mod
 use rtmg_lw_driver_mod, only: rrtmg_lw_driver

 implicit none

 type(configuration) :: config

 type(datafields) :: fields_traj
 type(datafields) :: fields_pert
 type(fluxfields) :: fluxes_traj
 type(fluxfields) :: fluxes_pert

 type(jacobianmat) :: jacobian

 integer, parameter :: nvars = 11
 integer :: n, k, vark(nvars)
 real :: pert, pertfac = 1000.0
 character(len=4) :: isc, jsc

 ! Some user input for the data to be used
 ! ---------------------------------------
 config%is = 1                  !Starting point i direction
 config%js = 1                  !End point i direction
 config%filename_in  = '/gpfsm/dnb31/drholdaw/Victor/IRRADTrainingData/f522_dh.trainingdata_in.lcv.20190401_0000z.nc4'  !Training data in
 config%filename_out = '/gpfsm/dnb31/drholdaw/Victor/IRRADTrainingData/f522_dh.trainingdata_out.lcv.20190401_0000z.nc4' !Training data out
 config%doy = 90                !Which day of the year is it?

 config%ie = config%is          !One profile at a time
 config%je = config%js          !One profile at a time
 config%im = (config%ie-config%is+1)
 config%jm = (config%je-config%js+1)
 config%lm = 72

 write(isc,"(I0.4)") config%is
 write(jsc,"(I0.4)") config%js
 jacobian%filename_out = '/gpfsm/dnb31/drholdaw/Victor/Jacobians/jacobian_i'//trim(isc)//'_j'//trim(jsc)//'.nc4'

 ! Run component
 ! -------------

 ! Allocate memory
 call allocate_fields(config,fields_traj)
 call allocate_fields(config,fields_pert)
 call allocate_fluxes(config,fluxes_traj)
 call allocate_fluxes(config,fluxes_pert)
 call allocate_jacobian(config,jacobian)

 ! Read fields from model output data
 call read_fields(config,fields_traj)

 ! Number of vertical layers for each variables
 vark = config%lm
 vark(10) = 1 !Ts
 vark(11) = 1 !Emis

 ! Call radiation scheme driver traj
 call rrtmg_lw_driver(config,fields_traj,fluxes_traj)

 ! Zero out the Jacobians
 jacobian%dflxdpl    = 0.0
 jacobian%dfdtsdpl   = 0.0
 jacobian%dflxdt     = 0.0
 jacobian%dfdtsdt    = 0.0
 jacobian%dflxdq     = 0.0
 jacobian%dfdtsdq    = 0.0
 jacobian%dflxdqi    = 0.0
 jacobian%dfdtsdqi   = 0.0
 jacobian%dflxdql    = 0.0
 jacobian%dfdtsdql   = 0.0
 jacobian%dflxdri    = 0.0
 jacobian%dfdtsdri   = 0.0
 jacobian%dflxdrl    = 0.0
 jacobian%dfdtsdrl   = 0.0
 jacobian%dflxdo3    = 0.0
 jacobian%dfdtsdo3   = 0.0
 jacobian%dflxdfcld  = 0.0
 jacobian%dfdtsdfcld = 0.0
 jacobian%dflxdts    = 0.0
 jacobian%dfdtsdts   = 0.0
 jacobian%dflxdemis  = 0.0
 jacobian%dfdtsdemis = 0.0

 ! Loop over variables
 do n = 1,nvars

   do k = 1,vark(n)

     !Overwrite the inputs
     call copy_fields(fields_pert,fields_traj)

     if (n==1) then
       pert = fields_pert%pl(config%is,config%js,k)/pertfac
       fields_pert%pl(config%is,config%js,k) = fields_pert%pl(config%is,config%js,k) + pert
     elseif (n==2) then
       pert = fields_pert%t(config%is,config%js,k)/pertfac
       fields_pert%t(config%is,config%js,k) = fields_pert%t(config%is,config%js,k) + pert
     elseif (n==3) then
       pert = fields_pert%q(config%is,config%js,k)/pertfac
       fields_pert%q(config%is,config%js,k) = fields_pert%q(config%is,config%js,k) + pert
     elseif (n==4) then
       pert = fields_pert%qi(config%is,config%js,k)/pertfac
       fields_pert%qi(config%is,config%js,k) = fields_pert%qi(config%is,config%js,k) + pert
     elseif (n==5) then
       pert = fields_pert%ql(config%is,config%js,k)/pertfac
       fields_pert%ql(config%is,config%js,k) = fields_pert%ql(config%is,config%js,k) + pert
     elseif (n==6) then
       pert = fields_pert%ri(config%is,config%js,k)/pertfac
       fields_pert%ri(config%is,config%js,k) = fields_pert%ri(config%is,config%js,k) + pert
     elseif (n==7) then
       pert = fields_pert%rl(config%is,config%js,k)/pertfac
       fields_pert%rl(config%is,config%js,k) = fields_pert%rl(config%is,config%js,k) + pert
     elseif (n==8) then
       pert = fields_pert%o3(config%is,config%js,k)/pertfac
       fields_pert%o3(config%is,config%js,k) = fields_pert%o3(config%is,config%js,k) + pert
     elseif (n==9) then
       pert = fields_pert%fcld(config%is,config%js,k)/pertfac
       fields_pert%fcld(config%is,config%js,k) = fields_pert%fcld(config%is,config%js,k) + pert
     elseif (n==10) then
       pert = fields_pert%ts(config%is,config%js)/pertfac
       fields_pert%ts(config%is,config%js) = fields_pert%ts(config%is,config%js) + pert
     elseif (n==11) then
       pert = fields_pert%emis(config%is,config%js)/pertfac
       fields_pert%emis(config%is,config%js) = fields_pert%emis(config%is,config%js) + pert
     endif

     ! Call radiation scheme driver
     call rrtmg_lw_driver(config,fields_pert,fluxes_pert)

     if (pert>0.0) then

       if (n==1) then
         jacobian%dflxdpl (:,k) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdpl(:,k) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==2) then
         jacobian%dflxdt (:,k) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdt(:,k) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==3) then
         jacobian%dflxdq (:,k) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdq(:,k) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==4) then
         jacobian%dflxdqi (:,k) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdqi(:,k) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==5) then
         jacobian%dflxdql (:,k) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdql(:,k) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==6) then
         jacobian%dflxdri (:,k) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdri(:,k) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==7) then
         jacobian%dflxdrl (:,k) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdrl(:,k) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==8) then
         jacobian%dflxdo3 (:,k) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdo3(:,k) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==9) then
         jacobian%dflxdfcld (:,k) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdfcld(:,k) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==10) then
         jacobian%dflxdts (:) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdts(:) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       elseif (n==11) then
         jacobian%dflxdemis (:) = (fluxes_pert%flx(config%is,config%js,:)   - fluxes_traj%flx(config%is,config%js,:)  )/(pert)
         jacobian%dfdtsdemis(:) = (fluxes_pert%dfdts(config%is,config%js,:) - fluxes_traj%dfdts(config%is,config%js,:))/(pert)
       endif

     endif

   enddo

 enddo

 ! Write the output comparison
 call write_jacobian(config,fields_traj,jacobian)

 ! Deallocate memory
 call deallocate_fields(fields_traj)
 call deallocate_fields(fields_pert)
 call deallocate_fluxes(fluxes_traj)
 call deallocate_fluxes(fluxes_pert)
 call deallocate_jacobian(jacobian)

end program
