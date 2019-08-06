program lw_driver

 use rtmg_lw_tools_mod
 use rtmg_lw_driver_mod, only: rrtmg_lw_driver

 use mpi

 implicit none

 type(configuration) :: configtraj
 type(configuration) :: configpert

 type(datafields) :: fields_traj
 type(datafields) :: fields_pert
 type(fluxfields) :: fluxes_traj
 type(fluxfields) :: fluxes_pert

 character(len=2048) :: jacobian_filename
 type(jacobianmat4d) :: jacobian
 type(jacobianmatwrite) :: jacobianwrite

 character(len=13) :: datetime

 character(len=4) :: stris, strjs

 integer, parameter :: nvars = 6
 integer :: n, k

 real :: pertfac
 real :: pertfac_t
 real :: pertfac_pl
 real :: pertfac_q
 real :: pertfac_qi
 real :: pertfac_ql
 real :: pertfac_o3

 integer :: i,j,count

 real(4) :: start, now, elapsed, percentdone, percentdoneincrt, percentdoneprint

 integer, parameter :: cube = 720 !  Number of grid points on cube face
 integer :: grank, gsize, csize, ierr
 real :: rtsize
 integer :: tsize, tile, nx, ny, npx, npy


 !--------------------------------------------------------------------------------------------------


 ! Domain decomposition
 ! --------------------
 call mpi_init ( ierr )
 call mpi_comm_rank (mpi_comm_world, grank, ierr)
 call mpi_comm_size (mpi_comm_world, gsize, ierr)

 csize = gsize

 rtsize = real(csize)/6.0
 tsize = int(rtsize)

 ! Check number of procs is multiple of 6
 if ( rtsize - tsize .ne. 0.0 ) then
   if (grank==0) print*, 'Number of mpi tasks must be divisisble by 6 (number of cube faces)'
   if (grank==0) print*, rtsize, tsize, rtsize - tsize
   stop
 endif

 ! Check procs on face are squarerootable
 if ( sqrt(real(tsize)) - int(sqrt(real(tsize))) .ne. 0.0 ) then
   if (grank==0) print*, 'sqrt of mpi task / should be integer'
   stop
 endif

 ! Layout on cube face
 nx = int(sqrt(real(tsize)))
 ny = nx

 ! Check number of grid point divisible by number of procs
 if ( real(cube)/real(nx) - int(real(cube)/real(nx)) .ne. 0 ) then
   if (grank==0) print*, 'Cube size not divisible by number of tasks along face'
   stop
 endif

 npx = cube/nx
 npy = npx

 tile = grank / tsize


 ! Config for trajecotory
 ! ----------------------
 datetime = '20190401_0000'
 configtraj%filename_in  = '/gpfsm/dnb31/drholdaw/Victor/IRRADTrainingData/f522_dh.trainingdata_in.lcv.'//datetime//'z.nc4'  !Training data in
 configtraj%filename_out = '/gpfsm/dnb31/drholdaw/Victor/IRRADTrainingData/f522_dh.trainingdata_out.lcv.'//datetime//'z.nc4' !Training data out
 configtraj%doy          = 90  !Which day of the year is it?
 configtraj%lm           = 72  !Number of model levels

 configtraj%is = (npx*mod(grank,nx))+1             !Starting point i direction
 configtraj%js = (grank/ny)*npy+1 !+ tile*cube      !End point i direction

 configtraj%ie = configtraj%is + npx - 1
 configtraj%je = configtraj%js + npy - 1


 ! Config for perturbations
 ! ------------------------
 configpert%filename_in  = configtraj%filename_in
 configpert%filename_out = configtraj%filename_out
 configpert%doy          = configtraj%doy
 configpert%lm           = configtraj%lm


 ! Jacobian output file
 ! --------------------
 write(stris,'(I4.4)') configtraj%is
 write(strjs,'(I4.4)') configtraj%js
 jacobian_filename = '/gpfsm/dnb31/drholdaw/Victor/Jacobians/jacobian_'//&
                     datetime//'z_is'//stris//'_js'//strjs//'.nc4'


 ! Allocate field and flux traj and pert data
 ! ------------------------------------------
 call allocate_fields(configtraj,fields_traj)
 call allocate_fluxes(configtraj,fluxes_traj)
 call allocate_fields(configtraj,fields_pert)
 call allocate_fluxes(configtraj,fluxes_pert)


 ! Read fields from model output data
 ! ----------------------------------
 if (grank==0) print*, 'Reading the trajectory'
 call read_fields(configtraj,fields_traj)
 if (grank==0) print*, 'Done reading the trajectory'


 ! Call radiation scheme driver traj
 ! ---------------------------------
 if (grank==0) print*, 'Generating flux trajectory'
 call rrtmg_lw_driver(configtraj,fields_traj,fluxes_traj)
 if (grank==0) print*, 'Done generating flux trajectory'


 ! Prepare perturbations
 ! ---------------------
 pertfac    = 1.0
 pertfac_pl = pertfac*1.0e-6
 pertfac_t  = pertfac*1.0e-6
 pertfac_q  = pertfac*1.0e-9
 pertfac_qi = pertfac*1.0e-9
 pertfac_ql = pertfac*1.0e-9
 pertfac_o3 = pertfac*1.0e-10


 ! Allocate Jacobian structure
 ! ---------------------------
 call allocate_jacobian4d(jacobian,configtraj)


 ! Timing diagnostics
 ! ------------------
 if (grank==0) then
   count = 0
   call cpu_time(start)
   percentdoneincrt = 1.0
   percentdoneprint = percentdoneincrt
   print*, 'Start Jacobian loop'
 endif

 ! Set Jacobian to zero
 jacobian%dflxdpl    = 0.0_4
 jacobian%dflxdt     = 0.0_4
 jacobian%dflxdq     = 0.0_4
 jacobian%dflxdqi    = 0.0_4
 jacobian%dflxdql    = 0.0_4
 jacobian%dflxdo3    = 0.0_4

 do j = configtraj%js,configtraj%je
   do i = configtraj%is,configtraj%ie

     configpert%is = i                  !Starting point i direction
     configpert%js = j                  !End point i direction
     configpert%ie = i                  !One profile at a time
     configpert%je = j                  !One profile at a time

     ! Loop over variables
     do n = 1,nvars
       do k = 1,configpert%lm

         !Overwrite the inputs
         call copy_fields(configpert,fields_pert,fields_traj,i,i,j,j)
         if (n==1) then
           fields_pert%pl(i,j,k) = fields_pert%pl(i,j,k) + pertfac_pl
         elseif (n==2) then
           fields_pert%t(i,j,k)  = fields_pert%t(i,j,k)  + pertfac_t
         elseif (n==3) then
           fields_pert%q(i,j,k)  = fields_pert%q(i,j,k)  + pertfac_q
         elseif (n==4) then
           fields_pert%qi(i,j,k) = fields_pert%qi(i,j,k) + pertfac_qi
         elseif (n==5) then
           fields_pert%ql(i,j,k) = fields_pert%ql(i,j,k) + pertfac_ql
         elseif (n==6) then
           fields_pert%o3(i,j,k) = fields_pert%o3(i,j,k) + pertfac_o3
         else
           if (grank==0) print*, 'Unrecognised variable when setting perts'
           stop
         endif

         ! Call radiation scheme driver
         call rrtmg_lw_driver(configpert,fields_pert,fluxes_pert)

         if (n==1) then
           jacobian%dflxdpl(i,j,:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_pl, 4)
         elseif (n==2) then
           jacobian%dflxdt (i,j,:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_t , 4)
         elseif (n==3) then
           jacobian%dflxdq (i,j,:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_q , 4)
         elseif (n==4) then
           jacobian%dflxdqi(i,j,:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_qi, 4)
         elseif (n==5) then
           jacobian%dflxdql(i,j,:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_ql, 4)
         elseif (n==6) then
           jacobian%dflxdo3(i,j,:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_o3, 4)
         else
           if (grank==0) print*, 'Unrecognised variable when building Jacobian'
           stop
         endif

         enddo
       enddo

     ! Timing diagnostics
     ! ------------------
     if (grank == 0) then

       count = count + 1
       percentdone = 100.0*count/((configtraj%ie-configtraj%is+1)*(configtraj%je-configtraj%js+1))

       if (percentdone > percentdoneprint) then
         print*, percentdone,'% done'
         call cpu_time(now)
         elapsed = now-start
         print*, '  Elapsed time:         ', (real((1.0/60.0)*elapsed,4)), 'minutes'
         if (percentdone<100.0) then
           print*, '  Remainig time (est.): ', (real((1.0/60.0)*elapsed*((100.0-percentdone)/percentdone),4)), 'minutes'
         endif
         print*, ' '

         percentdoneprint = percentdoneprint + percentdoneincrt

       endif


     endif

   enddo
 enddo

 if (grank==0) print*, 'Done Jacobian loop'

 ! Clean up Jacobian file
 ! ----------------------
 if (grank==0) print*, ' (IO) Start write Jacobian to file'

 call write_jacobian_setup(jacobianwrite,configtraj%ie-configtraj%is+1,configtraj%je-configtraj%js+1,&
                           configtraj,jacobian_filename)
 call write_jacobian_variables(jacobianwrite,configtraj,fields_traj,jacobian)
 call write_jacobian_final(jacobianwrite)

 if (grank==0) print*, ' (IO) Done write Jacobian to file'

 ! Finalize everything
 ! -------------------
 if (grank==0) print*, 'Start finalize'

 call deallocate_fields(fields_pert)
 call deallocate_fluxes(fluxes_pert)
 call deallocate_fields(fields_traj)
 call deallocate_fluxes(fluxes_traj)

 call deallocate_jacobian4d(jacobian)

 call MPI_Finalize( ierr )
 if (grank==0) print*, 'done finalize'

end program
