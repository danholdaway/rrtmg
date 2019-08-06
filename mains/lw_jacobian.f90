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

 type(jacobianmat)      :: jacobian
 type(jacobianmatwrite) :: jacobianwrite

 integer :: jmjm, jactotal
 real(4), allocatable :: JacSend(:)
 real(4), allocatable :: JacRecv(:)
 integer, allocatable :: iinRecv(:)
 integer, allocatable :: jinRecv(:)

 character(len=13) :: datetime

 integer, parameter :: nvars = 6
 integer :: n, k, vark(nvars)

 real :: pertfac
 real :: pertfac_t
 real :: pertfac_pl
 real :: pertfac_q
 real :: pertfac_qi
 real :: pertfac_ql
 real :: pertfac_o3

 integer :: i,j,count

 real(4) :: start, now, elapsed, percentdone, percentdoneincrt, percentdoneprint

 integer, parameter :: cube = 180 !  Number of grid points on cube face
 integer :: grank, gsize, csize, cerrr, iorank
 logical :: ioproc
 real :: rtsize
 integer :: tsize, tile, nx, ny, npx, npy
 integer :: kk, k1, k2, r

 !--------------------------------------------------------------------------------------------------

 ! MPI
 ! ---
 call mpi_init ( cerrr )
 call mpi_comm_rank (mpi_comm_world, grank, cerrr)
 call mpi_comm_size (mpi_comm_world, gsize, cerrr)
 call mpi_barrier(mpi_comm_world, cerrr)

 csize = gsize - 1 !Compute size, -1 IO processor

 ! IO processor
 iorank = csize
 ioproc = .false.
 if (grank == iorank) ioproc = .true.

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
 configtraj%doy = 90                !Which day of the year is it?
 configtraj%lm = 72

 configtraj%is = (npx*mod(grank,nx))+1                  !Starting point i direction
 configtraj%js = (npy*mod(grank,ny))+1 + tile*cube      !End point i direction
 configtraj%ie = configtraj%is + npx - 1
 configtraj%je = configtraj%js + npy - 1


 ! Config for perturbations
 ! ------------------------
 configpert%filename_in  = configtraj%filename_in
 configpert%filename_out = configtraj%filename_out
 configpert%doy = configtraj%doy
 configpert%lm = 72


 ! Perform trajectory calculations
 ! -------------------------------

 call allocate_fields(configtraj,fields_traj)
 call allocate_fluxes(configtraj,fluxes_traj)

 ! Read fields from model output data
 if (grank==0) print*, 'Reading the trajectory'
 if (.not.ioproc) call read_fields(configtraj,fields_traj)
 if (grank==0) print*, 'Done reading the trajectory'

 ! Call radiation scheme driver traj
 if (grank==0) print*, 'Generating flux trajectory'
 if (.not.ioproc) call rrtmg_lw_driver(configtraj,fields_traj,fluxes_traj)
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

 ! Number of vertical layers for each variables
 vark = configtraj%lm


 !Allocate Jacobian structure
 jmjm = configtraj%lm*configtraj%lm
 jactotal = nvars*configtraj%lm*configtraj%lm

 if (.not.ioproc) then
   !call allocate_jacobian(Jacobian,configtraj%is,configtraj%ie,configtraj%js,configtraj%je)
   call allocate_jacobian(Jacobian,configtraj,1,1,1,1)
 else
   call allocate_jacobian(Jacobian,configtraj,1,jmjm,1,1)
 endif

 ! Filename for Jacobian output
 jacobian%filename_out = '/gpfsm/dnb31/drholdaw/Victor/Jacobians/jacobian_'//datetime//'z.nc4'

 ! Allocate loop wise pert fields
 call allocate_fields(configtraj,fields_pert)
 call allocate_fluxes(configtraj,fluxes_pert)

 ! Jacobian to be communicated
 if (.not.ioproc == 0) then
   allocate(JacSend(jactotal))
 else
   allocate(JacRecv(grank*jactotal))
   allocate(iinRecv(grank))
   allocate(jinRecv(grank))
 endif

 ! Jacobian calculations
 ! ---------------------
 if (grank==0) print*, 'Start prepare Jacobian file'
 !call write_jacobian_setup(jacobianwrite,cube,configtraj,fields_traj,jacobian,mpi_comm_world)
 if (grank==0) print*, 'Done prepare Jacobian file'

 call mpi_barrier(mpi_comm_world, cerrr)
 count = 0
 call cpu_time(start)

 percentdoneincrt = 1.0
 percentdoneprint = percentdoneincrt

 if (grank==0) print*, 'Start Jacobian loop'

 do j = configtraj%js,configtraj%js
   do i = configtraj%is,configtraj%is

     configpert%is = i                  !Starting point i direction
     configpert%js = j                  !End point i direction
     configpert%ie = i                  !One profile at a time
     configpert%je = j                  !One profile at a time

     ! Loop over variables
      do n = 1,nvars

        do k = 1,vark(n)

        if (.not.ioproc) then

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

           ! Set Jacobian to zero
           jacobian%dflxdpl    = 0.0_4
           jacobian%dflxdt     = 0.0_4
           jacobian%dflxdq     = 0.0_4
           jacobian%dflxdqi    = 0.0_4
           jacobian%dflxdql    = 0.0_4
           jacobian%dflxdo3    = 0.0_4

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

         endif

       enddo

     enddo

     if (.not.ioproc) then

       ! Gather the Jacobian on write PE
       kk = 0
       do k1 = 1,config%lm
         do k2 = 1,config%lm
           kk = kk + 1
           JacSend(kk + 0*jactotal) = jacobian%dflxdpl(i,j,k1,k2)
           JacSend(kk + 1*jactotal) = jacobian%dflxdt (i,j,k1,k2)
           JacSend(kk + 2*jactotal) = jacobian%dflxdq (i,j,k1,k2)
           JacSend(kk + 3*jactotal) = jacobian%dflxdqi(i,j,k1,k2)
           JacSend(kk + 4*jactotal) = jacobian%dflxdql(i,j,k1,k2)
           JacSend(kk + 5*jactotal) = jacobian%dflxdo3(i,j,k1,k2)
         enddo
       enddo

     endif

     ! Everybody
     call mpi_gather(JacSend, jactotal, mpi_float, JacRecv, jactotal, mpi_float, iorank, mpi_comm_world)
     call mpi_gather(i, 1, mpi_int, iinRecv, 1, mpi_int, iorank, mpi_comm_world)
     call mpi_gather(j, 1, mpi_int, jinRecv, 1, mpi_int, iorank, mpi_comm_world)

     if (iproc == 1) then

       do r = 1,grank
         kk = grank*jactotal
         do k1 = 1,config%lm
           do k2 = 1,config%lm
             kk = kk + 1
             jacobian%dflxdpl(r,1,k1,k2) = JacRecv(kk + 0*jactotal)
             jacobian%dflxdt (r,1,k1,k2) = JacRecv(kk + 1*jactotal)
             jacobian%dflxdq (r,1,k1,k2) = JacRecv(kk + 2*jactotal)
             jacobian%dflxdqi(r,1,k1,k2) = JacRecv(kk + 3*jactotal)
             jacobian%dflxdql(r,1,k1,k2) = JacRecv(kk + 4*jactotal)
             jacobian%dflxdo3(r,1,k1,k2) = JacRecv(kk + 5*jactotal)
           enddo
         enddo
       enddo

       call write_jacobian_writegathered(jacobianwrite,configpert,jacobian,csize,iinRecv,jinRecv)

     endif

     count = count + 1
     percentdone = 100.0*count/((configtraj%ie-configtraj%is+1)*(configtraj%je-configtraj%js+1))

     if (percentdone > percentdoneprint) then

       if (grank==0) print*, percentdone,'% done'
       call cpu_time(now)
       elapsed = now-start         if (grank==0) print*, '  Elapsed time:         ', (real((1.0/60.0)*elapsed,4)), 'minutes'
       if (percentdone<100.0) then
       if (grank==0) print*, '  Remainig time (est.): ', (real((1.0/60.0)*elapsed*((100.0-percentdone)/percentdone),4)), 'minutes'
       endif
       if (grank==0) print*, ' '

       percentdoneprint = percentdoneprint + percentdoneincrt

     endif

   enddo
 enddo
 if (grank==0) print*, 'Done Jacobian loop'

 ! Deallocate
 if (allocated(JacSend)) deallocate(JacSend)
 if (allocated(JacRecv)) deallocate(JacRecv)
 if (allocated(iinRecv)) deallocate(iinRecv)
 if (allocated(jinRecv)) deallocate(jinRecv)

 ! Write Jacobian to file
 ! ----------------------
 !if (.not.ioproc) then
   !if (grank==0) print*, 'Start write Jacobian to file'
   !call mpi_barrier(mpi_comm_world, cerrr)
   !call write_jacobian_write(jacobianwrite,configpert,jacobian)
   !if (grank==0) print*, 'Done write Jacobian to file'
 !endif

 if (grank==0) print*, 'Start close Jacobian file'
 call mpi_barrier(mpi_comm_world, cerrr)
 call write_jacobian_final(jacobianwrite)
 if (grank==0) print*, 'Done close Jacobian to file'

 ! Finalize everything
 ! -------------------
 if (.not.ioproc) then
   if (grank==0) print*, 'Start finalize'
   call mpi_barrier(mpi_comm_world, cerrr)
   call deallocate_jacobian(jacobian)
   call deallocate_fields(fields_pert)
   call deallocate_fluxes(fluxes_pert)
   call deallocate_fields(fields_traj)
   call deallocate_fluxes(fluxes_traj)
   call mpi_barrier(mpi_comm_world, cerrr)
   call MPI_Finalize( cerrr )
   if (grank==0) print*, 'done finalize'
 endif

end program
