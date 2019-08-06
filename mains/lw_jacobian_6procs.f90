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
 type(jacobianmat) :: jacobian
 type(jacobianmat), allocatable :: jacobianIO(:)
 type(jacobianmatwrite) :: jacobianwrite

 integer :: lmlm, jactotal
 real(4) :: latsend, lonsend
 real(4), allocatable :: JacSend(:)
 real(4), allocatable :: JacRecv(:)
 integer, allocatable :: iinRecv(:)
 integer, allocatable :: jinRecv(:)
 real(4), allocatable :: latRecv(:)
 real(4), allocatable :: lonRecv(:)

 character(len=13) :: datetime

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
 integer :: grank, gsize, csize, ierr, iorank
 logical :: ioproc
 real :: rtsize
 integer :: tsize, tile, nx, ny, npx, npy
 integer :: kk, k1, k2
 real(4) :: iostart, iofinal

 !--------------------------------------------------------------------------------------------------

 ! Domain decomposition
 ! --------------------
 call mpi_init ( ierr )
 call mpi_comm_rank (mpi_comm_world, grank, ierr)
 call mpi_comm_size (mpi_comm_world, gsize, ierr)

 csize = gsize - 1 !Compute size (-1 IO processor)

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


 ! IO rank and flag
 ! ----------------
 iorank = gsize - 1
 ioproc = .false.
 if (grank == iorank) ioproc = .true.

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
 jacobian_filename = '/gpfsm/dnb31/drholdaw/Victor/Jacobians/jacobian_'//datetime//'z.nc4'
 if (grank==0) print*, ' (IO) Start prepare Jacobian file'
 !call write_jacobian_setup(jacobianwrite,cube,configtraj,jacobian_filename,mpi_comm_world,mpi_info_null)
 if (ioproc) call write_jacobian_setup(jacobianwrite,cube,6*cube,configtraj,jacobian_filename)
 if (grank==0) print*, ' (IO) Done prepare Jacobian file'


 if (.not.ioproc) then

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

 endif

 ! Allocate Jacobian structure
 ! ---------------------------
 if (.not.ioproc) then
   call allocate_jacobian(jacobian,configtraj)
 else
   allocate(jacobianIO(csize))
   do n = 1,csize
     call allocate_jacobian(jacobianIO(n),configtraj)
   enddo
 endif

! Jacobian to be communicated
! ---------------------------
lmlm = configtraj%lm*configtraj%lm
jactotal = nvars*lmlm

 allocate(JacSend(jactotal))
 if (ioproc) then
   allocate(JacRecv(gsize*jactotal))
   allocate(iinRecv(gsize))
   allocate(jinRecv(gsize))
   allocate(latRecv(gsize))
   allocate(lonRecv(gsize))
 endif

 ! Timing diagnostics
 ! ------------------
 if (grank==0) then

   count = 0
   call cpu_time(start)
   percentdoneincrt = 0.010
   percentdoneprint = percentdoneincrt

   print*, 'Start Jacobian loop'

 endif

 do j = configtraj%js,configtraj%je
   do i = configtraj%is,configtraj%ie

     if (.not.ioproc) then

       configpert%is = i                  !Starting point i direction
       configpert%js = j                  !End point i direction
       configpert%ie = i                  !One profile at a time
       configpert%je = j                  !One profile at a time

       ! Set Jacobian to zero
       jacobian%dflxdpl    = 0.0_4
       jacobian%dflxdt     = 0.0_4
       jacobian%dflxdq     = 0.0_4
       jacobian%dflxdqi    = 0.0_4
       jacobian%dflxdql    = 0.0_4
       jacobian%dflxdo3    = 0.0_4

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
             jacobian%dflxdpl(:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_pl, 4)
           elseif (n==2) then
             jacobian%dflxdt (:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_t , 4)
           elseif (n==3) then
             jacobian%dflxdq (:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_q , 4)
           elseif (n==4) then
             jacobian%dflxdqi(:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_qi, 4)
           elseif (n==5) then
             jacobian%dflxdql(:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_ql, 4)
           elseif (n==6) then
             jacobian%dflxdo3(:,k) = real((fluxes_pert%flx(i,j,:) - fluxes_traj%flx(i,j,:)) / pertfac_o3, 4)
           else
             if (grank==0) print*, 'Unrecognised variable when building Jacobian'
             stop
           endif

         enddo

       enddo

       ! Reshape Jacobian ready to send
       ! ------------------------------
       kk = 0
       do k1 = 1,configpert%lm
         do k2 = 1,configpert%lm
           kk = kk + 1
           JacSend(kk + 0*lmlm) = jacobian%dflxdpl(k1,k2)
           JacSend(kk + 1*lmlm) = jacobian%dflxdt (k1,k2)
           JacSend(kk + 2*lmlm) = jacobian%dflxdq (k1,k2)
           JacSend(kk + 3*lmlm) = jacobian%dflxdqi(k1,k2)
           JacSend(kk + 4*lmlm) = jacobian%dflxdql(k1,k2)
           JacSend(kk + 5*lmlm) = jacobian%dflxdo3(k1,k2)
         enddo
       enddo
       if (kk.ne.lmlm) then
         print*, 'kk not equal to jactotal when building JacSend'
         stop
       endif

     else

       JacSend = 0.0_4

     endif


     ! Gather Jacobian and grid information
     ! ------------------------------------
     if (.not.ioproc) then
       latsend = real(fields_traj%lats(i,j),4)
       lonsend = real(fields_traj%lons(i,j),4)
     else
       latsend = 0.0_4
       lonsend = 0.0_4
     endif

     call mpi_gather(JacSend, jactotal, mpi_float, JacRecv, jactotal, mpi_float, iorank, mpi_comm_world, ierr)

     call mpi_gather(i,       1, mpi_int,   iinRecv, 1, mpi_int,   iorank, mpi_comm_world, ierr)
     call mpi_gather(j,       1, mpi_int,   jinRecv, 1, mpi_int,   iorank, mpi_comm_world, ierr)
     call mpi_gather(latsend, 1, mpi_float, latRecv, 1, mpi_float, iorank, mpi_comm_world, ierr)
     call mpi_gather(lonsend, 1, mpi_float, lonRecv, 1, mpi_float, iorank, mpi_comm_world, ierr)

     ! Reshape back to Jacobian arrays
     ! -------------------------------
     if (ioproc) then

      call cpu_time(iostart)

       do n = 1,csize
         kk = (n-1)*jactotal
         do k1 = 1,configpert%lm
           do k2 = 1,configpert%lm
             kk = kk + 1
             jacobianIO(n)%dflxdpl(k1,k2) = JacRecv(kk + 0*lmlm)
             jacobianIO(n)%dflxdt (k1,k2) = JacRecv(kk + 1*lmlm)
             jacobianIO(n)%dflxdq (k1,k2) = JacRecv(kk + 2*lmlm)
             jacobianIO(n)%dflxdqi(k1,k2) = JacRecv(kk + 3*lmlm)
             jacobianIO(n)%dflxdql(k1,k2) = JacRecv(kk + 4*lmlm)
             jacobianIO(n)%dflxdo3(k1,k2) = JacRecv(kk + 5*lmlm)
           enddo
         enddo
       enddo

       call cpu_time(iofinal)
       print*, 'Time taken to reshape Jacobians = ', (real(iofinal-iostart,4)), 'seconds'
       call cpu_time(iostart)

       ! Write gathered Jacobians to disk
       ! --------------------------------
       !print*, " Start jacobian write gathered"
       !call write_jacobian_writegathered(jacobianwrite,configpert,csize,jacobianIO,&
       !                                  iinRecv(1:csize),jinRecv(1:csize),latRecv(1:csize),lonRecv(1:csize))
       !print*, " Done jacobian write gathered"

       call cpu_time(iofinal)
       print*, 'Time taken to write global Jacobians = ', (real(iofinal-iostart,4)), 'seconds'

     endif

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
       endif

       percentdoneprint = percentdoneprint + percentdoneincrt

     endif

   enddo
 enddo
 if (grank==0) print*, 'Done Jacobian loop'


 ! Let everyone catch up before finalize stages
 ! --------------------------------------------
 call mpi_barrier(mpi_comm_world, ierr)


 ! Clean up Jacobian file
 ! ----------------------
 if (grank==0) print*, ' (IO) Start close Jacobian file'
 call write_jacobian_final(jacobianwrite)
 if (grank==0) print*, ' (IO) Done close Jacobian to file'


 ! Finalize everything
 ! -------------------
 if (grank==0) print*, 'Start finalize'

 if (allocated(JacSend)) deallocate(JacSend)
 if (allocated(JacRecv)) deallocate(JacRecv)
 if (allocated(iinRecv)) deallocate(iinRecv)
 if (allocated(jinRecv)) deallocate(jinRecv)
 if (allocated(latRecv)) deallocate(latRecv)
 if (allocated(lonRecv)) deallocate(lonRecv)

 if (.not. ioproc) then
   call deallocate_fields(fields_pert)
   call deallocate_fluxes(fluxes_pert)
   call deallocate_fields(fields_traj)
   call deallocate_fluxes(fluxes_traj)
 endif

 if (.not.ioproc) then
   call deallocate_jacobian(jacobian)
 else
   do n = 1,csize
     call deallocate_jacobian(jacobianIO(n))
   enddo
   deallocate(jacobianIO)
 endif

 call MPI_Finalize( ierr )
 if (grank==0) print*, 'done finalize'

end program
