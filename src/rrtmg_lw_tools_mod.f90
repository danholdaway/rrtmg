module rtmg_lw_tools_mod

use netcdf

implicit none

public

! Configuration settings
type configuration
 character(len=2048) :: filename_in
 character(len=2048) :: filename_out
 integer :: is,ie,js,je,lm
 integer :: doy
end type configuration

! Fields that are read from file (training data)
type datafields
 real, allocatable, dimension(:,:,:) :: pl
 real, allocatable, dimension(:,:,:) :: t
 real, allocatable, dimension(:,:,:) :: q
 real, allocatable, dimension(:,:,:) :: qi
 real, allocatable, dimension(:,:,:) :: ql
 real, allocatable, dimension(:,:,:) :: ri
 real, allocatable, dimension(:,:,:) :: rl
 real, allocatable, dimension(:,:,:) :: o3
 real, allocatable, dimension(:,:,:) :: fcld
 real, allocatable, dimension(:,:)   :: ts
 real, allocatable, dimension(:,:)   :: emis
 real, allocatable, dimension(:,:)   :: lats
 real, allocatable, dimension(:,:)   :: lons
 real, allocatable, dimension(:,:,:) :: flx
 real, allocatable, dimension(:,:,:) :: dfdts
end type datafields

! Fields that are computed by the scheme
type fluxfields
 real, allocatable, dimension(:,:,:) :: flx
 real, allocatable, dimension(:,:,:) :: dfdts
end type fluxfields

! Jacobian matrix
type jacobianmat
 real(4), allocatable, dimension(:,:) :: dflxdpl
 real(4), allocatable, dimension(:,:) :: dflxdt
 real(4), allocatable, dimension(:,:) :: dflxdq
 real(4), allocatable, dimension(:,:) :: dflxdqi
 real(4), allocatable, dimension(:,:) :: dflxdql
 real(4), allocatable, dimension(:,:) :: dflxdo3
end type jacobianmat

! Jacobian matrix
type jacobianmat4d
 real(4), allocatable, dimension(:,:,:,:) :: dflxdpl
 real(4), allocatable, dimension(:,:,:,:) :: dflxdt
 real(4), allocatable, dimension(:,:,:,:) :: dflxdq
 real(4), allocatable, dimension(:,:,:,:) :: dflxdqi
 real(4), allocatable, dimension(:,:,:,:) :: dflxdql
 real(4), allocatable, dimension(:,:,:,:) :: dflxdo3
end type jacobianmat4d

type jacobianmatwrite

  integer :: ncid, varid(500)

end type jacobianmatwrite


! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine allocate_fields(config,fields)

 implicit none

 type(configuration), intent(in)    :: config
 type(datafields),    intent(inout) :: fields

 allocate(  fields%pl(config%is:config%ie,config%js:config%je,config%lm))
 allocate(   fields%t(config%is:config%ie,config%js:config%je,config%lm))
 allocate(   fields%q(config%is:config%ie,config%js:config%je,config%lm))
 allocate(  fields%qi(config%is:config%ie,config%js:config%je,config%lm))
 allocate(  fields%ql(config%is:config%ie,config%js:config%je,config%lm))
 allocate(  fields%ri(config%is:config%ie,config%js:config%je,config%lm))
 allocate(  fields%rl(config%is:config%ie,config%js:config%je,config%lm))
 allocate(  fields%o3(config%is:config%ie,config%js:config%je,config%lm))
 allocate(fields%fcld(config%is:config%ie,config%js:config%je,config%lm))

 allocate(  fields%ts(config%is:config%ie,config%js:config%je))
 allocate(fields%emis(config%is:config%ie,config%js:config%je))
 allocate(fields%lats(config%is:config%ie,config%js:config%je))
 allocate(fields%lons(config%is:config%ie,config%js:config%je))

 allocate(fields%flx  (config%is:config%ie,config%js:config%je,0:config%lm))
 allocate(fields%dfdts(config%is:config%ie,config%js:config%je,0:config%lm))

end subroutine allocate_fields

! ------------------------------------------------------------------------------

subroutine deallocate_fields(fields)

 implicit none

 type(datafields), intent(inout) :: fields

 deallocate(  fields%pl)
 deallocate(   fields%t)
 deallocate(   fields%q)
 deallocate(  fields%qi)
 deallocate(  fields%ql)
 deallocate(  fields%ri)
 deallocate(  fields%rl)
 deallocate(  fields%o3)
 deallocate(fields%fcld)

 deallocate(  fields%ts)
 deallocate(fields%emis)
 deallocate(fields%lats)
 deallocate(fields%lons)

 deallocate(fields%flx)
 deallocate(fields%dfdts)

end subroutine deallocate_fields

! ------------------------------------------------------------------------------

subroutine copy_fields(config,fields_lhs,fields_rhs,is_in,ie_in,js_in,je_in)

 implicit none

 type(configuration), intent(in) :: config
 type(datafields), intent(inout) :: fields_lhs
 type(datafields), intent(in)    :: fields_rhs
 integer, optional, intent(in)   :: is_in,ie_in,js_in,je_in

 integer :: is,ie,js,je

 is = config%is
 ie = config%ie
 js = config%js
 je = config%je

 if (present(is_in)) is = is_in
 if (present(ie_in)) ie = ie_in
 if (present(js_in)) js = js_in
 if (present(je_in)) je = je_in

 fields_lhs%pl   (is:ie,js:je,:)  = fields_rhs%pl   (is:ie,js:je,:)
 fields_lhs%t    (is:ie,js:je,:)  = fields_rhs%t    (is:ie,js:je,:)
 fields_lhs%q    (is:ie,js:je,:)  = fields_rhs%q    (is:ie,js:je,:)
 fields_lhs%qi   (is:ie,js:je,:)  = fields_rhs%qi   (is:ie,js:je,:)
 fields_lhs%ql   (is:ie,js:je,:)  = fields_rhs%ql   (is:ie,js:je,:)
 fields_lhs%ri   (is:ie,js:je,:)  = fields_rhs%ri   (is:ie,js:je,:)
 fields_lhs%rl   (is:ie,js:je,:)  = fields_rhs%rl   (is:ie,js:je,:)
 fields_lhs%o3   (is:ie,js:je,:)  = fields_rhs%o3   (is:ie,js:je,:)
 fields_lhs%fcld (is:ie,js:je,:)  = fields_rhs%fcld (is:ie,js:je,:)
 fields_lhs%ts   (is:ie,js:je)    = fields_rhs%ts   (is:ie,js:je)
 fields_lhs%emis (is:ie,js:je)    = fields_rhs%emis (is:ie,js:je)
 fields_lhs%lats (is:ie,js:je)    = fields_rhs%lats (is:ie,js:je)
 fields_lhs%lons (is:ie,js:je)    = fields_rhs%lons (is:ie,js:je)
 fields_lhs%flx  (is:ie,js:je,:)  = fields_rhs%flx  (is:ie,js:je,:)
 fields_lhs%dfdts(is:ie,js:je,:)  = fields_rhs%dfdts(is:ie,js:je,:)

end subroutine copy_fields

! ------------------------------------------------------------------------------

subroutine allocate_fluxes(config,fluxes)

 implicit none

 type(configuration), intent(in)    :: config
 type(fluxfields),    intent(inout) :: fluxes

 allocate(fluxes%flx  (config%is:config%ie,config%js:config%je,0:config%lm))
 allocate(fluxes%dfdts(config%is:config%ie,config%js:config%je,0:config%lm))

end subroutine allocate_fluxes

! ------------------------------------------------------------------------------

subroutine deallocate_fluxes(fluxes)

 implicit none

 type(fluxfields), intent(inout) :: fluxes

 deallocate(fluxes%flx)
 deallocate(fluxes%dfdts)

end subroutine deallocate_fluxes

! ------------------------------------------------------------------------------

subroutine allocate_jacobian(jacobian,config)

 implicit none

 type(configuration), intent(in)    :: config
 type(jacobianmat),   intent(inout) :: jacobian

 allocate(jacobian%dflxdpl  (config%lm,config%lm))
 allocate(jacobian%dflxdt   (config%lm,config%lm))
 allocate(jacobian%dflxdq   (config%lm,config%lm))
 allocate(jacobian%dflxdqi  (config%lm,config%lm))
 allocate(jacobian%dflxdql  (config%lm,config%lm))
 allocate(jacobian%dflxdo3  (config%lm,config%lm))

end subroutine allocate_jacobian

! ------------------------------------------------------------------------------

subroutine allocate_jacobian4d(jacobian,config)

 implicit none

 type(configuration), intent(in)    :: config
 type(jacobianmat4d), intent(inout) :: jacobian

 allocate(jacobian%dflxdpl(config%is:config%ie,config%js:config%je,config%lm,config%lm))
 allocate(jacobian%dflxdt (config%is:config%ie,config%js:config%je,config%lm,config%lm))
 allocate(jacobian%dflxdq (config%is:config%ie,config%js:config%je,config%lm,config%lm))
 allocate(jacobian%dflxdqi(config%is:config%ie,config%js:config%je,config%lm,config%lm))
 allocate(jacobian%dflxdql(config%is:config%ie,config%js:config%je,config%lm,config%lm))
 allocate(jacobian%dflxdo3(config%is:config%ie,config%js:config%je,config%lm,config%lm))

end subroutine allocate_jacobian4d

! ------------------------------------------------------------------------------

subroutine deallocate_jacobian(jacobian)

 implicit none

 type(jacobianmat), intent(inout) :: jacobian

 deallocate(jacobian%dflxdpl)
 deallocate(jacobian%dflxdt)
 deallocate(jacobian%dflxdq)
 deallocate(jacobian%dflxdqi)
 deallocate(jacobian%dflxdql)
 deallocate(jacobian%dflxdo3)

end subroutine deallocate_jacobian

! ------------------------------------------------------------------------------

subroutine deallocate_jacobian4d(jacobian)

 implicit none

 type(jacobianmat4d), intent(inout) :: jacobian

 deallocate(jacobian%dflxdpl)
 deallocate(jacobian%dflxdt)
 deallocate(jacobian%dflxdq)
 deallocate(jacobian%dflxdqi)
 deallocate(jacobian%dflxdql)
 deallocate(jacobian%dflxdo3)

end subroutine deallocate_jacobian4d

! ------------------------------------------------------------------------------

subroutine write_jacobian_setup(self,xsize,ysize,config,filename,comm,info)

 implicit none

 type(jacobianmatwrite), intent(inout) :: self
 integer,                intent(in)    :: xsize
 integer,                intent(in)    :: ysize
 type(configuration),    intent(in)    :: config
 character(len=2048),    intent(in)    :: filename
 integer, optional,      intent(in)    :: comm
 integer, optional,      intent(in)    :: info

 integer :: x_dimid, y_dimid, z_dimid, vc


  ! Create file
  ! -----------
  if (present(comm) .and. present(info)) then
    call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_MPIIO), self%ncid, &
                               comm = comm, info = info), "nf90_create" )
  else
    call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_CLOBBER), self%ncid), "nf90_create" )
  endif

  ! Define dimensions
  ! -----------------
  call nccheck ( nf90_def_dim(self%ncid, "Xdim", xsize,     x_dimid), "nf90_def_dim Xdim" )
  call nccheck ( nf90_def_dim(self%ncid, "Ydim", ysize,     y_dimid), "nf90_def_dim Ydim" )
  call nccheck ( nf90_def_dim(self%ncid, "lev",  config%lm, z_dimid), "nf90_def_dim lev" )


  ! Define variables
  ! ----------------
  vc = 0

  !Grid
  vc = vc + 1
  call nccheck( nf90_def_var(self%ncid, "lons", NF90_FLOAT, (/ x_dimid, y_dimid /), self%varid(vc)), "nf90_def_var lons" )
  vc = vc + 1
  call nccheck( nf90_def_var(self%ncid, "lats", NF90_FLOAT, (/ x_dimid, y_dimid /), self%varid(vc)), "nf90_def_var lats" )
  vc = vc + 1
  call nccheck( nf90_def_var(self%ncid, "levs", NF90_FLOAT, (/ z_dimid /),          self%varid(vc)), "nf90_def_var levs" )

  ! Variables
  vc = vc + 1
  call nccheck( nf90_def_var(self%ncid, "dflxdpl", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, z_dimid /), &
                self%varid(vc)), "nf90_def_var dflxdpl" )
  vc = vc + 1
  call nccheck( nf90_def_var(self%ncid, "dflxdt",  NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, z_dimid /), &
                self%varid(vc)), "nf90_def_var dflxdt" )
  vc = vc + 1
  call nccheck( nf90_def_var(self%ncid, "dflxdq",  NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, z_dimid /), &
                 self%varid(vc)), "nf90_def_var dflxdq" )
  vc = vc + 1
  call nccheck( nf90_def_var(self%ncid, "dflxdqi", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, z_dimid /), &
                 self%varid(vc)), "nf90_def_var dflxdqi" )
  vc = vc + 1
  call nccheck( nf90_def_var(self%ncid, "dflxdql", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, z_dimid /), &
                 self%varid(vc)), "nf90_def_var dflxdql" )
  vc = vc + 1
  call nccheck( nf90_def_var(self%ncid, "dflxdo3", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, z_dimid /), &
                 self%varid(vc)), "nf90_def_var dflxdo3" )

  ! End define
  call nccheck( nf90_enddef(self%ncid), "nf90_enddef" )

end subroutine write_jacobian_setup

! ------------------------------------------------------------------------------

subroutine write_jacobian_variables(self,config,fields,jacobian)

 implicit none

 type(jacobianmatwrite), intent(inout) :: self
 type(configuration),    intent(in)    :: config
 type(datafields),       intent(in)    :: fields
 type(jacobianmat4d),    intent(in)    :: jacobian

 integer :: istart(4), icount(4)
 integer :: vc
 integer :: k
 integer :: levs(config%lm)
 integer :: is,ie,js,je

 ! Write Jacobian arrays
 ! ---------------------

 is = config%is
 ie = config%ie
 js = config%js
 je = config%je

 ! Starts and counts for writing
 istart(1) = 1 !config%is
 istart(2) = 1 !config%js
 istart(3) = 1
 istart(4) = 1
 icount(1) = config%ie - config%is + 1
 icount(2) = config%je - config%js + 1
 icount(3) = config%lm
 icount(4) = config%lm

 do k = 1, config%lm
   levs(k) = k
 enddo

 vc = 0

 !Grid
 vc = vc + 1
 call nccheck( nf90_put_var( self%ncid, self%varid(vc), fields%lons(is:ie,js:je), &
                             start = istart(1:2), count = icount(1:2) ), "nf90_put_var lons" )
 vc = vc + 1
 call nccheck( nf90_put_var( self%ncid, self%varid(vc), fields%lats(is:ie,js:je), &
                             start = istart(1:2), count = icount(1:2) ), "nf90_put_var lats" )
 vc = vc + 1
 call nccheck( nf90_put_var( self%ncid, self%varid(vc), levs ), "nf90_put_var levs" )

 !Jacobian
 vc = vc + 1
 call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian%dflxdpl(is:ie,js:je,:,:), &
                             start = istart, count = icount ), "nf90_put_var dflxdpl" )
 vc = vc + 1
 call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian%dflxdt (is:ie,js:je,:,:), &
                             start = istart, count = icount ), "nf90_put_var dflxdt"  )
 vc = vc + 1
 call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian%dflxdq (is:ie,js:je,:,:), &
                             start = istart, count = icount ), "nf90_put_var dflxdq"  )
 vc = vc + 1
 call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian%dflxdqi(is:ie,js:je,:,:), &
                             start = istart, count = icount ), "nf90_put_var dflxdqi" )
 vc = vc + 1
 call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian%dflxdql(is:ie,js:je,:,:), &
                             start = istart, count = icount ), "nf90_put_var dflxdql" )
 vc = vc + 1
 call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian%dflxdo3(is:ie,js:je,:,:), &
                             start = istart, count = icount ), "nf90_put_var dflxdo3" )

end subroutine write_jacobian_variables

! ------------------------------------------------------------------------------

subroutine write_jacobian_final(self)

implicit none
type(jacobianmatwrite), intent(in) :: self

! Close file
! ----------
call nccheck ( nf90_close(self%ncid), "nf90_close" )

end subroutine write_jacobian_final

! ------------------------------------------------------------------------------

subroutine write_jacobian_writegathered(self,config,commsize,jacobian,iinRecv,jinRecv,latRecv,lonRecv)

 implicit none

 type(jacobianmatwrite), intent(inout) :: self
 type(configuration),    intent(in)    :: config
 integer,                intent(in)    :: commsize
 type(jacobianmat),      intent(in)    :: jacobian(commsize)
 integer,                intent(in)    :: iinRecv(commsize)
 integer,                intent(in)    :: jinRecv(commsize)
 real(4),                intent(in)    :: latRecv(commsize)
 real(4),                intent(in)    :: lonRecv(commsize)

 integer :: istart(4), icount(4)
 integer :: vc, n
 integer :: k
 integer :: levs(config%lm)
 logical, save :: donelevs = .false.

 ! Write Jacobian arrays
 ! ---------------------

 do n = 1,commsize

   ! Starts and counts for writing
   istart(1) = iinRecv(n)
   istart(2) = jinRecv(n)
   istart(3) = 1
   istart(4) = 1
   icount(1) = 1
   icount(2) = 1
   icount(3) = config%lm
   icount(4) = config%lm

!   istart(1) = 1
!   istart(2) = 1
!   istart(3) = jinRecv(n)
!   istart(4) = iinRecv(n)
!   icount(1) = config%lm
!   icount(2) = config%lm
!   icount(3) = 1
!   icount(4) = 1

   vc = 0

   !Grid
   vc = vc + 1
   call nccheck( nf90_put_var( self%ncid, self%varid(vc), lonRecv(n:n), &
                               start = istart(1:2), count = icount(1:2) ), "nf90_put_var lons" )

   vc = vc + 1
   call nccheck( nf90_put_var( self%ncid, self%varid(vc), latRecv(n:n), &
                               start = istart(1:2), count = icount(1:2) ), "nf90_put_var lats" )

   vc = vc + 1
   if (.not. donelevs) then
     do k = 1, config%lm
       levs(k) = k
     enddo
     call nccheck( nf90_put_var( self%ncid, self%varid(vc), levs ), "nf90_put_var levs" )
     donelevs = .true.
   endif

   vc = vc + 1
   call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian(n)%dflxdpl(:,:), &
                               start = istart, count = icount ), "nf90_put_var dflxdpl" )
   vc = vc + 1
   call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian(n)%dflxdt (:,:), &
                               start = istart, count = icount ), "nf90_put_var dflxdt"  )
   vc = vc + 1
   call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian(n)%dflxdq (:,:), &
                               start = istart, count = icount ), "nf90_put_var dflxdq"  )
   vc = vc + 1
   call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian(n)%dflxdqi(:,:), &
                               start = istart, count = icount ), "nf90_put_var dflxdqi" )
   vc = vc + 1
   call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian(n)%dflxdql(:,:), &
                               start = istart, count = icount ), "nf90_put_var dflxdql" )
   vc = vc + 1
   call nccheck( nf90_put_var( self%ncid, self%varid(vc), jacobian(n)%dflxdo3(:,:), &
                               start = istart, count = icount ), "nf90_put_var dflxdo3" )

 enddo



end subroutine write_jacobian_writegathered

! ------------------------------------------------------------------------------

subroutine read_fields(config,fields)

 implicit none

 type(configuration), intent(in)    :: config
 type(datafields),    intent(inout) :: fields

 integer :: ncid, varid
 integer :: istart3(4), icount3(4)
 integer :: istart2(3), icount2(3)
 character(len=20) :: field

 real, allocatable, dimension(:,:,:) :: flxu_int
 real, allocatable, dimension(:,:,:) :: flxd_int

 !Read ranges
 istart3(1) = config%is
 istart3(2) = config%js
 istart3(3) = 1
 istart3(4) = 1
 icount3(1) = config%ie - config%is + 1
 icount3(2) = config%je - config%js + 1
 icount3(3) = config%lm
 icount3(4) = 1

 istart2(1) = istart3(1)
 istart2(2) = istart3(2)
 istart2(3) = istart3(4)
 icount2(1) = icount3(1)
 icount2(2) = icount3(2)
 icount2(3) = icount3(4)

 ! Input file
 ! ----------

 call nccheck ( nf90_open(trim(config%filename_in), NF90_NOWRITE, ncid), "nf90_open"//trim(config%filename_in) )

 field = "pl"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%pl(config%is:config%ie,config%js:config%je,:), istart3, icount3 ),&
                "nf90_get_var "//trim(field) )

 field = "t"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%t (config%is:config%ie,config%js:config%je,:), istart3, icount3 ),&
                "nf90_get_var "//trim(field) )

 field = "q"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%q (config%is:config%ie,config%js:config%je,:), istart3, icount3 ),&
                "nf90_get_var "//trim(field) )

 field = "qi"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%qi(config%is:config%ie,config%js:config%je,:), istart3, icount3 ),&
               "nf90_get_var "//trim(field) )

 field = "ql"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%ql(config%is:config%ie,config%js:config%je,:), istart3, icount3 ),&
               "nf90_get_var "//trim(field) )

 field = "ri"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%ri(config%is:config%ie,config%js:config%je,:), istart3, icount3 ),&
               "nf90_get_var "//trim(field) )

 field = "rl"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%rl(config%is:config%ie,config%js:config%je,:), istart3, icount3 ),&
               "nf90_get_var "//trim(field) )

 field = "o3"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%o3(config%is:config%ie,config%js:config%je,:), istart3, icount3 ),&
               "nf90_get_var "//trim(field) )

 field = "fcld"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%fcld(config%is:config%ie,config%js:config%je,:), istart3, icount3 ),&
             "nf90_get_var "//trim(field) )

 field = "ts"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%ts  (config%is:config%ie,config%js:config%je), istart2, icount2 ),&
               "nf90_get_var "//trim(field) )

 field = "emis"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%emis(config%is:config%ie,config%js:config%je), istart2, icount2 ),&
             "nf90_get_var "//trim(field) )

 field = "lats"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%lats(config%is:config%ie,config%js:config%je), istart2, icount2 ),&
             "nf90_get_var "//trim(field) )

 field = "lons"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%lons(config%is:config%ie,config%js:config%je), istart2, icount2 ),&
             "nf90_get_var "//trim(field) )

 call nccheck ( nf90_close(ncid), "nf90_close" )


 ! Output file
 ! -----------

 allocate(flxu_int(config%is:config%ie,config%js:config%je,0:config%lm))
 allocate(flxd_int(config%is:config%ie,config%js:config%je,0:config%lm))

 call nccheck ( nf90_open(trim(config%filename_out), NF90_NOWRITE, ncid), "nf90_open"//trim(config%filename_out) )

 icount3(3) = config%lm+1

 field = "flxu"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, flxu_int(config%is:config%ie,config%js:config%je,:), istart3, icount3 ), &
                "nf90_get_var "//trim(field) )

 field = "flxd"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, flxd_int(config%is:config%ie,config%js:config%je,:), istart3, icount3 ), &
                "nf90_get_var "//trim(field) )

 field = "dfdts"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%dfdts(config%is:config%ie,config%js:config%je,:), istart3, icount3 ), &
                "nf90_get_var "//trim(field) )

 call nccheck ( nf90_close(ncid), "nf90_close" )

 ! Sum the fluxes
 fields%flx(config%is:config%ie,config%js:config%je,:) = flxd_int(config%is:config%ie,config%js:config%je,:) +&
                                                         flxu_int(config%is:config%ie,config%js:config%je,:)

 deallocate(flxu_int)
 deallocate(flxd_int)

end subroutine read_fields

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

subroutine compare_fluxes(config,fields,fluxes)

 implicit none
 type(configuration), intent(in) :: config
 type(datafields),    intent(in) :: fields
 type(fluxfields),    intent(in) :: fluxes

 integer :: i,j,k
 real, allocatable, dimension(:,:) :: flx_err
 real, allocatable, dimension(:,:) :: dfdts_err

 !Write the result with what GEOS produced
 open (unit = 101, file = "output.txt")
 do j=config%js,config%je
   do i=config%is,config%ie
     write(101,*) 'Profile print versus GEOS'
     write(101,*) 'Latitude', fields%lats(i,j), ',  longitude', fields%lons(i,j)
     write(101,*) ' '
     write(101,*) ' flx'
     do k=0,config%lm
       write(101,*) fields%flx(i,j,k), fluxes%flx(i,j,k), &
                    100*abs(fields%flx(i,j,k)-fluxes%flx(i,j,k))/(fields%flx(i,j,k)+1e-6)
     enddo
     write(101,*) ' dfdts '
     do k=0,config%lm
       write(101,*) fields%dfdts(i,j,k), fluxes%dfdts(i,j,k), &
                    100*abs(fields%dfdts(i,j,k)-fluxes%dfdts(i,j,k))/(fields%dfdts(i,j,k)+1e-6)
     enddo
   enddo
 enddo
 close(101)

 allocate(flx_err(config%is:config%ie,config%js:config%je))
 allocate(dfdts_err(config%is:config%ie,config%js:config%je))

 flx_err = sqrt(sum((fields%flx-fluxes%flx)**2,3)/config%lm)
 dfdts_err = sqrt(sum((fields%dfdts-fluxes%dfdts)**2,3)/config%lm)

 do j=config%js,config%je
   do i=config%is,config%ie
     print*, i,j
     print*, 'flx rmse: ', flx_err(i,j)
     print*, 'dfdts rmse: ', dfdts_err(i,j)
   enddo
 enddo

 deallocate(flx_err,dfdts_err)

end subroutine compare_fluxes

! ------------------------------------------------------------------------------

end module rtmg_lw_tools_mod
