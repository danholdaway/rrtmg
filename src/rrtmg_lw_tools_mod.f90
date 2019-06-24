module rtmg_lw_tools_mod

use netcdf

implicit none

private

public :: configuration
public :: datafields, allocate_fields, deallocate_fields, copy_fields, read_fields
public :: fluxfields, allocate_fluxes, deallocate_fluxes, compare_fluxes
public :: jacobianmat, allocate_jacobian, deallocate_jacobian, write_jacobian

! Configuration settings
type configuration
 character(len=2048) :: filename_in
 character(len=2048) :: filename_out
 integer :: im,jm,lm
 integer :: is,ie,js,je
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
 character(len=2048) :: filename_out
 real, allocatable, dimension(:,:) :: dflxdpl
 real, allocatable, dimension(:,:) :: dflxdt
 real, allocatable, dimension(:,:) :: dflxdq
 real, allocatable, dimension(:,:) :: dflxdqi
 real, allocatable, dimension(:,:) :: dflxdql
 real, allocatable, dimension(:,:) :: dflxdri
 real, allocatable, dimension(:,:) :: dflxdrl
 real, allocatable, dimension(:,:) :: dflxdo3
 real, allocatable, dimension(:,:) :: dflxdfcld
 real, allocatable, dimension(:)   :: dflxdts
 real, allocatable, dimension(:)   :: dflxdemis
 real, allocatable, dimension(:,:) :: dfdtsdpl
 real, allocatable, dimension(:,:) :: dfdtsdt
 real, allocatable, dimension(:,:) :: dfdtsdq
 real, allocatable, dimension(:,:) :: dfdtsdqi
 real, allocatable, dimension(:,:) :: dfdtsdql
 real, allocatable, dimension(:,:) :: dfdtsdri
 real, allocatable, dimension(:,:) :: dfdtsdrl
 real, allocatable, dimension(:,:) :: dfdtsdo3
 real, allocatable, dimension(:,:) :: dfdtsdfcld
 real, allocatable, dimension(:)   :: dfdtsdts
 real, allocatable, dimension(:)   :: dfdtsdemis
end type jacobianmat


! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine allocate_fields(config,fields)

 implicit none

 type(configuration), intent(in)    :: config
 type(datafields),    intent(inout) :: fields

 allocate(  fields%pl(config%im,config%jm,config%lm))
 allocate(   fields%t(config%im,config%jm,config%lm))
 allocate(   fields%q(config%im,config%jm,config%lm))
 allocate(  fields%qi(config%im,config%jm,config%lm))
 allocate(  fields%ql(config%im,config%jm,config%lm))
 allocate(  fields%ri(config%im,config%jm,config%lm))
 allocate(  fields%rl(config%im,config%jm,config%lm))
 allocate(  fields%o3(config%im,config%jm,config%lm))
 allocate(fields%fcld(config%im,config%jm,config%lm))

 allocate(  fields%ts(config%im,config%jm))
 allocate(fields%emis(config%im,config%jm))
 allocate(fields%lats(config%im,config%jm))
 allocate(fields%lons(config%im,config%jm))

 allocate(fields%flx  (config%im,config%jm,0:config%lm))
 allocate(fields%dfdts(config%im,config%jm,0:config%lm))

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

subroutine copy_fields(fields_lhs,fields_rhs)

 implicit none

 type(datafields), intent(inout) :: fields_lhs
 type(datafields), intent(in)    :: fields_rhs

 fields_lhs%pl    = fields_rhs%pl
 fields_lhs%t     = fields_rhs%t
 fields_lhs%q     = fields_rhs%q
 fields_lhs%qi    = fields_rhs%qi
 fields_lhs%ql    = fields_rhs%ql
 fields_lhs%ri    = fields_rhs%ri
 fields_lhs%rl    = fields_rhs%rl
 fields_lhs%o3    = fields_rhs%o3
 fields_lhs%fcld  = fields_rhs%fcld
 fields_lhs%ts    = fields_rhs%ts
 fields_lhs%emis  = fields_rhs%emis
 fields_lhs%lats  = fields_rhs%lats
 fields_lhs%lons  = fields_rhs%lons
 fields_lhs%flx   = fields_rhs%flx
 fields_lhs%dfdts = fields_rhs%dfdts

end subroutine copy_fields

! ------------------------------------------------------------------------------

subroutine allocate_fluxes(config,fluxes)

 implicit none

 type(configuration), intent(in)    :: config
 type(fluxfields),    intent(inout) :: fluxes

 allocate(fluxes%flx  (config%im,config%jm,0:config%lm))
 allocate(fluxes%dfdts(config%im,config%jm,0:config%lm))

end subroutine allocate_fluxes

! ------------------------------------------------------------------------------

subroutine deallocate_fluxes(fluxes)

 implicit none

 type(fluxfields), intent(inout) :: fluxes

 deallocate(fluxes%flx)
 deallocate(fluxes%dfdts)

end subroutine deallocate_fluxes

! ------------------------------------------------------------------------------

subroutine allocate_jacobian(config,jacobian)

 implicit none

 type(configuration), intent(in)    :: config
 type(jacobianmat),   intent(inout) :: jacobian

 allocate(jacobian%dflxdpl  (config%lm,config%lm))
 allocate(jacobian%dflxdt   (config%lm,config%lm))
 allocate(jacobian%dflxdq   (config%lm,config%lm))
 allocate(jacobian%dflxdqi  (config%lm,config%lm))
 allocate(jacobian%dflxdql  (config%lm,config%lm))
 allocate(jacobian%dflxdri  (config%lm,config%lm))
 allocate(jacobian%dflxdrl  (config%lm,config%lm))
 allocate(jacobian%dflxdo3  (config%lm,config%lm))
 allocate(jacobian%dflxdfcld(config%lm,config%lm))
 allocate(jacobian%dflxdts  (config%lm))
 allocate(jacobian%dflxdemis(config%lm))

 allocate(jacobian%dfdtsdpl  (config%lm,config%lm))
 allocate(jacobian%dfdtsdt   (config%lm,config%lm))
 allocate(jacobian%dfdtsdq   (config%lm,config%lm))
 allocate(jacobian%dfdtsdqi  (config%lm,config%lm))
 allocate(jacobian%dfdtsdql  (config%lm,config%lm))
 allocate(jacobian%dfdtsdri  (config%lm,config%lm))
 allocate(jacobian%dfdtsdrl  (config%lm,config%lm))
 allocate(jacobian%dfdtsdo3  (config%lm,config%lm))
 allocate(jacobian%dfdtsdfcld(config%lm,config%lm))
 allocate(jacobian%dfdtsdts  (config%lm))
 allocate(jacobian%dfdtsdemis(config%lm))

end subroutine allocate_jacobian

! ------------------------------------------------------------------------------

subroutine deallocate_jacobian(jacobian)

 implicit none

 type(jacobianmat), intent(inout) :: jacobian

 deallocate(jacobian%dflxdpl)
 deallocate(jacobian%dflxdt)
 deallocate(jacobian%dflxdq)
 deallocate(jacobian%dflxdqi)
 deallocate(jacobian%dflxdql)
 deallocate(jacobian%dflxdri)
 deallocate(jacobian%dflxdrl)
 deallocate(jacobian%dflxdo3)
 deallocate(jacobian%dflxdfcld)
 deallocate(jacobian%dflxdts)
 deallocate(jacobian%dflxdemis)

 deallocate(jacobian%dfdtsdpl)
 deallocate(jacobian%dfdtsdt)
 deallocate(jacobian%dfdtsdq)
 deallocate(jacobian%dfdtsdqi)
 deallocate(jacobian%dfdtsdql)
 deallocate(jacobian%dfdtsdri)
 deallocate(jacobian%dfdtsdrl)
 deallocate(jacobian%dfdtsdo3)
 deallocate(jacobian%dfdtsdfcld)
 deallocate(jacobian%dfdtsdts)
 deallocate(jacobian%dfdtsdemis)

end subroutine deallocate_jacobian

! ------------------------------------------------------------------------------

subroutine write_jacobian(config,fields,jacobian)

 implicit none

 type(configuration), intent(in) :: config
 type(datafields),    intent(in) :: fields
 type(jacobianmat),   intent(in) :: jacobian

 integer :: ncid, z_dimid, n_dimid, vc, varid(500)

 call nccheck( nf90_create( trim(jacobian%filename_out), NF90_NETCDF4, ncid), "nf90_create" )

 call nccheck ( nf90_def_dim(ncid, "lev", config%lm, z_dimid), "nf90_def_dim lev" )
 call nccheck ( nf90_def_dim(ncid, "n",           1, n_dimid), "nf90_def_dim n"   )


 ! Define mode
 vc = 0

 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "lat", NF90_DOUBLE, n_dimid, varid(vc)), "nf90_def_var lat" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "lon", NF90_DOUBLE, n_dimid, varid(vc)), "nf90_def_var lon" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdpl", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dflxdpl" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdt", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dflxdt" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdq", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dflxdq" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdqi", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dflxdqi" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdql", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dflxdql" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdri", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dflxdri" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdrl", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dflxdrl" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdo3", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dflxdo3" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdfcld", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dflxdfcld" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdts", NF90_DOUBLE, (/ z_dimid /), varid(vc)), "nf90_def_var dflxdts" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dflxdemis", NF90_DOUBLE, (/ z_dimid /), varid(vc)), "nf90_def_var dflxdemis" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdpl", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dfdtsdpl" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdt", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dfdtsdt" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdq", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dfdtsdq" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdqi", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dfdtsdqi" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdql", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dfdtsdql" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdri", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dfdtsdri" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdrl", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dfdtsdrl" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdo3", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dfdtsdo3" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdfcld", NF90_DOUBLE, (/ z_dimid, z_dimid /), varid(vc)), "nf90_def_var dfdtsdfcld" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdts", NF90_DOUBLE, (/ z_dimid /), varid(vc)), "nf90_def_var dfdtsdts" )
 vc = vc + 1
 call nccheck( nf90_def_var(ncid, "dfdtsdemis", NF90_DOUBLE, (/ z_dimid /), varid(vc)), "nf90_def_var dfdtsdemis" )

 call nccheck( nf90_enddef(ncid), "nf90_enddef" )

 ! Write mode
 vc = 0

 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), fields%lats(config%is,config%js) ), "nf90_put_var lat" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), fields%lons(config%is,config%js) ), "nf90_put_var lon" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdpl ), "nf90_put_var dflxdpl" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdt ), "nf90_put_var dflxdt" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdq ), "nf90_put_var dflxdq" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdqi ), "nf90_put_var dflxdqi" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdql ), "nf90_put_var dflxdql" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdri ), "nf90_put_var dflxdri" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdrl ), "nf90_put_var dflxdrl" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdo3 ), "nf90_put_var dflxdo3" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdfcld ), "nf90_put_var dflxdfcld" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdts ), "nf90_put_var dflxdts" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dflxdemis ), "nf90_put_var dflxdemis" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdpl ), "nf90_put_var dfdtsdpl" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdt ), "nf90_put_var dfdtsdt" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdq ), "nf90_put_var dfdtsdq" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdqi ), "nf90_put_var dfdtsdqi" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdql ), "nf90_put_var dfdtsdql" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdri ), "nf90_put_var dfdtsdri" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdrl ), "nf90_put_var dfdtsdrl" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdo3 ), "nf90_put_var dfdtsdo3" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdfcld ), "nf90_put_var dfdtsdfcld" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdts ), "nf90_put_var dfdtsdts" )
 vc = vc + 1
 call nccheck( nf90_put_var( ncid, varid(vc), jacobian%dfdtsdemis ), "nf90_put_var dfdtsdemis" )

 ! Close file
 call nccheck ( nf90_close(ncid), "nf90_close" )

end subroutine write_jacobian

! ------------------------------------------------------------------------------

subroutine read_fields(config,fields)

 implicit none

 type(configuration), intent(in)    :: config
 type(datafields),    intent(inout) :: fields

 integer :: ncid, varid
 integer, allocatable, dimension(:)    :: istart3, icount3
 integer, allocatable, dimension(:)    :: istart2, icount2
 character(len=20) :: field

 real, allocatable, dimension(:,:,:) :: flxu_int
 real, allocatable, dimension(:,:,:) :: flxd_int

 !Read ranges
 allocate(istart3(4), icount3(4))
 allocate(istart2(3), icount2(3))

 istart3(1) = config%is; icount3(1) = config%ie - config%is + 1
 istart3(2) = config%js; icount3(2) = config%je - config%js + 1
 istart3(3) = 1 ; icount3(3) = config%lm
 istart3(4) = 1 ; icount3(4) = 1

 istart2(1) = config%is; icount2(1) = config%ie - config%is + 1
 istart2(2) = config%js; icount2(2) = config%je - config%js + 1
 istart2(3) = 1 ; icount2(3) = 1


 ! Input file
 ! ----------

 call nccheck ( nf90_open(trim(config%filename_in), NF90_NOWRITE, ncid), "nf90_open"//trim(config%filename_in) )

 field = "pl"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%pl, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "t"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%t, istart3, icount3 ),    "nf90_get_var "//trim(field) )

 field = "q"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%q, istart3, icount3 ),    "nf90_get_var "//trim(field) )

 field = "qi"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%qi, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "ql"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%ql, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "ri"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%ri, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "rl"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%rl, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "o3"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%o3, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "fcld"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%fcld, istart3, icount3 ), "nf90_get_var "//trim(field) )

 field = "ts"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%ts, istart2, icount2 ),   "nf90_get_var "//trim(field) )

 field = "emis"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%emis, istart2, icount2 ), "nf90_get_var "//trim(field) )

 field = "lats"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%lats, istart2, icount2 ), "nf90_get_var "//trim(field) )

 field = "lons"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),          "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%lons, istart2, icount2 ), "nf90_get_var "//trim(field) )

 call nccheck ( nf90_close(ncid), "nf90_close" )


 ! Output file
 ! -----------

 allocate(flxu_int(config%im,config%jm,0:config%lm))
 allocate(flxd_int(config%im,config%jm,0:config%lm))

 call nccheck ( nf90_open(trim(config%filename_out), NF90_NOWRITE, ncid), "nf90_open"//trim(config%filename_out) )

 icount3(3) = config%lm+1

 field = "flxu"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, flxu_int, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "flxd"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, flxd_int, istart3, icount3 ),   "nf90_get_var "//trim(field) )

 field = "dfdts"
 call nccheck ( nf90_inq_varid( ncid, trim(field), varid ),                 "nf90_inq_var "//trim(field) )
 call nccheck ( nf90_get_var( ncid, varid, fields%dfdts, istart3, icount3 ),  "nf90_get_var "//trim(field) )

 call nccheck ( nf90_close(ncid), "nf90_close" )

 ! Sum the fluxes
 fields%flx = flxd_int + flxu_int

 deallocate(flxu_int)
 deallocate(flxd_int)

 deallocate(istart3, istart2)
 deallocate(icount3, icount2)

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
       write(101,*) fields%flx(i,j,k), fluxes%flx(i,j,k), 100*abs(fields%flx(i,j,k)-fluxes%flx(i,j,k))/(fields%flx(i,j,k)+1e-6)
     enddo
     write(101,*) ' dfdts '
     do k=0,config%lm
       write(101,*) fields%dfdts(i,j,k), fluxes%dfdts(i,j,k), 100*abs(fields%dfdts(i,j,k)-fluxes%dfdts(i,j,k))/(fields%dfdts(i,j,k)+1e-6)
     enddo
   enddo
 enddo
 close(101)

 allocate(flx_err(config%im,config%jm))
 allocate(dfdts_err(config%im,config%jm))

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
