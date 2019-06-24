module rtmg_lw_tools_mod

use netcdf

implicit none

private

public :: datafields, allocate_fields, deallocate_fields, copy_fields, read_fields
public :: fluxfields, allocate_fluxes, deallocate_fluxes, compare_fluxes
public :: jacobianmat, allocate_jacobian, deallocate_jacobian, write_jacobian


type datafields
 character(len=2048) :: filename_in
 character(len=2048) :: filename_out
 integer :: im,jm,lm
 integer :: is,ie,js,je
 integer :: doy
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

type fluxfields
 real, allocatable, dimension(:,:,:) :: flx
 real, allocatable, dimension(:,:,:) :: dfdts
end type fluxfields

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
end type jacobianmat


! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine allocate_fields(fields)

 implicit none

 type(datafields), intent(inout) :: fields

 fields%im = (fields%ie-fields%is+1)
 fields%jm = (fields%je-fields%js+1)
 fields%lm = 72

 allocate(  fields%pl(fields%im,fields%jm,fields%lm))
 allocate(   fields%t(fields%im,fields%jm,fields%lm))
 allocate(   fields%q(fields%im,fields%jm,fields%lm))
 allocate(  fields%qi(fields%im,fields%jm,fields%lm))
 allocate(  fields%ql(fields%im,fields%jm,fields%lm))
 allocate(  fields%ri(fields%im,fields%jm,fields%lm))
 allocate(  fields%rl(fields%im,fields%jm,fields%lm))
 allocate(  fields%o3(fields%im,fields%jm,fields%lm))
 allocate(fields%fcld(fields%im,fields%jm,fields%lm))

 allocate(  fields%ts(fields%im,fields%jm))
 allocate(fields%emis(fields%im,fields%jm))
 allocate(fields%lats(fields%im,fields%jm))
 allocate(fields%lons(fields%im,fields%jm))

 allocate(fields%flx  (fields%im,fields%jm,0:fields%lm))
 allocate(fields%dfdts(fields%im,fields%jm,0:fields%lm))

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

subroutine allocate_fluxes(fluxes,fields)

 implicit none

 type(fluxfields), intent(inout) :: fluxes
 type(datafields), intent(in)    :: fields

 allocate(fluxes%flx  (fields%im,fields%jm,0:fields%lm))
 allocate(fluxes%dfdts(fields%im,fields%jm,0:fields%lm))

end subroutine allocate_fluxes

! ------------------------------------------------------------------------------

subroutine deallocate_fluxes(fluxes)

 implicit none

 type(fluxfields), intent(inout) :: fluxes

 deallocate(fluxes%flx)
 deallocate(fluxes%dfdts)

end subroutine deallocate_fluxes

! ------------------------------------------------------------------------------

subroutine allocate_jacobian(jacobian,fields)

 implicit none

 type(jacobianmat), intent(inout) :: jacobian
 type(datafields), intent(in)    :: fields

 allocate(jacobian%dflxdpl  (fields%lm,fields%lm))
 allocate(jacobian%dflxdt   (fields%lm,fields%lm))
 allocate(jacobian%dflxdq   (fields%lm,fields%lm))
 allocate(jacobian%dflxdqi  (fields%lm,fields%lm))
 allocate(jacobian%dflxdql  (fields%lm,fields%lm))
 allocate(jacobian%dflxdri  (fields%lm,fields%lm))
 allocate(jacobian%dflxdrl  (fields%lm,fields%lm))
 allocate(jacobian%dflxdo3  (fields%lm,fields%lm))
 allocate(jacobian%dflxdfcld(fields%lm,fields%lm))
 allocate(jacobian%dflxdts  (fields%lm))
 allocate(jacobian%dflxdemis(fields%lm))

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

end subroutine deallocate_jacobian

! ------------------------------------------------------------------------------

subroutine write_jacobian(jacobian)

 implicit none

 type(jacobianmat), intent(inout) :: jacobian


end subroutine write_jacobian

! ------------------------------------------------------------------------------

subroutine read_fields(fields)

 implicit none

 type(datafields), intent(inout) :: fields

 integer :: ncid, varid
 integer, allocatable, dimension(:)    :: istart3, icount3
 integer, allocatable, dimension(:)    :: istart2, icount2
 character(len=20) :: field

 real, allocatable, dimension(:,:,:) :: flxu_int
 real, allocatable, dimension(:,:,:) :: flxd_int

 !Read ranges
 allocate(istart3(4), icount3(4))
 allocate(istart2(3), icount2(3))

 istart3(1) = fields%is; icount3(1) = fields%ie - fields%is + 1
 istart3(2) = fields%js; icount3(2) = fields%je - fields%js + 1
 istart3(3) = 1 ; icount3(3) = fields%lm
 istart3(4) = 1 ; icount3(4) = 1

 istart2(1) = fields%is; icount2(1) = fields%ie - fields%is + 1
 istart2(2) = fields%js; icount2(2) = fields%je - fields%js + 1
 istart2(3) = 1 ; icount2(3) = 1


 ! Input file
 ! ----------

 call nccheck ( nf90_open(trim(fields%filename_in), NF90_NOWRITE, ncid), "nf90_open"//trim(fields%filename_in) )

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

 allocate(flxu_int(fields%im,fields%jm,0:fields%lm))
 allocate(flxd_int(fields%im,fields%jm,0:fields%lm))

 call nccheck ( nf90_open(trim(fields%filename_out), NF90_NOWRITE, ncid), "nf90_open"//trim(fields%filename_out) )

 icount3(3) = fields%lm+1

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

subroutine compare_fluxes(fields,fluxes)

 implicit none
 type(datafields), intent(in)    :: fields
 type(fluxfields), intent(inout) :: fluxes

 integer :: i,j,k
 real, allocatable, dimension(:,:) :: flx_err
 real, allocatable, dimension(:,:) :: dfdts_err

 !Write the result with what GEOS produced
 open (unit = 101, file = "output.txt")
 do j=fields%js,fields%je
   do i=fields%is,fields%ie
     write(101,*) 'Profile print versus GEOS'
     write(101,*) 'Latitude', fields%lats(i,j), ',  longitude', fields%lons(i,j)
     write(101,*) ' '
     write(101,*) ' flx'
     do k=0,fields%lm
       write(101,*) fields%flx(i,j,k), fluxes%flx(i,j,k), 100*abs(fields%flx(i,j,k)-fluxes%flx(i,j,k))/(fields%flx(i,j,k)+1e-6)
     enddo
     write(101,*) ' dfdts '
     do k=0,fields%lm
       write(101,*) fields%dfdts(i,j,k), fluxes%dfdts(i,j,k), 100*abs(fields%dfdts(i,j,k)-fluxes%dfdts(i,j,k))/(fields%dfdts(i,j,k)+1e-6)
     enddo
   enddo
 enddo
 close(101)

 allocate(flx_err(fields%im,fields%jm))
 allocate(dfdts_err(fields%im,fields%jm))

 flx_err = sqrt(sum((fields%flx-fluxes%flx)**2,3)/fields%lm)
 dfdts_err = sqrt(sum((fields%dfdts-fluxes%dfdts)**2,3)/fields%lm)

 do j=fields%js,fields%je
   do i=fields%is,fields%ie
     print*, i,j
     print*, 'flx rmse: ', flx_err(i,j)
     print*, 'dfdts rmse: ', dfdts_err(i,j)
   enddo
 enddo

 deallocate(flx_err,dfdts_err)

end subroutine compare_fluxes

! ------------------------------------------------------------------------------

end module rtmg_lw_tools_mod
