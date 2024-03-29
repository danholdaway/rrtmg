####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

if( HAVE_OMP )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
else( )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp-stubs")
endif( )

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -traceback")

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -qopt-report0 -ftz -align all -fno-alias -traceback -assume realloc_lhs -convert big_endian -fPIC -fpe0 -fp-model source -heap-arrays 32 -assume noold_maxminloc -fimf-arch-consistency=true -align dcommons" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0 -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -warn unused" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 -ip -ipo -unroll -inline -no-heap-arrays" )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

####################################################################

# Meaning of flags
# ----------------
# todo

