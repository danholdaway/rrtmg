####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

if( HAVE_OMP )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
else( )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-openmp")
endif( )

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check -fdefault-real-8")

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -mtune=generic -funroll-loops -g -D__GFORTRAN__ -ffree-line-length-none -Wno-missing-include-dirs -fPIC -ffpe-trap=zero,overflow -fbacktrace -falign-commons")

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -mtune=generic -funroll-loops -g -D__GFORTRAN__ -ffree-line-length-none -Wno-missing-include-dirs -fPIC -ffpe-trap=zero,overflow -fbacktrace -falign-commons" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 -funroll-all-loops -finline-functions" )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

####################################################################

# Meaning of flags
# ----------------
# -fstack-arrays     : Allocate automatic arrays on the stack (needs large stacksize!!!)
# -funroll-all-loops : Unroll all loops
# -fcheck=bounds     : Bounds checking

