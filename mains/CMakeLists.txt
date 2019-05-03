ecbuild_add_executable( TARGET  lw_driver.x
                        SOURCES lw_driver.f90
                        LIBS	${NETCDF_LIBRARIES} rrtmg
                      )

