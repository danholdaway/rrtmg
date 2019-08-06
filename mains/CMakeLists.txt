ecbuild_add_executable( TARGET  lw_driver.x
                        SOURCES lw_driver.f90
                        LIBS	${NETCDF_LIBRARIES} rrtmg
                      )

#ecbuild_add_executable( TARGET  lw_jacobian.x
#                        SOURCES lw_jacobian.f90
#                        LIBS	${NETCDF_LIBRARIES} rrtmg
#                      )

ecbuild_add_executable( TARGET  lw_jacobian_6procs.x
                        SOURCES lw_jacobian_6procs.f90
                        LIBS	${NETCDF_LIBRARIES} rrtmg
                      )

ecbuild_add_executable( TARGET  lw_jacobian_allwrite.x
                        SOURCES lw_jacobian_allwrite.f90
                        LIBS	${NETCDF_LIBRARIES} rrtmg
                      )
