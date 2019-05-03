list( APPEND rrtmg_src_files
rrtmg_lw/gcm_model/modules/dbmemory.F90
rrtmg_lw/gcm_model/modules/_gpudef.inc
rrtmg_lw/gcm_model/modules/_memdef.inc
rrtmg_lw/gcm_model/modules/parrrtm.F90
rrtmg_lw/gcm_model/modules/rrlw_cld.F90
rrtmg_lw/gcm_model/modules/rrlw_con.F90
rrtmg_lw/gcm_model/modules/rrlw_kg01.F90
rrtmg_lw/gcm_model/modules/rrlw_kg02.F90
rrtmg_lw/gcm_model/modules/rrlw_kg03.F90
rrtmg_lw/gcm_model/modules/rrlw_kg04.F90
rrtmg_lw/gcm_model/modules/rrlw_kg05.F90
rrtmg_lw/gcm_model/modules/rrlw_kg06.F90
rrtmg_lw/gcm_model/modules/rrlw_kg07.F90
rrtmg_lw/gcm_model/modules/rrlw_kg08.F90
rrtmg_lw/gcm_model/modules/rrlw_kg09.F90
rrtmg_lw/gcm_model/modules/rrlw_kg10.F90
rrtmg_lw/gcm_model/modules/rrlw_kg11.F90
rrtmg_lw/gcm_model/modules/rrlw_kg12.F90
rrtmg_lw/gcm_model/modules/rrlw_kg13.F90
rrtmg_lw/gcm_model/modules/rrlw_kg14.F90
rrtmg_lw/gcm_model/modules/rrlw_kg15.F90
rrtmg_lw/gcm_model/modules/rrlw_kg16.F90
rrtmg_lw/gcm_model/modules/rrlw_ncpar.F90
rrtmg_lw/gcm_model/modules/rrlw_ref.F90
rrtmg_lw/gcm_model/modules/rrlw_tbl.F90
rrtmg_lw/gcm_model/modules/rrlw_vsn.F90
rrtmg_lw/gcm_model/modules/rrlw_wvn.F90
rrtmg_lw/gcm_model/src/gmcica_subcol_gen_lw.F90
rrtmg_lw/gcm_model/src/gpu_rrtmg_lw_cldprmc.F90
rrtmg_lw/gcm_model/src/gpu_rrtmg_lw_rtrnmc.F90
rrtmg_lw/gcm_model/src/gpu_rrtmg_lw_setcoef.F90
rrtmg_lw/gcm_model/src/gpu_rrtmg_lw_taumol.F90
rrtmg_lw/gcm_model/src/mcica_random_numbers.F90
rrtmg_lw/gcm_model/src/rrtmg_lw_init.F90
rrtmg_lw/gcm_model/src/rrtmg_lw_k_g.F90
rrtmg_lw/gcm_model/src/rrtmg_lw_rad.F90
rrtmg_lw/gcm_model/src/rrtmg_lw_rtrn.F90
rrtmg_lw/gcm_model/src/rrtmg_lw_rtrnmr.F90
rrtmg_lw/gcm_model/src/rrtmg_lw_setcoef.F90
rrtmg_lw/gcm_model/src/WaterDistributionMod.F90
)

#Append includes and libraries with particular model 
list (APPEND RRTMG_LIB_DEP ${NETCDF_LIBRARIES})

#Append src with particular model 
list (APPEND RRTMG_SRC_DEP ${rrtmg_src_files} )

#Add the library
ecbuild_add_library( TARGET          rrtmg
                     SOURCES         ${rrtmg_src_files}
                     LIBS            ${NETCDF_LIBRARIES}
                     LINKER_LANGUAGE ${RRTMG_LINKER_LANGUAGE}
                   )

