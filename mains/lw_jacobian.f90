program lw_driver

 use rtmg_lw_tools_mod
 use rtmg_lw_driver_mod, only: rrtmg_lw_driver

 implicit none

 type(datafields) :: fields
 type(fluxfields) :: fluxes

 type(jacobianmat) :: jacobian

 ! Some user input for the data to be used
 ! ---------------------------------------
 fields%is = 1                   !Starting point i direction
 fields%js = 1                   !Starting point j direction
 fields%ie = 1                  !End point i direction
 fields%je = 1                  !End point j direction

 fields%filename_in  = '/gpfsm/dnb31/drholdaw/Victor/IRRADTrainingData/f522_dh.trainingdata_in.lcv.20190401_0000z.nc4'  !Training data in
 fields%filename_out = '/gpfsm/dnb31/drholdaw/Victor/IRRADTrainingData/f522_dh.trainingdata_out.lcv.20190401_0000z.nc4' !Training data out

 fields%doy = 90                 !Which day of the year is it?


 ! Allocate memory
 call allocate_fields(fields)
 call allocate_fluxes(fluxes,fields)
 call allocate_jacobian(jacobian,fields)

 ! Read fields from model output data
 call read_fields(fields)

 ! Call radiation scheme driver
 call rrtmg_lw_driver(fields,fluxes)

 ! Write the output comparison
 call write_jacobian(jacobian)

 ! Deallocate memory
 call deallocate_fields(fields)
 call deallocate_fluxes(fluxes)
 call deallocate_jacobian(jacobian)

end program
