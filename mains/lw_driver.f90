program lw_driver

 use rtmg_lw_tools_mod
 use rtmg_lw_driver_mod, only: rrtmg_lw_driver

 implicit none

 type(configuration) :: config
 type(datafields) :: fields
 type(fluxfields) :: fluxes


 ! Some user input for the data to be used
 ! ---------------------------------------
 config%is = 1                  !Starting point i direction
 config%js = 1                  !Starting point j direction
 config%ie = 1                  !End point i direction
 config%je = 1                  !End point j direction

 config%filename_in  = '/gpfsm/dnb31/drholdaw/Victor/IRRADTrainingData/f522_dh.trainingdata_in.lcv.20190401_0000z.nc4'  !Training data in
 config%filename_out = '/gpfsm/dnb31/drholdaw/Victor/IRRADTrainingData/f522_dh.trainingdata_out.lcv.20190401_0000z.nc4' !Training data out

 config%doy = 90                 !Which day of the year is it?

 config%im = (config%ie-config%is+1)
 config%jm = (config%je-config%js+1)
 config%lm = 72


 ! Run component
 ! -------------

 ! Allocate memory
 call allocate_fields(config,fields)
 call allocate_fluxes(config,fluxes)

 ! Read fields from model output data
 call read_fields(config,fields)

 ! Call radiation scheme driver
 call rrtmg_lw_driver(config,fields,fluxes)

 ! Write the output comparison
 call compare_fluxes(config,fields,fluxes)

 ! Deallocate memory
 call deallocate_fields(fields)
 call deallocate_fluxes(fluxes)

end program
