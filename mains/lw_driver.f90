program lw_driver

  use netcdf
  use rrtmg_lw_rad, only: rrtmg_lw  !  RRTMG Code
  use rrtmg_lw_init, only: rrtmg_lw_ini

  implicit none

  call RRTMG_LW_INI(1.004e3)

end program
