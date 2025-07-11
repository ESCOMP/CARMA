! Copyright 2025 University Corporation for Atmospheric Research
!
! Simple test to check the CARMA version
program test_carma_version
  use carma_version, only: get_carma_version

  implicit none

  ! Print the CARMA version
  write(*,*) "CARMA version: ", get_carma_version()

end program test_carma_version