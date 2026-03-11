
!----------------------
!s_wf(spin[down down, up down, down up, up up], ispin[down up, up down], particle number[1,2],coordinate[1,2,3])
!s_wf(spin[down down, up down, down up, up up], ispin[down up, up down],
! sigma_1 coordinate[1,2,3], sigma_2 coordinate[1,2,3] ) 
!
!
!
module spin_calc
  implicit none
  complex*16::s_wf(4,2,2,3),ss_wf(4,2,3,3)
contains
  subroutine run_spin(n,niso,wf_in)
    implicit none
    integer intent(in) n,niso
    
