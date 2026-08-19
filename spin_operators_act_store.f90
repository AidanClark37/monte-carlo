!-----------------------------------------------
!s_wf - nspin, niso, particle number, coordinate index
!calculate simga_1,2 x y and z operator actions
!ss_wf - nspin, niso, coordinate index for sigma_1 and sigma_2
!calculates sigma_1^i sigma_2^j for all i,j<=3
!t_wf - nspin, niso, particle number
!only calculates z component of tau_1,2
!tt_wf nspin, niso, coordinate index for tau_1 and tau_2
!only calculates tau_1^z * tau_2^z, tau_1^x * tau_2^y, tau_1^y * tau_2^x
!st_wf nspin, niso, particle number for sigma, coordinate index for sigma, particle number for tau
!sigma_n^i tau_m^z for n,m=1,2 i=1,2,3
!stt_wf niso, nspin particle number for sigma,coordinate index for sigma, coordinate index for tau_1 nad tau_2
!sigma_n^i for n=1,2 i=1,2,3, only calculates tau_1^z * tau_2^z, tau_1^x * tau_2^y, tau_1^y * tau_2^x
!sstt_wf nspin, niso, coordinate index for sigma_1 and sigma_2, coordinate index for tau_1 and tau_2
!calculates sigma_1^i sigma_2^j for all i,j<=3,only calculates tau_1^z * tau_2^z, tau_1^x * tau_2^y, tau_1^y * tau_2^x
!----------------------------------------------------
module spin_ops_act_store
  use mpi_modules
  use param_calc
  use spin_ops_act
  implicit none
  complex*16,allocatable::s_wf(:,:,:,:),ss_wf(:,:,:,:,:,:)!,t_wf(4,2,2),tt_wf(4,2,3,3),st_wf(4,2,2,3,2),stt_wf(4,2,2,3,3,3)
!  complex*16::sst_wf(4,2,3,3,2),sstt_wf(4,2,3,3,3,3)
contains
  subroutine all_operators(wf_in)
    use param_calc
    use mpi_modules
    use spin_ops_act
    integer::p1,p2,a1,a2,a3
    complex*16,intent(in)::wf_in(:,:)
    complex*16,allocatable::wf1(:,:),wf2(:,:),wf3(:,:),wf4(:,:)
    if (.not. allocated(s_wf)) allocate(s_wf(nspin,niso,msg%npart,3))

    if (.not. allocated(ss_wf)) allocate(ss_wf(nspin,niso,msg%npart,msg%npart,3,3))
    allocate(wf1(nspin,niso))

    allocate(wf2(nspin,niso))
allocate(wf3(nspin,niso))
allocate(wf4(nspin,niso))

    do p1 = 1, msg%npart
       do a1= 1,3

          call spin(wf_in,p1,a1,wf2)
          s_wf(:,:,p1,a1)=wf2(:,:)
          do p2 = 1,msg%npart
             do a2=1,3
                call spin(wf_in,p1,a1,wf1)
                call spin(wf1,p2,a2,wf2)
                ss_wf(:,:,p1,p2,a1,a2)=wf2(:,:)
             enddo
          enddo
       enddo
       
    enddo
    
             
return
  end subroutine all_operators
end module spin_ops_act_store






! spin(wf,p,b,N,niso,sigma_wf)
! isospin(iarray,wf,p,b,N,niso,tau_wf)
