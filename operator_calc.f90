!-----------------------------------------------
!s_wf - nspin, niso, particle number, coordinate index
!calculate simga_1,2 x y and z operator actions
!ss_wf - nspin, niso, coordinate index for sigma_2 and sigma_2
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
module operator_calc
  implicit none
  complex*16::s_wf(4,2,2,3),ss_wf(4,2,3,3),t_wf(4,2,2),tt_wf(4,2,3,3),st_wf(4,2,2,3,2),stt_wf(4,2,2,3,3,3)
  complex*16::sst_wf(4,2,3,3,2),sstt_wf(4,2,3,3,3,3)
contains
  subroutine all_operators(iarray,wf_in,N,niso)
    integer,intent(in)::n,niso,iarray(niso)
    integer::p1,p2,a1,a2,a3,nspin
    complex*16,intent(in)::wf_in(2**N,niso)
    complex*16,allocatable::wf1(:,:),wf2(:,:),wf3(:,:),wf4(:,:)
    nspin=2**N
    allocate(wf1(nspin,niso))

    allocate(wf2(nspin,niso))
allocate(wf3(nspin,niso))
allocate(wf4(nspin,niso))

    do p1 = 1, N
      
       call isospin(iarray,wf_in,p1,3,N,niso,wf1)
       
       t_wf(:,:,p1)=wf1(:,:)
       
       
      
       
       do a1= 1,3
          call spin(wf_in,p1,a1,N,niso,wf2)
          s_wf(:,:,p1,a1)=wf2(:,:)
          do p2 =1, N
             call isospin(iarray,wf2,p2,3,N,niso,wf3)
             st_wf(:,:,p1,a1,p2)=wf3(:,:)
          enddo
          
             do a2=1,3
       call isospin(iarray,wf2,1,a2,N,niso,wf3)
       call isospin(iarray,wf3,2,a2,N,niso,wf4)
       stt_wf(:,:,p1,a1,a2,a2)=wf4(:,:)
    enddo
    call isospin(iarray,wf2,1,1,N,niso,wf3)
    call isospin(iarray,wf3,2,2,N,niso,wf4)
    stt_wf(:,:,p1,a1,1,2)=wf4(:,:)
    call isospin(iarray,wf2,1,2,N,niso,wf3)
    call isospin(iarray,wf3,2,1,N,niso,wf4)
    stt_wf(:,:,p1,a1,2,1)=wf4(:,:)
   
    
 enddo
enddo
    do a1=1,3
       call isospin(iarray,wf_in,1,a1,N,niso,wf1)
       call isospin(iarray,wf1,2,a1,N,niso,wf2)
       tt_wf(:,:,a1,a1)=wf2(:,:)
    enddo
    call isospin(iarray,wf_in,1,1,N,niso,wf1)
    call isospin(iarray,wf1,2,2,N,niso,wf2)
    tt_wf(:,:,1,2)=wf2(:,:)
    call isospin(iarray,wf_in,1,2,N,niso,wf1)
    call isospin(iarray,wf1,2,1,N,niso,wf2)
    tt_wf(:,:,2,1)=wf2(:,:)
    do a1=1,3
       do a2=1,3
          call spin(wf_in,a1,1,N,niso,wf1)
          call spin(wf1,a2,2,N,niso,wf2)
          ss_wf(:,:,a1,a2)=wf2(:,:)
          do p1 = 1,N
             call isospin(iarray,wf2,p1,3,N,niso,wf3)
             sst_wf(:,:,a1,a2,p1)=wf3(:,:)
          enddo
          do a3=1,3
             call isospin(iarray,wf2,1,a3,N,niso,wf3)
             call isospin(iarray,wf3,2,a3,N,niso,wf4)
             sstt_wf(:,:,a1,a2,a3,a3)=wf4(:,:)
          enddo
          call isospin(iarray,wf2,1,1,N,niso,wf3)
          call isospin(iarray,wf3,2,2,N,niso,wf4)
          sstt_wf(:,:,a1,a2,1,2)=wf4(:,:)
          call isospin(iarray,wf2,1,2,N,niso,wf3)
          call isospin(iarray,wf3,2,1,N,niso,wf4)
          sstt_wf(:,:,a1,a2,2,1)=wf4(:,:)
       enddo

    enddo
    
             

  end subroutine all_operators
  end module operator_calc






! spin(wf,p,b,N,niso,sigma_wf)
! isospin(iarray,wf,p,b,N,niso,tau_wf)
