module spin_ops_act
  implicit none
  contains
subroutine  spin2(wfin,p,b,sigma_wf)   !wf - matrix input, wavefunction
                 
  !p - particle number, 1,2, up to the number of particles N
  !b - coordinate index
  !N - total number of particles
  !niso
  use mpi_modules
  use param_calc
  use sim_ops
  implicit none
  integer,intent(in)::p,b
  integer::i,j
  !integer,external::flip,state_sign
  complex*16,intent(in)::wfin(4,2)
  complex*16,intent(out)::sigma_wf(4,2)
  !write(*,*)'niso:',niso
 ! write(*,*)shape(sigma_wf)
 ! writE(*,*)nspin,niso
  select case(b)
     case(1)  !sigma_x
        do i = 1,nspin
           do j = 1,niso





              sigma_wf(i,j)=wfin(flip(p,i),j)
!              write(*,*)i,j,' | ',flip(p,i),j
!              write(*,*)wf(i,j),'| ',wf(flip(p,i),j)
              !write(*,*)sigma_wf(i,j)
           enddo
        enddo
     case(2)   !sigma_y
!        write(*,*)'case 2'
        do i = 1,nspin
           do j=1,niso
              
             sigma_wf(i,j)=-cmplx(0.d0,1.d0)*state_sign(p,i)*wfin(flip(p,i),j)
          enddo
       enddo
       !write(*,*)'end case 2'
    case(3)  !sigma_z
       do i = 1,nspin
          do j=1,niso
 
             sigma_wf(i,j)=state_sign(p,i)*wfin(i,j)
          enddo
       enddo
    case default
       write(*,*)'fail'   
    end select
    return
  end subroutine spin2

subroutine spin(wf_in,p,b,wf_out)
use mpi_modules
  use param_calc
  use sim_ops
  implicit none
  integer,intent(in)::p,b
  integer::i,j,iout
  !integer,external::flip,state_sign
  complex*16::cout
  complex*16,intent(in)::wf_in(:,:)
  complex*16,intent(out)::wf_out(:,:)
  iout=0
  cout=dcmplx(0.d0,0.d0)
  !allocate(wf_out(nspin,niso))
     !spin x
     do i = 1,nspin
        call pauli(p,b,i,iout,cout)
        wf_out(i,:)=cout*wf_in(iout,:)
        ! write(*,*)shape(sigma_wf)
     enddo
     
 ! writE(*,*)nspin,niso
  
end subroutine spin

end module spin_ops_act

