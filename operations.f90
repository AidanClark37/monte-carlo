
  integer function state_sign(p,i)
  implicit none
  integer,intent(in)::i,p
  
  !sign=1-2*MOD(CEILING(MOD(DBLE(i/(2.d0**(p-1)),2.d0)),2) !given row index i and particle number p,
  !write(*,*) 'start sign'                                                        !return the sign of z component spin of particle p
  !write(*,*)'here'
  state_sign = 2*mod((i-1)/(2**(p-1)),2)-1

  !ss=i
  !write(*,*) '2^(p-1)',2**(p-1)
  !write(*,*) '(i-1)/2^(p-1)',(i-1)/(2**(p-1))
  !write(*,*) 'mod((i-1)/2^(p-1),2)',mod((i-1)/(2**(p-1)),2)
 ! write(*,*) '
  return
end function state_sign

integer function flip(p,i)    ! given row index i and particle number p, return the index of the row
  implicit none
  integer,external::state_sign                            ! that results from flipping the spin of particle p
  integer::i,p                
  !write(*,*)'start flip'

  flip=i-state_sign(p,i)*2**(p-1)
  !write(*,*)'end flip'
  return
end function flip

integer function iso_index(i,niso,iso_array)
  integer::i,niso,iso_array(niso)
  integer,external::state_sign
  iso_index = iso_array(i)
  return
end function iso_index

!subroutine param_calc(npart,Tz)
!  implicit none
!  integer,intent(in)::npart,Tz
!  integer,intent(out)::niso,nspin
!  integer,allocatable::iarray(:),inviarray(:)
!  integer::i,p,j,state_sign,nu,itot,count 
!  real*8::gamma
!  nspin=2**npart
!  allocate(itot(nspin))
!  allocate(inviarray(npsin))
!  nu = (Tz+npart)/2
  
!  niso=gamma(real(npart+1))/(int(dgamma(dble(nu+1)))*int(dgamma(dble(npart-nu+1))))
!  allocate(iarray(niso))
!  count=1
!  do i=1,nspin

!     do p = 1, npart
!        itot=0
!        itot = itot + state_sign(p,i)
!     enddo
!     if(itot.eq.Tz)then
!        iso_states(count)=itot
!        count=count+1
!     endif
!  enddo
!  do i =1,niso
!     inviarray(iarray(i)) = i
     
 ! write(*,*)"nspin!",nspin,gamma(real(nspin+1))
 ! write(*,*)"nu!",nu,(gamma(real(nu+1)))
 ! write(*,*)"(nspin-nu)!",(gamma(real(nspin-nu+1)))
 ! write(*,*)"niso",niso
  




 

 
 

  !write(*,*)"test"
!  do i = 1,nspin
!    j=1
!    if(itot(i)==Tz) then
!       write(*,*)"one"
!       iarray(j)=i
!       write(*,*)"two"
!       j=j+1
!    endif
! enddo
! write(*,*)iarray(1)


!end subroutine param_calc

  
  

subroutine  spin(wfin,p,b,sigma_wf)   !wf - matrix input, wavefunction
                   
  !p - particle number, 1,2, up to the number of particles N
  !b - coordinate index
  !N - total number of particles
  !niso
  use mpi_modules
  use param_calc
  implicit none
  integer,intent(in)::p,b
  integer::i,j
  integer,external::flip,state_sign
  complex*16,intent(in)::wfin(4,2)
  complex*16,intent(out)::sigma_wf(4,2)
  !write(*,*)'niso:',niso
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
              
             sigma_wf(i,j)=-dcmplx(0,1)*state_sign(p,i)*wfin(flip(p,i),j)
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
  end subroutine spin

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !input for pauli:
  !p: particle number, which paricle the pauli matrix is operating on
  !values: 1 through the number of particles (n)
  !b: coordinate index, which pauli matrix is operating
  !values: 1,2,3
  !i: row index, which row of the state vector is being moved
  !values: 1 through 2**n
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !ouput for pauli:
  !iout: where the row indicated by index i would be swapped with
  !values: 1 through 2**n
  !cout: the constant this row would be multiplied by for a specific pauli matrix
  !values: -1, 1, -i, or i

  subroutine pauli(p,b,i,iout,cout)
    implicit none
    integer,intent(in)::p,b,i
    integer,intent(out)::iout
    complex*16,intent(out)::cout
    integer,external::state_sign,flip
    
    select case(b)
    case(1)
       iout = flip(p,i)
       cout=1
    case(2)
       iout= flip(p,i)
       cout= -cmplx(0,1)*state_sign(p,i)
    case(3)
       iout=i
       cout=state_sign(p,i)
    case default
       write(*,*) "fail"
    end select
    
  end subroutine pauli

  
  subroutine isospin(wf_in,p,b,tau_wf)
    use mpi_modules
    use param_calc
    implicit none
    integer,intent(in)::p,b
    integer::i,j,jj,k
    integer,external::flip,state_sign
  complex*16,intent(in)::wf_in(4,2)
  complex*16,intent(out)::tau_wf(4,2)
  !allocate(iarray(niso))
  select case(b)
     case(1)  !tau_x                                                                                                  
        do i = 1,nspin
           do j = 1,niso
              !write(*,*)'for j=',j, iarray(j),flip(p,iarray(j))
              jj = flip(p,iarray(j))
              if (ANY(iarray==jj)) then
                 tau_wf(i,j)=wf_in(i,jj)
              else
                 tau_wf(i,j)=cmplx(0,0)

              endif
              !write(*,*)i,j,' | ',flip(p,i),j                                                                          
!              write(*,*)wf(i,j),'| ',wf(flip(p,i),j)                                                                   
              !write(*,*)sigma_wf(i,j)                                                                                  
           enddo
        enddo
     case(2)   !tau_y                                                                                                 
!        write(*,*)'case 2'                                                                                             
        do i = 1,nspin
           do j=1,niso
              jj=flip(p,iarray(j))
              if (ANY(iarray==jj)) then
                 tau_wf(i,j)=-dcmplx(0,1)*state_sign(p,iarray(j))*wf_in(i,jj)
              else
                 tau_wf(i,j)=cmplx(0,0)

              endif

             !tau_wf(i,j)=-dcmplx(0,1)*state_sign(p,iarray(j))*wf(i,flip(p,iarray(j)))
          enddo
       enddo
       !write(*,*)'end case 2'                                                                                          
    case(3)  !tau_z                                                                                                  
       do i = 1,nspin
          do j=1,niso
             tau_wf(i,j)=state_sign(p,iarray(j))*wf_in(i,j)
             !write(*,*)i,j,tau_wf(i,j),state_sign(p,iarray(j)),wf(i,j)
          enddo
       enddo
    case default
       write(*,*)'fail'
    end select
       !do i = 1,nspin
          !do j=1,niso
             !write(*,*)i,j,wf(i,j),state_sign(p,iarray(j)),tau_wf(i,j)
          !enddo
       !enddo
    return
  end subroutine isospin

  !subroutine rho_NNg(iarray,wf_in,N,niso,pm,exp)
  !  implicit none
  !  integer,intent(in)::iarray(niso),N,niso
  !  integer::p,nspin,i,j,pm,state_sign
  !  complex*16,intent(in)::wf_in(4,niso)
  !!  real*8,intent(out)::exp
  !  complex*16,allocatable::cc_wf(:,:),tau_wf(:,:),spin_tau_wf(:,:),wf_out(:,:),tau_one_wf(:,:)
  !  nspin=2**N
    
 !   allocate(cc_wf(nspin,niso))
 !   allocate(tau_wf(nspin,niso))
 !   allocate(spin_tau_wf(nspin,niso))
 !   allocate(wf_out(nspin,niso))
 !   allocate(tau_one_wf(nspin,niso))
 !   exp=0.d0
 !   wf_out=dcmplx(0.d0,0.d0)
 !   do p = 1,N
 !   call isospin(wf_in,p,3,tau_wf)
    !tau_one_wf(:,:)=tau_wf(:,:) + wf_in(:,:)
 !         tau_one_wf = 0.5d0*(tau_wf + wf_in)
    !   end do
    !end do
    !write(*,*)"particle = ",p
    !do j=1,2
    !   write(*,*)j,state_sign(p,iarray(j))
    !end do
!    if(proc_rank.eq.3)then
!    do i=1,4
!       write(*,*)"spin = ",i
!       write(*,*)"wf ",(wf_in(i,j),j=1,2)
!       write(*,*)"tau_z",(tau_wf(i,j),j=1,2)
!       write(*,*)"1+tau_z wf ",(tau_one_wf(i,j),j=1,2)
!    end do
!    end if
    
 !   call spin(tau_one_wf,p,3,2,2,spin_tau_wf)
 !         wf_out = wf_out + spin_tau_wf
    
 !enddo

! cc_wf = conjg(wf_in)
 
!    do i = 1,nspin
!       do j = 1,niso
          !if(proc_rank.eq.1)write(*,*)wf_in(i,j),cc_wf(i,j),wf_out(i,j)
!          exp =	exp + cc_wf(i,j)*wf_out(i,j)
          !write(*,*)proc_rank,i,j,exp
!       enddo
       
!    enddo
!    write(*,*) proc_rank,exp
    !stop
!return

!end subroutine rho_NNg

  subroutine tau_exp_val(iarray,cwf,b,N,niso,nspin,tau_exp)
    implicit none
    integer,intent(in)::iarray(niso),b,N,niso,nspin
    complex*16,intent(in)::cwf(nspin,niso)
    integer::p,i,j
    real*8::tau_exp
    complex*16,allocatable::cc_wf(:,:),tau_wf(:,:)
    allocate(cc_wf(nspin,niso))
    allocate(tau_wf(nspin,niso))
    do p = 1,N
       call isospin(cwf,p,3,tau_wf)
    enddo
    cc_wf=conjg(cwf)
    tau_exp=0.d0
    do i = 1,nspin
       do j = 1,niso
          tau_exp=tau_exp+cc_wf(i,j)*tau_wf(i,j)
          !write(*,*)tau_exp
       enddo
    enddo
    return
  end subroutine tau_exp_val
  





