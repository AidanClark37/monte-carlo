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

integer function flip(p,i)    ! give n row index i and particle number p, return the index of the row
  implicit none
  ! that results from flipping the spin of particle p
  integer::i,p,state_sign                
  !write(*,*)'start flip'

  flip=i-state_sign(p,i)*2**(p-1)
  !write(*,*)'end flip'
  return
end function flip

integer function iso_index(i,niso,iso_array)
  integer::i,niso,iso_array(niso),state_sign
  iso_index = iso_array(i)
  return
end function iso_index

subroutine  spin(wf,p,b,N,niso,sigma_wf)   !wf - matrix input, wavefunction
                                    !N - number of particles, matrix has 2^N rows 
  !p - particle number, 1,2, up to the number of particles N
  implicit none
  integer,intent(in)::p,b,N,niso
  integer::i,j,flip,nspin,state_sign
  complex*16,intent(in)::wf(4,2)
  complex*16,intent(out)::sigma_wf(4,2)
  nspin=2**N
!  write(*,*)'niso:',niso
  select case(b)
     case(1)  !sigma_x
        do i = 1,nspin
           do j = 1,niso              
              sigma_wf(i,j)=wf(flip(p,i),j)
!              write(*,*)i,j,' | ',flip(p,i),j
!              write(*,*)wf(i,j),'| ',wf(flip(p,i),j)
              !write(*,*)sigma_wf(i,j)
           enddo
        enddo
     case(2)   !sigma_y
!        write(*,*)'case 2'
        do i = 1,nspin
           do j=1,niso
              
             sigma_wf(i,j)=-dcmplx(0,1)*state_sign(p,i)*wf(flip(p,i),j)
          enddo
       enddo
       !write(*,*)'end case 2'
    case(3)  !sigma_z
       do i = 1,nspin
          do j=1,niso
 
             sigma_wf(i,j)=state_sign(p,i)*wf(i,j)
          enddo
       enddo
    case default
       write(*,*)'fail'   
    end select
    return
  end subroutine spin
  
  subroutine spin_exp_val(cwf,b,N,niso,nspin,spin_exp)
    implicit none
    integer,intent(in)::b,N,niso,nspin
    complex*16,intent(in)::cwf(nspin,niso)
    complex*16::spin_cwf(nspin,niso),cc_cwf(nspin,niso)
    real*8,intent(out)::spin_exp
    integer::p,i,j
    spin_exp=0.d0
    cc_cwf=conjg(cwf)
    do p = 1,N
       call spin(cwf,p,b,N,niso,spin_cwf)
       do i = 1,nspin
          do j=1,niso
             spin_exp=spin_exp+cc_cwf(i,j)*spin_cwf(i,j)
          enddo
       enddo
    enddo
  end subroutine spin_exp_val
