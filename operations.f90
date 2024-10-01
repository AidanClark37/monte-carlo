integer function state_sign(p,i)
  implicit none
  integer,intent(in)::i,p
  
  
  !sign=1-2*MOD(CEILING(MOD(DBLE(i/(2.d0**(p-1)),2.d0)),2) !given row index i and particle number p,
                                                          !return the sign of z component spin of particle p
  state_sign = 2*mod(floor((i-1)/(2.d0**(p-1))),2)-1
  return
end function state_sign

integer function flip(p,i)    ! give n row index i and particle number p, return the index of the row
  implicit none
  ! that results from flipping the spin of particle p
  integer::flip,i,p,state_sign                

  flip=i-state_sign(p,i)*2**(p-1)
  return
end function flip

! intege function iso_index(i,niso,iso_array)
!  integer::i,niso,iso_array(niso)
!  iso_index = iso_array(i)
!  return
!end function iso_index

subroutine  spin(A,p,b,N,sigma)   !A - matrix input, wavefunction
                                    !N - number of particles, matrix has 2^N rows 
  !p - particle number, 1,2, up to the number of particles N
  implicit none
  integer::b,p,i,N,state_sign,flip!i - pauli index, 1,2,3
  complex*16::A(4,2),sigma(4,2)
  select case(b)
     case(1)  !sigma_x
       do i = 1,2**N
          sigma(i,:)=A(flip(p,i),:)
       enddo
    case(2)   !sigma_y
       do i = 1,2**N
          sigma(i,:)=dcmplx(0,1)*state_sign(p,i)*A(flip(p,i),:)
       enddo
    case(3)  !sigma_z
       do i = 1,2**N
          sigma(i,:)=state_sign(p,i)*A(i,:)
       enddo
    end select
    return
  end subroutine spin
  
       
