module sim_ops
  !contains flip and statesign functions, and pauli subroutine
  implicit none
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !input: wavefunction array index i corresonding to one spin combination                                           
  !       particle index p looking at the state of                                                                  
  !output: sign of spin of particle p in spin combination i
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer function state_sign(p,i)
    implicit none
    integer,intent(in)::i,p
    state_sign = 2*mod((i-1)/(2**(p-1)),2)-1
  return
end function state_sign


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!input: wavefunction array index i corresonding to one spin combination                                           
!       particle index p looking at the state of                                                                  
!output: wavefunction index corresponding to the state of i with the spin of particle p flipped
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer function flip(p,i)    
  implicit none
  integer::i,p                
  flip=i-state_sign(p,i)*2**(p-1)
  return
end function flip

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !input for pauli:
  !p: particle number, which paricle the pauli matrix is operating on
  !values: 1 through the number of particles (n)
  !b: coordinate index, which pauli matrix is operating
  !values: 1,2,3
  !i: row index, index of output array
  !values: 1 through 2**n
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !ouput for pauli:
  !iout: row index, index of input array (yes confusing, easier to change definition here than everywhere else)
  !values: 1 through 2**n
  !cout: the constant that i row of output array would be multiplied by for a specific pauli matrix
  !values: -1, 1, -i, or i
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine pauli(p,b,i,iout,cout)
    implicit none
    integer,intent(in)::p,b,i
    integer,intent(out)::iout
    complex*16,intent(out)::cout
    !integer,external::state_sign,flip
    
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
       write(*,*) "coordinate index of pauli matrix is out of bounds"
       stop
    end select
    return
    
  end subroutine pauli

  
  end module sim_ops
