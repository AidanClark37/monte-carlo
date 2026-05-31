module param_calc
implicit none
  integer::nspin,niso
  integer,allocatable::iarray(:),inviarray(:)
  real*8,parameter::gA=1.267,nmass=938,fpi=92.4
contains
  
  subroutine parameter_calc(npart,Tz)
    implicit none
    integer,intent(in)::npart,Tz

    

    integer::i,p,j,state_sign,nu,count,itot 
    real*8::gamma
  nspin=2**npart

  allocate(inviarray(nspin))
  nu = (Tz+npart)/2
  
  niso=int(dgamma(dble(npart+1)))/(int(dgamma(dble(nu+1)))*int(dgamma(dble(npart-nu+1))))
!  write(*,*)"npart!:",dgamma(dble(npart+1))
!  write(*,*)"nu!:",dgamma(dble(nu+1))
!  write(*,*)"nisoinfile",niso
  allocate(iarray(niso))
  if (allocated(iarray))then
!     write(*,*)"allocated"
     endif

 ! write(*,*)"[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
  count=1
  do i=1,nspin
     itot=0
     do p = 1, npart
        
        itot = itot + state_sign(p,i)
  !      write(*,*)"--------------------"
  !      write(*,*)"i=",i,"p=",p
  !      write(*,*)"--------------------"
  !      write(*,*)"sign:",state_sign(p,i)
  !      write(*,*)"itot:",itot
     enddo
  !   write(*,*)"========================="
  !   write(*,*)"itot for state",i,itot
  !   write(*,*)"=========================="
     if(itot.eq.Tz)then
  !      write(*,*)"adding to iarray"
  !      write(*,*)"itot:",itot
  !      write(*,*)"count:",count
  !      write(*,*)"Tz:",Tz
        iarray(count)=i
        count=count+1
     endif
  enddo
  !write(*,*)"[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["
  do i =1,niso
     inviarray(iarray(i)) = i
!     write(*,*)"inviarry",iarray(i),inviarray(iarray(i))
  end do
!  write(*,*)"iarray",iarray(1),iarray(2)
  return
end subroutine parameter_calc
end module param_calc
