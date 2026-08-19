module diff
  implicit none
  contains
  subroutine derivative(wf_in,rr,i,wf_out)
  use mpi_modules
  use param_calc
  use deuteron
  implicit none
  integer,intent(in)::i
  integer::n
  real*8,intent(in)::rr(3,2)
  real*8::dr1(3),r1,rr1(3,2),h
  complex*16,intent(in)::wf_in(4,2)
  complex*16,intent(inout)::wf_out(4,2)
  complex*16,allocatable::wf1(:,:)
  allocate(wf1(nspin,niso))
  write(*,*)"spin",nspin,niso
  h=0.000d0
  rr1(:,:)=rr(:,:)
    rr1(i,2)=rr(i,2)+h
 !   write(*,*)"rr",rr(i,:)
 !   write(*,*)"rr1",rr1(i,:)
    
  call deut_wave(rr1,wf1,msg%wf,dr1,r1)


  wf_out(:,:) = (wf1(:,:)-wf_in(:,:))/h
 
 do n = 1,4
    write(*,*)"in",wf_in(n,1),wf_in(n,2)
    enddo
 do n = 1,4
    write(*,*)"shift",wf1(n,1),wf1(n,2)
 enddo
 stop
    return
end subroutine derivative
end module diff

