module tau_operator
contains
subroutine tau0(wf_in, p1,p2,wf_out)
  use mpi_modules
  use isospin_operators
  use param_calc
  implicit none
  complex*16,intent(in)::wf_in(:,:)
  complex*16,intent(out)::wf_out(:,:)
  complex*16,allocatable::temp_wf(:,:)
  integer,intent(in)::p1,p2
  integer::i,j,b
  allocate(temp_wf(nspin,niso))
  !wf_out(:,:)=0
  do b = 1,3
     do i = 1, nspin
        do j = 1, niso
           wf_out(i,j)=wf_out(i,j)+c0arr(p1,p2,b,iarray(j))*wf_in(i,inviarray(p0arr(p1,p2,b,iarray(j))))
        enddo
     enddo
   !  wf_out(:,:)=wf_out(:,:)+temp_wf(:,:)
  enddo
return
end subroutine tau0

subroutine tau1(wf_in,p,wf_out)
  use mpi_modules
  use isospin_operators
  use param_calc
  implicit none
  complex*16,intent(in)::wf_in(:,:)
  complex*16,intent(out)::wf_out(:,:)
  complex*16,allocatable::temp_wf(:,:)
  integer,intent(in)::p
  integer::i,j
allocate(temp_wf(nspin,niso))
  do i = 1,nspin
     do j = 1,niso
!        write(*,*)"iarray,p1arr",iarray(j),p1arr(p,iarray(j))
!        write(*,*)"c1",c1arr(p,iarray(j))
        wf_out(i,j)=c1arr(p,iarray(j))*wf_in(i,inviarray(p1arr(p,iarray(j))))
!        write(*,*)"wf",i,j,wf_in(i,inviarray(p1arr(p,iarray(j)))),wf_out(i,j)
     enddo
  enddo
  return
end subroutine tau1

subroutine tau2(wf_in,p1,p2,wf_out)
  use mpi_modules
  use isospin_operators
  use param_calc
  implicit none
  complex*16,intent(in)::wf_in(:,:)
  complex*16,intent(out)::wf_out(:,:)
  integer,intent(in)::p1,p2
  integer::i,j
  do i = 1, nspin
     do j = 1, niso
        wf_out(i,j)=c2arr(p1,p2,iarray(j))*wf_in(i,inviarray(p2arr(p1,p2,iarray(j))))
     enddo
  enddo
  return
end subroutine tau2

subroutine taux(wf_in,p1,p2,wf_out)
  use mpi_modules
  use isospin_operators
  use param_calc
  implicit none
  complex*16,intent(in)::wf_in(:,:)
  complex*16,intent(out)::wf_out(:,:)
  complex*16,allocatable::temp_wf(:,:)
  integer,intent(in)::p1,p2
  integer::i,j
  allocate(temp_wf(nspin,niso))
  do i = 1,nspin
     do j = 1,niso
        temp_wf(i,j)=cxarr(p1,p2,iarray(j))*wf_in(i,inviarray(pxarr(p1,p2,iarray(j))))
        temp_wf(i,j)=temp_wf(i,j)-cxarr(p2,p1,iarray(j))*wf_in(i,inviarray(pxarr(p2,p1,iarray(j))))
        
     enddo
  enddo
  return
end subroutine taux
end module tau_operator
