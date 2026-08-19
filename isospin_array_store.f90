module iso_arr_store
  implicit none
  integer,allocatable::p0arr(:,:,:,:),p1arr(:,:)
  integer,allocatable::p2arr(:,:,:),pxarr(:,:,:)
  complex*16,allocatable::c0arr(:,:,:,:),c1arr(:,:)
  complex*16,allocatable::c2arr(:,:,:),cxarr(:,:,:)

contains

  subroutine isospin_ops()
    use mpi_modules
    use param_calc
    use sim_ops
  implicit none

  complex*16::intcarr(2)
  integer::nstates,i,b,ncomb,p1,p2,intp
  ncomb=msg%npart*(msg%npart-1)/2
  intcarr(:)=0
  intp=0
!    write(*,*)"starthere"
    allocate(p0arr(msg%npart,msg%npart,3,nspin))
!  write(*,*)"p0"
  allocate(p1arr(msg%npart,nspin))
  allocate(p2arr(msg%npart,msg%npart,nspin))
  allocate(pxarr(msg%npart,msg%npart,nspin))
!  write(*,*)"swaparrays"
  allocate(c0arr(msg%npart,msg%npart,3,nspin))
  allocate(c1arr(msg%npart,nspin))
  allocate(c2arr(msg%npart,msg%npart,nspin))
  allocate(cxarr(msg%npart,msg%npart,nspin))
!  write(*,*)"dome with allocation"
  do p1=1,msg%npart-1
     do p2=p1+1,msg%npart
        do i = 1,nspin
           call pauli(p2,2,i,intp,intcarr(1))
           call pauli(p1,1,intp,pxarr(p1,p2,i),intcarr(2))
           cxarr(p1,p2,i)=intcarr(1)*intcarr(2)
           call pauli(p2,1,i,intp,intcarr(1))
           call pauli(p1,2,intp,pxarr(p2,p1,i),intcarr(2))
!           write(*,*)"dot and cross"
           cxarr(p2,p1,i)=intcarr(1)*intcarr(2)
           if (p1.eq.1)then
              if (p2.eq.2)then
                 call pauli(1,3,i,p1arr(1,i),c1arr(1,i))
 !                write(*,*)"p2:",1
 !                write(*,*)"vec"
              endif
 !             write(*,*)"p2:", p2
              call pauli(p2,3,i,p1arr(p2,i),c1arr(p2,i))
              
              
  !            write(*,*)"iniso:",p1arr(p2,i)
           endif
           call pauli(p2,3,i,intp,intcarr(1))
           call pauli(p1,3,intp,p2arr(p1,p2,i),intcarr(2))
           c2arr(p1,p2,i)=intcarr(1)*intcarr(2)
           do b = 1,3
              call pauli(p2,b,i,intp,intcarr(1))
              call pauli(p1,b,intp,p0arr(p1,p2,b,i),intcarr(2))
              c0arr(p1,p2,b,i)=intcarr(1)*intcarr(2)
           enddo
        enddo
     enddo
  enddo
 ! write(*,*)"finishing i array calc"
  return
end subroutine isospin_ops
end module iso_arr_store
