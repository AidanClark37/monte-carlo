subroutine step(rnd,rpart_o,jz,cwf,norm0,rr,dr,r,acc)
  use mpi_modules
  use param_calc
  use deuteron
  implicit none


  real*8             ,intent(inout)::rnd
  real*8             ,intent(inout)::norm0
  integer            ,intent(in)::jz
  real*8             ,intent(inout)::rpart_o(3,2)
  complex*16         ,intent(out)::cwf(4,2)
  real*8::sig
  complex*16::warray(4,2)
  complex*16,allocatable::cwf0(:,:)
  real*8::ysol(2,40),dr(3),r
  logical,intent(out)::acc
  
  real*8::rr(3,2)
  real*8,allocatable::random(:)
  real*8::norm
  integer::k,j,i
allocate(random(3*msg%npart+2))
  sig=msg%sigma
  ysol=msg%wf
  call rndnb(3*msg%npart+2 , rnd , random )!call the random number
!  write(*,*)random(1)
!  stop
  !selection of the new position configuration
  do k=1,msg%npart
     rr(1,k)=rpart_o(1,k)+sig*(random(3*k-2)-0.5d0)
     rr(2,k)=rpart_o(2,k)+sig*(random(3*k-1)-0.5d0)
     rr(3,k)=rpart_o(3,k)+sig*(random(3*k  )-0.5d0)
!     write(*,*)k,rr(1,k),rr(2,k),rr(3,k)
  end do

  allocate(cwf0(nspin,niso))


  call deut_wave(rr,cwf0,ysol,dr,r)

  !write(*,*)'here'
!!!Here you have to insert the subroutine that compute
  !the wave function (output cwf)
  !write(*,*) 'got to wave function calc'
  !computation of |psi|^2
  norm=0.d0
  do k=1,nspin
     do j=1,niso
        norm=norm+conjg(cwf0(k,j))*cwf0(k,j)
     end do
  end do

  !Metropolis step
  acc=.false.
  if(norm.gt.norm0*random(3*msg%npart+2))then
     norm0=norm
     rpart_o=rr
     cwf=cwf0
     acc=.true.
  end if


  return
end subroutine step


subroutine rndnb( nrnd , rnd0 , rnd )
!_______________________________________________________________________
!     This subroutine computes random number between 0 and 1
!     Input:
!     nrnd = number of random position
!     rnd0 = seed
!     Output:
!     rnd (nrnd)= random number generated      

  implicit real*8(a-h,o-z)
  save
  
  dimension rnd(nrnd) 
  
  data rndl / 16807.d0 / , rndm / 2147483647.d0 / &
       rndp / 2147483648.d0 /
  
  do  irnd=1,nrnd
     rnd0=mod( rndl*rnd0 , rndm )
     rnd(irnd)=rnd0/rndp
  end do
   
  return
end subroutine rndnb
