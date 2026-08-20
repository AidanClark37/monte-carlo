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
  real*8             ,intent(inout)::dr(3),r
  real*8::sig
  complex*16::warray(4,2)
  complex*16,allocatable::cwf0(:,:)
  real*8::ysol(2,40),trial_dr(3),trial_r
  logical,intent(out)::acc

  real*8,intent(inout)::rr(3,2)
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


  call deut_wave(rr,cwf0,ysol,trial_dr,trial_r)

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
  !on rejection, rr/dr/r/cwf must all stay tied to the last
  !accepted configuration (rpart_o) rather than the rejected trial
  acc=.false.
  if(norm.gt.norm0*random(3*msg%npart+2))then
     norm0=norm
     rpart_o=rr
     cwf=cwf0
     dr=trial_dr
     r=trial_r
     acc=.true.
  end if
  rr=rpart_o

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

subroutine gauss_init(rnd,rpart_o,npart,sig)
!_______________________________________________________________________
!     Draws an initial particle configuration from a Gaussian,
!     mean 0, standard deviation sig, independently for each
!     coordinate of each particle (Box-Muller on top of rndnb).
!     rnd is the walker's own RNG seed, threaded through so each
!     walker gets an independent draw.

  implicit none
  real*8,intent(inout)::rnd
  integer,intent(in)::npart
  real*8,intent(in)::sig
  real*8,intent(out)::rpart_o(3,npart)
  real*8,allocatable::u(:),z(:)
  real*8::pi,rad
  integer::n,m

  pi=dacos(-1.d0)
  n=3*npart
  if(mod(n,2).ne.0) n=n+1
  allocate(u(n),z(n))

  call rndnb(n,rnd,u)

  do m=1,n/2
     rad=dsqrt(-2.d0*dlog(u(2*m-1)))
     z(2*m-1)=rad*dcos(2.d0*pi*u(2*m))
     z(2*m)  =rad*dsin(2.d0*pi*u(2*m))
  end do

  rpart_o=sig*reshape(z(1:3*npart),(/3,npart/))

  return
end subroutine gauss_init
