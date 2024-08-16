subroutine step(rnd,sigma,rpart_o,npart,jz,cwf,norm0,acc)
  implicit none

  real*8             ,intent(inout)::rnd
  real*8             ,intent(inout)::norm0
  real*8             ,intent(in)::sigma
  integer            ,intent(in)::npart
  integer            ,intent(in)::jz
  real*8             ,intent(inout)::rpart_o(3,npart)
  complex*16,intent(out)::cwf(nspin,niso)
  logical,intent(out)::acc
  
  real*8::rr(3,npart)
  real*8::random(3*npart+2)
  real*8::norm
  integer::k,j
  
  call rndnb(3*npart+2 , rnd , random )!call the random number

  !selection of the new position configuration
  do k=1,npart
     rr(1,k)=rpart_o(1,k)+sigma*(random(3*k-2)-0.5d0)
     rr(2,k)=rpart_o(2,k)+sigma*(random(3*k-1)-0.5d0)
     rr(3,k)=rpart_o(3,k)+sigma*(random(3*k  )-0.5d0)
!     write(*,*)k,rr(1,k),rr(2,k),rr(3,k)
  end do

!!!Here you have to insert the subroutine that compute
  !the wave function (output cwf)

  !computation of |psi|^2
  norm=0.d0
  do k=1,nspin
     do j=1,niso
        norm=norm+conjg(cwf(k,j))*cwf(k,j)
     end do
  end do

  !Metropolis step
  acc=.false.
  if(norm.gt.norm0*random(3*wt%npart+2))then
     norm0=norm
     rpart_o=rr
     acc=.true.
  end if


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
