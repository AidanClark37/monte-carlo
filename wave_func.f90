subroutine wave(func0)
 
  use mpi_modules
  use struct_funcs
  use interpolation
  use param_calc
  implicit none
  real*8,intent(out)::func0
  real*8::dr(3),rr(3,2),h,b5,cns,f_lambda
  real*8,parameter::gamma=4.5d0
  integer::i,j,k,li,lj,M,M0,ri,npoints,ias
  real*8::asum(4,2,2),r,plaguer,thetax,phix,x,prod,lr,lsum(2),f_prime_lambda
  real*8,allocatable::af(:)
  complex*16::ysum(4,2,2),ylm
  !number of points for b5
  npoints=2000
  h=60.d0/(npoints-1)
  allocate(af(npoints))
  
  ias=1
  do i = 1,npoints
     !lr - r value 
     af(i)=0.d0
     lr=real(i-1)*h

     !scale, dimensionless integral
     x=gamma*lr
     !constructing wavefunction
     ! k - L = 0,2
     !li - laguerre polynomial basis, corresponding index
   do k=1,2
     lsum(k)=0.d0
     do li=1,40
        lj=li-1
        lsum(k)=lsum(k)+dsqrt((gamma**3)/((li)*(li+1)))*plaguer(lj,2.d0,x)*msg%wf(k,li)*EXP(-x/2.d0)
      enddo
   enddo
   ! (2/3)r^3(f_0 + 1/sqrt2*f_2)^2
      af(i)=(2.d0/3.d0)*((lr**3)/msg%hbarc)*(lsum(1)+(1/dsqrt(2.d0))*lsum(2))**2.d0
enddo

func0= b5(1,npoints,0.d0,0.d0,h,0.d0,0.d0,npoints,af,ias)

return
end subroutine wave
