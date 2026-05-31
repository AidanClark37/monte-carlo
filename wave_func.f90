subroutine wave(func0)
 
  use mpi_modules
  use struct_funcs
  use interpolation
  use param_calc
  implicit none
  real*8,intent(out)::func0
  real*8::dr(3),rr(3,2),h,b5,cns
  real*8,parameter::gamma=4.5d0
  integer::i,j,k,li,lj,M,M0,ri,npoints,ias
  real*8::asum(4,2,2),r,plaguer,thetax,phix,x,prod,lr,lsum(2),f_prime_lambda
  real*8,allocatable::af(:)
  complex*16::ysum(4,2,2),ylm

  npoints=500
   h=10.d0/(npoints-1)
   allocate(af(npoints))

  ias=1
  do i = 1,npoints

     af(i)=0
     lr=real(i-1)*h


     x=gamma*lr

open(unit=25,file="func0.txt")
    
   do k=1,2
     lsum(k)=0.d0
     do li=1,40
        lj=li-1
        lsum(k)=lsum(k)+dsqrt((gamma**3)/((li)*(li+1)))*plaguer(lj,2.d0,x)*msg%wf(k,li)*EXP(-x/2.d0)
        
     enddo
     !lsum(k)=lsum(k)*f_lambda(r/msg%hbarc)*msg%lambda
     write(*,*)lr,gamma,plaguer(lj,2.d0,x),msg%wf(k,1),EXP(-x/2.d0),lsum(k)
  enddo
  call interpolate(fl1arr,lr/msg%hbarc,cns)
     af(i)=0.25d0*(lr**3)*(cns/msg%hbarc)*msg%lambda*(lsum(2)-dsqrt(2.d0)*lsum(1))**2.d0
     write(*,*)"f",cns,"l",msg%lambda,"1,2",lsum(1),lsum(2),"af",af(i)
     write(25,*)af(i)
enddo

write(*,*)"lr","gamma", "plaguer","wf","exp"
!ba(:)=af(:,1)
func0= b5(1,npoints,0.d0,0.d0,h,0.d0,0.d0,npoints,af,ias)
!ba(:)=af(:,2)
!func2=b5(1,npoints,0.d0,0.d0,h,0.d0,0.d0,npoints,ba,ias)

write(*,*)"func0",func0
close(25)

return
end subroutine wave
