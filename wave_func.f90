
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
  npoints=10000
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
   ! r^3(f_0 + 1/sqrt2*f_2)^2
      af(i)=((lr**3)/msg%hbarc)*(lsum(1)+(1/dsqrt(2.d0))*lsum(2))**2.d0
enddo

func0= b5(1,npoints,0.d0,0.d0,h,0.d0,0.d0,npoints,af,ias)

return
end subroutine wave


subroutine wave2(func0)
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
  npoints=10000
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
      af(i)=(lsum(1)**2+1/(5*dsqrt(2.d0))*lsum(1)*lsum(2)-(2.d0/5.d0)*lsum(2)**2)
enddo

func0= b5(1,npoints,0.d0,0.d0,h,0.d0,0.d0,npoints,af,ias)

return
end subroutine wave2


subroutine wave3()
  use wigner
  use sim_ops
  use pre_deut
  use param_calc
  use spin_ops_act
  implicit none
  real*8::ysol(2,40)
  real*8,parameter::gamma=4.5d0
  integer::p,si,sj,i,j,k,li,lj,M,M0,ri,l1,l2,npoints,ias
  real*8::sum(4,2,2),r,plaguer,thetax,phix,x,h,lt,lf,b5
  complex*16,allocatable::af(:,:,:,:),af2(:,:,:),yaf(:,:),yaf2(:)
  complex*16::ysum(4,2,2),ylm,prod(4,2,2),sigma_wf(4,2,2),cc_wf(4,2,2),yint
  complex*16::func0(2,2)
  ! r = r_1-r_2 vector, labelled dr                                             npoints=10000
  
  
         npoints=1000
  h=2.d0/(npoints-1)
  allocate(af(npoints,npoints,2,2))
  allocate(af2(npoints,2,2))
  allocate(yaf(npoints,npoints))
  allocate(yaf2(npoints))
  ias=1
  do i = 1,npoints
     do j = 1,npoints
     !lr - r value                                                              
        af(i,j,:,:)=dcmplx(0.d0,0.d0)
        af2(j,:,:)=dcmplx(0.d0,0.d0)
        lt=real(i-1)*h-1.d0
     lf = real(j-1)*h*dacos(-1.d0)
  ! |r|, labeled r     
     yaf(i,j)=1!conjg(ylm(0,0,lt,lf))*ylm(0,0,lt,lf)
     !write(*,*)ylm(0,0,lt,lf)
  !radial wavefunction, f_0 and f_2 part                                        
  

     ! dependance on sperical harmonics                                            
  prod=0.d0   
     
        
do k = 1,2
           do si=1,4

     M0=2-CEILING((si-1.d0)/2)
     
	
           ysum(si,:,k) = warray(si,:,k)*ylm(2*k-2,M0,lt,lf)
        
        
                  
        enddo
!        write(*,*)"------------"
!do si=1,4
!   write(*,*)warray(si,1,1),warray(si,2,1)
!enddo
!write(*,*)"--------------"
!do si=1,4
do p = 1,2


   call spin(ysum(:,:,k),p,1,sigma_wf(:,:,k))
   !write(*,*)"zero!",sigma_wf(4,1,2)
   
   prod(:,:,k)=prod(:,:,k)+dsqrt(1-lt**2)*cos(lf)*sigma_wf(:,:,k)
   
   call spin(ysum(:,:,k),p,2,sigma_wf(:,:,k))

   prod(:,:,k)=prod(:,:,k)+dsqrt(1-lt**2)*sin(lf)*sigma_wf(:,:,k)

   call spin(ysum(:,:,k),p,3,sigma_wf(:,:,k))

   prod(:,:,k)=prod(:,:,k)+lt*sigma_wf(:,:,k)
  
  enddo
enddo
  prod = prod*lt

  cc_wf=dconjg(ysum)
!  write(*,*)prod(1,1,1)
  do l1=1,2
     do l2=1,2

        

        do si = 1,4
           do sj=1,2
              af(i,j,l1,l2)=af(i,j,l1,l2)+cc_wf(si,sj,l1)*ysum(si,sj,l2)
           enddo
        enddo
        
     enddo
  enddo
  
enddo
enddo

do j = 1,npoints
   yaf2(j)=b5(1,npoints,0.d0,0.d0,h,0.d0,0.d0,npoints,yaf(j,:),ias)
   write(*,*) yaf2(j)
enddo

yint = b5(1,npoints,0.d0,0.d0,h,0.d0,0.d0,npoints,yaf2,ias)
write(*,*)"y00", yint
stop
do l1=1,2
   do l2=1,2
      do j=1,npoints
      af2(j,l1,l2)= b5(1,npoints,0.d0,0.d0,dacos(-1.d0)*h,0.d0,0.d0,npoints,af(:,j,l1,l2),ias)
enddo
func0(l1,l2)=       b5(1,npoints,0.d0,0.d0,h,0.d0,0.d0,npoints,af2(:,l1,l2),ias)
write(*,*) func0(l1,l2)
   enddo
enddo
end subroutine wave3
