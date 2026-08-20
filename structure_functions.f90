module struct_funcs
  implicit none
  real*8,allocatable::flarr(:,:),Clarr(:,:),glarr(:,:),Llarr(:,:),Hlarr(:,:),G1larr(:,:)
  !real*8,allocatable::fl1arr(:),Cl1arr(:),gl1arr(:),Ll1arr(:),Hl1arr(:),G1l1arr(:)
  !real*8,allocatable::gl2arr(:),G1l2arr(:)
  real*8,allocatable::xarray(:)
contains
  subroutine struct_func_calc(xmin,xmax)
    use mpi_modules
    use param_calc
    implicit none
    real*8,intent(in)::xmin,xmax
    real*8::f_lambda,f_prime_lambda,C_lambda,C_prime_lambda,g_lambda,g_prime_lambda
    real*8::g_2prime_lambda,L_lambda,L_prime_lambda,H_lambda,H_prime_lambda
    real*8::G1_lambda,G1_prime_lambda,G1_2prime_lambda
    real*8::xrange,x
    integer::i,p
    xrange=xmax-xmin
    allocate(flarr(msg%npoint,4))
    !allocate(fl1arr(msg%npoint,4))
    allocate(Clarr(msg%npoint,4))
    !allocate(Cl1arr(msg%npoint,4))
    allocate(glarr(msg%npoint,4))
    !allocate(gl1arr(msg%npoint))
    !allocate(gl2arr(msg%npoint))
    allocate(Llarr(msg%npoint,4))
!    allocate(Ll1arr(msg%npoint))
    allocate(Hlarr(msg%npoint,4))
!    allocate(Hl1arr(msg%npoint))
    allocate(G1larr(msg%npoint,4))
!    allocate(G1l2arr(msg%npoint))
    
    
    allocate(xarray(msg%npoint))
    do p = 1,4
    do i=1,msg%npoint
       if(p.eq.1) xarray(i)=real(i-1)*xrange/(msg%npoint-1)+xmin
       flarr(i,p)  =f_lambda        (p-1,xarray(i))
       !fl1arr(i) =f_prime_lambda  (xarray(i))
       Clarr(i,p)  =C_lambda        (p-1,xarray(i))
       !Cl1arr(i) =C_prime_lambda  (xarray(i))
       glarr(i,p)  =g_lambda        (p-1,xarray(i))
       !gl1arr(i) =g_prime_lambda  (xarray(i))
       !gl2arr(i) =g_2prime_lambda (xarray(i))
       Llarr(i,p)  =L_lambda        (p-1,xarray(i))
       !Ll1arr(i) =L_prime_lambda  (xarray(i))
       Hlarr(i,p)  =H_lambda        (p-1,xarray(i))
       !Hl1arr(i) =H_prime_lambda  (xarray(i))
       G1larr(i,p)=G1_lambda (p-1,xarray(i))
       !G1l2arr(i)=G1_2prime_lambda(xarray(i))
    enddo
    enddo
    write(*,*)"module:fl1 filled",flarr(1,2),flarr(msg%npoint,2)
    return
   end subroutine struct_func_calc

 end module struct_funcs
 
