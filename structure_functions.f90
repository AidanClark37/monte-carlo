module struct_funcs
  implicit none
  real*8,allocatable::flarr(:),Clarr(:),glarr(:),Llarr(:),Hlarr(:)
  real*8,allocatable::fl1arr(:),Cl1arr(:),gl1arr(:),Ll1arr(:),Hl1arr(:),G1l1arr(:)
  real*8,allocatable::gl2arr(:),G1l2arr(:)
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
    integer::i
    xrange=xmax-xmin
    allocate(flarr(msg%npoint))
    allocate(fl1arr(msg%npoint))
    allocate(Clarr(msg%npoint))
    allocate(Cl1arr(msg%npoint))
    allocate(glarr(msg%npoint))
    allocate(gl1arr(msg%npoint))
    allocate(gl2arr(msg%npoint))
    allocate(Llarr(msg%npoint))
    allocate(Ll1arr(msg%npoint))
    allocate(Hlarr(msg%npoint))
    allocate(Hl1arr(msg%npoint))
    allocate(G1l1arr(msg%npoint))
    allocate(G1l2arr(msg%npoint))
    
    
    allocate(xarray(msg%npoint))
    do i=1,msg%npoint
       xarray(i)=real(i-1)*xrange/(msg%npoint-1)+xmin
       flarr(i)  =f_lambda        (xarray(i))
       fl1arr(i) =f_prime_lambda  (xarray(i))
       Clarr(i)  =C_lambda        (xarray(i))
       Cl1arr(i) =C_prime_lambda  (xarray(i))
       glarr(i)  =g_lambda        (xarray(i))
       gl1arr(i) =g_prime_lambda  (xarray(i))
       gl2arr(i) =g_2prime_lambda (xarray(i))
       Llarr(i)  =L_lambda        (xarray(i))
       Ll1arr(i) =L_prime_lambda  (xarray(i))
       Hlarr(i)  =H_lambda        (xarray(i))
       Hl1arr(i) =H_prime_lambda  (xarray(i))
       G1l1arr(i)=G1_prime_lambda (xarray(i))
       G1l2arr(i)=G1_2prime_lambda(xarray(i))
    enddo
    write(*,*)"module:fl1 filled",fl1arr(1),fl1arr(595)
    return
   end subroutine struct_func_calc

 end module struct_funcs
 
