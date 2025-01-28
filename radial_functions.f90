  real*8 function C_lambda(k,lambda)
    implicit none
    real*8::k,lambda
    C_lambda = exp(-(k/lambda)**4)
  return
end function C_lambda


  
  real*8 function sbessel(x)
     implicit none
  real*8::x
  if( x.eq.0.d0) then
  sbessel=1
  else
     sbessel=sin(x)/x
     endif
  return

  end function sbessel

  real*8 function f_0(k,R,lambda,m)
  implicit none
  real*8::k,r,m,sbessel,lambda,C_lambda
  f_0 = 1/(2.d0*(3.14159265358979d0**2)*lambda)
  f_0 = C_lambda(k,lambda)*f_0*(sbessel(k*r/197.326))/(k**2+m**2)
  return
end function f_0

  
  real*8 function f_lambda(R,m,lambda,nla)
implicit none
  integer::i,nla
  real*8,allocatable::wei(:),xx(:)
  real*8::R,m,lambda,f_0
  allocate(xx(nla))
  allocate(wei(nla))
  f_lambda=0.d0
  call setgaulag(2.d0,nla,wei,xx)
  do i = 1,nla
     f_lambda = f_lambda+exp(xx(i))*f_0(xx(i),R,lambda,m)*wei(i)
  enddo
  return
end function f_lambda

real*8 function ff_0(k,R,q,y,lambda,m)
  implicit none
  real*8::k,R,q,m,y,sbessel,lambda,C_lambda
  ff_0 = 1/(2.d0*(3.14159265358979d0**2)*lambda)
  ff_0 = C_lambda(k,lambda)*ff_0*(sbessel(k*r/197.326))/(k**2+m**2+(q**2)*(y-1)/4)
  return
end function ff_0

real*8 function ff_lambda(R,q,y,m,lambda)
  implicit none
  integer::i,nla
  real*8,allocatable::wei(:),xx(:)
  real*8::R,q,m,y,lambda,ff_0
  allocate(xx(nla))
  allocate(wei(nla))
  ff_lambda=0.d0
  call setgaulag(2.d0,nla,wei,xx)
  do i = 1,nla
     ff_lambda = ff_lambda+exp(xx(i))*ff_0(xx(i),R,q,y,lambda,m)*wei(i)
  enddo



end function ff_lambda  !nla=181
!  m=138.03919333333d0
 ! lambda=600.d0
 









