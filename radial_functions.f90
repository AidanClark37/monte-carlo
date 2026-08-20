!====================================
!k-space functions
!====================================

real*8 function Ck(k)
  use mpi_modules
    implicit none
    real*8::k
    Ck = exp(-(k/msg%lambda)**4)
  return
end function Ck

real*8 function sk(k)
  use mpi_modules
  implicit none
  real*8,intent(in)::k
  sk=sqrt(k**2+4*msg%mass**2)
  return
end function sk

real*8 function Lk(k)
  implicit none
  real*8,intent(in)::k
  real*8::sk
  if (k.eq.0)then
     Lk=1
  else
  Lk=0.5d0*(sk(k)/k)*log((sk(k)+k)/(sk(k)-k))
  endif
  return
end function Lk

real*8 function Hk(k)
  use mpi_modules
  implicit none
  real*8,intent(in)::k
  real*8::sk,Lk
  Hk=msg%mass**2/sk(k)**2*Lk(k)
  return
end function Hk

!====================================
!r-space functions
!====================================

real*8 function jbessel(sel,x)
  implicit none
  integer,intent(in)::sel
  real*8,intent(in)::x
  real*8::bessel
  select case(sel)
  case(0)
     jbessel=bessel(0,x)
  case(1)
     jbessel=-bessel(1,x)
  case(2)
     jbessel=(1.d0/3.d0)*(2*bessel(2,x)-bessel(0,x))
  case(3)
     jbessel = (1.d0/5.d0)*(2*bessel(3,x)-3*bessel(1,x))
  end select
  return
end function jbessel

!------------------------------------
! f lambda
!------------------------------------

real*8 function f_lambda(p,r)
  use mpi_modules
  implicit none
  integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,jbessel,b5,Ck,h
  integer::ias
  integer,intent(in)::p
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla
       
     af(i)=0
     lkk=real(i)*h

     af(i)=(lkk**(2+p))*Ck(lkk)*jbessel(p,lkk*r)/(lkk**2+msg%mass**2)

     
  enddo
  f_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  f_lambda=f_lambda/(2.d0*msg%lambda*acos(-1.d0)**2)
  return
end function f_lambda



real*8 function f_prime_lambda(r)
use mpi_modules
  implicit none
integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,jbessel,b5,Ck,h
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla


     af(i)=0
     lkk=real(i)*h
     af(i)=(lkk**3)*Ck(lkk)*jbessel(1,lkk*r)/(lkk**2+msg%mass**2)


  enddo
  f_prime_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  f_prime_lambda=f_prime_lambda/(2.d0*msg%lambda*acos(-1.d0)**2)
  return
end function f_prime_lambda

real*8 function f_2prime_lambda(r)
use mpi_modules
  implicit none
integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,jbessel,b5,Ck,h
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla
     af(i)=0
     lkk=real(i)*h
     af(i)=(lkk**4)*Ck(lkk)*jbessel(2,lkk*r)/(lkk**2+msg%mass**2)
     
  enddo
  f_2prime_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  f_2prime_lambda=f_2prime_lambda/(2.d0*msg%lambda*acos(-1.d0)**2)
  return
end function f_2prime_lambda


!------------------------------------
!g lambda
!------------------------------------
real*8 function g_lambda(p,r)
use mpi_modules
  implicit none
integer::i
real*8,allocatable::af(:)
integer,intent(in)::p
  real*8::r,lkk,jbessel,b5,Ck,h
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla
       
     
     af(i)=0
     lkk=real(i)*h
     af(i)=(lkk**(2-p))*Ck(lkk)*jbessel(p,lkk*r)/(lkk**2+msg%mass**2)**2

     
  enddo
  g_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  g_lambda=g_lambda*msg%lambda/(2.d0*acos(-1.d0)**2)
  return
end function g_lambda



real*8 function g_prime_lambda(r)
use mpi_modules
  implicit none
integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,bessel,b5,Ck,h
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla


     af(i)=0
     lkk=real(i)*h
     af(i)=(lkk**3)*Ck(lkk)*bessel(1,lkk*r)/(lkk**2+msg%mass**2)**2


  enddo
  g_prime_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  g_prime_lambda=-g_prime_lambda*msg%lambda/(2.d0*acos(-1.d0)**2)
  return
end function g_prime_lambda

real*8 function g_2prime_lambda(r)
use mpi_modules
  implicit none
integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,bessel,b5,Ck,h
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla
     af(i)=0
     lkk=real(i)*h
     af(i)=(lkk**4)*Ck(lkk)*(2*bessel(2,lkk*r)-bessel(0,lkk*r))/(lkk**2+msg%mass**2)**2
  enddo
  g_2prime_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  g_2prime_lambda=g_2prime_lambda*msg%lambda/(6.d0*acos(-1.d0)**2)
  return
end function g_2prime_lambda

!-----------------------------------
!L_lambda
!-----------------------------------

real*8 function L_lambda(p,r)
  use mpi_modules
  implicit none
  integer::i,ias
  integer,intent(in)::p
  real*8,allocatable::af(:)
  real*8::r,lkk,b5,jbessel,Ck,h,Lk
  ias = 1
  allocate(af(msg%nla))
  h = 2.d0*msg%lambda/msg%nla
  do i=1,msg%nla
     af(i)=0
     lkk=h*real(i)
     af(i)=(lkk**(2+p))*Ck(lkk)*jbessel(p,lkk*r)/(lkk**2+msg%mass**2)*Lk(lkk)/msg%mass**2
  enddo
  L_lambda=b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  L_lambda=L_lambda*msg%lambda/(2.d0*acos(-1.d0)**2)
  return
end function L_lambda

real*8 function L_prime_lambda(r)
  use mpi_modules
  implicit none
  integer::i,ias
  real*8,allocatable::af(:)
  real*8::r,lkk,b5,bessel,Ck,h,Lk
  ias = 1
  allocate(af(msg%nla))
  h = 2.d0*msg%lambda/msg%nla
  do i=1,msg%nla
     af(i)=0
     lkk=h*real(i)
     af(i)=(lkk**3)*Ck(lkk)*bessel(1,lkk*r)/(lkk**2+msg%mass**2)*Lk(lkk)/msg%mass**2
  enddo
  L_prime_lambda=b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  L_prime_lambda=-L_prime_lambda*msg%lambda/(2.d0*acos(-1.d0)**2)
  return
end function L_prime_lambda

real*8 function L_2prime_lambda(r)
use mpi_modules
  implicit none
integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,bessel,b5,Ck,h,Lk
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla
     af(i)=0
     lkk=real(i)*h
     af(i)=(lkk**4)*Ck(lkk)*(2*bessel(2,lkk*r)-bessel(0,lkk*r))/(lkk**2+msg%mass**2)*Lk(lkk)/msg%mass**2
  enddo
  L_2prime_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  L_2prime_lambda=L_2prime_lambda*msg%lambda/(6.d0*acos(-1.d0)**2)
  return
end function L_2prime_lambda



!------------------------------------
!C lambda
!------------------------------------

real*8 function C_lambda(p,r)
  use mpi_modules
  implicit none
  integer::i
  integer,intent(in)::p
  real*8,allocatable::af(:)
  real*8::r,lkk,jbessel,b5,Ck,h
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla

     af(i)=0
     lkk=real(i)*h

     af(i)=(lkk**(2+p))*Ck(lkk)*jbessel(p,lkk*r)


  enddo
  C_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  C_lambda=C_lambda/((msg%lambda**3)*2.d0*acos(-1.d0)**2)
  return
end function C_lambda


real*8 function C_prime_lambda(r)
  use mpi_modules
  implicit none
  integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,bessel,b5,Ck,h
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla

     af(i)=0
     lkk=real(i)*h

     af(i)=-(lkk**3)*Ck(lkk)*bessel(1,lkk*r)


  enddo
  C_prime_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  C_prime_lambda=C_prime_lambda/((msg%lambda**3)*2.d0*acos(-1.d0)**2)
  return
end function C_prime_lambda


real*8 function C_2prime_lambda(r)
use mpi_modules
  implicit none
integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,bessel,b5,Ck,h,Lk
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla
     af(i)=0
     lkk=real(i)*h
     af(i)=(lkk**4)*Ck(lkk)*(2*bessel(2,lkk*r)-bessel(0,lkk*r))
  enddo
  C_2prime_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  C_2prime_lambda=C_2prime_lambda/((msg%lambda**3)*6.d0*acos(-1.d0)**2)
  return
end function C_2prime_lambda


!------------------------------------ 
!H lambda
!------------------------------------


real*8 function H_lambda(p,r)
  use mpi_modules
  implicit none
  integer::i,ias
  integer,intent(in)::p
  real*8,allocatable::af(:)
  real*8::r,lkk,b5,jbessel,Hk,h,Lk,Ck
  ias = 1
  allocate(af(msg%nla))
  h = 2.d0*msg%lambda/msg%nla
  do i=1,msg%nla
     af(i)=0
     lkk=h*real(i)
     af(i)=(lkk**(2+p))*Ck(lkk)*jbessel(p,lkk*r)/(lkk**2+msg%mass**2)*Hk(lkk)/msg%mass**2
  enddo
  H_lambda=b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  H_lambda=H_lambda*msg%lambda/(2.d0*acos(-1.d0)**2)
  return
end function H_lambda


real*8 function H_prime_lambda(r)
  use mpi_modules
  implicit none
  integer::i,ias
  real*8,allocatable::af(:)
  real*8::r,lkk,b5,bessel,Ck,h,Hk
  ias = 1
  allocate(af(msg%nla))
  h = 2.d0*msg%lambda/msg%nla
  do i=1,msg%nla
     af(i)=0
     lkk=h*real(i)
     af(i)=(lkk**3)*Ck(lkk)*(bessel(1,lkk*r)/(lkk**2+msg%mass**2))*Hk(lkk)/(msg%mass**2)
  enddo
  H_prime_lambda=b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  H_prime_lambda=-H_prime_lambda*msg%lambda/(2.d0*acos(-1.d0)**2)
  return
end function H_prime_lambda


real*8 function H_2prime_lambda(r)
use mpi_modules
  implicit none
integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,bessel,b5,Ck,h,Hk
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla
     af(i)=0
     lkk=real(i)*h
     af(i)=(lkk**4)*Ck(lkk)*(2*bessel(2,lkk*r)-bessel(0,lkk*r))/(lkk**2+msg%mass**2)*Hk(lkk)/msg%mass**2
  enddo
  H_2prime_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  H_2prime_lambda=H_2prime_lambda*msg%lambda/(6.d0*acos(-1.d0)**2)
  return
end function H_2prime_lambda

!------------------------------------
!G lambda
!------------------------------------

real*8 function G1_lambda(p,r)
  use mpi_modules
  implicit none
  integer::i,ias
  integer,intent(in)::p
  real*8,allocatable::af(:)
  real*8::r,lkk,b5,jbessel,Hk,h,Ck,sk
  ias = 1
  allocate(af(msg%nla))
  h = 2.d0*msg%lambda/msg%nla
  do i=1,msg%nla
     af(i)=0
     lkk=h*real(i)
     af(i)=(lkk**(2+p))*Ck(lkk)*jbessel(p,lkk*r)/(lkk**2+msg%mass**2)
     af(i)=af(i)*(4*Hk(lkk)*(1/msg%mass**2+1/sk(lkk)**2)+1/sk(lkk)**2)
  enddo
  G1_lambda=b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  G1_lambda=G1_lambda*msg%lambda/(2.d0*acos(-1.d0)**2)
  return
end function G1_lambda

real*8 function G1_prime_lambda(r)
  use mpi_modules
  implicit none
  integer::i,ias
  real*8,allocatable::af(:)
  real*8::r,lkk,b5,bessel,Ck,h,Hk,sk
  ias = 1
  allocate(af(msg%nla))
  h = 2.d0*msg%lambda/msg%nla
  do i=1,msg%nla
     af(i)=0
     lkk=h*real(i)
     af(i)=(lkk**3)*Ck(lkk)*bessel(1,lkk*r)/(lkk**2+msg%mass**2)
     af(i)=af(i)*(4*Hk(lkk)*(1/msg%mass**2+1/sk(lkk)**2)+1/sk(lkk)**2)
  enddo
  G1_prime_lambda=b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  G1_prime_lambda=-G1_prime_lambda*msg%lambda/(2.d0*acos(-1.d0)**2)
  return
end function G1_prime_lambda


real*8 function G1_2prime_lambda(r)
use mpi_modules
  implicit none
integer::i
  real*8,allocatable::af(:)
  real*8::r,lkk,bessel,b5,Ck,h,Hk,sk
  integer::ias
  h=2.d0*msg%lambda/msg%nla
  allocate(af(msg%nla))
  ias=1
  do i = 1,msg%nla
     af(i)=0
     lkk=real(i)*h
     af(i)=(lkk**4)*Ck(lkk)*(2*bessel(2,lkk*r)-bessel(0,lkk*r))/(lkk**2+msg%mass**2)
     af(i)=af(i)*(4*Hk(lkk)*(1/msg%mass**2+1/sk(lkk)**2)+1/sk(lkk)**2)
  enddo
  G1_2prime_lambda = b5(1,msg%nla,0.d0,0.d0,h,0.d0,0.d0,msg%nla,af,ias)
  G1_2prime_lambda=G1_2prime_lambda*msg%lambda/(6.d0*acos(-1.d0)**2)
  return
end function G1_2prime_lambda
