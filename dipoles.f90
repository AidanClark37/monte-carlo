!-------------
!dipoles b1, c
!-------------

subroutine  g1V(r,dr,wf_in,index,exp)
  use operator_calc
  use mpi_modules
  use param_calc
  use tau_operator
  use struct_funcs
  use interpolation
  implicit none
  complex*16,intent(in)::wf_in(4,2)
  integer,intent(in)::index
  real*8,intent(in)::r,dr(3)
  real*8,intent(out)::exp
  complex*16::t1wf(nspin,niso),t2wf(nspin,niso),t3wf(nspin,niso),spinwf(nspin,niso),ccwf(nspin,niso),wfout(nspin,niso)
  integer::i,p,j,si
  real*8::f_lambda,cns,x,dx(3)
  exp=0.d0
  x=r/msg%hbarc
  dx(:)=dr(:)/msg%hbarc
  spinwf(:,:)=dcmplx(0.d0,0.d0)
  !calc S_+ . rhat
  do i=1,3
     do p=1,2
        spinwf(:,:) = spinwf(:,:)+s_wf(:,:,p,i)*dx(i)/x
     enddo
  enddo
  
  !multiply by z
  spinwf(:,:)=spinwf(:,:)*dx(index)
  
  
  
  
  
  
  !calc tau operators tau_1 . tau_2 and tau_1^z tau_2^z
  ! call tau0(spinwf,1,2,t1wf)
  ! call tau2(spinwf,1,2,t2wf)
  
  
  !calc structure function
  
  !     call interpolate(flarr,1,x,cns)
  !     cns=msg%lambda*cns/4.d0
  
  !     cns=msg%lambda*f_lambda(1,x)/4.d0
  
  wfout(:,:)=spinwf(:,:)
  !wfout(i,j)=(t1wf(:,:)-t2wf(:,:))*cns
  ccwf=conjg(wf_in)
  exp=0.d0
  ! inner product
  do i = 1,nspin
     do j=1,niso
        exp=exp+wfout(i,j)*ccwf(i,j)
     enddo
  enddo
  return
end subroutine g1V

!--------
!dipole d
!--------

subroutine d0(wf_in,index,exp)
  use operator_calc
  use mpi_modules
  use param_calc
  implicit none
  complex*16,intent(in)::wf_in(4,2)
  integer,intent(in)::index
  real*8,intent(out)::exp
  complex*16::spinwf(4,2),ccwf(4,2)
  integer::i,p,j
  
  spinwf(:,:)=dcmplx(0.d0,0.d0)
  exp=0.d0
  !calc S_+^i 
  do p=1,2
     spinwf(:,:) = spinwf(:,:)+s_wf(:,:,p,index)
  enddo
  
  spinwf(:,:)=-0.5d0*spinwf(:,:)
  ccwf=conjg(wf_in)
  !inner product
  do i = 1,nspin
     do j=1,niso
        
	exp=exp+spinwf(i,j)*ccwf(i,j)
 
enddo
enddo

return
end subroutine d0

!--------
!dipole d
!--------

subroutine d1(wf_in,index,exp)
  use operator_calc
  use mpi_modules
  use param_calc
  use tau_operator
  implicit none
  complex*16,intent(in)::wf_in(4,2)
  integer,intent(in)::index
  real*8,intent(out)::exp
  complex*16::spinwf(4,2),ccwf(4,2),twf(4,2),wf_out(4,2)
  integer::i,p,j
  wf_out(:,:)=dcmplx(0.d0,0.d0)
  twf(:,:)=dcmplx(0.d0,0.d0)
  spinwf(:,:)=dcmplx(0.d0,0.d0)
  exp=0.d0
  
  do p=1,2
!    write(*,*)"---------------------------"
!    write(*,*)"p=",p
!    write(*,*)"---------------------------"
!    write(*,*)"==========================="
!    do i=1,nspin
!       write(*,*)"wf", wf_in(i,1),wf_in(i,2)
!    enddo
!    write(*,*)"==========================="
     !operate sigma_p
     spinwf(:,:) = s_wf(:,:,p,index)
!    do i =1,nspin
!       write(*,*) "s1",spinwf(i,1),spinwf(i,2)
!    enddo
!    write(*,*)"==========================="
     !operate tau_p
     call tau1(spinwf,p,twf)
!    do i = 1,nspin
!       write(*,*) "t1",twf(i,1),twf(i,2)
!    enddo
!    write(*,*)"==========================="
     
     wf_out = wf_out - twf
  enddo

  
  ccwf=conjg(wf_in)
  !inner product                                                                                                        
  do i = 1,nspin
     do j=1,niso

        exp=exp+wf_out(i,j)*ccwf(i,j)

enddo
enddo
!write(*,*)"exp",exp

return
end subroutine d1

!----------------------
!dipoles a1, a2, b2, b3
!----------------------

subroutine g0(r,dr,wf_in,dwf,index,exp)
  use operator_calc
  use tau_operator
  use mpi_modules
  use param_calc
  implicit none
  complex*16,intent(in)::wf_in(:,:),dwf(:,:)
  integer,intent(in)::index
  real*8,intent(in)::r,dr(:)
  real*8,intent(out)::exp
  complex*16::swf_1(nspin,niso),swf_2(nspin,niso),swf0_1(nspin,niso)
  complex*16::swf0_2(nspin,niso),swf1_1(nspin,niso),swf1_2(niso,nspin)
  complex*16::owf(nspin,niso),ccwf(nspin,niso)
  integer::i,j
  real*8::cst,f_lambda
  swf_1(:,:)=s_wf(:,:,1,index)
  swf_2(:,:)=s_wf(:,:,2,index)
  cst=-gA*msg%lambda/(2*fpi*nmass)*f_lambda(0,r)
  call tau0(swf_1,1,2,swf0_1)
  call tau1(swf_1,2,swf1_1)
  call tau0(swf_2,1,2,swf0_2)
  call tau1(swf_2,1,swf1_2)
  owf(:,:)=cst*(swf0_1(:,:)+swf0_2(:,:)+swf1_1(:,:)+swf1_2(:,:))
  ccwf=conjg(wf_in)
  exp=0.d0
  do i = 1,nspin
     do j=1,niso
        exp=exp+ccwf(i,j)*owf(i,j)
     enddo
  enddo
  return
end subroutine g0

!---------
!dipole a2
!---------

subroutine g1(r,dr,wf_in,index,exp)
  use param_calc
  use tau_operator
  use operator_calc
  use mpi_modules
  use interpolation
  use struct_funcs
  implicit none
  real*8,intent(in)::r,dr(3)
  real*8,intent(out)::exp
  real*8::x,cns,f_lambda
  integer,intent(in)::index
  integer::p,i,j
  complex*16,intent(in)::wf_in(4,2)
  complex*16::spinwf(4,2),twf(4,2),wf_out(4,2),cc_wf(4,2)
  wf_out(:,:)=dcmplx(0.d0,0.d0)
x=r/msg%hbarc

do p = 1,msg%npart
     !operate sigma_p
     spinwf=s_wf(:,:,p,index)
     !operate tau_p
     !do i = 1,nspin
        
     !      write(*,*) "spin", spinwf(i,1),spinwf(i,2)
        
     !enddo
     
     call tau1(spinwf,p,twf)
     !do i = 1,nspin
        
     !      write(*,*) "isospin",twf(i,1),twf(i,2)
        
     !enddo


     ! add sigma_p(tau_p+1) to total acted wavefunction
     wf_out = wf_out + twf+spinwf
  enddo
  
  call interpolate(flarr,0,x,cns)                                                                                 
  cns=msg%lambda*cns
  !cns=f_lambda(0,x)*msg%lambda
  wf_out=wf_out*cns
  cc_wf=conjg(wf_in)
exp=0.d0
  do i = 1,nspin
     do j = 1,niso
        exp=exp+wf_out(i,j)*cc_wf(i,j)
     enddo
  enddo
  return

  
end  subroutine g1



!---------                                    
!dipole a2                                                                                                            
!---------                                                                                  

subroutine g2(r,dr,wf_in,index,exp)
  use param_calc
  use mpi_modules
  use tau_operator
  use operator_calc
  use interpolation
  use struct_funcs
  implicit none
  real*8,intent(in)::r,dr(3)
  real*8,intent(out)::exp
  real*8::cns,x,f_lambda
  integer,intent(in)::index
  integer::p,i,j
  complex*16,intent(in)::wf_in(4,2)
  complex*16::spinwf(4,2),t1wf(4,2),t2wf(4,2),wf_out(4,2),cc_wf(4,2)
  wf_out(:,:)=dcmplx(0.d0,0.d0)
  x=r/msg%hbarc
  do p = 1,msg%npart
     spinwf=s_wf(:,:,p,index)
     call tau1(spinwf,p,t1wf)
     call tau2(spinwf,1,2,t2wf)
     
     
     wf_out = wf_out + t1wf+t2wf
  
  enddo
  call interpolate(flarr,0,x,cns)
  cns=msg%lambda*cns
  !cns=f_lambda(0,x)*msg%lambda
  wf_out=wf_out*cns

  cc_wf=conjg(wf_in)
exp=0.d0
  do i = 1,nspin
     do	j = 1,niso
        exp=exp+wf_out(i,j)*cc_wf(i,j)
     enddo
  enddo

return
end  subroutine g2



!------------
!dipoles e, f
!------------

subroutine Delta(r,dr,wf_in,index,exp)
  implicit none
  real*8,intent(in)::r,dr(3)
  real*8,intent(out)::exp
  real*8::cns,x,dx(3)
  integer,intent(in)::index
  integer::p,i,j
  complex*16,intent(in)::wf_in(4,2)

  
end subroutine Delta

!------------------------------
!test operators
!------------------------------

!------------------------------
!S_+^i r^i
!------------------------------
subroutine  Srop(r,dr,index,wf_in,exp)
  use operator_calc
  use mpi_modules
  use param_calc
  use tau_operator
  use interpolation
  use struct_funcs
  implicit none
  complex*16,intent(in)::wf_in(4,2)
  integer,intent(in)::index
  real*8,intent(in)::r,dr(3)
  real*8,intent(out)::exp
  complex*16::spinwf(nspin,niso),ccwf(nspin,niso),wfout(nspin,niso)
  integer::i,p,j,si
  real*8::x,dx(3)
  exp=0
  x=r/msg%hbarc
  dx(:)=dr(:)/msg%hbarc
  spinwf(:,:)=dcmplx(0.d0,0.d0)
  
  do p=1,2
     spinwf(:,:) = spinwf(:,:)+s_wf(:,:,p,index)*dx(index)/x
     
  enddo
  
  ccwf=conjg(wf_in)
  exp=0.d0
  !inner product
  do i = 1,nspin
     do j=1,niso
        exp=exp+spinwf(i,j)*ccwf(i,j)
     enddo
  enddo
  return
end subroutine Srop

!-------------------------------
!z S_+^i r^i
!-------------------------------
subroutine  rSrop(r,dr,index,wf_in,exp)
  use operator_calc
  use mpi_modules
  use param_calc
  use tau_operator
  use struct_funcs
  use interpolation
  implicit none
  complex*16,intent(in)::wf_in(4,2)
  integer,intent(in)::index
  real*8,intent(in)::r,dr(3)
  real*8,intent(out)::exp
  complex*16::spinwf(nspin,niso),ccwf(nspin,niso),wfout(nspin,niso)
  integer::i,p,j,si
  real*8::x,dx(3)
  x=r/msg%hbarc
  dx(:)=dr(:)/msg%hbarc
  spinwf(:,:)=dcmplx(0.d0,0.d0)
  ! calc S_+^i r^i
  do p=1,2
     spinwf(:,:) = spinwf(:,:)+s_wf(:,:,p,index)*dx(index)/x
  enddo
  ccwf=conjg(wf_in)
  exp=0.d0
  !multiply by z
  spinwf(:,:)=spinwf(:,:)*dx(3)
  !inner product
  do i = 1,nspin
     do j=1,niso
        
        exp=exp+spinwf(i,j)*ccwf(i,j)
     enddo
  enddo
  return
end subroutine rSrop

