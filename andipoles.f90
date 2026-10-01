
!================================
!diagram a
!================================

!--------------------------------
!a1
!--------------------------------

subroutine  an_NLO_diagram_a1_g0(exp)
    use mpi_modules
  use param_calc
    use struct_funcs
    use const_param
    
  implicit none
  real*8,intent(out)::exp
  
  
    exp = 0.d0
  return
end subroutine an_NLO_diagram_a1_g0

!-----------------------------------------
!a2
!-----------------------------------------

subroutine an_NLO_diagram_a2_g0(exp)
  
  use mpi_modules
  use param_calc
  use const_param
  
  use struct_funcs
  
  implicit none
  
  real*8,intent(out)::exp
  real*8::cst,x,t12wf,spwf
  integer::i,j

  exp=0.d0
  
  call sz_op(.true.,flarr,0,spwf)
  t12wf=-3.d0*spwf

  cst = -msg%hbarc*gA*msg%lambda/(4*fpi*nmass)
  exp=cst*t12wf
  
  return
end subroutine an_NLO_diagram_a2_g0


subroutine an_NLO_diagram_a2_g1(exp)
 
  use mpi_modules
  use param_calc
  use const_param
  
  use struct_funcs
 
  implicit none
  
  real*8,intent(out)::exp
  real*8::cst,x,wfout
  
 
  
  integer::i,j

  exp=0.d0

  
  
  
  
  cst = -gA*msg%hbarc*msg%lambda/(4*fpi*nmass)
  !write(*,*)cst
  !stop
  call sz_op(.true.,flarr,0,wfout)
   exp = cst*wfout
  return
end subroutine an_NLO_diagram_a2_g1

subroutine an_NLO_diagram_a2_g2(exp)

  use mpi_modules
  use param_calc
  use const_param
  
  use struct_funcs
  
  implicit none

  real*8,intent(out)::exp
  real*8::cst,x,spwf,t12wf,wfout

 

  integer::i,j
 
  exp=0.d0
 
  call sz_op(.true.,flarr,0,spwf)
  t12wf=-spwf
  
  cst = -msg%hbarc*gA*msg%lambda/(4*fpi*nmass  )
  wfout=cst*t12wf
  exp=wfout
  return
end subroutine an_NLO_diagram_a2_g2


!============================================
!diagram b
!============================================

!---------------------------------------------
!daigram b1
!---------------------------------------------

subroutine an_NLO_diagram_b1_g1V(exp)
  
  use mpi_modules
  use param_calc
  use const_param
  
  use struct_funcs
  
  implicit none
  
  real*8,intent(out)::exp
  real*8::cst,x,spinwf,t12wf,wfout
  
  
  
  integer::i,j
  
  exp=0.d0
  
  
  
  call zsr_op(.true.,flarr,1,spinwf)
  
  t12wf=-2.d0*spinwf

  
  cst = -msg%hbarc*gA*msg%lambda/(4*fpi**2)
  wfout=cst*t12wf
  
  exp=wfout
  return
end subroutine an_NLO_diagram_b1_g1V

!---------------------------------------------
!diagram b2
!---------------------------------------------
subroutine an_NLO_diagram_b2_g0(exp)
  
  use mpi_modules
  use param_calc
  use const_param
  
  use struct_funcs
  
  implicit none
  
  real*8,intent(out)::exp
  real*8::cst,spinwf,wfout
  
  integer::i,j,p
 
   exp=0.d0
 
   call sz_op(.true.,glarr,0,spinwf)
       
   
  cst = gA*msg%hbarc*dmass2/(4*fpi**2*msg%lambda)
  
  
  
  wfout=-2.d0*cst*spinwf
  exp=wfout
  return
end subroutine an_NLO_diagram_b2_g0

!---------------------------------------------
!diagram b3
!---------------------------------------------
subroutine an_NLO_diagram_b3_g0(exp)
  
  use mpi_modules
  use param_calc
  use const_param
  
  use struct_funcs
  
  implicit none
  
  real*8,intent(out)::exp
  real*8::cst(2),x,sr(3)
  
  
  
  integer::i,j,p
  exp = 0.d0
  return
end subroutine an_NLO_diagram_b3_g0
!=============================================
!diagram c
!=============================================

subroutine an_NLO_diagram_c_g1V(exp)
 
  use mpi_modules
  use param_calc
  use const_param
 
  use struct_funcs
 
  implicit none
  
  real*8,intent(out)::exp
  real*8::cst,x,spinwf,t12wf,wfout
  
  
  
  integer::i,j
  
  exp=0.d0
  
  
  
  call zsr_op(.true.,flarr,1,spinwf)
  
  t12wf=-2.d0*spinwf

  
  cst = -gA*msg%hbarc*msg%lambda/(4*fpi**2)
  exp=cst*(t12wf)
return  
  
  
end subroutine an_NLO_diagram_c_g1V

!=============================================
!diagram d
!=============================================

subroutine an_LO_diagram_d_d0(exp)
 
  use mpi_modules
  use param_calc
  use struct_funcs
  implicit none
  
  real*8::spinwf
  real*8,intent(out)::exp
  
  integer::i,p,j
  
  
  exp=0.d0
  !calc S_+^i 
  call sz_op(.false.,flarr,0,spinwf)
  
  spinwf=-0.5d0*spinwf
  exp = spinwf

return
end subroutine an_LO_diagram_d_d0



subroutine an_LO_diagram_d_d1(exp)
  
  use mpi_modules
  use param_calc
  use struct_funcs
  implicit none
  
  
  real*8,intent(out)::exp
  
  
 exp=0.d0

return
end subroutine an_LO_diagram_d_d1

subroutine an_NLO_diagram_e_Delta(exp)
   use spin_ops_act_store
  use mpi_modules
  use param_calc
  use const_param
  use iso_ops_act
  use struct_funcs
  use interpolation
    implicit none
  real*8,intent(out)::exp
  real*8:: cst

  exp=0.d0  
  

  return
   end subroutine an_NLO_diagram_e_Delta

   
subroutine an_NLO_diagram_f_Delta(exp)
  
  use mpi_modules
  use param_calc
  use const_param
  
  use struct_funcs
  
  implicit none
 
  real*8,intent(out)::exp
  real*8::cst,x,lcst,hcst,fcst
 
 

 
  exp=0.d0
  cst= 3*gA**3*nmass/(32*msg%lambda*pi**2*fpi**3)
  call zsr_op(.true.,Hlarr,1,hcst)
  call zsr_op(.true.,Llarr,1,lcst)
  call zsr_op(.true.,flarr,1,fcst)
  exp= cst*(pimass**2*(3*lcst-4*hcst)-msg%lambda**2*fcst)
end subroutine an_NLO_diagram_f_Delta

   subroutine an_N3LO_diagram_CT_CS(rr,dr,r,wf_in,index,exp)
     use spin_ops_act_store
     use mpi_modules
     use param_calc
     use const_param
     use iso_ops_act
     use struct_funcs
     use interpolation
     use diff
     implicit none
     real*8,intent(in)::rr(3,2),dr(3),r
     real*8,intent(out)::exp
     real*8::cst(2),x,sr(3)
     complex*16,intent(in)::wf_in(nspin,niso)
     complex*16::ccwf(nspin,niso),spin_plus(nspin,niso,3),spin_minus(nspin,niso,3),spin_x(nspin,niso,3),smdotr(nspin,niso),spdotdr(nspin,niso),sxdotr(nspin,niso),smdotdr(nspin,niso),wfout(nspin,niso),tempdrwf(2,nspin,niso)
     integer,intent(in)::index
     integer::i,j
     x=r/msg%hbarc
     sr = (rr(:,1)+rr(:,2))/2
     spin_plus = s_wf(:,:,1,:)+s_wf(:,:,2,:)
  !   do i = 1,4
  !      write(*,*)spin_plus(i,1,3),spin_plus(i,2,3)
  !   enddo
!     stop
     spin_minus= s_wf(:,:,1,:)-s_wf(:,:,2,:)
     spin_x(:,:,1) =ss_wf(:,:,1,2,2,3)-ss_wf(:,:,1,2,3,2)
     spin_x(:,:,2) =ss_wf(:,:,1,2,3,1)-ss_wf(:,:,1,2,1,3)
     spin_x(:,:,3) =ss_wf(:,:,1,2,1,2)-ss_wf(:,:,1,2,2,1)
     smdotr=dcmplx(0.d0,0.d0)
     spdotdr=dcmplx(0.d0,0.d0)
     smdotdr=dcmplx(0.d0,0.d0)
     sxdotr=dcmplx(0.d0,0.d0)
     do i = 1,3
        smdotr = smdotr + spin_minus(:,:,i)*dr(i)
        call derivative(spin_minus(:,:,i),rr,1,i,tempdrwf(1,:,:))
        call derivative(spin_minus(:,:,i),rr,2,i,tempdrwf(2,:,:))
        smdotdr = smdotdr + tempdrwf(1,:,:)-tempdrwf(2,:,:)
        call derivative(spin_plus(:,:,i),rr,1,i,tempdrwf(1,:,:))
        call derivative(spin_plus(:,:,i),rr,1,i,tempdrwf(2,:,:))
        spdotdr=spdotdr+tempdrwf(1,:,:)-tempdrwf(2,:,:)
        sxdotr=sxdotr+spin_x(:,:,i)*dr(i)
     enddo
     smdotr=smdotr/r
     sxdotr=sxdotr/r
     call interpolate(Clarr,0,x,cst(1))
     call interpolate(Clarr,1,x,cst(2))
     wfout=dcmplx(0.d0,0.d0)
     cst=-1*dcmplx(0.d0,1.d0)*cst/(64*nmass*(fpi*msg%lambda)**3)
     wfout=cst(1)*(0.5d0*spin_plus(:,:,index)+sr(index)*(spdotdr-2*smdotdr))
     wfout=wfout-cst(2)*sr(index)*(dcmplx(0.d0,1.d0)*smdotr+4*sxdotr)
 !    write(*,*)spin_plus(1,1,index)
     exp = 0.d0
     ccwf=dconjg(wf_in)
     do i = 1,nspin
!        write(*,*)spin_plus(i,1,3),spin_plus(i,2,3)
        do j = 1,niso
           exp=wfout(i,j)*wfout(i,j)
        end do
     end do
     !stop
     exp = 0.d0
   end subroutine an_N3LO_diagram_CT_CS
