program monte_carlo
  use mpi_modules
  use pre_deut
  use operator_calc
  implicit none

  integer,parameter::npart = 2
  real*8,allocatable::rpart_o(:,:)
  real*8::norm,f_lambda
  logical::acc
  complex*16::cwf(4,2)
  integer::i,j,ii,jj,iq,ir
  real*8::rnd
  integer::jz,n
  integer::nq,si,sj
  integer::nspin,niso
  integer::accp,acc_move
  real*8::r,rr(3,2),rcm(3),dr(3),q(3)
  real*8,allocatable::obs(:)
  real*8,allocatable::obs_av(:)
  real*8,allocatable::obs_av_w(:)
  real*8,allocatable::obs_sg_w(:)
  real*8,allocatable::mean_obs(:)
  real*8,allocatable::sigma_obs(:)
  
  real*8::wfa(2,40)
  real*8::obs0,rr2
  !write(*,*)'test'
  
  real*8::b_int,b5,h,bounds(2),f_0
  integer::ias,ndim
  real*8,allocatable::a(:)
  
  





 ! call operator_calc()
   call start_mpi()
  if(proc_rank.eq.0)then
     do i = 1,80
       ! write(*,*)(i-1)/100+1,mod(i,100)
        read(5,*)wfa((i-1)/40+1,mod(i,41)+i/41) !<-- This is used to read the file of the wave function change as you need
     enddo
     msg%wf=wfa
     read(5,*)msg%nwalk
     read(5,*)msg%neq
     read(5,*)msg%nav
     read(5,*)msg%ncorr
     read(5,*)msg%sigma
     read(5,*)msg%nla
     read(5,*)msg%iarray(1)
     read(5,*)msg%iarray(2)
     read(5,*)msg%mass
     read(5,*)msg%lambda
  end if
!  do i=1,2
 !    do j=1,40
  !      write(*,*)'wfa(',i,j,')=', msg%wf(i,j)
   !  enddo
   ! enddo
  !enddo
  !write(*,*)'neq=', msg%neq
  !write(*,*)'nav=', msg%nav
  !write(*,*)'ncorr=',msg%ncorr
  !write(*,*)'sigma=',msg%sigma
  call mpi_broadcast_input()
  call general_setting()
 
  !Initialize the wave function here
  call pre_deut_wave()

  ndim=1000
  ias=1
  bounds(1)=0.d0
  bounds(2)=2.d0*msg%lambda
  h=(bounds(2)-bounds(1))/real(ndim)
  allocate(a(ndim))
  b_int = 0.d0
  do i = 1,ndim



     a(i)=((real(i)*h)**2)*f_0(real(i)*h,2.d0,msg%lambda,msg%mass)


  enddo
b_int = b5(1,ndim,0.d0,0.d0,h,0.d0,0.d0,ndim,a,ias)
  
  write(*,*)"gaulag f(2)",f_lambda(2.d0,msg%mass,msg%lambda,msg%nla)
  write(*,*)"b5 f(2)",b_int
  stop

!  write(*,*)proc_rank!quantity you want to check
  allocate(rpart_o(3,npart))
  jz= 1;!select the jz (typically jz=tot j
  nq=1

  !initialization of the observables vectors
  allocate(obs(nq)) !<-selcect nq based on what you need
  allocate(obs_av(nq))
  allocate(obs_av_w(nq))
  allocate(obs_sg_w(nq))
  allocate(mean_obs(nq))
  allocate(sigma_obs(nq))
  !allocate(cwf(nspin,niso))!initialize matrix wave function
  
    
  obs_av_w=0.d0
  obs_sg_w=0.d0
  acc_move=0
  do i=1,nwalks_for_proc

     !Thermalization-----------------------------------------
       rnd=proc_rank*1047.d0+i*353.d0
     !rnd=1047.d0+i*353.d0
     rpart_o=0.d0
     norm=0.d0

     do j=1,neq
        call step(rnd,rpart_o,npart,jz,cwf,norm,rr,dr,r,acc)
        ! write(*,*)proc_rank,j,acc,norm
        
     end do
     !write(*,*)'here'
     !-----------------------------------------------------

     !Sampling phase------------------------------------------
     obs_av=0.d0
     do ii=1,nav
     do jj=1,ncorr
        call step(rnd,rpart_o,npart,jz,cwf,norm,rr,dr,r,acc)
       ! write(*,*)ii,jj,acc,norm
        if(acc)acc_move=acc_move+1
     end do
!     write(*,*)proc_rank,acc_move
     rcm=0.d0
     !remove center of mass
     do jj=1,npart
        rcm(1)=rcm(1)+rpart_o(1,jj)
        rcm(2)=rcm(2)+rpart_o(2,jj)
        rcm(3)=rcm(3)+rpart_o(3,jj)
     end do
     rcm=rcm/dfloat(npart)
     !Here remove the center of mass of the particle 
!!!!! Here you have to call your subroutine that compute the observable (obs as output)
     !     stop
     !call spin_exp_val(cwf,3,2,2,4,obs(1))
     !call tau_exp_val(iarray,cwf,3,2,2,4,obs(1))
     !stop
     
     
     call all_operators(iarray,cwf,2,2)
     !write(*,*) st_wf(1,1,1,3,1)

     !call rho_NNg_other(cwf,2,2,1,obs(1))
     call rho_NNpTRV_NNpgPC_1(rr,dr,r,q,msg%iarray,cwf,2,2,obs(1))
     
     
     
     
     
     
           !     do si = 1, 4
!        write(*,*)"wf:",cwf(si,:)
!     enddo
!     do si = 1,4
!        write(*,*)"twf:",tt_wf(si,:,1,2)
!        enddo
!     stop
     !call rho_NNg(iarray,cwf,2,2,1,obs(1))
     !call spin_exp_val(cwf,3,2,2,4,obs(2))
     ! write(*,*)'here1'
     !call tau_exp_val(iarray,cwf,3,2,2,4,obs(1))
     
     !write(*,*)'here2'
     !write(*,*)obs(1)
     !,obs(2)!  if(proc_rank.eq.1)write(*,*)'obs:',obs(1),norm
     !obs(1)=obs(1)/norm
     obs=obs/norm
     !write(*,*)obs(1)
     
     
     obs_av=obs_av+obs
    
!     write(*,*)obs_av(1),rho(1)
     end do
!     write(*,*)proc_rank,obs_av,obs_av/nav
     obs_av_w=obs_av_w+obs_av/nav
     obs_sg_w=obs_sg_w+(obs_av/nav)**2
     ! write(*,*)1.d0*acc_move(i)/nav/ncorr/i
     !stop
  end do

!  write(*,*)proc_rank,obs_av_w
!  call spin_exp_val(cwf,3,2,2,4,obs0)
!  write(*,*) obs0
  
  call op_reduction(obs_av_w,obs_sg_w,nq,mean_obs,sigma_obs)! Average on the walkers
  call int_reduction(acc_move,accp)

  if(proc_rank.eq.0)then
     
     write(*,*)"Acceptance = ",dfloat(accp)/dfloat(nwalk*nav*ncorr)
     write(*,*)"Final results"
     write(*,*)"----------------------------"
     do iq=1,nq
        write(*,*)iq,mean_obs(iq),sigma_obs(iq)
     end do
  end if

  call end_mpi()  
end program monte_carlo

