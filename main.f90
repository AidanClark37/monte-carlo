program monte_carlo
  use mpi_modules
  use pre_deut 
  implicit none

  integer,parameter::npart = 2 
  real*8,allocatable::rpart_o(:,:)
  real*8::norm
  logical::acc
  complex*16::cwf(4,2)
  integer::i,j,ii,jj,iq
  real*8::rnd
  integer::jz
  integer::nq
  integer::nspin,niso
  integer::accp,acc_move
  real*8::r,rr(3),rcm(3)
  real*8,allocatable::obs(:)
  real*8,allocatable::obs_av(:)
  real*8,allocatable::obs_av_w(:)
  real*8,allocatable::obs_sg_w(:)
  real*8,allocatable::mean_obs(:)
  real*8,allocatable::sigma_obs(:)

  real*8::wfa(2,40)
  real*8::obs0
  nspin=4
  niso=2
  !call operator_calc()
  !stop
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

  write(*,*)proc_rank!quantity you want to check
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
  do i=1,nwalks_for_proc

     !Thermalization-----------------------------------------
     rnd=proc_rank*1047.d0+i*356.d0
     rpart_o=0.d0
     norm=0.d0

     do j=1,neq
        call step(rnd,rpart_o,npart,jz,cwf,norm,acc)
        !write(*,*)j,acc,norm
     end do
     !write(*,*)'here'
     !-----------------------------------------------------

     !Sampling phase------------------------------------------
     obs_av=0.d0
     do ii=1,nav
     do jj=1,ncorr
        call step(rnd,rpart_o,npart,jz,cwf,norm,acc)
       ! write(*,*)ii,jj,acc,norm
        if(acc)acc_move=acc_move+1
     end do
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
     call spin_exp_val(cwf,2,2,2,4,obs(1))
     !write(*,*)'obs:',obs_av(1)
     obs_av=obs_av+obs
!     write(*,*)obs_av(1),rho(1)
     end do
     obs_av_w=obs_av_w+obs_av/nav
     obs_sg_w=obs_sg_w+(obs_av/nav)**2
     ! write(*,*)1.d0*acc_move(i)/nav/ncorr/i
  end do

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

