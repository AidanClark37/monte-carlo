program monte_carlo
  use mpi_modules
  implicit none

  integer,parameter::npart = 2 
  real*8,allocatable::rpart_o(:,:)
  real*8::norm
  logical::acc
  complex*16,allocatable::cwf(:,:)
x  integer::i,j,ii,jj,iq
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
  real*8,allocatable::warray(:,:,:)


  call start_mpi()

  if(proc_rank.eq.0)then
     read(5,*)msg%wf_file !<-- This is used to read the file of the wave function change as you need
     read(5,*)msg%nwalk
     read(5,*)msg%neq
     read(5,*)msg%nav
     read(5,*)msg%ncorr
     read(5,*)msg%sigma
  end if

  call mpi_broadcast_input()
  call general_setting()
  allocate(warray(4,2,2))
  !Initialize the wave function here
  call pre_deut_wave(warray)

  write(*,*)proc_rank!quantity you want to check
  allocate(rpart_o(3,npart))
  jz= 1;!select the jz (typically jz=tot j


  !initialization of the observables vectors
  allocate(obs(nq)) !<-selcect nq based on what you need
  allocate(obs_av(nq))
  allocate(obs_av_w(nq))
  allocate(obs_sg_w(nq))
  allocate(mean_obs(nq))
  allocate(sigma_obs(nq))
  allocate(cwf(nspin,niso))!initialize matrix wave function
  
    
  obs_av_w=0.d0
  obs_sg_w=0.d0
  do i=1,nwalks_for_proc

     !Thermalization-----------------------------------------
     rnd=proc_rank*1047.d0+i*356.d0
     rpart_o=0.d0
     norm=0.d0
     do j=1,neq
        call step(rnd,sigma,rpart_o,npart,jz,cwf,norm,acc)
        !write(*,*)j,acc,norm
     end do
     !-----------------------------------------------------

     !Sampling phase------------------------------------------
     obs_av=0.d0
     do ii=1,nav
     do jj=1,ncorr
        call step(rnd,sigma,rpart_o,npart,jz,cwf,norm,acc)
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
     
     obs_av=obs_av+obs
!     write(*,*)obs_av(1),rho(1)
     end do
     obs_av_w=obs_av_w+obs_av/nav
     obs_sg_w=obs_sg_w+(obs_av/nav)**2
    ! write(*,*)1.d0*acc_move(i)/nav/ncorr/i     
  end do

  
  
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

