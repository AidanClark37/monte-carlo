program monte_carlo
  use mpi_modules
  use pre_deut
  use operator_calc
  use wigner
  use isospin_operators
 use tau_operator
 use param_calc
 use struct_funcs
 use interpolation
  implicit none


  real*8,allocatable::rpart_o(:,:)
  real*8::norm,f_prime_lambda,xxxxxx,f_lambda,f_2prime_lambda,g_lambda,g_prime_lambda,g_2prime_lambda
  logical::acc
  complex*16::testcwf2(4,2),testcwf3(4,2),cout
 complex*16,allocatable::cwf(:,:),testcwf(:,:),twf(:,:)
  integer::is,i,j,ii,jj,iq,ir,p
  real*8::rnd,l_lambda,l_prime_lambda,l_2prime_lambda,c_lambda,c_prime_lambda,c_2prime_lambda
  real*8::h_lambda,h_prime_lambda,h_2prime_lambda,g1_lambda,g1_prime_lambda,g1_2prime_lambda
  integer::jz,n,iout
  integer::nq,si,sj,sp,sb,fi
  integer::accp,acc_move
  real*8::r,rr(3,2),rcm(3),dr(3),q(3),func0,func2
  real*8::bessel,f0sum,f2sum,ftotal
  real*8,allocatable::obs(:),out(:)
  real*8,allocatable::obs_av(:)
  real*8,allocatable::obs_av_w(:)
  real*8,allocatable::obs_sg_w(:)
  real*8,allocatable::mean_obs(:)
  real*8,allocatable::sigma_obs(:)
  real*8::wfa(2,40),xi
  real*8::start,finish
  real*8::obs0,rr2,x,ytest
  character(len=5),allocatable::names(:)
  real*8::b_int,b5,h,bounds(2),f_0
  integer::ias,ndim

  real*8,allocatable::a(:)
call cpu_time(start)
  call pre_deut_wave()
!  do si = 1,4
!     write(*,*) warray(si,1,2),warray(si,2,2)
!  enddo
  !stop
  call start_mpi()

  if(proc_rank.eq.0)then
     open(unit=100,file="input.txt")
     do i = 1,80
        read(100,*)wfa((i-1)/40+1,mod(i,41)+i/41) !<-- This is used to read the file of the wave function change as you need
     enddo
    
     msg%wf=wfa
     read(100,*)msg%nwalk
     read(100,*)msg%neq
     read(100,*)msg%nav
     read(100,*)msg%ncorr
     read(100,*)msg%sigma
     read(100,*)msg%nla
     
     
     read(100,*)msg%npart
     read(100,*)msg%Tz
     read(100,*)msg%mass
     read(100,*)msg%lambda
     read(100,*)msg%hbarc
     read(100,*)msg%npoint
    ! write(*,*)"npoint",msg%hbarc

!     open(unit=79,file="sqif595.txt")
!     write(*,*)"eh?"
!     call struct_func_calc(0.d0,10.d0/msg%hbarc)
  !   write(*,*)"eh"
 !    write(*,*)xarray(2)
!     write(*,*)"eh."
 !    do si=0,99
!        call interpolate(flarr,(dble(si)+0.5)*10/(99*msg%hbarc),ytest)
        !write(*,*)xarray(si+1)*msg%hbarc
        !write(*,*)dble(si)*10/(99)
!        write(79,*)ytest
    !call interpolate(fl1arr,0.3d0/msg%hbarc,ytest)
!     enddo
! close(79)    
!     stop
!     xxxxxx=f_prime_lambda(0.3d0/msg%hbarc)
!     write(*,*)"real:",xxxxxx,"inter:",ytest
!     open(unit=99,file="sbf.txt")
!     open(unit=98,file="bf1.txt")
!     open(unit=97,file="bf2.txt")
!     open(unit=96,file="bC.txt")
!     open(unit=95,file="bC1.txt")
!     open(unit=94,file="bC2.txt")
!     open(unit=93,file="bg.txt")
!     open(unit=92,file="bg1.txt")
!     open(unit=91,file="bg2.txt")
!     open(unit=90,file="bL.txt")
!     open(unit=89,file="bL1.txt")
!     open(unit=88,file="bL2.txt")
!     open(unit=87,file="bH.txt")
!     open(unit=86,file="bH1.txt")
!     open(unit=85,file="bH2.txt")
!     open(unit=84,file="bG.txt")
!     open(unit=83,file="bG1.txt")
!     open(unit=82,file="bG2.txt")
!  do si = 0,99
!     xi=(dble(si)+0.5)*10/(99*msg%hbarc)
!         write(*,*)xi*msg%hbarc
!     xxxxxx=f_lambda(xi)
    ! write(99,*)xxxxxx
     !write(*,*)"flam:",si,xi*msg%hbarc
!     xxxxxx=f_prime_lambda(xi)
!     write(*,*)xxxxxx
!     write(98,*)xxxxxx
!     xxxxxx=f_2prime_lambda(xi)
!     write(97,*)xxxxxx
!     xxxxxx=C_lambda(xi)
!     write(96,*)xxxxxx
!     xxxxxx=C_prime_lambda(xi)
!     write(95,*)xxxxxx
!     xxxxxx=C_2prime_lambda(xi)
!     write(94,*)xxxxxx
!     xxxxxx=g_lambda(xi)
!     write(93,*)xxxxxx
!     xxxxxx=g_prime_lambda(xi)
!     write(92,*)xxxxxx
!     xxxxxx=g_2prime_lambda(xi)
!     write(91,*)xxxxxx
!     xxxxxx=L_lambda(xi)
!     write(90,*)xxxxxx
!     xxxxxx=L_prime_lambda(xi)
!     write(89,*)xxxxxx
!     xxxxxx=L_2prime_lambda(xi)
!     write(88,*)xxxxxx
!     xxxxxx=H_lambda(xi)
     !     write(*,*)xxxxxx
!     write(87,*)xxxxxx
!     xxxxxx=H_prime_lambda(xi)
!     write(86,*)xxxxxx
!     xxxxxx=H_2prime_lambda(xi)
!     write(85,*)xxxxxx
!     xxxxxx=G1_lambda(xi)
!     write(84,*)xxxxxx
!     xxxxxx=G1_prime_lambda(xi)
!     write(83,*)xxxxxx
!     xxxxxx=G1_2prime_lambda(xi)
!     write(82,*)xxxxxx
!  enddo
!  close(99)
     !  stop
     
call parameter_calc(msg%npart,msg%Tz)
call isospin_ops()
  call struct_func_calc(0.d0,10.d0/msg%hbarc)
end if
  call mpi_broadcast_input()
  call general_setting()
 
!  write(*,*)"niso in main",niso
!  write(*,*)"iarray:",iarray(1),iarray(2)
  
  !Initialize the wave function here
  allocate(rpart_o(3,msg%npart))
  jz= 1;!select the jz (typically jz=tot j
  nq=2

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
  acc_move=0
 
  
  
        open(unit=15,file="ervals.txt")
!        open(unit=16,file="qifvals.txt")
  !        open(unit=17,file="qbfvals.txt")
  !obs(3)=0.d0
  open(unit=20,file="est.txt")
  do i=1,nwalks_for_proc

     !Thermalization-----------------------------------------
       rnd=proc_rank*1047.d0+i*353.d0
     !rnd=1047.d0+i*353.d0
     rpart_o=0.d0
     norm=0.d0

 do j=1,neq
        call step(rnd,rpart_o,jz,cwf,norm,rr,dr,r,acc)
        ! write(*,*)proc_rank,j,acc,norm                                                                             

     end do
     
     !Sampling phase------------------------------------------
     obs_av=0.d0
     do ii=1,nav
        do jj=1,ncorr
           
           call step(rnd,rpart_o,jz,cwf,norm,rr,dr,r,acc)
           
       ! write(*,*)ii,jj,acc,norm
        if(acc)acc_move=acc_move+1
     end do
     !write(*,*)proc_rank,acc_move
     rcm=0.d0
     !remove center of mass
     do jj=1,msg%npart
        rcm(1)=rcm(1)+rpart_o(1,jj)
        rcm(2)=rcm(2)+rpart_o(2,jj)
        rcm(3)=rcm(3)+rpart_o(3,jj)
     end do
     rcm=rcm/dfloat(msg%npart)
     !Here remove the center of mass of the particle 
!!!!! Here you have to call your subroutine that compute the observable (obs as output)


        
 
 call all_operators(cwf)
     
!call struct_func_calc(0.d0,10.d0/msg%hbarc)     
 
 if (proc_rank.eq.0)then
    !write(*,*)"rmax::",xarray(msg%npoint)
!    call interpolate(flarr,r/msg%hbarc,ytest)
!    xxxxxx=f_lambda(r/msg%hbarc)
    write(15,*)r
!    write(16,*)ytest
!    write(17,*)xxxxxx
    !write(*,*)flarr(msg%npoint)

        if (ii.eq.1)then
           if (i.eq.1)then
              

 !    allocate(testcwf(4,2))                                                                                            
!     allocate(twf(4,2))
!     call tau2(cwf,1,2,testcwf)
!     ytest=0.d0
!     twf = conjg(cwf)
!     do si  =1,4
!        do sj=1,2
!           ytest=ytest+twf(si,sj)*testcwf(si,sj)
!        enddo
!     enddo
!     do si = 1,4
 !       write(*,*) cwf(si,1),cwf(si,2)
  !   enddo
     
   !     write(*,*)"-----------------------"
!do si=1,4
 !       write(*,*) s_wf(si,1,1,3),s_wf(si,2,1,3)
  !   enddo
     !stop
   !     write(*,*)"psi* t1 dot t2 psi",ytest
              f0sum=0.d0
              f2sum=0.d0
              do si=1,40
                 f0sum=f0sum + wfa(1,si)**2
                 f2sum =f2sum + wfa(2,si)**2
              enddo
              func0=0.d0
              func2=0.d0
              ftotal = 0.5*(f0sum-0.5*f2sum)
              call wave(func0)
              write(*,*)"0",func0
           endif
        endif
     endif

     call g1V(r,dr,cwf,3,obs(1))
     call d0(cwf,3,obs(2))
     !obs(3)=obs(3)+r*norm
     !write(*,*)"a!",r,obs(2)
     !stop
!     call g0(r,dr,cwf,3,obs(3))
     !write(*,*)obs(1)
     !,obs(2)!  if(proc_rank.eq.1)write(*,*)'obs:',obs(1),norm
     !obs(1)=obs(1)/norm
     obs=obs/norm
     !write(*,*)obs(1)
     
     write(20,*)obs(1)
     obs_av=obs_av+obs
    
!     write(*,*)obs_av(1),rho(1)



  end do
!     write(*,*)proc_rank,obs_av,obs_av/nav
     obs_av_w=obs_av_w+obs_av/nav
     obs_sg_w=obs_sg_w+(obs_av/nav)**2
     ! write(*,*)1.d0*acc_move(i)/nav/ncorr/i

  end do
  close(20)
  close(15)
 ! close(16)
 ! close(17)
!  write(*,*)proc_rank,obs_av_w
!  call spin_exp_val(cwf,3,2,2,4,obs0)
!  write(*,*) obs0
  
  call op_reduction(obs_av_w,obs_sg_w,nq,mean_obs,sigma_obs)! Average on the walkers
  call int_reduction(acc_move,accp)

  if(proc_rank.eq.0)then
     allocate(names(2))
     names(1) = "g1V"
     names(2) = "d0"
     write(*,*)"Acceptance = ",dfloat(accp)/dfloat(nwalk*nav*ncorr)
     write(*,*)"Final results"
     write(*,*)"----------------------------"
     open(unit=21,position="append",file="g1v.txt")
!     write(21,*)mean_obs(1)
     close(21)
     do iq=1,nq
        write(*,*)names(iq),mean_obs(iq),sigma_obs(iq)
        write(*,*)fl1arr(1)
     end do
     write(*,*) "expected result for 2:",-2*ftotal
     write(*,*) "expected result for 1:",(1.d0/3.d0)*func0
     !write(*,*)-func0,-func2,sqrt(2.d0)
  end if

  deallocate(iarray)
  deallocate(inviarray)
  deallocate(p0arr)
  deallocate(p1arr)
  deallocate(p2arr)
  deallocate(pxarr)
  deallocate(c0arr)
  deallocate(c1arr)
  deallocate(c2arr)
  deallocate(cxarr)
  deallocate(flarr)
  deallocate(Clarr)
  deallocate(glarr)
  deallocate(Llarr)
  deallocate(Hlarr)
  deallocate(fl1arr)
  deallocate(Cl1arr)
  deallocate(gl1arr)
  deallocate(Ll1arr)
  deallocate(Hl1arr)
  deallocate(G1l1arr)
  deallocate(gl2arr)
  deallocate(G1l2arr)
call end_mpi()
call cpu_time(finish)
write(*,*)"time",finish-start
end program monte_carlo

