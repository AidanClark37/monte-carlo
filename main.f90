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
  real*8::norm,f_prime_lambda,xxxxxx,f_lambda,f_2prime_lambda,g_lambda,g_prime_lambda,g_2prime_lambda,exp
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
!call cpu_time(start)
  call pre_deut_wave()
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
 
     
call parameter_calc(msg%npart,msg%Tz) !calculate constants - niso,nspin,
call isospin_ops() ! call arrays for spin and isopin operations
  call struct_func_calc(0.d0,10.d0/msg%hbarc) !calculate interpolation poins for structure funcitons
end if
  call mpi_broadcast_input()
  call general_setting()
 
  
  !Initialize the wave function here
  allocate(rpart_o(3,msg%npart))
  jz= 1;!select the jz (typically jz=tot j
  nq=10

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


        
 
 call all_operators(cwf)!precalculate spin operated wavefuncitons
     

 
 if (proc_rank.eq.0)then


        if (ii.eq.1)then
           if (i.eq.1)then
              !-------------------------
              !calculating spin operator example expectation vaules
              !-------------------------
              
              allocate(testcwf(4,2))
              testcwf2(:,:)=conjg(cwf)
     write(*,*)"------------------------"
     write(*,*)"base wavefunction"
     write(*,*)"--------------------------"
     do si = 1,4
        write(*,*)cwf(si,1),cwf(si,2)
     enddo
     write(*,*)"-------------------"
     !sigma 1 
     do sp =1,3
        exp=dcmplx(0.d0,0.d0)
        testcwf(:,:)=(s_wf(:,:,1,sp))
       write(*,*) "S1",sp
       write(*,*)"-------------------"
       do si  =1,4
          write(*,*) testcwf(si,1),testcwf(si,2)
          ! do sj=1,2
         !    exp = exp+testcwf(si,sj)*testcwf2(si,sj)
         !    enddo
       enddo
       call d0(cwf,sp,exp)
       write(*,*)"--------------------"
       write(*,*)"exp:",exp
       write(*,*)"------------------"
    enddo

    !sigma 2
    do sp =1,3
	exp=dcmplx(0.d0,0.d0)
        testcwf(:,:)=(s_wf(:,:,2,sp))
       write(*,*) "S2",sp
       write(*,*)"-------------------"
       do si  =1,4
          write(*,*) testcwf(si,1),testcwf(si,2)
          do sj=1,2
             exp = exp+testcwf(si,sj)*testcwf2(si,sj)
             enddo
       enddo
       write(*,*)"--------------------"
       write(*,*)"exp:",exp
       write(*,*)"------------------"
    enddo

    
    !S_+^i r^i
     do sp =1,3
	exp=0.d0
       testcwf(:,:)=(s_wf(:,:,1,sp)+s_wf(:,:,2,sp))*dr(sp)/r
       write(*,*) "S_+r",sp
       write(*,*)"-------------------"
       do si  =1,4
          write(*,*) testcwf(si,1),testcwf(si,2)
          !do sj=1,2
          !   exp = exp+testcwf(si,sj)*testcwf2(si,sj)
          !   enddo
       enddo
       call Srop(r,dr,sp,cwf,exp)
       write(*,*)"--------------------"
       write(*,*)"exp:",exp
       write(*,*)"------------------"
    enddo
    !z S_+^i r^i
    do sp =1,3
        exp=dcmplx(0.d0,0.d0)
       testcwf(:,:)=(dr(3)/msg%hbarc)*(s_wf(:,:,1,sp)+s_wf(:,:,2,sp))*dr(sp)/r
       write(*,*) "rS_+r",sp
       write(*,*)"-------------------"
       do si  =1,4
          write(*,*) testcwf(si,1),testcwf(si,2)
          !do sj=1,2
          !   exp = exp+testcwf(si,sj)*testcwf2(si,sj)
          !   enddo
       enddo
       call rSrop(r,dr,sp,cwf,exp)
       write(*,*)"--------------------"
       write(*,*)"exp:",exp
       write(*,*)"------------------"
    enddo
    !calculating analytic for d0 (S_+^z)
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
             ! write(*,*)"0",func0
           endif
        endif
     endif

     !call g1V operator (~  z S_+ . rhat)
     call g1V(r,dr,cwf,3,obs(1))
     !call S_+^i for i=x,y,z
     call d0(cwf,1,obs(2))
     call d0(cwf,2,obs(3))
     call d0(cwf,3,obs(4))
     !call S_+^i r^i for i=x,y,z
     call Srop(r,dr,1,cwf,obs(5))
     call Srop(r,dr,2,cwf,obs(6))
     call Srop(r,dr,3,cwf,obs(7))
     !call z S_+^i r^i for i=x,y,z
     call rSrop(r,dr,1,cwf,obs(8))
     call rSrop(r,dr,2,cwf,obs(9))
     call rSrop(r,dr,3,cwf,obs(10))
     obs=obs/norm
     obs_av=obs_av+obs
  end do
  !     write(*,*)proc_rank,obs_av,obs_av/nav
  !a
     obs_av_w=obs_av_w+obs_av/nav
     obs_sg_w=obs_sg_w+(obs_av/nav)**2
     ! write(*,*)1.d0*acc_move(i)/nav/ncorr/i

  end do
  !close(20)
 ! close(15)
 ! close(16)
 ! close(17)
!  write(*,*)proc_rank,obs_av_w
!  call spin_exp_val(cwf,3,2,2,4,obs0)
!  write(*,*) obs0
  
  call op_reduction(obs_av_w,obs_sg_w,nq,mean_obs,sigma_obs)! Average on the walkers
  call int_reduction(acc_move,accp)

  if(proc_rank.eq.0)then
     allocate(names(nq))
     names(1) = "g1V"
     names(2) = "s+x"
     names(3) = "s+y"
     names(4) = "s+z"
     names(5) = "sx rx"
     names(6) = "sy ry"
     names(7) = "sz rz"
     names(8) = "r sx rx"
     names(9) = "r sy ry"
     names(10)= "r sz rz"
     write(*,*)"Acceptance = ",dfloat(accp)/dfloat(nwalk*nav*ncorr)
     write(*,*)"Final results"
     write(*,*)"----------------------------"
     open(unit=21,position="append",file="g1v.txt")
!     write(21,*)mean_obs(1)
     close(21)
     do iq=1,nq
        write(*,*)names(iq),mean_obs(iq),sigma_obs(iq)
        !write(*,*)fl1arr(1)
     end do
     write(*,*)"expected for g1V",mean_obs(8)+mean_obs(9)+mean_obs(10)
     write(*,*) "expected result for 2:",-2*ftotal
     write(*,*) "expected result for 1:",func0
     !write(*,*)-func0,-func2,sqrt(2.d0)
  end if

  !deallocate(iarray)
  !deallocate(inviarray)
  !deallocate(p0arr)
  !deallocate(p1arr)
  !deallocate(p2arr)
  !deallocate(pxarr)
  !deallocate(c0arr)
  !deallocate(c1arr)
  !deallocate(c2arr)
  !deallocate(cxarr)
  !deallocate(flarr)
  !deallocate(Clarr)
  !deallocate(glarr)
  !deallocate(Llarr)
  !deallocate(Hlarr)
  !deallocate(G1larr)
  !deallocate(fl1arr)
  !deallocate(Cl1arr)
  !deallocate(gl1arr)
  !deallocate(Ll1arr)
  !deallocate(Hl1arr)
  !deallocate(G1l1arr)
  !deallocate(gl2arr)
  !deallocate(G1l2arr)
call end_mpi()
end program monte_carlo

