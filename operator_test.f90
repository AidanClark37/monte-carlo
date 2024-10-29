subroutine operator_calc()
  use pre_deut
  integer::a,p,state_sign,iarray(2)

  real*8::rr(3,2),ysol(2,40),obs
  complex*16::wfa(4,2),tau_wfa(4,2)
  iarray(1)=2
  iarray(2)=3
  call pre_deut_wave()
  do i = 1,80
!     write(*,*) (i-1)/40+1,mod(i,41)+i/41
     read(5,*)ysol((i-1)/40+1,mod(i,41)+i/41)
  enddo
 ! do i = 1,40
 !    write(*,*)ysol(1,i)
 ! enddo
  
     
  rr(1,1)=0.3
  rr(1,2)=-0.3
  rr(2,1)=-0.2
  rr(2,2)=0.2
  rr(3,1)=1.04
  rr(3,2)=0.5
  call deut_wave(rr,wfa,ysol)
 ! write(*,*)wfa(1,1)
 a=1
 p=1
 write(*,*)'here'
!  do p = 1,2
 !    do a = 1,3
        !call tau_exp_val(iarray,cwf,b,N,niso,nspin,tau_exp)
        !write(*,*) obs

        write(*,*)'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        write(*,*)"wavefunction:"
        write(*,*) ' '
        do i = 1,4
           write(*,*) wfa(i,1),wfa(i,2)
           !if (a==1) then

              !write(*,*) 'particle',p,'index',i
              !write(*,*)'here'
              !write(*,*)'spin:',state_sign(p,i)
            !  endif
        enddo
        call isospin(iarray,wfa,p,a,4,2,tau_wfa)
        write(*,*)'--------------------------------------'
        write(*,*) ' '
        write(*,*)"spin acted wavefunction"
        !write(*,*)"p=",p,"a=",a
        
        do i = 1,4
           write(*,*) tau_wfa(i,1),tau_wfa(i,2)
           !write(*,*)'here'
        enddo
!        write(*,*)'no here'
!     enddo
!     write(*,*)'actually here'
!  enddo
  

end subroutine operator_calc

