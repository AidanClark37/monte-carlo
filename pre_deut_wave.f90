module pre_deut
  implicit none




  real*8::warray(4,2,2)
  
contains

  subroutine pre_deut_wave()
  use wigner
  implicit none

  integer::i,j,k,m1,m2,m10,m20,t1,t2
  
  call start_wigner()
!  write(*,*) 'starting pre deut wave'
  do j =1,2
     do i=1,4
        do k=1,2
           m1=1-2*MOD(CEILING(MOD(DBLE(i),2.d0)),2)
           m2=1-2*MOD(CEILING(MOD(i/2.d0,2.d0)),2)
           !m1=sign(1,i)
           !m2=sign(2,i)
           !t1=sign(1,iso_array(j,2,(2,3)))
           !t2=sign(1,iso_array(j,2,(2,3)))
           t1=2*j-3
           t2=3-2*j
           !m1=1-2*FLOOR(0.5*(i-0.1))
           !m2=2*MOD(i,2)-1
           !write(*,*) "m1: ",m1,"m10 ",m10,"m2: ",m2,"m20: ",m20 
           warray(i,j,k) =  cgor(4*k-4,2-m1-m2,2,m1+m2,2,2)*cgor(1,m1,1,m2,2,m1+m2)*cgor(1,t1,1,t2,0,0)
           !write(*,*)warray(i,j,k)
           !write(*,'(A,1i3,1i3,1i3,A,f8.4)')'SO(',i,j,k,') = ',cgor(4*k-4,2-m1-m2,2,m1+m2,2,2)
           !write(*,'(A,1i3,1i3,1i3,A,f8.4)')'SS(',i,j,k,') = ',cgor(1,m1,1,m2,2,m1+m2)
           !write(*,'(A,1i3,1i3,1i3,A,f8.4)')'TT(',i,j,k,') = ',cgor(1,3-2*j,1,2*j-3,0,0)
        enddo
     enddo
  enddo

  

  return

endsubroutine pre_deut_wave
end module pre_deut
