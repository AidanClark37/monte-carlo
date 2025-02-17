subroutine deut_wave(rr,cwf,ysol,dr,r)
  use wigner
  use mpi_modules
  use pre_deut
  implicit none
  real*8::dr(3),ysol(2,40),rr(3,2)
  real*8,parameter::gamma=4.5
  integer::i,j,k,li,ll,M,M0,npart,ri
  real*8::asum(4,2,2),lsum(2),r,plaguer,thetax,phix,x,prod
  complex*16::ysum(4,2,2),ylm,cwf(4,2)
  !write(*,*)'starting  deut wave'
  !write(*,*)'here'
  !write(*,*)'wf1:',msg%wf(2,40)
!  do i = 1,2
!     do j = 1,40
!        write(*,*)'wf:',ysol(i,j)
        !ysol(i,j)=msg%wf(i,j)
!     enddo
 !    enddo
  !write(*,*)'here'
  !write(*,*)wf
  !gamma=deut_par%gamma 
  !warray=deut_wave%pre_wave !angular component of wavefunction
  !nla=deut_par%nla    !number of laguerre polynomials used
  r=0.d0
  do ri = 1,3
     dr(ri)=rr(ri,1)-rr(ri,2)
     !write(*,*)'rr:',rr(ri,1),rr(ri,2)
     !write(*,*)'dr:',dr(ri)
     r=r+dr(ri)**2
  enddo
  r=SQRT(r)
  thetax=ACOS(dr(3)/r)
  phix=ACOS(dr(1)/(dsqrt(dr(1)**2+dr(2)**2)))
  if (dr(2).lt.0)then
     phix=2.d0*ACOS(-1.d0)-phix
  end if
  x=gamma*r
 ! write(*,*)"r: ",r,'theta:',thetax,'phi:',phix
  
  do k=1,2
     lsum(k)=0.d0
     do li=1,40
        ll=li-1
        lsum(k)=lsum(k)+dsqrt((gamma**3)/((li)*(li+1)))*plaguer(ll,2.d0,x)*ysol(k,li)*EXP(-x/2.d0)
        !        write(*,*)"L:",2*k-2,"ll=",ll,"N:",dsqrt((gamma**3)/((li)*(li+1))),"plaguer:",plaguer(ll,2.d0,x),"wf:",ysol(k,li)
     enddo
     asum(:,:,k)=lsum(k)*warray(:,:,k)
  enddo
  !write(*,*)"L=0:",lsum(1),"L=2",lsum(2)
  !write(*,*)'asum:'
  do i=1,4
!     M=CEILING((i-1.d0)/2)
     M0=2-CEILING((i-1.d0)/2)
!     write(*,*)"M: ",M,"M1: ",M0
     do j=1,2
        do k=1,2
           !m1=1-2*FLOOR(0.5*(i-0.1))
           !m2=2*MOD(i,2)-1
           !write(*,'(A,1i3,1i3,A,f10.5)')'ylm(',2*k-2,M,')',AIMAG(ylm(2*k-2,M,thetax,phix)
           !write(*,*)i,j,k
           !write(*,*)asum(i,j,k)
           ysum(i,j,k) = asum(i,j,k)*ylm(2*k-2,M0,thetax,phix)
           !write(*,*)i,j,k,ysum(i,j,k)
!           write(*,*)"M:",M0,"L:",2*k-2,ylm(2*k-2,M0,thetax,phix)
        enddo  
     enddo
  enddo
  
!  do i = 1,4
!     do j=1,2
!        write(*,'(2i3)')1-2*MOD(CEILING(MOD(DBLE(i),2.d0)),2),1-2*MOD(CEILING(MOD(i/2.d0,2.d0)),2)
!        write(*,*)"real: ",REAL(ysum(i,j,1)+ysum(i,j,2)),"imaginary: ",AIMAG(ysum(i,j,1)+ysum(i,j,2))
!     enddo
  !  enddo
  do i = 1,4
     do j=1,2
        cwf(i,j)=ysum(i,j,1)+ysum(i,j,2)
        !write(*,*)cwf(i,j)
     enddo
  enddo
  
  !deut_wave%wave_out=cwf
  !write(*,*)'got to end of deut wave'
return
end subroutine deut_wave

subroutine derivative(rr,ysol,h,dcwf)
  implicit none
  real*8,intent(in)::rr(3,2),ysol(2,40),h
  complex*16,intent(out)::dcwf(4,2,a,p)
  integer::a,p,i,j
  real*8::dr(3),r,h_add(3,2),rr_new(3,2)
  complex*16::cwf_plus(4,2),cwf_minus(4,2)
  do a = 1,3
     do p = 1,2
        h_add=0.d0
        h_add(a,p)=h
        rr_new = rr+h
        call deut_wave(rr_new,cwf_plus,ysol,dr,r)
        rr_new=rr-h
        call deut_wave(rr_new,cwf_minus,ysol,dr,r)
        do i = 1,4
           do j = 1,2
              dcwf(i,j,a,p)=(cwf_plus(i,j)-cwf_minus(i,j))/(2*h)
           enddo
        enddo
     enddo
  enddo

  
  return
end subroutine derivative
