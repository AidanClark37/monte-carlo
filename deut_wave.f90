subroutine deut_wave(rr,cwf)
  use input
  use storing
  use wigner
  implicit none

  real*8::dr(3)
  real*8::rr(2,3),cwf
  integer::i,j,k,li,ll,M,nla,M0,npart,ri
  real*8::warray(4,2,2),asum(4,2,2),ysol(2,100),gamma,lsum(2),r,plaguer,thetax,phix,x,prod
  complex*16::ysum(4,2,2),ylm

  gamma=deut_par%gamma 
  warray=deut_wave%pre_wave !angular component of wavefunction
  nla=deut_par%nla    !number of laguerre polynomials used
  ysol=deut_wave%ysol !wavefunction coefficients
  r=0
  do ri = 1,3
     dr(ri)=rr(1,ri)-rr(2,ri)
     r=r+dr**2
  enddo
  r=SQRT(r)
  thetax=ACOS(dr(3)/r)
  phix=ACOS(dr(1)/(dsqrt(dr(1)**2+dr(2)**2)))
  if (dr(2).lt.0)then
     phix=2.d0*ACOS(-1.d0)-phix
  end if
  x=gamma*r
  !write(*,*)"r: ",r
  
  lsum=0.d0
  do k=1,2
     do li=1,nla
        ll=li-1
        lsum(k)=lsum(k)+dsqrt((gamma**3)/((li)*(li+1)))*plaguer(ll,2.d0,x)*ysol(k,li)*EXP(-x/2.d0)
!        write(*,*)"L:",2*k-2,"ll=",ll,"N:",dsqrt((gamma**3)/((li)*(li+1))),"plaguer:",plaguer(ll,2.d0,x),"ysol:",ysol(k,li)
     enddo
     asum(:,:,k)=lsum(k)*warray(:,:,k)
  enddo
  write(*,*)"L=0:",lsum(1),"L=2",lsum(2)

  do i=1,4
!     M=CEILING((i-1.d0)/2)
     M0=2-CEILING((i-1.d0)/2)
!     write(*,*)"M: ",M,"M1: ",M0
     do j=1,2
        do k=1,2
           !m1=1-2*FLOOR(0.5*(i-0.1))
           !m2=2*MOD(i,2)-1
           !write(*,'(A,1i3,1i3,A,f10.5)')'ylm(',2*k-2,M,')',AIMAG(ylm(2*k-2,M,thetax,phix))
           ysum(i,j,k) = asum(i,j,k)*ylm(2*k-2,M0,thetax,phix)
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

  cwf=ysum(:,:,1)+ysum(:,:,2)
  deut_wave%wave_out=cwf
return
end subroutine deut_wave
