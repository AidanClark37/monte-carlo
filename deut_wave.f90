module deuteron
  implicit none
  real*8::lsum(2)
  contains
subroutine deut_wave(rr,cwf,ysol,dr,r)
  use wigner
  use mpi_modules
  use pre_deut
  use param_calc
  implicit none
  real*8::dr(3),ysol(2,40),rr(3,2)
  real*8,parameter::gamma=4.5d0
  integer::i,j,k,li,lj,M,M0,ri
  real*8::asum(4,2,2),r,plaguer,thetax,phix,x,prod
  complex*16::ysum(4,2,2),ylm
  complex*16,intent(inout)::cwf(:,:)
  ! r = r_1-r_2 vector, labelled dr
  cwf=dcmplx(0.d0,0.d0)
  ! |r|, labeled r
  r=0.d0
  do ri = 1,3
     dr(ri)=rr(ri,1)-rr(ri,2)
     r=r+dr(ri)**2
  enddo
  r=SQRT(r)
  !calculate theta and phi for r vector input
  thetax=ACOS(dr(3)/r)
  phix=ACOS(dr(1)/(dsqrt(dr(1)**2+dr(2)**2)))
  if (dr(2).lt.0)then
     phix=2.d0*ACOS(-1.d0)-phix
  end if
  x=gamma*r

  !radial wavefunction, f_0 and f_2 part
  do k=1,2
     lsum(k)=0.d0
     do li=1,40
        lj=li-1
        lsum(k)=lsum(k)+dsqrt((gamma**3)/((li)*(li+1)))*plaguer(lj,2.d0,x)*ysol(k,li)*EXP(-x/2.d0)
     enddo
     asum(:,:,k)=lsum(k)*warray(:,:,k)

  enddo

  ! dependance on sperical harmonics
  do i=1,4

     M0=2-CEILING((i-1.d0)/2)
     do j=1,2
        do k=1,2
           ysum(i,j,k) = asum(i,j,k)*ylm(2*k-2,M0,thetax,phix)
        enddo  
     enddo
  enddo
  !sum over L
  do i = 1,4
     do j=1,2
        cwf(i,j)=ysum(i,j,1)+ysum(i,j,2)
       
     enddo
  enddo

return
end subroutine deut_wave

subroutine derivative_old(rr,ysol,h,dcwf)
  implicit none
  real*8,intent(in)::rr(3,2),ysol(2,40),h
  complex*16,intent(out)::dcwf(4,2,3,2)
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
end subroutine derivative_old
end module deuteron
