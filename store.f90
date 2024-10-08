module storing
  implicit none

  !kinematic table
  integer::n_ene
  real*8,allocatable::ene(:)
  real*8,allocatable::w_ene(:)
  real*8,allocatable::p_rel(:)
  real*8,allocatable::v_rel(:)
  real*8,allocatable::q_ex(:)
  real*8,allocatable::q4_ex(:)
  
  !deuteron wave function
  type d_wave
     integer::neq
     integer::alf0
     integer::nnl
     integer::l(2),s(2),t(2),j(2)
     real*8::ysol(2,100)
     real*8::pre_wave(4,2,2)
     complex*16::wave_out(4,2)
     real*8::tot_en
     real*8::kin_en
     real*8::pot_en
     real*8::pws
     real*8::pwd
     real*8::edm
  end type d_wave
  type(d_wave)::deut_wave

  type wave_data
     integer::l(2)
     integer::s(2)
     integer::t(2)
     integer::j
     integer::jp
     integer::ich
  end type wave_data
  type(wave_data)::wave_qn(20)
  integer::num_waves
  
  type wave_scatt_out
     integer::l
     integer::s
     integer::t
     integer::j
     integer::jz
     integer::st
     integer::nwqn
  end type wave_scatt_out
  type(wave_scatt_out),allocatable::waves_out(:)
  integer::num_waves_out
  
  type nn_wave
     integer::neq
     real*8 ::gamma
     integer::nnl
     real*8 ::eps
     integer::l(2),s(2),t(2),j(2)
     real*8 ::ac(2,200)
     real*8 ::rd(2,2)
     real*8 ::mix,d1,d2
  end type nn_wave
  type(nn_wave),allocatable::wave_proc(:)
  type(nn_wave),allocatable::scatt_wave(:)

  type wksi_grid
     integer::ir
     integer::ithe
     integer::iphi
  end type wksi_grid
  type(wksi_grid),allocatable::grid_proc(:)

  type wksi_grid_2
     integer::jz
     integer::iw
     integer::iq
     integer::i
     integer::j
  end type wksi_grid_2
  type(wksi_grid_2),allocatable::grid_proc_2(:)

  !MATRIX ELEMENTS
  complex*16,allocatable::cww(:,:,:,:,:)
  complex*16,allocatable::cww0(:,:,:,:,:)

  !FINAL RESULTS
  real*8,allocatable::dgdp0(:,:)
  real*8,allocatable::dgdp1(:,:)
  real*8,allocatable::tot_dgdp(:,:)
  real*8,allocatable::err2_dgdp(:,:,:)
  real*8,allocatable::tot_err_dgdp(:,:)
  real*8,allocatable::xsec0(:)
  real*8,allocatable::xsec1(:)

  real*8,allocatable::tot_xsec(:)
  real*8,allocatable::err2_xsec(:,:)
  real*8,allocatable::tot_err_xsec(:)
  !SENSITIVITY STUDIES
  real*8::ra_2(100),c_DD(100)
  real*8,allocatable::sense_xsec(:,:)
  
contains
    subroutine set_waves(jmax,tz,num_waves_comp)
      implicit none
      integer,intent(in)::jmax
      integer,intent(in)::tz
      integer,intent(in)::num_waves_comp
      integer::ii,j,l,s,t
      integer::j0,jp0,jp
      integer::jz,i,ich
      
      ii=0
      j0=-1;jp0=-1
      do j=0,jmax
         do s=0,1
         do t=abs(tz),1
         do l=abs(j-s),abs(j+s)      
            if((-1)**(l+s+t).eq.-1)then
               jp=((-1)**l+1)/2
               if(j.ne.j0.or.jp.ne.jp0)then
                  ii=ii+1
                  wave_qn(ii)%j=j
                  wave_qn(ii)%jp=jp
                  wave_qn(ii)%ich =1
                  wave_qn(ii)%s(1)=s
                  wave_qn(ii)%l(1)=l
                  wave_qn(ii)%t(1)=t
                  j0=j;jp0=jp
               else
                  wave_qn(ii)%ich =2
                  wave_qn(ii)%s(2)=s
                  wave_qn(ii)%l(2)=l
                  wave_qn(ii)%t(2)=t
                  j0=-1;jp0=-1
               end if
            end if
         end do
         end do
         end do
      end do
      if(ii< num_waves_comp)num_waves=ii
      if(ii>=num_waves_comp)num_waves=num_waves_comp

      num_waves_out=0
      do i=1,num_waves
         num_waves_out=num_waves_out+(2*wave_qn(i)%j+1)*wave_qn(i)%ich
      end do
      allocate(waves_out(num_waves_out))
      ii=0
      do i=1,num_waves
      do j=1,wave_qn(i)%ich
      do jz=-wave_qn(i)%j,wave_qn(i)%j
         ii=ii+1
         waves_out(ii)%l=wave_qn(i)%l(j)
         waves_out(ii)%s=wave_qn(i)%s(j)
         waves_out(ii)%t=wave_qn(i)%t(j)
         waves_out(ii)%j=wave_qn(i)%j
         waves_out(ii)%jz=jz
         waves_out(ii)%st=j
         waves_out(ii)%nwqn=i
      end do
      end do
      end do
      return
  end subroutine set_waves

  subroutine allocate_grid(ndim)
    implicit none
    integer,intent(in)::ndim
    allocate(grid_proc(ndim))
  end subroutine allocate_grid

  subroutine allocate_grid_2(ndim)
    implicit none
    integer,intent(in)::ndim
    allocate(grid_proc_2(ndim))
  end subroutine allocate_grid_2

  subroutine allocate_cww(nd1,nd2,nd3,nd4,nd5)
    implicit none
    integer,intent(in)::nd1,nd2,nd3,nd4,nd5
    allocate(cww(nd1,nd2,nd3,nd4,nd5))
  end subroutine allocate_cww
end module storing
  
