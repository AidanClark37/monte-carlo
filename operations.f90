
  subroutine rho_NNg(iarray,wf_in,N,niso,pm,exp)
    implicit none
    integer,intent(in)::iarray(niso),N,niso
    integer::p,nspin,i,j,pm,state_sign
    complex*16,intent(in)::wf_in(4,niso)
    real*8,intent(out)::exp
    complex*16,allocatable::cc_wf(:,:),tau_wf(:,:),spin_tau_wf(:,:),wf_out(:,:),tau_one_wf(:,:)
    nspin=2**N
    
    allocate(cc_wf(nspin,niso))
    allocate(tau_wf(nspin,niso))
    allocate(spin_tau_wf(nspin,niso))
    allocate(wf_out(nspin,niso))
    allocate(tau_one_wf(nspin,niso))
    exp=0.d0
    wf_out=dcmplx(0.d0,0.d0)
    do p = 1,N
    call isospin(iarray,wf_in,p,3,2,2,tau_wf)
    !tau_one_wf(:,:)=tau_wf(:,:) + wf_in(:,:)
    tau_one_wf = 0.5d0*(tau_wf + wf_in)
    !   end do
    !end do
    !write(*,*)"particle = ",p
    !do j=1,2
    !   write(*,*)j,state_sign(p,iarray(j))
    !end do
!    if(proc_rank.eq.3)then
!    do i=1,4
!       write(*,*)"spin = ",i
!       write(*,*)"wf ",(wf_in(i,j),j=1,2)
!       write(*,*)"tau_z",(tau_wf(i,j),j=1,2)
!       write(*,*)"1+tau_z wf ",(tau_one_wf(i,j),j=1,2)
!    end do
!    end if
    
    call spin(tau_one_wf,p,3,2,2,spin_tau_wf)
          wf_out = wf_out + spin_tau_wf
    
 enddo

 cc_wf = conjg(wf_in)
 
    do i = 1,nspin
       do j = 1,niso
          !if(proc_rank.eq.1)write(*,*)wf_in(i,j),cc_wf(i,j),wf_out(i,j)
          exp =	exp + cc_wf(i,j)*wf_out(i,j)
          !write(*,*)proc_rank,i,j,exp
       enddo
       
    enddo
!    write(*,*) proc_rank,exp
return

end subroutine rho_NNg

subroutine rho_NNg_other(wf_in,N,niso,pm,exp)
  use operator_calc
    implicit none
    integer::p,nspin,i,j,pm,state_sign
    integer,intent(in)::N,niso
    complex*16,intent(in)::wf_in(2**N,niso)
    real*8,intent(out)::exp
    complex*16,allocatable::cc_wf(:,:),tau_wf(:,:),spin_tau_wf(:,:),wf_out(:,:),tau_one_wf(:,:)
    nspin=2**N
    exp=0.d0
     allocate(cc_wf(nspin,niso))
    allocate(wf_out(nspin,niso))
   

    wf_out = dcmplx(0,0)
    do p = 1, N
       wf_out(:,:) = wf_out(:,:) +pm*0.5*st_wf(1,1,p,3,p) + 0.5*s_wf(:,:,p,3)
    enddo
    cc_wf = conjg(wf_in)
    do i = 1, nspin
       do j  = 1, niso
          exp = exp + wf_out(i,j)*cc_wf(i,j)
       enddo
       enddo
          return
  end subroutine rho_NNg_other

  subroutine rho_NNpTRV_NNpgPC_1(rr,dr,r,q,iarray,cwf,N,niso,exp)
    use operator_calc
    implicit none
   
    real*8,intent(in)::rr(3,2)
    complex*16,intent(in)::cwf(2**N,niso)
    integer,intent(in)::iarray(niso),N,niso
    real*8,intent(out)::exp
    integer::p,i,j,a,ii,jj,nspin
    real*8::dr(3),rhat(3),r,qr,q(3)
    complex*16::total(2**N,niso),ccwf(2**N,niso)
    nspin=2**N
    !do i = 1,3
    !   dr(i) = rr(i,1)-rr(i,2)
    !   r = r + dr(i)**2
    !enddo
    !r = dsqrt(r)
    total = dcmplx(0,0)
    do i=1,3
       rhat(i)=dr(i)/r
    enddo
    do p=1,2
       qr=0.d0
       do i = 1,3
          qr = qr + q(i)*rr(i,p)
       enddo
       
       do a = 1,3
          do ii = 1,nspin
             do jj = 1,niso
                total(ii,jj) = total(ii,jj) +((-1)**(p-1))*stt_wf(ii,jj,p,a,1,2)*dcmplx(COS(qr),SIN(qr))*rhat(a)
                
                
                
                total(ii,jj) = total(ii,jj) +((-1)**p)*stt_wf(ii,jj,p,a,2,1)*dcmplx(COS(qr),SIN(qr))*rhat(a)
           
             enddo
             
       enddo
    enddo
    
    enddo
    
    ccwf=conjg(cwf)
    exp=0.d0
    do i = 1,nspin
       do j = 1,niso
          exp=exp+ccwf(i,j)*total(i,j)
       enddo
    enddo
    
       
    
    
  end subroutine rho_NNpTRV_NNpgPC_1

  subroutine rho_NNpTRV_NNpgPC_2(rr,dr,r,q,iarray,cwf,N,niso,exp)
    
  end subroutine rho_NNpTRV_NNpgPC_2

  subroutine rho_NNpTRV_NNpgHB_0(rr,dr,r,q,iarray,cwf,N,niso,exp)
    use operator_calc
    implicit none

    real*8,intent(in)::rr(3,2)
    complex*16,intent(in)::cwf(2**N,niso)
    integer,intent(in)::iarray(niso),N,niso
    real*8,intent(out)::exp
    integer::p,i,j,a,ii,jj,nspin
    real*8::dr(3),rhat(3),r,qr,q(3)
    complex*16::total(2**N,niso),ccwf(2**N,niso)
    nspin=2**N
    total = dcmplx(0.d0,0.d0)
    do p = 1,2
       do a = 1,3
          do ii = 1,nspin
             do jj= 1,niso
                do i=1,3
                   total(ii,jj)=total(ii,jj)+stt_wf(ii,jj,p,a,i,i)*q(a)
                enddo
                total(ii,jj)=total(ii,jj)+st_wf(ii,jj,p,a,p)*q(a)
             enddo
          enddo
       enddo
    enddo
    
    ccwf=conjg(cwf)
    exp=0.d0
    do i = 1,nspin
       do j = 1,niso
          exp=exp+ccwf(i,j)*total(i,j)
       enddo
    enddo

          return      
       
  end subroutine rho_NNpTRV_NNpgHB_0

  subroutine rho_NNpTRV_NNpgHB_1(rr,dr,r,q,iarray,cwf,N,niso,exp)
    use operator_calc
    implicit none

    real*8,intent(in)::rr(3,2)
    complex*16,intent(in)::cwf(2**N,niso)
    integer,intent(in)::iarray(niso),N,niso
    real*8,intent(out)::exp
    integer::p,i,j,a,ii,jj,nspin
    real*8::dr(3),rhat(3),r,qr,q(3)
    complex*16::total(2**N,niso),ccwf(2**N,niso)
    nspin=2**N
    total = dcmplx(0.d0,0.d0)
    do p = 1,2
       do a = 1,2
          do ii=1,nspin
             do jj=1,niso
                total(ii,jj)=total(ii,jj) + (st_wf(ii,jj,p,a,p) + s_wf(ii,jj,p,a))*q(a)
             enddo
          enddo
       enddo
    enddo
    ccwf=conjg(cwf)
    exp=0.d0
    do i = 1,nspin
       do j = 1,niso
          exp=exp+ccwf(i,j)*total(i,j)
       enddo
    enddo
    return
  end subroutine rho_NNpTRV_NNpgHB_1

  subroutine rho_NNpTRV_NNpgHB_2(rr,dr,r,q,iarray,cwf,N,niso,exp)
    use operator_calc
    implicit none

    real*8,intent(in)::rr(3,2)
    complex*16,intent(in)::cwf(2**N,niso)
    integer,intent(in)::iarray(niso),N,niso
    real*8,intent(out)::exp
    integer::p,i,j,a,ii,jj,nspin
    real*8::dr(3),rhat(3),r,qr,q(3)
    complex*16::total(2**N,niso),ccwf(2**N,niso)
    nspin=2**N
    total = dcmplx(0.d0,0.d0)
    do p = 1,2
       do a = 1,3
          do ii=1,nspin
             do jj= 1, niso
                total(ii,jj) = total(ii,jj) + (stt_wf(ii,jj,p,a,3,3) + st_wf(ii,jj,p,a,p))*q(a)
             enddo
          enddo
       enddo
    enddo
    
    ccwf=conjg(cwf)
    exp=0.d0
    do i = 1,nspin
       do j = 1,niso
          exp=exp+ccwf(i,j)*total(i,j)
       enddo
    enddo
    return
  end subroutine rho_NNpTRV_NNpgHB_2
  
  subroutine rho_NNpTRV_ppg_PC(rr,dr,r,q,iarray,cwf,N,niso,exp)
    use operator_calc
    implicit none

    real*8,intent(in)::rr(3,2)
    complex*16,intent(in)::cwf(2**N,niso)
    integer,intent(in)::iarray(niso),N,niso
    real*8,intent(out)::exp
    integer::p,i,j,a,ii,jj,nspin
    real*8::dr(3),rhat(3),r,qr,q(3)
    complex*16::total(2**N,niso),ccwf(2**N,niso)
    nspin=2**N
    total = dcmplx(0.d0,0.d0)
    do i=1,3
       rhat(i)=dr(i)/r
    enddo
    do p = 1,2
       do a = 1, 3
          do ii = 1,nspin
             do jj=1,niso
                do i = 1,2
                   total(ii,jj) = total(ii,jj) + stt_wf(ii,jj,p,a,i,i)*rhat(a)
                enddo
             enddo
          enddo
       enddo
    enddo
    exp=0.d0
    do i = 1,nspin
       do j = 1,niso
          exp=exp+ccwf(i,j)*total(i,j)
       enddo
    enddo
    return

  end subroutine rho_NNpTRV_ppg_PC

  subroutine tau_exp_val(iarray,cwf,b,N,niso,nspin,tau_exp)
    implicit none
    integer,intent(in)::iarray(niso),b,N,niso,nspin
    complex*16,intent(in)::cwf(nspin,niso)
    integer::p,i,j
    real*8::tau_exp
    complex*16,allocatable::cc_wf(:,:),tau_wf(:,:)
    allocate(cc_wf(nspin,niso))
    allocate(tau_wf(nspin,niso))
    do p = 1,N
       call isospin(iarray,cwf,p,3,N,niso,tau_wf)
    enddo
    cc_wf=conjg(cwf)
    tau_exp=0.d0
    do i = 1,nspin
       do j = 1,niso
          tau_exp=tau_exp+cc_wf(i,j)*tau_wf(i,j)
          !write(*,*)tau_exp
       enddo
    enddo
    return
  end subroutine tau_exp_val
  


  subroutine spin_exp_val(cwf,b,N,niso,nspin,spin_exp)
    implicit none
    integer,intent(in)::b,N,niso,nspin
    complex*16,intent(in)::cwf(nspin,niso)
    complex*16::spin_cwf(nspin,niso),cc_cwf(nspin,niso)
    real*8,intent(out)::spin_exp
    integer::p,i,j
    spin_exp=0
    do p = 1,N
       call spin(cwf,p,b,N,niso,spin_cwf)
       cc_cwf=conjg(cwf)
       do i = 1,nspin
          do j=1,niso
             spin_exp=spin_exp+cc_cwf(i,j)*spin_cwf(i,j)
 !            write(*,*)p,i,j,cwf(i,j),spin_cwf(i,j)
          enddo
       enddo
    enddo
!    stop
  end subroutine spin_exp_val


