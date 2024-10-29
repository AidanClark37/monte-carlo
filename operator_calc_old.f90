module operator_calc
  implicit none
  complex*16::operator_array(4,2,3,3,3,3,3,3,3,3)
contains
  subroutine start_operator(iarray,wf_in,N,niso)
    integer,intent(in)::n,niso,iarray(niso)
    integer::,tp1,tp2,sp1,sp2,ta1,ta2,sa1,sa2,nspin
    complex*16,intent(in)::wf_in
    complex*16,allocatable::wf1(:,:),wf2(:,:),wf3(:,:),wf4(:,:)
    nspin=2**N
    allocate(wf1(nspin,niso))
    allocate(wf2(nspin,niso))
    allocate(wf3(nspin,niso))
    allocate(wf4(nspin,niso))
    operator_array(:,:,1,1,1,1,1,1,1,1)=wf_in(:,:)
    do tp1 = 1,N
       do ta1= 1,3
          call isospin(iarray,wf_in,tp1,ta1,N,niso,wf1)
          operator_array(:,:,tp1+1,ta1,1,1,1,1,1,1)=wf1(:,:)
          do tp2=2,N+1
             do ta2=1,N
                call isospin(iarray,wf1,tp2,ta2,N,niso,wf2)
                operator_array(:,:,tp1+1,ta1,tp2+1,ta2,1,1,1,1)=wf2(:,:)
                do sp1=1,N
                   do sa1=1,3
                      call spin(wf2,sp1,sa1,N,niso,wf3)
                      operator_array(:,:,tp1+1,ta1,tp2+1,ta2,sp1+1,sa1,1,1)=wf3(:,:)
                      do sp2=1,N
                         do sa2=1,3
                            call spin(wf3,sp2,sa2,N,niso,wf4)
                            operator_array(:,:,tp1+1,ta1,tp2+1,ta2,sp1+1,sa1,sp2+1,sa2)=wf4(:,:)
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          do sp1=1,N
             do sa1=1,3
                call spin(wf1,sp1,sa1,N,niso,wf2)
                operator_array(:,:,tp1+1,ta1,1,1,sp1+1,sa1,1,1)=wf2(:,:)
                do sp2=1,N
                   do  sa2=1,3
                      call spin(wf2,sp2,sa2,N,niso,wf3)
                      operator_array(:,:,tp1+1,ta1,1,1,sp1+1,sa1,sp2+1,sa2)=wf3(:,:)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    do sp1=1,N
       do sa1= 1,3
          call spin(wf_in,sp1,sa1,N,niso,wf1)
          operator_array(:,:,1,1,1,1,sp1+1,sa1,1,1)=wf1(:,:)
          do sp1=1,N
             do sa1=1,3
                call spin(wf1,sp2,sa2,N,niso,wf1)
                operator_array(:,:,1,1,1,1,sp1+1,sa1,sp2+1,sa2)=wf2(:,:)
             enddo
          enddo
       enddo
    enddo
    
  end subroutine start_operator
  
  !spin(wf,p,b,N,niso,sigma_wf)
!isospin(iarray,wf,p,b,N,niso,tau_wf)
  
  subroutine operator_call(t,s,N,niso,wf_out)
    integer,intent(in)::t(2,2),s(2,2)
    integer::nspin,i,j
    complex*16,intent(out)::wf_out(nspin,niso)
    nspin=2**N
    do i= 1, nspin
       do j = 1, niso
          if(tn==0) then
             if(sn==1) then
                wf_out(i,j)=operator_array(i,j,t(1,1)+1,t(1,2)+1,t(2,1)+1,t(2,2)+1,s(1,1)+1,s(1,2)+1,s(2,1)+1,s(2,2)+1)               
       enddo
    enddo
    
end module operator_calc







