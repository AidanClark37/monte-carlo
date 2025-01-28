module mpi_modules
  use mpi
  implicit none


  type input
     real*8            ::wf(2,40)
     integer           ::nwalk
     integer           ::neq
     integer           ::nav
     integer           ::ncorr
     real*8            ::sigma
     integer           ::nla
     integer           ::iarray(2)
     real*8            ::mass
     real*8            ::lambda

  end type input

  type(input)::msg

  real*8::wf(2,40)   !wave function file input
  integer           ::nwalk     !number of walkers
  integer           ::neq       !number of points for thermalization (each walker)
  integer           ::nav       !number of point for averaging
  integer           ::ncorr     !number of point to reduce correlation
  real*8            ::sigma     !length step of the move
  integer           ::nla
  integer           ::iarray(2)
  real*8            ::lambda
  real*8            ::mass

  !general settings
  integer::nwalks_for_proc
  

  !mpi variables
  integer::mpi_procs
  integer::proc_rank
  integer::ierr
  integer::mpi_world

  !mpi type 
  integer::num_blk
  integer::v_len_blk(6)
  integer::v_head(6)
  integer::v_el_typ(6)
  integer::new_type

contains
  subroutine start_mpi()
    implicit none

    call mpi_init(ierr)
    mpi_world = MPI_COMM_WORLD
    call mpi_comm_rank(mpi_world, proc_rank, ierr)
    call mpi_comm_size(mpi_world, mpi_procs, ierr)
    return
  end subroutine start_mpi

  subroutine end_mpi()
    implicit none
   
    call mpi_finalize(ierr)
    return
  end subroutine end_mpi

  subroutine general_setting()
    implicit none

    if(mod(nwalk,mpi_procs).ne.0)then
       write(*,*)'The number of walkers is not a multiple of mpi_procs'
       write(*,*)'Program terminated'
       write(*,*)'Please select a multiple of = ',mpi_procs
       stop
    end if

    nwalks_for_proc=nwalk/mpi_procs
    if(proc_rank.eq.1)then
       write(*,*)"Number of walkers for processor: ",nwalks_for_proc
    end if
    return
  end subroutine general_setting
  
  
  subroutine mpi_broadcast_input()
    implicit none
    integer::i
    
    num_blk=6
    v_len_blk=[80,1,1,1,1,1]
    v_head   =[0,640,644,648,652,656]
    v_el_typ=[MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_INTEGER,&
         MPI_INTEGER,MPI_INTEGER,MPI_DOUBLE_PRECISION]

    call mpi_type_struct(num_blk,v_len_blk,v_head,v_el_typ,new_type,ierr)

    call mpi_type_commit(new_type,ierr)
    
    call mpi_bcast(msg,1,new_type,0,mpi_world,ierr)

!    call mpi_scatter(rnd_mpi,1,MPI_DOUBLE_PRECISION,&
!         rnd0,1,MPI_DOUBLE_PRECISION,0,mpi_world,ierr)

   ! wf_file =msg%wf_file
    wf      =msg%wf
    nwalk   =msg%nwalk      
    neq     =msg%neq           
    nav     =msg%nav       
    ncorr   =msg%ncorr     
    sigma   =msg%sigma
 
    if(proc_rank.eq.1)then
       write(*,*)wf(1,1)
       write(*,*)nwalk
       write(*,*)neq
       write(*,*)nav
       write(*,*)ncorr
       write(*,*)sigma
    end if
    
    return
  end subroutine mpi_broadcast_input

   subroutine op_reduction(op,op2,nqq,mean,st_dev)
     implicit none     

     integer,intent(in) ::nqq
     real(8),intent(in) ::op (nqq)
     real(8),intent(in) ::op2(nqq)
     real(8)::sum_op (nqq)
     real(8)::sum_op2(nqq)
     real(8),intent(out)::mean (nqq)
     real(8),intent(out)::st_dev(nqq)

     integer::ndim

     ndim=nqq
     
     call mpi_reduce(op ,sum_op ,ndim,MPI_DOUBLE_PRECISION,&
          MPI_SUM,0,mpi_world,ierr)
     call mpi_reduce(op2,sum_op2,ndim,MPI_DOUBLE_PRECISION,&
          MPI_SUM,0,mpi_world,ierr)

     
     mean  =sum_op/nwalk
     sum_op2=sum_op2/nwalk
     st_dev=sqrt((sum_op2-mean**2)/nwalk)
     
     return
   end subroutine op_reduction

   subroutine int_reduction(ii,sum_ii)
     implicit none
     
     integer,intent(in) ::ii
     integer,intent(out)::sum_ii

     call mpi_reduce(ii,sum_ii,1,MPI_INTEGER,MPI_SUM,0,mpi_world,ierr)
     
     return
   end subroutine int_reduction
   
 end module mpi_modules
    
