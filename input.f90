module input
  implicit none

  character(len=10)::process
  integer::proc_id
  
  type current
     integer::ilem
     integer::ilcc
     integer::ireg
     integer::iord
     integer::icrg
     integer::iewk
     real*8 ::r_l
  end type current

  type potential
     integer::ipot
     integer::ipte
     integer::ipls
     integer::ipbb
     integer::ipqq
     integer::ipsb
     integer::ipco
     integer::lemp
     integer::ilb
     character(len=10)::name
  end type potential
  
  type in_deut
     integer::nla
     integer::nn
     real*8::gamma
  end type in_deut

  type in_scatt
     integer::tz
     integer::jmax
     real*8::eps
     real*8::gamma
     integer::nla
     integer::nn
     integer::num_waves_comp
  end type in_scatt

  type in_wksi
     real*8 ::alf
     integer::nrr
     integer::nthe
     integer::nphi
  end type in_wksi

  type in_rate
     integer::n_phat
  end type in_rate

  logical::ierr
  logical::isense
  logical::variational
  
  type(current)::cur_par
  type(potential)::pot_par
  type(in_deut)::deut_par
  type(in_scatt)::scatt_par
  type(in_wksi)::wksi_par
  type(in_rate)::rate_par

contains

  subroutine read_input()
    implicit none

    read(5,*)process
    read(5,*)variational
    
    !current_parameters
    read(5,*)cur_par%ilcc
    read(5,*)cur_par%iord
!    read(5,*)cur_par%r_l
    
    !potential parameters
    read(5,*)pot_par%ipot
    read(5,*)pot_par%ipte
    read(5,*)pot_par%ipls
    read(5,*)pot_par%ipqq
    read(5,*)pot_par%ipbb
    read(5,*)pot_par%ipsb
    read(5,*)pot_par%ipco
    read(5,*)pot_par%lemp       

    !deuteron parameters
    read(5,*)deut_par%nla
    read(5,*)deut_par%nn
    read(5,*)deut_par%gamma

    !scattering wave parameters
    read(5,*)scatt_par%jmax
    read(5,*)scatt_par%eps
    read(5,*)scatt_par%gamma
    read(5,*)scatt_par%nla
    read(5,*)scatt_par%nn
    read(5,*)scatt_par%num_waves_comp
    
    !wksi parameters
    read(5,*)wksi_par%nrr
    read(5,*)wksi_par%alf
    read(5,*)wksi_par%nthe
    read(5,*)wksi_par%nphi

    !integration on theta_p
    read(5,*)rate_par%n_phat

    !for computing error
    read(5,*)ierr

    !for sensitivity studies
    read(5,*)isense
    write(*,*)isense
    call set_par()
    
  end subroutine read_input


  subroutine set_par()
    use constant_par
    implicit none

    if(trim(process).eq."mucap")then
       scatt_par%tz=-1
       cur_par%icrg=2
       cur_par%iewk=2
       proc_id=1
    else if(trim(process).eq."pp")then
       scatt_par%tz=1
       cur_par%icrg=2
       cur_par%iewk=2
       proc_id=2
    else if(trim(process).eq."npdg")then
       scatt_par%tz=0
       cur_par%icrg=0
       cur_par%iewk=0
       proc_id=3
    else
       proc_id=0
       write(*,*)"The process is not yet considered"
       stop
    end if

    if(scatt_par%tz.eq.-1)htm=41.44252d0
    if(scatt_par%tz.eq. 0)htm=41.47108d0
    if(scatt_par%tz.eq.+1)htm=41.49964d0

    if(pot_par%ipot<100)then
       if ( pot_par%ipot == 31) then
          cur_par%ilem=0
          pot_par%ilb=0
          pot_par%name='AV18'
          cur_par%ireg=0
       end if
       if ( pot_par%ipot == 61) then
          cur_par%ilem=16
          pot_par%ilb=6
          pot_par%name='NvIa'
          cur_par%ireg=1
       end if
       if ( pot_par%ipot == 62) then
          cur_par%ilem=19
          pot_par%ilb=10
          pot_par%name='NvIIa'
          cur_par%ireg=1
       end if
       if ( pot_par%ipot == 63) then
          cur_par%ilem=17
          pot_par%ilb=5
          pot_par%name='NvIb'
          cur_par%ireg=2
       end if
       if ( pot_par%ipot == 64) then
          cur_par%ilem=20
          pot_par%ilb=9
          pot_par%name='NvIIb'
          cur_par%ireg=2
       end if
    else if ( pot_par%ipot > 100 .and. pot_par%ipot<1000 ) then
       pot_par%ilb=0
       cur_par%ireg=0
       if(pot_par%ipot.eq.101)then
          pot_par%name='EM LO450'
          cur_par%ilem=22
       end if
       if(pot_par%ipot.eq.102)then
          pot_par%name='EM NLO450'
          cur_par%ilem=25
       end if
       if(pot_par%ipot.eq.103)then
          pot_par%name='EM N2LO450'
          cur_par%ilem=28
       end if
       if(pot_par%ipot.eq.104)then
          pot_par%name='EM N3LO450'
          cur_par%ilem=31
       end if
       if(pot_par%ipot.eq.105)then
          pot_par%name='EM N4LO450'
          cur_par%ilem=34
       end if
       if(pot_par%ipot.eq.106)then
          pot_par%name='EM LO500'
          cur_par%ilem=23
       end if
       if(pot_par%ipot.eq.107)then
          pot_par%name='EM NLO500'
          cur_par%ilem=26
       end if
       if(pot_par%ipot.eq.108)then
          pot_par%name='EM N2LO500'
          cur_par%ilem=29
       end if
       if(pot_par%ipot.eq.109)then
          pot_par%name='EM N3LO500'
          cur_par%ilem=32
       end if
       if(pot_par%ipot.eq.110)then
          pot_par%name='EM N4LO500'
          cur_par%ilem=35
       end if
       if(pot_par%ipot.eq.111)then
          pot_par%name='EM LO550'
          cur_par%ilem=24
       end if
       if(pot_par%ipot.eq.112)then
          pot_par%name='EM NLO550'
          cur_par%ilem=27
       end if
       if(pot_par%ipot.eq.113)then
          pot_par%name='EM N2LO550'
          cur_par%ilem=30
       end if
       if(pot_par%ipot.eq.114)then
          pot_par%name='EM N3LO550'
          cur_par%ilem=33
       end if
       if(pot_par%ipot.eq.115)then
          pot_par%name='EM N4LO550'
          cur_par%ilem=36
       end if
       if(pot_par%ipot.eq.121)then
          pot_par%name='EM N3LO500(2003)'
          cur_par%ilem=8
       end if
    else if(pot_par%ipot>1000)then
          pot_par%name='ChiEFT LosAlamos'
          write(*,*)"Setting LO of the currents"
          cur_par%iord=0
          cur_par%ilem=1000
          cur_par%ireg=1
          cur_par%ilem=16
    else
       write(*,*)"Currents are not available"
       write(*,*)"Stopping the program"
       stop
    end if
  end subroutine set_par
  
end module input
