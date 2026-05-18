subroutine map65_mmdec(nutc,id2,nqd,nsubmode,nfa,nfb,nfqso,ntol,newdat,   &
     nagain,max_drift,ndepth,mycall,hiscall,hisgrid)

  use prog_args
  use timer_module, only: timer
  use q65_decode
  use decode_callbacks, only: map65_q65_decoder, map65_q65_decoded_callback

  include 'jt9com.f90'
  include 'timer_common.inc'

  logical single_decode,bVHF,lnewdat,lagain,lclearave,lapcqonly
  integer*2 id2(300*12000)
  integer nqf(20)
!  type(params_block) :: params
  character(len=12) :: mycall, hiscall
  character(len=6) :: hisgrid
  data ntr0/-1/
  save
  type(map65_q65_decoder) :: my_q65

! Cast C character arrays to Fortran character strings
!  datetime=transfer(params%datetime, datetime)
!  mycall=transfer(params%mycall,mycall)
!  hiscall=transfer(params%hiscall,hiscall)
!  mygrid=transfer(params%mygrid,mygrid)
!  hisgrid=transfer(params%hisgrid,hisgrid)

  my_q65%decoded = 0
  ncontest=0
  nQSOprogress=0
  lclearave=.false.
  single_decode=.false.
  lapcqonly=.false.
  lnewdat=(newdat.ne.0)
  lagain=(nagain.ne.0)
  bVHF=.true.
  emedelay=2.5
  ntrperiod=60

  call timer('dec_q65 ',0)
  call my_q65%decode(map65_q65_decoded_callback,id2,nqd,nutc,ntrperiod,nsubmode,nfqso,       &
       ntol,ndepth,nfa,nfb,lclearave,single_decode,lagain,max_drift,lnewdat,  &
       emedelay,mycall,hiscall,hisgrid,nQSOProgress,ncontest,lapcqonly,navg0,nqf)
  call timer('dec_q65 ',1)

  return

end subroutine map65_mmdec
