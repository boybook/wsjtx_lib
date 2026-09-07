subroutine multimode_decoder(ss,id2,params,nfsample)

!$ use omp_lib
  use prog_args
  use timer_module, only: timer
  use ft8_decode, only: c2fox, g2fox, nsnrfox, nfreqfox, n30fox, &
       n30z, nfox
  use decode_callbacks, only: counting_jt4_decoder, counting_jt65_decoder, &
       counting_jt9_decoder, counting_ft8_decoder, counting_ft4_decoder, &
       counting_fst4_decoder, counting_q65_decoder, &
       jt4_decoded_callback, jt4_average_callback, jt65_decoded_callback, &
       jt9_decoded_callback, ft8_decoded_callback, ft4_decoded_callback, &
       fst4_decoded_callback, q65_decoded_callback

  include 'jt9com.f90'

  real ss(184,NSMAX)
  logical baddata,newdat65,newdat9,single_decode,bVHF,bad0,newdat,ex,ltry_a8
  logical lprinthash22
  integer*2 id2(NTMAX*12000)
  integer nqf(20)
  type(params_block) :: params
  real*4 dd(NTMAX*12000)
  character(len=20) :: datetime
  character(len=12) :: mycall, hiscall
  character(len=6) :: mygrid, hisgrid
  character*60 line
  external wsjtx_decode_stats
  data ndec8/0/,ntr0/-1/,ndec41/0/,ndec47/0/
  save
  type(counting_jt4_decoder) :: my_jt4
  type(counting_jt65_decoder) :: my_jt65
  type(counting_jt9_decoder) :: my_jt9
  type(counting_ft8_decoder) :: my_ft8
  type(counting_ft4_decoder) :: my_ft4
  type(counting_fst4_decoder) :: my_fst4
  type(counting_q65_decoder) :: my_q65

  my_jt4%decoded = 0
  my_jt65%decoded = 0
  my_jt9%decoded = 0
  my_ft8%decoded = 0
  my_ft4%decoded = 0
  my_fst4%decoded = 0
  my_q65%decoded = 0
  nsynced=0
  navg0=0
  if(params%nmode.eq.8 .and. params%nzhsym.eq.41) then
     ndec41=0
     ndec47=0
  endif

  if(.not.params%newdat .and. params%ntr.gt.ntr0) go to 800
  ntr0=params%ntr
  rms=sqrt(dot_product(float(id2(1:180000)),                         &
       float(id2(1:180000)))/180000.0)
  ! Match WSJT-X 3.0.2: the input gate is 0.5 counts. Higher thresholds
  ! discard low-gain recordings before FT8/FT4 can attempt a decode.
  if(rms.lt.0.5) go to 800

  !cast C character arrays to Fortran character strings
  datetime=transfer(params%datetime, datetime)
  mycall=transfer(params%mycall,mycall)
  hiscall=transfer(params%hiscall,hiscall)
  mygrid=transfer(params%mygrid,mygrid)
  hisgrid=transfer(params%hisgrid,hisgrid)

  ! initialize decode counts
  my_jt4%decoded = 0
  my_jt65%decoded = 0
  my_jt9%decoded = 0
  my_ft8%decoded = 0
  my_ft4%decoded = 0
  my_fst4%decoded = 0
  my_q65%decoded = 0
  
! For testing only: return Rx messages stored in a file as decodes
  inquire(file='rx_messages.txt',exist=ex)
  if(ex) then
     if(params%nzhsym.eq.41) then
        open(39,file='rx_messages.txt',status='old')
        do i=1,9999
           read(39,'(a60)',end=5) line
           if(line(1:1).eq.' ' .or. line(1:1).eq.'-') go to 800
           write(*,'(a)') trim(line)
        enddo
5       close(39)
     endif
     go to 800
  endif

  ncontest=iand(params%nexp_decode,7)
  single_decode=iand(params%nexp_decode,32).ne.0
  bVHF=iand(params%nexp_decode,64).ne.0
  if(mod(params%nranera,2).eq.0) ntrials=10**(params%nranera/2)
  if(mod(params%nranera,2).eq.1) ntrials=3*10**(params%nranera/2)
  if(params%nranera.eq.0) ntrials=0

  my_jt4%nutc = params%nutc
  my_jt65%nutc = params%nutc
  my_jt65%bVHF = bVHF
  my_jt9%nutc = params%nutc
  my_ft8%nutc = params%nutc
  my_ft8%ncontest = ncontest
  my_ft8%mycall = mycall
  my_ft4%nutc = params%nutc
  
  nfail=0
!10 if (params%nagain) then
!     open(13,file=trim(temp_dir)//'/decoded.txt',status='unknown',            &
!          position='append',iostat=ios13)
!  else
!     open(13,file=trim(temp_dir)//'/decoded.txt',status='unknown',iostat=ios13)
!  endif
!  if(ios13.ne.0) then
!     nfail=nfail+1
!     if(nfail.le.3) then
!        call sleep_msec(10)
!        go to 10
!     endif
!  endif

  if(params%nmode.eq.8) then
! We're in FT8 mode
     
     if(ncontest.eq.6) then
! Fox mode: initialize and open houndcallers.txt     
        inquire(file=trim(temp_dir)//'/houndcallers.txt',exist=ex)
        if(.not.ex) then
           c2fox='            '
           g2fox='    '
           nsnrfox=-99
           nfreqfox=-99
           n30z=0
           my_ft8%nwrap=0
           my_ft8%fox_initialized=.true.
           nfox=0
        endif
        open(19,file=trim(temp_dir)//'/houndcallers.txt',status='unknown')
     endif

     call timer('decft8  ',0)
     ltry_a8=.true.
     newdat=params%newdat
     if(params%emedelay.ne.0.0) then
        id2(1:156000)=id2(24001:180000)  ! Drop the first 2 seconds of data
        id2(156001:180000)=0
     endif
     call my_ft8%decode(ft8_decoded_callback,id2,params%nQSOProgress,params%nfqso,    &
          params%nftx,newdat,params%nutc,params%nfa,params%nfb,              &
          params%nzhsym,params%ndepth,params%emedelay,ncontest,              &
          logical(params%nagain),logical(params%lft8apon),                   &
          ltry_a8,logical(params%lapcqonly),params%napwid,mycall,hiscall,hisgrid,            &
          params%ndiskdat)
     call timer('decft8  ',1)
     if(nfox.gt.0) then
        n30min=minval(n30fox(1:nfox))
        n30max=maxval(n30fox(1:nfox))
     endif
     j=0

     if(ncontest.eq.6) then
! Fox mode: save decoded Hound calls for possible selection by FoxOp
        rewind 19
        if(nfox.eq.0) then
           endfile 19
           rewind 19
        else
           do i=1,nfox
              n=n30fox(i)
              if(n30max-n30fox(i).le.4) then
                 j=j+1
                 c2fox(j)=c2fox(i)
                 g2fox(j)=g2fox(i)
                 nsnrfox(j)=nsnrfox(i)
                 nfreqfox(j)=nfreqfox(i)
                 n30fox(j)=n
                 m=n30max-n
                 if(len(trim(g2fox(j))).eq.4) then
                    call azdist(mygrid,g2fox(j)//'  ',0.d0,nAz,nEl,nDmiles, &
                         nDkm,nHotAz,nHotABetter)
                 else
                    nDkm=9999
                 endif
                 write(19,1004) c2fox(j),g2fox(j),nsnrfox(j),nfreqfox(j),nDkm,m
1004             format(a12,1x,a4,i5,i6,i7,i3)
              endif
           enddo
           nfox=j
           flush(19)
        endif
     endif
     go to 800
  endif

  if(params%nmode.eq.5) then
     call timer('decft4  ',0)
     call my_ft4%decode(ft4_decoded_callback,id2,params%nQSOProgress,params%nfqso,    &
          params%nfa,params%nfb,params%ndepth,                               &
          logical(params%lapcqonly),ncontest,mycall,hiscall)
     call timer('decft4  ',1)
     go to 800
  endif

  if(params%nmode.eq.66) then        !NB: JT65 = 65, Q65 = 66.
! We're in Q65 mode
     open(17,file=trim(temp_dir)//'/red.dat',status='unknown')
     open(14,file=trim(temp_dir)//'/avemsg.txt',status='unknown')
     call timer('dec_q65 ',0)
     nqd=1
     call my_q65%decode(q65_decoded_callback,id2,nqd,params%nutc,params%ntr,      &
          params%nsubmode,params%nfqso,params%ntol,params%ndepth,        &
          params%nfa,params%nfb,logical(params%nclearave),               &
          single_decode,logical(params%nagain),params%max_drift,         &
          logical(params%newdat),params%emedelay,mycall,hiscall,hisgrid, &
          params%nQSOProgress,ncontest,logical(params%lapcqonly),navg0,nqf)
     params%nclearave=.false.

     if(.not.params%nagain) then
! Go through identified candidates again, treating each as if it had been
! double-clicked on the waterfall.
        do k=1,20
           if(nqf(k).eq.0) exit
           if(params%nagain .and. abs(nqf(k)-params%nfqso).gt.params%ntol) cycle
           nqd=1
           navg0=0
           ntol=5
           call my_q65%decode(q65_decoded_callback,id2,nqd,params%nutc,params%ntr,    &
                params%nsubmode,nqf(k),ntol,params%ndepth,                   &
                params%nfa,params%nfb,logical(params%nclearave),             &
                .true.,.true.,params%max_drift,                              &
                .false.,params%emedelay,mycall,hiscall,hisgrid,              &
                params%nQSOProgress,ncontest,logical(params%lapcqonly),      &
                navg0,nqf)
        enddo
     endif

     call timer('dec_q65 ',1)
     close(17)
     go to 800
  endif

  if(params%nmode.eq.240) then
! We're in FST4 mode
     ndepth=iand(params%ndepth,3)
     iwspr=0
     lprinthash22=.false.
     params%nsubmode=0
     call timer('dec_fst4',0)
     call my_fst4%decode(fst4_decoded_callback,id2,params%nutc,                &
          params%nQSOProgress,params%nfa,params%nfb,                  &
          params%nfqso,ndepth,params%ntr,params%nexp_decode,          &
          params%ntol,params%emedelay,logical(params%nagain),         &
          logical(params%lapcqonly),mycall,hiscall,iwspr,lprinthash22)
     call timer('dec_fst4',1)
     go to 800
  endif

    if(params%nmode.eq.241 .or. params%nmode.eq.242) then
! We're in FST4W mode
     ndepth=iand(params%ndepth,3)
     iwspr=1
     lprinthash22=.false.
     if(params%nmode.eq.242) lprinthash22=.true. 
     call timer('dec_fst4',0)
     call my_fst4%decode(fst4_decoded_callback,id2,params%nutc,                &
          params%nQSOProgress,params%nfa,params%nfb,                  &
          params%nfqso,ndepth,params%ntr,params%nexp_decode,          &
          params%ntol,params%emedelay,logical(params%nagain),         &
          logical(params%lapcqonly),mycall,hiscall,iwspr,lprinthash22)
     call timer('dec_fst4',1)
     go to 800
  endif

! Zap data at start that might come from T/R switching transient?
  nadd=100
  k=0
  bad0=.false.
  do i=1,240
     sq=0.
     do n=1,nadd
        k=k+1
        sq=sq + float(id2(k))**2
     enddo
     rms=sqrt(sq/nadd)
     if(rms.gt.10000.0) then
        bad0=.true.
        kbad=k
        rmsbad=rms
     endif
  enddo
  if(bad0) then
     nz=min(NTMAX*12000,kbad+100)
!     id2(1:nz)=0                ! temporarily disabled as it can breaak the JT9 decoder, maybe others
  endif
  
  if(params%nmode.eq.4 .or. params%nmode.eq.65) open(14,file=trim(temp_dir)// &
       '/avemsg.txt',status='unknown')

  if(params%nmode.eq.4) then
     jz=52*nfsample
     if(params%newdat) then
        if(nfsample.eq.12000) call wav11(id2,jz,dd)
        if(nfsample.eq.11025) dd(1:jz)=id2(1:jz)
     else
        jz=52*11025
     endif
     call my_jt4%decode(jt4_decoded_callback,dd,jz,params%nutc,params%nfqso,         &
          params%ntol,params%emedelay,params%dttol,logical(params%nagain),  &
          params%ndepth,logical(params%nclearave),params%minsync,           &
          params%minw,params%nsubmode,mycall,hiscall,         &
          hisgrid,params%nlist,params%listutc,jt4_average_callback)
     go to 800
  endif

  npts65=52*12000
  if(baddata(id2,npts65)) then
     nsynced=0
     ndecoded=0
     go to 800
  endif
 
  ntol65=params%ntol              !### is this OK? ###
  newdat65=params%newdat
  newdat9=params%newdat

!$call omp_set_dynamic(.true.)
!$omp parallel sections num_threads(2) copyin(/timer_private/) shared(ndecoded) if(.true.) !iif() needed on Mac

!$omp section
  if(params%nmode.eq.65) then
! We're in JT65 mode

     if(newdat65) dd(1:npts65)=id2(1:npts65)
     nf1=params%nfa
     nf2=params%nfb
     call timer('jt65a   ',0)
     call my_jt65%decode(jt65_decoded_callback,dd,npts65,newdat65,params%nutc,      &
          nf1,nf2,params%nfqso,ntol65,params%nsubmode,params%minsync,      &
          logical(params%nagain),params%n2pass,logical(params%nrobust),    &
          ntrials,params%naggressive,params%ndepth,params%emedelay,        &
          logical(params%nclearave),mycall,hiscall,          &
          hisgrid,params%nexp_decode,params%nQSOProgress,           &
          logical(params%ljt65apon))
     call timer('jt65a   ',1)

  else if(params%nmode.eq.9 .or. (params%nmode.eq.(65+9) .and. params%ntxmode.eq.9)) then
! We're in JT9 mode, or should do JT9 first
     call timer('decjt9  ',0)
     call my_jt9%decode(jt9_decoded_callback,ss,id2,params%nfqso,       &
          newdat9,params%npts8,params%nfa,params%nfsplit,params%nfb,       &
          params%ntol,params%nzhsym,logical(params%nagain),params%ndepth,  &
          params%nmode,params%nsubmode,params%nexp_decode)
     call timer('decjt9  ',1)
  endif

!$omp section
  if(params%nmode.eq.(65+9)) then       !Do the other mode (we're in dual mode)
     if (params%ntxmode.eq.9) then
        if(newdat65) dd(1:npts65)=id2(1:npts65)
        nf1=params%nfa
        nf2=params%nfb
        call timer('jt65a   ',0)
        call my_jt65%decode(jt65_decoded_callback,dd,npts65,newdat65,params%nutc,   &
             nf1,nf2,params%nfqso,ntol65,params%nsubmode,params%minsync,   &
             logical(params%nagain),params%n2pass,logical(params%nrobust), &
             ntrials,params%naggressive,params%ndepth,params%emedelay,     &
             logical(params%nclearave),mycall,hiscall,       &
             hisgrid,params%nexp_decode,params%nQSOProgress,        &
             logical(params%ljt65apon))
        call timer('jt65a   ',1)
     else
        call timer('decjt9  ',0)
        call my_jt9%decode(jt9_decoded_callback,ss,id2,params%nfqso,                &
             newdat9,params%npts8,params%nfa,params%nfsplit,params%nfb,    &
             params%ntol,params%nzhsym,logical(params%nagain),             &
             params%ndepth,params%nmode,params%nsubmode,params%nexp_decode)
        call timer('decjt9  ',1)
     end if
  endif

!$omp end parallel sections

! JT65 is not yet producing info for nsynced, ndecoded.
800 ndecoded = my_jt4%decoded + my_jt65%decoded + my_jt9%decoded +       &
         my_ft8%decoded + my_ft4%decoded + my_fst4%decoded +             &
         my_q65%decoded
  if(params%nmode.eq.8 .and. params%nzhsym.eq.41) ndec41=ndecoded
  if(params%nmode.eq.8 .and. params%nzhsym.eq.47) ndec47=ndecoded
  if(params%nmode.eq.8 .and. params%nzhsym.eq.50) then
     ndecoded=ndec41+ndec47+ndecoded
  endif
  call wsjtx_decode_stats(params%nzhsym,nsynced,ndecoded,navg0)
  if(params%nmode.ne.8 .or. params%nzhsym.eq.50 .or.                     &
       .not.params%ndiskdat) then

     write(*,1010) nsynced,ndecoded,navg0
1010 format('<DecodeFinished>',2i4,i9)
     call flush(6)
  endif
  close(13)
  if(ncontest.eq.6) close(19)
  if(params%nmode.eq.4 .or. params%nmode.eq.65 .or. params%nmode.eq.66) close(14)
  return
end subroutine multimode_decoder
