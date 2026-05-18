module decode_callbacks

  use jt4_decode, only: jt4_decoder
  use jt65_decode, only: jt65_decoder
  use jt9_decode, only: jt9_decoder
  use ft8_decode, only: ft8_decoder, MAXFOX, c2fox, g2fox, nsnrfox, &
       nfreqfox, n30fox, n30z, nfox
  use ft4_decode, only: ft4_decoder
  use fst4_decode, only: fst4_decoder
  use q65_decode, only: q65_decoder, cq0, msg0, nsnr0, xdt0, nfreq0

  implicit none

  private

  public :: counting_jt4_decoder, counting_jt65_decoder
  public :: counting_jt9_decoder, counting_ft8_decoder
  public :: counting_ft4_decoder, counting_fst4_decoder
  public :: counting_q65_decoder, map65_q65_decoder
  public :: jt4_decoded_callback, jt4_average_callback
  public :: jt65_decoded_callback, jt9_decoded_callback
  public :: ft8_decoded_callback, ft4_decoded_callback
  public :: fst4_decoded_callback, q65_decoded_callback
  public :: map65_q65_decoded_callback

  type, extends(jt4_decoder) :: counting_jt4_decoder
     integer :: decoded = 0
     integer :: nutc = 0
  end type counting_jt4_decoder

  type, extends(jt65_decoder) :: counting_jt65_decoder
     integer :: decoded = 0
     integer :: nutc = 0
     logical :: bVHF = .false.
  end type counting_jt65_decoder

  type, extends(jt9_decoder) :: counting_jt9_decoder
     integer :: decoded = 0
     integer :: nutc = 0
  end type counting_jt9_decoder

  type, extends(ft8_decoder) :: counting_ft8_decoder
     integer :: decoded = 0
     integer :: nutc = 0
     integer :: ncontest = 0
     integer :: nwrap = 0
     logical :: fox_initialized = .false.
     character(len=12) :: mycall = '            '
  end type counting_ft8_decoder

  type, extends(ft4_decoder) :: counting_ft4_decoder
     integer :: decoded = 0
     integer :: nutc = 0
  end type counting_ft4_decoder

  type, extends(fst4_decoder) :: counting_fst4_decoder
     integer :: decoded = 0
  end type counting_fst4_decoder

  type, extends(q65_decoder) :: counting_q65_decoder
     integer :: decoded = 0
  end type counting_q65_decoder

  type, extends(q65_decoder) :: map65_q65_decoder
     integer :: decoded = 0
  end type map65_q65_decoder

  external :: wsjtx_decoded

contains

  subroutine jt4_decoded_callback(this,snr,dt,freq,have_sync,sync,is_deep, &
       decoded0,qual,ich,is_average,ave)
    implicit none
    class(jt4_decoder), intent(inout) :: this
    integer, intent(in) :: snr
    real, intent(in) :: dt
    integer, intent(in) :: freq
    logical, intent(in) :: have_sync
    logical, intent(in) :: is_deep
    character(len=1), intent(in) :: sync
    character(len=22), intent(in) :: decoded0
    real, intent(in) :: qual
    integer, intent(in) :: ich
    logical, intent(in) :: is_average
    integer, intent(in) :: ave

    character*22 decoded
    character*3 cflags
    integer :: nutc

    nutc = 0
    select type(ctx => this)
    type is (counting_jt4_decoder)
       nutc = ctx%nutc
    end select

    if(ich.eq.-99) stop
    if (have_sync) then
       decoded=decoded0
       cflags='   '
       if(decoded.ne.'                      ') then
          cflags='f  '
          if(is_deep) then
             cflags='d  '
             write(cflags(2:2),'(i1)') min(int(qual),9)
             if(qual.ge.10.0) cflags(2:2)='*'
             if(qual.lt.3.0) decoded(22:22)='?'
          endif
          if(is_average) then
             write(cflags(3:3),'(i1)') min(ave,9)
             if(ave.ge.10) cflags(3:3)='*'
             if(cflags(1:1).eq.'f') cflags=cflags(1:1)//cflags(3:3)//' '
          endif
       endif
       write(*,1000) nutc,snr,dt,freq,sync,decoded,cflags
1000   format(i4.4,i4,f5.1,i5,1x,'$',a1,1x,a22,1x,a3)
    else
       write(*,1000) nutc,snr,dt,freq
    end if

    select type(ctx => this)
    type is (counting_jt4_decoder)
       ctx%decoded = ctx%decoded + 1
    end select
  end subroutine jt4_decoded_callback

  subroutine jt4_average_callback(this, used, utc, sync, dt, freq, flip)
    implicit none
    class(jt4_decoder), intent(inout) :: this
    logical, intent(in) :: used
    integer, intent(in) :: utc
    real, intent(in) :: sync
    real, intent(in) :: dt
    integer, intent(in) :: freq
    logical, intent(in) :: flip
    character(len=1) :: cused, csync

    cused = '.'
    csync = '*'
    if (used) cused = '$'
    if (flip) csync = '$'
    write(14,1000) cused,utc,sync,dt,freq,csync
1000 format(a1,i5.4,f6.1,f6.2,i6,1x,a1)
  end subroutine jt4_average_callback

  subroutine jt65_decoded_callback(this,sync,snr,dt,freq,drift,nflip,width, &
       decoded0,ft,qual,nsmo,nsum,minsync)
    implicit none

    class(jt65_decoder), intent(inout) :: this
    real, intent(in) :: sync
    integer, intent(in) :: snr
    real, intent(in) :: dt
    integer, intent(in) :: freq
    integer, intent(in) :: drift
    integer, intent(in) :: nflip
    real, intent(in) :: width
    character(len=22), intent(in) :: decoded0
    integer, intent(in) :: ft
    integer, intent(in) :: qual
    integer, intent(in) :: nsmo
    integer, intent(in) :: nsum
    integer, intent(in) :: minsync

    integer i,n,nap,nutc
    logical is_deep,is_average,bVHF
    character decoded*22,csync*2,cflags*3

    nutc = 0
    bVHF = .false.
    select type(ctx => this)
    type is (counting_jt65_decoder)
       nutc = ctx%nutc
       bVHF = ctx%bVHF
    end select

!$omp critical(decode_results)
    decoded=decoded0
    cflags='   '
    is_deep=ft.eq.2

    if(ft.eq.0 .and. minsync.ge.0 .and. int(sync).lt.minsync) then
       write(*,1010) nutc,snr,dt,freq
    else
       is_average=nsum.ge.2
       if(bVHF .and. ft.gt.0) then
          cflags='f  '
          if(is_deep) then
             cflags='d  '
             write(cflags(2:2),'(i1)') min(qual,9)
             if(qual.ge.10) cflags(2:2)='*'
             if(qual.lt.3) decoded(22:22)='?'
          endif
          if(is_average) then
             write(cflags(3:3),'(i1)') min(nsum,9)
             if(nsum.ge.10) cflags(3:3)='*'
          endif
          nap=ishft(ft,-2)
          if(nap.ne.0) then
             if(nsum.lt.2) write(cflags(1:3),'(a1,i1," ")') 'a',nap
             if(nsum.ge.2) write(cflags(1:3),'(a1,2i1)') 'a',nap,min(nsum,9)
          endif
       endif
       csync='# '
       i=0
       if(bVHF .and. nflip.ne.0 .and. sync.ge.max(0.0,float(minsync))) then
          csync='#*'
          if(nflip.eq.-1) then
             csync='##'
             if(decoded.ne.'                      ') then
                do i=22,1,-1
                   if(decoded(i:i).ne.' ') exit
                enddo
                if(i.gt.18) i=18
                decoded(i+2:i+4)='OOO'
             endif
          endif
       endif
       n=len(trim(decoded))
       if(n.eq.2 .or. n.eq.3) csync='# '
       if(cflags(1:1).eq.'f') then
          cflags(2:2)=cflags(3:3)
          cflags(3:3)=' '
       endif
       write(*,1010) nutc,snr,dt,freq,csync,decoded,cflags
1010   format(i4.4,i4,f5.1,i5,1x,a2,1x,a22,1x,a3)
    endif
    call wsjtx_decoded(nutc,snr,dt,freq,decoded)
    call flush(6)

!$omp end critical(decode_results)
    select type(ctx => this)
    type is (counting_jt65_decoder)
       ctx%decoded = ctx%decoded + 1
    end select
  end subroutine jt65_decoded_callback

  subroutine jt9_decoded_callback(this, sync, snr, dt, freq, drift, decoded)
    implicit none

    class(jt9_decoder), intent(inout) :: this
    real, intent(in) :: sync
    integer, intent(in) :: snr
    real, intent(in) :: dt
    real, intent(in) :: freq
    integer, intent(in) :: drift
    character(len=22), intent(in) :: decoded
    integer :: nutc

    nutc = 0
    select type(ctx => this)
    type is (counting_jt9_decoder)
       nutc = ctx%nutc
    end select

!$omp critical(decode_results)
    write(*,1000) nutc,snr,dt,nint(freq),decoded
1000 format(i4.4,i4,f5.1,i5,1x,'@ ',1x,a22)
    call wsjtx_decoded(nutc,snr,dt,nint(freq),decoded)
    call flush(6)
!$omp end critical(decode_results)
    select type(ctx => this)
    type is (counting_jt9_decoder)
       ctx%decoded = ctx%decoded + 1
    end select
  end subroutine jt9_decoded_callback

  subroutine ft8_decoded_callback(this,sync,snr,dt,freq,decoded,nap,qual)
    implicit none

    class(ft8_decoder), intent(inout) :: this
    real, intent(in) :: sync
    integer, intent(in) :: snr
    real, intent(in) :: dt
    real, intent(in) :: freq
    character(len=37), intent(in) :: decoded
    integer, intent(in) :: nap
    real, intent(in) :: qual
    character c1*12,c2*12,g2*4,w*4
    integer i0,i1,i2,i3,i4,i5,n30,n
    integer ncontest,nutc
    character(len=12) :: mycall
    character*2 annot
    character*37 decoded0
    logical isgrid4,b0,b1,b2

    isgrid4(w)=(len_trim(w).eq.4 .and.                                        &
         ichar(w(1:1)).ge.ichar('A') .and. ichar(w(1:1)).le.ichar('R') .and.  &
         ichar(w(2:2)).ge.ichar('A') .and. ichar(w(2:2)).le.ichar('R') .and.  &
         ichar(w(3:3)).ge.ichar('0') .and. ichar(w(3:3)).le.ichar('9') .and.  &
         ichar(w(4:4)).ge.ichar('0') .and. ichar(w(4:4)).le.ichar('9'))

    nutc = 0
    ncontest = 0
    mycall = '            '
    select type(ctx => this)
    type is (counting_ft8_decoder)
       nutc = ctx%nutc
       ncontest = ctx%ncontest
       mycall = ctx%mycall
       if(.not.ctx%fox_initialized) then
          c2fox='            '
          g2fox='    '
          nsnrfox=-99
          nfreqfox=-99
          n30z=0
          ctx%nwrap=0
          nfox=0
          ctx%fox_initialized=.true.
       endif
    end select

    decoded0=decoded

    annot='  '
    if(nap.ne.0) then
       write(annot,'(a1,i1)') 'a',nap
       if(qual.lt.0.17) decoded0(37:37)='?'
    endif

    i0=1
    if(i0.le.0) write(*,1000) nutc,snr,dt,nint(freq),decoded0(1:22),annot
1000 format(i6.6,i4,f5.1,i5,' ~ ',1x,a22,1x,a2)
    if(i0.gt.0) write(*,1001) nutc,snr,dt,nint(freq),decoded0,annot
1001 format(i6.6,i4,f5.1,i5,' ~ ',1x,a37,1x,a2)

    call wsjtx_decoded(nutc,snr,dt,nint(freq),decoded0)
    if(ncontest.eq.6) then
       i1=index(decoded0,' ')
       i2=i1 + index(decoded0(i1+1:),' ')
       i3=i2 + index(decoded0(i2+1:),' ')
       if(i1.ge.3 .and. i2.ge.7 .and. i3.ge.10) then
          c1=decoded0(1:i1-1)//'            '
          c2=decoded0(i1+1:i2-1)
          g2=decoded0(i2+1:i3-1)
          b0=c1.eq.mycall
          if(c1(1:3).eq.'DE ' .and. index(c2,'/').ge.2) b0=.true.
          if(len(trim(c1)).ne.len(trim(mycall))) then
             i4=index(trim(c1),trim(mycall))
             i5=index(trim(mycall),trim(c1))
             if(i4.ge.1 .or. i5.ge.1) b0=.true.
          endif
          b1=i3-i2.eq.5 .and. isgrid4(g2)
          b2=i3-i2.eq.1
          if(b0 .and. (b1.or.b2) .and. nint(freq).ge.1000) then
             n=nutc
             n30=(3600*(n/10000) + 60*mod((n/100),100) + mod(n,100))/30
             select type(ctx => this)
             type is (counting_ft8_decoder)
                if(n30.lt.n30z) ctx%nwrap=ctx%nwrap+5760
                n30z=n30
                n30=n30+ctx%nwrap
             end select
             if(nfox.lt.MAXFOX) nfox=nfox+1
             c2fox(nfox)=c2
             g2fox(nfox)=g2
             nsnrfox(nfox)=snr
             nfreqfox(nfox)=nint(freq)
             n30fox(nfox)=n30
          endif
       endif
    endif

    call flush(6)

    select type(ctx => this)
    type is (counting_ft8_decoder)
       ctx%decoded = ctx%decoded + 1
    end select
  end subroutine ft8_decoded_callback

  subroutine ft4_decoded_callback(this,sync,snr,dt,freq,decoded,nap,qual)
    implicit none

    class(ft4_decoder), intent(inout) :: this
    real, intent(in) :: sync
    integer, intent(in) :: snr
    real, intent(in) :: dt
    real, intent(in) :: freq
    character(len=37), intent(in) :: decoded
    integer, intent(in) :: nap
    real, intent(in) :: qual
    integer :: nutc
    character*2 annot
    character*37 decoded0

    nutc = 0
    select type(ctx => this)
    type is (counting_ft4_decoder)
       nutc = ctx%nutc
    end select

    decoded0=decoded

    annot='  '
    if(nap.ne.0) then
       write(annot,'(a1,i1)') 'a',nap
       if(qual.lt.0.17) decoded0(37:37)='?'
    endif

    write(*,1001) nutc,snr,dt,nint(freq),decoded0,annot
1001 format(i6.6,i4,f5.1,i5,' + ',1x,a37,1x,a2)
    call wsjtx_decoded(nutc,snr,dt,nint(freq),decoded0)
10  call flush(6)

    select type(ctx => this)
    type is (counting_ft4_decoder)
       ctx%decoded = ctx%decoded + 1
    end select
  end subroutine ft4_decoded_callback

  subroutine fst4_decoded_callback(this,nutc,sync,nsnr,dt,freq,decoded,nap, &
       qual,ntrperiod,fmid,w50)
    implicit none

    class(fst4_decoder), intent(inout) :: this
    integer, intent(in) :: nutc
    real, intent(in) :: sync
    integer, intent(in) :: nsnr
    real, intent(in) :: dt
    real, intent(in) :: freq
    character(len=37), intent(in) :: decoded
    integer, intent(in) :: nap
    real, intent(in) :: qual
    integer, intent(in) :: ntrperiod
    real, intent(in) :: fmid
    real, intent(in) :: w50

    character*2 annot
    character*37 decoded0
    character*70 line

    decoded0=decoded
    annot='  '
    if(nap.ne.0) then
       write(annot,'(a1,i1)') 'a',nap
       if(qual.lt.0.17) decoded0(37:37)='?'
    endif

    if(ntrperiod.lt.60) then
       write(line,1001) nutc,nsnr,dt,nint(freq),decoded0,annot
1001   format(i6.6,i4,f5.1,i5,' ` ',1x,a37,1x,a2)
       call wsjtx_decoded(nutc,nsnr,dt,nint(freq),decoded0)
    else
       write(line,1003) nutc,nsnr,dt,nint(freq),decoded0,annot
1003   format(i4.4,i4,f5.1,i5,' ` ',1x,a37,1x,a2,2f7.3)
       call wsjtx_decoded(nutc,nsnr,dt,nint(freq),decoded0)
    endif

    if(fmid.ne.-999.0) then
       if(w50.lt.0.95) write(line(65:70),'(f6.3)') w50
       if(w50.ge.0.95) write(line(65:70),'(f6.2)') w50
    endif

    write(*,1005) line
1005 format(a70)

    call flush(6)

    select type(ctx => this)
    type is (counting_fst4_decoder)
       ctx%decoded = ctx%decoded + 1
    end select
  end subroutine fst4_decoded_callback

  subroutine q65_decoded_callback(this,nutc,snr1,nsnr,dt,freq,decoded,idec, &
       nused,ntrperiod)
    implicit none

    class(q65_decoder), intent(inout) :: this
    integer, intent(in) :: nutc
    real, intent(in) :: snr1
    integer, intent(in) :: nsnr
    real, intent(in) :: dt
    real, intent(in) :: freq
    character(len=37), intent(in) :: decoded
    integer, intent(in) :: idec
    integer, intent(in) :: nused
    integer, intent(in) :: ntrperiod
    character*3 cflags

    cflags='   '
    if(idec.ge.0) then
       cflags='q  '
       write(cflags(2:2),'(i1)') idec
       if(nused.ge.2) write(cflags(3:3),'(i1)') nused
    endif

    if(ntrperiod.lt.60) then
       write(*,1001) nutc,nsnr,dt,nint(freq),decoded,cflags
1001   format(i6.6,i4,f5.1,i5,' : ',1x,a37,1x,a3)
    else
       write(*,1003) nutc,nsnr,dt,nint(freq),decoded,cflags
1003   format(i4.4,i4,f5.1,i5,' : ',1x,a37,1x,a3)
    endif
    call flush(6)

    select type(ctx => this)
    type is (counting_q65_decoder)
       if(idec.ge.0) ctx%decoded = ctx%decoded + 1
    end select
  end subroutine q65_decoded_callback

  subroutine map65_q65_decoded_callback(this,nutc,snr1,nsnr,dt,freq,decoded, &
       idec,nused,ntrperiod)
    implicit none

    class(q65_decoder), intent(inout) :: this
    integer, intent(in) :: nutc
    real, intent(in) :: snr1
    integer, intent(in) :: nsnr
    real, intent(in) :: dt
    real, intent(in) :: freq
    character(len=37), intent(in) :: decoded
    integer, intent(in) :: idec
    integer, intent(in) :: nused
    integer, intent(in) :: ntrperiod

    if(nutc+snr1+nsnr+dt+freq+idec+nused+ntrperiod.eq.-999) stop
    if(decoded.eq.'-999') stop

    cq0='q  '
    write(cq0(2:2),'(i1)') idec
    if(nused.ge.2) write(cq0(3:3),'(i1)') nused
    nsnr0=nsnr
    xdt0=dt
    nfreq0=nint(freq)
    msg0=decoded

    select type(ctx => this)
    type is (map65_q65_decoder)
       if(idec.ge.0) ctx%decoded = ctx%decoded + 1
    end select
  end subroutine map65_q65_decoded_callback

end module decode_callbacks
