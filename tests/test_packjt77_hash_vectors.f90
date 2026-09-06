program test_packjt77_hash_vectors
  use packjt77, only: ihashcall
  implicit none

  integer :: ntests

  ntests=0
  call expect_hash_vector('K1ABC',       712, 2851, 2920267)
  call expect_hash_vector('W9XYZ',       972, 3889, 3982604)
  call expect_hash_vector('PJ4/K1ABC',   346, 1387, 1420834)
  call expect_hash_vector('KH1/KH7Z',    201,  806,  825805)
  call expect_hash_vector('K1JT/P',      959, 3839, 3931158)
  call expect_hash_vector('W3CCX',        665, 2662, 2726483)
  call expect_hash_vector('3DA0ABC',       48,  192,  196695)
  call expect_hash_vector('3XABC',        563, 2254, 2308721)
  call expect_hash_vector('///////////',  671, 2685, 2749801)
  call expect_hash_vector('ZZZZZZZZZZZ',  902, 3609, 3695718)
  call expect_hash_vector('99999999999',  762, 3050, 3123740)
  call expect_hash_vector('A23456789Z/',  122,  488,  499902)
  call expect_hash_vector(' ',              0,    0,       0)

  write(*,'(a,i0)') 'packjt77 hash vector tests passed: ', ntests

contains

  subroutine expect_hash_vector(callsign,want10,want12,want22)
    character(len=*), intent(in) :: callsign
    integer, intent(in) :: want10,want12,want22
    character(len=13) :: c13
    integer :: got10,got12,got22

    c13=' '
    c13=callsign
    got10=ihashcall(c13,10)
    got12=ihashcall(c13,12)
    got22=ihashcall(c13,22)

    if(got10.ne.want10 .or. got12.ne.want12 .or. got22.ne.want22) then
       write(*,'(a,a,a,3(1x,i0),a,3(1x,i0))') &
            'Hash vector failure for "',trim(callsign),'"; got', &
            got10,got12,got22,' wanted',want10,want12,want22
       error stop 1
    endif

    if(got10.lt.0 .or. got10.gt.1023 .or. got12.lt.0 .or. got12.gt.4095 .or. &
       got22.lt.0 .or. got22.gt.4194303) then
       write(*,'(a,a,a,3(1x,i0))') &
            'Hash range failure for "',trim(callsign),'"; got',got10,got12,got22
       error stop 1
    endif

    ntests=ntests+1
  end subroutine expect_hash_vector

end program test_packjt77_hash_vectors
