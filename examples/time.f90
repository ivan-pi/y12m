
! Returns the time in millisecons
subroutine time(imsec)
  integer, intent(out) :: imsec
  integer(8) :: count, rate
  call system_clock(count,rate)
  imsec = real(count) / real(rate) * 1.0e-6
end subroutine