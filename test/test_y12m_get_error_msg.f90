! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:gpt-5
!
! test_y12m_get_error_msg.f90 - Tests for y12m_get_error_msg.
!
program test_y12m_get_error_msg
   use y12m, only: y12m_get_error_msg
   implicit none

   integer :: nfail
   integer :: length
   character(len=256) :: message
   character(len=20) :: short_message

   nfail = 0

   call y12m_get_error_msg(0, message, length)
   call check_equal('ifail=0 message', trim(message), 'Success: no error detected.', nfail)
   call check_int('ifail=0 length', length, len('Success: no error detected.'), nfail)

   call y12m_get_error_msg(1, message, length)
   call check_equal('ifail=1 message', trim(message), &
      'The coefficient matrix A is not factorized. Call Y12MC before Y12MD, or ensure IFLAG(1) is 0 or greater.', &
      nfail)
   call check_int('ifail=1 length', length, len_trim(message), nfail)

   call y12m_get_error_msg(25, message, length)
   call check_equal('ifail=25 message', trim(message), &
      'Row index error. Found an element with a row index outside the range 1 to N.', nfail)
   call check_int('ifail=25 length', length, len_trim(message), nfail)

   call y12m_get_error_msg(999, message, length)
   call check_equal('ifail=999 message', trim(message), 'Unknown error code.', nfail)
   call check_int('ifail=999 length', length, len('Unknown error code.'), nfail)

   call y12m_get_error_msg(1, short_message, length)
   call check_equal('truncated output', short_message, 'The coefficient matr', nfail)
   call check_int('truncated length', length, &
      len_trim('The coefficient matrix A is not factorized. Call Y12MC before Y12MD, or ensure IFLAG(1) is 0 or greater.'), &
      nfail)

   call y12m_get_error_msg(2, message)
   call check_equal('without optional length', trim(message), &
      'The coefficient matrix A is not ordered. Ensure Y12MB is called successfully before calling Y12MC.', nfail)

   if (nfail /= 0) then
      write(*,'(i0,a)') nfail, ' test(s) FAILED'
      stop 1
   end if

   write(*,'(a)') 'All test_y12m_get_error_msg tests PASSED'

contains

   subroutine check_equal(label, got, expected, nfail)
      character(len=*), intent(in) :: label, got, expected
      integer, intent(inout) :: nfail

      if (got /= expected) then
         write(*,'(a)') 'FAIL '//trim(label)
         write(*,'(a)') '  got:      ['//trim(got)//']'
         write(*,'(a)') '  expected: ['//trim(expected)//']'
         nfail = nfail + 1
      else
         write(*,'(a)') 'PASS '//trim(label)
      end if
   end subroutine check_equal

   subroutine check_int(label, got, expected, nfail)
      character(len=*), intent(in) :: label
      integer, intent(in) :: got, expected
      integer, intent(inout) :: nfail

      if (got /= expected) then
         write(*,'(a,i0,a,i0)') 'FAIL '//trim(label)//' got=', got, ' expected=', expected
         nfail = nfail + 1
      else
         write(*,'(a,i0)') 'PASS '//trim(label)//' value=', got
      end if
   end subroutine check_int

end program test_y12m_get_error_msg
