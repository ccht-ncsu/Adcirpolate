!> @author Ali Samii - Institute for Computational Engineering and Sciences (ICES), UT Austin - 2018
!! @brief This module handles the errors and messages that we send to the user.
!! This module handles the errors and messages that we send to the user.
module errors_and_msgs

   use MPI

   interface throw_fatal_error
      procedure throw_error_info_and_stop
      procedure throw_error_and_stop
   end interface

   interface show_message
      procedure just_show_message
   end interface

contains

   subroutine check_error(line_number, file_name, rc)
      integer, intent(in)              :: line_number, rc
      character(len=*), intent(in)     :: file_name
      integer                          :: localPet, ierr

      call MPI_Comm_rank(MPI_COMM_WORLD, localPet, ierr)
      if (rc .NE. 0) then
         write (*, "(A, I4, A, I6, A, A, A, I4)") "Processor ", localPet, &
            " exited line: ", line_number, " in ", file_name, " with error: ", rc
      end if
   end subroutine

   subroutine throw_error_info_and_stop(line_number, file_name, error_message)
      integer, intent(in)               :: line_number
      character(len=*), intent(in)      :: file_name, error_message
      call MPI_Comm_rank(MPI_COMM_WORLD, localPet, ierr)
      write (*, "(A, I4, A, I6, A, A, A, A)") "Processor ", localPet, &
         " exited line: ", line_number, " in: ", file_name, &
         " with error: ", error_message
      stop
   end subroutine

   subroutine throw_error_and_stop(error_message)
      character(len=*), intent(in)      :: error_message
      write (*, "(A, I4, A, I6, A, A, A, A)") error_message
      stop
   end subroutine

   subroutine just_show_message(in_message)
      character(len=*), intent(in)      :: in_message
      write (*, "(A)") in_message
   end subroutine

end module
