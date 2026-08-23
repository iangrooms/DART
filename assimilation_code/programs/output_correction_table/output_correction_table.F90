program output_correction_table
  use types_mod,       only : r8
  use utilities_mod,   only : initialize_utilities
  use assim_tools_mod, only : get_correction_from_table, assim_tools_init
  implicit none

  integer, parameter :: num_correl = 1000
  integer, parameter :: ens_sizes(3) = [20, 40, 80]
  integer :: i, e
  integer :: ens_size, iunit, io_status

  real(r8) :: scorrel, mean_factor, expected_true_correl, sec

  call initialize_utilities()
  call assim_tools_init()

  ! Open ASCII output file
  open(newunit=iunit, file='correction_results.txt', status='replace', &
       action='write', iostat=io_status)
  if (io_status /= 0) then
     print *, 'Error opening output file.'
     stop 1
  end if

  ! Write header for straightforward MATLAB parsing
  write(iunit, '(A)') 'ens_size scorrel mean_factor expected_true_correl sec'

  ! Loop over ensemble sizes and sample correlations
  do e = 1, size(ens_sizes)
     ens_size = ens_sizes(e)

     do i = 1, num_correl
        scorrel = real(i, r8) / real(num_correl, r8)

        call get_correction_from_table(scorrel, mean_factor, &
                                       expected_true_correl, ens_size)

        if (scorrel /= 0.0_r8) then
           sec = (expected_true_correl / scorrel) * mean_factor
        else
           sec = 0.0_r8
        end if

        write(iunit, '(I6, 1X, ES22.14, 1X, ES22.14, 1X, ES22.14, 1X, ES22.14)') &
             ens_size, scorrel, mean_factor, expected_true_correl, sec
     end do
  end do

  close(iunit)
  print *, 'Data successfully written to correction_results.txt'

end program output_correction_table