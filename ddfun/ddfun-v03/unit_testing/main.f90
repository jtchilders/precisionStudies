program main
   use test_ddabs
   use test_ddacosh
   use test_ddadd
   use test_ddasinh
   use test_ddatanh
   use test_dddiv
   use test_dddivd
   use test_ddexp
   use test_ddlog
   use test_ddmul
   use test_ddmuld
   use test_ddmuldd
   use test_ddneg
   use test_ddnpwr
   use test_ddpower
   use test_ddsqrt
   use test_ddsub
   implicit none

   ! variables for input file path
   character(len=100) :: input_files_path

   ! read path to input files from command line
   call get_command_argument(1, input_files_path)

   write(*,*) "Starting unit tests for DDFUN library..."

   ! Call unit tests with respective input files
   call unittest_ddabs(trim(input_files_path) // "/ddabs_test_cases.bin")
   call unittest_ddacosh(trim(input_files_path) // "/ddacosh_test_cases.bin")
   call unittest_ddadd(trim(input_files_path) // "/ddadd_test_cases.bin")
   call unittest_ddasinh(trim(input_files_path) // "/ddasinh_test_cases.bin")
   call unittest_ddatanh(trim(input_files_path) // "/ddatanh_test_cases.bin")
   call unittest_dddiv(trim(input_files_path) // "/dddiv_test_cases.bin")
   call unittest_dddivd(trim(input_files_path) // "/dddivd_test_cases.bin")
   call unittest_ddexp(trim(input_files_path) // "/ddexp_test_cases.bin")
   call unittest_ddlog(trim(input_files_path) // "/ddlog_test_cases.bin")
   call unittest_ddmul(trim(input_files_path) // "/ddmul_test_cases.bin")
   call unittest_ddmuld(trim(input_files_path) // "/ddmuld_test_cases.bin")
   call unittest_ddmuldd(trim(input_files_path) // "/ddmuldd_test_cases.bin")
   call unittest_ddneg(trim(input_files_path) // "/ddneg_test_cases.bin")
   call unittest_ddnpwr(trim(input_files_path) // "/ddnpwr_test_cases.bin")
   call unittest_ddpower(trim(input_files_path) // "/ddpower_test_cases.bin")
   call unittest_ddsqrt(trim(input_files_path) // "/ddsqrt_test_cases.bin")
   call unittest_ddsub(trim(input_files_path) // "/ddsub_test_cases.bin")

   write(*,*) "All tests completed."

end program main
