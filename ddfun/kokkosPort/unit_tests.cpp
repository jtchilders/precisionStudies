#include <string>
#include "Kokkos_Core.hpp"
#include "test_ddabs.h"
#include "test_ddacosh.h"
#include "test_ddadd.h"
#include "test_ddasinh.h"
#include "test_ddatanh.h"
#include "test_dddiv.h"
#include "test_dddivd.h"
#include "test_ddexp.h"
#include "test_ddlog.h"
#include "test_ddmul.h"
#include "test_ddmuld.h"
// #include "test_ddmuldd.h"
// #include "test_ddneg.h"
// #include "test_ddnpwr.h"
// #include "test_ddpolyr.h"
// #include "test_ddpower.h"
// #include "test_ddsqrt.h"
#include "test_ddsub.h"

int main(int argc, char* argv[]) {

   Kokkos::ScopeGuard scope_guard;

   // user passes path to input files as first parameter on command line
   // extract the path
   std::string path = std::string(argv[1]);

   unittest_ddabs(path + "/ddabs_test_cases.bin");
   unittest_ddacosh(path + "/ddacosh_test_cases.bin");
   unittest_ddadd(path + "/ddadd_test_cases.bin");
   unittest_ddasinh(path + "/ddasinh_test_cases.bin");
   unittest_ddatanh(path + "/ddatanh_test_cases.bin");
   unittest_dddiv(path + "/dddiv_test_cases.bin");
   unittest_dddivd(path + "/dddivd_test_cases.bin");
   unittest_ddexp(path + "/ddexp_test_cases.bin");
   unittest_ddlog(path + "/ddlog_test_cases.bin");
   unittest_ddmul(path + "/ddmul_test_cases.bin");
   unittest_ddmuld(path + "/ddmuld_test_cases.bin");
   // unittest_ddmuldd(path + "/ddmuldd_test_cases.bin");
   // unittest_ddneg(path + "/ddneg_test_cases.bin");
   // unittest_ddnpwr(path + "/ddnpwr_test_cases.bin");
   // unittest_ddpolyr(path + "/ddpolyr_test_cases.bin");
   // unittest_ddpower(path + "/ddpower_test_cases.bin");
   // unittest_ddsqrt(path + "/ddsqrt_test_cases.bin");
   unittest_ddsub(path + "/ddsub_test_cases.bin");

   return 0;
}
