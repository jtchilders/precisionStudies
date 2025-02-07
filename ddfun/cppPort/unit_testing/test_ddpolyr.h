#ifndef TEST_DDPOLYR_H
#define TEST_DDPOLYR_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

void unittest_ddpolyr(const std::string &filename) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    std::cout << "Running tests for ddpolyr from: " << filename << std::endl;

    int degree, placeholder;
    int total_tests = 0;
    int passed_tests = 0;
    const double tolerance = 30;

    while (infile.read(reinterpret_cast<char *>(&degree), sizeof(int))) {
        
        infile.read(reinterpret_cast<char *>(&placeholder), sizeof(int));

        std::vector<ddouble> coefficients(degree + 1);
        for (int i = 0; i <= degree; ++i) {
            infile.read(reinterpret_cast<char *>(&coefficients[i].hi), sizeof(double));
            infile.read(reinterpret_cast<char *>(&coefficients[i].lo), sizeof(double));
        }

        double hi_x0, lo_x0, expected_hi, expected_lo;
        infile.read(reinterpret_cast<char *>(&hi_x0), sizeof(double));
        infile.read(reinterpret_cast<char *>(&lo_x0), sizeof(double));
        infile.read(reinterpret_cast<char *>(&expected_hi), sizeof(double));
        infile.read(reinterpret_cast<char *>(&expected_lo), sizeof(double));

        ddouble x0{hi_x0, lo_x0};
        ddouble expected{expected_hi, expected_lo};

        total_tests++;

        // Perform ddpolyr operation
        ddouble result = ddpolyr(degree, coefficients, x0);

        // Use scale difference for comparison
        int scale_diff = calculate_scale_difference(result, expected);
        bool test_passed = (scale_diff >= tolerance or scale_diff == 0);

        if (test_passed) {
            passed_tests++;
        } else {
            std::cout << "Test Failed: " << std::setprecision(16)
                      << "degree: " << degree
                      << " initial x0: [" << x0.hi << ", " << x0.lo << "] "
                      << "result: [" << result.hi << ", " << result.lo << "] "
                      << "expected: [" << expected.hi << ", " << expected.lo << "] "
                      << "error: [" << std::abs(result.hi - expected.hi) << ", " << std::abs(result.lo - expected.lo) << "]"
                   << "scale difference: [" << scale_diff << "]\n";
        }
    }

    infile.close();

    std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDPOLYR_H