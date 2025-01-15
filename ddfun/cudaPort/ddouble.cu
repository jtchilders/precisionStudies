#include "ddouble.h"


// Assume ddmul, ddsub, ddadd, dddiv, etc. are implemented as operators on ddouble.
// Assume that ddouble::to_double() returns a double approximation and that
// we have a working constructor ddouble(double) and ddouble(int).
std::string ddouble::str() const {
    if (hi == 0.0 && lo == 0.0) {
        return "0.0";
    }

    // Determine sign and make absolute value
    bool negative = hi < 0 || (hi == 0 && lo < 0);
    ddouble abs_val = negative ? ddouble(-hi, -lo) : *this;

    // Scale abs_val to [1, 10) using npwr and calculate exponent
    const ddouble ten(10.0);
    int exponent = 0;
    if (abs_val.hi >= 10.0) {
        exponent = static_cast<int>(std::log10(abs_val.hi));
        abs_val = abs_val / ten.npwr(exponent);
    } else if (abs_val.hi < 1.0) {
        exponent = static_cast<int>(std::log10(abs_val.hi)) - 1;
        abs_val = abs_val / ten.npwr(exponent);
    }

    // Adjust exponent if necessary
    while (abs_val.hi >= 10.0) {
        abs_val = abs_val / ten;
        ++exponent;
    }
    while (abs_val.hi < 1.0) {
        abs_val = abs_val * ten;
        --exponent;
    }

    // Extract integer part of mantissa
    int first_digit = static_cast<int>(abs_val.hi);
    ddouble fraction = abs_val - ddouble(first_digit);

    // Build string representation
    std::ostringstream oss;
    if (negative) {
        oss << "-";
    }
    oss << first_digit << ".";

    // Extract fractional digits using chunks
    const int DIGITS = 31; // Total digits after the decimal point
    for (int i = 0; i < DIGITS; ++i) {
        fraction = fraction * ten;
        int digit = static_cast<int>(fraction.hi);
        oss << digit;
        fraction = fraction - ddouble(digit);
    }

    // Round last digit if necessary
    if (fraction.hi >= 0.5) {
        std::string result = oss.str();
        for (int i = result.size() - 1; i >= 0; --i) {
            if (result[i] == '.') continue;
            if (result[i] == '9') {
                result[i] = '0';
            } else {
                result[i] += 1;
                break;
            }
        }
        oss.str("");
        oss << result;
    }

    // Append exponent
    oss << "e" << (exponent >= 0 ? "+" : "") << exponent;
    return oss.str();
}


std::string ddouble::strs() const{
    std::stringstream ret;
    // convert hi to string
    ret << "(" << std::setprecision(16) << hi << ",";
    // convert lo to string
    ret << lo << ")";
    return ret.str();
}
