#pragma once
#include <locale>
class comma : public std::numpunct<char> {
public:
    comma() : std::numpunct<char>() {}
protected:
    char do_decimal_point() const {
        return ',';
    }
};