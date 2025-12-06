#pragma once
#include <stdexcept>

class InvalidMatrixSize : public std::domain_error {
public:
	InvalidMatrixSize(const std::string& what_arg) : std::domain_error(what_arg) 
    {}
};