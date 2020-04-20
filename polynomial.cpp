#include "polynomial.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cfloat>
#include <cmath>

// Name: Nghi Vo (Ivy)
// Lab 1: Polynominal
// CSC 2431
// Description: practice dynamic memory with a simple class.
// Given the polynomial and with functions provide by Dr.Arias. Finishing the sum, subtract, multiply and derive.

using std::istream;
using std::ostream;
using std::string;
using std::stringstream;
using std::fixed;
using std::setprecision;
using std::showpos;

Polynomial::Polynomial(size_t degree) : _degree(degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = 0.0;
	}
}
Polynomial::Polynomial(size_t degree, const float* coefficients): _degree(degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = coefficients[i];
	}
}
Polynomial::Polynomial(const Polynomial& polynomial): _degree(polynomial._degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = polynomial._coefficients[i];
	}
}
Polynomial::~Polynomial(){
	// DO THIS FIRST TO PREVENT MEMORY LEAKS!
	if (_coefficients != nullptr)
	    delete [] _coefficients;
}
const Polynomial Polynomial::Sum(const Polynomial& rhs)const{
    // If else to check which degree to use in the return polynomial.
    // Using degree of larger one
    if (this->_degree < rhs._degree){
        Polynomial retVal(rhs._degree);
        // Adding the with the same degree
        for (size_t i = 0; i <= this->_degree; i++){
            retVal._coefficients[i] = this->_coefficients[i] + rhs._coefficients[i];
        }
        // Adding the degree higher (which do not have the same
        for (size_t j = rhs._degree; j <= rhs._degree - this->_degree; j--){
            retVal._coefficients[j] = rhs._coefficients[j];
        }
        return retVal;
    }
    else {
        Polynomial retVal(this->_degree);
        for (size_t i = 0; i <= this->_degree; i++) {
            retVal._coefficients[i] = this->_coefficients[i] + rhs._coefficients[i];
        }
        for (size_t j = this->_degree; j <= this->_degree - rhs._degree; j--) {
            retVal._coefficients[j] = this->_coefficients[j];
        }
        return retVal;
    }
}
const Polynomial Polynomial::Subtract(const Polynomial& rhs)const{
    // If else to check which degree to use in the return polynomial.
    // Using degree of larger one
    if (this->_degree < rhs._degree){
        Polynomial retVal(rhs._degree);
        // Subtract the with the same degree
        for (size_t i = 0; i <= rhs._degree; i++){
            retVal._coefficients[i] =  rhs._coefficients[i] - this->_coefficients[i];
        }
        // Adding the rest the with different degree
        for (size_t j = rhs._degree; j <= rhs._degree - this->_degree; j--){
            retVal._coefficients[j] = rhs._coefficients[j];
        }
        return retVal;
    }
    else {
        Polynomial retVal(this->_degree);
        for (size_t i = 0; i <= this->_degree; i++) {
            retVal._coefficients[i] = this->_coefficients[i] - rhs._coefficients[i];
        }
        for (size_t j = this->_degree; j <= this->_degree - rhs._degree; j--) {
            retVal._coefficients[j] = this->_coefficients[j];
        }
        return retVal;
    }
}
const Polynomial Polynomial::Minus()const{
	Polynomial retVal(this->_degree);
	for (size_t i = 0; i < _degree + 1; i++) {
		retVal._coefficients[i] *= -1;
	}
	return retVal;
}
const Polynomial Polynomial::Multiply(const Polynomial& rhs)const{
    // Using the degree of sum of both polynomial
	Polynomial retVal(this->_degree + rhs._degree);
	for (size_t i = 0; i <= this->_degree; i++){
        for (size_t j = 0; j <= rhs._degree; j++){
            retVal._coefficients[i+j] += this->_coefficients[i] * rhs._coefficients[j];
        }
	}
	return retVal;
}
const Polynomial Polynomial::Divide(const Polynomial& rhs)const{
    return Polynomial(0);

}
const Polynomial Polynomial::Derive()const{
    Polynomial retVal (this->_degree - 1);
    for (size_t i = 1; i <= this->_degree; i++){
        retVal._coefficients[i-1] = this->_coefficients[i] * i;
    }
    return retVal;
}
float Polynomial::Evaluate(float x)const{
	float eVal = 0.0;
	for (size_t i = 0; i <= this->_degree; i++){
	    eVal += this->_coefficients[i] * powf(x, i);
	}
	return eVal;
}
float Polynomial::Integrate(float start, float end)const{
    Polynomial retVal (this->_degree  + 1);
    for (size_t i = 0; i <= this->_degree; i++){
        retVal._coefficients[i+1] = this->_coefficients[i] / (i+1);
    }
    float result = retVal.Evaluate(end) - retVal.Evaluate(start);

	return result;

}
const Polynomial& Polynomial::operator=(const Polynomial& rhs){
	if (&rhs == this){
		return *this;
	}
	if (_degree != rhs._degree){
		if (_coefficients){
			delete[] _coefficients;
		}
		_degree = rhs._degree;
		_coefficients = new float[_degree + 1];
	}
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = rhs._coefficients[i];
	}
	return *this;
}
bool Polynomial::Equals(const Polynomial& rhs)const{
	if (_degree != rhs._degree){
		return false;
	}
	for (size_t i=0; i < _degree; i++){
		if (abs(_coefficients[i] - rhs._coefficients[i]) > 0.0001){
			return false;
		}
	}
	return true;
}
string Polynomial::ToString()const{
	stringstream ss;
	for (size_t i = _degree; i > 0; i--) {
		ss << showpos << fixed << setprecision(2) << _coefficients[i] << "x^" << i << " ";
	}
	ss << showpos << fixed << setprecision(2) << _coefficients[0];
	return ss.str();
}
ostream& Polynomial::Write(ostream& output)const{
	output << _degree << " ";
	for (size_t i = 0; i < _degree + 1; i++) {
		output << _coefficients[i] << " ";
	}
	return output;
}
istream& Polynomial::Read(istream& input){
	size_t degree;
	input >> degree;
	if (input.fail()){
		return input;
	}
	float* coefficients = new float[degree + 1];
	for (size_t i = 0; i < degree + 1; i++) {
		input >> coefficients[i];
		if (input.fail()){
			delete[] coefficients;
			return input;
		}
	}

	if (degree != _degree){
		if (_coefficients){
			delete[] _coefficients;
		}
		_degree = degree;
		_coefficients = coefficients;
	}else{
		for (size_t i = 0; i < _degree + 1; i++) {
			_coefficients[i] = coefficients[i];
		}
		delete[] coefficients;
	}
	return input;
}
