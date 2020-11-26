#pragma once
#include <cmath>

class Complex
{
public:
	double re, im;

	Complex(float _re = 0, float _im = 0) :re(_re), im(_im)
	{
	}

	float abs()
	{
		return sqrt(re*re + im*im);
	}

	Complex operator+(const Complex &c)
	{
		return Complex(re + c.re, im + c.im);
	}

	Complex operator*(const Complex &c)
	{
		return Complex(re*c.re - im*c.im, re*c.im + im*c.re);
	}

	Complex operator*(float a)
	{
		return Complex(re*a, im*a);
	}
	Complex operator*(double a)
	{
		return Complex(re*a, im*a);
	}

};