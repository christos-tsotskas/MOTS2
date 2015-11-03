/*
 * objective_function_formulae.cpp
 *
 *  Created on: Jan 20, 2014
 *      Author: ctsotskas
 */
#include <cassert>
#include <cstdlib>
#include <cmath>

#include "objective_function_formulae.h"

#define pi 3.14159

extern "C" {
int ffd_ (const double&, const int*, const double&, const double&, const int&, const int&);
//      int ffd_ (const double&, const int*, const char*, const int&, const int&);
int xfoil_ (const double&, const double&, const char*, const int&, const double&, const double&, const int&);
//      int xfoil_ (const double&, const double&, const char*, const int&);
}

ObjFunction2 eval_ZDT1(const Point2 &input)
{

	int i = 0;
	int n = input.size();
	double f1 = 0;
	double g = 0;
	double h = 0;
	ObjFunction2 output(2,0.0);

#ifdef CONTROL_TEST_FUNCTIONS
	if( input.size() != 30){
		std::cout << " ZDT1 input has" << input.size() << "variables, instead of 30" << std::endl;
		exit(-100);
	}
	assert(n == 30); //variable vector length
#endif //CONTROL_TEST_FUNCTIONS
	assert(output.size() == 2);

	f1 = input[0];

	for (i = 1; i < n; i++)
	{
		g += input[i];
	}
	g = 1 + 9 * g / (n-1);
	h = 1 - sqrt(f1 / g);

	output[0]=f1;
	output[1] = g * h;

	return output;
}


ObjFunction2 eval_ZDT2(const Point2 &input)
{
	int i = 0;
	int n = input.size();
	double f1 = 0;
	double g = 0;
	double h = 0;
	ObjFunction2 output(2,0.0);

#ifdef CONTROL_TEST_FUNCTIONS
	if( input.size() != 30){
		std::cout << " ZDT3 input has" << input.size() << "variables, instead of 30" << std::endl;
		exit(-100);
	}
	assert(n == 30); //variable vector length
#endif //CONTROL_TEST_FUNCTIONS

	assert(output.size() == 2);

	f1 = input[0];

	for (i = 1; i < n; i++)
	{
		g += input[i];
	}
	g = 1 + 9 * g / (n-1);
	h = 1 - pow(f1 / g, 2);

	output[0]=f1;
	output[1] = g * h;

	return output;
}


ObjFunction2 eval_ZDT3(const Point2 &input)
{

	int i = 0;
	int n = input.size();
	double f1 = 0.0;
	double g = 0.0;
	double h = 0.0;


	ObjFunction2 output(2,0.0);

#ifdef CONTROL_TEST_FUNCTIONS
	if( input.size() != 30){
		std::cout << " ZDT3 input has" << input.size() << "variables, instead of 30" << std::endl;
		exit(-100);
	}
	assert(n == 30); //variable vector length
#endif //CONTROL_TEST_FUNCTIONS
	assert(output.size() == 2);

	f1 = input[0];

	for (i = 1; i < n; i++)
	{
		g += input[i];
	}
	g = 1 + 9 * g / (n-1);
	h = 1 - std::sqrt(f1 / g) - (f1 / g) * std::sin(10 * pi * f1);

	output[0]=f1;
	output[1] = g * h;

	return output;
}


ObjFunction2 eval_ZDT4(const Point2 &input)
{

	int i = 0;
	int n = input.size();
	double f1 = 0;
	double g = 0;
	double h = 0;
	ObjFunction2 output(2,0.0);

#ifdef CONTROL_TEST_FUNCTIONS
	if( input.size() != 10){
		std::cout << " ZDT4 input has" << input.size() << "variables, instead of 10" << std::endl;
		exit(-100);
	}
	assert(n == 10); //variable vector length
#endif //CONTROL_TEST_FUNCTIONS
	assert(output.size() == 2);

	f1 = input[0];

	for (i = 1; i < n; i++)
	{
		double x = input[i];
		g += x * x - 10 * cos(4 * pi * x);
	}
	g = 1 + 10 * (n - 1) + g;
	h = 1 - sqrt(f1 / g);

	output[0]=f1;
	output[1] = g * h;

	return output;
}


ObjFunction2  eval_ZDT6(const Point2 &input)
{
	int i = 0;
	int n = input.size();
	double f1 = 0;
	double g = 0;
	double h = 0;
	ObjFunction2 output(2,0.0);

#ifdef CONTROL_TEST_FUNCTIONS
	if( input.size() != 10){
		std::cout << " ZDT6 input has" << input.size() << "variables, instead of 10" << std::endl;
		exit(-100);
	}
	assert(n == 10); //variable vector length
#endif //CONTROL_TEST_FUNCTIONS

	assert(output.size() == 2);

	f1 = 1 - exp(-4 * input[0]) * pow(sin(6 * pi * input[0]), 6);

	for (i = 1; i < n; i++)
	{
		g += input[i];
	}
	g = 1 + 9 * pow(g / (n-1), 0.25);
	h = 1 - pow(f1 / g, 2);



	output[0]=f1;
	output[1] = g * h;

	return output;
}


ObjFunction2  eval_DTLZ1(Point2 &input, unsigned int const nObj_size)
{
	ObjFunction2 output(nObj_size,0.0);
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int n = input.size();
	unsigned int k = n - output.size() + 1;

	double g = 0;
	for (i = n - k + 1; i <= n; i++)
	{
		g += pow(input[i-1]-0.5,2) - cos(20 * pi * (input[i-1]-0.5));
	}
	g = 100 * (k + g);

	for (i = 1; i <= output.size(); i++)
	{
		double f = 0.5 * (1 + g);
		for (j = output.size() - i; j >= 1; j--)
		{
			f *= input[j-1];
		}
		if (i > 1)
		{
			f *= 1 - input[(output.size() - i + 1) - 1];
		}

		output[i-1] = f;
	}



	return output;
}


ObjFunction2  eval_DTLZ2(Point2 &input, unsigned int const nObj_size)
{
	ObjFunction2 output(nObj_size,0.0);
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int n = input.size();
	unsigned int k = n - output.size() + 1;

	double g = 0;
	for (i = n - k + 1; i <= n; i++)
	{
		g += pow(input[i-1]-0.5,2);
	}

	for (i = 1; i <= output.size(); i++)
	{
		double f = (1 + g);
		for (j = output.size() - i; j >= 1; j--)
		{
			f *= cos(input[j-1] * pi / 2);
		}
		if (i > 1)
		{
			f *= sin(input[(output.size() - i + 1) - 1] * pi / 2);
		}

		output[i-1] = f;
	}



	return output;
}




ObjFunction2  eval_DTLZ3(Point2 &input, unsigned int const nObj_size)
{
	ObjFunction2 output(nObj_size,0.0);
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int n = input.size();
	unsigned int k = n - output.size() + 1;

	double g = 0;
	for (i = n - k + 1; i <= n; i++)
	{
		g += pow(input[i-1]-0.5,2) - cos(20 * pi * (input[i-1]-0.5));
	}
	g = 100 * (k + g);

	for (i = 1; i <= output.size(); i++)
	{
		double f = (1 + g);
		for (j = output.size() - i; j >= 1; j--)
		{
			f *= cos(input[j-1] * pi / 2);
		}
		if (i > 1)
		{
			f *= sin(input[(output.size() - i + 1) - 1] * pi / 2);
		}

		output[i-1] = f;
	}



	return output;
}


ObjFunction2  eval_DTLZ4(Point2 &input, unsigned int const nObj_size)
{
	ObjFunction2 output(nObj_size,0.0);
	unsigned int i = 0;
	unsigned int j = 0;
	double alpha = 100;
	unsigned int n = input.size();
	unsigned int k = n - output.size() + 1;

	double g = 0;
	for (i = n - k + 1; i <= n; i++)
	{
		g += pow(input[i-1]-0.5,2);
	}

	for (i = 1; i <= output.size(); i++)
	{
		double f = (1 + g);
		for (j = output.size() - i; j >= 1; j--)
		{
			f *= cos(pow(input[j-1],alpha) * pi / 2);
		}
		if (i > 1)
		{
			f *= sin(pow(input[(output.size() - i + 1) - 1],alpha) * pi / 2);
		}

		output[i-1] = f;
	}



	return output;
}


ObjFunction2  eval_DTLZ5(Point2 &input, unsigned int const nObj_size)
{
	ObjFunction2 output(nObj_size,0.0);
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int n = input.size();
	unsigned int k = n - output.size() + 1;
	double *theta = new double[output.size()];
	double t = 0;
	double g = 0;

	for (i = n - k + 1; i <= n; i++)
	{
		g += pow(input[i-1] - 0.5, 2);
	}

	t = pi / (4 * (1 + g));
	theta[0] = input[0] * pi / 2;
	for (i = 2; i <= output.size() - 1; i++)
	{
		theta[i-1] = t * (1 + 2 * g * input[i-1]);
	}

	for (i = 1; i <= output.size(); i++)
	{
		double f = (1 + g);
		for (j = output.size() - i; j >= 1; j--)
		{
			f *= cos(theta[j-1]);
		}
		if (i > 1)
		{
			f *= sin(theta[(output.size() - i + 1) - 1]);
		}

		output[i-1] = f;
	}

	free(theta);


	return output;
}


ObjFunction2  eval_DTLZ6(Point2 &input, unsigned int const nObj_size)
{
	ObjFunction2 output(nObj_size,0.0);
	unsigned int i = 0;
	int j = 0;
	unsigned int n = input.size();
	unsigned int k = n - output.size() + 1;
	double *theta = new double[output.size()] ;
	double t = 0;
	double g = 0;

	for (i = n - k + 1; i <= n; i++)
	{
		g += pow(input[i-1], 0.1);
	}

	t = pi / (4 * (1 + g));
	theta[0] = input[0] * pi / 2;
	for (i = 2; i <= output.size() - 1; i++)
	{
		theta[i-1] = t * (1 + 2 * g * input[i-1]);
	}

	for (i = 1; i <= output.size(); i++)
	{
		double f = (1 + g);
		for (j = output.size() - i; j >= 1; j--)
		{
			f *= cos(theta[j-1]);
		}
		if (i > 1)
		{
			f *= sin(theta[(output.size() - i + 1) - 1]);
		}

		output[i-1] = f;
	}

	free(theta);


	return output;
}


ObjFunction2  eval_DTLZ7(Point2 &input, unsigned int const nObj_size)
{
	ObjFunction2 output(nObj_size,0.0);
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int n = input.size();
	unsigned int k = n - output.size() + 1;
	double g = 0;
	double h = 0;

	for (i = n - k + 1; i <= n; i++)
	{
		g += input[i-1];
	}
	g = 1 + 9 * g / k;

	for (i = 1; i <= output.size() - 1; i++)
	{
		output[i-1] = input[i-1];
	}

	for (j = 1; j <= output.size() - 1; j++)
	{
		h += input[j-1] / (1 + g) * (1 + sin(3 * pi * input[j-1]));
	}
	h = output.size() - h;
	output[output.size() - 1] = (1 + g) * h;



	return output;
}



ObjFunction2  eval_COMET(Point2 &input)
{
	double x1;
	double x2;
	double x3;
	double g;
	ObjFunction2 output(3,0.0);

	assert(input.size() == 3);
	assert(output.size() == 3);

	x1 = 1 + (input[0] * 2.5);
	x2 = -2 + (input[1] * 4);
	x3 = input[2];

	g = x3;

	output[0] = (1 + g) * (pow(x1,3) * pow(x2,2) - 10 * x1 - 4 * x2);
	output[1] = (1 + g) * (pow(x1,3) * pow(x2,2) - 10 * x1 + 4 * x2);
	output[2] = 3 * (1 + g) * pow(x1,2);

	output[0] = output[0] + 100;
	output[1] = output[1] + 100;
	//output[2] = output[2];


	return output;
}

objective_function_formulae::objective_function_formulae(unsigned int nOfObjectives): __nObj(nOfObjectives) {}

ObjFunction2 objective_function_formulae::operator() (const Point2 &input){


	return eval_ZDT2(input);
}
