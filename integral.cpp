/*
 * =====================================================================================
 *
 *       Filename:  tema5.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/15/14 12:29:18
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;

const double epsilon = 0.0000000001;

FILE *in, *out;

double f(double x)
{
	return (pow(x, 7) * sqrt(1 - x * x)) / (pow(2 - x, (double)13 / 2));
}


double integrateRectangle(int a, int b, int n)
{
	double h, result, t; // result este A din pdf, t este B din pdf
	h = (b - a) / n;

	for (int i = 0; i < n; i++)
		result = result + h * f(a + i * h);

	do
	{
		t = result;
		n = 2 * n;
		h = (double)(b - a) / n;
		result = 0;
		
		for (int i = 0; i < n; i++)
			result = result + h * f(a + i * h);
	}
	while (abs(t - result) < epsilon);

	return result;
}


double integrateTrapeze(int a, int b, int n)
{
	double h, result, t;
	
	h = (double)(b - a) / n;
	result = 0;
	result = result + (h / 2) * f(a);

	for (int i = 1; i < n; i++)
		result = result + h * f(a + i * h);

	result = result + (h / 2) * f(b);

	do
	{
		t = result;
		n = 2 * n;
		h = (double)(b - a) / n;
		result = 0;
		result = result + (h / 2) * f(a);

		for (int i = 1; i < n ; i++)
			result = result + h * f(a + i * h);

		result = result + (h / 2) * f(b);
	}
	while (abs(t - result) < epsilon);

	return result;
}


double integrateSimpson(int a, int b, int n) 
{
	double rez = 0, d = (double)(b - a) / n;
	double xi, xi1;

	for (int i = 0; i < n; i++) 
	{
		xi = a + i * d;
		xi1 = a + (i + 1) * d;
		rez += ((b - a) / (6.0 * n)) * (f(xi) + 4.0 * f((xi + xi1) / 2.0) + f(xi1));
	}
	return rez;

}

double newtonCotes2(double a, double b, double n) 
{
	double rez = 0, d = (b - a) / n;
	double xi, xi1;

	for (int i = 0; i < n; i++) 
	{
		xi = a + i * d;
		xi1 = a + (i + 1) * d;
		rez += ((b - a) / (2.0 * n)) *
		       (f((2.0 * xi + xi1) / 3.0) + f((xi + 2 * xi1) / 3.0));
	}
	return rez;
}

double newtonCotes3(double a, double b, double n) 
{
	double rez = 0, d = (b - a) / n;
	double xi, xi1;

	for (int i = 0; i < n; i++) 
	{
		xi = a + i * d;
		xi1 = a + (i + 1) * d;
		rez += ((b - a) / (3.0 * n)) * 
			(2.0 * f((3.0 * xi + xi1) / 4.0) - 
			f((xi + xi1) / 2.0) + 2.0*f((3 * xi + xi1) / 4.0));
	}
	return rez;
}

double gauss2(double a, double b, double n) 
{
	double rez = 0, d = (b - a)/n;
	double yi, yi1, ai, bi;
	
	for (int i = 1; i <= n; i++) 
	{
		yi1 = a + (i - 1)*d;
		yi = a + i * d;
		ai = (yi1 + yi) / 2.0 - (b - a) / (2.0 * n * sqrt(3));
		bi = (yi1 + yi) / 2.0 + (b - a) / (2.0 * n * sqrt(3));

		rez += (d / 2.0) * (f(ai) + f(bi));
	}
	return rez;
}

double gauss3(double a, double b, double n) 
{
	double rez = 0, d = (b - a) / n;
	double yi, yi1, ai, bi, ci;
	
	for (int i = 1; i <= n; i++) 
	{
		yi1 = a + (i - 1)*d;
		yi = a + i * d;
		ai = (yi1 + yi) / 2.0 - ((b - a) * sqrt(15.0)) / (10.0 * n);
		bi = (yi1 + yi) / 2.0; 
		ci = (yi1 + yi) / 2.0 + ((b - a) * sqrt(15.0)) / (10.0 * n);

		rez += (d / 18.0) * (5*f(ai) + 8*f(bi) + 5*f(ci));
	}
	return rez;
}



double integrateBoole(double StartPoint, double EndPoint, int n) 
{
	vector<double> X(n + 1, 0.0);
	vector<double> Y(n + 1, 0.0);
	double delta_x = (EndPoint - StartPoint) / double(n);
	for (int i=0; i<=n; i++) 
	{
		X[i] = StartPoint + i*delta_x;
		Y[i] = f(X[i]);
	}
	double sum = 0;
	for (int t = 0; t <= (n - 1) / 4; t++) 
	{
		int ind = 4 * t;
	    sum += (1 / 45.0) * (14 * Y[ind] + 64 * Y[ind + 1] + 
			    24 * Y[ind + 2] + 64 * Y[ind + 3] + 14 * Y[ind + 4]) * delta_x;
	}
	return sum;
}



int main()
{
	int a, b, n;
	
	in = fopen("integral.in", "r");
	out = fopen("integral.out", "w");
	
		// integral from a to b, with n iterations	
	fscanf(in, "%d %d %d", &a, &b, &n);	
	
	fprintf(out, "--------- n = %d --------\n", n);
	
	fprintf(out, "\nRectangle Method\n");
	fprintf(out, "%.10lf\n", integrateRectangle(a, b, n));
		
	fprintf(out, "\nTrapeze Method\n");
	fprintf(out, "%.10lf\n", integrateTrapeze(a, b, n));
	
	fprintf(out, "\nSimpson's Method\n");
	fprintf(out, "%.10lf\n", integrateSimpson(a, b, n));
	
	fprintf(out, "\nBoole's Method\n");
	fprintf(out, "%.10lf\n", integrateBoole(a, b, n));
	
	fprintf(out, "\nNewton-Cotes with 2 nodes method\n");
	fprintf(out, "%.10lf\n", newtonCotes2(a, b, n));
	
	fprintf(out, "\nNewton-Cotes with 3 nodes method\n");
	fprintf(out, "%.10lf\n", newtonCotes3(a, b, n));	
	
	fprintf(out, "\nGauss with 2 nodes method\n");
	fprintf(out, "%.10lf\n", gauss2(a, b, n));
	
	fprintf(out, "\nGauss with 3 nodes method\n");
	fprintf(out, "%.10lf\n", gauss3(a, b, n));
		

	fclose(in);
	fclose(out);
 
	return 0;
}
