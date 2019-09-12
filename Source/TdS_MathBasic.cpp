#include <math.h>
#include <cmath>
#include "TdS_Math.h"

const double eps = 1e-14;


double TdS_RandomNumber(double start, double end)
{
	return (double)rand()*(end - start) / RAND_MAX + start;
}

void TdS_RandomNumbers(double start, double end, int n_v, double *v)
{
	for (int i = 0; i < n_v; i++)
		v[i] = TdS_RandomNumber(start, end);
}

void TdS_BubbleSort(double *x, int n) // Sort: x[0] <= x[1] <= ...
{
	for (int i = 0; i<n - 1; i++)
		for (int j = n - 1; j>i; j--)
			if (x[j]<x[j - 1])
			{
				double t = x[j];
				x[j] = x[j - 1];
				x[j - 1] = t;
			}
}

int	TdS_Min(int num_v, int *v, int *i_min)
{
	int temp_min = v[0];
	int temp_i_min = 0;

	for (int i = 1; i < num_v; i++)
		if (v[i] < temp_min)
		{
			temp_min = v[i];
			temp_i_min = i;
		}
	if (i_min != NULL) *i_min = temp_i_min;
	return temp_min;
}

double	TdS_Min(int num_v, double *v, int *i_min)
{
	double temp_min = v[0];
	int temp_i_min = 0;

	for (int i = 1; i < num_v; i++)
		if (v[i] < temp_min)
		{
			temp_min = v[i];
			temp_i_min = i;
		}
	if (i_min != NULL) *i_min = temp_i_min;
	return temp_min;
}

int	TdS_Max(int num_v, int *v, int *i_max)
{
	int temp_max = v[0];
	int temp_i_max = 0;

	for (int i = 1; i < num_v; i++)
		if (v[i] > temp_max)
		{
			temp_max = v[i];
			temp_i_max = i;
		}
	if (i_max != NULL) *i_max = temp_i_max;
	return temp_max;
}

double	TdS_Max(int num_v, double *v, int *i_max)
{
	double temp_max = v[0];
	int temp_i_max = 0;

	for (int i = 1; i < num_v; i++)
		if (v[i] > temp_max)
		{
			temp_max = v[i];
			temp_i_max = i;
		}
	if (i_max != NULL) *i_max = temp_i_max;
	return temp_max;
}

int TdS_SimultaneousEquations3(double *x, double a1, double a2, double a3, double a0, double b1, double b2, double b3, double b0, double c1, double c2, double c3, double c0)
{
	//After http://cplusplus.happycodings.com/mathematics/code5.html

	//a1*x + a2*y + a3*z + a0 = 0
	//b1*x + b2*y + b3*z + b0 = 0
	//c1*x + c2*y + c3*z + c0 = 0

	double denominator = (a1*b2*c3 + a2*c1*b3 + a3*b1*c2) - (a1*b3*c2 + a2*b1*c3 + a3*b2*c1);
	if (fabs(denominator) < eps)
		return -1; // error

	x[0] = ((a2*c3*b0 + a3*b2*c0 + a0*b3*c2) - (a2*b3*c0 + a3*c2*b0 + a0*b2*c3)) / denominator;
	x[1] = ((a1*b3*c0 + a3*c1*b0 + a0*b1*c3) - (a1*c3*b0 + a3*b1*c0 + a0*b3*c1)) / denominator;
	x[2] = ((a1*c2*b0 + a2*b1*c0 + a0*b2*c1) - (a1*b2*c0 + a2*c1*b0 + a0*b1*c2)) / denominator;

	return 0; // fine
}