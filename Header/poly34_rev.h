// poly34.h : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com

// Jun/16/2018: Retrieved from http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// Jun/17/2018: Revised by Tadamasa Sawada and Vasily Minkov
// Aug/11/2019 (ver2): Minor revisions by Tadamasa Sawada and Vasily Minkov


#ifndef poly34
#define poly34

int   SolveP2(double *x, double a, double b);			// solve equation x^2 + a*x + b = 0
//int   SolveP3(double *x, double a, double b, double c);			// solve cubic equation x^3 + a*x^2 + b*x + c = 0
//int   SolveP4(double *x, double a, double b, double c, double d);	// solve equation x^4 + a*x^3 + b*x^2 + c*x + d = 0 by Dekart-Euler method

int   SolveP3_ver2(double *x,double a,double b,double c);			// solve cubic equation x^3 + a*x^2 + b*x + c = 0
int   SolveP4_ver2(double *x,double a,double b,double c,double d);	// solve equation x^4 + a*x^3 + b*x^2 + c*x + d = 0 by Dekart-Euler method

//Shared by ver1 and ver2
void  CSqrt(double x, double y, double &a, double &b); // returns:  a+i*s = sqrt(x+i*y)
double root3(double x);
int   SolveP4Bi(double *x, double b, double d);

#endif