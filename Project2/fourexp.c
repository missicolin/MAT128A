#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

const double pi = 3.14159265358979323846264338327950288;


void fourcoefs(int N, double complex (*funptr)(double), double complex *coefs)
/*
 *  Approximate the Fourier coefficients
 *
 *    a_{-N}, a_{-N+1}, ..., a_{-1}, a_0, a_1, ... a_{N-1}, a_N
 *
 *  of a user-specified function f(x) using the (2N+1)-point periodic
 *  trapezoidal rule on the interval [-\pi,\pi].  The function is specified
 *  by a pointer to a function (see below).
 *
 *  Input parameters:
 *    N - the order of the Fourier series used to approximate f
 *    funptr - a pointer to a function with calling syntax
 *
 *      double complex fun(double complex x);
 *
 *    which returns the value of the function f at the point x.
 *
 *  Output parameters:
 *    coefs - this array of length 2N+1, which must be preallocated by the
 *      caller, will contain the desired approximations upon return.  More
 *      explicitly, coefs[j] should be equal to a_{-N+j} for each j=0,1,..,2N.
 *
 *  Return value:
 *     None
 *
 */
{

/* Define number of nodes, step size, quadrature weight, and other variables
Also define "double complex y
You might find the line " y = funptr(x); " useful, where the variable x is the quadrature node. */

int j, m;
double complex x,y,sum;

for(m=-N;m <= N; m++) {

  sum = 0;
  for(j=0;j<=2*N;j++) {
    x   = -pi + 2*pi/(2*N+1)*j;
    y   = funptr(x);
    sum = sum + y*cexp(-I*m*x);
  }

  coefs[m+N]=sum/(2*N+1);
}

}

double complex foureval(int N, double complex *coefs, double t)
/*
 *  Evaluate the Fourier series
 *
 *             N
 *            sum    a_n \exp(int)                                 (1)
 *            n=-N
 *
 *  at a user-specified point on the interval [-pi,pi].
 *
 *  Input parameters:
 *    N - the integer N which specifies the order of the series
 *    coefs - an array of length 2*N+1 such that coefs[j] = a_{j-N}
 *    t - the point in the interval [-pi,pi] at which to evaluate the
 *      series
 *
 *  Output parameters:
 *    None
 *
 *  Return value:
 *    The value of the expansion (1) at the point t.
 *
 */
{
int n, j;
double complex res;
res = 0.0;
for (j = 0; j <= 2*N; j++){
    n = j-N;
    res = res + coefs[j]*cexp(I*n*t);
}
return res;
}
