#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

const double pi = 3.14159265358979323846264338327950288;


void chebpolys(int N, double x, double *pols)
/*
 *  Return the values of the Chebyshev polynomials of degrees 0 through N
 *  at the point x.
 *
 *  Input parameters:
 *    N - an integer > 1 specifying the Chebyshev polynomial of largest degree
 *      which is to be evaluated
 *    x - the point on the interval [-1,1] at which to evaluate the Chebyshev
 *      polynomials
 *
 *  Output parameters:
 *    pols - an array of length N+1, which must be allocated by the caller,
 *     which will contain the values of the Chebyshev polynomials on return.
 *     More explicitly, vals[j] = T_j(x).
 *
 *  Return values:
 *    None
 *
 */
{
    int n;
    for (n=0; n <= N; n++){
        if (n == 0) {
            pols[n] = 1;
        }
        else if (n == 1) {
            pols[n] = x;
        }
        else {
            pols[n] = 2*x*pols[n-1] - pols[n-2];
        }
    }
}

void chebcoefs(int N, double a, double b, double  (*funptr)(double), double  *coefs)
/*
 *  Approximate the Chebyshev coefficients
 *
 *    a_0, a_1, ..., a_N
 *
 *  of a user-specified function f(x) given on the interval [a,b] using the
 *  (N+1)-point Chebyshev extrema quadrature rule.
 *
 *  Input parameters:
 *    N - the order of the Fourier series used to approximate f
 *    a - the left endpoint of the interval under consideration
 *    b - the right endpoint of hte interval under consideration
 *    funptr - a pointer to a function with calling syntax
 *
 *      double fun(double x);
 *
 *    which returns the value of the function f at the point x.
 *
 *  Output parameters:
 *    coefs - this array of length N+1, which must be preallocated by the
 *      caller, will contain the desired approximations upon return.  More
 *      explicitly, coefs[j] should be equal to a_j.
 *
 *  Return value:
 *     None
 *
 */
{
   double *pols;
   double sum, x_j, y, wha;
   pols = (double *) malloc( sizeof(double) * (N+1));
   int n, j;
   wha = 2.0/(N+1);
   for (int j=0; j<= N; j++)
   coefs[j] = 0.0;

   for ( n = 0; n <= N; n++){
        x_j = cos((n+0.5)*(pi)/(N+1));
        chebpolys(N,x_j,pols);
        y   = funptr(((b-a)/2.0)*x_j + ((b+a)/2.0));

    for (j = 0; j <= N; j++){
        coefs[j] = coefs[j] + y*pols[j]*wha;
    }
   }

   free(pols);

}

double chebeval(int N, double a, double b, double *coefs, double y)
/*
 *  Evaluate the Chebyshev series
 *
 *             N
 *            sum    a_n T_n(2/(b-a) y - (b+a)/(b-a))                      (1)
 *            n=0
 *
 *  at a user-specified point on the interval [a,b].
 *
 *  Input parameters:
 *    N - the integer N which specifies the order of the series
 *    a - the left endpoint of the interval
 *    b - the right endpoint of the interval
 *    coefs - an array of length N+1 such that coefs[j] = a_j
 *    y - the point on the interval [a,b] at which to evaluate (1)
 *
 *  Output parameters:
 *    None
 *
 *  Return value:
 *    The value of the expansion (1) at the point y.
 *
 */
{
   double val, x_2;
   double *pols;
   int n;

   pols = (double *) malloc( sizeof(double) * (N+1) );
   val = 0;
   x_2 = 2*y/(b-a)- (b+a)/(b-a);
   chebpolys(N,x_2,pols);
   for(n = 0; n <= N; n++){
        if (n ==0){
        val = val + coefs[n]*(0.5)*pols[n];
        }
        else{
        val = val + coefs[n]*pols[n];
        }
   }
   return val;
}
