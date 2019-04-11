#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

//const double pi = 3.14159265358979323846264338327950288;

void chebadap(int *ier, double eps, int maxints, int N, double a, double b,
		  double (*funptr)(double), int *m, double *as);
void chebadap_coefs(int N, int m, double *as, double (*funptr)(double), double *coefs);
double chebadap_eval(int N, int m, double *as, double *coefs, double x);
double chebadap_int(int N, int m, double *as, double *coefs);



double chebint(int N, double a, double b, double *coefs)
{
/*
 *  Evaluate the definite integral
 *
 *         b
 *    \int    p(x) dx
 *         a
 *
 *  where
 *
 )*               N
 *      p(x) = \sum  coefs[n]  T_n(2/(b-a) y - (b+a)/(b-a))
 *              n=0
 *
 *
 *
 *  Input parameters:
 *    N - the integer N which specifies the order of the series
 *    a - the left endpoint of the interval
 *    b - the right endpoint of the interval
 *    coefs - an array of length N+1 such that coefs[j] = a_j
 *
 *  Output parameters:
 *    None
 *
 *  Return value:
 *    The value of the integral of the
 *
 */
int n;
double sum= 0.0, it;
it = 0.5*(b-a);
for (n=0; n <= N/2; n++){
	if (n == 0){
		sum += (0.5)*it*(2*coefs[2*n])/(1-4*(n*n));
			}
	else{
		sum += it*(2*coefs[2*n])/(1-4*(n*n));
	}
}
return sum;
}

double chebadap_int(int N, int m, double *as, double *coefs)
{
/*
 *  Compute the coefficients in a piecewise Chebyshev expansion of a user-specified
 *  function f.  The collection of intervals on which the Chebyshev expansions
 *  are to be computed is specified by the input array as.  The endpoints of the
 *  jth interval are as[j] and as[j+1].
 *
 *  The coefficients for the j^th interval should be stored in the positions
 *
 *    j*(N+1), ..., (j+1)*(N+1)
 *
 *  in the coefs array, which is a return parameter for this function.
 *
 *  Input parameters:
 *    N - the order of the Chebyshev expansions on each interval
 *    m - the number of intervals
 *    as - an array of length m+1 specifying the endpoints of the intervals
 *    funptr - a pointer to the input function
 *
 *  Output parameters:
 *    None
 *
 *  Return value:
 *
 */
int j,n, ind;

double sum = 0.0, con;

for (j=0; j <= (m-1); j++){
	con = 0.5*(as[j+1] - as[j]);
	ind = j*(N+1);
	for(n=0 ; n <= N/2; n++){
		if (n == 0){
			sum += 0.5*con*(2*coefs[2*n + ind]/(1-4*(n*n)));
		}
		else{
			sum += con*(2*coefs[2*n + ind]/(1-4*(n*n)));
		}
	}
}
return sum;
}
