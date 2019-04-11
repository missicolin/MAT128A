#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* function prototypes */

void chebcoefs(int N, double a, double b, double  (*funptr)(double), double  *coefs);
double chebeval(int N, double a, double b, double *coefs, double y);


int ifsplit(double eps, int N, double *coefs);


void chebadap(int *ier, double eps, int maxints, int N, double a, double b,
  double (*funptr)(double), int *m, double *as)
/*
 *  Use an adaptive discretization procedure to find a partition
 *
 *     a = a_0 < a_1 < ... < a_m                                         (1)
 *
 *  such that the restriction of a user-specified function f(x) to each of
 *  the intervals [a_j, a_{j+1}] is represented with precision eps by a
 *  (N+1)-term Chebyshev expansion.
 *
 *  Input parameters:
 *    eps - the desired precision for the Chebyshev expansions
 *    maxints - the maximum possible number of intervals
 *    N - the order of the Chebyshev expansions used to represent f on each
 *      subinterval
 *    (a,b) - the extents of the interval on which f is given
 *    funptr - a pointer to the function f
 *
 *  Output parameters:
 *    ier - an error return code;
 *        ier = 0    means that the procedure was successful
 *        ier = 4    means that the maximum number of intervals was exceeded
 *        ier = 8    means that an interval of length less than 2^(-15) was encountered
 *    m - the number of intervals, provided ier = 0
 *    as -  this user-allocated array, which must be of length maxints+1,
 *     will contain the points in the partition (1) --- more explicitly,
 *     it will be the case that as[j] = a_j for j= 0,1,...,m
 *
 */
{
  double *as0, *coefs;
  double a0,b0,c0;
  int nas0, nas1;

  *ier = 0;
  *m   = 0;

  /* as0 is the list of intervals which must still be considered */

  as0   = (double *)malloc( sizeof(double) * 2*maxints );
  coefs = (double *)malloc( sizeof(double) * (N+1) );

  nas0  = 1;
  nas1  = 0;

  as0[0] = a;
  as0[1] = b;

  while( nas0 > 0)
  {
      a0 = as0[2*nas0-2];
      b0 = as0[2*nas0-1];
      c0 = (a0+b0)/2.0;

      if ( (b0-a0) < 1.0e-15) {
        *ier = 8;
	*m   = -1;
	return;
      }

      nas0 = nas0 - 1;
//   print_array("a0 = ",1,&a0);
//      print_array("b0 = ",1,&b0);
//      print_array("c0 = ",1,&c0);

      chebcoefs(N, a0, b0, funptr, coefs);

      if ( ifsplit(eps,N,coefs)  ) {

	if (nas0 +2 > maxints) {
	  *ier = 4;
	  *m   = -1;
	  return;
	}

	nas0 += 1;
        as0[2*nas0-2] = c0;
        as0[2*nas0-1] = b0;

	nas0 += 1;
        as0[2*nas0-2] = a0;
        as0[2*nas0-1] = c0;
      } else {
	*m = *m +1;
	if ( (*m) > maxints ) {
	  *ier = 4;
	   *m  = -1;
	   return;
	}

	if (*m == 1) {
	  as[0] = a0;
	}
	as[*m] = b0;
      }
  }

  free(coefs);
  free(as0);
}


int ifsplit(double eps, int N, double *coefs)
/*
 *  Determine whether chebadap should split a specified interval or not
 *  given the coefficients in the local Chebyshev expansion.
 *
 *  Input parameters:
 *    eps - the precision parameter passed to chebadap
 *    N - the order of the Chebyshev expansion
 *    coefs - an array of length (N+1) giving the coefficients in the
 *     Chebyshev expansion
 *
 *  Output parameters:
 *    None
 *
 *  Return value:
 *    0 if the interval is not to be split and 1 if it is to be split
 */
{
  double dd1, dd2, dd;
  int i;

  dd1 = 0.0;
  dd2 = 0.0;


  for(i=0; i<= N/2; i++) {
    dd1 = dd1 + coefs[i]*coefs[i];
  }

  for(; i<= N; i++) {
    dd1 = dd1 + coefs[i]*coefs[i];
    dd2 = dd2 + coefs[i]*coefs[i];
  }

  dd = dd2/dd1;

  if(dd < eps*eps)
    return 0;


  return 1;

}



void chebadap_coefs(int N, int m, double *as, double (*funptr)(double), double *coefs)
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
 *    coefs - an (N+1) * m array which must be allocated by the caller and
 *      which, upon return, contains the coefficients in the Chebyshev expansions;
 *      more explicitly, the entries j*(N+1) through (j+1)*(N+1) will contain the
 *      coefficients of the expansion on the j^th  interval [as[j], as[j+1]].
 *
 */
{

int ind;
int i;
for (i = 0; i <= (m-1); i++){
    ind = i*(N+1);
    chebcoefs(N,as[i],as[i+1],funptr,coefs + ind);
}
}

double chebadap_eval(int N, int m, double *as, double *coefs, double x)
/*
 *  Evaluate a piecewise Chebyshev expansion at a specified point.
 *
 *  Input parameters:
 *    N - the order of the Chebyshev expansions on each interval
 *    m - the number of intervals
 *    as - an array of length m+1 specifying the endpoints of the intervals
 *    coefs - the array returned by chebadap_coefs which contains the
 *      coefficients of the local Chebyshev expansions
 *
 *  Output parameters:
 *    x - the point at which to evaluate the piecewise Chebyshev expansion
 *
 *  Return values:
 *    the value of the expansion at the point x
 *
 */
{
int ind;
int i;
double val;
for (i = 0; i <= (m-1); i++){
	if ( x >= as[i] && x <= as[i+1]) {
	ind = i*(N+1);
	val = chebeval(N,as[i],as[i+1],coefs + ind,x);
	}
}
	return val;
}
