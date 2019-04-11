#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* function prototypes */
const double pi = 3.14159265358979323846264338327950288;

void chebadap(int *ier, double eps, int maxints, int N, double a, double b,
  double (*funptr)(double), int *m, double *as);
void chebadap_coefs(int N, int m, double *as, double (*funptr)(double), double *coefs);
double chebadap_eval(int N, int m, double *as, double *coefs, double x);
double chebadap_int(int N, int m, double *as, double *coefs);

void print_array(char *str, int n, double *vals);
void print_intarray(char *str, int n, int *vals);


double fun1(double x)
{
  return cos(x*x);
}


double fun2(double x)
{
  return exp(cos(x*x));
}

double fun3(double x)
{
  return 1.0/(1.0e-7+x*x);
}

double fun4(double x)
{
  return 2+sin(27*x*x);
}

double fun5(double x)
{
  return 1.0/(25*x*x+1);
}


void print_array(char *str, int n, double *vals)
{
/*
 *  Utility routine which prints an array of doubles in a formatted column.
 *
 *  Input parameters:
 *    str - a string to print
 *    n - the length of the array of doubles to print
 *    vals - the array of length n containing the values to print
 *
 *  Output parameters:
 *    None
 *
 *  Return value:
 *    None
 *
 */

  int nlines,i,idx,nn;

  nlines = n/3;
  nn     = n - 3*nlines;
  idx    = 0;

  printf("%s\n",str);

  for(i =1; i <= nlines; i++) {
    printf("%24.15e   %24.15e   %24.15e\n",vals[idx],vals[idx+1],vals[idx+2]);
    idx+=3;
  }

  for(i=1; i<= nn; i++) {
    printf("%24.15e   ",vals[idx]);
    idx+=1;
    printf("\n");
  }

}


void print_intarray(char *str, int n, int *vals)
{
/*
 *  Utility routine which prints an array of integers in a formatted column.
 *
 *  Input parameters:
 *    str - a string to print
 *    n - the length of the array of integers to print
 *    vals - the array of length n containing the values to print
 *
 *  Output parameters:
 *    None
 *
 *  Return value:
 *    None
 *
 */

  int nlines,i,idx,nn;

  nlines = n/3;
  nn     = n - 3*nlines;
  idx    = 0;

  printf("%s\n",str);

  for(i =1; i <= nlines; i++) {
    printf("   %08d   %08d   %08d\n",vals[idx],vals[idx+1],vals[idx+2]);
    idx+=3;
  }

  for(i=1; i<= nn; i++) {
    printf("   %08d   ",vals[idx]);
    idx+=1;
    printf("\n");
  }

}

double adap_error_check(int N, double a, double b, int m, double *as,
  double *coefs, double (*funptr)(double))
{
/*
 *  Test the error in an adaptive discretization of a specified function.
 *
 *  Input parameters:
 *    N - the order of the  local Chebyshev expansions
 *    m - the number of interval
 *    as - an array of length m+1 giving the partition points
 *    coefs - an  m by (N+1) matrix whose jth column gives the coefficients
 *      for the jth interval as[j], as[j+1]
 *    funptr - a pointer to the function
 *
 *  Output parameters:
 *     None
 *
 *  Return value:
 *    The largest observed relative error.
 *
 */

  int i,nn;
  double x, y, y0, errmax;

  errmax = 0.0;
  nn     = 1000;

  for (i=0; i <= nn; i++) {
    x      = a + (b-a) * i/nn;
    y0     = funptr(x);
    y      = chebadap_eval(N, m, as, coefs, x);
    errmax = fmax(errmax,fabs(y-y0)/(fabs(y0)));

  }

  return errmax;

}

int main(int argc, char **argv)
{
  int maxints, ier, m, N;
  double eps, *as, *coefs, a, b, errmax, rint;
  int nscore;

  eps     = 1.0e-13;
  maxints = 10000;
  nscore  = 0;



  /* Adaptively discretize f(x) = cos(x^2) */

  N      = 20;
  a      = -1.0;
  b      =  1.0;

  printf("===[ f(x) = cos(x^2)           ]====================================================\n");

  as    = (double *) malloc( sizeof(double) * (maxints+1) );
  coefs = (double *) malloc( sizeof(double) * (N+1)*(maxints) );

  chebadap(&ier,eps,maxints,N,a,b,fun1,&m,as);
  print_intarray("after chebadap, ier = ",1,&ier);
  if (ier != 0) exit(0);

  print_intarray("after chebadap, m = ",1,&m);
  print_array("after chebadap, as = ",m+1,as);
  chebadap_coefs(N, m, as, fun1, coefs);
  errmax = adap_error_check(N,a,b,m,as,coefs,fun1);
  print_array("errmax = ",1, &errmax);

  if (errmax < 1.0e-12) {
    nscore = nscore+1;
    printf("===[ PASSED                   ]====================================================\n");
  }
  else {
    printf("===[ FAILED                   ]====================================================\n");
  }
  printf("\n");

  free(as);
  free(coefs);



  /* Adaptively discretize f(x) = cos(x^2) */

  N      = 12;
  a      =  2.0;
  b      =  3.0;

  printf("===[ f(x) = exp(cos(x)^2)       ]====================================================\n");

  as    = (double *) malloc( sizeof(double) * (maxints+1) );
  coefs = (double *) malloc( sizeof(double) * (N+1)*(maxints) );

  chebadap(&ier,eps,maxints,N,a,b,fun2,&m,as);
  print_intarray("after chebadap, ier = ",1,&ier);
  if (ier != 0) exit(0);

  print_intarray("after chebadap, m = ",1,&m);
  print_array("after chebadap, as = ",m+1,as);
  chebadap_coefs(N, m, as, fun2, coefs);
  errmax = adap_error_check(N,a,b,m,as,coefs,fun2);
  print_array("errmax = ",1, &errmax);

  if (errmax < 1.0e-12) {
    nscore = nscore+1;
    printf("===[ PASSED                   ]====================================================\n");
  }
  else {
    printf("===[ FAILED                   ]====================================================\n");
  }
  printf("\n");

  free(as);
  free(coefs);




  /* Adaptively discretize f(x) = 1/(10^(-7)+x^2) */

  N      =  30;
  a      = -1.0;
  b      =  1.0;

  printf("===[ f(x) = 1/(10^(-7)+x^2)     ]====================================================\n");

  as    = (double *) malloc( sizeof(double) * (maxints+1) );
  coefs = (double *) malloc( sizeof(double) * (N+1)*(maxints) );

  chebadap(&ier,eps,maxints,N,a,b,fun3,&m,as);
  print_intarray("after chebadap, ier = ",1,&ier);
  if (ier != 0) exit(0);

  print_intarray("after chebadap, m = ",1,&m);
  print_array("after chebadap, as = ",m+1,as);
  chebadap_coefs(N, m, as, fun3, coefs);
  errmax = adap_error_check(N,a,b,m,as,coefs,fun3);
  print_array("errmax = ",1, &errmax);

  if (errmax < 1.0e-12) {
    nscore = nscore+1;
    printf("===[ PASSED                   ]====================================================\n");
  }
  else {
    printf("===[ FAILED                   ]====================================================\n");
  }
  printf("\n");

  free(as);
  free(coefs);



  /* Adaptively discretize f(x) = sin(27*x) */

  N      =  30;
  a      =  0.0;
  b      =  3.0;

  printf("===[ f(x) = 2+ sin(27*x^2)       ]====================================================\n");

  as    = (double *) malloc( sizeof(double) * (maxints+1) );
  coefs = (double *) malloc( sizeof(double) * (N+1)*(maxints) );

  chebadap(&ier,eps,maxints,N,a,b,fun4,&m,as);
  print_intarray("after chebadap, ier = ",1,&ier);
  if (ier != 0) exit(0);

  print_intarray("after chebadap, m = ",1,&m);
  print_array("after chebadap, as = ",m+1,as);
  chebadap_coefs(N, m, as, fun4, coefs);
  errmax = adap_error_check(N,a,b,m,as,coefs,fun4);
  print_array("errmax = ",1, &errmax);

  if (errmax < 1.0e-12) {
    nscore = nscore+1;
    printf("===[ PASSED                   ]====================================================\n");
  }
  else {
    printf("===[ FAILED                   ]====================================================\n");
  }
  printf("\n");

  free(as);
  free(coefs);



  /* Adaptively discretize f(x) = Runge Example */

  N      =  30;
  a      = -1.0;
  b      =  1.0;

  printf("===[ f(x) = Runge's Example     ]====================================================\n");

  as    = (double *) malloc( sizeof(double) * (maxints+1) );
  coefs = (double *) malloc( sizeof(double) * (N+1)*(maxints) );

  chebadap(&ier,eps,maxints,N,a,b,fun5,&m,as);
  print_intarray("after chebadap, ier = ",1,&ier);
  if (ier != 0) exit(0);

  print_intarray("after chebadap, m = ",1,&m);
  print_array("after chebadap, as = ",m+1,as);
  chebadap_coefs(N, m, as, fun5, coefs);
  errmax = adap_error_check(N,a,b,m,as,coefs,fun5);
  print_array("errmax = ",1, &errmax);

  if (errmax < 1.0e-12) {
    nscore = nscore+1;
    printf("===[ PASSED                   ]====================================================\n");
  }
  else {
    printf("===[ FAILED                   ]====================================================\n");
  }
  printf("\n");

  free(as);
  free(coefs);

printf("\n\n");
printf("SCORE = %d / %d\n",nscore,5);
return nscore;


}
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
