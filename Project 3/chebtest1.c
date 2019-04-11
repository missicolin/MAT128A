#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void chebpolys(int N, double x, double *pols);
double chebeval(int N, double a, double b, double *coefs, double y);
void chebcoefs(int N, double a, double b, double  (*funptr)(double), double  *coefs);

const double pi = 3.14159265358979323846264338327950288;

double relative_error(int N, double  *coefs, double  *coefs0)
{
  double derr,dnorm,dd;
  int i;
  derr  = 0.0;
  dnorm = 0.0;

  for(i=0; i<= N; i++) {
    dd    = fabs(coefs[i]-coefs0[i]);
    derr  = derr + dd*dd;
    dd    =  fabs(coefs0[i]);
    dnorm = dnorm + dd*dd;
  }
  derr = sqrt(derr / dnorm);
  return derr;
}

void print_coefs(int N, double  *coefs)
{
  int i;
  for(i=0; i <= N; i++) {
     printf("a_{%4.3d} = %24.16f\n",i,coefs[i]);
  }
}


double fun0(double x)
{
  return x*x;
}



double fun1(double x)
{
  return 2.0*(-2.0+x) / (-5 + 4*x);
}


double fun2(double x)
{
  return exp(x);
}


double fun3(double x)
{
  return cos(13/(0.01+x*x));
}


int main(int argc, char **argv)
{
int nscore;
int N, i;
double *coefs, *coefs0, val, val0;
double y,errrel,a,b;

nscore = 0;


/* Approximate the Chebyshev coefficients of f(x) = x^2 on [-1,1] */

N      = 4;
a      = -1.0;
b      =  1.0;
coefs  = (double *)malloc(sizeof(double)*(N+1));
coefs0 = (double *)malloc(sizeof(double)*(N+1));


for(i=0; i<= N; i++) {
  coefs0[i] = 0;
}

coefs0[0] = 1.0;
coefs0[2] = 0.5;

chebcoefs(N,a,b,fun0,coefs);
 print_coefs(N,coefs);


printf("\n");
printf("===[ f(x) = x^2               ]====================================================\n");
errrel = relative_error(N,coefs,coefs0);

printf("N                              = %7.4d\n",N);
printf("a                              = %7.4f\n",a);
printf("b                              = %7.4f\n",b);

printf("relative error in coefficients = %24.16g\n",errrel);

if (errrel < 1.0e-13)  {
  printf("===[ PASSED                   ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                   ]====================================================\n");
printf("\n");


free(coefs);
free(coefs0);



/* Approximate the Chebyshev coefficients of f(x) = x^2 on [3,4] */

N      = 4;
a      = 3.0;
b      = 4.0;
coefs  = (double *)malloc(sizeof(double)*(N+1));
coefs0 = (double *)malloc(sizeof(double)*(N+1));


for(i=0; i<= N; i++) {
  coefs0[i] = 0;
}

coefs0[0] = 24.75;
coefs0[1] = 3.5;
coefs0[2] = 0.125;

chebcoefs(N,a,b,fun0,coefs);

printf("\n");
printf("===[ f(x) = x^2               ]====================================================\n");
errrel = relative_error(N,coefs,coefs0);

printf("N                              = %7.4d\n",N);
printf("a                              = %7.4f\n",a);
printf("b                              = %7.4f\n",b);
printf("relative error in coefficients = %24.16g\n",errrel);

if (errrel < 1.0e-13)  {
  printf("===[ PASSED                   ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                   ]====================================================\n");
printf("\n");


free(coefs);
free(coefs0);




/* Approximate the Chebyshev coefficients of f(x) = 2(-2+x)/(-5+4x) on [-1,1] */

N      = 60;
a      =-1.0;
b      = 1.0;
coefs  = (double *)malloc(sizeof(double)*(N+1));
coefs0 = (double *)malloc(sizeof(double)*(N+1));


for(i=0; i<= N; i++) {
  coefs0[i] = pow(2.0,-i);
}
coefs0[0]=coefs0[0]*2;

chebcoefs(N,a,b,fun1,coefs);
//print_coefs(N,coefs);

printf("\n");
printf("===[ f(x) = 2(-2+x)/(-5+4x)   ]====================================================\n");
errrel = relative_error(N,coefs,coefs0);

printf("N                              = %7.4d\n",N);
printf("a                              = %7.4f\n",a);
printf("b                              = %7.4f\n",b);
printf("relative error in coefficients = %24.16g\n",errrel);

if (errrel < 1.0e-13)  {
  printf("===[ PASSED                   ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                   ]====================================================\n");
printf("\n");


free(coefs);
free(coefs0);



/* Evaluate the approximate series for  f(x) = exp(x) when x = 1 */

N      = 30;
a      = 0.0;
b      = 1.0;
y      = 1.0;
val0   = fun2(y);

coefs  = (double *)malloc(sizeof(double)*(N+1));

chebcoefs(N,a,b,fun2,coefs);
val = chebeval(N,a,b,coefs,y);

printf("\n");
printf("===[ f(x) = exp(x)            ]====================================================\n");
errrel = fabs(val-val0)/fabs(val0);

printf("N                              = %7.4d\n",N);
printf("a                              = %7.4f\n",a);
printf("b                              = %7.4f\n",b);
printf("relative error in evaluation   = %24.16g\n",errrel);

if (errrel < 1.0e-13)  {
  printf("===[ PASSED                ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                ]====================================================\n");
printf("\n");

//print_coefs(N,coefs);

free(coefs);




/* Evaluate the approximate series for  f(x) = cos(13/(.01+x**2)) when x = 1/10 */

N      = 3000;
a      = 0.0;
b      = 1.0;
y      = 0.1;
val0   = fun3(y);

coefs  = (double *)malloc(sizeof(double)*(N+1));

chebcoefs(N,a,b,fun3,coefs);
val = chebeval(N,a,b,coefs,y);

printf("\n");
printf("===[ f(x) = cos(13/(.01+x**2)) ]====================================================\n");
errrel = fabs(val-val0)/fabs(val0);

printf("N                              = %7.4d\n",N);
printf("a                              = %7.4f\n",a);
printf("b                              = %7.4f\n",b);
printf("relative error in evaluation   = %24.16g\n",errrel);

if (errrel < 1.0e-13)  {
  printf("===[ PASSED                ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                ]====================================================\n");
printf("\n");


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
            pols[n] = 2.0*x*pols[n-1] - pols[n-2];
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
   val = 0.0;
   x_2 = 2.0*y/(b-a)- (b+a)/(b-a);
   chebpolys(N,x_2,pols);
   for(n = 0; n <= N; n++ ){
        if (n ==0){
        val = val + coefs[n]*(0.5)*pols[n];
        }
        else{
        val = val + coefs[n]*pols[n];
        }
   }
   return val;
}
