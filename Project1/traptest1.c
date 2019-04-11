#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Declaration of the external traprule function */
void traprule(unsigned int n, double a, double b, double *xs, double *whts);


/* The test functions */
double fun1(double x) {
return exp(sin(x)-x*x);
}

double fun2(double x) {
return exp(cos(x));
}

double fun3(double x) {
return x*x*x-x+1.0;
}

double fun4(double x) {
return cos(10.0*x)*cos(10.0*x);
}

double trapint(unsigned int n,double a, double b, double (*fun)(double))

/*  Use the (n+1) trapezoidal rule to approximate the integral of a function
 *  f(x) on the interval (a,b).
 *
 *  Input parameters:
 *    n - indicates that the (n+1)-point trapezoidal rule should be used
 *    (a,b) - the interval which is the domain of integration
 *    fun - a pointer to the function which evaluates the integrand
 *
 *  Output parameters:
 *    None
 *
 *  Return value:
 *    the obtained approximation of the integral
 */

{
  double *xs, *whts;
  double val;
  unsigned int i;

  xs   = (double *)malloc( sizeof(double) * (n+1) );
  whts = (double *)malloc( sizeof(double) * (n+1) );

  if (xs == NULL || whts == NULL) {
    printf("int trapint, memory allocation failed!\n");
    exit(-1);
  }

  traprule(n,a,b,xs,whts);

  val = 0;
  for(i=0;i<=n;i++)
    val = val + whts[i]*fun(xs[i]);

  free((void*)xs);
  free((void*)whts);


  return val;
}



double trapadapint(double eps, double a, double b,double (*fun)(double))
/*
 *  Adaptively integrate a user-supplied function over the interval (a,b)
 *  using trapezoidal rules of varying lengths.
 *
 *  Input parameters:
 *    eps - the desired precision for the result
 *    (a,b) - the integration domain
 *    fun - pointer to a user-supplied function which returns the values
 *      of the integrand
 *
 *  Output parameters:
 *    None
 *
 *  Return value:
 *    The obtained approximation of the integral
 *
 */
{
  unsigned int i;
  unsigned int n    = 16;
  unsigned int maxn = 2<<20;
  double val, val2, errrel;

  val = trapint(n,a,b,fun);

  do {
  n    = n*2;
  if (n > maxn) {
    printf("integrate: maximum quadrature size exceeded");
    exit(-1);
  }

  val2   = trapint(n,a,b,fun);
  errrel = fabs(val2-val)/ (fabs(val)+1.0);
  val     = val2;

//  printf("n = %8d    val  = %24.16f   val2  = %24.16f derr =   %24.16f \n",n,val,val2,errrel);


  } while (errrel > eps);

  exit(0);
}





int main(int argc, char **argv)
{
  double a,b,val,val0,errrel,pi;
  unsigned int n,i,nscore;

  pi = acos(-1.0);
  nscore =  0;

  n      =  10;
  a      =  0.0;
  b      =  2*pi;
  val    =  trapint(n,a,b,fun1);
  val0   =  1.49447794313489562752224044227326427;
  errrel =  fabs(val-val0)/fabs(val0);
  printf("\n");
  printf("-------------------------------------------------\n");
  printf("f(x) = exp(sin(x)-x^2)\n");
  printf("-------------------------------------------------\n");
  printf("    n  = %24d\n",n);
  printf("    a  = %24.15f\n",a);
  printf("    b  = %24.15f\n",b);
  printf("  val  = %24.16f\n",val);
  printf("  val0 = %24.16f\n",val0);
  printf("  derr = %24.16f\n",errrel);
  printf("-------------------------------------------------\n");

  if (errrel < .05) {
    nscore = nscore+1;
    printf("TEST SUCCEEDED -- RELATIVE ERROR < 0.05\n");
  }else {
    printf("TEST FAILED    -- RELATIVE ERROR > 0.05\n");
  }
  printf("-------------------------------------------------\n");
  printf("\n");

  n      =  100;
  a      =  0.0;
  b      =  2*pi;
  val    =  trapint(n,a,b,fun1);
  val0   =  1.49447794313489562752224044227326427;
  errrel =  fabs(val-val0)/fabs(val0);
  printf("\n");
  printf("-------------------------------------------------\n");
  printf("f(x) = exp(sin(x)-x^2)\n");
  printf("-------------------------------------------------\n");
  printf("    n  = %24d\n",n);
  printf("    a  = %24.15f\n",a);
  printf("    b  = %24.15f\n",b);
  printf("  val  = %24.16f\n",val);
  printf("  val0 = %24.16f\n",val0);
  printf("  derr = %24.16f\n",errrel);
  printf("-------------------------------------------------\n");

  if (errrel < .0005) {
    nscore = nscore+1;
    printf("TEST SUCCEEDED -- RELATIVE ERROR < 0.0005\n");
  }else {
    printf("TEST FAILED    -- RELATIVE ERROR > 0.0005\n");
  }
  printf("-------------------------------------------------\n");
  printf("\n");


  n      = 16;
  a      = 0.0;
  b      = 2*pi;
  val    =  trapint(n,a,b,fun2);
  val0   =  7.95492652101284527451321966532939433;
  errrel =  fabs(val-val0)/fabs(val0);
  printf("\n");
  printf("-------------------------------------------------\n");
  printf("f(x) = exp(cos(x))\n");
  printf("-------------------------------------------------\n");
  printf("    n  = %24d\n",n);
  printf("    a  = %24.15f\n",a);
  printf("    b  = %24.15f\n",b);
  printf("  val  = %24.16f\n",val);
  printf("  val0 = %24.16f\n",val0);
  printf("  derr = %24.16f\n",errrel);
  printf("-------------------------------------------------\n");

  if (errrel < 1.0e-12) {
    nscore = nscore+1;
    printf("TEST SUCCEEDED -- RELATIVE ERROR < 1.0e-12\n");
  }else {
    printf("TEST FAILED    -- RELATIVE ERROR > 1.0e-12\n");
  }
  printf("-------------------------------------------------\n");
  printf("\n");

  n      =  10000;
  a      =  0.0;
  b      =  1.0;
  val    =  trapint(n,a,b,fun3);
  val0   =  0.750000000000000000000000000000000000;
  errrel =  fabs(val-val0)/fabs(val0);
  printf("\n");
  printf("-------------------------------------------------\n");
  printf("f(x) = x^3 - x + 1\n");
  printf("-------------------------------------------------\n");
  printf("    n  = %24d\n",n);
  printf("    a  = %24.15f\n",a);
  printf("    b  = %24.15f\n",b);
  printf("  val  = %24.16f\n",val);
  printf("  val0 = %24.16f\n",val0);
  printf("  derr = %24.16f\n",errrel);
  printf("-------------------------------------------------\n");

  if (errrel < 1.0e-7) {
    nscore = nscore+1;
    printf("TEST SUCCEEDED -- RELATIVE ERROR < 1.0e-7\n");
  }else {
    printf("TEST FAILED    -- RELATIVE ERROR > 1.0e-7\n");
  }
  printf("-------------------------------------------------\n");
  printf("\n");


  n      =  11;
  a      =  0.0;
  b      =  2*pi;
  val    =  trapint(n,a,b,fun4);
  val0   =  pi;
  errrel =  fabs(val-val0)/fabs(val0);
  printf("\n");
  printf("-------------------------------------------------\n");
  printf("f(x) = cos(x)\n");
  printf("-------------------------------------------------\n");
  printf("    n  = %24d\n",n);
  printf("    a  = %24.15f\n",a);
  printf("    b  = %24.15f\n",b);
  printf("  val  = %24.16f\n",val);
  printf("  val0 = %24.16f\n",val0);
  printf("  derr = %24.16f\n",errrel);
  printf("-------------------------------------------------\n");

  if (errrel < 1.0e-12) {
    nscore = nscore+1;
    printf("TEST SUCCEEDED -- RELATIVE ERROR < 1.0e-12\n");
  }else {
    printf("TEST FAILED    -- RELATIVE ERROR > 1.0e-12\n");
  }
  printf("-------------------------------------------------\n");
  printf("\n");


  printf("SCORE = %d / 5 \n",nscore);

  return 0;
}

