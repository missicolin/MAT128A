#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void fourcoefs(int N, double complex (*funptr)(double), double complex *coefs);
double complex foureval(int N, double complex *coefs, double t);


double relative_error(int N, double complex *coefs, double complex *coefs0)
{
  double derr,dnorm;
  int i;
  derr = 0.0;
  dnorm = 0.0;

  for(i=0; i<= 2*N; i++) {
    derr  = derr + cabs(coefs[i]-coefs0[i])*cabs(coefs[i]-coefs0[i]);
    dnorm = dnorm + cabs(coefs0[i])*cabs(coefs0[i]);
  }
  derr = sqrt(derr / dnorm);
  return derr;  
}

void print_coefs(int N, double complex *coefs)
{
  int i;
  for(i=0; i <= 2*N; i++) {
     printf("a_{%4.3d} = (%24.16f,%24.16f)\n",i-N,creal(coefs[i]),cimag(coefs[i]));
  }
}


double complex fun0(double t) 
{
  return (double complex)cos(t);
}

double complex fun1(double t) 
{
  return 3.0 / (5.0 - 4*cos(t)) ;
}

double complex fun2(double t) 
{
  return 2.0 / (2-cexp(I*t));
}

double complex fun3(double t)
{
  return cexp(pow(fabs(sin(t)),3));
}

double complex fun4(double t)
{
  return cexp(pow(fabs(sin(t)),9));
}



int main(int argc, char **argv)
{
int nscore;
int N, i;
double complex *coefs, *coefs0, val, val0;
double t,errrel;

nscore = 0;


/* Approximate the Fourier coefficients of f(t) = cos(t) with  N = 4  */

N      = 4;

coefs  = (double complex *)malloc(sizeof(double complex)*(2*N+1));
coefs0 = (double complex *)malloc(sizeof(double complex)*(2*N+1));


for(i=0; i<= 2*N; i++) {
  coefs0[i] = 0;
}
coefs0[N-1]=0.5;
coefs0[N+1]=0.5;

fourcoefs(N,fun0,coefs);
val = foureval(N,coefs,t);

printf("\n");
printf("===[ f(x) = cos(t)           ]====================================================\n");
errrel = relative_error(N,coefs,coefs0);

printf("N                              = %7.4d\n",N);
printf("relative error in coefficients = %24.16g\n",errrel);

if (errrel < 1.0e-13)  {
  printf("===[ PASSED                ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                ]====================================================\n");
printf("\n");

//print_coefs(N,coefs0);
//print_coefs(N,coefs);

free(coefs);
free(coefs0);



/* Approximate the Fourier coefficients of f(t) = 3/(5-4*cos(t)) with  N = 20  */

N      = 20;

coefs  = (double complex *)malloc(sizeof(double complex)*(2*N+1));
coefs0 = (double complex *)malloc(sizeof(double complex)*(2*N+1));

for(i=0; i<= 2*N; i++) {
  coefs0[i] = pow(0.5,abs(i-N));
}

fourcoefs(N,fun1,coefs);
val = foureval(N,coefs,t);

printf("\n");
printf("===[ f(x) = 3/(5-4*cos(t)) ]====================================================\n");
errrel = relative_error(N,coefs,coefs0);

printf("N                              = %7.4d\n",N);
printf("relative error in coefficients = %24.16g\n",errrel);

if (errrel < 2.0e-6)  {
  printf("===[ PASSED                ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                ]====================================================\n");
printf("\n");

//print_coefs(N,coefs0);
//print_coefs(N,coefs);

free(coefs);
free(coefs0);



/* Approximate the Fourier coefficients of f(t) = 2/(2-exp(I*t)), N = 60 */

N      = 60;

coefs  = (double complex *)malloc(sizeof(double complex)*(2*N+1));
coefs0 = (double complex *)malloc(sizeof(double complex)*(2*N+1));

for(i=0; i< N; i++) {
  coefs0[i] = 0.0;
}

for(i=N; i<= 2*N; i++) {
  coefs0[i] = pow(0.5,i-N);
}

fourcoefs(N,fun2,coefs);

printf("\n");
printf("===[ f(x) = 2/(2-exp(it))    ]====================================================\n");
errrel = relative_error(N,coefs,coefs0);

printf("N                              = %7.4d\n",N);
printf("relative error in coefficients = %24.16g\n",errrel);

if (errrel < 2.0e-14)  {
  printf("===[ PASSED                ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                ]====================================================\n");
printf("\n");

/* //print_coefs(N,coefs0); */
/* //print_coefs(N,coefs); */

free(coefs);
free(coefs0);


/* Evaluate the approximate series for  f(t) = exp(|sin(t)|^3) when t = 1 */

N      = 100;
t      = 1.0;
val0   = fun3(t);

coefs  = (double complex *)malloc(sizeof(double complex)*(2*N+1));

fourcoefs(N,fun3,coefs);
val = foureval(N,coefs,t);

printf("\n");
printf("===[ f(x) = exp(|sin(t)|^3)  ]====================================================\n");
errrel = cabs(val-val0)/cabs(val0);

printf("N                              = %7.4d\n",N);
printf("relative error in evaluation   = %24.16g\n",errrel);

if (errrel < 1.0e-7)  {
  printf("===[ PASSED                ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                ]====================================================\n");
printf("\n");

//print_coefs(N,coefs);

free(coefs);


/* Evaluate the approximate series for  f(t) = exp(|sin(t)|^9) when t = 1/2 */

N      = 100;
t      = 0.5;
val0   = fun4(t);

coefs  = (double complex *)malloc(sizeof(double complex)*(2*N+1));

fourcoefs(N,fun4,coefs);
val = foureval(N,coefs,t);

printf("\n");
printf("===[  f(x) = exp(|sin(t)|^9) ]====================================================\n");
errrel = cabs(val-val0)/cabs(val0);

printf("N                              = %7.4d\n",N);
printf("relative error in evaluation   = %24.16g\n",errrel);

if (errrel < 1.0e-13)  {
  printf("===[ PASSED                ]====================================================\n");
   nscore = nscore+1;
} else
  printf("===[ FAILED                ]====================================================\n");
printf("\n");

//print_coefs(N,coefs);

free(coefs);


printf("\n\n");
printf("SCORE = %d / %d\n",nscore,5);
return nscore;

}
