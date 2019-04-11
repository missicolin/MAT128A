#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void traprule(unsigned int n, double a, double b, double *xs, double *whts)

/*
 *  Return the nodes
 *
 *    x0, x1, ..., xn
 *
 *  and weights
 *
 *    w0, w1, ..., wn
 *
 *  of the (n+1)-point trapezoidal rule on the interval [a,b].  Storage for the
 *  output arrays xs and whts is supplied by the user.  They must both be of
 *  length greater than or equal to the size of (n+1) double precision numbers.
 *
 *  Input parameters:
 *    n - a positive integer specifying the number of quadrature nodes in the
 *        desired rule
 *    a - the left-hand endpoint of the interval on which the quadrature rule
 *        is defined
 *    b - the right-hand endpoint of the interval on which the quadrature rule
 *        is defined
 *
 *  Output parameters:
 *    xs - an array whose j^th  element xs[j] will be the value of the j^{th}
 *         quadrature node x_j
 *    whts - an array whose j^th element whts[j] will be the value of j^{th}
 *         quadrature weight
 *
 *
 */

{
    unsigned int j;
    int lambda;
    double h = (b-a)/n;
    for (j = 0; j < n+1; j++ ){
        xs[j] = a + h*j;
        if (j == 0 || j == n){
            lambda = 1;
            whts[j] = (h/2.0)*(lambda);
        }
        else{
            lambda = 2;
            whts[j] = (h/2.0)*(lambda);
        }
    }
}


