#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "evauate_quad.h"

main()
{
	double a,b;			// a = left end point, b = right end point
	int n;				// n = total number of gauss points
	double *w,*x;		// w = weights , x = locations
	int i;				// loop counter

	a = -1;
	b = 1;
	n = 500;

	w = calloc(n,sizeof(double));			// w = weights for the gauss points
	x = calloc(n,sizeof(double));			// x = location of the gauss points

	legendre_compute_glr ( n, x, w );


	for (i = 0 ; i < n ; i ++)
		printf("weights:%f location:%f\n",w[i],x[i]);
}