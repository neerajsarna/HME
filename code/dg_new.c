/*The code performs the simulation for a shock tube considering the fluid behaviour to be governed by Euler equations*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>

#define NDEBUG		//uncomment the code when not using assert
#include <assert.h>
#include "time.h"


#define p 0	 				//degree of the polynomial
#define var 5				//total convective variables, set for euler equations + stiffened as law
#define ngp 4				// total number of gauss points being used in the code
#define gamma 5/3.0

/*
fluid-1 = fluid on the left of the initial shock
fluid-2 = fluid on the right of the initial shock
*/

double w[ngp];				//values of the quadrature
double x_quad[ngp];			//location of the quadrature on a 0 to 1 scale

int N;						//total number of intervals between 0 and 1
int ngc;					//total number of ghost cells
double xl = -2;					//left hand boundary
double xr = 2;					//right hand boundary
double Kn = 0.5;				// the Knudsen number

#include "tools.h"
#include "evauate_quad.h"		// evalutes the quadrature location and the weights
#include "quad_values.h"
#include "set_boundary.h"
#include "set_initial.h"
#include "pres_sound.h"			// evaluation of the speed of sound
#include "develop_G.h"			// develops the non-conservative matrix G
#include "defs.h"				//defining the legendre polynomials and fluxes
#include "fluxes_cons.h"		//flux evaluation
#include "quad_cons.h"			//gaussian quadrature
#include "leg.h"				//development of legendre polynomials
#include "grid_initial.h"		//generates the points on the grid and sets the initial condition
#include "shu.h"				// shu limiter to be used 
//#include "limiter.h"				// file for characteristic limiting 
//#include "limiter_multispecies.h"	// limiter for the multispecies
//#include "limiter_moments.h"		// limiter for the moment equations
#include "limiter_moments_pressure.h"
#include "production_term.h"
/*The following two limiters are available
1. limiter_primitive: provides limiting on the primitive variables
2. limiter_euler: provides limiting on the characteristic variables
*/

#include "time_stepping.h"					//evaluates the maximum propogation speed depending upon the average
#include "fluxes_noncons.h" 				//evalutes the flux part coming from the non conservative part
#include "quad_noncons.h"
#include "update.h"							// performs the time update of the solution

main(int agrc, char **argv)
{
 legendre_compute_glr ( ngp,x_quad, w );

 N = atoi(argv[1]);

 ngc = 2;									//set for the sod's problem

 
 double time_factor = 2;					//adjusts the time stepping for stabilization
 double delta_x,delta_t;					//time and space discretization
 int i,j,k;	
 double CFL = 0.3;
 int limiting = 1;							//sets the limiting


 /*
 List of test cases:
 0. sods shock tube problem
 1. periodic data for limiter testing
 */

 int test_case = 0;				

 char opts_out[30] = "mean";

/****memory allocation**********/
 
 polynomial l;							//contains all the properties of the legendre polynomial
 space grid;							//contains the midpoint and the interfaces of each cell
 data solution;							//contains the current coefficients, the temporary coefficients and the fluxes 
 stiff noncons;

		

 delta_x = (xr-xl)/(N*1.0);				//grid spacing


 leg(l.legendre,l.der_legendre);					//develops the legendre polynomials, ref. leg.h
 quadrature_values(l.legendre,l.der_legendre,l.values,l.der_values);	//evaluates the value of the legendre at the quadrature points,  ref. leg.h
 boundary_values(l.legendre,l.b_values);				//evaluates the value of the legendre at the boundaries {0,1}, ref. leg.h

 #include "mem_allocate.h"
 
 //setting of the initial conditions and the grid points
 set_grid(grid.x,grid.x_mid,delta_x,xl);					//evaluates the location of the mid points and the cell interfaces, ref. gridinitial.h

 
 // setting of the initial conditions 
 set_initial(solution.coeff,l.values,delta_x,grid.x_mid,grid.x,test_case);		//evaluates the coefficients for t=0, ref.gridinitial.h

 // setting of the boundary conditions
 set_boundary(solution.coeff);

 /*********testing of the limiter ****/
 /*set_initial(solution.coeff,l.values,delta_x,grid.x_mid,grid.x,test_case);	
 printf("LIMITING PRESSURE \n");
 limiter_moment(solution.coeff,l.b_values);
 set_initial(solution.coeff,l.values,delta_x,grid.x_mid,grid.x,test_case);	
 printf("LIMITING THETA\n");
 limiter_moment_prim1(solution.coeff,l.b_values);*/
 /********* end of testing of the limiter ********/
 

 double t=0;
 double T= 0.3;
 int counter = 0;		//counts the number of iterations
 unsigned int flag = 0 ; 

 //update of the solution
 while (t < T)
	{

		delta_t = CFL*delta_x/(time_factor*time_stepping(solution.coeff));		

		if (t + delta_t > T && flag == 0)
		{
			delta_t = T - t;
			flag = 1;    
		}

		t = t + delta_t;


		double tau;

		tau = delta_t/2;
		time_update(&solution,&noncons,l,delta_x,tau);
		
		tau = delta_t;
		production_term(tau,solution.coeff,l.values);

		tau = delta_t / 2;
		time_update(&solution,&noncons,l,delta_x,tau);

		counter++;
		if (counter % 100 ==0 )
			printf("time :%f\n",t);

	}

 //plotting of the solution obtained
 output(solution.coeff,grid.x,l.b_values);			//ref: grid_initial.h



 #include "mem_free.h"

}

   
