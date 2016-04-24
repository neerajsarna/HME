

//evaluates the maxium characteristic speed at any interfaces using the values from both the sides
/* the input argument to this function uleft, is the value of the solution u to the left of the interface i, i.e uright[i-1] and uright is the value of the solution u to the right of the interface i i.e uleft[i]*/
double max_char_speed(double *uleft, double *uright)
{
	int i;
	const double c_left,c_right; //speed of sound at left and right interface
	double lambda_left,lambda_right;	//characteristic speed to left and right
	const double max_root_hermite = 2.7;
	const double theta_left = uleft[2];
	const double theta_right = uright[2];

	lambda_left = fabs(uleft[1]) + 2.7 * sqrt(theta_left);
	lambda_right = fabs(uright[1]) + 2.7 * sqrt(theta_right);

	return(fmax(lambda_left,lambda_right));
	

}

// computes the original flux in the system
/* description of the input arguments 
1. u = array of the solution
2. f = array for the original flux value
*/
void actual_flux(double *u,double *f)
{
	for (unsigned int i = 0 ; i < var ; i ++)
		f[i] = 0;

    return;
}
//evaluates the numerical flux at each interface using local lax friedrich
/*Could be seen from the numerical flux that the uleft corresponding to first cell is not being used and the uright corresponding to the last cell is not being used*/
/*the value of uleft for first cell and the value of u right for the last cell will be evaulated depending upon the boundary conditions to be used. Instead of evaluating the convective variables, directly the fluxes have been evaulated*/

/*At the interface i, the left of the interface has i-1 cell and the right of the interface has i cell.*/



void numerical_flux(double **fleft,double **fright,double **uleft,double **uright,double **q)
{
	double beta;							//maximum of the characteristic speed

	// the following routine update the value of fluxes for all the interfaces
	// no distinction has been made regarding the test case: CAUTION.

	for (unsigned int i = 0 ; i < N + 1 ; i++)					
	{
		beta = max_char_speed(uright[i-1 + ngc/2],uleft[i + ngc/2]);		//maximum characteristic speed at booth the interfaces
		for (unsigned int j = 0 ; j < var ; j++)
			q[i][j] = 0.5*(fright[i-1 + ngc/2][j] + fleft[i + ngc/2][j]) - 0.5*beta*(uleft[i + ngc/2][j] - uright[i-1 + ngc/2][j]);
	}

}


//evaluates the conservative part of the flux and the solution value
void fluxes_cons(double **coeff,double (*legendre)[p+1],double (*b_values)[2],double **q,double delta_x) 
{
	double **uleft;		//value of the solution at the interface left to the cell,x-1/2	, including ghost cells
	double **uright;	//value of the solution at the interface right to the cell,x+1/2, including ghost cells
	double **fleft;		//actual flux at the interface left to the cell,x-1/2
	double **fright;	//actual flux at the interface right to the cell,x+1/2
	

/*For N cells we have total N+1 interfaces but only N x+1/2 interfaces and N x-1/2 interfaces*/

	uleft = calloc(N+ngc,sizeof(double));				
	uright = calloc(N+ngc,sizeof(double));
	fleft = calloc(N+ngc,sizeof(double));
	fright = calloc(N+ngc,sizeof(double));

	for (unsigned int i = 0 ; i < N + ngc ; i ++)
	{
		uleft[i] = calloc(var,sizeof(double));
		uright[i] = calloc(var,sizeof(double));
		fleft[i] = calloc(var,sizeof(double));
		fright[i] = calloc(var,sizeof(double));
	}
	

	//evaluation of the solution and the fluxes at the cell interfaces
	for (unsigned int i = 0 ;i < N + ngc ; i++)						//loop over each cell
		for (unsigned int j = 0 ; j < var ; j++)					//loop over each convective variable		
		  for (unsigned int k = 0 ; k < p+1 ; k++)					//loop over all the legendre polynomials
				{
					uleft[i][j] +=	b_values[k][0]*coeff[i][k + j*(p+1)];		//x-1/2	
					uright[i][j] += b_values[k][1]*coeff[i][k + j*(p+1)];		//x+1/2
				}


	for (unsigned int i = 0 ; i < N+ngc ; i++)
		{
			actual_flux(uleft[i],fleft[i]);
			actual_flux(uright[i],fright[i]);
		} 

	numerical_flux(fleft,fright,uleft,uright,q);


	for (unsigned int i = 0 ; i < N+ngc ; i++)
	{
		free(fleft[i]);
		free(fright[i]);
		free(uright[i]);
		free(uleft[i]);
	}

	free(fleft);
	free(fright);
	free(uleft);
	free(uright);	
}


		

