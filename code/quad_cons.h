
/***the following function evaluates the value of the fluxes at the quadrature points,No ghost cells**/
void flux_quad(double **sol_quad,double **quad_flux)
{

	for (unsigned int i = 0 ; i < N ; i++)							//loop over cells
	{
		for (unsigned int j = 0 ; j < ngp ; j++)						//loop over quadrature points
			{

				quad_flux[i][j] = 0;		
				quad_flux[i][j + 1*ngp] = 0;
				quad_flux[i][j + 2*ngp] = 0;
				quad_flux[i][j + 3*ngp] = 0;
	
			}	

	}

	
}

//the following function evaluates the value of the integral of flux*derivative of the quadrature,NO ghost cells included
double quad_cons(double **quadrature, double **coeff,double (*legendre)[p+1],double (*values)[ngp],double (*der_values)[ngp])
{

 double delta_x = (xr-xl)/(N*1.0);
 int i,j,k,m;
 double **sol_quad;
 double **quad_flux;

 quad_flux = calloc(N,sizeof(double));		//stores the value of the actual flux at the four quadrature points, no ghost cells
 sol_quad = calloc(N,sizeof(double));		//stores the value of the solution at the four quadrature points

 for (i = 0 ; i < N ; i++)
	{
		sol_quad[i] = calloc(var*ngp, sizeof(double));	//value of 3 variables at the four quadrature points
		quad_flux[i] = calloc(var*ngp,sizeof(double));	//vaulue of the 3 corresponding fluxes at the four quadrature points
	}


 values_quad(coeff,values,sol_quad);
 flux_quad(sol_quad,quad_flux);

 for (i = 0 ; i < N ; i ++)
		for (j = 0 ; j < var  ; j ++)
			for (k = 0 ; k < p+1 ; k ++)
				quadrature[i][k + j *(p+1)] = 0;



 
  for (i = 0 ; i < N; i++)					//for number of cells
	   for (k = 0 ; k < var ; k++)
		for (j = 0 ; j < p+1 ; j++)			//for the legendre	
			for (m = 0 ; m < ngp; m++)		// for number of quadrature points
			     quadrature[i][j + k*(p+1)] += 0.5*w[m]*der_values[j][m]*(quad_flux[i][m + k*ngp])/delta_x;	//delta_x comes in because of the transformation

 for (i = 0 ; i < N ; i++)
	{
		free(sol_quad[i]);
		free(quad_flux[i]);
	}

	free(sol_quad);
	free(quad_flux);
			

}


