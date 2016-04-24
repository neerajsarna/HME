/****evaluates the value of the solution at the four quadrature points (No ghost cells)****/
void values_quadrature(double **coeff, double (*values)[ngp],double **sol_quad)
{
	int i,j,k,m;

	for (i = 0 ; i < N ; i++)							//loop over all interior cells
	   for (m = 0 ; m < var ; m++)						//loop over all available variables
		for (j = 0 ; j < ngp ; j++)						//loop over all quadrature points
			for (k = 0 ; k < p+1 ;k++)					//loop over all legendre polynomials
				sol_quad[ngp*i + j][m] += coeff[i+ ngc/2][k+m*(p+1)] * values[k][j];


}

// the following routine calculates the value of the derivative of the solution at the four quadrature points
void der_values_quad(double **coeff,double (*der_values)[ngp],double **der_sol_quad)
{
	int i,j,k,m;

	for (i = 0 ; i < N ; i++)							//loop over all interior cells
	   for (m = 0 ; m < var ; m++)
		for (j = 0 ; j < ngp ; j++)						//loop over all quadrature points
			for (k = 0 ; k < p+1 ;k++)					//loop over all legendre polynomials
				der_sol_quad[ngp*i + j][m] += coeff[i+ ngc/2][k+m*(p+1)]*der_values[k][j];
}


// the following routine evaluates the value obtained of the quadrature obtained from the non conservative part
void quad_noncons(double **coeff,double **quad,double (*values)[ngp],double (*der_values)[ngp])
{
	int i,j,k,l;
	double *result;						//stores the result of the matrix vector multiplication
	double **sol_quad;					//value of the solution at the quadrature points
	double **der_sol_quad;				//value of the derivative of solution at the quadrature points
	double **G;							//value of nonconservative matrix
	double delta_x = (xr-xl)/N;
	
	sol_quad = calloc(ngp*N,sizeof(double));				//ghost cells have been excluded, for every interior cell four quadrature points
	der_sol_quad = calloc(ngp*N,sizeof(double));			
	result = calloc(var,sizeof(double));
	G = calloc(var,sizeof(double));


	for(i = 0 ; i < ngp * N ; i ++)
	{
		sol_quad[i] = calloc(var,sizeof(double));
		der_sol_quad[i] = calloc(var,sizeof(double));
	}

	for(i = 0 ; i < var ; i ++)
		G[i] = calloc(var,sizeof(double));

	values_quadrature(coeff,values,sol_quad);				//evaluates the value of the solution at the quadrature points
	der_values_quad(coeff,der_values,der_sol_quad);			//evaluates the derivative of the solution at the quadrature points

	for(i = 0 ; i < N ; i ++)
		for(k = 0 ; k < p + 1 ; k ++)
			for(l = 0 ; l < var ; l++)
				quad[i][k + l*(p+1)] = 0;


	for(i = 0 ; i < N ; i ++)
	   for(k = 0 ; k < p + 1 ; k ++)
		for(j = 0 ; j < ngp ; j ++)
		{
			develop_G(G,sol_quad[ngp * i + j]);
			matrix_multiply(G,der_sol_quad[ngp * i + j],result, var, var);	

			for(l = 0 ; l < var ; l ++)
				quad[i][k + l*(p+1)] += 0.5 * w[j] * result[l] * values[k][j]/delta_x;

		}



// freeing the 1-D pointers
for(i = 0 ; i < var ; i ++)
	free(G[i]);

free(result);	

for(i = 0 ; i < ngp * N ; i ++)
{
	free(sol_quad[i]);
	free(der_sol_quad[i]);
}

// freeing the 2-D pointers
free(sol_quad);
free(der_sol_quad);
free(G);

}