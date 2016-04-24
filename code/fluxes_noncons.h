// evaluates a matrix vector dot product
// m = number of rows
// n = number of coloums
// has been checked
void matrix_multiply(const double **G,const double *path_der,double *result,const  int m,const int n)	
{
	int i,j;
	double sum = 0;

	for(i = 0 ; i < m ; i ++ )			//loop over rows
	{
		sum = 0;
		for(j = 0 ; j < n ; j ++)		//loop over coloumns
			sum = G[i][j] * path_der[j] + sum;

		result[i] = sum;
	}
		
		
}	

//evaluates the value of the path and the derivative of the path wrt tau
// u right = value at the right of the interface
//u left = value at the left of the interface

void path_evaluate(double *uright,double *uleft,double *path,double *path_der,double tau)
{
	double alpha = 1;		//can be changed for testing purposes
	int i;

	for(i = 0 ; i < var ; i ++)
	{
		path[i] = (uright[i] - uleft[i]) * pow(tau,alpha) + uleft[i];
		path_der[i] = (uright[i] - uleft[i]) * alpha * pow(tau,alpha-1);
	}
}


// evaluates the contribution to the fluxes coming from the non conservative kind
void fluxes_noncons(double **coeff,double **qr,double **ql,double (*b_values)[2])
{

	double *path;				//value of the path at different tau
	double *path_der;			//evaluates the derivative of the path with respect to tau
	double **G;					//evaluates the matrix of the non-conservative part
	double **uleft;				//value of solution at x-1/2
	double **uright;				//value of solution at x+1/2
	double *result;				//result obtained after multiplication of path and G
	double *sum;				//result obtained after integration
	double loc_quad;			//value of the quadrature point on a zero to one scale
	int i,j,k;					//loop index

	uleft = calloc(N+ngc,sizeof(double));			//evaluation at all the cells			
	uright = calloc(N+ngc,sizeof(double));
	path = calloc(var,sizeof(double));
	G = calloc(var,sizeof(double));
	sum = calloc(var,sizeof(double));
	

	for (i = 0 ; i < N + ngc ; i ++)
	{
		uleft[i] = calloc(var,sizeof(double));
		uright[i] = calloc(var,sizeof(double));
	}

	for(i = 0 ; i < var ; i ++)
		G[i] = calloc(var,sizeof(double));

	path = calloc(var,sizeof(double));
	path_der = calloc(var,sizeof(double));
	result = calloc(var,sizeof(double));
	

	//evaluation of the solution and the fluxes at the cell interfaces for every cell
	for (i = 0 ;i < N + ngc ; i++)						//loop over each cell
		for (j = 0 ; j < var ; j++)						//loop over each convective variable		
		  for (k = 0 ; k < p+1 ; k++)					//loop over all the legendre polynomials
				{
					uleft[i][j] += b_values[k][0]*coeff[i][k + j*(p+1)];		//x-1/2	
					uright[i][j] += + b_values[k][1]*coeff[i][k + j*(p+1)];		//x+1/2
				}


	for(i = 0 ; i < N + 1; i ++)			// loop over all the interior interfaces
	{
		for(k = 0 ; k < var ; k ++)
			sum[k] = 0;

		for( j = 0 ; j < ngp ; j ++)			// loop over the quadrature points
		{
			loc_quad = (x_quad[j] + 1)/2;				//checked, mapped the quadrature points from zero to one
			
			path_evaluate(uleft[ngc/2 + i],uright[ngc/2-1 + i],path,path_der,loc_quad);	//evaluation of derivative of path and path derivative at gauss points


			develop_G(G,path);						// develop the G matrix for the sytem
			
			
			matrix_multiply(G,path_der,result,var,var);

			for(k = 0 ; k < var ; k ++)		// loop over the number of variables
				sum[k] = 0.5 * w[j] *result[k] + sum[k];			//flux to the right of the interface
			
		}

		for(k = 0 ; k < var ; k ++)		// loop over the number of variables
				qr[i][k] = sum[k];			//flux to the right of the interface


	}

	

	free(path);
	free(path_der);
	free(result);
	free(sum);

	for(i = 0 ; i < var ; i ++)
		free(G[i]);
	
	for (i = 0 ; i < N + ngc ; i++)
	{
		free(uleft[i]);
		free(uright[i]);
	}	

	free(G);
	free(uleft);
	free(uright);

}