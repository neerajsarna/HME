// evaluates the coefficients of the solution from the solution at the quadrature point
void inv_values_quad(double **u_quad,double *coeff,double (*values)[ngp])
{
	int i,j,k;

	empty(coeff,var * (p+1));


	for (i = 0 ; i < p + 1 ; i++)
		for (j = 0 ; j < var ; j++)
			for (k = 0 ; k < ngp ; k++)
				coeff[i + j*(p+1)] += 0.5 * w[k] 
									  * u_quad[k][j] * values[i][k]; 


	return;
}

// the following function evaluates the value of the solution at the quadrature points 
void values_quad_per_cell(double *coeff,double (*values)[ngp],double **sol_quad)
{
	empty_2d(sol_quad,ngp,var);

	for (unsigned int i = 0 ; i < ngp ; i ++)
		for (unsigned int j = 0 ; j < var ; j ++)
			for (unsigned int k = 0 ; k < p + 1 ; k ++)
			sol_quad[i][j] += coeff[k + j * (p + 1)] * values[k][i];

}
/****evaluates the value of the solution at the four quadrature points (No ghost cells)****/
void values_quad(double **coeff, double (*values)[ngp],double **sol_quad)
{
	int i,j,k,m;
	empty_2d(sol_quad,ngp,var);

	for (i = 0 ; i < N ; i++)							//loop over all interior cells
	   for (m = 0 ; m < var ; m++)
		for (j = 0 ; j < ngp ; j++)						//loop over all quadrature points
			for (k = 0 ; k < p+1 ;k++)					//loop over all legendre polynomials
				sol_quad[i][j + m * ngp] += coeff[i+ ngc/2][k+m*(p+1)]*values[k][j];


}