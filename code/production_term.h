void production_term(const double delta_t,double **coeff,(*values)[ngp])
{
	double **u_quad;						// store the value of the solution at the quadrature points

	u_quad = calloc(ngp,sizeof(double));

	for (unsigned int i = 0 ; i < ngp ; i ++)
		u_quad[i] = calloc(var,sizeof(double));

	for (unsigned int i = 0 ; i < N ; i ++) // loop over all the cells
	{
		values_quad_per_cell(coeff[i + ngc/2],values,u_quad);

		for (unsigned int j = 0 ; j < ngp ; j ++)
			for (unsigned int k = 3 ; k < var ; k ++)
			{
				double rho = u_quad[j][0];
				u_quad[j][k] *= exp(-rho * delta_t / Kn);
			}

		inv_values_quad(u_quad,coeff[i+ngc/2],values);

	}

    set_boundary(coeff);
    
    for(unsigned int i = 0 ; i < ngp ; i ++)
    	free(u_quad[i]);

    free(u_quad);
}