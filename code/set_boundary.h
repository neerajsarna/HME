set_boundary(double **coeff)
{

	for (unsigned int i = 0 ; i < ngc/2 ; i++)									
		for (unsigned int j = 0 ; j < var*(p+1) ; j++)
			{
					coeff[ngc/2-1-i][j] = coeff[ngc/2 + i][j];			// value from first cell in the domain
					coeff[N+ngc/2+i][j] = coeff[N+ngc/2-1-i][j];		// value from the last cell in the domain
			}
	

}

