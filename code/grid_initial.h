
/***sets the points on the grid*******/
set_grid(double *x,double *x_mid,double delta_x, double xl)
{
  int i;
  for (i = 0 ; i < N+1+ngc ; i ++)
	{
	x[i] = xl + (i-ngc/2)*delta_x;		//considering equal number of gauss points at the start and end of the domain

	if (i < N+ngc)				//total number of mid points are 1 less than total number of interfaces
	x_mid[i] = xl + ( i-ngc/2 + 0.5) * delta_x;
	}
}



void output(double **coeff,double *x,double (*b_values)[2])	
{

	int i,j,k,l;

	unsigned int sample_point = 3;		// points where the solution is needed, only the points on the cell and quadrature have been implemented
	double **u;							// value of the solution at the sample points
	double *location;					// location of the sample points on the global space
	double delta_x = (xr-xl)/N;			//spacing in the x-direction
	FILE *fp;
	FILE *fp2;							// store the location of the points on the grid

	fp = fopen("result_high_resolution/solutionKn0p5order1","w+");
	
	assert(fp != NULL);

	fp2 = fopen("result_high_resolution/locationKn0p5order1","w+");

	assert(fp2 != NULL);

	u = calloc(sample_point * (N + ngc),sizeof(double));
	location = calloc(sample_point * (N+ngc),sizeof(double));

	for(i = 0 ; i < sample_point * (N + ngc) ; i ++)
		u[i] = calloc(var,sizeof(double));

	for (i = 0 ; i < N + ngc ; i ++)								//loop over cells
		  for(l = 0 ; l < var ; l++)								// loop over variables
		  {
		  	for (k = 0 ; k < p +1 ; k ++)							// loop over total polynomials
		  	{
		  		u[3 * i][l] += coeff[i][k + l*(p+1)]*b_values[k][0];
				u[3 * i + 2][l] += coeff[i][k + l*(p+1)]*b_values[k][1];
		  	}
		  		u[3 * i + 1][l] = coeff[i][l*(p+1)];
		  }
			
  for(i = 0 ; i < N + ngc ; i ++)
  	{
  		location[3* i] = x[i];

  		location[3 * i + 1] = x[i] + delta_x/2;

  		location[3 * i + 2] = x[i] + delta_x; 
  	}


  fprintf(fp, "# all the primitive variables\n");
  fprintf(fp2,"#location\n");

  for (i = 0 ; i < (N + ngc) ; i ++)
  {
  	fprintf(fp2, "%f\n",location[3 * i]);

  	fprintf(fp2, "%f\n",location[3 * i + 2]);

    for (l = 0 ; l < var ; l++)
  		fprintf(fp, " %f ",u[3 * i][l]);

  	fprintf(fp, "\n");

  	for (l = 0 ; l < var ; l++)
  	  fprintf(fp, " %f ",u[3 * i + 2][l]);

  	fprintf(fp, "\n" );

    
  }	

  free(location);
  fclose(fp2);
  fclose(fp);

  for (i = 0 ; i < sample_point * (N + ngc) ; i++)
  	free(u[i]);

  free(u);
  

}









