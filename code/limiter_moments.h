// evaluatio of the value of the solution to the right boundary
void boundary_right(double *coeffs,double *result,double (*b_values)[2])
{
	for (unsigned int i = 0 ; i < var ; i ++)
	{
		result[i]  = 0;
		for (unsigned int j = 0 ; j < p + 1 ; j++)
			result[i] += coeffs[j + i * (p + 1 )] * b_values[j][1];
	}


}

// evaluation of the value of the solution to the right boundary
void boundary_left(double *coeffs,double *result,double (*b_values)[2])
{
	for (unsigned int i = 0 ; i < var ; i ++)
	{
		result[i]  = 0;
		for (unsigned int j = 0 ; j < p + 1 ; j++)
			result[i] += coeffs[j + i * (p + 1 )] * b_values[j][0];
	}


}


// limiting of the solution using shu limiter
void limiter_moment_prim1(double **coeff,double (*b_values)[2])
{

	double mean_middle, mean_right, mean_left;	// mean values of the solution
	double *right_middle, *left_middle;			// value of the solution to the right and the left
	double modified_high_order_left,modified_high_order_right;	//value obtained from the limiting process
	bool flag = false;

	right_middle = calloc(var,sizeof(double));
	left_middle = calloc(var,sizeof(double));

	for (unsigned int i = 0 ; i < N ; i ++)
	{
		flag = false;

		boundary_left(coeff[i + ngc/2],left_middle,b_values);					// value of the solution to the left boundary
		boundary_right(coeff[i + ngc/2],right_middle,b_values);					// value of the solution to the right boudnary
	
		for (unsigned int j = 0 ; j < var ; j ++)
			{
				mean_middle = coeff[i + ngc/2][j * ( p + 1)];
				mean_left = coeff[i + ngc/2 - 1][j * ( p + 1 )];
				mean_right = coeff[i + ngc/2 + 1][j * ( p + 1 )];

				modified_high_order_right = shu(right_middle[j]-mean_middle,mean_right-mean_middle,mean_middle-mean_left); 
				modified_high_order_left = shu (mean_middle-left_middle[j],mean_right-mean_middle,mean_middle-mean_left);
			

				if (j == 0)
				{
					printf("slopeMR %f slopeHR %f corrected_slopeR %f\n",
							 mean_right-mean_middle,
							 right_middle[j]-mean_middle,
							 modified_high_order_right);
					
					printf("slopeML %f slopeHL %f corrected_slopeL %f\n",
							 mean_middle-mean_left,
							 mean_middle-left_middle[j],
							 modified_high_order_left);
										
				}

				coeff[i + ngc/2][1 + j * (p + 1)] = (modified_high_order_left
													+modified_high_order_right)/(2 * sqrt(3));


/*				if ((fabs(modified_high_order_left - (mean_middle-left_middle[j])) > 1e-5)
						 || (fabs(modified_high_order_right - (right_middle[j]-mean_middle)) > 1e-5))
				{
					flag = true;
					break;
				}*/

			}

		if (flag)
		for (unsigned int j = 0 ; j < var ; j ++)
			for (unsigned int k = 1 ; k < p + 1 ; k++)
				coeff[i + ngc/2][k + j * (p + 1)] = 0;


	}

		
}
