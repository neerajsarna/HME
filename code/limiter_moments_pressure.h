// the following functions implement the same limiter as discussed in limiter_moments it just
// uses pressure instead of theta
// the primitve variables used in the equations
// CAUTION: only implemented for p = 1
void mean_prim1(const double *coeff,double *mean)
{

	for (unsigned int i = 0 ; i < var ; i ++)
		mean[i] = coeff[i * (p + 1)];
		
}

// the following function replace theta by pressure
void prim1_prim2(const double *prim1,double *prim2)
{
	for (unsigned int i = 0 ; i < var ; i ++ )
		prim2[i] = prim1[i];

	double density = prim2[0];
	double theta = prim2[2];
	double pressure = density * theta;

	prim2[2] = density * theta;				// pressure = theat * density
}

void prim2_prim1(const double *prim2,double *prim1)
{
	for (unsigned int i = 0 ; i < var ; i ++)
		prim1[i] = prim2[i];

	double density = prim1[0];
	double pressure = prim1[2];
	double theta = pressure / density;

	prim1[2] = theta;
}

// evaluation of the value of the solution to the right boundary
void boundary_right_prim1(double *coeffs,double *result,double (*b_values)[2])
{
	for (unsigned int i = 0 ; i < var ; i ++)
	{
		result[i]  = 0;
		for (unsigned int j = 0 ; j < p + 1 ; j++)
			result[i] += coeffs[j + i * (p + 1 )] * b_values[j][1];
	}


}

// evaluation of the value of the solution to the right boundary
void boundary_left_prim1(double *coeffs,double *result,double (*b_values)[2])
{
	for (unsigned int i = 0 ; i < var ; i ++)
	{
		result[i]  = 0;
		for (unsigned int j = 0 ; j < p + 1 ; j++)
			result[i] += coeffs[j + i * (p + 1 )] * b_values[j][0];
	}
}

void limiter_moment(double **coeff,double (*b_values)[2])
{
	double *meanM_p1, *meanL_p1, *meanR_p1;
	double *meanM_p2, *meanL_p2, *meanR_p2;

	double *right_p1, *left_p1;
	double *right_p2, *left_p2;

	double *limited_right_p1, *limited_left_p1;
	double *limited_right_p2, *limited_left_p2;

	bool flag = false;

	meanM_p1 = calloc(var,sizeof(double));
	meanL_p1 = calloc(var,sizeof(double));
	meanR_p1 = calloc(var,sizeof(double));

	meanM_p2 = calloc(var,sizeof(double));
	meanL_p2 = calloc(var,sizeof(double));
	meanR_p2 = calloc(var,sizeof(double));

	right_p1 = calloc(var,sizeof(double));
	left_p1 = calloc(var,sizeof(double));

	right_p2 = calloc(var,sizeof(double));
	left_p2 = calloc(var,sizeof(double));

	limited_right_p1 = calloc(var,sizeof(double));
	limited_left_p1 = calloc(var,sizeof(double));

	limited_right_p2 = calloc(var,sizeof(double));
	limited_left_p2 = calloc(var,sizeof(double));

	for (unsigned int i = 0 ; i < N ; i ++)
	{
		flag = false;
		mean_prim1(coeff[i+ngc/2],meanM_p1);
		mean_prim1(coeff[i + ngc/2 - 1],meanL_p1);
		mean_prim1(coeff[i + ngc/2 + 1],meanR_p1);

		boundary_left_prim1(coeff[i + ngc / 2], left_p1,b_values );
		boundary_right_prim1(coeff[i + ngc / 2] , right_p1, b_values );

		prim1_prim2(meanM_p1,meanM_p2);
		prim1_prim2(meanL_p1,meanL_p2);
		prim1_prim2(meanR_p1,meanR_p2);

		prim1_prim2(left_p1,left_p2);
		prim1_prim2(right_p1,right_p2);

		for (unsigned int j = 0 ; j < var ; j ++)
		{
			limited_right_p2[j] = shu(right_p2[j]-meanM_p2[j],meanR_p2[j] - meanM_p2[j],
								meanM_p2[j]-meanL_p2[j]); 


			limited_left_p2[j] = shu(meanM_p2[j]- left_p2[j],meanR_p2[j] - meanM_p2[j],
								meanM_p2[j]-meanL_p2[j]); 

			if (fabs(limited_right_p2[j]-(right_p2[j]-meanM_p2[j])) > 1e-5
				|| fabs(limited_left_p2[j]- (meanM_p2[j] - left_p2[j])) > 1e-5)
			{
				flag = true;
				break;
			}

		/*	// converting back to values inside the edge
			limited_right_p2[j] = limited_right_p2[j] + meanM_p2[j];
			limited_left_p2[j] = -limited_left_p2[j] + meanM_p2[j];*/
		}

		/*prim2_prim1(limited_right_p2, limited_right_p1);
		prim2_prim1(limited_left_p2, limited_left_p1);

		for (unsigned int j = 0 ; j < var ; j ++)
		{
			limited_left_p1[j] = -limited_left_p1[j] + meanM_p1[j];
			limited_right_p1[j] = limited_right_p1[j] - meanM_p1[j];

		}*/



	   if (flag)
		for (unsigned int j = 0 ; j <var ; j ++)
			//coeff[i + ngc/2][1 + j * (p + 1)] = (limited_left_p1[j]+limited_right_p1[j])/(2 * sqrt(3));
			coeff[i + ngc/2][1 + j * (p + 1)] = 0;

	}

	free(meanM_p1);
	free(meanL_p1);
	free(meanR_p1);
	free(meanM_p2);
	free(meanL_p2);
	free(meanR_p2);
	
	free(right_p1);
	free(left_p1);
	free(right_p2);
	free(left_p2);

	free(limited_right_p1);
	free(limited_left_p1);
	free(limited_right_p2);
	free(limited_left_p2);
}
