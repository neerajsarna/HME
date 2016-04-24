

// the following function evaluates the maximum characteristic speed
double time_stepping(double **coeff)
{
	double *sound;		//sound speed for every cell, boundary included
	double max = 0;			//maximum of the charecteristic speed
    double *u_avg;				//stores the average value of variables for every cell

	sound = calloc(N + ngc,sizeof(double));
	u_avg = calloc(var,sizeof(double));


//evaluates the pressure in every cell
	for (unsigned int i = 0 ; i < N + ngc ;i++ )
	{
		for (unsigned int j = 0 ; j < var ; j++)
			u_avg[j] = coeff[i][j * ( p + 1 )]; 
	 
	 	sound[i] = c_sound(u_avg);
	 	
	 	if (fabs(u_avg[1]) + 3 * sound[i] > max )
	 		max = fabs(u_avg[1]) + 3 * sound[i]; 		// a very crude approximation of the characteristic speed

	}

	return(max);
	
}


