// the following file contains a limiter on primitive variables for the euler equations. Only a single gas case has been implemented
void mean_convective(double *mean,double *coeff)
{
	double density;
	double momentum;
	double energy;

	density = coeff[0];
	momentum = coeff[0 + 1 * (p + 1)];
	energy = coeff[0 + 2 * ( p + 1)];

	mean[0] = density;
	mean[1] = momentum;
	mean[2] = energy;
}


void mean_primtive(double *mean_prim,double *mean_conv)
{
	double density;
	double velocity;
	double energy;
	double pressure;

	density = mean_conv[0];

	velocity = mean_conv[1]/mean_conv[0];

	energy = mean_conv[2];

	pressure = (gamma_1 - 1)*(energy - 0.5 * density * pow(velocity,2));

	mean_prim[0] = density;

	mean_prim[1] = velocity;

	mean_prim[2] = pressure;
}


// higher order of the convective variables
void interface_convective(double *right_con,double *coeff,double (*b_values)[2])
{
	double density = 0;
	double momentum = 0;
	double energy = 0;
	int i;

	for(i = 0 ; i < p + 1 ; i ++)
	{
		density = b_values[i][1] * coeff[i] + density;
		momentum = b_values[i][1] * coeff[i + 1 * (p+1)] + momentum;
		energy = b_values[i][1] * coeff[i + 2 * (p+1) ] + energy;
	}

	right_con[0] = density;
	right_con[1] = momentum;
	right_con[2] = energy;
}


// higher order of the primitve variables
void interface_primitive(double *right_prim,double *right_con,double *mean_con)
{
	// mean values of the convective variables
	double m_density;
	double m_momentum;
	double m_velocity;
	double m_energy;

	//high order values of the convective variables
	double right_density;
	double right_momentum;
	double right_energy;

	//high order values of the primitive variables
	double right_velocity;
	double right_pressure;

	m_density = mean_con[0];
	m_momentum = mean_con[1];
	m_energy = mean_con[2];
	m_velocity = m_momentum/m_density;

	right_density = right_con[0];
	right_momentum = right_con[1];
	right_energy = right_con[2];

	// computation of the higher order of primitive variables
	right_prim[0] = right_density;

	right_velocity = right_momentum/right_density;

	right_pressure = (gamma_1-1) * (right_energy - 0.5 * right_density * pow(right_velocity,2));

	right_prim[1] = right_velocity;

	right_prim[2] = right_pressure; 


}


//the following routine converts back to the convective variables
void inv_interface_conv(double *tilda_conv,double *tilda_prim,double *mean_prim)
{
	double m_density;
	double m_velocity;
	double m_pressure;

	double right_density;
	double right_velocity;
	double right_pressure;

	double right_momentum;
	double right_energy;

	m_density = mean_prim[0];
	m_velocity = mean_prim[1];
	m_pressure = mean_prim[2];

	right_density = tilda_prim[0];
	right_velocity = tilda_prim[1];
	right_pressure = tilda_prim[2];

	right_momentum = right_density * right_velocity;

	right_energy  = right_pressure/(gamma_1-1) +  0.5 * right_density * pow(right_velocity,2);

	tilda_conv[0] = right_density;

	tilda_conv[1] = right_momentum;

	tilda_conv[2] = right_energy;


}

//evaluates the minium out of three positive numbers

double min(double a, double b, double c)
{
	double minimum;
	minimum =a ;
	
	if (minimum > b)
		minimum = b;

	if (minimum >c)
		minimum = c;
	
	return(minimum);
}

//sign of a double, 0 has no sign 

int sign(double a)
{
	return (a > 0 ? 1 : ( a < 0 ? -1 : 0));
}

/****Shu limiter*****/
double shu(double a,double b,double c)
{
	double M=50;		//the shu constant (how to select M2?),presently working on the van leer limiter
    double alpha;
    double delta_x = (xr-xl)/N;



		
	

	if (fabs(a) <= M*delta_x*delta_x)
		return(a);
	
		

	if (fabs(a) > M*delta_x*delta_x)
		{
			if (sign(a) == sign(b) && sign(b) == sign(c))
				return (sign(a)*min(fabs(a),fabs(b),fabs(c)));
			else 
				return (0);
		}

}






void limiter_singlegas(double **coeff,double (*b_values)[2])
{
	double **mean_con;			//mean values of the convective variables
	double **mean_prim;			//mean values of the primitive variables
	double **right_con;			//value of the conservative variables at the right of the cell 
	double **right_prim;		//value of the primitive variables at the right of the cell
	double **tilda_prim;		//modified higher order of the primitive variables
	double **tilda_conv;		//modified value of the convective variables
	double meanR,meanM,meanL,barM;			//variables for limiter application
	int i,j,k;


	mean_con = calloc(N + ngc,sizeof(double));
	mean_prim = calloc(N + ngc,sizeof(double));
	right_con = calloc(N + ngc,sizeof(double));
	right_prim = calloc(N + ngc,sizeof(double));
	tilda_prim = calloc(N,sizeof(double));
	tilda_conv = calloc(N,sizeof(double));

	//printf("size of mean %lu\n",sizeof(*mean_prim));
	for(i = 0 ; i < N + ngc ; i ++)
	{
		mean_con[i] = calloc(var-1,sizeof(double));			//phi not considered
		mean_prim[i] = calloc(var-1,sizeof(double));

		right_con[i] = calloc(var-1,sizeof(double));
		right_prim[i] = calloc(var-1,sizeof(double));

	
	}

	for(i = 0 ; i < N ; i ++)
	{
		tilda_prim[i] = calloc(var-1,sizeof(double));
		tilda_conv[i] = calloc(var-1,sizeof(double));
	}


	
	for(i = 0 ; i < N + ngc ; i ++)
		mean_convective(mean_con[i],coeff[i]);			//routine to evaluate the mean of the convective variables
	

	for(i = 0 ; i < N + ngc ; i ++)
		mean_primtive(mean_prim[i],mean_con[i]);		//routine to evaluate the mean of the primitive variables


	for(i = 0 ; i < N + ngc ; i ++)
		interface_convective(right_con[i],coeff[i],b_values);

	for(i = 0 ; i < N + ngc ; i ++)
		interface_primitive(right_prim[i],right_con[i],mean_con[i]);



	// application of the limiter
	for ( i = 0 ; i < N ; i ++ )					//loop over the interior cells
		for(j = 0 ; j < var - 1 ; j ++)				//loop over all the possible variables
		{
			meanR = mean_prim[i + ngc/2 + 1][j];	//mean value to the right of the current cell
			meanM = mean_prim[i + ngc/2][j];		//mean value in the current cell
			meanL = mean_prim[i + ngc/2 - 1][j];	//mean value to the left of the current cell

			barM = right_prim[i + ngc/2][j]-meanM;			//higher order of the current cell

			tilda_prim[i][j] = shu(barM,meanR-meanM,meanM-meanL);

			tilda_prim[i][j] = tilda_prim[i][j] + meanM;	// converting to value to the right of the cell

			//printf("%f %f %d\n",tilda_prim[i][j],barM,j);
		}


	// converting back to the convective variables
		for(i = 0 ; i < N ; i ++)
			inv_interface_conv(tilda_conv[i],tilda_prim[i],mean_prim[i + ngc/2]);			//converting back to convective variables


 /*
	// for testing
		for(i = 0 ; i < N ; i ++)
			inv_high_conv(high_con[i + ngc/2],high_prim[i + ngc/2],mean_prim[i + ngc/2]);			//converting back to convective variables*/

	//evaluation of the higher order coeffs, only for p = 1
		for(i = 0 ; i < N ; i ++)
			for(j = 0 ; j < var - 1 ; j ++)
				for(k = 1 ; k < p + 1 ; k ++)
				    coeff[i + ngc/2][k + j *(p+1)] = (tilda_conv[i][j]-mean_con[i + ngc/2][j])/sqrt(3);



	for(i = 0 ; i < N + ngc ; i ++)
	{
		free(mean_con[i]);
		free(mean_prim[i]);
		free(right_prim[i]);
		free(right_con[i]);
	}

	for(i = 0 ; i < N ; i ++)
	{
		free(tilda_prim[i]);
		free(tilda_conv[i]);
	}
	

}