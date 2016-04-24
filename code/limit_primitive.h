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
    double delta_x = 1/N;

		

	

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


// evaluates the higher order for pressure
void higher_order_pressure(double *bar,double *mean,double *mean_phy_var,double *high_phy_var,double *test)
{
	#include "variable.h"
	int i;
	int j_press = 4;		//index for pressure in bar
	int j_density = 0;					//index for density in bar
	int j_velocity = 1;					//index for velocity in bar
	int j_energy = 2;					//index for energy in bar
	int j_phi = 3;						//index for phi in bar
	int j_gamma = 0;
	int j_pi = 1;

    for(i = 0 ; i < 2 ; i ++)
    {
    	m_gamma = mean_phy_var[j_gamma];
    	m_pi = mean_phy_var[j_pi];
    	m_density = mean[j_density];
    	m_velocity = mean[j_velocity];
    	m_energy = mean[j_energy];

    	h_gamma = high_phy_var[i + 2 * j_gamma];
    	h_pi = high_phy_var[i + 2 * j_pi];
    	h_density = bar[i + 2 * j_density];
    	h_velocity = bar[i + 2 * j_velocity];
    	h_energy = bar[i + 2 * j_energy];

    	test[0] = bar[0 + 2 * j_energy];
    	test[1] = bar[1 + 2 * j_energy];

    	bar[i + 2 * j_press] = (m_gamma-1)*(h_energy-m_density*m_velocity*h_velocity-0.5*pow(m_velocity,2)*h_density) + h_gamma*(m_energy-0.5*m_density*pow(m_velocity,2))-m_pi*h_gamma-m_gamma*h_pi;
    	//bar[i + 2 * j_press] = (mean_phy_var[j_gamma] - 1)*(bar[i + 2 * j_energy] - mean[j_density] * mean[j_velocity] * bar[i + 2 * j_velocity] - bar[i + 2 * j_density] * pow(mean[j_velocity],2)/2) + high_phy_var[i + 2 * j_gamma]*(mean[j_energy] - 0.5 * mean[j_density] * pow(mean[j_velocity],2)-mean_phy_var[j_pi]) - high_phy_var[i + 2 * j_pi] * mean_phy_var[j_gamma];
    	//h_energy = (bar[i + 2 * j_press] - (h_gamma*(m_energy-0.5*m_density*pow(m_velocity,2))-m_pi*h_gamma-m_gamma*h_pi))/(m_gamma-1)+(m_density*m_velocity*h_velocity+0.5*pow(m_velocity,2)*h_density);

    }

		


   
}


// evaluates the high order for velocity
void higher_order_velocity(double *bar,double *mean)
{
	int j_momentum = 1;		//index for momentum 
	int j_density = 0;		//index for density

	bar[0 + 2*j_momentum] = (bar[0 + 2 * j_momentum] - bar[0 + 2 * j_density] * mean[j_momentum])/(mean[j_density]);	//update for x + 1/2 i.e minus
	bar[1 + 2 * j_momentum] = (bar[1 + 2 * j_density] * mean[j_momentum] - bar[1 + 2 * j_momentum])/(-mean[j_density]);
}






//evaluates the higher order for the convective variables
void higher_order(double *coeff,double *mean,double *bar, double (*b_values)[2], int prim_variables)
{
	int i,j;

	/*Order of variables in bar
	1. density	
	2. velocity
	3. energy
	4. phi
	5. pressure
	*/
	int j_press = prim_variables-1;		//index for pressure in bar
	int j_density = 0;					//index for density in bar
	int j_velocity = 1;					//index for velocity in bar
	int j_energy = 2;					//index for energy in bar
	int j_phi = 3;						//index for phi in bar


	for(i = 0 ; i < prim_variables-1 ; i ++)	//evaluation of higher order coeffs except pressure
	{
		for(j = 1 ; j < p + 1 ; j ++)
		{
			bar[0 + 2 * i] = coeff[j + i * (p + 1)] * b_values[j][1] + bar[0 + 2 * i];
			bar[1 + 2 * i] = -coeff[j + i * (p + 1)] * b_values[j][0] + bar[1 + 2 * i];
		}

	}
		
	higher_order_velocity(bar,mean);	//evaluation of the higher order coefficients for velocity

}




// the main function for applying the limiter to the values
void limiter_primitive(double **coeff, double (*b_values)[2])
{
	int i,j,k;
	int prim_variables = 5; 			// including pressure
	int j_press = prim_variables-1;		//index for pressure in bar
	int j_density = 0;					//index for density in bar
	int j_velocity = 1;					//index for velocity in bar
	int j_energy = 2;					//index for energy in bar
	int j_phi = 3;						//index for phi in bar
	double **mean_phy_var;		//stores the mean of the physical variables i.e gamma and PI for every cell

	double **high_phy_var;		//stores the higher order coefficients of the physical variables, first two coloumns for gamma then for PI
	double *tilda_high_phy_var;	//stores the value of the physical variables obtained by altering the value using Shu limiter

	/*The following variable stores the mean values of the following variables (in order):
	1. density
	2. velocity
	3. energy
	4. phi
	5. pressure
	for every cell
	*/
	double **mean;

	/*stores the higher order contribution from every convective variable
	first coloumn of every variable for minus
	second coloumn of every variable for plus (x + 1/2)
	the order of storage is same as the mean
	refer notes for plus and minus
	*/
	double **bar;

	/*Data structure similar to bar, contains altered values from the shu limiter
	*/
	double **tilda;

	// allocation of memory including ghost cells
	mean = calloc(N + ngc,sizeof(double));
	bar = calloc(N + ngc,sizeof(double));
	tilda = calloc(N + ngc,sizeof(double));

	// variables for the physical variables
	mean_phy_var = calloc(N + ngc,sizeof(double));
	high_phy_var = calloc(N + ngc,sizeof(double));
	tilda_high_phy_var = calloc(4,sizeof(double));

	for(i = 0 ; i < N + ngc ; i ++)
	{
		mean[i] = calloc(prim_variables,sizeof(double));				// one value per variable per cell
		bar[i] = calloc(2 * prim_variables,sizeof(double));				// two values per variable per cell
		tilda[i] = calloc(2 * prim_variables,sizeof(double));			// two values per variable per cell
		mean_phy_var[i] = calloc(2,sizeof(double));						// two means, one for gamma and one for PI
		high_phy_var[i] = calloc(4,sizeof(double));						// two higher order value for each gamma and PI
	}
	
		double test1,test2;
		double *test;
		test = calloc(2,sizeof(double));

	//routine to evaluate means and higher order coefficients
	for(i = 0 ; i < N + ngc ; i ++)
	{
		for(j = 0 ; j < prim_variables-1 ; j ++)				//evaluation of all the means except pressure
			mean[i][j] = coeff[i][j*(p+1)];
		
		// original mean[i][1] = mean[i][1]/mean[i][0];, now for the test case
		mean[i][1] = mean[i][1]/mean[i][0];						//conversion of momentum to velocity
			

		// mean value for the pressure, send in order: density, velocity, energy, phi
	    mean[i][prim_variables-1] = p_pointwise(mean[i][0],mean[i][1],mean[i][2],mean[i][3]);	


	    // mean value of gamma and PI
	    physical_variable(mean[i][j_phi],mean_phy_var[i]);



	    higher_order(coeff[i],mean[i],bar[i],b_values,prim_variables);					//evaluation of higher order coefficiens (except pressure)

	    high_physical_variable(bar[i],mean_phy_var[i],high_phy_var[i]);	//evaluates the higher order of the physical coefficients //ref pres_sound.h


	    higher_order_pressure(bar[i],mean[i],mean_phy_var[i],high_phy_var[i],test);
		

		
	  /* test1 = bar[i][0 + 2 * j_energy];
	    test2 = bar[i][1 + 2 * j_energy];
		
		// evalution of the value of the energy using altered value of pressure from above
		e_evaluation(mean[i],bar[i],mean_phy_var[i],high_phy_var[i],test);

		/*if(test1 != bar[i][0 + 2 * j_energy])
			printf("idiot\n");


		if(test2 != bar[i][1 + 2 * j_energy])
			printf("idiot\n");*/

	}



	
	int index;

	//application of the limiter to the interior cells, excluding the ghost cells
	for (i = ngc/2; i < N + ngc/2 ; i++)
	{
		for(j = 0 ; j < prim_variables ; j ++)
		{
			//if ( j == 2)		// no limiting the value of energy, evaluated with pressure latter
			//	continue;

			// CAUTION: value form ghost cells used, apply boundary conditions first
			for(k = 0 ; k < 2 ; k ++)
				 tilda[i][k + 2*j] = shu(bar[i][k + 2 * j], mean[i+1][j]-mean[i][j],mean[i][j]-mean[i - 1][j]);

			
		}

		index = j_press;
		// evauation of the altered value of the higher order coefficients of the physical variables
	   printf("%e %e  final \n",tilda[i][0 + 2 * index]-bar[i][0 + 2 * index],tilda[i][1 + 2 * index]-bar[i][1 + 2 * index]);	

		high_physical_variable(tilda[i],mean_phy_var[i],tilda_high_phy_var);
		
		

		// evalution of the value of the energy using altered value of pressure from above
		e_evaluation(mean[i],tilda[i],mean_phy_var[i],tilda_high_phy_var,test);


		// evaluation of the new coeffs using the 
		if ( p == 1)
		{
			for(j = 0 ; j < var ; j ++)
			{
				if (j == 1)			//update for momentum
					coeff[i][p + j*(p+1)] = (tilda[i][0 + 2 * j_density] * mean[i][j_velocity] + tilda[i][0 + 2 * j_velocity] * mean[i][j_density])/sqrt(3);


				else				//update for other variables
					coeff[i][p + j*(p+1)] = (tilda[i][0 + 2 * j])/(sqrt(3));

			}
		}

		// CAUTION: update for momentum has not been corrected yet, update from the top
		if (p == 2)
		{
			for (j = 0; j < var; j++)
			{
				if ( j == 1 )		//update for momentum
				{
					coeff[i][p + j*(p + 1)] = (tilda[i][0 + j * 2] * tilda[i][0 + (j-1) * 2] - tilda[i][1 + j * 2] * tilda[i][1 + (j-1) * 2])/(2 * sqrt(5));
					coeff[i][p-1 + j * (p+1)] = (tilda[i][ 0 + j * 2] * tilda[i][0 + (j-1) * 2] + tilda[i][1 + j * 2] * tilda[i][1 + (j-1) * 2])/(2 * sqrt(3));
				}
				else				//update for other variables
				{
					coeff[i][p + j*(p + 1)] = (tilda[i][0 + j * 2] - tilda[i][1 + j * 2])/(2 * sqrt(5));
					coeff[i][p-1 + j * (p+1)] = (tilda[i][ 0 + j * 2] + tilda[i][1 + j * 2])/(2 * sqrt(3));
				}
				
			}
		}

	}
	 
	 

	for(i = 0 ; i < N + ngc ; i ++)
	{
		free(mean[i]);
		free(bar[i]);
		free(tilda[i]);
	}
}