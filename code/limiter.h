
/***********develops the eigen vector matrix for one cell********/
void Char_matrix(double **A,double *coeff, int s)	//A-matrix to be returned, coeff-coefficients of the solution, s-switch for left or right eigen vectors
{
	int i;
	double *u;						//stores the cell average of the three convective variables
	double *phy_var;				//stores the mean value of the physical variables
	double *beta_value;				//stores the value of the beta 

	phy_var = calloc(2,sizeof(double));
	u = calloc(var,sizeof(double));
	beta_value = calloc(2,sizeof(double));

	for(i = 0 ; i < var ; i ++)
		u[i] = coeff[i*(p+1)];

	

	//if ( s == 0)
	//printf ("loop4 %f %f %f \n",coeff[0],coeff[1*(p+1)],coeff[2*(p+1)]);

	

	double rho, velocity, pressure, energy, c, enthalpy;		//primitive variables at the cell averages (c = speed of sound)

	double tau;			//tau = (gamma-1)/c^2
   
    beta(u[3],beta_value);

    physical_variable(beta_value,phy_var);				//evaluation of gamma and PI

	rho = u[0];
	velocity = u[1]/u[0];
	pressure = p_pointwise(u[0],u[1]/u[0],u[2],beta_value);
	energy = u[2];
	c = c_sound(u);
	tau = (phy_var[0]-1)/pow(c,2);
	enthalpy = (energy + pressure)/rho;

	//development of R
	if (s == 0)	
	{
		for (i = 0 ; i < 3 ; i ++)
			A[0][i] = 1;

		A[1][0] = velocity - c;
		A[1][1] = velocity;
		A[1][2] = velocity + c;

		A[2][0] = enthalpy - velocity*c;
		A[2][1] = 0.5*pow(velocity,2);
		A[2][2] = enthalpy + velocity*c;

		A[3][3] = 1;
		//printf ("h %f v %f c %f \n",enthalpy,velocity,c);
	}


	//development of R_inv
	if (s == 1)
	{
		A[0][0] = 0.5*(1 + tau*(pow(velocity,2) - enthalpy) + velocity/c);
		A[0][1] = -0.5*(tau * velocity + 1/c);
		A[0][2] = 0.5 * tau;

		A[1][0] = -tau * (pow(velocity,2) - enthalpy);
		A[1][1] = tau * velocity;
		A[1][2] = -tau;

		A[2][0] = 0.5* (1 + tau * (pow(velocity,2) - enthalpy) - velocity/c);
		A[2][1] = -0.5 * (tau * velocity - 1/c);
		A[2][2] = 0.5 * tau;

		A[3][3] = 1;


	}

}

/*******converts from convective to characteristic variable (one cell at a time)*******/
void Con_Char(double *coeff,double *w)
{
	double **R_inv;
	int i,j,k;
	double sum,sum1;		// sum for cell average, sum0 for first order coeff


	R_inv = calloc(var,sizeof(double));
	
	for (i = 0 ; i < var ; i ++)
		R_inv[i] = calloc(var,sizeof(double));
	
	Char_matrix(R_inv,coeff,1);

	
	
	for (i = 0 ; i < var ; i ++)
	{
		sum = 0;
		sum1 = 0;
			for (j = 0 ; j < var ; j ++)
			{
				sum = sum + R_inv[i][j] * coeff[j*(p+1)];		//zeroth order coeff
			  if (p)							//first order does not exist for p = 0
				sum1 = sum1 + R_inv[i][j] * coeff[j * (p+1) + 1];	//first order coeff
			}
		w[i*2] = sum;
		w[i*2 + 1] = sqrt(3.0)*sum1;			//wx = sqrt(3)*w1
	}

	for (i = 0 ; i < var ; i ++)
		free(R_inv[i]);

}

/******converts from characteristic variables to convective variables**********/
void Char_Con(double *w,double *coeff)
{
	double sum,sum1;
	double **R;
	int i,j;

	R = calloc(var,sizeof(double));
	
	for (i = 0 ; i < var ; i ++)
		R[i] = calloc(var,sizeof(double));
	
	Char_matrix(R,coeff,0);

      //printf ("loop3 %f %f %f \n",coeff[0],coeff[1*(p+1)],coeff[2*(p+1)]);

	for (i = 0 ; i < var ; i ++)
	{
		sum = 0;
		sum1 = 0;
			for (j = 0 ; j < var ; j ++)
				{
					sum = sum + R[i][j] * w[j*2];			//zeroth order coeff
					sum1 = sum1 + R[i][j] * w[j*2 + 1];		//first order coeff
					//printf ("R %f",R[i][j]);
				}
		//printf ("\n");
		//printf ("value of sum1 %f \n",sum);
		coeff[i*(p+1)] = sum;
		if (p)
		coeff[i*(p+1) + 1] = sum1/sqrt(3);
	}
	
	for (i = 0 ; i < var ; i ++)
		free (R[i]);
	

}

/****function to perform the limiting on the characteristic variables using the shu limiter*****/
/**CAUTION: limiting has only been applied to the cells in the interior of the domain********/
void limiter_euler(double **coeff)
{

	double **w;			//characteristic variable for every cell
	int i,j,k,m,n;
	double *wxtilda;			//modified values for the gradient from the shu limiter

	w = calloc(N + ngc,sizeof(double));
	wxtilda = calloc(var,sizeof(double));		//modified gradient values for one cell

	for (i = 0 ; i < N + ngc ; i++)
		w[i] = calloc(var*2,sizeof(double));			//3 characteristic variable for every cell, each with 2 coeffs

	for (i = 0 ; i < N + ngc ; i ++)
		Con_Char(coeff[i],w[i]);						//converts from convective to characteristic variables


	for (i = 1 + ngc/2 ; i < N-1 + ngc/2 ; i++)					//dont limit the values at the boundaries
	{
		//values from shu limiter
		for (j = 0 ; j < var ; j++)
			wxtilda[j] =  shu(w[i][j*2 + 1], w[i+1][j*2] - w[i][j*2], w[i][j*2] - w[i-1][j*2]);

		//comparison of the values
		for (j = 0 ; j < var ; j ++)
		{
		
			if (wxtilda[j] != w[i][2*j + 1])
			{

				for (k = 0 ; k < var ; k ++)					
						for (m = 2 ; m < p+1; m ++)
							coeff[i][k*(p+1) + m] = 0;		//setting all the higher order coeffs to zero

				for ( k = 0 ; k < var ; k ++)
					w[i][2*k + 1] = wxtilda[k];

				Char_Con(w[i],coeff[i]);
			
				break;
			}
		}
	}

		

	for (i = 0 ; i < N + ngc; i ++)
		free(w[i]);

	free(wxtilda);

}

	
	


					








	
				



	
	
