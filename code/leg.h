/*********develops the legendre polynomials for 1-d********************/
/*Note:the legendre polynomials have not been shifted to every cell, x belongs from [0,1],divide by delta_x for tilda and shift the coordinates for the normal one*/

void leg(double (*legendre)[p+1],double (*der_legendre)[p+1])
{
 int i,j;
 for (i = 0 ; i < p+1 ; i++)
	for (j = 0 ; j < p +1 ; j++)
		{
			legendre[i][j] =0;
			der_legendre[i][j] =0;
		}

 legendre[0][p] = 1;

 
 if (p)
 {
 	legendre[1][p-1]= 2;	//removed for p=0;
	 legendre[1][p] = -1;		//shifting the coordinates to 2*x-1
 }
 

 for (i = 2 ; i <p+1; i ++)
	{
		for (j = 0; j < p+1; j++)
			{
					legendre[i][j-1] = legendre[i][j-1]+2*((2*(i-1)+1)*legendre[i-1][j])/(i);

					legendre[i][j] = legendre[i][j]-((i-1)*legendre[i-2][j])/i-((2*(i-1)+1)*legendre[i-1][j])/(i);
			}

	}


 
 /************scaling for normalization*************/
 for (i = 0 ; i < p+1; i++)
   for (j = 0 ; j < p+1; j++)
		legendre[i][j] = sqrt(2*i+1)*legendre[i][j];


 /************derivatives of legendre polynomials******/

 for (i = 0 ; i < p+1; i++)
		for (j = 0 ; j < p; j++)
			der_legendre[i][j+1] = legendre[i][j]*(p-j);
			

 
}

/*****evaluates the value of the legendre polynomials at the quadrature points*****/
void quadrature_values(double (*legendre)[p+1],double (*der_legendre)[p+1], double (*values)[ngp], double (*der_values)[ngp])
{
 int i,k,m;


 for (i = 0 ; i < p+1 ;i++)
	{
		for (k = 0 ; k < ngp ; k++)
		{
			values[i][k] = 0 ;
			der_values[i][k] = 0;
		}
	}


 for (i = 0 ; i < p+1 ;i++)
		for (m = 0 ; m < ngp ; m ++)
		 for (k = 0 ; k < p+1; k++)
		{
			
				
			values[i][m] = values[i][m] + legendre[i][k]*pow((x_quad[m]+1)/2,p-k);


			der_values[i][m] = der_values[i][m] + der_legendre[i][k]*pow((x_quad[m]+1)/2,p-k);
			
		}

	

}

/*********evaluates the value of the polynomials at the boundaries**************/
void boundary_values(double (*legendre)[p+1], double (*b_values)[2])
{

 int i,j;

 for (i = 0 ; i < p+1 ; i++)
	{
		for (j = 0 ; j < 2 ; j++)
			b_values[i][j] = 0;
	}

 for (i = 0 ; i < p+1 ; i++)
	{
		  for (j = 0 ; j < p+1 ; j++)
		     {
			b_values[i][0] = b_values[i][0] + legendre[i][j]*pow(0,p-j);
			b_values[i][1] = b_values[i][1] + legendre[i][j]*pow(1,p-j);
		     }	
	}
}

/*****evaluates the maximum of an array****/
double array_max(double *v,int size)
{
	int i;
	
	double max = fabs(v[0]);

	
	
	for (i = 1; i <size ; i++)
		if (fabs(v[i]) > max)
			max = fabs(v[i]);

	for (i = 0 ; i < size ; i++)
		v[i] = 0;


	return(max);
}

void norms(double (*legendre)[p+1],double (*values)[4],double *ci1, double *ci_infi,double *cp)
{
	/**division of interval from 0 to 1 to find the maximum**/
	int spacing = 1000;		 //number of intervals	
	double interval = 1*1.0/spacing; //spacing between nodes
	double x[spacing+1];		 //location of nodes
	int i,j,k;
	*cp = 0;

	for (i = 0 ; i < p+1 ; i++)
		ci1[i] = 0;


	for (i = 0 ; i < spacing+1 ; i++)	//local discretization for evaluating the l-1 norm of the legendre
		x[i] = i*interval;

	/****evaluation of the l_infinity norm****/
	
	
	double *legendre_value;
	legendre_value = calloc(spacing+1,sizeof(double));
	

	
	   for (j = 0 ; j < p+1 ; j++)		//for legendre
	   {
		
		for (i = 0 ; i < spacing+1; i++)	//for grid points
			for (k = 0 ; k<p+1; k++)
				legendre_value[i] = legendre[j][k]*pow(x[i],p-k)+legendre_value[i];
				
		ci_infi[j] = array_max(legendre_value,spacing+1);
	   }

	/**evalution of the l-1 norm**/

	  for (i = 0; i < p+1 ; i++)
			for (j = 0 ; j < 4; j++)
				ci1[i] = 0.5*w[j]*fabs(values[i][j]) + ci1[i];


	 for (i = 0 ; i < p+1 ; i++)
		*cp = ci1[i]*ci_infi[i];		
}

// the following function evaluates the average of the basis functions for every cell, tilda basis functions = basis/delta_x
void basis_average(double (*b_values)[2],double **avg_basis)
{
	int count = N;
	int i,j;
	double delta_x = 1.0/N;

	for (i = 0; i < N; i++)
		for (j = 0; j < p+1; j++)
			avg_basis[i][j] = (b_values[j][0] + b_values[j][1])/(2 * delta_x); 
}

