//Note: uses two extra ghost cells, the ghost cells have only been included. The ghost cells only appear in the definition of the grid and in the definition of the coefficients of the solution.

typedef struct 
{
	double legendre[p+1][p+1];		//coeffs of the legendre
	double der_legendre[p+1][p+1];	//coeffs of the derivatives of legendre
	double values[p+1][ngp];			//values at the quad point 
	double b_values[p+1][2];		//values at the boundary of the cells
	double der_values[p+1][ngp];		//derivative values at the quad points
	double ci1[p+1];				//l-1 norm (see skript)
	double ci_infi[p+1];			//l-inifnity norm	(see skript)
	double cp;
}polynomial;

/***information about the grid******/
typedef struct 
{
	double *x_mid;	//x value at the nodes, no ghost cells
	double *x;	//x value at the interfaces, no ghost cells
}space;


/**variables for the solution****/
/*description of storage of the solution:
for every cell, the first p+1 coefficients are that of rho, the next p+1 are for rho*velocity and the next p+1 are for Energy*/

/*for the update of the solution, ghost cells have not been used rather the boundary conditions will be set on the fluxes */

typedef struct 
{
	//variables for the coefficients of the basis expansion
	double **coeff;			//coefficients of the legendre polynomials for every cell, each set of p+1 elements correspond to one convective variable
    double **temp;			//temporary coefficients for the second order RK scheme

    //variables for the fluxes at the boundary of the interfaces
	double **q;					//value of  numerical flux, from convective part
	double **qtemp;				//value of  numerical flux with RK, from convective part

	//variables for the quadrature
	double **quadrature;		//integral of flux and der_legendre for every cell
	double **quadtemp;			//temporary integral variable for second order RK scheme


	//stores the value of the average of legendre polynomials  for every cell
	double **avg_legendre;	//average value of the legendre polynomials in every cell 
	
}data;

typedef struct 
{
	double **qr;				//value of numerical flux from non conservative part for right side of the interface
	double **ql;				//value of numerical flux from non conservative part for left side of the interface
	double **qrtemp;			//value of numerical flux for right side of the interface for RK stage
	double **qltemp;			//value of numerical flux for left side of the interface for RK stage
	double **quad;				//value of the integral from non conservative part
	double **quadtemp;
	double **avg_basis;		//average of the values of the test functions at the boundary
}stiff;

