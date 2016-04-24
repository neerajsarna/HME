void develop_G(double **G,const double *solution)
{
	const double density = solution[0];
	const double velocity = solution[1];
	const double theta = solution[2];
	const double f3 = solution[3];
	const double f4 = solution[4];

	/**********/
	/* (velocity, density, 0, 0, 0) */
	/*(theta/density,velocity,1,0,0)*/
	/*(0,2theta,velocity,6/density,0)*/
	/*(0,4 * f3,density*theta/2,velocity,4)*/
	/*(-f3 * theta/density,0,-f3,theta,velocity)*/
	/*********/

	G[0][0] = velocity;
	G[0][1] = density;
	
	G[1][0] = theta/density;
	G[1][1] = velocity;
	G[1][2] = 1;

	G[2][1] = 2 * theta;
	G[2][2] =  velocity;
	G[2][3] = 6 / density;

	G[3][1] = 4 * f3;
	G[3][2] = density * theta / 2;
	G[3][3] = velocity;
	G[3][4] = 4.0;

	G[4][0] = -f3 * theta / density;
	G[4][2] = -f3;
	G[4][3] = theta;
	G[4][4] = velocity;

}

