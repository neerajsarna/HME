//evaluates the speed of sound, the variable u contains (density, momentum, energy, phi)
double c_sound(double *u_average)
{
	double theta = u_average[2];

	assert(theta > 0 );

	return(sqrt(gamma * 3 * theta / 2));
}
