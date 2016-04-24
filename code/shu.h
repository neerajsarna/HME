// contains the functions for the shu limiter

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
	double M=200.0;		
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
