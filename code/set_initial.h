set_initial(double **coeff,double (*values)[ngp],double delta_x,double *x_mid,double *x,const int test_case)
{

  assert(var == 5 );

 const double shock_location = 0.0;
 const double state_left_shock[var] = {7,0,1,0,0};
 const double state_right_shock[var] = {1,0,1,0,0};
 const unsigned int alpha = 0;

 double **u_quad;
 u_quad = calloc(ngp,sizeof(double));
 		
 for (unsigned int j = 0 ; j < ngp ; j ++)
 	u_quad[j] = calloc(var,sizeof(double));

 // Riemann problem for the Grad's-13 case
 if (test_case == 0)
 for (unsigned int i = ngc/2 ; i < N + ngc/2 ; i++)
	{
				if (x_mid[i] < shock_location)
				 for (unsigned int j = 0 ; j < var ; j++)
					coeff[i][j*(p+1)] = state_left_shock[j];		//density
					

				if (x_mid[i] >= shock_location)
 					for (unsigned int j = 0 ; j < var ; j++)
						coeff[i][j*(p+1)] = state_right_shock[j];		//density
						
	}

 if (test_case == 1)
 	for (unsigned int i = 0 ; i < N + ngc ; i ++)
 	{

			for (unsigned int k = 0 ; k < ngp ; k ++)	
			{
			 //needed to be check, the index for x
				double loc_x = (delta_x*(0.5*(x_quad[k]+1)) + x[i]);

				u_quad[k][0] = sin(2 * M_PI * loc_x)+2;		// number density
				u_quad[k][1] = sin(2 * M_PI * loc_x);		// velocity
				u_quad[k][2] = sin(2 * M_PI * loc_x) + 2;	// theta
				u_quad[k][3] = alpha*sin(2 * M_PI * loc_x);		// f3
				u_quad[k][4] = alpha*sin(2 * M_PI * loc_x);		// f4
			}

		inv_values_quad(u_quad,coeff[i],values);
 	}
 

}