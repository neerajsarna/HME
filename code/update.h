// performs the time update for the solution
void time_update(data *solution,stiff *noncons,polynomial l,double delta_x,double delta_t)
{


	int flag = 1;		// 0 = characteristic limiting, 1 = primitive limiting
	int i,j,k;


 
 	 fluxes_cons(solution->coeff,l.legendre,l.b_values,solution->q,delta_x);		//evaluates the fluxes at each cell interface,ref. fluxes.h
 
 	 quad_cons(solution->quadrature, solution->coeff,l.legendre,l.values,l.der_values);	//evaluates the value of flux*derivative of legendre in every cell,ref. quad.h



	 fluxes_noncons(solution->coeff,noncons->qr,noncons->ql,l.b_values);
 

	 quad_noncons(solution->coeff,noncons->quad,l.values,l.der_values);

		for (i = 0 ; i < N ; i ++)
			for (k = 0 ; k < var ; k++)
				for (j = 0 ; j < p+1 ;  j++)			
					solution->temp[i + ngc/2][j + k*(p+1)] = solution->coeff[i + ngc/2][j + k*(p+1)] + delta_t*(-l.b_values[j][1]*solution->q[i+1][k]/delta_x + l.b_values[j][0]*solution->q[i][k]/delta_x+solution->quadrature[i][j + k*(p+1)]  - noncons->qr[i+1][k]*l.b_values[j][1]/(2 * delta_x) - noncons->qr[i][k] * l.b_values[j][0] / (2 * delta_x)  - noncons->quad[i][j + k *(p+1)]);

	

	set_boundary(solution->temp);

	limiter_moment(solution->temp,l.b_values);
	
		fluxes_cons(solution->temp,l.legendre,l.b_values,solution->qtemp,delta_x);			
 		
		quad_cons(solution->quadtemp, solution->temp,l.legendre,l.values,l.der_values);

		
		fluxes_noncons(solution->temp,noncons->qrtemp,noncons->qltemp,l.b_values);

		quad_noncons(solution->temp,noncons->quadtemp,l.values,l.der_values);
		
		

		//second update in time
		for (i = 0 ; i < N ; i ++)
			for (k = 0 ; k < var ; k++)
				for (j = 0 ; j < p+1 ;  j++)
					solution->coeff[i + ngc/2][j + k*(p+1)] = solution->coeff[i + ngc/2][j + k*(p+1)]*0.5 + solution->temp[i + ngc/2][j + k*(p+1)]*0.5 + delta_t*(-l.b_values[j][1]*solution->qtemp[i+1][k]/delta_x + l.b_values[j][0]*solution->qtemp[i][k]/delta_x+solution->quadtemp[i][j + k*(p+1)]  - noncons->qrtemp[i+1][k] * l.b_values[j][1]/(2 * delta_x) - noncons->qrtemp[i][k] * l.b_values[j][0]/(2 * delta_x) - noncons->quadtemp[i][j + k *(p+1)])*0.5; 

		
		set_boundary(solution->coeff);

		limiter_moment(solution->coeff,l.b_values);

}
