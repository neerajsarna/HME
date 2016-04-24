// memory allocation for all the pointers
  grid.x_mid = calloc(N+ngc,sizeof(double));
  grid.x = calloc(N+1+ngc,sizeof(double));

  solution.coeff = calloc(N+ngc,sizeof(double));			//coefficients for the ghost cells
  solution.temp = calloc(N+ngc,sizeof(double));
  solution.q = calloc(N+1,sizeof(double));
  solution.qtemp = calloc(N+1,sizeof(double));

  noncons.qr = calloc(N+1,sizeof(double));
  noncons.qrtemp = calloc(N+1,sizeof(double));

  noncons.ql = calloc(N+1,sizeof(double));
  noncons.qltemp = calloc(N+1,sizeof(double));

  noncons.quad = calloc(N,sizeof(double));
  noncons.quadtemp = calloc(N,sizeof(double));

  solution.quadrature = calloc(N,sizeof(double));
  solution.quadtemp = calloc(N,sizeof(double));

  noncons.avg_basis = calloc(N,sizeof(double));

  for (i = 0 ; i < N + ngc ; i ++)
	{
		solution.coeff[i] = calloc(var*(p+1),sizeof(double));
		solution.temp[i] =  calloc(var*(p+1),sizeof(double));
	}


// the following quantities are only for the interior of the cells
 for (i = 0 ; i < N ; i++)
	{
		solution.quadrature[i] = calloc(var*(p+1),sizeof(double));
		solution.quadtemp[i] = calloc(var*(p+1),sizeof(double));


		noncons.quad[i] = calloc(var*(p+1),sizeof(double));
		noncons.quadtemp[i] = calloc(var*(p+1),sizeof(double));

		noncons.avg_basis[i] = calloc(p+1,sizeof(double));		
	}

  for (i = 0 ; i < N+1 ; i ++)
	{
		solution.q[i] = calloc(var,sizeof(double));
		solution.qtemp[i] = calloc(var,sizeof(double));

  		noncons.qr[i] = calloc(var,sizeof(double));			//left and right at every interface 
  		noncons.qrtemp[i] = calloc(var,sizeof(double));

  		noncons.ql[i] = calloc(var,sizeof(double));			//left and right at every interface 
  		noncons.qltemp[i] = calloc(var,sizeof(double));
	}