
 //deallocation
  free(grid.x_mid);
  free(grid.x);
 
 for (i = 0 ; i < N + ngc ; i ++)
	{
		free(solution.coeff[i]);
		free(solution.temp[i]);
	}

 for (i = 0 ; i < N ; i ++)
 {
 		free(solution.quadrature[i]);
		free(solution.quadtemp[i]);
 }

 for(i = 0 ; i < N+1 ; i ++)
	{
		free(solution.q[i]);
		free(solution.qtemp[i]);
	}


  for (i = 0 ; i < N ; i ++)
 {
 		free(noncons.quad[i]);
		free(noncons.quadtemp[i]);
 }

 for(i = 0 ; i < N+1 ; i ++)
	{
		free(noncons.qr[i]);
		free(noncons.ql[i]);

		free(noncons.qrtemp[i]);
		free(noncons.qltemp[i]);
	}





// freeing the 2-d pointer
	free(solution.coeff);
	free(solution.temp);

	free(solution.quadrature);
	free(solution.quadtemp);
	free(solution.q);
	free(solution.qtemp);


 	free(noncons.quad);
	free(noncons.quadtemp);
 
	free(noncons.qr);
	free(noncons.ql);

	free(noncons.qrtemp);
	free(noncons.qltemp);
	

