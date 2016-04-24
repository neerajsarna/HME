void empty(double *a,const int size)
{
	for (unsigned int i = 0 ; i < size ; i ++)
		a[i] = 0;
}

void empty_2d(double **a , const unsigned int m,const unsigned int n)
{
	for (unsigned int i = 0 ; i < m ; i ++)
		for (unsigned int j = 0 ; j < n ; j ++)
			a[i][j] = 0;
}