double g11(double theta) {
    return (theta + sin(theta)*cos(theta));
}

double g22(double theta) {
    return (theta - sin(theta)*cos(theta));
}

double g33(double theta) {
    return (theta);
}

double g23(double theta) {
    return (-cos(theta));
}

double g13(double theta) {
    return (sin(theta));
}

double g12(double theta) {
    double st = sin(theta);
    return (st*st);
}

typedef double (*g)(double);
g gxx[7] = {&g11, &g22, &g33, &g23, &g13, &g12, &g33};  

double gr(double r, void * params)
{
	Conic ep = *(Conic *) params;
	
	double F = ep.F - r*r;
	int vindex = ep.index;
  
	double Mp11 = gsl_matrix_get(ep.M1,0,0);
	double Mp22 = gsl_matrix_get(ep.M1,1,1);
	double Mp33 = gsl_matrix_get(ep.M1,2,2);
	
	double a = sqrt((r*r-Mp33)/Mp11);
	double b = sqrt((r*r-Mp33)/Mp22);
	double k = sqrt(1-(b*b)/(a*a));
	
	double alpha[7] = { f(r)*0.5*r*(r*r-Mp33)/(Mp11*sqrt(Mp11*Mp22)),
                      f(r)*0.5*r*(r*r-Mp33)/(Mp22*sqrt(Mp11*Mp22)),
                      f(r)*r/sqrt(Mp11*Mp22),
                      f(r)*r*sqrt(r*r-Mp33)/(Mp22*sqrt(Mp11)),
                      f(r)*r*sqrt(r*r-Mp33)/(Mp11*sqrt(Mp22)),
                      f(r)*0.5*r*(r*r-Mp33)/(Mp11*Mp22),
                      V(r)*r/sqrt(Mp11*Mp22) };

	//printf("\na = %lf, b = %lf, k = %lf, alpha = %lf \n",a,b,k,alpha);
	
	if(r<=(ep.rmin+EPS) || r>=(ep.rmax-EPS)) return 0;
	
	//printf("The modified conic is (%g)x^2+(%g)y^2+(%g)=0\n",Mp11,Mp22,Mp33);
	
	gsl_complex csol[8];
	/* A * x^2 + D * x + F = 0 */
	gsl_poly_complex_solve_quadratic(ep.A, ep.D, F, &csol[0], &csol[1]);	
	gsl_poly_complex_solve_quadratic(ep.C, ep.B+ep.E, ep.A+ep.D+F, &csol[2], &csol[3]);	
	gsl_poly_complex_solve_quadratic(ep.A, ep.B+ep.D, ep.C+ep.E+F, &csol[4], &csol[5]);	
	gsl_poly_complex_solve_quadratic(ep.C, ep.E, F, &csol[6], &csol[7]);	
	//print_roots8(csol);
	
  for (int i = 0; i < 8; i++)
	{
    if(abs(GSL_REAL(csol[i]))<EPS) GSL_SET_REAL(&csol[i], 0);  
  }
  //print_roots8(csol);

	double rsol[24],drsol[24]; int count=0;
	for (int i = 0; i < 4; i++)
	{
		int loc = i%2;
		double val = 1.5-fabs(1.5-i);
		if(fabs(GSL_IMAG(csol[i*2]))<EPS) //real roots
		{
			if( abs(GSL_REAL(csol[i*2])-GSL_REAL(csol[i*2+1])) > EPS_ROOT ) // roots are not identical
			{
				for (int j = 0; j < 2; j++)
				{
					if(abs(GSL_REAL(csol[i*2+j])-0.5)<(0.5+EPS)) //root is within the square
					{
						rsol[count*3+loc]=GSL_REAL(csol[i*2+j]); 
						rsol[count*3+1-loc]=val; 
						rsol[count*3+2]=1;
						
						drsol[count*3+loc]=1-2*(i/2); 
						drsol[count*3+1-loc]=0; 
						drsol[count*3+2]=0;
												
						count++;
					}
				}		
			}
		}
	}
	
	//printf("count = %d\n",count);
	//mat_print1("rsol ",rsol,count,3);
	//if(count==1) count=0; //accounting for corner touch
  
	if(count==0) 
	//ellipse does not intersect square. In this case 
	// (1) ellipse could be within the square, or 
	// (2) square could be within the ellipse, or
	// (3) ellipse and square are separate from each other
	// For (2) and (3) g33 should be 0. For (1) g33 should be related to perimeter of ellipse/circle
	{
		int in_square = (fabs(ep.xc-0.5)<0.5) && (fabs(ep.yc-0.5)<0.5);
    
    if(F>0 && (ep.A+ep.D+F)>0 && (ep.A + ep.B + ep.C + ep.D + ep.E + F) > 0 && (ep.C + ep.E + F) > 0 && in_square)
		// condition (1), as all square corners are outside ellipse 
		{
			//printf("\n ellipse (A,B,C,D,E,F) = (%.3f,%.3f,%.3f,%.3f,%.3f,%.3f) with center (%.3f,%.3f) in_square %d\n",ep.A,ep.B,ep.C,ep.D,ep.E,F,ep.xc,ep.yc,in_square);
      double theta[2][2];
			
			theta[0][0] = -PI; 	theta[0][1] = -1;
			theta[1][0] = PI;	theta[1][1] = 1;
			
			double output=0;
			for (int i = 0; i < 2; i++)
				output += (gxx[vindex](theta[i][0]))*theta[i][1];
			
			return alpha[vindex]*output;
		}
		else return 0; //condition (2) or (3)
	}
	
	gsl_matrix * zsolT = gsl_matrix_alloc (count, 3);
	gsl_matrix * zsol = gsl_matrix_alloc (3,count);
	gsl_matrix * dzsolT = gsl_matrix_alloc (3,count);
	
	for (int i = 0; i < count; i++)
		for (int j = 0; j < 3; j++)
		{
			gsl_matrix_set (zsolT, i, j, rsol[3*i+j]);
			gsl_matrix_set (zsol, j, i, rsol[3*i+j]);
			gsl_matrix_set (dzsolT, j, i, drsol[3*i+j]);
		}
	
	//gsl_mat_print(zsolT,count,3);
	//gsl_mat_print(dzsolT,3,count);
	
	gsl_matrix * dG = gsl_matrix_alloc (count, 1);
	
	for (int i = 0; i < count; i++)
	{
		gsl_matrix_view row = gsl_matrix_submatrix(zsolT,i,0,1,3);
		gsl_matrix_view drow = gsl_matrix_submatrix(dzsolT,0,i,3,1);
		
		gsl_matrix * Mdrow = gsl_matrix_alloc (3, 1);
		gsl_matrix * Gi = gsl_matrix_alloc (1, 1);
		
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
			  1.0, ep.M, &drow.matrix,
			  0.0, Mdrow);
			  
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
			  1.0, &row.matrix, Mdrow,
			  0.0, Gi);
		
		double sign = (gsl_matrix_get(Gi,0,0)>0)-(gsl_matrix_get(Gi,0,0)<0);
		//change// gsl_matrix_set (dG, i, 0, sign);
    gsl_matrix_set (dG, i, 0, -sign);
	}
	
	////gsl_mat_print(dG,count,1);
	
	gsl_matrix * Z = gsl_matrix_alloc (3, count);
	
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		1.0, ep.T, zsol,
		0.0, Z);
	
	//gsl_mat_print(Z,3,count);
			
	double **theta = (double **)malloc(count * sizeof(double *));
	for (int i = 0; i < count; i++) theta[i] = (double *)malloc(2 * sizeof(double));
	
	for (int i = 0; i < count; i++)
	{
		//change// theta[i][0] = atan2( gsl_matrix_get(Z,0,i)/a,gsl_matrix_get(Z,1,i)/b );
    theta[i][0] = atan2( gsl_matrix_get(Z,1,i)/b,gsl_matrix_get(Z,0,i)/a );
		theta[i][1] = gsl_matrix_get(dG,i,0);
	}
	
	//mat_print("theta ",theta,count,2);
	qsort(theta, count, sizeof theta[0], compare); //how am i sorting this?
	//mat_print("theta ",theta,count,2);
	
	//remove dulicates in theta
	for(int i=1; i<count; i++)
	{
		if(fabs(theta[i][0]-theta[i-1][0])<EPS) theta[i][1] = 0; 
	}
		
	if(theta[0][1]==1) theta[count-1][0]-=2*PI; //the is to ensure arcs are not discontinuous
	//mat_print("theta: ",theta,count,2);
	
	//for testing perpendicular fibers (0.5<a<1)
	double output1=0;
	if(vindex==6 && ep.err_test == 1) {
		double a1 = -ep.D/2;
		double b1 = ep.E/2;
				
		double f1 = sqrt(b1*b1+(1-a1)*(1-a1));
		double f2 = sqrt(b1*b1+a1*a1);
		double f3 = sqrt((1+b1)*(1+b1)+(1-a1)*(1-a1));
		double f4 = sqrt((1+b1)*(1+b1)+a1*a1);
				
		output1 = ((r-b1)>0)*((f1-r)>0)*(2*acos(b1/r)) + ((r-f1)>0)*((f2-r)>0)*(acos(b1/r)+asin((1-a1)/r)) + ((r-f2)>0)*((1+b1-r)>0)*(asin(a1/r)+asin((1-a1)/r));
		output1 += ((r-1-b1)>0)*((f3-r)>0)*(asin(a1/r)+asin((1-a1)/r)-2*acos((1+b1)/r)) + ((r-f3)>0)*((f4-r)>0)*( asin(a1/r)-acos((1+b1)/r) );
		output1 = output1*r;
	
		//printf("a b %.18f %.18f\n",a1,b1);
	}
	////////////////
	
	
	//printf("k = %lf \n",k);
	double output=0;
	for (int i = 0; i < count; i++)
		output += (gxx[vindex](theta[i][0]))*theta[i][1];
		
	//printf("output output1 %.18f %.18f\n",alpha[vindex]*output,output1);
	
	if(vindex==6 && ep.err_test == 1) return (output1-alpha[vindex]*output)*(output1-alpha[vindex]*output);
	
	return alpha[vindex]*output;
}
