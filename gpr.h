double gp11(double x, double y) {
    return ( x*x*y );
}

double gp22(double x, double y) {
    return ( y*y*y/3.0 );
}

double gp12(double x, double y) {
    return ( x*y*y/2.0 );
}

double gp13(double x, double y) {
    return ( x*y );
}

double gp23(double x, double y) {
    return ( y*y/2 );
}

double gp33(double x, double y) {
    return ( y );
}

typedef double (*gp)(double, double);
gp gpxx[7] = {&gp11, &gp22, &gp33, &gp23, &gp13, &gp12, &gp33};  

double gpr(double r, void * params)
{
	Conic ep = *(Conic *) params;
	
	double F = ep.F - r*r;
	int vindex = ep.index;
  
	double Mp11 = gsl_matrix_get(ep.M1,0,0);
	double Mp13 = gsl_matrix_get(ep.M1,0,2);
		
	//double a = sqrt((r*r-Mp33)/Mp11);
	//double b = sqrt((r*r-Mp33)/Mp22);
	//double k = sqrt(1-(b*b)/(a*a));
	
	double XC = -Mp13/Mp11;
	double XD2 = (XC*XC - F/Mp11);
	double XD = sqrt(XD2);
	double XD20 = (XC*XC - ep.F/Mp11);
	double prefac = r/(Mp11*XD); //careful when XD -> 0
	
	if(ep.debug == 1) printf("XD2=%.18f eps=%.18f F=%0.18f \n XD20 = %.18f prefac=%.18f\n",XD2,EPS,F,XD20,prefac);
	
	if (abs(XD2)<=EPS) 
	{
		//prefac = 1.0/sqrt(Mp11);
		//return INFINITY;
		return 0;
	}
		
  double alpha[7] = { f(r)*prefac,
                      f(r)*prefac,
                      f(r)*prefac,
                      f(r)*prefac,
                      f(r)*prefac,
                      f(r)*prefac,
                      V(r)*prefac };

	//printf("\na = %lf, b = %lf, k = %lf, alpha = %lf \n",a,b,k,alpha);
	
	//printf("Hello\n");
	
	if(r<=ep.rmin || r>=ep.rmax) return 0;
	
	//printf("Hello\n");
	
	//printf("The modified conic is (%g)x^2+(%g)y^2+(%g)=0\n",Mp11,Mp22,Mp33);
	
	gsl_complex csol[8];
	/* A * x^2 + D * x + F = 0 */
	gsl_poly_complex_solve_quadratic(ep.A, ep.D, F, &csol[0], &csol[1]);	
	gsl_poly_complex_solve_quadratic(ep.C, ep.B+ep.E, ep.A+ep.D+F, &csol[2], &csol[3]);	
	gsl_poly_complex_solve_quadratic(ep.A, ep.B+ep.D, ep.C+ep.E+F, &csol[4], &csol[5]);	
	gsl_poly_complex_solve_quadratic(ep.C, ep.E, F, &csol[6], &csol[7]);	
	if(ep.debug == 1) print_roots8(csol);
	
  for (int i = 0; i < 8; i++)
	{
    if(abs(GSL_REAL(csol[i]))<EPS) GSL_SET_REAL(&csol[i], 0);  
  }
	if(ep.debug == 1) print_roots8(csol);

	double rsol[24],drsol[24]; int count=0;
	for (int i = 0; i < 4; i++)
	{
		int loc = i%2;
		double val = 1.5-fabs(1.5-i);
		if(fabs(GSL_IMAG(csol[i*2]))<EPS) //real roots
		{
			//if( abs(GSL_REAL(csol[i*2])-GSL_REAL(csol[i*2+1])) > EPS_ROOT ) // roots are not identical
			{
				for (int j = 0; j < 2; j++)
				{
					if(abs(GSL_REAL(csol[i*2+j])-0.5)<(0.5+EPS)) //root is within the square
					{
						double xval = GSL_REAL(csol[i*2+j])*(loc==0) + val*(loc==1);
						double yval = GSL_REAL(csol[i*2+j])*(loc==1) + val*(loc==0);
						
						double sim_sum=0;
						for( int k=0; k<count; k++) //only for parabolic
						{
							double d[] = {xval-rsol[k*3],yval-rsol[k*3+1],0};
							sim_sum += (fabs(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])<EPS);
						}
						
						if(ep.debug == 1) printf("simsum = %lf\n", sim_sum);
						
						if(sim_sum<1)
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
	}

	//printf("count = %d\n",count);
	if(ep.debug == 1) mat_print1("rsol: \n",rsol,count,3);
	//if(count==1) count=0; //accounting for corner touch
  
	if(count==0) 
	//degenrate parabola (pair of lines) does not intersect square. 
	// In this case gij = 0
	{
		return 0; 
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
	
	//mat_print_gsl("dG: ",dG,count,1);
	
	gsl_matrix * Z = gsl_matrix_alloc (count, 3);
	
	gsl_blas_dgemm (CblasTrans, CblasTrans,
		1.0, zsol, ep.T,
		0.0, Z);
	
	if(ep.debug == 1) mat_print_gsl("Z: ",Z,count,3);
	
	double **Xsol = (double **)malloc(count * sizeof(double *));
	for (int i = 0; i < count; i++) Xsol[i] = (double *)malloc(3 * sizeof(double));
	
	for (int i = 0; i < count; i++)
	{
	  Xsol[i][0] = gsl_matrix_get(Z,i,0);
    Xsol[i][1] = gsl_matrix_get(Z,i,1);
  }
  
  if(ep.debug == 1) mat_print("Xsol (before sorting): ",Xsol,count,3);
  
  qsort(Xsol, count, sizeof Xsol[0], compare); //how am i sorting this?
  
  if(ep.debug == 1) mat_print("Xsol (before fixing column 3): ",Xsol,count,3);
  
  //if(( abs(XD20)<EPS || abs(F) < EPS )  && count==4) //exception when c(X,Y) is a perfect square.
  if( abs(XD20)<EPS && count==4) //exception when c(X,Y) is a perfect square.
  //if(count==4) //exception when c(X,Y) is a perfect square.
  {
			double d1 = (Xsol[1][0]-Xsol[0][0]), d2 = (Xsol[1][1]-Xsol[0][1]);
			double d3 = (Xsol[3][0]-Xsol[2][0]), d4 = (Xsol[3][1]-Xsol[2][1]);
			
			if( (d1*d1+d2*d2)<EPS ) Xsol[1][2] = 10;// 1st two rows are identical
			if( (d3*d3+d4*d4)<EPS ) Xsol[3][2] = 10;// last two rows are identical
			
			qsort(Xsol, count, sizeof Xsol[0], compare3);
			count-=2;
	}
	
	if(ep.debug == 1) mat_print("Xsol (after fixing repeat rows): ",Xsol,count,3);
    
  for (int i = 1; i < count; i+=2)
	{
		if( Xsol[i][1] >= Xsol[i-1][1] ) 
		{
			 Xsol[i][2] = 1; Xsol[i-1][2] = -1;
		}
		else
		{
			Xsol[i][2] = -1; Xsol[i-1][2] = 1;
		}
	}
		
	if(ep.debug == 1) mat_print("Final Xsol: ",Xsol,count,3);
	
	//mat_print(theta,count,2);
	//qsort(theta, count, sizeof theta[0], compare);
	////mat_print(theta,count,2);
	//if(theta[0][1]==1) theta[count-1][0]-=2*PI; //the is to ensure arcs are not discontinuous
	//mat_print(theta,count,2);
	
	//for testing parallel fibers
	double output1=0;
	if(vindex==6 && ep.err_test == 1) {
		double a = ep.E/2;
		double b = sqrt(ep.F - a*a);
		double f = sqrt(r*r-b*b);
		
		//printf("hello output 1");
		
		if((r*r-b*b)>0) {
		output1 = ((f-a)>0)*((1-f+a)>0)*(1-f+a) + ((a-f)>0)*((1+f-a)>0)*(1+f-a) + ((-f-a)>0)*((1+f+a)>0)*(1+f+a) + ((f+a)>0)*((1-f-a)>0)*(1-f-a);
		output1 = output1*r/f;
		
		//printf("hello output 2");
		}
	}
	////////////////
	
	//printf("k = %lf \n",k);
	double output=0;
		
	for (int i = 0; i < count; i++)
		output += (gpxx[vindex](Xsol[i][0],Xsol[i][1]))*Xsol[i][2];
	
	if(vindex==6 && ep.err_test == 1) 
	{
		double diff = 1.0*(output1-alpha[vindex]*output);
		//printf("ep_err_test r output-output1 %d %.18f %.18f %.18f %.18f \n",ep.err_test,r,alpha[vindex]*output,output1,diff*diff);
		return diff*diff;
	}
		
	return alpha[vindex]*output;
	
	//return output;
		 
	//return 0;
}
