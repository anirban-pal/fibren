void fibren(Fibren *params)
{
  double pq[9],P0[3][4];
    
  for(int i=0; i<9; i++) pq[i] = params->pq[i];  
  
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
    {
      P0[i][j] = params->P0[i][j];
    }
  
  double pmag = sqrt(pq[0]*pq[0]+pq[1]*pq[1]+pq[2]*pq[2]);
  double qmag = sqrt(pq[3]*pq[3]+pq[4]*pq[4]+pq[5]*pq[5]);
  //vec_print("pq: ",pq,9);
  
  gsl_matrix * PQ = gsl_matrix_alloc (3, 3);
	gsl_matrix * M = gsl_matrix_alloc (3, 3);
	
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
		gsl_matrix_set (PQ, i, j, pq[3*i+j]);
	}
	
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,
                  1.0, PQ, PQ,
                  0.0, M);
    
  double A = gsl_matrix_get (M, 0, 0); double C = gsl_matrix_get (M, 1, 1); double F = gsl_matrix_get (M, 2, 2);
	double B = 2*gsl_matrix_get (M, 0, 1); double D = 2*gsl_matrix_get (M, 0, 2); double E = 2*gsl_matrix_get (M, 1, 2);
  
  double det1 = B*B-4*A*C;
	double det2 = sqrt(B*B+(A-C)*(A-C));
	double det0 = -0.25*(C*D*D+A*E*E-B*E*D+F*det1);
	
	if( params->debug == 1) mat_print_gsl("M:",M,3,3);
	if( params->debug == 1) printf("det0 det1 det2 %.18f %.18f %.18f\n",det0,det1,det2);
	//if(det0==0) printf("\n Conic is degenerate \n"); //handle carefully
  
	if(abs(det1)>=EPS3) //regular ellipse, maybe relax this criterion to abs(det1)>EPS
	{
    if( params->debug == 1) printf("CASE: elliptic %.18f (det1) > %.18f(EPS3) \n",abs(det1),EPS3);
    
    double perp = C-A-det2;
		double COS_THETA = B/sqrt(B*B+perp*perp);
		double SIN_THETA = perp/sqrt(B*B+perp*perp);
		if(COS_THETA<0) {COS_THETA=-COS_THETA; SIN_THETA = -SIN_THETA;}
			
		if(B==0 && A<=C) {COS_THETA = 1.0; SIN_THETA = 0.0;}
    
		double xc = (2*C*D-B*E)/det1, yc = (2*A*E-B*D)/det1;
    
    double T0[3][3] = {
			{COS_THETA, SIN_THETA, (-xc*COS_THETA-yc*SIN_THETA)},
			{-SIN_THETA, COS_THETA, (xc*SIN_THETA-yc*COS_THETA)},
			{0, 0, 1} };
			
		double T0i[3][3] = {
			{COS_THETA, -SIN_THETA, xc},
			{SIN_THETA, COS_THETA, yc},
			{0, 0, 1} };
				
		gsl_matrix * T = gsl_matrix_alloc (3, 3);
		gsl_matrix * Ti = gsl_matrix_alloc (3, 3);
		gsl_matrix * TiT = gsl_matrix_alloc (3, 3);

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				gsl_matrix_set (T, i, j, T0[i][j]);
				gsl_matrix_set (Ti, i, j, T0i[i][j]);
			}
		
		if( params->debug == 1) mat_print_gsl("T:",T,3,3);
		
		gsl_matrix * MTi = gsl_matrix_alloc (3, 3);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
			  1.0, M, Ti,
			  0.0, MTi);
			  
		gsl_matrix * M1 = gsl_matrix_alloc (3, 3);
		gsl_blas_dgemm (CblasTrans, CblasNoTrans,
			  1.0, Ti, MTi,
			  0.0, M1);
			  
		if( params->debug == 1) mat_print_gsl("\n ellipse M1: ",M1,3,3);
		
		double Omega0[3][4] = {
			{0, 1, 1, 0},
			{0, 0, 1, 1},
			{1, 1, 1, 1} };
			
		gsl_matrix * Omega = gsl_matrix_alloc (3, 4);
		gsl_matrix * TOmega = gsl_matrix_alloc (3, 4);
		
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 4; j++)
			{
				gsl_matrix_set (Omega, i, j, Omega0[i][j]);
			}
		
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
			  1.0, T, Omega,
			  0.0, TOmega);	
		
		if( params->debug == 1)
		{
			printf("\nThe original conic is (%g)x^2+(%g)xy+(%g)y^2+(%g)x+(%g)y+(%g - r^2)=0\n",A,B,C,D,E,F);
			printf("The original domain is:"); mat_print_gsl("",Omega,3,4);
			
			printf("\nThe modified conic is (%g)x^2+(%g)y^2+(%g - r^2)=0\n",gsl_matrix_get(M1, 0, 0),gsl_matrix_get(M1, 1, 1),gsl_matrix_get(M1, 2, 2));
	    printf("The modified domain is:"); mat_print_gsl("",TOmega,3,4);
		}
		//printf("\n rmin given through %lf %lf \\",	-4*det0/det1, gsl_matrix_get(M1,2,2));
	  	  
	  gsl_complex rminc = gsl_complex_sqrt_real(gsl_matrix_get(M1,2,2)); //how about sqrt(Mp33) instead? same thing
		
		double p0[] = {P0[0][0],P0[1][0],P0[2][0]}; 
		double p1[] = {P0[0][1],P0[1][1],P0[2][1]};
		double q0[] = {P0[0][2],P0[1][2],P0[2][2]};
		double q1[] = {P0[0][3],P0[1][3],P0[2][3]};  
		
		double vminmags[4];
		vminmags[0] = pDistance(p0,q0,q1);
		vminmags[1] = pDistance(p1,q0,q1);
		vminmags[2] = pDistance(q0,p0,p1);
		vminmags[3] = pDistance(q1,p0,p1);
		
		double v1[] = {pq[6], pq[7], pq[8]};
		double v2[] = {pq[6]+pq[0], pq[7]+pq[1], pq[8]+pq[2]};
		double v3[] = {pq[6]+pq[3], pq[7]+pq[4], pq[8]+pq[5]};
		double v4[] = {pq[6]+pq[0]+pq[3], pq[7]+pq[1]+pq[4], pq[8]+pq[2]+pq[5]};
		
		double vmags[4];
		vmags[0] = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
		vmags[1] = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
		vmags[2] = sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
		vmags[3] = sqrt(v4[0]*v4[0] + v4[1]*v4[1] + v4[2]*v4[2]);
		
		double rmax = max(max(vmags[0],vmags[1]),max(vmags[2],vmags[3]));
		double rmin0 = min(min(vminmags[0],vminmags[1]),min(vminmags[2],vminmags[3]));
    
    double rmin;
    
    if( params->debug == 1) printf("xc = %.18f yc = %.18f\n",xc,yc);
    
    if ( xc*(xc-1) < 0 && yc*(yc-1) < 0 ) //points of nearest approach lie on the segments p(t1) and q(t2)
			rmin = GSL_REAL(rminc);
		else 
			rmin = rmin0;
        
    //double rmin = max(rmin0,GSL_REAL(rminc)); // this might need revision, rmin is quite complicated
    
		if( params->debug == 1) printf("rmin = %.18f rmax = %.18f\n",rmin,rmax);
		if (rmin<EPS2) printf("Error rmin = %lf is very close to zero, fibers intersect!! \n",rmin);
			  
    Conic ellipse_params;
    ellipse_params.A = A;    		ellipse_params.B = B;    		ellipse_params.C = C;
    ellipse_params.D = D;    		ellipse_params.E = E;    		ellipse_params.F = F;
    ellipse_params.det0 = det0; 	ellipse_params.det1 = det1; 	ellipse_params.det2 = det2;
    ellipse_params.xc = xc; 		ellipse_params.yc = yc;
    ellipse_params.rmin = rmin; 		ellipse_params.rmax = rmax;
    
    
    ellipse_params.M = gsl_matrix_alloc (3, 3);	gsl_matrix_memcpy(	ellipse_params.M, M );
    ellipse_params.T = gsl_matrix_alloc (3, 3);	gsl_matrix_memcpy(	ellipse_params.T, T );
    ellipse_params.M1 = gsl_matrix_alloc (3, 3);	gsl_matrix_memcpy(	ellipse_params.M1, M1 );
    
    ellipse_params.index = 6; //what does this do? it specifies which gij component to use
    
    //for(double r = (0.1-500e-3); r<=(0.1+500e-3); r+=1e-3)    
    //for(double r = (0.812-1e-3); r<=(1.2+1e-3); r+=1e-3)
    if( params->debug == 1)
    {
	    double spc = 1e-3;
	    //double r=0.4; printf("r gr %.18f %.18f\n",r,gr(r,&ellipse_params));
	    for(double r = (rmin-spc); r<=(rmax+spc); r+=spc)
	      printf("r gr %.18f %.18f\n",r,gr(r,&ellipse_params));
	  }
	  
	  //if(1>2)
	  {
	  ////////////
	  char * line = NULL;
    size_t len = 0;
    ssize_t read;
    fprintf(params->outfile,"rmin = %.6f rmax = %.6f\n",rmin,rmax);
    fprintf(params->outfile,"\t\t\t\t r \t\t\t\t\t\t\t\t g(r)\n");
    
	  while ((read = getline(&line, &len, params->grid)) != -1) {
        double r = atof(line);
        fprintf(params->outfile,"%.8f %.18f\n",r,gr(r,&ellipse_params));
    }
    /////////////
    
    gsl_function G;
    G.function = &gr;
    G.params = &ellipse_params;
    double result, error;
    
    //gsl_integration_workspace * w 
    //  = gsl_integration_workspace_alloc (100000);
      
    gsl_integration_cquad_workspace * w1
			= gsl_integration_cquad_workspace_alloc(100);
		
		size_t nev;
    
    double correct[7] = {1.0/3, 1.0/3, 1.0, 0.5, 0.5, 0.25, 1.0};
    double fgl[7];
    for(int i = 0 ; i < 7; i++) {
      clock_t begin = clock();
      
      ellipse_params.index = i;
      //gsl_integration_qags (&G, rmin, rmax+1.0, 1e-8, 1e-8, 100000,
      //            w, &result, &error);
      
      double pts[2]={rmin, rmax};
      //gsl_integration_qagp (&G, pts, 2, ABS_ERR, REL_ERR, 100000,
      //            w, &result, &error);
                  
      gsl_integration_cquad (&G, rmin, rmax, ABS_ERR, params->rel_err,
                  w1, &result, &error, &nev);
        
      clock_t end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        
      //printf ("result = % .18f error = %.10e Time taken = %lf sec\n", result,(result-correct[i]),time_spent);
      if( params->debug == 1) printf ("result = % .18f Time taken = %lf sec\n", result,time_spent);
      fgl[i] = result;
    }
    
    /////////////////error calculation
    if(params->errL2_flag == 1)
    {
			int i = 6;
			clock_t begin = clock();
      
      ellipse_params.index = i;
      ellipse_params.err_test = params->errL2_flag; //error testing flag
      
      //gsl_integration_qags (&G, rmin, rmax+1.0, 1e-8, 1e-8, 100000,
      //            w, &result, &error);
            
      double pts[2]={rmin, rmax};
      //gsl_integration_qagp (&G, pts, 2, ABS_ERR, REL_ERR, 100000,
        //          w, &result, &error);
                  
      gsl_integration_cquad (&G, rmin, rmax, 1e-12, REL_ERR,
                  w1, &result, &error, &nev);
              
      clock_t end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        
      //printf ("result = % .18f error = %.10e Time taken = %lf sec\n", result,(result-correct[i]),time_spent);
      if( params->debug == 1) printf ("result = % .18f Time taken = %lf sec\n", result,time_spent);
      params->errL2 = sqrt(result);
		}
		ellipse_params.err_test = 0;
		///////////////////////////////
     
    params->energy = fgl[6]*pmag*qmag; 
    ///////
    double H0[4][3] = {
      {-1, 0, 1},
      {1, 0, 0},
      {0, 1, -1},
      {0, -1, 0}, };
    
    gsl_matrix * H = gsl_matrix_alloc (4, 3);
    //gsl_matrix * HT = gsl_matrix_alloc (3, 4);
    
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++)
      {
        gsl_matrix_set (H, i, j, H0[i][j]);
        //gsl_matrix_set (HT, j, i, H0[i][j]);
      }
    
    //gsl_mat_print(H,4,3);
    //gsl_mat_print(HT,3,4);
    
    //double P0[3][4] = { {p0[0], p0[0]+p[0], q0[0], q0[0]+q[0]},
              //{p0[1], p0[1]+p[1], q0[1], q0[1]+q[1]},
              //{p0[2], p0[2]+p[2], q0[2], q0[2]+q[2]} };
    
    gsl_matrix * P = gsl_matrix_alloc (3, 4);
    
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++)
      {
        gsl_matrix_set (P, i, j, P0[i][j]);
      }
    
    //gsl_mat_print(P,3,4);	
      
      //////
      
    gsl_matrix * fg = gsl_matrix_alloc (3, 3);
    double fg0[3][3] = {	{fgl[0], fgl[5], fgl[4]},
              {fgl[5], fgl[1], fgl[3]},
              {fgl[4], fgl[3], fgl[2]}	};
    
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      gsl_matrix_set (fg, i, j, fg0[i][j]);
    
    gsl_matrix * HTi = gsl_matrix_alloc (4, 3);
    gsl_matrix * TiTHT = gsl_matrix_alloc (3, 4);
    
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
      1.0, H, Ti,
      0.0, HTi);
      
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++)
      {
        gsl_matrix_set (TiTHT, i, j, gsl_matrix_get(HTi,j,i) );
      }
      
    //gsl_mat_print(HTi,4,3);	
    //gsl_mat_print(TiTHT,3,4);	
      
    gsl_matrix * PHTi = gsl_matrix_alloc (3, 3);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
      1.0, P, HTi,
      0.0, PHTi);
      
    gsl_matrix * fgTiTHT = gsl_matrix_alloc (3, 4);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
      1.0, fg, TiTHT,
      0.0, fgTiTHT);
    
    gsl_matrix * F = gsl_matrix_alloc (3, 4);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
      -1.0, PHTi, fgTiTHT,
      0.0, F);
    
    if( params->debug == 1) mat_print_gsl("F: ",F,3,4);	
    
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++)
      {
        params->force[i+3*j] = gsl_matrix_get(F,i,j)*pmag*qmag;
      }
		}
  }
  //////////////////////////////////////////////////////////////////////
  else //degenerate parabola
  {
    //printf("CASE: degenerate parabola \n");
    if( params->debug == 1) printf("CASE: degenerate parabola %.18f (det1) < %.18f(EPS3) \n",abs(det1),EPS3);
    double k = -0.5*B/C;
    double kd = sqrt(k*k+1);
    double COS_PHI = k/kd;
		double SIN_PHI = -1/kd;
	
		double T0[3][3] = {
			{COS_PHI, SIN_PHI, 0},
			{-SIN_PHI, COS_PHI, 0},
			{0, 0, 1} };
		
		double T0i[3][3] = {
			{COS_PHI, -SIN_PHI, 0},
			{SIN_PHI, COS_PHI, 0},
			{0, 0, 1} };
			
		gsl_matrix * T = gsl_matrix_alloc (3, 3);
		gsl_matrix * Ti = gsl_matrix_alloc (3, 3);
	
		
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				gsl_matrix_set (T, i, j, T0[i][j]);
				gsl_matrix_set (Ti, i, j, T0i[i][j]);
			}
			
		if( params->debug == 1) mat_print_gsl("T: ",T,3,3);
			
		gsl_matrix * MTi = gsl_matrix_alloc (3, 3);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
			  1.0, M, Ti,
			  0.0, MTi);
				  
		gsl_matrix * M1 = gsl_matrix_alloc (3, 3);
		gsl_blas_dgemm (CblasTrans, CblasNoTrans,
			  1.0, Ti, MTi,
			  0.0, M1);
				  
		if( params->debug == 1) mat_print_gsl("\n parabola M1: ",M1,3,3);
			
		double Omega0[3][4] = {
			{0, 1, 1, 0},
			{0, 0, 1, 1},
			{1, 1, 1, 1} };
			
		gsl_matrix * Omega = gsl_matrix_alloc (3, 4);
		gsl_matrix * TOmega = gsl_matrix_alloc (3, 4);
		
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 4; j++)
			{
				gsl_matrix_set (Omega, i, j, Omega0[i][j]);
			}
		
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
			  1.0, T, Omega,
			  0.0, TOmega);	
			
		if( params->debug == 1)
		{
			printf("\nThe original conic is (%g)x^2+(%g)xy+(%g)y^2+(%g)x+(%g)y+(%g - r^2)=0\n",A,B,C,D,E,F);
			printf("The original domain is:"); mat_print_gsl("",Omega,3,4);
			
			printf("\nThe modified conic is (%g)x^2+(%g)x+(%g - r^2)=0\n",gsl_matrix_get(M1, 0, 0),gsl_matrix_get(M1, 0, 2),gsl_matrix_get(M1, 2, 2));
	    printf("The modified domain is:"); mat_print_gsl("",TOmega,3,4);
	  }
    
    double Mp11 = gsl_matrix_get(M1,0,0);
		double Mp13 = gsl_matrix_get(M1,0,2);
		
		//gsl_complex rminc = gsl_complex_sqrt_real(gsl_matrix_get(M1,2,2));
    gsl_complex rminc = gsl_complex_sqrt_real(F-Mp13*Mp13/Mp11);
			
		double v1[3] = {pq[6], pq[7], pq[8]};
		double v2[3] = {pq[6]+pq[0], pq[7]+pq[1], pq[8]+pq[2]};
		double v3[3] = {pq[6]+pq[3], pq[7]+pq[4], pq[8]+pq[5]};
		double v4[3] = {pq[6]+pq[0]+pq[3], pq[7]+pq[1]+pq[4], pq[8]+pq[2]+pq[5]};
		
		double vmags[4];
		vmags[0] = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
		vmags[1] = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
		vmags[2] = sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
		vmags[3] = sqrt(v4[0]*v4[0] + v4[1]*v4[1] + v4[2]*v4[2]);
		
		double rmax = max(max(vmags[0],vmags[1]),max(vmags[2],vmags[3]));
		double rmin0 = min(min(vmags[0],vmags[1]),min(vmags[2],vmags[3]));
    
    double rmin;
    
    double lineside[4] = {-0.5*E/C, k-0.5*E/C, k-1-0.5*E/C, -1-0.5*E/C};	
    bool check = 1;
    if( params->debug == 1) printf("closest approach equation (%lf)t1 - t2 - (%lf) = 0",k,0.5*E/C);
    if( params->debug == 1) vec_print("\nlineside: ",lineside,4);
    for(int i=1; i<4; i++)
    {
			check = check && ( lineside[i]*lineside[i-1]>0 );
		} 
    if (check == 1) //no intersection for closet approach calculation 
    {
			if( params->debug == 1) printf("No intersection\n");
			rmin = rmin0;
		}
		else
			rmin = GSL_REAL(rminc);
    
    //double rmin = max(rmin0,GSL_REAL(rminc)); // this might need revision, rmin is quite complicated
    
		//printf("rmin = %lf + (%lf)i\n",GSL_REAL(rminc),GSL_IMAG(rminc));
		if( params->debug == 1) printf("rmin = %lf + (%lf)i\n",rmin,0.0);
		if( params->debug == 1) printf("rmax = %lf + (%lf)i\n",rmax,0.0);
    
    Conic parabola_params;
    parabola_params.A = A;    		parabola_params.B = B;    		parabola_params.C = C;
    parabola_params.D = D;    		parabola_params.E = E;    		parabola_params.F = F;
    parabola_params.det0 = det0; 	parabola_params.det1 = det1; 	parabola_params.det2 = det2;
    //parabola_params.xc = xc; 		parabola_params.yc = yc;
    parabola_params.rmin = rmin; 		parabola_params.rmax = rmax;
    
    parabola_params.debug = params->debug;
    
    parabola_params.M = gsl_matrix_alloc (3, 3);	gsl_matrix_memcpy(	parabola_params.M, M );
    parabola_params.T = gsl_matrix_alloc (3, 3);	gsl_matrix_memcpy(	parabola_params.T, T );
    parabola_params.M1 = gsl_matrix_alloc (3, 3);	gsl_matrix_memcpy(	parabola_params.M1, M1 );
    
    parabola_params.index = 6; //what does this do? it specifies which gij component to use
    
    //rmin = 0.01;
    if( params->debug == 1)
    {
      //double r = 1.3; printf("r gpr %lf %lf\n",r,gpr(r,&parabola_params));
      for(double r = (rmin-1e-3); r<=(rmax+1e-3); r+=1e-3)
				printf("r gpr %lf %.18f\n",r,gpr(r,&parabola_params));
    }
    
    //if(1>2)  
    if(params->debug != 1)
    {
		////////////
	  char * line = NULL;
    size_t len = 0;
    ssize_t read;
    fprintf(params->outfile,"rmin = %.6f rmax = %.6f\n",rmin,rmax);
    fprintf(params->outfile,"\t\t\t\t r \t\t\t\t\t\t\t\t g(r)\n");
    
	  while ((read = getline(&line, &len, params->grid)) != -1) {
        double r = atof(line);
        fprintf(params->outfile,"%.8f %.18f\n",r,gpr(r,&parabola_params));
    }
    /////////////
    
		gsl_function G;
    G.function = &gpr;
    G.params = &parabola_params;
    double result, error;
    
    //gsl_integration_workspace * w 
    //  = gsl_integration_workspace_alloc (100000);
    
    gsl_integration_cquad_workspace * w1
			= gsl_integration_cquad_workspace_alloc(100);
    
    size_t nev;
    
    //double correct[7] = {1.0/3, 1.0/3, 1.0, 0.5, 0.5, 0.25, 1.0};
    double fgl[7];
    for(int i = 0 ; i < 7; i++) {
      clock_t begin = clock();
      
      parabola_params.index = i;
      //gsl_integration_qags (&G, rmin, rmax+1.0, 1e-8, 1e-8, 100000,
      //            w, &result, &error);
            
      double pts[2]={rmin, rmax};
      //gsl_integration_qagp (&G, pts, 2, ABS_ERR, REL_ERR, 100000,
        //          w, &result, &error);
                  
      gsl_integration_cquad (&G, rmin, rmax, ABS_ERR, params->rel_err,
                  w1, &result, &error, &nev);
              
      clock_t end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        
      //printf ("result = % .18f error = %.10e Time taken = %lf sec\n", result,(result-correct[i]),time_spent);
      if( params->debug == 1) printf ("result = % .18f Time taken = %lf sec\n", result,time_spent);
      fgl[i] = result;
    }
    
    ///////////////////error calculation
    if(params->errL2_flag == 1)
    {
			//gsl_integration_workspace * w 
				//= gsl_integration_workspace_alloc (100000);
						
			int i = 6;
			clock_t begin = clock();
      
      parabola_params.index = i;
      parabola_params.err_test = params->errL2_flag; //error testing flag
      
      //gsl_integration_qags (&G, rmin, rmax+1.0, 1e-8, 1e-8, 100000,
      //            w, &result, &error);
            
      double pts[2]={rmin, rmax};
      //gsl_integration_qagp (&G, pts, 2, ABS_ERR, REL_ERR, 100000,
        //         w, &result, &error);
                  
      gsl_integration_cquad (&G, rmin, rmax, 1e-12, REL_ERR,
                  w1, &result, &error, &nev);
              
      clock_t end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        
      //printf("rmin,rmax = %lf %lf error %.4e result = %.4e Time taken = %lf sec\n",rmin,rmax,error, result,time_spent);
      if( params->debug == 1) printf ("result = % .18f Time taken = %lf sec\n", result,time_spent);
      params->errL2 = sqrt(result);
		}
		parabola_params.err_test = 0;
		/////////////////////////////////
		    	
    params->energy = fgl[6]*pmag*qmag;
    ////////
    double H0[4][3] = {
      {-1, 0, 1},
      {1, 0, 0},
      {0, 1, -1},
      {0, -1, 0}, };
    
    gsl_matrix * H = gsl_matrix_alloc (4, 3);
    //gsl_matrix * HT = gsl_matrix_alloc (3, 4);
    
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++)
      {
        gsl_matrix_set (H, i, j, H0[i][j]);
        //gsl_matrix_set (HT, j, i, H0[i][j]);
      }
    
    //gsl_mat_print(H,4,3);
    //gsl_mat_print(HT,3,4);
    
    //double P0[3][4] = { {p0[0], p0[0]+p[0], q0[0], q0[0]+q[0]},
              //{p0[1], p0[1]+p[1], q0[1], q0[1]+q[1]},
              //{p0[2], p0[2]+p[2], q0[2], q0[2]+q[2]} };
    
    gsl_matrix * P = gsl_matrix_alloc (3, 4);
    
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++)
      {
        gsl_matrix_set (P, i, j, P0[i][j]);
      }
    
    //gsl_mat_print(P,3,4);	
      
      //////
      
    gsl_matrix * fg = gsl_matrix_alloc (3, 3);
    double fg0[3][3] = {	{fgl[0], fgl[5], fgl[4]},
              {fgl[5], fgl[1], fgl[3]},
              {fgl[4], fgl[3], fgl[2]}	};
    
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      gsl_matrix_set (fg, i, j, fg0[i][j]);
    
    gsl_matrix * HTi = gsl_matrix_alloc (4, 3);
    gsl_matrix * TiTHT = gsl_matrix_alloc (3, 4);
    
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
      1.0, H, Ti,
      0.0, HTi);
      
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++)
      {
        gsl_matrix_set (TiTHT, i, j, gsl_matrix_get(HTi,j,i) );
      }
      
    //gsl_mat_print(HTi,4,3);	
    //gsl_mat_print(TiTHT,3,4);	
      
    gsl_matrix * PHTi = gsl_matrix_alloc (3, 3);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
      1.0, P, HTi,
      0.0, PHTi);
      
    gsl_matrix * fgTiTHT = gsl_matrix_alloc (3, 4);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
      1.0, fg, TiTHT,
      0.0, fgTiTHT);
    
    gsl_matrix * F = gsl_matrix_alloc (3, 4);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
      -1.0, PHTi, fgTiTHT,
      0.0, F);
    
    if( params->debug == 1) mat_print_gsl("F: ",F,3,4);	
    
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++)
      {
        params->force[i+3*j] = gsl_matrix_get(F,i,j)*pmag*qmag;
      }
    }
  }
}


void fibren_discrete(Fibren *params, int n)
{
  double pq[9],P0[3][4];
    
  for(int i=0; i<9; i++) pq[i] = params->pq[i];  
  
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
    {
      P0[i][j] = params->P0[i][j];
    }
  
  double pmag = sqrt(pq[0]*pq[0]+pq[1]*pq[1]+pq[2]*pq[2]);
  double qmag = sqrt(pq[3]*pq[3]+pq[4]*pq[4]+pq[5]*pq[5]);
  
  gsl_vector * t = gsl_vector_alloc (3);
  
  gsl_vector * p0 = gsl_vector_alloc (3);
	gsl_vector * p = gsl_vector_alloc (3);
	gsl_vector * q0 = gsl_vector_alloc (3);
	gsl_vector * q = gsl_vector_alloc (3);
		
	gsl_vector * pi = gsl_vector_alloc (3);
	gsl_vector * qj = gsl_vector_alloc (3);
		
	for (int i = 0; i < 3; i++) {
		gsl_vector_set (p0, i, params->P0[i][0]);
		gsl_vector_set (p, i, params->P0[i][1]-params->P0[i][0]);
				
		gsl_vector_set (q0, i, params->P0[i][2]);
		gsl_vector_set (q, i, params->P0[i][3]-params->P0[i][2]);
	}
	
	gsl_vector_set_all(t,1);
	
	double H0[4][3] = {
		{-1, 0, 1},
		{1, 0, 0},
		{0, 1, -1},
		{0, -1, 0}, };
	
	gsl_matrix * H = gsl_matrix_alloc (4, 3);
		
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 3; j++)
		{
			gsl_matrix_set (H, i, j, H0[i][j]);
		}
	
  gsl_matrix * P = gsl_matrix_alloc (3, 4);
	
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 4; j++)
		{
			gsl_matrix_set (P, i, j, P0[i][j]);
		}
	
	double ti,tj;
	double U=0,r,dA = 1.0/(n*n);
	gsl_matrix * F = gsl_matrix_alloc (3, 4);
	gsl_matrix_set_zero(F);
	
	gsl_matrix * fg = gsl_matrix_alloc (3, 3);
	gsl_matrix_set_zero(fg);
	
	for(int i=0; i<=n; i++)
		for(int j=0; j<=n; j++)	{
				
				ti = (i*1.0-0.5)/n;
				tj = (j*1.0-0.5)/n;
			
				gsl_blas_dcopy(p0,pi); gsl_blas_daxpy(ti,p,pi);
				gsl_blas_dcopy(q0,qj); gsl_blas_daxpy(tj,q,qj);
				
				gsl_blas_daxpy(-1,qj,pi);
				r = gsl_blas_dnrm2(pi);
		
				U += V(r)/(n*n);
				
				gsl_vector_set(t,0,ti); gsl_vector_set(t,1,tj);
				
				gsl_blas_dger(f(r),t,t,fg);
		}
		
  gsl_matrix * PH = gsl_matrix_alloc (3, 3);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		1.0, P, H,
		0.0, PH);
		
	gsl_matrix * fgHT = gsl_matrix_alloc (3, 4);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,
		1.0, fg, H,
		0.0, fgHT);
	
	gsl_blas_dgemm( CblasNoTrans, CblasNoTrans,
		-1.0, PH, fgHT,
		1.0, F);
		      
	params->energy = U*pmag*qmag;
	for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      {
        params->force[i+3*j] = gsl_matrix_get(F,i,j)*dA*pmag*qmag;
      }
  
  
}
