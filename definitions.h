#define EPS 1e-12
#define EPS2 1e-8
#define EPS3 1e-10
#define PI 3.141592653589793238462643383279502884197169399375105820974
#define EPS_ROOT 1e-6  // This number is very important, // roots are not identical
#define max(x,y) ((x) >= (y)) ? (x) : (y)
#define min(x,y) ((x) <= (y)) ? (x) : (y)
#define ABS_ERR 0
#define REL_ERR 1e-12

//int compare ( const void *pa, const void *pb ) {
    //const double *a = *(const double **)pa;
    //const double *b = *(const double **)pb;
    //return a[0] < b[0] ? -1 : a[0] == b[0] ? 0 : 1;
//}

int compare ( const void *pa, const void *pb ) {
    const double *a = *(const double **)pa;
    const double *b = *(const double **)pb;
    if( fabs(a[0] - b[0])<EPS )
        return a[1]>b[1];
    else
        return a[0]>b[0];
}

int compare3 ( const void *pa, const void *pb ) {
    const double *a = *(const double **)pa;
    const double *b = *(const double **)pb;
	return a[2]>b[2];
}

void mat_print_gsl(const char *text, gsl_matrix * M,int nrow, int ncol)
{
	printf ("%s",text);
  for (int i = 0; i < nrow; i++)	/* OUT OF RANGE ERROR */
	{
		printf ("\n");
		for (int j = 0; j < ncol; j++)
			printf ("%lf ", gsl_matrix_get (M, i, j));
		
	}
	printf ("\n");
}

void vec_print(const char *text,double *V,int nel)
{
	printf ("%s",text);
  for (int i = 0; i < nel; i++)	/* OUT OF RANGE ERROR */
	{
		printf ("%lf ", V[i]);
  }
	printf ("\n");
}

void mat_print(const char *text,double ** M,int nrow, int ncol)
{
	printf ("%s",text);
	for (int i = 0; i < nrow; i++)	/* OUT OF RANGE ERROR */
	{
		printf ("\n");
		for (int j = 0; j < ncol; j++)
			printf ("%.18f ", M[i][j]);
		
	}
	printf ("\n");
}

void mat_print1(const char *text, double * M,int nrow, int ncol)
{
	printf ("%s",text);
  for (int i = 0; i < nrow; i++)	/* OUT OF RANGE ERROR */
	{
		printf ("\n");
		for (int j = 0; j < ncol; j++)
			printf ("%lf ", M[i*ncol+j]);
		
	}
	printf ("\n");
}

void print_roots8(gsl_complex r[])
{
	for(int i=0;i<8;i++)
		printf("%g + (%g)i\n", GSL_REAL(r[i]), GSL_IMAG(r[i]));
	printf ("\n");
}

double pDistance(double x[], double x1[], double x2[]) 
{
  double A[] = {x[0]-x1[0], x[1]-x1[1], x[2]-x1[2]};
  double C[] = {x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2]};
  
  //var A = x - x1;
  //var B = y - y1;
  //var C = x2 - x1;
  //var D = y2 - y1;

  double dot = A[0]*C[0] + A[1]*C[1] + A[2]*C[2];
  double len_sq = C[0]*C[0] + C[1]*C[1] + C[2]*C[2];
  double param = -1.0;
  
  //var dot = A * C + B * D;
  //var len_sq = C * C + D * D;
  //var param = -1;
    
  if (len_sq != 0) //in case of 0 length line
      param = dot / len_sq;

  double xx[3];
  //var xx, yy;

  if (param < 0) {
    for (int i=0; i<3; i++) xx[i] = x1[i];
  }
  else if (param > 1) {
    for (int i=0; i<3; i++) xx[i] = x2[i];
  }
  else {
    for (int i=0; i<3; i++) xx[i] = x1[i]+param*C[i];
  }

  double dx[] = {x[0]-xx[0],x[1]-xx[1],x[2]-xx[2]}; 
  //var dx = x - xx;
  //var dy = y - yy;
  return sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
}

typedef struct Conicparams {
	
  double A,B,C,D,E,F;
	double det0,det1,det2;
	double xc,yc;
	double rmin,rmax;
	
	gsl_matrix * M;
	gsl_matrix * T;
	gsl_matrix * M1;
	
	int debug = 0;
	int index;
	int err_test = 0;
	
} Conic;

typedef struct Fibrenparams {
	
  double pq[9],P0[3][4];
	double energy = 0, errL2 = 0;
  double force[12] = { 0 };
  int debug = 0;
  FILE *grid,*outfile;
	
	double rel_err=1e-8;
  int errL2_flag = 0;
    
} Fibren;
