//INCLUDE REGULAR C HEADERS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

//INCLUDE GSL HEADERS
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>

//INCLUDE LOCAL HEADERS
#include "definitions.h"
#include "potential.h"
#include "gr.h"
#include "gpr.h"
#include "fibren.h"
