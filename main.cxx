// COMPILE WITH: gcc main.cxx -o main -lm -lgsl -lgslcblas
// RUN WITH: ./main <qtheta> <qphi> <qmidx> <qmidy> <qmidz> <gridfile>
/*
 * main.cxx
 * 
 * Copyright 2020 Anirban Pal, apal@wtamu.edu
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#include "headers.h"

int main(int argc, char **argv)
{
	double qtheta = atof(argv[1])*PI/180;  
	double qphi = atof(argv[2])*PI/180;  
	
	double qmid[3] = { atof(argv[3]), atof(argv[4]), atof(argv[5]) };
	double hdir[3] = { 0.5*cos(qphi)*sin(qtheta), 0.5*sin(qphi)*sin(qtheta), 0.5*cos(qtheta) };
  
  double p0[] = {-0.5, 0.0, 0.0};	double p[] = {1.0, 0.0, 0.0};
	double q0[] = {qmid[0] - hdir[0], qmid[1] - hdir[1], qmid[2] - hdir[2]};	double q[] = {2*hdir[0], 2*hdir[1], 2*hdir[2]};
	
  double pq[] = { p[0], p[1], p[2], -q[0], -q[1], -q[2], p0[0]-q0[0], p0[1]-q0[1], p0[2]-q0[2] };
  double P0[3][4] = { {p0[0], p0[0]+p[0], q0[0], q0[0]+q[0]},
                      {p0[1], p0[1]+p[1], q0[1], q0[1]+q[1]},
                      {p0[2], p0[2]+p[2], q0[2], q0[2]+q[2]}  };
  //vec_print("p0: ",p0,3);   vec_print("q0: ",q0,3); vec_print("p: ",p,3);   vec_print("q: ",q,3);   
  //vec_print("pq: ",pq,9);
  
  Fibren params; 
  for(int i=0; i<9; i++) 
    params.pq[i]=pq[i];
    
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
    {
      params.P0[i][j] = P0[i][j];
    }
  //vec_print("pq: ",params.pq,9);
  //params.debug = 1;
  //{vec_print("p0: ",p0,3);   vec_print("q0: ",q0,3); vec_print("p: ",p,3);   vec_print("q: ",q,3);   }
  
  if (params.debug == 1) {vec_print("p0: ",p0,3);   vec_print("q0: ",q0,3); vec_print("p: ",p,3);   vec_print("q: ",q,3);   }
  
  char outfilename[100]="outfile.gr";
  params.grid = fopen(argv[6],"r");
  params.outfile = fopen(outfilename,"w");  
	
	if(1) { //default option 
		params.rel_err = REL_ERR;
			
		clock_t begin1 = clock();
		fibren(&params); //calling the function fibren to compute the energy and forces using pq
	  clock_t end1 = clock();
	  
	  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
	  printf ("Energy_adaptive: %.18f rel_error: %.2e Time taken = %lf sec\n",params.energy, fabs(params.energy-1), time_spent1); //printing computed energy and forces 
	}
	
	
	if(0) { //option for computing errL2 for parallel and perpendicular cases. This will only work for cases referred to in the paper. Make sure V(r) = f(r) = 1 in potential.h
		// for parallel fibers, run ./main 90 0 0.0 0.6 0.0 gridfile.dat
		// for coplanar perpendicular fibers, run ./main 0 0 0.0 0.0 1.0 gridfile.dat
		params.errL2_flag = 1;		
		params.rel_err = REL_ERR;
			
		clock_t begin1 = clock();
		fibren(&params); //calling the function fibren to compute the energy and forces using pq
	  clock_t end1 = clock();
	  
	  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
	  printf ("Energy_adaptive: %.18f rel_error: %.2e errorL2: %.2e Time taken = %lf sec\n",params.energy, fabs(params.energy-1), params.errL2, time_spent1); //printing computed energy and forces 
	}
	
	if(0) { //option to run performance comparison with discrete scheme
		for (double e=2; e<=15; e++)
		{	
			params.rel_err = pow(10,-e);
			
			clock_t begin1 = clock();
			fibren(&params); //calling the function fibren to compute the energy and forces using pq
		  clock_t end1 = clock();
		  
		  double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
		  printf ("Energy_adaptive: %.18f rel_error(assigned): %.2e Time taken = %lf sec\n",params.energy,params.rel_err,time_spent1); //printing computed energy and forces 
		}
	  //printf ("Energy: %.18f\n",params.energy); //printing computed energy and forces 
	  //mat_print1("Forces: ",params.force,4,3);
	  
	  
	  for(double e = 0; e < 5; e+=0.5)
	  {
		  int n = (int) pow(10,e);
		  
		  clock_t begin2 = clock();
		  fibren_discrete(&params,n); //calling the function fibren discrete 
		  
		  clock_t end2 = clock();
		  double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
		  
		  printf ("Energy_discrete: %.18f n: %d Time taken = %lf sec\n",params.energy,n,time_spent2); //printing computed energy and forces 
		  //mat_print1("Forces: ",params.force,4,3);
		}
	}
	
  fclose(params.grid);
  fclose(params.outfile);
  
	return 0;
}

