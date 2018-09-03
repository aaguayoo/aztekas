/*
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"./Headers/main.h"

void RESTART()
{
   FILE *file;
   int i, j, k, idum;
   double dum;
   char line[100];
   
   //Initialize dt
   dt = 0.1; 
   
   file = fopen(restartfile,"r");
   
   // Skip first line
   idum = fscanf(file,"%s\n",line) ;

   // Read time
   idum = fscanf(file,"%lf\n",&time) ;

   printf("%i %f\n",idum,time) ;

   // Skip third line   
   idum = fscanf(file,"%s\n",line) ;

   // Read rest of file an initialize variables      
   if(dim == 1)
   {
     for(i = 3; i <= Nx1-3; i++)
     {
        idum = fscanf(file,"%lf %lf %lf %lf\n",&dum,\
               &U[c1(0,i)],&U[c1(1,i)],&U[c1(2,i)]);
     }   
   }
   if(dim == 2)
   {   
     for(i = 3; i <= Nx1-3; i++)
  	 {
        for(j = 3; j <= Nx2-3; j++)
        {
            idum = fscanf(file,"%lf %lf %lf %lf %lf %lf\n",&dum,&dum,\
                &U[c2(0,i,j)],&U[c2(1,i,j)],&U[c2(2,i,j)],&U[c2(3,i,j)]);
        }
     }
   }
   if(dim == 3)
   {
     for(i = 3; i <= Nx1-3; i++)
  	 {
        for(j = 3; j <= Nx2-3; j++)
        {
           for(k = 3; k <= Nx3-3; k++)
           {
               idum = fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf\n",\
                      &dum,&dum,&dum,\
                      &U[c3(0,i,j,k)],&U[c3(1,i,j,k)],\
                      &U[c3(2,i,j,k)],&U[c3(3,i,j,k)],&U[c3(4,i,j,k)]);
           }    
        }
     }   
   }

   // CALLING BOUNDARIES TO GET GHOST CELLS RIGHT
   BOUNDARIES(U);

   fclose(file);
}

void INITIAL()
{
   int n, i, j, k;
   
   //Initialize time
   time    = 0.0; 
   
   //Initialize dt
   dt      = 0.1; 
   

   if(dim == 1)
   {
/*      for(i = 0; i <= Nx1; i++)
      {
         if(X1[i] < 0.5)
         {
            U[c1(0,i)] = nl;
            U[c1(1,i)] = pl;
            U[c1(2,i)] = vx1l;
         }
         else
         {
            U[c1(0,i)] = nr;
            U[c1(1,i)] = pr;
            U[c1(2,i)] = vx1r;
         }
      }
*/   }
   else if(dim == 2)
   {
      ///////////////////////////
      //--------Jet------------//
      ///////////////////////////
      /*
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            if(fabs(X1[i]) <= r_jet)
            {
               if(fabs(X2[j]) <= z_jet)
               {
                  U[c2(0,i,j)] = n_jet;
                  U[c2(1,i,j)] = p_jet;
                  U[c2(2,i,j)] = vx1_jet;
                  U[c2(3,i,j)] = vx2_jet;
               }
               else
               {
                  U[c2(0,i,j)] = n_atm;
                  U[c2(1,i,j)] = p_atm;
                  U[c2(2,i,j)] = vx1_atm;
                  U[c2(3,i,j)] = vx2_atm;
               }
            }
            else
            {
               U[c2(0,i,j)] = n_atm;
               U[c2(1,i,j)] = p_atm;
               U[c2(2,i,j)] = vx1_atm;
               U[c2(3,i,j)] = vx2_atm;
            }
         }
      }
      */
      /////////////////////////////

      ///////////////////////////
      //-----Acretion----------//
      ///////////////////////////

		for(i = 0; i <= Nx1; i++)
		{	
			for(j = 0; j <= Nx2; j++)
			{
					U[c2(0,i,j)] = density_0;
					U[c2(1,i,j)] = pressure_0;
					U[c2(2,i,j)] = 0.0;
					U[c2(3,i,j)] = 0.0;
			}
		}
		
    int jmin;
    double r, R, z, dummy;

		for(i = 3; i <= Nx1; i++)
		{	
		  dummy = r_bou*r_bou - X1[i]*X1[i];
		  if (dummy >=0 ) 
		  {
		    jmin = sqrt(dummy)/dx2 + 4;
		  }  
		  else
		  {
		    jmin = 3;
		  }
		  		  
			for(j = jmin; j <= Nx2; j++)
			{
			  R = X1[i];
			  z = X2[j];
			  r = sqrt(R*R + z*z);
			  
				U[c2(0,i,j)] = density_0;
				U[c2(1,i,j)] = pressure_0;
				U[c2(2,i,j)] = velocity_0*(R/r);
				U[c2(3,i,j)] = velocity_0*(z/r);
			}
		}

      /////////////////////////////

      ///////////////////////////
      //--------Sedov----------//
      ///////////////////////////
      /*
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            if(pow(X1[i] - 0.5,2.0) + pow(X2[j] - 0.5,2.0) < pow(0.15,2.0))
            {
               U[c2(0,i,j)] = 1.0;
               U[c2(1,i,j)] = 1.0;
               U[c2(2,i,j)] = 0.0;
               U[c2(3,i,j)] = 0.0;
            }
            else
            {
               U[c2(0,i,j)] = 1.0;
               U[c2(1,i,j)] = 0.1;
               U[c2(2,i,j)] = 0.0;
               U[c2(3,i,j)] = 0.0;
            }
         }
      }
      */
      /////////////////////////////

      ///////////////////////////
      //-------Riemann---------//
      ///////////////////////////
      /*
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            if(X1[i] < 0.5)
            {
               U[c2(0,i,j)] = nl;
               U[c2(1,i,j)] = pl;
               U[c2(2,i,j)] = vx1l;
               U[c2(3,i,j)] = vx2l;
            }
            else
            {
               U[c2(0,i,j)] = nr;
               U[c2(1,i,j)] = pr;
               U[c2(2,i,j)] = vx1r;
               U[c2(3,i,j)] = vx2r;
            }
         }
      }
   */
   /////////////////////////////
   }
   else if(dim == 3)
   {
/*      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            for(k = 0; k <= Nx3; k++)
            {
               if(X1[i] < (x1max - x1min)/2.0)
               {
                  U[c3(0,i,j,k)] = nl;
                  U[c3(1,i,j,k)] = pl;
                  U[c3(2,i,j,k)] = vx1l;
                  U[c3(3,i,j,k)] = vx2l;
                  U[c3(4,i,j,k)] = vx3l;
               }
               else
               {
                  U[c3(0,i,j,k)] = nr;
                  U[c3(1,i,j,k)] = pr;
                  U[c3(2,i,j,k)] = vx1r;
                  U[c3(3,i,j,k)] = vx2r;
                  U[c3(4,i,j,k)] = vx3r;
               }
            }
         }
      }
*/   }
   else
   {
      printf("ERROR: Variable dim distinta de {1,2,3}");
   }
}
