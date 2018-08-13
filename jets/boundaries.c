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
#include<stdlib.h>
#include<string.h>
#include"./Headers/main.h"

double gtheta(double R, double z) ;

int BOUNDARIES(double *B)
{
   int n, i, j, k;

   if(dim == 1)
   {
      for(n = 0; n < eq; n++)
      {
         B[c1(n,3)] = B[c1(n,4)];
         B[c1(n,2)] = B[c1(n,4)];
         B[c1(n,1)] = B[c1(n,4)];
         B[c1(n,0)] = B[c1(n,4)];

         B[c1(n,Nx1-3)] = B[c1(n,Nx1-4)];
         B[c1(n,Nx1-2)] = B[c1(n,Nx1-4)];
         B[c1(n,Nx1-1)] = B[c1(n,Nx1-4)];
         B[c1(n,Nx1  )] = B[c1(n,Nx1-4)];
      }
   }
   else if(dim == 2)
   {
   ////////////////////////////////////////////////
   ///*---------INJECTION-BOUNDARY-------------*///
   ////////////////////////////////////////////////
   int imin, jmin;
	 double r, R, z, dummy;
   double pre, inv_r3;
   
   //imin = r_bou*sin(theta_f)/dx1 + 4 ;
   imin = 3;
	    	    
   for(i = imin; i <= Nx1; i++)
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
			  
         //pre = gtheta(R,z);
         inv_r3 = 1./pow(r,3) ;
			  
         U[c2(0,i,j)] = density_0;
         U[c2(1,i,j)] = pressure_0;
         U[c2(2,i,j)] = -alpha_B*velocity_0*R*(inv_r3 + 0.5*inv_s3);
         U[c2(3,i,j)] = -alpha_B*velocity_0*z*(inv_r3 - inv_s3);
			}
		}
   
	  // take care of ghost cells up, down, right and left
      for(n = 0; n < eq; n++)
      {
         for(i = 0; i <= Nx1-0; i++)
         {
       	  // set reflexion symmetry along y = 0 (x axis)
            B[c2(n,i,2)] =  B[c2(n,i,3)];
            B[c2(n,i,1)] =  B[c2(n,i,4)];
            B[c2(n,i,0)] =  B[c2(n,i,5)];

            B[c2(n,i,3)] = 0.5*(B[c2(n,i,2)] + B[c2(n,i,4)])

		  // copy values of Nx2-4 cell to all ghost cells above
            B[c2(n,i,Nx2-3)] = B[c2(n,i,Nx2-4)];
            B[c2(n,i,Nx2-2)] = B[c2(n,i,Nx2-4)];
            B[c2(n,i,Nx2-1)] = B[c2(n,i,Nx2-4)];
            B[c2(n,i,Nx2)]   = B[c2(n,i,Nx2-4)];
         }

         for(j = 0; j <= Nx2; j++)
         {
       	  // set reflexion symmetry along x = 0 (y axis)
            B[c2(n,2,j)] = B[c2(n,4,j)];
            B[c2(n,1,j)] = B[c2(n,5,j)];
            B[c2(n,0,j)] = B[c2(n,6,j)];

            B[c2(n,3,j)] = 0.5*(B[c2(n,2,j)] + B[c2(n,4,j)])

		  // copy values of Nx1-4 cell to all ghost cells to the right	
            B[c2(n,Nx1-3,j)] = B[c2(n,Nx1-4,j)];
            B[c2(n,Nx1-2,j)] = B[c2(n,Nx1-4,j)];
            B[c2(n,Nx1-1,j)] = B[c2(n,Nx1-4,j)];
            B[c2(n,Nx1  ,j)] = B[c2(n,Nx1-4,j)];
         }
      }

    // change sign for vz for ghost cells below x axis
      for(i = 0; i <= Nx1; i++)
      {
        B[c2(3,i,3)] = 0.0;
        B[c2(3,i,2)] = -B[c2(3,i,2)];
        B[c2(3,i,1)] = -B[c2(3,i,1)];
        B[c2(3,i,0)] = -B[c2(3,i,0)];
      }

    // change sign for vR for ghost cells left of y axis
      for(j = 0; j <= Nx2; j++)
      {
        B[c2(2,3,j)] = 0.0;
        B[c2(2,2,j)] = -B[c2(2,2,j)];
        B[c2(2,1,j)] = -B[c2(2,1,j)];
        B[c2(2,0,j)] = -B[c2(2,0,j)];
      }
   }

   else if(dim == 3)
   {
      for(n = 0; n < eq; n++)
      {
         for(i = 0; i <= Nx1-0; i++)
         {
            for(k = 0; k <= Nx3-0; k++)
            {
               B[c3(n,i,2,k)] =  B[c3(n,i,3,k)];
               B[c3(n,i,1,k)] =  B[c3(n,i,4,k)];
               B[c3(n,i,0,k)] =  B[c3(n,i,5,k)];

               B[c3(n,i,Nx2-2,k)] = B[c3(n,i,Nx2-3,k)];
               B[c3(n,i,Nx2-1,k)] = B[c3(n,i,Nx2-4,k)];
               B[c3(n,i,Nx2  ,k)] = B[c3(n,i,Nx2-5,k)];
            }
         }

         for(j = 0; j <= Nx2-0; j++)
         {
            for(k = 0; k <= Nx3-0; k++)
            {
               B[c3(n,2,j,k)] = B[c3(n,3,j,k)];
               B[c3(n,1,j,k)] = B[c3(n,4,j,k)];
               B[c3(n,0,j,k)] = B[c3(n,5,j,k)];

               B[c3(n,Nx1-2,j,k)] = B[c3(n,Nx1-3,j,k)];
               B[c3(n,Nx1-1,j,k)] = B[c3(n,Nx1-4,j,k)];
               B[c3(n,Nx1  ,j,k)] = B[c3(n,Nx1-5,j,k)];
            }
         }

         for(i = 0; i <= Nx1-0; i++)
         {
            for(j = 0; j <= Nx2-0; j++)
            {
               B[c3(n,i,j,2)] = B[c3(n,i,j,3)];
               B[c3(n,i,j,1)] = B[c3(n,i,j,4)];
               B[c3(n,i,j,0)] = B[c3(n,i,j,5)];

               B[c3(n,i,j,Nx3-2)] = B[c3(n,i,j,Nx3-3)];
               B[c3(n,i,j,Nx3-1)] = B[c3(n,i,j,Nx3-4)];
               B[c3(n,i,j,Nx3  )] = B[c3(n,i,j,Nx3-5)];
            }
         }
      }
   }

   return 0;
}

double gtheta(double R, double z) 
{
   double theta;   
   
   if( rho_boundary == 1)
   // constant profile
   {
      return 1.0;
   }
   
   if(R > 0 && z >= 0)
   {
      theta = atan(z/R);
   }
   else if(R == 0 && z > 0)
   {
      theta = M_PI_2;
   }
   else if(R < 0)
   {
      theta = atan(z/R) + M_PI;
   }
   else if(R == 0 && z < 0)
   {
      theta = 3.*M_PI_2;
   }
   else if(R >  0 && z < 0)
   {
      theta = atan(z/R) + 2.*M_PI;
   }			  

   if( rho_boundary == 2)
   // gaussian profile for the density at the injection boundary
   {
      return exp(-pow(theta/(sqrt(2.0)*theta_0),2.0));
   }
   else if( rho_boundary == 3)
   // step profile for the density at the injection boundary
   {
      return 1.0/(exp((theta - theta_0)/delta_theta) + 1.);
   }   
   
   return 1.0;
}
