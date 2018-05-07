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

int BOUNDARIES(double *B)
{
	int n, i, j, k;

	if(dim == 1)
	{	
///////////////////////////
//-----Shock Tube--------//
///////////////////////////
/*		for(n = 0; n < eq; n++)
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
*/

///////////////////////////
//-----Bondi-------------//
///////////////////////////
        int bound = 5;
        
        // INNER BOUNDARY
// constant		
		for(i = 0; i < bound; i++)
		{
			B[c1(0,i)] = B[c1(0,10)];
			B[c1(1,i)] = B[c1(1,10)];
			B[c1(2,i)] = B[c1(2,10)];
		}

// second order extrapolation
/*
        double r, r1, r2, r3, rho1, rho2, rho3, P1, P2, P3, vel1, vel2, vel3;
        double r1r2, r1r3, r2r3, fac1, fac2, fac3;

        r1 = X1[bound]; 
        r2 = X1[bound + 1];
        r3 = X1[bound + 2];
        
        r1r2 = r1 - r2;
        r1r3 = r1 - r3;
        r2r3 = r2 - r3;
                        
        rho1 = B[c1(0,bound)]; 
        rho2 = B[c1(0,bound + 1)];
        rho3 = B[c1(0,bound + 2)];
                        
        P1 = B[c1(1,bound)]; 
        P2 = B[c1(1,bound + 1)]; 
        P3 = B[c1(1,bound + 2)]; 
                        
        vel1 = B[c1(2,bound)];
        vel2 = B[c1(2,bound + 1)];
        vel3 = B[c1(2,bound + 2)];
                		
		for(i = 0; i < bound; i++)
		{
		    r = X1[i];
		    
		    fac1 = (r - r2)*(r - r3)/(r1r2*r1r3);
		    fac2 = -(r - r1)*(r - r3)/(r1r2*r2r3);
		    fac3 = (r - r1)*(r - r2)/(r2r3*r1r3);
		    		    		    
			B[c1(0,i)] = fac1*rho1 + fac2*rho2 + fac3*rho3;
  			B[c1(1,i)] = fac1*P1 + fac2*P2 + fac3*P3;
	  		B[c1(2,i)] = fac1*vel1 + fac2*vel2 + fac3*vel3;
		}
*/
        // OUTER BOUNDARY
/*		for(i = 0; i < 4; i++)
		{
			B[c1(0,Nx1-i)] = density_0;
			B[c1(1,Nx1-i)] = pressure_0;
			B[c1(2,Nx1-i)] = velocity_0;
		}
*/
	}
	else if(dim == 2)
	{
		for(n = 0; n < eq; n++)
		{
			for(i = 0; i <= Nx1-0; i++)
			{
				B[c2(n,i,2)] =  B[c2(n,i,4)];
				B[c2(n,i,1)] =  B[c2(n,i,5)];
 				B[c2(n,i,0)] =  B[c2(n,i,6)];
	
				B[c2(n,i,Nx2-3)] = B[c2(n,i,Nx2-4)];
				B[c2(n,i,Nx2-2)] = B[c2(n,i,Nx2-4)];
				B[c2(n,i,Nx2-1)] = B[c2(n,i,Nx2-4)];
				B[c2(n,i,Nx2)]   = B[c2(n,i,Nx2-4)];

				B[c2(n,i,3)] = 0.5*(B[c2(n,i,2)] + B[c2(n,i,4)]);
			}
		
			for(j = 0; j <= Nx2; j++)
			{
				B[c2(n,2,j)] = B[c2(n,4,j)];
				B[c2(n,1,j)] = B[c2(n,5,j)];
 				B[c2(n,0,j)] = B[c2(n,6,j)];
	
				B[c2(n,Nx1-3,j)] = B[c2(n,Nx1-4,j)];
				B[c2(n,Nx1-2,j)] = B[c2(n,Nx1-4,j)];
				B[c2(n,Nx1-1,j)] = B[c2(n,Nx1-4,j)];
				B[c2(n,Nx1  ,j)] = B[c2(n,Nx1-4,j)];

				B[c2(n,3,j)] = 0.5*(B[c2(n,2,j)] + B[c2(n,4,j)]);
			}
		}

		////////////////////////////////////////////////
		///*---------JET-BOUNDARY-------------------*///
		////////////////////////////////////////////////
/*
		for(j = 0; j <= Nx2; j++)
		{
			B[c2(2,3,j)] = -B[c2(2,4,j)];
			B[c2(2,2,j)] = -B[c2(2,5,j)];
			B[c2(2,1,j)] = -B[c2(2,6,j)];
 			B[c2(2,0,j)] = -B[c2(2,7,j)];
		}
		
		for(i = 0; i <= Nx1; i++)
		{
			for(j = 0; j <= Nx2; j++)
			{
				if(fabs(X1[i]) <= r_jet)
				{
					if(fabs(X2[j]) <= z_jet)
					{
						B[c2(0,i,j)] = n_jet;
						B[c2(1,i,j)] = p_jet;
						B[c2(2,i,j)] = vx1_jet;
						B[c2(3,i,j)] = vx2_jet + 0.29*fabs(sin(2*PI*time/5.0));
					}
				}	
			}
		}
*/		
		////////////////////////////////////////////////
		////////////////////////////////////////////////

		////////////////////////////////////////////////
		///*---------ACRETION-BOUNDARY--------------*///
		////////////////////////////////////////////////

		for(i = 0; i <= Nx1; i++)
		{
			B[c2(3,i,2)] = -B[c2(3,i,4)];
			B[c2(3,i,1)] = -B[c2(3,i,5)];
 			B[c2(3,i,0)] = -B[c2(3,i,6)];

			B[c2(3,i,3)] = 0.5*(B[c2(3,i,2)] + B[c2(3,i,4)]);
		}

		for(j = 0; j <= Nx2; j++)
		{
			B[c2(2,2,j)] = -B[c2(2,4,j)];
			B[c2(2,1,j)] = -B[c2(2,5,j)];
 			B[c2(2,0,j)] = -B[c2(2,6,j)];

			B[c2(2,3,j)] = 0.5*(B[c2(2,2,j)] + B[c2(2,4,j)]);
		}

		for(i = 0; i <= Nx1; i++)
		{	
			for(j = 0; j <= Nx2; j++)
			{
				double R = sqrt(X1[i]*X1[i] + X2[j]*X2[j]);
				double r  = X1[i];
				double z  = X2[j];
				double rho, pre, vr, vz;
				double theta, th_0, gtheta;
				
				th_0 = PI*theta_0; // rescale theta_0 in units of pi

				if(z >= r_int && z < (r_int + dx2) && r == 0.0)
				{
					rho = U[c2(0,i,j)];
					pre = U[c2(1,i,j)];
					vr  = U[c2(2,i,j)];
					vz  = U[c2(3,i,j)];
				}

				if(R < r_int)
				{
					U[c2(0,i,j)] = rho;
					U[c2(1,i,j)] = pre;
					U[c2(2,i,j)] = vr;
					U[c2(3,i,j)] = vz;
				}
				else if(R > r_bou)
				{
					if(R > 0 && z >= 0)
					{
						theta = atan(z/R);
					}
					else if(R == 0 && z > 0)
					{
						theta = 0.5*PI;
					}
					else if(R < 0)
					{
						theta = atan(z/R) + PI;
					}
					else if(R == 0 && z < 0)
					{
						theta = 1.5*PI;
					}
					else if(R >  0 && z < 0)
					{
						theta = atan(z/R) + 2.*PI;
					}
					
					//gtheta = 1.*exp(-pow(theta/(sqrt(2.0)*th_0),2.0));
					gtheta = 1./(exp((theta - th_0)/delta_theta) + 1.);

					U[c2(0,i,j)] = n_1*gtheta;
					U[c2(1,i,j)] = pow(U[c2(0,i,j)],K)/K;
					U[c2(2,i,j)] = -v_1*(r/R);
					U[c2(3,i,j)] = -v_1*(z/R);
				}
			}
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

