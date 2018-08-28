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
#include"./Headers/limiters.h"
#include"./Headers/vector.h"

int INTEGRATION()
{
   int n, i, j, k;
   double uu[eq+1];
   double qq[eq+1];

   //Runge-Kutta 2th-Order and Piecewie Polynomial Reconstruction
   if(dim == 1)
   {
      funct_U2Q(Q,U);
      RK1D(U,Q,Q1,Q2,1);
      funct_Q2U(U1,Q1);
      BOUNDARIES(U1);
      RK1D(U1,Q,Q1,Q2,2);
      funct_Q2U(U,Q2);
      BOUNDARIES(U);
   }
   else if(dim == 2)
   {
      funct_U2Q(Q,U);
      RK2D(U,Q,Q1,Q2,1);
      funct_Q2U(U1,Q1);
      BOUNDARIES(U1);
      RK2D(U1,Q,Q1,Q2,2);
      funct_Q2U(U,Q2);
      BOUNDARIES(U);
   }
   else if(dim == 3)
   {
      funct_U2Q(Q,U);
      RK3D(U,Q,Q1,Q2,1);
      funct_Q2U(U1,Q1);
      BOUNDARIES(U1);
      RK3D(U1,Q,Q1,Q2,2);
      funct_Q2U(U,Q2);
      BOUNDARIES(U);
   }

   return 0;
}
