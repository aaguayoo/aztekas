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

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

int AMATRIX1D(double *u, vec_ *v, int *I)
{
   int m, n;
   int i;

   x1  = X1[I[0]];
   x2  = 0;
   x3  = 0;

   funct_S(v->S,u);

   return 0;
}

int AMATRIX2D(double *u, vec_ *v, int *I)
{
   int m, n;

   x1  = X1[I[0]];
   x2  = X2[I[1]];
   x3  = 0;

   funct_S(v->S,u);

   return 0;
}

int AMATRIX3D(double *u, vec_ *v, int *I)
{
   int m, n;

   x1  = X1[I[0]];
   x2  = X2[I[1]];
   x3  = X3[I[2]];

   funct_Q(v->Q,u);
   funct_S(v->S,u);

   return 0;
}

///////////////////////////////////////////////////////////////////////////

int VECTOR(int pm, char flux, lim_ *l, flx_ *f, int *I)
{
   int n;
   double *u, lr, ll;
   double up[eq+1];
   double um[eq+1];
   double Np, Pp, vx1p=0, vx2p=0, vx3p=0;
   double Nm, Pm, vx1m=0, vx2m=0, vx3m=0;
   double x1p=0, x2p=0, x3p=0;
   double x1m=0, x2m=0, x3m=0;
   double dp[3];
   double dm[3];

   x1 = X1[I[0]];

   if(dim == 2)
   {
      x2 = X2[I[1]];
   }
   else if(dim == 3)
   {
      x2 = X2[I[1]];
      x3 = X3[I[2]];
   }

   if(pm == 1)
   {
      switch(flux)
      {
         case 'f':
            u  = l->ux1p;
            x1 = X1p[I[0]];
         break;

         case 'g':
            u = l->ux2p;
            x2 = X2p[I[1]];
         break;

         case 'h':
            u = l->ux3p;
            x3 = X3p[I[2]];
         break;
      }
   }
   else
   {
      switch(flux)
      {
         case 'f':
            u = l->ux1m;
            x1 = X1m[I[0]];
         break;

         case 'g':
            u = l->ux2m;
            x2 = X2m[I[1]];
         break;

         case 'h':
            u = l->ux3m;
            x3 = X3m[I[2]];
         break;
      }
   }

   for(n = 0; n < eq; n++)
   {
      up[n] = u[1*eq + n];
      um[n] = u[0*eq + n];
   }

   funct_Q(f->qp,up);
   funct_Q(f->qm,um);

   switch(flux)
   {
      case 'f':
         funct_F(f->fp,up);
         funct_Dm(dp,up);

         funct_F(f->fm,um);
         funct_Dm(dm,um);
      break;

      case 'g':
         funct_G(f->fp,up);
         funct_Dn(dp,up);

         funct_G(f->fm,um);
         funct_Dn(dm,um);
      break;

      case 'h':
         funct_H(f->fp,up);
         funct_Do(dp,up);

         funct_H(f->fm,um);
         funct_Do(dm,um);
      break;
   }

   lr = max(dp[0],dp[1]);
   lr = max(lr,dp[2]);
   ll = max(dm[0],dm[1]);
   ll = max(ll,dm[2]);

   f->lp = max(lr,ll);

   lr = min(dp[0],dp[1]);
   lr = min(lr,dp[2]);
   ll = min(dm[0],dm[1]);
   ll = min(ll,dm[2]);

   f->lm = min(lr,ll);

   return 0;
}

///////////////////////////////////////////////////////////////////////////
