#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_S(double *a, double *uu)
{
   int i;
   double r;
   double n, p, u=0, v=0, w=0;
   double E;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   r = sqrt(x1*x1 + x2*x2);
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
   E = 0.5*n*(u*u + v*v + w*w) + p/(K-1);
    
   for(i = 0; i <= eq; i++)
   {
      if(i == 0)
      {
         a[i] = -(n*u)/x1;
      }
      else if(i == 1)
      {
         a[i] = -(x1*n*v*x2+((p+E)*pow(r,3.0)+pow(x1,2.0)*n)*u)/(x1*pow(r,3.0));
      }
      else if(i == 2)
      {
         a[i] = -(n*pow(r,3.0)*pow(u,2.0)+pow(x1,2.0)*n)/(x1*pow(r,3.0));
      }
      else if(i == 3)
      {
         a[i] = -(x1*n*x2+n*pow(r,3.0)*u*v)/(x1*pow(r,3.0));
      }
      else if(i == 4)
      {
         a[i] = -(n*u*w)/x1;
      }
   }

   return 0;
}
