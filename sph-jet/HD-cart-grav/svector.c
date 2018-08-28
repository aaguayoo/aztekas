#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_S(double *a, double *uu)
{
   int i;
   double r;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   r = sqrt(x1*x1 + x2*x2);
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
    
   for(i = 0; i <= eq; i++)
   {
      if(i == 0)
      {
         a[i] = 0;
      }
      else if(i == 1)
      {
         a[i] = -(n*w*x3+n*v*x2+n*u*x1)/pow(r,3.0);
      }
      else if(i == 2)
      {
         a[i] = -(n*x1)/pow(r,3.0);
      }
      else if(i == 3)
      {
         a[i] = -(n*x2)/pow(r,3.0);
      }
      else if(i == 4)
      {
         a[i] = -(n*x3)/pow(r,3.0);
      }
   }

   return 0;
}