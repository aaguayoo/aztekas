#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_F(double *a, double *uu)
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
         a[i] = g11*n*u;
      }
      else if(i == 1)
      {
         a[i] = ((g11*K-g11)*n*u*pow(w,2.0)+(g11*K-g11)*n*u*pow(v,2.0)+(g11*K-g11)*n*pow(u,3.0)+2*g11*K*p*u)/(2*K-2);
      }
      else if(i == 2)
      {
         a[i] = g11*n*pow(u,2.0)+g11*p;
      }
      else if(i == 3)
      {
         a[i] = g11*n*u*v;
      }
      else if(i == 4)
      {
         a[i] = g11*n*u*w;
      }
   }

   return 0;
}
