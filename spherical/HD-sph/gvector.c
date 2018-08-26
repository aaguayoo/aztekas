#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_G(double *a, double *uu)
{
   int i;
   double r;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   r = x1;
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
    
   for(i = 0; i <= eq; i++)
   {
      if(i == 0)
      {
         a[i] = n*(x1)*(sin(x2))*v;
      }
      else if(i == 1)
      {
         a[i] = ((K-1)*n*(x1)*(sin(x2))*v*pow(w,2.0)+(K-1)*n*(x1)*(sin(x2))*pow(v,3.0)+((K-1)*n*(x1)*(sin(x2))*pow(u,2.0)+2*K*p*(x1)*(sin(x2)))*v)/(2*K-2);
      }
      else if(i == 2)
      {
         a[i] = n*(x1)*(sin(x2))*u*v;
      }
      else if(i == 3)
      {
         a[i] = n*(x1)*(sin(x2))*pow(v,2.0)+p*(x1)*(sin(x2));
      }
      else if(i == 4)
      {
         a[i] = n*(x1)*(sin(x2))*v*w;
      }
   }

   return 0;
}
