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
   if(dim == 1){x2 = M_PI_2;};
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
   E = 0.5*n*(u*u + v*v + w*w) + p/(K-1);
    
   for(i = 0; i <= eq; i++)
   {
      if(i == 0)
      {
         a[i] = -(2*n*sin(x2)*u*x1+n*x1*cos(x2)*v)/(x1*sin(x2)*x1);
      }
      else if(i == 1)
      {
         a[i] = -((2*p+2*E)*sin(x2)*u*x1+(p+E)*x1*cos(x2)*v)/(x1*sin(x2)*x1);
      }
      else if(i == 2)
      {
         a[i] = ((n*sin(x2)*pow(w,2.0)+n*sin(x2)*pow(v,2.0)-2*n*sin(x2)*pow(u,2.0))*x1-n*x1*cos(x2)*u*v)/(x1*sin(x2)*x1);
      }
      else if(i == 3)
      {
         a[i] = ((n*cos(x2)*pow(w,2.0)-3*n*sin(x2)*u*v)*x1-n*x1*cos(x2)*pow(v,2.0))/(x1*sin(x2)*x1);
      }
      else if(i == 4)
      {
         a[i] = -((n*cos(x2)*v+3*n*sin(x2)*u)*w*x1+n*x1*cos(x2)*v*w)/(x1*sin(x2)*x1);
      }
   }

   return 0;
}
