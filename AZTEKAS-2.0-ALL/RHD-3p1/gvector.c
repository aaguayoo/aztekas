#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
double funct_G(double *a, double *uu)
{
   int i;
   double n, p, u=0, v=0, w=0;
   double r, W, O;
   double dWu, dWv, dWw;
   double dOn, dOp;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
    
   r = x1;
   W = 1/sqrt((-pow(x1,2.0)*pow(w,2.0))-pow(v,2.0)-pow(u,2.0)+1);
   O = (K*p+(K-1)*n)/(K-1);
    
   if(u >= 1)
   {
      u = 0.9;
   }
   else if(v >= 1)
   {
      v = 0.9;
   }
   else if(w >= 1)
   {
      w = 0.9;
   }
    
   for(i = 0; i <= eq; i++)
   {
      if(i == 0)
      {
         a[i] = W*n*v;
      }
      else if(i == 1)
      {
         a[i] = O*pow(W,2.0)*v;
      }
      else if(i == 2)
      {
         a[i] = O*pow(W,2.0)*u*v;
      }
      else if(i == 3)
      {
         a[i] = O*pow(W,2.0)*pow(v,2.0)+p;
      }
      else if(i == 4)
      {
         a[i] = O*pow(W,2.0)*pow(x1,2.0)*v*w;
      }
   }
   
   return 0;
}
