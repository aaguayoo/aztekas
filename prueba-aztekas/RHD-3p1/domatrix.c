#include<stdio.h>
#include<math.h>
#include"../Headers/matrix.h"
#include"../Headers/main.h"
    
double funct_Do(double *a, double *uu)
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
   W = 1/sqrt(-pow(x1,2.0)*pow(w,2.0)-pow(v,2.0)-pow(u,2.0)+1);
   O = (K*p+(K-1)*n)/(K-1);
    
   dWu = -u/(sqrt(-pow(x1,2.0)*pow(w,2.0)-pow(v,2.0)-pow(u,2.0)+1)*(pow(x1,2.0)*pow(w,2.0)+pow(v,2.0)+pow(u,2.0)-1));
   dWv = -v/(sqrt(-pow(x1,2.0)*pow(w,2.0)-pow(v,2.0)-pow(u,2.0)+1)*(pow(x1,2.0)*pow(w,2.0)+pow(v,2.0)+pow(u,2.0)-1));
   dWw = -pow(x1,2.0)*w/(sqrt(-pow(x1,2.0)*pow(w,2.0)-pow(v,2.0)-pow(u,2.0)+1)*(pow(x1,2.0)*pow(w,2.0)+pow(v,2.0)+pow(u,2.0)-1));
   dOn = 1;
   dOp = K/(K-1);
    
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
    
   for(i = 0; i <= 2; i++)
   {
      if(i == 0)
      {
         a[i] = ((x1*w*O-K*p*x1*w)*fabs(O)+O*sqrt((K*p*pow(x1,4.0)*pow(w,4.0)+(K*p*pow(x1,2.0)*pow(v,2.0)+K*p*pow(x1,2.0)*pow(u,2.0)-2*K*p*pow(x1,2.0))*pow(w,2.0)-K*p*pow(v,2.0)-K*p*pow(u,2.0)+K*p)*O+(pow(K,2.0)*pow(p,2.0)*pow(x1,2.0)*pow(v,2.0)+pow(K,2.0)*pow(p,2.0)*pow(x1,2.0)*pow(u,2.0))*pow(w,2.0)+pow(K,2.0)*pow(p,2.0)*pow(v,4.0)+(2*pow(K,2.0)*pow(p,2.0)*pow(u,2.0)-pow(K,2.0)*pow(p,2.0))*pow(v,2.0)+pow(K,2.0)*pow(p,2.0)*pow(u,4.0)-pow(K,2.0)*pow(p,2.0)*pow(u,2.0)))/((x1*O-K*p*pow(x1,3.0)*pow(w,2.0)-K*p*x1*pow(v,2.0)-K*p*x1*pow(u,2.0))*fabs(O));
      }
      else if(i == 1)
      {
         a[i] = -((K*p*x1*w-x1*w*O)*fabs(O)+O*sqrt((K*p*pow(x1,4.0)*pow(w,4.0)+(K*p*pow(x1,2.0)*pow(v,2.0)+K*p*pow(x1,2.0)*pow(u,2.0)-2*K*p*pow(x1,2.0))*pow(w,2.0)-K*p*pow(v,2.0)-K*p*pow(u,2.0)+K*p)*O+(pow(K,2.0)*pow(p,2.0)*pow(x1,2.0)*pow(v,2.0)+pow(K,2.0)*pow(p,2.0)*pow(x1,2.0)*pow(u,2.0))*pow(w,2.0)+pow(K,2.0)*pow(p,2.0)*pow(v,4.0)+(2*pow(K,2.0)*pow(p,2.0)*pow(u,2.0)-pow(K,2.0)*pow(p,2.0))*pow(v,2.0)+pow(K,2.0)*pow(p,2.0)*pow(u,4.0)-pow(K,2.0)*pow(p,2.0)*pow(u,2.0)))/((x1*O-K*p*pow(x1,3.0)*pow(w,2.0)-K*p*x1*pow(v,2.0)-K*p*x1*pow(u,2.0))*fabs(O));
      }
      else if(i == 2)
      {
         a[i] = w;
      }
   }

   return 0;
}
