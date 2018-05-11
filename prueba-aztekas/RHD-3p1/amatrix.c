#include<stdio.h>
#include<math.h>
#include"../Headers/matrix.h"
#include"../Headers/main.h"
    
double funct_A(double *a, double *uu)
{
   int i, j;
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
    
   for(i = 0; i < eq; i++)
   {
      for(j = 0; j < eq; j++)
      {
         if(i == 0 && j == 0)
         {
            a[i*eq + j] = -(2*w*O*dWw+2*v*O*dWv+2*u*O*dWu-O*dOp*pow(W,3.0)+O*W)/((n*w*dOn-2*w*O)*W*dWw+(n*v*dOn-2*v*O)*W*dWv+(n*u*dOn-2*u*O)*W*dWu+O*dOp*pow(W,4.0)-O*pow(W,2.0));
         }
         else if(i == 0 && j == 1)
         {
            a[i*eq + j] = (n*w*dOp*dWw+n*v*dOp*dWv+n*u*dOp*dWu)/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 0 && j == 2)
         {
            a[i*eq + j] = -(n*dOp*pow(W,2.0)-n)*dWu/((n*w*dOn-2*w*O)*pow(W,2.0)*dWw+(n*v*dOn-2*v*O)*pow(W,2.0)*dWv+(n*u*dOn-2*u*O)*pow(W,2.0)*dWu+O*dOp*pow(W,5.0)-O*pow(W,3.0));
         }
         else if(i == 0 && j == 3)
         {
            a[i*eq + j] = -(n*dOp*pow(W,2.0)-n)*dWv/((n*w*dOn-2*w*O)*pow(W,2.0)*dWw+(n*v*dOn-2*v*O)*pow(W,2.0)*dWv+(n*u*dOn-2*u*O)*pow(W,2.0)*dWu+O*dOp*pow(W,5.0)-O*pow(W,3.0));
         }
         else if(i == 0 && j == 4)
         {
            a[i*eq + j] = -(n*dOp*pow(W,2.0)-n)*dWw/((n*pow(x1,2.0)*w*dOn-2*pow(x1,2.0)*w*O)*pow(W,2.0)*dWw+(n*pow(x1,2.0)*v*dOn-2*pow(x1,2.0)*v*O)*pow(W,2.0)*dWv+(n*pow(x1,2.0)*u*dOn-2*pow(x1,2.0)*u*O)*pow(W,2.0)*dWu+pow(x1,2.0)*O*dOp*pow(W,5.0)-pow(x1,2.0)*O*pow(W,3.0));
         }
         else if(i == 1 && j == 0)
         {
            a[i*eq + j] = -O*dOn*pow(W,2.0)/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 1 && j == 1)
         {
            a[i*eq + j] = -((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu-O*W)/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 1 && j == 2)
         {
            a[i*eq + j] = (n*dOn-2*O)*dWu/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 1 && j == 3)
         {
            a[i*eq + j] = (n*dOn-2*O)*dWv/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 1 && j == 4)
         {
            a[i*eq + j] = (n*dOn-2*O)*dWw/((n*pow(x1,2.0)*w*dOn-2*pow(x1,2.0)*w*O)*dWw+(n*pow(x1,2.0)*v*dOn-2*pow(x1,2.0)*v*O)*dWv+(n*pow(x1,2.0)*u*dOn-2*pow(x1,2.0)*u*O)*dWu+pow(x1,2.0)*O*dOp*pow(W,3.0)-pow(x1,2.0)*O*W);
         }
         else if(i == 2 && j == 0)
         {
            a[i*eq + j] = u*dOn/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 2 && j == 1)
         {
            a[i*eq + j] = -u*dOp*W/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 2 && j == 2)
         {
            a[i*eq + j] = ((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+O*dOp*pow(W,3.0)-O*W)/((n*w*O*dOn-2*w*pow(O,2.0))*pow(W,2.0)*dWw+(n*v*O*dOn-2*v*pow(O,2.0))*pow(W,2.0)*dWv+(n*u*O*dOn-2*u*pow(O,2.0))*pow(W,2.0)*dWu+pow(O,2.0)*dOp*pow(W,5.0)-pow(O,2.0)*pow(W,3.0));
         }
         else if(i == 2 && j == 3)
         {
            a[i*eq + j] = -(n*u*dOn-2*u*O)*dWv/((n*w*O*dOn-2*w*pow(O,2.0))*pow(W,2.0)*dWw+(n*v*O*dOn-2*v*pow(O,2.0))*pow(W,2.0)*dWv+(n*u*O*dOn-2*u*pow(O,2.0))*pow(W,2.0)*dWu+pow(O,2.0)*dOp*pow(W,5.0)-pow(O,2.0)*pow(W,3.0));
         }
         else if(i == 2 && j == 4)
         {
            a[i*eq + j] = -(n*u*dOn-2*u*O)*dWw/((n*pow(x1,2.0)*w*O*dOn-2*pow(x1,2.0)*w*pow(O,2.0))*pow(W,2.0)*dWw+(n*pow(x1,2.0)*v*O*dOn-2*pow(x1,2.0)*v*pow(O,2.0))*pow(W,2.0)*dWv+(n*pow(x1,2.0)*u*O*dOn-2*pow(x1,2.0)*u*pow(O,2.0))*pow(W,2.0)*dWu+pow(x1,2.0)*pow(O,2.0)*dOp*pow(W,5.0)-pow(x1,2.0)*pow(O,2.0)*pow(W,3.0));
         }
         else if(i == 3 && j == 0)
         {
            a[i*eq + j] = v*dOn/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 3 && j == 1)
         {
            a[i*eq + j] = -v*dOp*W/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 3 && j == 2)
         {
            a[i*eq + j] = -(n*v*dOn-2*v*O)*dWu/((n*w*O*dOn-2*w*pow(O,2.0))*pow(W,2.0)*dWw+(n*v*O*dOn-2*v*pow(O,2.0))*pow(W,2.0)*dWv+(n*u*O*dOn-2*u*pow(O,2.0))*pow(W,2.0)*dWu+pow(O,2.0)*dOp*pow(W,5.0)-pow(O,2.0)*pow(W,3.0));
         }
         else if(i == 3 && j == 3)
         {
            a[i*eq + j] = ((n*w*dOn-2*w*O)*dWw+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W)/((n*w*O*dOn-2*w*pow(O,2.0))*pow(W,2.0)*dWw+(n*v*O*dOn-2*v*pow(O,2.0))*pow(W,2.0)*dWv+(n*u*O*dOn-2*u*pow(O,2.0))*pow(W,2.0)*dWu+pow(O,2.0)*dOp*pow(W,5.0)-pow(O,2.0)*pow(W,3.0));
         }
         else if(i == 3 && j == 4)
         {
            a[i*eq + j] = -(n*v*dOn-2*v*O)*dWw/((n*pow(x1,2.0)*w*O*dOn-2*pow(x1,2.0)*w*pow(O,2.0))*pow(W,2.0)*dWw+(n*pow(x1,2.0)*v*O*dOn-2*pow(x1,2.0)*v*pow(O,2.0))*pow(W,2.0)*dWv+(n*pow(x1,2.0)*u*O*dOn-2*pow(x1,2.0)*u*pow(O,2.0))*pow(W,2.0)*dWu+pow(x1,2.0)*pow(O,2.0)*dOp*pow(W,5.0)-pow(x1,2.0)*pow(O,2.0)*pow(W,3.0));
         }
         else if(i == 4 && j == 0)
         {
            a[i*eq + j] = w*dOn/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 4 && j == 1)
         {
            a[i*eq + j] = -w*dOp*W/((n*w*dOn-2*w*O)*dWw+(n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W);
         }
         else if(i == 4 && j == 2)
         {
            a[i*eq + j] = -(n*w*dOn-2*w*O)*dWu/((n*w*O*dOn-2*w*pow(O,2.0))*pow(W,2.0)*dWw+(n*v*O*dOn-2*v*pow(O,2.0))*pow(W,2.0)*dWv+(n*u*O*dOn-2*u*pow(O,2.0))*pow(W,2.0)*dWu+pow(O,2.0)*dOp*pow(W,5.0)-pow(O,2.0)*pow(W,3.0));
         }
         else if(i == 4 && j == 3)
         {
            a[i*eq + j] = -(n*w*dOn-2*w*O)*dWv/((n*w*O*dOn-2*w*pow(O,2.0))*pow(W,2.0)*dWw+(n*v*O*dOn-2*v*pow(O,2.0))*pow(W,2.0)*dWv+(n*u*O*dOn-2*u*pow(O,2.0))*pow(W,2.0)*dWu+pow(O,2.0)*dOp*pow(W,5.0)-pow(O,2.0)*pow(W,3.0));
         }
         else if(i == 4 && j == 4)
         {
            a[i*eq + j] = ((n*v*dOn-2*v*O)*dWv+(n*u*dOn-2*u*O)*dWu+O*dOp*pow(W,3.0)-O*W)/((n*pow(x1,2.0)*w*O*dOn-2*pow(x1,2.0)*w*pow(O,2.0))*pow(W,2.0)*dWw+(n*pow(x1,2.0)*v*O*dOn-2*pow(x1,2.0)*v*pow(O,2.0))*pow(W,2.0)*dWv+(n*pow(x1,2.0)*u*O*dOn-2*pow(x1,2.0)*u*pow(O,2.0))*pow(W,2.0)*dWu+pow(x1,2.0)*pow(O,2.0)*dOp*pow(W,5.0)-pow(x1,2.0)*pow(O,2.0)*pow(W,3.0));
         }
     
      }
     
   }
   return 0;
}
