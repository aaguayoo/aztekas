#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int funct_A(double *a, double *uu)
{
   int i, j;
   double n, p, u=0, v=0, w=0;
   double R, W, h;
   double dWu, dWv, dWw;
   double dhn, dhp;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}

   R = x1;
   W = (x1*fabs(sin(x2)))/sqrt(-((x1+2*MM)*pow(w,2.0)+(x1+2*MM)*pow(sin(x2),2.0)*pow(v,2.0)+pow(x1,3.0)*pow(sin(x2),2.0)*pow(u,2.0)+((-pow(x1,3.0))-2*MM*pow(x1,2.0))*pow(sin(x2),2.0))/(x1+2*MM));
   h = 1;

   dWu = (pow(W,3.0)*x1*u)/(x1+2*MM);
   dWv = (pow(W,3.0)*v)/pow(x1,2.0);
   dWw = (pow(W,3.0)*w)/(pow(x1,2.0)*pow(sin(x2),2.0));
   dhn = 0;
   dhp = 0;

   a[0*eq + 0] = (dWw*w+dWv*v+dWu*u+pow(W,2.0)*h)/(pow(W,3.0)*h);
   a[0*eq + 1] = (dWw*w+dWv*v+dWu*u)/(pow(W,3.0)*h);
   a[0*eq + 2] = -dWu/(pow(W,3.0)*h);
   a[0*eq + 3] = -dWv/(pow(W,3.0)*h);
   a[0*eq + 4] = -dWw/(pow(W,3.0)*h);
   a[1*eq + 0] = -((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,2.0)*h*dhn*n+pow(W,2.0)*pow(h,2.0)-W*h)/(pow(W,3.0)*h*dhp*n);
   a[1*eq + 1] = -((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u-W*h)/(pow(W,3.0)*h*dhp*n);
   a[1*eq + 2] = (dWu*dhn*n-dWu*h)/(pow(W,3.0)*h*dhp*n);
   a[1*eq + 3] = (dWv*dhn*n-dWv*h)/(pow(W,3.0)*h*dhp*n);
   a[1*eq + 4] = (dWw*dhn*n-dWw*h)/(pow(W,3.0)*h*dhp*n);
   a[2*eq + 0] = -u/(pow(W,2.0)*h*n);
   a[2*eq + 1] = -u/(pow(W,2.0)*h*n);
   a[2*eq + 2] = 1/(pow(W,2.0)*h*n);
   a[2*eq + 3] = 0;
   a[2*eq + 4] = 0;
   a[3*eq + 0] = -v/(pow(W,2.0)*h*n);
   a[3*eq + 1] = -v/(pow(W,2.0)*h*n);
   a[3*eq + 2] = 0;
   a[3*eq + 3] = 1/(pow(W,2.0)*h*n);
   a[3*eq + 4] = 0;
   a[4*eq + 0] = -w/(pow(W,2.0)*h*n);
   a[4*eq + 1] = -w/(pow(W,2.0)*h*n);
   a[4*eq + 2] = 0;
   a[4*eq + 3] = 0;
   a[4*eq + 4] = 1/(pow(W,2.0)*h*n);

   return 0;
}
