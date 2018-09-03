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
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "./Headers/main.h"

double gtheta(double th);

void OUTFLOW(double *B)
{
  int i, j, k, n, cell;

  if (dim == 1)
  {
    for (n = 0; n < eq; n++)
    {
      for (cell = 0; cell < gc; cell++)
      {
        if (alfa > 0)
        {
          B[c1(n, gc)] = B[c1(n, gc + 1)];
        }

        B[c1(n, cell)] = B[c1(n, gc)];
        B[c1(n, Nx1 - cell)] = B[c1(n, Nx1 - gc)];
      }
    }
  }
  else if (dim == 2)
  {
    for (n = 0; n < eq; n++)
    {
      for (i = 0; i <= Nx1 - 0; i++)
      {
        if (alfa == 2)
        {
          B[c2(n, i, gc)] = B[c2(n, i, gc + 1)];
          B[c2(n, i, Nx2 - gc)] = B[c2(n, i, Nx2 - gc - 1)];
        }
        for (cell = 0; cell < gc; cell++)
        {
          B[c2(n, i, cell)] = B[c2(n, i, gc)];
          B[c2(n, i, Nx2 - cell)] = B[c2(n, i, Nx2 - gc)];
        }
      }

      for (j = 0; j <= Nx2; j++)
      {
        for (cell = 0; cell < gc; cell++)
        {
          if (alfa > 0)
          {
            B[c2(n, gc, j)] = B[c2(n, gc + 1, j)];
          }

          B[c2(n, cell, j)] = B[c2(n, gc, j)];
          B[c2(n, Nx1 - cell, j)] = B[c2(n, Nx1 - gc, j)];
        }
      }
    }
  }
}

void PERIODIC(double *B, int r, int l, int u, int d, int f, int b)
{
  int i, j, k, n, cell;
  if (dim == 1)
  {
    for (n = 0; n < eq; n++)
    {
      for (i = 0; i <= Nx1 - 0; i++)
      {
        for (cell = 0; cell < gc; cell++)
        {
          if (l == 1)
          {
            B[c1(n, cell)] = B[c1(n, Nx1 - 2 * gc + cell + 1)];
          }
          if (r == 1)
          {
            B[c1(n, Nx1 - cell)] = B[c1(n, 2 * gc - cell - 1)];
          }
        }
      }
    }
  }
  else if (dim == 2)
  {
    for (n = 0; n < eq; n++)
    {
      for (i = 0; i <= Nx1 - 0; i++)
      {
        for (cell = 0; cell < gc; cell++)
        {
          if (d == 1)
          {
            B[c2(n, i, cell)] = B[c2(n, i, Nx2 - 2 * gc + cell + 1)];
          }
          if (u == 1)
          {
            B[c2(n, i, Nx2 - cell)] = B[c2(n, i, 2 * gc - cell - 1)];
          }
        }
      }

      for (j = 0; j <= Nx2; j++)
      {
        for (cell = 0; cell < gc; cell++)
        {
          if (l == 1)
          {
            B[c2(n, cell, j)] = B[c2(n, Nx1 - 2 * gc + cell + 1, j)];
          }
          if (r == 1)
          {
            B[c2(n, Nx1 - cell, j)] = B[c2(n, 2 * gc - cell - 1, j)];
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////
// REFLECTION:
//    This function stablish a reflective condition in the x1, -x1,
//    x2, -x2, x3, -x3 boundary, respectibly.
//
//    The conditions fix the point gc and Nxi-gc, i.e., the first
//    and last points of the integration domain, to fulfil the
//    condition: Tangent velocity = 0.0 and V[gc-1] = - V[gc+1]
////////////////////////////////////////////////////////////////////

void REFLECTION(double *B, int r, int l, int u, int d, int f, int b)
{
  int i, j, k, n, cell;
  if (dim == 1)
  {
    for (i = 0; i <= Nx1 - 0; i++)
    {
      for (cell = 0; cell < gc; cell++)
      {
        if (l == 1)
        {
          for (n = 0; n < eq; n++)
          {
            B[c1(n, cell)] = B[c1(n, 2 * gc - cell - 1)];
          }

          B[c1(2, cell)] = -B[c1(2, 2 * gc - cell - 1)];
        }
        if (r == 1)
        {
          for (n = 0; n < eq; n++)
          {
            B[c1(n, Nx1 - cell)] = B[c1(n, Nx1 - 2 * gc + cell + 1)];
          }

          B[c1(2, Nx1 - cell)] = -B[c1(2, Nx1 - 2 * gc + cell + 1)];
        }
      }
    }
  }
  else if (dim == 2)
  {
    for (i = 0; i <= Nx1 - 0; i++)
    {
      for (cell = 0; cell < gc; cell++)
      {
        if (d == 1)
        {
          for (n = 0; n < eq; n++)
          {
            B[c2(n, i, cell)] = B[c2(n, i, 2 * gc - cell)];
          }

          B[c2(3, i, cell)] = -B[c2(3, i, 2 * gc - cell)];
          B[c2(3, i, gc)] = 0.0;
        }

        if (u == 1)
        {
          for (n = 0; n < eq; n++)
          {
            B[c2(n, i, Nx2 - cell)] = B[c2(n, i, Nx2 - 2 * gc + cell)];
          }

          B[c2(3, i, Nx2 - cell)] = -B[c2(3, i, Nx2 - 2 * gc + cell)];
          B[c2(3, i, Nx2 - gc)] = 0.0;
        }
      }
    }

    for (j = 0; j <= Nx2; j++)
    {
      for (cell = 0; cell < gc; cell++)
      {
        if (l == 1)
        {
          for (n = 0; n < eq; n++)
          {
            B[c2(n, cell, j)] = B[c2(n, 2 * gc - cell, j)];
          }

          B[c2(2, cell, j)] = -B[c2(2, 2 * gc - cell, j)];
          B[c2(2, gc, j)] = 0.0;
        }

        if (r == 1)
        {
          for (n = 0; n < eq; n++)
          {
            B[c2(n, Nx1 - cell, j)] = B[c2(n, Nx1 - 2 * gc + cell, j)];
          }

          B[c2(2, Nx1 - cell, j)] = -B[c2(2, Nx1 - 2 * gc + cell, j)];
          B[c2(2, Nx1 - gc, j)] = 0.0;
        }
      }
    }
  }
}

void JET_LAUNCH(double *B)
{
  int i, j, k, n;
  if (dim == 1)
  {
    for (i = 0; i <= Nx1; i++)
    {
      if (fabs(X1[i]) <= r_jet)
      {
        B[c1(0, i)] = n_jet;
        B[c1(1, i)] = p_jet;
        B[c1(2, i)] = vx1_jet;
      }
    }
  }
  else if (dim == 2)
  {
    for (i = 0; i <= Nx1; i++)
    {
      for (j = 0; j <= Nx2; j++)
      {
        if (fabs(X1[i]) <= 1.0)
        {
          if (fabs(X2[j]) <= 1.0)
          {
            B[c2(0, i, j)] = n_jet;
            B[c2(1, i, j)] = p_jet;
            B[c2(2, i, j)] = vx1_jet;
            B[c2(3, i, j)] = vx2_jet;
          }
        }
      }
    }
  }
}

void IN_OUT_BOUND(double *B)
{
  int i, j, k, n;
  double x1, x2, x3, r;
  double dummy;

  if (dim == 1)
  {
    for (i = 0; i <= Nx1; i++)
    {
      r = X1[i];
      r_in = X1[gc + 2];

      if (r > r_out)
      {
        B[c1(0, i)] = density_0;
        B[c1(1, i)] = pressure_0;
        B[c1(2, i)] = velocity_0;
      }
      else if (r <= r_in)
      {
        B[c1(0, i)] = density_0;
        B[c1(1, i)] = pressure_0;
        B[c1(2, i)] = 0.0;
      }
    }
  }
  else if (dim == 2)
  {
    int i1, i2, i3;
    double r, r1, r2, r3, rho, rho1, rho2, rho3;
    double vr1, vr2, vr3, vth1, vth2, vth3, th;
    double r1r2, r1r3, r2r3, fac1, fac2, fac3;

    // set analytic values at injection boundary
    i = Nx1 - gc;
    for (j = 0; j <= Nx2; j++)
    {
      th = X2[j];
      rho = gtheta(th) * density_0;

      B[c2(0, i, j)] = rho;
      B[c2(1, i, j)] = pow(rho, K) / K;
      // B[c2(2, i, j)] = velocity_0;
      // B[c2(3, i, j)] = 0.0;
    }

    // Inner boundary
    switch (extrapolation)
    {
    case 0:
      // copy constant values
      for (i = gc; i >= 0; i--)
      {
        for (j = 0; j <= Nx2; j++)
        {
          B[c2(0, i, j)] = B[c2(0, gc + 1, j)];
          B[c2(1, i, j)] = B[c2(1, gc + 1, j)];
          B[c2(2, i, j)] = B[c2(2, gc + 1, j)];
          B[c2(3, i, j)] = B[c2(3, gc + 1, j)];
        }
      }
      break;

    case 1:
      // first order extrapolation to ghost cells
      i1 = gc + 1; 
      i2 = gc + 2;
      r1 = X1[i1];
      r2 = X1[i2];
      r1r2 = r1 - r2;
      for (i = gc; i >= 0; i--)
      {
        r = X1[i];
        for (j = 0; j <= Nx2; j++)
        {

          rho1 = B[c2(0, i1, j)];
          rho2 = B[c2(0, i2, j)];

          vr1 = B[c2(2, i1, j)];
          vr2 = B[c2(2, i2, j)];

          vth1 = B[c2(3, i1, j)];
          vth2 = B[c2(3, i2, j)];

          fac1 = (r - r2)/r1r2;
          fac2 = -(r - r1)/r1r2;

          rho = fac1 * rho1 + fac2 * rho2 ;

          B[c2(0, i, j)] = rho;
          B[c2(1, i, j)] = pow(rho, K) / K;
          B[c2(2, i, j)] = fac1 * vr1 + fac2 * vr2 ;
          B[c2(3, i, j)] = fac1 * vth1 + fac2 * vth2 ;
        }
      }
      break;

    case 2:
      // second order extrapolation to ghost cells
      r1 = X1[gc + 1];
      r2 = X1[gc + 2];
      r3 = X1[gc + 3];
      r1r2 = r1 - r2;
      r1r3 = r1 - r3;
      r2r3 = r2 - r3;
      for (i = gc; i >= 0; i--)
      {
        r = X1[i];

        for (j = 0; j <= Nx2; j++)
        {
          rho1 = B[c2(0, gc + 1, j)];
          rho2 = B[c2(0, gc + 2, j)];
          rho3 = B[c2(0, gc + 3, j)];

          vr1 = B[c2(2, gc + 1, j)];
          vr2 = B[c2(2, gc + 2, j)];
          vr3 = B[c2(2, gc + 3, j)];

          vth1 = B[c2(3, gc + 1, j)];
          vth2 = B[c2(3, gc + 2, j)];
          vth3 = B[c2(3, gc + 3, j)];

          fac1 = (r - r2) * (r - r3) / (r1r2 * r1r3);
          fac2 = -(r - r1) * (r - r3) / (r1r2 * r2r3);
          fac3 = (r - r1) * (r - r2) / (r2r3 * r1r3);

          rho = fac1 * rho1 + fac2 * rho2 + fac3 * rho3;

          B[c2(0, i, j)] = rho;
          B[c2(1, i, j)] = pow(rho, K) / K;
          B[c2(2, i, j)] = fac1 * vr1 + fac2 * vr2 + fac3 * vr3;
          B[c2(3, i, j)] = fac1 * vth1 + fac2 * vth2 + fac3 * vth3;
        }
      }
      break;
    }
    // Outer boundary
    switch (extrapolation)
    {
    case 0:
      // copy constant values
      for (i = Nx1 - gc + 1; i <= Nx1; i++)
      {
        for (j = 0; j <= Nx2; j++)
        {
          B[c2(0, i, j)] = B[c2(0, Nx1 - gc, j)];
          B[c2(1, i, j)] = B[c2(1, Nx1 - gc, j)];
          B[c2(2, i, j)] = B[c2(2, Nx1 - gc, j)];
          B[c2(3, i, j)] = B[c2(3, Nx1 - gc, j)];
        }
      }
      break;

    case 1:
      // first order extrapolation to ghost cells
      i1 = Nx1 - gc ;
      i2 = Nx1 - gc - 1;
      r1 = X1[i1];
      r2 = X1[i2];
      r1r2 = r1 - r2;
      for (i = Nx1 - gc + 1; i <= Nx1; i++)      
      {
        r = X1[i];
        for (j = 0; j <= Nx2; j++)
        {
          rho1 = B[c2(0, i1, j)];
          rho2 = B[c2(0, i2, j)];

          vr1 = B[c2(2, i1, j)];
          vr2 = B[c2(2, i2, j)];

          vth1 = B[c2(3, i1, j)];
          vth2 = B[c2(3, i2, j)];

          fac1 = (r - r2)/r1r2;
          fac2 = -(r - r1)/r1r2;

          rho = fac1 * rho1 + fac2 * rho2 ;

          B[c2(0, i, j)] = rho;
          B[c2(1, i, j)] = pow(rho, K) / K;
          B[c2(2, i, j)] = fac1 * vr1 + fac2 * vr2 ;
          B[c2(3, i, j)] = fac1 * vth1 + fac2 * vth2 ;
        }
      }
      break;

    case 2:
      // second order extrapolation to ghost cells above
      i1 = Nx1 - gc ;
      i2 = Nx1 - gc - 1;
      i3 = Nx1 - gc - 2;
      r1 = X1[i1];
      r2 = X1[i2];
      r3 = X1[i3];
      r1r2 = r1 - r2;
      r1r3 = r1 - r3;
      r2r3 = r2 - r3;
      for (i = Nx1 - gc + 1; i <= Nx1; i++)
      {
        r = X1[i];
        for (j = 0; j <= Nx2; j++)
        {
          rho1 = B[c2(0, i1, j)];
          rho2 = B[c2(0, i2, j)];
          rho3 = B[c2(0, i3, j)];

          vr1 = B[c2(2, i1, j)];
          vr2 = B[c2(2, i2, j)];
          vr3 = B[c2(2, i3, j)];

          vth1 = B[c2(3, i1, j)];
          vth2 = B[c2(3, i2, j)];
          vth3 = B[c2(3, i3, j)];

          fac1 = (r - r2) * (r - r3) / (r1r2 * r1r3);
          fac2 = -(r - r1) * (r - r3) / (r1r2 * r2r3);
          fac3 = (r - r1) * (r - r2) / (r2r3 * r1r3);

          rho = fac1 * rho1 + fac2 * rho2 + fac3 * rho3;

          B[c2(0, i, j)] = rho;
          B[c2(1, i, j)] = pow(rho, K) / K;
          B[c2(2, i, j)] = fac1 * vr1 + fac2 * vr2 + fac3 * vr3;
          B[c2(3, i, j)] = fac1 * vth1 + fac2 * vth2 + fac3 * vth3;
        }
      }
      break;
    }
  }
}

double gtheta(double th)
{
  if (rho_boundary == 1)
  // constant profile
  {
    return 1.0;
  }

  if (rho_boundary == 2)
  // gaussian profile for the density at the injection boundary
  {
    return 1.0; //exp(-pow(th / (sqrt(2.0) * th_0), 2.0));
  }
  else if (rho_boundary == 3)
  // sin**2 profile for the density at the injection boundary
  {
    return (1.0 - rho_atm) * pow(sin(th), n_pow) + rho_atm;
  }

  return 1.0;
}