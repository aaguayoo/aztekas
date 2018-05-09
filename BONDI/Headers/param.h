//DOMAIN//
#define x1max 11.0
#define x2max 11.0
#define x3max 1.0
#define x1min 0.0
#define x2min 0.0
#define x3min 0.0

//MESH//
#define Nx1 106
#define Nx2 106
#define Nx3 0

//MAX TIME & PRINT TIME//
#define tmax      51.
#define timefile  10.

//DIMENSION PARAMETERS (num of equations, dim of graf, physical dimension)//
#define eq    4
#define graf  2
#define dim   2
#define alfa  1

//COURANT NUMBER//
#define cou   0.5

//POLYTROPIC INDEX//
#define K  1.333333333333

//ACCRETION DOMAIN//
#define r_int 0.10
#define r_bou 10.0

//DENSITY PROFILE
// (theta_0 should be in units of pi, delta_theta is the width of the step function)
#define theta_0 0.6
#define delta_theta 0.02

//OUTPUT FILE AND DIRECTORY
#define outputdirectory "pi_6/"
#define outputfile  "back_"

