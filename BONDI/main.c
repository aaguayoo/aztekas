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
#include<stdio.h>
#include<omp.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"./Headers/main.h"

int main(int argc, char* argv[])
{
   int itprint;
   double dtprint, tprint;

   if(argc != 2) 
   {
      printf("%s\n","Wrong number of arguments") ;
      printf("%s\n","Execute as:") ;    
      printf("%s\n","./aztekas paramfile") ;  
      exit(EXIT_FAILURE);
   }

   strcpy(paramfile_name, argv[1]);
   
   read_parameters_file(paramfile_name);   
      
	 // create output directory
   char create_dir[] = "mkdir -p ";	
	 strcat(create_dir,outputdirectory);	
	 int sysret = system(create_dir);

   allocateArray();

   //We set the mesh func_planarMESH.c
   MESH();
   
   
   // time interval between data dumps
   dtprint = timefile;

   //We set the initial parameters func_planarINITIAL.c
   if ( restart_simulation ) 
   {
      RESTART();
      
      //Initialize printing parameter   
      tprint  = time; 
      //Initialize file numeration
      itprint = restart_filecount;
   }
   else
   {
      INITIAL();
      
      //Initialize printing parameter   
      tprint  = 0.0; 
      //Initialize file numeration
      itprint = 0;   
   }
   
   printf("dump : %i\n", itprint);

   printf("time : %f\n",time) ;   
   
   start = omp_get_wtime();
   while(time <= tmax)
   {
      //In this part we compute the time step
      dt = TIMESTEP();

      //We print the values: file (DATOS*) and to terminal func_planarOUTPUT.c
      PrintValues(&tprint,&dtprint,&itprint);

      //In here we set the integration method (Finite volume method)
      INTEGRATION();
   }

   delta = omp_get_wtime() - start;
   printf("Delta %.4g seconds with %d threads\n",delta,4);

   free(U);
   free(U1);
   free(X1);
   free(X2);
   free(X3);

   return 0;
}
