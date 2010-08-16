/*
   File Name: rand.c
   System: hmmSNP
   Programmer: Glen Beane and Jin Szatkeiwicz
   Comments: 
   
   these functions were provided by Jin Szatkeiwicz from her dissertation, but
   they are pretty standard algorithms
   
   
 */

/* RANDOM NUBMER GENERATORS FOR RANDOM UNIFORM AND RANDOM EXPONENTIAL*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<limits.h>
#include<time.h>
#include<assert.h>

#define MIN(a,b) ( ((a)<(b)) ? (a):(b) )
#define MAX(a,b) ( ((a)>(b)) ? (a):(b) )


#include "rand.h"

/* Generate a standard uniform */
/* UINT_MAX should be defined in /usr/include/limits.h */

double uniformRandom(void) /* return a single Uniform(0,1) */
{

  long j;
  static unsigned int seed = 0;
  
  if (seed == 0)
  {
      srand(time(NULL));
      seed = (unsigned) rand();
      for (j = 0; j < 20000000; j++)
      {
              seed = 109441 * seed * seed + 850373 * seed + 1299827;
      } 
  }
 
  for (j = 0; j < 3; j++)
  {
     seed = 850373 * seed * seed + 1000003 * ((unsigned) rand()) * seed + 1299827;
  }

  return (double) seed / UINT_MAX;

}


/* return a deviate distributed as a gamma distribution of interger order iA */
/* this is a standard cookbook gamma deviate function */
double gammaDev(int shape)
{

   double am,e,s,v1,v2,y;
   
   double gamma;

   assert(shape >= 1);
   
   if (shape < 6)
   {
      gamma = 1.0;
      
      for (int i = 0; i <= shape; i++)
         gamma *= uniformRandom();
         
      gamma = -(log(gamma));  
   }
   else
   {   
      do
      {
         do
         {
            do
            {
               v1 = uniformRandom();
               v2 = 2.0 * uniformRandom() - 1.0;
            } while (v1 * v1 + v2 * v2 > 1.0);
            
            y = v2 / v1;
            am = shape - 1;
            s = sqrt(2.0 * am + 1.0);
            gamma = s * y + am;
            
         } while (gamma <= 0.0);
         
         e = (1.0 + y * y) * exp(am * log(gamma / am) - s * y);
         
      } while (uniformRandom() > e);
   }
      
   return gamma;
}

double beta(int a, int b)
{
   return gammaDev(a) / ( gammaDev(a) + gammaDev(b) );
}


/* Generate a single exponential variable with mean 1/lambda. */

double rexp(double lambda)

{ 
  return(-log(uniformRandom())/lambda);
}
