/*
   File Name: version.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
   Date developed: 3/07
   Purpose: contains function to print version information
   Comments:
 
   
  Copyright (c) 2010 The Jackson Laboratory
 
  This software was developed by Gary Churchill's Lab at The Jackson
  Laboratory (see http://research.jax.org/faculty/churchill).
 
  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This software is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this software. If not, see <http://www.gnu.org/licenses/>.
 */
 
 
#include <stdio.h> 
 
#include "constants.h" 
 
 
/*******************************************************************************
 * void version()
 *
 * INPUTS: void
 *
 * RETURNS: void
 *
 * ASSUMES: if HMM_VERSION and HMM_RELEASE_DATE are defined they will be
 *          strings
 *
 * EFFECTS: prints verison and disclaimer to stdout
 *
 * ERROR CONDITIONS: none
 *
 * COMMENTS: simple function prints a message with the version and release date
 *           preprocessor variables HMM_VERSION, HMM_RELEASE_DATE
 *
 ******************************************************************************/
void version()
{
   printf("\nhmmSNP version %s  release date %s\n", 
          HMM_VERSION, HMM_RELEASE_DATE);
          
#ifdef __DATE__ 
#ifdef  __TIME__
   printf("  compile date: %s,  %s\n", __DATE__, __TIME__);
#endif
#endif
  
   printf("\n  Algorithm designed by Jin Szatkeiwicz, PhD (jin.szatkiewicz@med.unc.edu)\n");
   printf("   under the supervision of Gary Churchill, PhD (gary.churchill@jax.org)\n\n");
     
   printf("  Programmer (C implementation): \n"
          "    Glen Beane, MS (glen.beane@jax.org), with assistance from Jin\n");
   printf("    A matlab implementation is available from Jin.\n\n");
   
   printf("  For more details see:\n"
          "    Szatkiewicz JP, Beane GL, Ding Y, Hutchins L, Pardo-Manuel de Villena F,\n"
          "      Churchill GA.  \"An imputed genotype resource for the laboratory mouse.\"\n"
          "      Mamm Genome 2008; 19(3):199-208.\n\n"); 



  printf("\n\n  Copyright (c) 2010 The Jackson Laboratory\n\n"
  
         "    This software was developed by Gary Churchill's Lab at The Jackson\n"
         "    Laboratory (see http://research.jax.org/faculty/churchill).\n\n"
 
         "    This is free software: you can redistribute it and/or modify\n"
         "    it under the terms of the GNU General Public License as published by\n"
         "    the Free Software Foundation, either version 3 of the License, or\n"
         "    (at your option) any later version.\n\n"
 
         "    This software is distributed in the hope that it will be useful,\n"
         "    but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
         "    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
         "    GNU General Public License for more details.\n\n"
 
         "    You should have received a copy of the GNU General Public License\n"
         "    along with this software. If not, see <http://www.gnu.org/licenses/>\n\n");


}
