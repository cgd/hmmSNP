/*
   File Name: temp_file_cleanup.c
   System: hmmSNP
   Programmer: Glen Beane
   Date developed: 2/13/08
   Purpose: track and clean temporary files
   Comments:
  
   I found that I had a case where I created a temporary file in one function
   but I don't want to clean that file up until just before the program exits,
   so I decided to make a general solution that I could use to clean up any 
   temporary files created anywhere in the program with one function call
   which is to be called before exiting.
   
   
   
   
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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "constants.h"

static char *g_filename_buffers[TMP_FILE_NAME_BUFFER];
static int g_count = 0;



/*******************************************************************************
 * register_temp_file - register a temporary file for deletion prior to 
 *    program termination
 *
 * INPUTS:    char *filename - path to temporary file
 *
 * RETURNS:   if file name buffer array is full return EHMM_ENOMEM, 0 if no error
 *
 * ASSUMES:   filename is a valid path
 *
 * EFFECTS:   allocates memory to store a copy of filename
 *
 * ERROR CONDITIONS:  returns error if buffer is full
 *
 * COMMENTS:
 *
 ******************************************************************************/
int register_temp_file(char *filename)
{
   if (g_count == TMP_FILE_NAME_BUFFER) {
      return EHMM_NOMEM;
   }
   
   g_filename_buffers[g_count++] = strdup(filename);
   
   return 0;


}

/*******************************************************************************
 * temp_file_cleanup - cleanup all temp files listed in g_filename_buffers
 *
 * INPUTS: void
 *
 * RETURNS: void
 *
 * ASSUMES: 
 *
 * EFFECTS:  will unlink all g_count files in g_filename_buffers
 *
 * ERROR CONDITIONS: does not report any unlink errors
 *
 * COMMENTS:
 *
 ******************************************************************************/
void temp_file_cleanup(void)
{
   /* call unlink on each file name in our list of file names,  free each 
      file name buffer */
   for (int i = 0; i < g_count; i++) {
      unlink(g_filename_buffers[i]);
      free(g_filename_buffers[i]);
   }
   g_count = 0;
}
