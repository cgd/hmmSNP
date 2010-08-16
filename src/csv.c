/*
   File Name: csv.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
   Date developed: March 2007
   Purpose: general purpose csv library
   Overview: this basic csv library is based on an example from "The Practice
    of Programming" by Brian W. Kernighan and Rob Pike. It may be extended and 
    improved, but their example provides us with enough functionality to start
    with.

   
   
   
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
 
 
/* Portions of this code Copyright (C) 1999 Lucent Technologies 
 * Adapted from 'The Practice of Programming' 
 * by Brian W. Kernighan and Rob Pike 

 * You may use this code for any purpose, as long as you leave the copyright 
 * notice and book citation attached.
 */

 
#include <ctype.h> 
#include <stdlib.h>
#include <string.h>

#include "csv.h"

static void reset(void);
static int endofline(FILE *fsIn, int c);
static int split(void);
static char *advquoted(char *p);


enum { NOMEM = -2};

static char *g_line = NULL;           /* input chars */
static char *g_split_line = NULL;     /* copy used by split() */
static int  g_max_line = 0;           /* size of pgsLine[] and pgsLineCopy[] */
static char **g_field_array = NULL;   /* array of field pointers */
static int  g_max_field = 0;          /* size of g_field_array[] */
static int  g_num_fields = 0;         /* number of feilds */ 

static char g_field_separators[] = ",";

/*******************************************************************************
 * char *csvGetLine(FILE *fsIn)
 *
 * INPUTS: open file stream to input csv file
 *
 * RETURNS: pointer to buffer containing line read from file 
 *
 * ASSUMES:  fsIn is opened for reading
 *
 * EFFECTS:  allocates buffer for line
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
char *csvGetLine(FILE *fsIn)
{
   int i;
   int c;
   
   char *newl;
   char *news;
   
   if (g_line == NULL)
   {
      g_max_line = 1;
      g_max_field = 1;
      
      g_line = malloc(g_max_line);
      g_split_line = malloc(g_max_line);
      g_field_array = malloc(g_max_field * sizeof(g_field_array[0]));
      
      if (g_line == NULL || g_split_line == NULL || g_field_array == NULL)
      {
         reset();
         return NULL;  /* out of memory */
      }
   }
   for (i = 0; (c = getc(fsIn)) != EOF && !endofline(fsIn, c); i++)
   {
      if (i >= g_max_line - 1)
      {
         /* grow line */
         g_max_line *= 2;
         newl = realloc(g_line, g_max_line);
         news = realloc(g_split_line, g_max_line);
         if (newl == NULL || news == NULL)
         {
            reset();
            return NULL;
         }
         
         g_line = newl;
         g_split_line = news;
      }
      g_line[i] = (char)c;
   }
   g_line[i] = '\0';
   if (split() == NOMEM)
   {
      reset();
      return NULL;
   }
   return (c == EOF && i ==0) ? NULL : g_line;
}

/*******************************************************************************
 * char *csvField(int n)
 *
 * INPUTS: int n:  index of the field we want to retrieve
 *
 * RETURNS: pointer to nth field of last line ready by csvGetLine
 *
 * ASSUMES: csvGetLine has been called sucessfully
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
char *csvField(int n)
{
   if (n < 0 || n >= g_num_fields)
   {
      return NULL;
   }
   else
   {
      return g_field_array[n];
   }
}

/*******************************************************************************
 * int   csvNumFields(void)
 *
 * INPUTS: void
 *
 * RETURNS: returns number of fields that were in the last line read 
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
int   csvNumFields(void)
{
   return g_num_fields;
}

/*******************************************************************************
 * int csvGetChrColumn(void) - special function to determine which field c
 *   contains the chromosome informaiton
 *
 * INPUTS:
 *
 * RETURNS: the column number for the column containing the chromosome data, 
 *   if this can not be determined then -1 will be returned
 *
 * ASSUMES:  assumes it is called after the first line is read (header) but 
 *   before the first data line is read
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:  none of the column headers start with "chr", -1 is returned
 *
 * COMMENTS:
 *
 ******************************************************************************/
int csvGetChrColumn(void)
{
   
   for (int i = 0; i < g_num_fields; i++)
   {
     /* check the first "word" to see if it is the chromosome column */  
      if ( (g_field_array[i][0] == 'C' || g_field_array[i][0] == 'c') &&
           (g_field_array[i][1] == 'H' || g_field_array[i][1] == 'h') && 
           (g_field_array[i][2] == 'R' || g_field_array[i][2] == 'r'))
      {
         return i;
      }
   }
 
   return -1;

}

/*******************************************************************************
 * int csvGetID_Column(void)
 *
 * INPUTS:  void
 *
 * RETURNS:  the column number for the column containing the SNP ID
 *
 * ASSUMES:  it is called after the header has been read but before the first 
 *  data line is read by csvGetLine
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:  if the SNP ID column can not be determined then -1 is 
 *    returned
 *
 * COMMENTS:
 *
 ******************************************************************************/
int csvGetID_Column(void)
{
   char *tmp_str;
   for (int i = 0; i < g_num_fields; i++)
   {
   
      tmp_str = strdup(g_field_array[i]);
      for (int j = 0; j < strlen(tmp_str); j++)
      {
         tmp_str[j] = toupper(tmp_str[j]);
      }
      if (strcmp(tmp_str, "SNPID") == 0 || strcmp(tmp_str, "SNP_ID") == 0 
          || strcmp(tmp_str, "SNP.ID") == 0 || strcmp(tmp_str, "SNP ID") == 0)
      {
         free(tmp_str);
         return i;
      }
      free(tmp_str);
   }
 
   return -1;

}

/* internal functions, from The Practice of Programming*/

/*******************************************************************************
 * static void reset(void) - free bufferes and reset sizes
 *
 * INPUTS: void
 *
 * RETURNS: void
 *
 * ASSUMES: nothing 
 *
 * EFFECTS:  frees g_line, g_split_line, g_field_array,  clears g_max_line, 
 *  g_max_field, g_num_fields
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:  anything ANSI-C or later should allow free(NULL), so this function
 *   shouldn't cause any problems if it is called before csvGetLine has 
 *   been called and returned a valid pointer
 *
 ******************************************************************************/
static void reset(void)
{
   free(g_line);
   free(g_split_line);
   free(g_field_array);
   g_line = NULL;
   g_split_line = NULL;
   g_field_array = NULL;
   g_max_line = 0;
   g_max_field = 0;
   g_num_fields = 0;
}

/*******************************************************************************
 * static int endofline(FILE *fsIn, int c) - handles different styles of EOL 
 *   delimiters
 *
 * INPUTS: FILE *fsIn - input file stream
 *         int c - last value returned from a getc call
 *
 * RETURNS: returns value of expression (c == '\r' || c == '\n')
 *
 * ASSUMES: fsIn is opened for input, c is the last value returned by getc
 *
 * EFFECTS:  consumes EOL (\r, \n, \r\n) or EOF from fsIn
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
static int endofline(FILE *fsIn, int c)
{
   int eol;
   
   eol = (c == '\r' || c == '\n');
   if ( c == '\r')
   {
      c = getc(fsIn);
      if (c != '\n' && c != EOF)
      {
         ungetc(c, fsIn); /* read too far; put c back */
      }
   }
   return eol;
}

/*******************************************************************************
 * static int split(void) - splits line on , and handles commans contained in 
 *   double quotes
 *
 * INPUTS: void
 *
 * RETURNS: number of fields
 *
 * ASSUMES: line has been read in with csvGetLine
 *
 * EFFECTS: replaces commas and field enclosing trailing " with null in 
 *   gpsSplitIine, updates g_field_array pointers to point to field strings in 
 *   g_split_line
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
static int split(void)
{
   char *p;
   char **newf;
   char *sep_ptr;
   int  separator;
   
   g_num_fields = 0;
   if (g_line[0] == '\0')
   {
      return 0;
   }
   strcpy(g_split_line, g_line);
   p = g_split_line;
   
   do
   {
      if (g_num_fields >= g_max_field)
      {  
         g_max_field *= 2;
         newf = realloc(g_field_array, g_max_field * sizeof(g_field_array[0]));
         if (newf == NULL)
         {
            return NOMEM;
         }

         g_field_array = newf;
      }
      if (*p == '"')
      {
         sep_ptr = advquoted(++p);
      }
      else
      {
         sep_ptr = p + strcspn(p, g_field_separators);
      }
      separator = sep_ptr[0];
      sep_ptr[0] = '\0';
      g_field_array[g_num_fields++] = p;
      p = sep_ptr + 1;
   
   } while (separator == ',');
   
   return g_num_fields;
}

/*******************************************************************************
 * static char *advquoted(char *p) - 
 *
 * INPUTS: char *p - pointer to quoted field
 *
 * RETURNS: pointer to next separator
 *
 * ASSUMES:
 *
 * EFFECTS: removes enclosing quotes, turns adjacent quotes ("") into "
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
static char *advquoted(char *p)
{
   int i, j;
   
   for (i = 0, j= 0; p[j] != '\0'; i++, j++)
   {
      if (p[j] == '"' && p[++j] != '"')
      {
         int k = strcspn(p+j, g_field_separators);
         memmove(p+i, p+j, k);
         i += k;
         j += k;
         break;
      }
      p[i] = p[j];
   }
   p[i] = '\0';
   return p + j;
}

