/*
   File Name: rand.h
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
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

#ifndef RAND_H
#define RAND_H

double uniformRandom(void);
double gammaDev(int iShape);
double beta(int a, int b);
double rexp(double lambda);

#endif
