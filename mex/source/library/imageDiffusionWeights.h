/*
Functions related to diffusion.

Author: Jarno Ralli
E-mail: jarno@ralli.fi
Web: www.jarnoralli.fi (or www.jarnoralli.com)

If you use this code, please reference (some) of my papers available at http://www.jarnoralli.fi

The program is delivered as it is and the author accepts no responsibility what so ever of its use and/or results.
Errors and suggestions are kindly to be communicated to the author.

(C) Copyright 2011-2014, Jarno Ralli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _DIFFUSION_H
#define _DIFFUSION_H

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "malloc.h"
#include <string.h>


/*
//-------------------------------------------------
//--- Matlab "matrix", which really is a vector ---
//-------------------------------------------------
*/
typedef struct matrixM
{
	unsigned int ndims;
	const int *dimElems;
	float *data;
	unsigned int npadding;
	unsigned int spadding;
	unsigned int epadding;
	unsigned int wpadding;
} matrixM;

typedef struct mparam
{
	float omega;	/*Relaxation parameter for SOR*/
	float alpha;
	float beta;
	float iter;	/*Number of iterations*/
	float fiter;	/*Number of fixed point iterations*/
	float b1;
	float b2;
	float ht;
	float hx;
	float hy;
	
} mparam;

/*
//--------------------------------
//---Calculates diffusion weights 
//--------------------------------
*/
void diffWeights6_2D_c(		struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct matrixM *D,
				float eps );

#endif
