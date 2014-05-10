/*
Functions related to RANSAC

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

#ifndef _RANSACC_H
#define _RANSACC_H

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "malloc.h"
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>


/*
//-------------------------------------------------
//--- Matlab "matrix", which really is a vector ---
//-------------------------------------------------
*/
typedef struct matrixM
{
	unsigned int ndims;
	/*unsigned int *dimElems;*/
	const mwSize *dimElems;
	float *data;
	unsigned int npadding;
	unsigned int spadding;
	unsigned int epadding;
	unsigned int wpadding;
} matrixM;

void RANSAC( struct matrixM *A_in,				/*Data to be fitted A_in * model_out = B_in*/
	     struct matrixM *B_in,				/*Target*/
	     struct matrixM *model_out,				/*The best model found by RANSAC*/
	     struct matrixM *model_in,				/*A model participating in the competition*/
	     struct matrixM *error_out,				/*Residual (or fitting error) of each data for the best model*/
	     /*function pointer as an argument*/
	     void (*createModel) (	struct matrixM *,	/* A */
					struct matrixM *,	/* B */
					struct matrixM *,	/* model_out */
					struct matrixM *,	/* model_in */
					struct matrixM *,	/* error_out */
					unsigned int *,		/* rand_set */
					unsigned int ),		/* n */
	     float err_thr, 					/*Error threshold for determining if the data fits the model*/
	     float min_set_size,	 			/*Minimum number of close data (percentage; 0.0-1.0) values required to assert that a model fits well to data*/
	     unsigned int iter, 				/*Number of iterations*/
	     unsigned int n );					/*Number of data used for generating the model*/

void randVect(unsigned int *vect, unsigned int min, unsigned int max, unsigned int length);

#endif
