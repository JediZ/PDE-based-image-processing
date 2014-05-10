/*
Functions related to image derivative calculation

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

#ifndef _IMAGEDERIVATIVESC_H
#define _IMAGEDERIVATIVESC_H

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

/*
//------------------------------------------------------------------
//--- First order image derivatives using Simoncelli-like operators
//------------------------------------------------------------------
*/
void fstSimoncelli_c(	struct matrixM *Idt,
			struct matrixM *Idx,
			struct matrixM *Idy,
			struct matrixM *It0,
			struct matrixM *It1,
			struct matrixM *temp,
			unsigned int size,	/*Size of the convolution operator, 3 or 5*/
			float *Smoother,
			float *SpatialDerivator,
			float *TemporalDerivator
		);

/*
//-------------------------------------------------------------------
//--- Second order image derivatives using Simoncelli-like operators
//-------------------------------------------------------------------
*/
void sndSimoncelli_c(	struct matrixM *Idxt,
			struct matrixM *Idyt,
			struct matrixM *Idxx,
			struct matrixM *Idyy,
			struct matrixM *Idxy,
			struct matrixM *It0,
			struct matrixM *It1,
			struct matrixM *temp1,
			struct matrixM *temp2,
			unsigned int size,	/*Size of the convolution operator => 5*/
			float *Smoother,
			float *fstSpatialDerivator,
			float *sndSpatialDerivator,
			float *TemporalDerivator
		);



#endif
