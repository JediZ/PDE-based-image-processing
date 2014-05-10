/*
Functions related to variational disparity calcuation

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

#ifndef _DISPEMINC_H
#define _DISPEMINC_H

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "malloc.h"
#include <string.h>

/*
	    This switch enables compilation of the GS-SOR solver programmed using FPU in GAS for IA32 (Intel)
	    architecture...will not work with X86-64 Intel architecture, for example!
*/
/*#define GASASM*/




/*
//-------------------------------------------------
//--- Matlab "matrix", which really is a vector ---
//-------------------------------------------------
*/
typedef struct matrixM
{
	unsigned int ndims;
	const unsigned int *dimElems;
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
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void	GS_SOR_llin4_2d(	struct matrixM *U,
				struct matrixM *dU,
				struct matrixM *Cu,
				struct matrixM *Du,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Alternating line relaxation Gauss-Seidel (Gauss-Seidel block variant)
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization.
//---Late linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void	GS_ALR_SOR_llin4_2d(	struct matrixM *U,
				struct matrixM *dU,
				struct matrixM *Cu,
				struct matrixM *Du,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions with symmetry.
//---------------------------------------------------------------------------------------
*/
void	GS_SOR_llinsym4_2d(	struct matrixM *U0, 
				struct matrixM *dU0,
				struct matrixM *Cu0,
				struct matrixM *Du0,
				struct matrixM *wW0,
				struct matrixM *wN0,
				struct matrixM *wE0,
				struct matrixM *wS0,
				struct matrixM *U1, 
				struct matrixM *dU1,
				struct matrixM *Cu1,
				struct matrixM *Du1,
				struct matrixM *wW1,
				struct matrixM *wN1,
				struct matrixM *wE1,
				struct matrixM *wS1,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions with symmetry.
//---------------------------------------------------------------------------------------
*/
void	GS_ALR_SOR_llinsym4_2d(	struct matrixM *U0, 
				struct matrixM *dU0,
				struct matrixM *Cu0,
				struct matrixM *Du0,
				struct matrixM *wW0,
				struct matrixM *wN0,
				struct matrixM *wE0,
				struct matrixM *wS0,
				struct matrixM *U1, 
				struct matrixM *dU1,
				struct matrixM *Cu1,
				struct matrixM *Du1,
				struct matrixM *wW1,
				struct matrixM *wN1,
				struct matrixM *wE1,
				struct matrixM *wS1,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------------------
// Calculate residuals r (Ax=b => r=b-Ax) for non-linear (late linearization) case with 4-neighbours
//---------------------------------------------------------------------------------------------------
*/
void Residuals_llin4_2d(	struct matrixM *RU,
				struct matrixM *U,
				struct matrixM *dU,
				struct matrixM *Cu,
				struct matrixM *Du,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin8_2d(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Alternating line relaxation Gauss-Seidel (Gauss-Seidel block variant)
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization.
//---Late linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void	GS_ALR_SOR_llin8_2d(	struct matrixM *U,
				struct matrixM *dU,
				struct matrixM *Cu,
				struct matrixM *Du,
				struct matrixM *wW,
				struct matrixM *wNW,
				struct matrixM *wN,
				struct matrixM *wNE,
				struct matrixM *wE,
				struct matrixM *wSE,
				struct matrixM *wS,
				struct matrixM *wSW,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions with symmetry.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llinsym8_2d(	struct matrixM *U0, 
				struct matrixM *dU0,
				struct matrixM *Cu0,
				struct matrixM *Du0,
				struct matrixM *wW0,
				struct matrixM *wNW0,
				struct matrixM *wN0,
				struct matrixM *wNE0,
				struct matrixM *wE0,
				struct matrixM *wSE0,
				struct matrixM *wS0,
				struct matrixM *wSW0,
				struct matrixM *U1, 
				struct matrixM *dU1,
				struct matrixM *Cu1,
				struct matrixM *Du1,
				struct matrixM *wW1,
				struct matrixM *wNW1,
				struct matrixM *wN1,
				struct matrixM *wNE1,
				struct matrixM *wE1,
				struct matrixM *wSE1,
				struct matrixM *wS1,
				struct matrixM *wSW1,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Succesive over relaxation using both brightness and gradient constraints
//---and 2D first and second order regularization. Non-lizearized constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin4_2d12(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW1,
			struct matrixM *wN1,
			struct matrixM *wE1,
			struct matrixM *wS1,
			struct matrixM *wW2,
			struct matrixM *wN2,
			struct matrixM *wE2,
			struct matrixM *wS2,
			struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions.
//---Inline assembler implementation (AT&T syntax for Linux GAS).
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin4ASM_2d(	float *U,
				float *dU,
				float *Cu,
				float *Du,
				float *wW,
				float *wN,
				float *wE,
				float *wS,
				unsigned int nrows,
				unsigned int ncols,
				unsigned int iterations,
				float omega);

/*
//------------------------------------------------------
// TDMA vertical line solving - east most column solver
//------------------------------------------------------
*/
void eastColumn4(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA vertical line solving - middle columns solver
//------------------------------------------------------
*/
void middleColumn4(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA vertical line solving - west most column solver
//------------------------------------------------------
*/
void westColumn4(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA horizontal line solving - top most row solver
//------------------------------------------------------
*/
void northRow4(		struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA horizontal line solving - middle rows solver
//------------------------------------------------------
*/
void middleRow4(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA horizontal line solving - bottom most row solver
//------------------------------------------------------
*/
void southRow4(		struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA vertical line solving - east most column solver
//------------------------------------------------------
*/
void eastColumn8(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA vertical line solving - middle columns solver
//------------------------------------------------------
*/
void middleColumn8(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA vertical line solving - west most column solver
//------------------------------------------------------
*/
void westColumn8(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA horizontal line solving - top most row solver
//------------------------------------------------------
*/
void northRow8(		struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA horizontal line solving - middle rows solver
//------------------------------------------------------
*/
void middleRow8(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);

/*
//------------------------------------------------------
// TDMA horizontal line solving - bottom most row solver
//------------------------------------------------------
*/
void southRow8(		struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows);


#endif

