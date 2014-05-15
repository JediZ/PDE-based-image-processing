 /*
Functions related to PDE solvers.
 
Solvers search for a solution to Ax = B, such that
-TRACE = TRACE(A)
-x is the solution

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

#ifndef _PDESOLVER_H
#define _PDESOLVER_H

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
//---Point-wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_4_2d(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam);
			
/*
//---------------------------------------------------------------------------------------
//---Point-wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_8_2d(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
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
//---Alternating Line Relaxation (ALR) scheme
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 4-neighbours
//---------------------------------------------------------------------------------------
*/
void GS_ALR_SOR_4_2d(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Alternating Line Relaxation (ALR) scheme
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 8-neighbours
//---------------------------------------------------------------------------------------
*/
void GS_ALR_SOR_8_2d(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
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
//--------------------------------------------
// TDMA vertical line solving - east column
//--------------------------------------------
*/
void TDMA_ecolumn_ALR_4(struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int frames,
			float omega);

/*
//--------------------------------------------
// TDMA vertical line solving - middle columns
//--------------------------------------------
*/
void TDMA_mcolumn_ALR_4(struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int frames,
			float omega);

/*
//------------------------------------------
// TDMA vertical line solving - west column
//------------------------------------------
*/
void TDMA_wcolumn_ALR_4(struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int frames,
			float omega);
/*
//------------------------------------------
// TDMA horizontal line solving - north row
//------------------------------------------
*/
void TDMA_nrow_ALR_4(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega);

/*
//-------------------------------------------
// TDMA horizontal line solving - middle rows
//-------------------------------------------
*/
void TDMA_mrow_ALR_4(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega);

/*
//------------------------------------------
// TDMA horizontal line solving - south row
//------------------------------------------
*/
void TDMA_srow_ALR_4(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega);
/*
//----------------------------
// TDMA vertical line solving
//----------------------------
*/
void TDMAcolumn_ALR_8(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega);
			
/*
//------------------------------
// TDMA horizontal line solving 
//------------------------------
*/
void TDMArow_ALR_8(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega);
			
/*
//------------------------------------------------------
// TDMA horizontal line solving - 8 neighbours
// !!!NOT IMPLEMENTED!!!
//------------------------------------------------------
*/
void TDMA_row8(		struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes);
#endif
