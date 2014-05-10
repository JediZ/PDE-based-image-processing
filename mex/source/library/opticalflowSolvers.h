/*
Functions related to variational optical-flow calculation

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

#ifndef _OPFLOWEMINC_H
#define _OPFLOWEMINC_H

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
	float alpha;	/*Smoothness (diffusion) factor*/
	float iter;	/*Number of iterations*/
	float fiter;	/*Number of fixed point iterations*/
	float b1;	/*Brightness constancy co-efficient*/
	float b2;	/*Gradient constancy co-efficient*/
	float ht;
	float hx;
	float hy;
	float beta;
	
} mparam;

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Early lizearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_elin4_2d(	struct matrixM *U, 
			struct matrixM *V, 
			struct matrixM *M,
			struct matrixM *Cu,
			struct matrixM *Cv,
			struct matrixM *Du,
			struct matrixM *Dv,
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
//---Early linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void	GS_ALR_SOR_elin4_2d(	struct matrixM *U, 
				struct matrixM *V, 
				struct matrixM *M,
				struct matrixM *Cu,
				struct matrixM *Cv,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------------------
// Calculate residuals r (Ax=b => r=b-Ax) for early linearization case with 4 neighbours.
//---------------------------------------------------------------------------------------------------
*/
void Residuals_elin4_2d(	struct matrixM *RU, 
				struct matrixM *RV, 
				struct matrixM *U,
				struct matrixM *V,
				struct matrixM *M,
				struct matrixM *Cu,
				struct matrixM *Cv,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam);

/*
//----------------------------------------------------------------------
// Calculates LHS (Left hand size of Ax=b) for early linearization case
//----------------------------------------------------------------------
*/
void LHS_elin4_2d(	struct matrixM *AU, 
			struct matrixM *AV, 
			struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
			struct matrixM *Du,
			struct matrixM *Dv,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel successive over relaxation with 2D (spatial) regularization
//---Late linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin4_2d(	struct matrixM *U, 
			struct matrixM *V, 
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
			struct matrixM *Cu,
			struct matrixM *Cv,
			struct matrixM *Du,
			struct matrixM *Dv,
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
				struct matrixM *V, 
				struct matrixM *dU,
				struct matrixM *dV,
				struct matrixM *M,
				struct matrixM *Cu,
				struct matrixM *Cv,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------------------
// Calculates residuals r (Ax=b => r=b-Ax) for late linearization case with 4-neighbours
//---------------------------------------------------------------------------------------------------
*/
void Residuals_llin4_2d(	struct matrixM *RU,
				struct matrixM *RV,
				struct matrixM *U, 
				struct matrixM *V, 
				struct matrixM *dU,
				struct matrixM *dV,
				struct matrixM *M,
				struct matrixM *Cu,
				struct matrixM *Cv,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------
// Calculates LHS (Left hand size of Ax=b) for late linearization case
//---------------------------------------------------------------------
*/
void LHS_llin4_2d(		struct matrixM *AU,
				struct matrixM *AV,
				struct matrixM *U, 
				struct matrixM *V, 
				struct matrixM *dU,
				struct matrixM *dV,
				struct matrixM *M,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam);

/*
//----------------------------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel successive over relaxation with 2D (spatial) regularization for symmetric case
//---Late linearization of constancy assumptions with symmetry.
//----------------------------------------------------------------------------------------------------------
*/
void GS_SOR_llinsym4_2d(	struct matrixM *U0, 
				struct matrixM *V0, 
				struct matrixM *dU0,
				struct matrixM *dV0,
				struct matrixM *M0,
				struct matrixM *Cu0,
				struct matrixM *Cv0,
				struct matrixM *Du0,
				struct matrixM *Dv0,
				struct matrixM *wW0,
				struct matrixM *wN0,
				struct matrixM *wE0,
				struct matrixM *wS0,
				struct matrixM *U1, 
				struct matrixM *V1, 
				struct matrixM *dU1,
				struct matrixM *dV1,
				struct matrixM *M1,
				struct matrixM *Cu1,
				struct matrixM *Cv1,
				struct matrixM *Du1,
				struct matrixM *Dv1,
				struct matrixM *wW1,
				struct matrixM *wN1,
				struct matrixM *wE1,
				struct matrixM *wS1,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Alternating line relaxation Gauss-Seidel (Gauss-Seidel block variant)
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization.
//---Late linearization of constancy assumptions imposing symmetry
//---------------------------------------------------------------------------------------
*/
void GS_ALR_SOR_llinsym4_2d(	struct matrixM *U0, 
				struct matrixM *V0, 
				struct matrixM *dU0,
				struct matrixM *dV0,
				struct matrixM *M0,
				struct matrixM *Cu0,
				struct matrixM *Cv0,
				struct matrixM *Du0,
				struct matrixM *Dv0,
				struct matrixM *wW0,
				struct matrixM *wN0,
				struct matrixM *wE0,
				struct matrixM *wS0,
				struct matrixM *U1, 
				struct matrixM *V1, 
				struct matrixM *dU1,
				struct matrixM *dV1,
				struct matrixM *M1,
				struct matrixM *Cu1,
				struct matrixM *Cv1,
				struct matrixM *Du1,
				struct matrixM *Dv1,
				struct matrixM *wW1,
				struct matrixM *wN1,
				struct matrixM *wE1,
				struct matrixM *wS1,
				struct mparam Mparam);

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel successive over relaxation with 2D (spatial) regularization
//---Late linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin8_2d(	struct matrixM *U, 
			struct matrixM *V, 
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
			struct matrixM *Cu,
			struct matrixM *Cv,
			struct matrixM *Du,
			struct matrixM *Dv,
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
				struct matrixM *V, 
				struct matrixM *dU,
				struct matrixM *dV,
				struct matrixM *M,
				struct matrixM *Cu,
				struct matrixM *Cv,
				struct matrixM *Du,
				struct matrixM *Dv,
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
//*********************************************************************************
// LINE SOLVERS FOR EARLY LINEARIZATION CASE
//*********************************************************************************
*/

/*
//------------------------------------------------------
// TDMA vertical line solving - west most columns solver
//------------------------------------------------------
*/
void westColumn_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
void middleColumn_elin4(struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
void eastColumn_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
// TDMA horizontal line solving - north most row solver
//------------------------------------------------------
*/
void northRow_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
void middleRow_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
// TDMA horizontal line solving - south most rows solver
//------------------------------------------------------
*/
void southRow_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
//*********************************************************************************
// LINE SOLVERS FOR LATE LINEARIZATION CASE
//*********************************************************************************
*/

/*
//------------------------------------------------------
// TDMA vertical line solving - west most columns solver
//------------------------------------------------------
*/
void westColumn_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
void middleColumn_llin4(struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
void eastColumn_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
// TDMA horizontal line solving - north most row solver
//------------------------------------------------------
*/
void northRow_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
void middleRow_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
// TDMA horizontal line solving - south most rows solver
//------------------------------------------------------
*/
void southRow_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
// TDMA vertical line solving - west most columns solver
//------------------------------------------------------
*/
void westColumn_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
void middleColumn_llin8(struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
// TDMA vertical line solving - east most column solver
//------------------------------------------------------
*/
void eastColumn_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
// TDMA horizontal line solving - north most row solver
//------------------------------------------------------
*/
void northRow_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
void middleRow_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
// TDMA horizontal line solving - south most rows solver
//------------------------------------------------------
*/
void southRow_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
