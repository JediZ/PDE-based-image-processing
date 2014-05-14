/*
Functions related to levelset:s

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

#ifndef _LEVEL_SETS_H
#define _LEVEL_SETS_H

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "malloc.h"
#include "float.h"
#include <string.h>
#include <xmmintrin.h>
#include <omp.h>
/*#include <pmmintrin.h>*/

/*Maximum buffer size*/
#define MAX_BUF_SIZE 2048

/*	SWITCHES ARE AS FOLLOWS:
	SSE and GASASM = assembler + SSE code for blurredSignFunction()
	SSE = C SSE code for blurredSignFunction()
	GASASM alone = C code for blurredSignFunction()
	without SSE or GASASM = C code for blurredSignFunction()
	
	A word of warning!
	SSE has been tested in both IA32 and X86-64-bit Intel architectures with GCC and should work O.K!
	SSE and GASASM works only in IA32!!!
*/
#define SSE
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
//---Additive Operator Splitting (AOS) scheme for Active Contour Models
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 4-neighbours
//---------------------------------------------------------------------------------------
*/
void	CV_AOS_4_2d(	struct matrixM *PHI_out,
			struct matrixM *PHI_in,
			struct matrixM *D_in,		/*Data*/
			struct matrixM *GradNorm_in,	/*Derived from N = (nabla O)/|nabla O|*/
			struct matrixM *Diff_in,	/*Diffusivity, inside DIV*/
			float tau,
			float nu );
			
/*
//---------------------------------------------------------------------------------------
//---Additive Operator Splitting (AOS) scheme for Active Contour Models
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 4-neighbours
//--- Slightly optimized version using Open MP when solving for more than 1 level-set equation
//---------------------------------------------------------------------------------------
*/
void	CV_AOSOMP_4_2d(	struct matrixM *PHI_out,
			struct matrixM *PHI_in,
			struct matrixM *D_in,		/*Data*/
			struct matrixM *GradNorm_in,	/*Derived from N = (nabla O)/|nabla O|*/
			struct matrixM *Diff_in,	/*Diffusivity, inside DIV*/
			float tau,
			float nu );

/*
//---------------------------------------------------------------------------------------
//---Additive Operator Splitting (AOS) scheme for Active Contour Models
//---With convective derivative term
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 4-neighbours
//---------------------------------------------------------------------------------------
*/
void	AC_AOS_4_2d(	struct matrixM *PHI_out,
			struct matrixM *PHI_in,
			struct matrixM *D_in,		/*Data*/
			struct matrixM *GradNorm_in,	/*Derived from N = (nabla PHI)/|nabla PHI|*/
			struct matrixM *Diff_in,	/*Diffusivity, inside DIV*/
			float tau,			/*tau = delta t = time step*/
			float nu);/*, 			/*diffusivity coefficient*/
/*TODO: Fix this!	float Treinit);*/			/*t=0:0.25:Treinit...reinitialisation steps*/

/*
//---------------------------------------------------------------------------------------
//---Reinitialization of the level-set function to a signed distance function
//---------------------------------------------------------------------------------------
*/
void	reinit(		struct matrixM *PHI_inout,
			float T);
			
/*
//-------------------------------------------------------------------
//--- TDMA vertical line solving - 4 neighbours - harmonic averaging
//--- Slightly optimized version using Open MP
//-------------------------------------------------------------------
*/
void CV_TDMA_Column4_omp(	struct matrixM *PHI_out, 
				struct matrixM *PHI_in, 
				struct matrixM *D_in, 
				struct matrixM *GradNorm_in, 
				struct matrixM *Diff_in, 
				float tau, 
				float nu, 
				int ncols, 
				int nrows, 
				int nframes);
				

/*
//---------------------------------------------------------------------
//--- TDMA horizontal line solving - 4 neighbours - harmonic averaging
//--- Slightly optimized version using Open MP
//---------------------------------------------------------------------
*/
void CV_TDMA_Row4_omp(	struct matrixM *PHI_out, 
			struct matrixM *PHI_in, 
			struct matrixM *D_in, 
			struct matrixM *GradNorm_in, 
			struct matrixM *Diff_in, 
			float tau, 
			float nu, 
			int ncols, 
			int nrows, 
			int nframes);

/*
//-------------------------------------------------------------------
//--- TDMA vertical line solving - 4 neighbours - harmonic averaging
//-------------------------------------------------------------------
*/
void CV_TDMA_Column4(		struct matrixM *PHI_out, 
				struct matrixM *PHI_in, 
				struct matrixM *D_in, 
				struct matrixM *GradNorm_in, 
				struct matrixM *Diff_in, 
				float tau, 
				float nu, 
				float *cp,
				float *dp,
				int ncols, 
				int nrows, 
				int nframes);
				

/*
//---------------------------------------------------------------------
//--- TDMA horizontal line solving - 4 neighbours - harmonic averaging
//---------------------------------------------------------------------
*/
void CV_TDMA_Row4(	struct matrixM *PHI_out, 
			struct matrixM *PHI_in, 
			struct matrixM *D_in, 
			struct matrixM *GradNorm_in, 
			struct matrixM *Diff_in, 
			float tau, 
			float nu, 
			float *cp,
			float *dp,
			int ncols, 
			int nrows, 
			int nframes);

/*
//------------------------------------------------------------------------------------------
//--- TDMA vertical line solving - 4 neighbours - harmonic averaging - active contour model
//------------------------------------------------------------------------------------------
*/
void AC_TDMA_column4(	struct matrixM *PHI_out, 
			struct matrixM *PHI_in, 
			struct matrixM *D_in, 
			struct matrixM *GradNorm_in, 
			struct matrixM *Diff_in, 
			float tau, 
			float nu, 
			float *cp,
			float *dp,
			int ncols, 
			int nrows, 
			int nframes);

/*
//--------------------------------------------------------------------------------------------
//--- TDMA horizontal line solving - 4 neighbours - harmonic averaging - active contour model
//--------------------------------------------------------------------------------------------
*/
void AC_TDMA_row4(	struct matrixM *PHI_out, 
			struct matrixM *PHI_in, 
			struct matrixM *D_in, 
			struct matrixM *GradNorm_in, 
			struct matrixM *Diff_in, 
			float tau, 
			float nu, 
			float *cp,
			float *dp,
			int ncols, 
			int nrows, 
			int nframes);

/*
//-------------------------
//--- Vertical convolution 
//-------------------------
*/
void VerticalConv(float *Result, float *I, float *operator, unsigned int rows, unsigned int cols, unsigned int frames);

/*
//----------------------------
//--- Horizontal convolution 
//----------------------------
*/
void HorizontalConv(float *Result, float *I, float *operator, unsigned int rows, unsigned int cols, unsigned int frames);

/*
//----------------------------------------------------------
//--- Blurred sign function
//--- Implementation depends on the switches SSE and GASASM
//----------------------------------------------------------
*/
void blurredSignFunction(float *S, float *PHI, float *PHIx, float *PHIy, unsigned int rows, unsigned int cols, unsigned int frames);

/*
//-------------------------------------------------------------------------------------
//--- Godunov's method for discretising the hyperbolic $S(\Phi_0) | \Delta x |$ term
//-------------------------------------------------------------------------------------
*/
void godunovUpwind(	float *PHI2dx_out,
			float *PHI2dy_out,
			float *a_in,
			float *PHI_in,
			unsigned int rows,
			unsigned int cols,
			unsigned int frames);

/*
//-----------------------------------------
//--- Fast 'unprecise' inverse square root
//-----------------------------------------
*/
float InvSqrt (float x);

#endif
