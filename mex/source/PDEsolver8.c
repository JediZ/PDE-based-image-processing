/*
function Xnew = PDEsolver( Xold, TRACE, B, wW, wNW, wN, wNE, wE, wSE, wS, wSW, iter, omega, solver )

A MEX c-program for Matlab.

INPUT
Xold	=	Old solution.
wW	=	West diffusion weights.
wN	=	North diffusion weights.
wE	=	East diffusion weights.
wS	=	South diffusion weights.
iter	=	Number of iteration cycles.
omega	=	Relaxation parameter.
solver		=		Type of solver used:
			      1 = Normal point-wise Gauss-Seidel relaxation.
			      2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.

OUTPUT
Xnew	=	New solution.


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

#include "pdeSolvers.h"

/*
//------------------
//--- MAIN (MEX) ---
//------------------
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
	/*Structs*/
	struct matrixM Xnew = {0,NULL,NULL,0,0,0,0};
	struct matrixM Xold = {0,NULL,NULL,0,0,0,0};
	struct matrixM TRACE = {0,NULL,NULL,0,0,0,0};
	struct matrixM B = {0,NULL,NULL,0,0,0,0};
	
	struct matrixM wW = {0,NULL,NULL,0,0,0,0};
	struct matrixM wN = {0,NULL,NULL,0,0,0,0};
	struct matrixM wNW = {0,NULL,NULL,0,0,0,0};
	struct matrixM wNE = {0,NULL,NULL,0,0,0,0};
	struct matrixM wE = {0,NULL,NULL,0,0,0,0};
	struct matrixM wS = {0,NULL,NULL,0,0,0,0};
	struct matrixM wSW = {0,NULL,NULL,0,0,0,0};
	struct matrixM wSE = {0,NULL,NULL,0,0,0,0};

	struct mparam Mparam = {0,0,0,0,0,0,0,0,0,0};
	int solver = 0;
	int elems = -1;
	void (*solveME)(struct matrixM *Xnew,
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
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 14) {
		mexErrMsgTxt("error: wrong number of input parameters!");
	} 
	/*
	//------------
	//--- Xold in ---
	//------------
	*/
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("error: 'Xold' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Xold.ndims = mxGetNumberOfDimensions(prhs[0]);
	Xold.dimElems = mxGetDimensions(prhs[0]);
	Xold.data = (float*)mxGetPr( prhs[0] );
	/*
	//------------
	//--- TRACE ---
	//------------
	*/
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("error: 'TRACE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	TRACE.ndims = mxGetNumberOfDimensions(prhs[1]);
	TRACE.dimElems = mxGetDimensions(prhs[1]);
	TRACE.data = (float*)mxGetPr( prhs[1] );
	/*
	//------------
	//--- B ---
	//------------
	*/
	if (!mxIsSingle(prhs[2])){
		mexErrMsgTxt("error: 'B' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	B.ndims = mxGetNumberOfDimensions(prhs[2]);
	B.dimElems = mxGetDimensions(prhs[2]);
	B.data = (float*)mxGetPr( prhs[2] );
	/*
	//------------
	//--- wW in ---
	//------------
	*/
	if (!mxIsSingle(prhs[3])){
		mexErrMsgTxt("error: 'wW' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wW.ndims = mxGetNumberOfDimensions(prhs[3]);
	wW.dimElems = mxGetDimensions(prhs[3]);
	wW.data = (float*)mxGetPr( prhs[3] );
	/*
	//------------
	//--- wNW in ---
	//------------
	*/
	if (!mxIsSingle(prhs[4])){
		mexErrMsgTxt("error: 'wNW' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wNW.ndims = mxGetNumberOfDimensions(prhs[4]);
	wNW.dimElems = mxGetDimensions(prhs[4]);
	wNW.data = (float*)mxGetPr( prhs[4] );
	/*
	//------------
	//--- wN in ---
	//------------
	*/
	if (!mxIsSingle(prhs[5])){
		mexErrMsgTxt("error: 'wN' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wN.ndims = mxGetNumberOfDimensions(prhs[5]);
	wN.dimElems = mxGetDimensions(prhs[5]);
	wN.data = (float*)mxGetPr( prhs[5] );
	/*
	//------------
	//--- wNE in ---
	//------------
	*/
	if (!mxIsSingle(prhs[6])){
		mexErrMsgTxt("error: 'wNE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wNE.ndims = mxGetNumberOfDimensions(prhs[6]);
	wNE.dimElems = mxGetDimensions(prhs[6]);
	wNE.data = (float*)mxGetPr( prhs[6] );
	/*
	//------------
	//--- wE in ---
	//------------
	*/
	if (!mxIsSingle(prhs[7])){
		mexErrMsgTxt("error: 'wE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wE.ndims = mxGetNumberOfDimensions(prhs[7]);
	wE.dimElems = mxGetDimensions(prhs[7]);
	wE.data = (float*)mxGetPr( prhs[7] );
	/*
	//------------
	//--- wSE in ---
	//------------
	*/
	if (!mxIsSingle(prhs[8])){
		mexErrMsgTxt("error: 'wSE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wSE.ndims = mxGetNumberOfDimensions(prhs[8]);
	wSE.dimElems = mxGetDimensions(prhs[8]);
	wSE.data = (float*)mxGetPr( prhs[8] );
	/*
	//------------
	//--- wS in ---
	//------------
	*/
	if (!mxIsSingle(prhs[9])){
		mexErrMsgTxt("error: 'wS' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wS.ndims = mxGetNumberOfDimensions(prhs[9]);
	wS.dimElems = mxGetDimensions(prhs[9]);
	wS.data = (float*)mxGetPr( prhs[9] );
	/*
	//------------
	//--- wSW in ---
	//------------
	*/
	if (!mxIsSingle(prhs[10])){
		mexErrMsgTxt("error: 'wS' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wSW.ndims = mxGetNumberOfDimensions(prhs[10]);
	wSW.dimElems = mxGetDimensions(prhs[10]);
	wSW.data = (float*)mxGetPr( prhs[10] );
	/*
	//-----------------------------------
	//--- Number of iterations 'iter' ---
	//-----------------------------------
	*/
	if (!mxIsSingle(prhs[11])){
		mexErrMsgTxt("error: 'iter' must be a noncomplex, single-type scalar");
	}
	Mparam.iter = *((float*)mxGetPr( prhs[11] ) );
	/*
	//-----------------------------------
	//--- Relaxation parameter 'omega'
	//-----------------------------------
	*/
	if (!mxIsSingle(prhs[12])){
		mexErrMsgTxt("error: 'omage' must be a noncomplex, single-type scalar");
	}
	Mparam.omega = *((float*)mxGetPr( prhs[12] ) );
	/*
	//-------------------
	//--- Solver type ---
	//-------------------
	*/
	if (!mxIsSingle(prhs[13])){
		mexErrMsgTxt("error: 'solver' must be a noncomplex, single-type scalar");
	}
	solver = (int)*((float*)mxGetPr( prhs[13] ) );

		
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	Xnew.ndims = Xold.ndims;
	Xnew.dimElems = Xold.dimElems;
	
	if( nlhs<1 ){
		mexErrMsgTxt("error insufficient number of outputs.");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(Xnew.ndims, Xnew.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("error: error reserving space for output variable 'Xnew'");
		if( (Xnew.data = (float *)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("error: error obtaining pointer to output variable 'Xnew'");
	}

	/*Choose solver*/
	switch (solver)
	{
		case 1:
			solveME = &GS_SOR_8_2d;
			break;
		case 2:	
			solveME = &GS_ALR_SOR_8_2d;
			break;
		case 3:	
			break;
		default:
			mexErrMsgTxt("error: no such solver");
	}
	
	elems = Xnew.dimElems[0]*Xnew.dimElems[1];
		
	if(Xnew.ndims>2)
	  elems *= Xnew.dimElems[2];

	memcpy( Xnew.data, Xold.data, elems*sizeof(float) );
	solveME(	&Xnew,
			&TRACE,
			&B,
			&wW,
			&wNW,
			&wN,
			&wNE,
			&wE,
			&wSE,
			&wS,
			&wSW,
			Mparam);

}
