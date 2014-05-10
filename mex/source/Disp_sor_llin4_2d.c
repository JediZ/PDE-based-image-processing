/*
function [dU_out varargout] = Disp_sor_llin4_2d(U_in, dU_in, Cu_in, Du, wW, wN, wE, wS, iter)

Calculates disparity using nonlinear diffusion as smoothness term and early linearization of the constancy terms. 

A MEX c-program for Matlab.

INPUT
U_in		=		Disparity.
dU_in		=		Change in disparity.
Cu_in		=		Right hand side for U (A*dU = Cu).
Du		=		Diagonal component of A from the data term.
wW		=		Diffusion weight (west).
wN		=		Diffusion weight (north).
wE		=		Diffusion weight (east).
wS		=		Diffusion weight (south).
iter		=		Number of iterations.
solver		=		Type of solver used:
				    1 = Normal point-wise Gauss-Seidel relaxation.
				    2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.

OUTPUT
dU_out		=		New estimation of dU after 'iter' iterations.
varargout{1}	=		Residual for U. A*dU = Cu => r = Cu - A*dU
varargout{2}	=		Residual for V. A*dU = Cv => r = Cv - A*dV

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

#include "disparitySolvers.h"

/*
//------------------
//--- MAIN (MEX) ---
//------------------
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
	/*Structs*/
	struct matrixM U_in = {0,NULL,NULL,0,0,0,0};	
	struct matrixM dU_in = {0,NULL,NULL,0,0,0,0};	
	struct matrixM dU_out = {0,NULL,NULL,0,0,0,0};
	struct matrixM RU = {0,NULL,NULL,0,0,0,0};	
	/*Constant terms*/
	struct matrixM Cu = {0,NULL,NULL,0,0,0,0};
	struct matrixM Du = {0,NULL,NULL,0,0,0,0};
	/*Regularization diffusion terms*/
	struct matrixM wW = {0,NULL,NULL,0,0,0,0};
	struct matrixM wN = {0,NULL,NULL,0,0,0,0};
	struct matrixM wE = {0,NULL,NULL,0,0,0,0};
	struct matrixM wS = {0,NULL,NULL,0,0,0,0};

	struct mparam Mparam = {0,0,0,0,0,0,0,0,0,0};

	int solver = 0;
	void (*solveME)(struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam);

	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 11) {
		mexErrMsgTxt("dispsor_llin4_2d parameter error: wrong number of input parameters!");
	} 
	/*
	//------------
	//--- U in ---
	//------------
	*/
	/* U_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'U_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	U_in.ndims = mxGetNumberOfDimensions(prhs[0]);
	U_in.dimElems = mxGetDimensions(prhs[0]);
	U_in.data = (float*)mxGetPr( prhs[0] );
	/*
	//------------
	//--- dU in ---
	//------------
	*/
	/* dU_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'dU_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	dU_in.ndims = mxGetNumberOfDimensions(prhs[1]);
	dU_in.dimElems = mxGetDimensions(prhs[1]);
	dU_in.data = (float*)mxGetPr( prhs[1] );
	/*
	//------------
	//--- Cu in ---
	//------------
	*/
	/* Cu must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[2])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'Cu' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Cu.ndims = mxGetNumberOfDimensions(prhs[2]);
	Cu.dimElems = mxGetDimensions(prhs[2]);
	Cu.data = (float*)mxGetPr( prhs[2] );
	/*
	//------------
	//--- Du in ---
	//------------
	*/
	/* Du must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[3])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'Du' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Du.ndims = mxGetNumberOfDimensions(prhs[3]);
	Du.dimElems = mxGetDimensions(prhs[3]);
	Du.data = (float*)mxGetPr( prhs[3] );
	/*
	//----------
	//--- wW ---
	//----------
	*/
	/* wW must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[4])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'wW' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wW.ndims = mxGetNumberOfDimensions(prhs[4]);
	wW.dimElems = mxGetDimensions(prhs[4]);
	wW.data = (float*)mxGetPr( prhs[4] );
	/*
	//----------
	//--- wN ---
	//----------
	*/
	/* wN must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[5])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'wN' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wN.ndims = mxGetNumberOfDimensions(prhs[5]);
	wN.dimElems = mxGetDimensions(prhs[5]);
	wN.data = (float*)mxGetPr( prhs[5] );
	/*
	//----------
	//--- wE ---
	//----------
	*/
	/* wE must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[6])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'wE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wE.ndims = mxGetNumberOfDimensions(prhs[6]);
	wE.dimElems = mxGetDimensions(prhs[6]);
	wE.data = (float*)mxGetPr( prhs[6] );
	/*
	//----------
	//--- wS ---
	//----------
	*/
	/* wS must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[7])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'wS' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wS.ndims = mxGetNumberOfDimensions(prhs[7]);
	wS.dimElems = mxGetDimensions(prhs[7]);
	wS.data = (float*)mxGetPr( prhs[7] );
	/*
	//-----------------------------------
	//--- Number of iterations 'iter' ---
	//-----------------------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[8])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'iter' must be a noncomplex, single-type scalar");
	}
	Mparam.iter = *((float*)mxGetPr( prhs[8] ) );
		/*
	//------------------------------------
	//--- Relaxation parameter 'omega' ---
	//------------------------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[9])){
		mexErrMsgTxt("dispsor_llin4_2d error: 'omega' must be a noncomplex, single-type scalar");
	}
	Mparam.omega = *((float*)mxGetPr( prhs[9] ) );
	/*
	//-------------------
	//--- Solver type ---
	//-------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[10])){
		mexErrMsgTxt("dispsor_elin_2d error: 'solver' must be a noncomplex, single-type scalar");
	}
	solver = (int)*((float*)mxGetPr( prhs[10] ) );
		
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	dU_out.ndims = dU_in.ndims;
	dU_out.dimElems = dU_in.dimElems;

		
	if( nlhs<1 ){
		mexErrMsgTxt("dispsor_llin4_2d insufficient number of outputs. Outputs from this function is 'dU'");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(dU_out.ndims,dU_out.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("dispsor_llin4_2d: error reserving space for output variable 'dU_out'");
		if( (dU_out.data = (float *)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("dispsor_llin4_2d: error obtaining pointer to output variable 'dU_out'");
	}

	if( nlhs>=2 )
	{
		RU.ndims = U_in.ndims;
		RU.dimElems = U_in.dimElems;

		if( (plhs[1] = mxCreateNumericArray(RU.ndims,RU.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("dispsor_llin4_2d: error reserving space for output variable 'RU'");
		if( (RU.data = (float *)mxGetPr( plhs[1] ))==NULL )
			mexErrMsgTxt("dispsor_llin4_2d: error obtaining pointer to output variable 'RU'");
	}

	/*Choose solver*/
	switch (solver)
	{
		case 1:
			solveME = &GS_SOR_llin4_2d;
			break;
		case 2:	
			solveME = &GS_ALR_SOR_llin4_2d;
			break;
		default:
			mexErrMsgTxt("dispsor_llin4_2d: no such solver");
	}

	/*Call solver*/		
	if( Mparam.iter>0 )
	{
		memcpy( dU_out.data, dU_in.data, dU_in.dimElems[0]*dU_in.dimElems[1]*sizeof(float) );
		solveME( &U_in, &dU_out, &Cu, &Du, &wW, &wN, &wE, &wS, Mparam );
	}

}
