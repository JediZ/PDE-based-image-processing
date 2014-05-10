/*
function [U_out V_out varargout] = oflow_sor_elin4_2d(U_in, V_in, M, Cu_in, Cv_in, Du, Dv, wW, wN, wE, wS, iter)

Calculates optical flow using nonlinear diffusion as smoothness term and early linearization of the constancy terms. 

A MEX c-program for Matlab.

INPUT
U_in		=		Optical flow, horizontal component.
V_in		=		Optical flow, vertical component.
M		=		Common multiplier
Cu_in		=		Right hand side for U (Au*U = Cu - M*V).
Cv_in		=		Right hand size for V (Av*V = Cv - M*U).
Du		=		Diagonal component of A from the data term.
Dv		=		Diagonal component of A from the data term.
wW		=		Diffusion weight (west).
wN		=		Diffusion weight (north).
wE		=		Diffusion weight (east).
wS		=		Diffusion weight (south).
iter		=		Number of iterations.
solver		=		Type of solver used:
			      1 = Normal point-wise Gauss-Seidel relaxation.
			      2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.

OUTPUT
U_out		=		New estimation of U after 'iter' iterations.
V_out		=		New estimation of V after 'iter' iterations.
varargout{1}	=		Residual for U. Au*U = Cu - M*V => r = Cu - M*V - Au*U
varargout{2}	=		Resitual for V. Av*V = Cv - M*U => r = Cv - M*U - Av*V


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

#include "opticalflowSolvers.h"

/*
//------------------
//--- MAIN (MEX) ---
//------------------
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
	/*Structs*/
	struct matrixM U_in = {0,NULL,NULL,0,0,0,0};	
	struct matrixM V_in = {0,NULL,NULL,0,0,0,0};	
	struct matrixM U_out = {0,NULL,NULL,0,0,0,0};	
	struct matrixM V_out = {0,NULL,NULL,0,0,0,0};
	struct matrixM RU = {0,NULL,NULL,0,0,0,0};	
	struct matrixM RV = {0,NULL,NULL,0,0,0,0};
	/*Constant terms*/
	struct matrixM M = {0,NULL,NULL,0,0,0,0};
	struct matrixM Cu = {0,NULL,NULL,0,0,0,0};
	struct matrixM Cv = {0,NULL,NULL,0,0,0,0};
	struct matrixM Du = {0,NULL,NULL,0,0,0,0};
	struct matrixM Dv = {0,NULL,NULL,0,0,0,0};
	/*Regularization diffusion terms*/
	struct matrixM wW = {0,NULL,NULL,0,0,0,0};
	struct matrixM wN = {0,NULL,NULL,0,0,0,0};
	struct matrixM wE = {0,NULL,NULL,0,0,0,0};
	struct matrixM wS = {0,NULL,NULL,0,0,0,0};

	struct mparam Mparam = {0,0,0,0,0,0,1,1,1};
	
	int solver = 0;
	void (*solveME)(struct matrixM *U, 
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
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 14) {
		mexErrMsgTxt("opsor_elin_2d parameter error: wrong number of input parameters!");
	} 
	/*
	//------------
	//--- U in ---
	//------------
	*/
	/* U_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("opsor_elin_2d error: 'U_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	U_in.ndims = mxGetNumberOfDimensions(prhs[0]);
	U_in.dimElems = mxGetDimensions(prhs[0]);
	U_in.data = (float*)mxGetPr( prhs[0] );
	/*
	//------------
	//--- V in ---
	//------------
	*/
	/* V_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("opsor_elin_2d error: 'V_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	V_in.ndims = mxGetNumberOfDimensions(prhs[1]);
	V_in.dimElems = mxGetDimensions(prhs[1]);
	V_in.data = (float*)mxGetPr( prhs[1] );
	/*
	//------------
	//--- M in ---
	//------------
	*/
	/* M must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[2])){
		mexErrMsgTxt("opsor_elin_2d error: 'M' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	M.ndims = mxGetNumberOfDimensions(prhs[2]);
	M.dimElems = mxGetDimensions(prhs[2]);
	M.data = (float*)mxGetPr( prhs[2] );
	/*
	//------------
	//--- Cu in ---
	//------------
	*/
	/* Cu must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[3])){
		mexErrMsgTxt("opsor_elin_2d error: 'Cu' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Cu.ndims = mxGetNumberOfDimensions(prhs[3]);
	Cu.dimElems = mxGetDimensions(prhs[3]);
	Cu.data = (float*)mxGetPr( prhs[3] );
	/*
	//------------
	//--- Cv in ---
	//------------
	*/
	/* Cv must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[4])){
		mexErrMsgTxt("opsor_elin_2d error: 'Cv' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Cv.ndims = mxGetNumberOfDimensions(prhs[4]);
	Cv.dimElems = mxGetDimensions(prhs[4]);
	Cv.data = (float*)mxGetPr( prhs[4] );
	/*
	//------------
	//--- Du in ---
	//------------
	*/
	/* Du must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[5])){
		mexErrMsgTxt("opsor_elin_2d error: 'Du' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Du.ndims = mxGetNumberOfDimensions(prhs[5]);
	Du.dimElems = mxGetDimensions(prhs[5]);
	Du.data = (float*)mxGetPr( prhs[5] );
	/*
	//------------
	//--- Dv in ---
	//------------
	*/
	/* Dv must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[6])){
		mexErrMsgTxt("opsor_elin_2d error: 'Dv' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Dv.ndims = mxGetNumberOfDimensions(prhs[6]);
	Dv.dimElems = mxGetDimensions(prhs[6]);
	Dv.data = (float*)mxGetPr( prhs[6] );
	/*
	//----------
	//--- wW ---
	//----------
	*/
	/* wW must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[7])){
		mexErrMsgTxt("opsor_elin_2d error: 'wW' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wW.ndims = mxGetNumberOfDimensions(prhs[7]);
	wW.dimElems = mxGetDimensions(prhs[7]);
	wW.data = (float*)mxGetPr( prhs[7] );
	/*
	//----------
	//--- wN ---
	//----------
	*/
	/* wN must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[8])){
		mexErrMsgTxt("opsor_elin_2d error: 'wN' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wN.ndims = mxGetNumberOfDimensions(prhs[8]);
	wN.dimElems = mxGetDimensions(prhs[8]);
	wN.data = (float*)mxGetPr( prhs[8] );
	/*
	//----------
	//--- wE ---
	//----------
	*/
	/* wE must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[9])){
		mexErrMsgTxt("opsor_elin_2d error: 'wE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wE.ndims = mxGetNumberOfDimensions(prhs[9]);
	wE.dimElems = mxGetDimensions(prhs[9]);
	wE.data = (float*)mxGetPr( prhs[9] );
	/*
	//----------
	//--- wS ---
	//----------
	*/
	/* wS must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[10])){
		mexErrMsgTxt("opsor_elin_2d error: 'wS' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wS.ndims = mxGetNumberOfDimensions(prhs[10]);
	wS.dimElems = mxGetDimensions(prhs[10]);
	wS.data = (float*)mxGetPr( prhs[10] );
	/*
	//-----------------------------------
	//--- Number of iterations 'iter' ---
	//-----------------------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[11])){
		mexErrMsgTxt("opsor_elin_2d error: 'iter' must be a noncomplex, single-type scalar");
	}
	Mparam.iter = *((float*)mxGetPr( prhs[11] ) );
	/*
	//------------------------------------
	//--- Relaxation parameter 'omega' ---
	//------------------------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[12])){
		mexErrMsgTxt("opsor_elin_2d error: 'omega' must be a noncomplex, single-type scalar");
	}
	Mparam.omega = *((float*)mxGetPr( prhs[12] ) );
	/*
	//-------------------
	//--- Solver type ---
	//-------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[13])){
		mexErrMsgTxt("opsor_elin_2d error: 'solver' must be a noncomplex, single-type scalar");
	}
	solver = (int)*((float*)mxGetPr( prhs[13] ) );
		
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	U_out.ndims = U_in.ndims;
	U_out.dimElems = U_in.dimElems;
	V_out.ndims = V_in.ndims;
	V_out.dimElems = V_in.dimElems;
	
	if( nlhs<2 ){
		mexErrMsgTxt("opsor_elin_2d insufficient number of outputs. Outputs from this function are 'U' and 'V'");
	}
	else{
		if( (plhs[0] = mxCreateNumericArray(U_out.ndims, V_out.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("opsor_elin_2d: error reserving space for output variable 'U_out'");
		if( (U_out.data = (float*)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("opsor_elin_2d: error obtaining pointer to output variable 'U_out'");
	
		if( (plhs[1] = mxCreateNumericArray(V_out.ndims, V_out.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("opsor_elin_2d: error reserving space for output variable 'V_out'");
		if( (V_out.data = (float*)mxGetPr( plhs[1] ))==NULL )
			mexErrMsgTxt("opsor_elin_2d: error obtaining pointer to output variable 'V_out'");
	}
	if( nlhs>=4 )
	{
		RU.ndims = M.ndims;
		RU.dimElems = M.dimElems;
		RV.ndims = M.ndims;
		RV.dimElems = M.dimElems;

		if( (plhs[2] = mxCreateNumericArray(RU.ndims, RU.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("opsor_elin_2d: error reserving space for output variable 'RU'");
		if( (RU.data = (float*)mxGetPr( plhs[2] ))==NULL )
			mexErrMsgTxt("opsor_elin_2d: error obtaining pointer to output variable 'RU'");
	
		if( (plhs[3] = mxCreateNumericArray(RV.ndims, RV.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("opsor_elin_2d: error reserving space for output variable 'RV'");
		if( (RV.data = (float*)mxGetPr( plhs[3] ))==NULL )
			mexErrMsgTxt("opsor_elin_2d: error obtaining pointer to output variable 'RV'");
	}
	
	/*Choose solver*/
	switch (solver)
	{
		case 1:
			solveME = &GS_SOR_elin4_2d;
			break;
		case 2:	
			solveME = &GS_ALR_SOR_elin4_2d;
			break;
		default:
			mexErrMsgTxt("opsor_elin4_2d: no such solver");
	}
		
	/*Call solver*/		
	if( Mparam.iter>0 )
	{
		memcpy( U_out.data, U_in.data, U_in.dimElems[0]*U_in.dimElems[1]*sizeof(float) );
		memcpy( V_out.data, V_in.data, V_in.dimElems[0]*V_in.dimElems[1]*sizeof(float) );
		solveME(&U_out, &V_out, &M, &Cu, &Cv, &Du, &Dv, &wW, &wN, &wE, &wS, Mparam);
	}

	/*Calculate residuals*/
	if( nlhs>=4 )
		Residuals_elin4_2d(&RU, &RV, &U_in, &V_in, &M, &Cu, &Cv, &Du, &Dv, &wW, &wN, &wE, &wS, Mparam);

}
