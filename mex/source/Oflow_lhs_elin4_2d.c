/*
function [AU UV] = oflow_lhs_elin4_2d(U_in, V_in, M, Du, Dv, wW, wN, wE, wS, iter)

Calculates LHS (left hand size of Ax = b) of the optical flow equations using nonlinear diffusion as smoothness term and 
early linearization of the constancy terms. 

A MEX c-program for Matlab.

INPUT
U_in		=		Optical flow, horizontal component.
V_in		=		Optical flow, vertical component.
M		=		Common multiplier
Du		=		Diagonal component of Au from the data term.
Dv		=		Diagonal component of Av from the data term.
wW		=		Diffusion weight (west).
wN		=		Diffusion weight (north).
wE		=		Diffusion weight (east).
wS		=		Diffusion weight (south).

OUTPUT
AU		=		LHS for U term.
AV		=		LHS for V term.

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
	struct matrixM AU = {0,NULL,NULL,0,0,0,0};	
	struct matrixM AV = {0,NULL,NULL,0,0,0,0};
	/*Constant terms*/
	struct matrixM M = {0,NULL,NULL,0,0,0,0};
	struct matrixM Du = {0,NULL,NULL,0,0,0,0};
	struct matrixM Dv = {0,NULL,NULL,0,0,0,0};
	/*Regularization diffusion terms*/
	struct matrixM wW = {0,NULL,NULL,0,0,0,0};
	struct matrixM wN = {0,NULL,NULL,0,0,0,0};
	struct matrixM wE = {0,NULL,NULL,0,0,0,0};
	struct matrixM wS = {0,NULL,NULL,0,0,0,0};

	struct mparam Mparam = {0,0,0,0,0,0,1,1,1};
  
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 9) {
		printf("%d\n",nrhs);
		mexErrMsgTxt("oplhs_lin_2d parameter error: wrong number of input parameters!");
	} 
	/*
	//------------
	//--- U in ---
	//------------
	/*
	/* U_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("oplhs_lin_2d error: 'U_in' must be a noncomplex single-valued matrix.");
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
		mexErrMsgTxt("oplhs_lin_2d error: 'V_in' must be a noncomplex single-valued matrix.");
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
		mexErrMsgTxt("oplhs_lin_2d error: 'M' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	M.ndims = mxGetNumberOfDimensions(prhs[2]);
	M.dimElems = mxGetDimensions(prhs[2]);
	M.data = (float*)mxGetPr( prhs[2] );
	/*
	//------------
	//--- Du in ---
	//------------
	*/
	/* Du must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[3])){
		mexErrMsgTxt("oplhs_lin_2d error: 'Du' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Du.ndims = mxGetNumberOfDimensions(prhs[3]);
	Du.dimElems = mxGetDimensions(prhs[3]);
	Du.data = (float*)mxGetPr( prhs[3] );
	/*
	//------------
	//--- Dv in ---
	//------------
	*/
	/* Dv must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[4])){
		mexErrMsgTxt("oplhs_lin_2d error: 'Dv' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Dv.ndims = mxGetNumberOfDimensions(prhs[4]);
	Dv.dimElems = mxGetDimensions(prhs[4]);
	Dv.data = (float*)mxGetPr( prhs[4] );
	/*
	//----------
	//--- wW ---
	//----------
	*/
	/* wW must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[5])){
		mexErrMsgTxt("oplhs_lin_2d error: 'wW' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wW.ndims = mxGetNumberOfDimensions(prhs[5]);
	wW.dimElems = mxGetDimensions(prhs[5]);
	wW.data = (float*)mxGetPr( prhs[5] );
	/*
	//----------
	//--- wN ---
	//----------
	*/
	/* wN must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[6])){
		mexErrMsgTxt("oplhs_lin_2d error: 'wN' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wN.ndims = mxGetNumberOfDimensions(prhs[6]);
	wN.dimElems = mxGetDimensions(prhs[6]);
	wN.data = (float*)mxGetPr( prhs[6] );
	/*
	//----------
	//--- wE ---
	//----------
	*/
	/* wE must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[7])){
		mexErrMsgTxt("oplhs_lin_2d error: 'wE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wE.ndims = mxGetNumberOfDimensions(prhs[7]);
	wE.dimElems = mxGetDimensions(prhs[7]);
	wE.data = (float*)mxGetPr( prhs[7] );
	/*
	//----------
	//--- wS ---
	//----------
	*/
	/* wS must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[8])){
		mexErrMsgTxt("oplhs_lin_2d error: 'wS' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wS.ndims = mxGetNumberOfDimensions(prhs[8]);
	wS.dimElems = mxGetDimensions(prhs[8]);
	wS.data = (float*)mxGetPr( prhs[8] );

	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	AU.ndims = M.ndims;
	AU.dimElems = M.dimElems;
	AV.ndims = M.ndims;
	AV.dimElems = M.dimElems;
	
	if( nlhs<2 ){
		mexErrMsgTxt("oplhs_lin_2d insufficient number of outputs. Outputs from this function are 'AU' and 'AV'");
	}
	else{
		if( (plhs[0] = mxCreateNumericArray(AU.ndims, AU.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("oplhs_lin_2d: error reserving space for output variable 'AU'");
		if( (AU.data = (float*)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("oplhs_lin_2d: error obtaining pointer to output variable 'AU'");
	
		if( (plhs[1] = mxCreateNumericArray(AV.ndims, AV.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("oplhs_lin_2d: error reserving space for output variable 'AV'");
		if( (AV.data = (float*)mxGetPr( plhs[1] ))==NULL )
			mexErrMsgTxt("oplhs_lin_2d: error obtaining pointer to output variable 'AV'");
	}

	LHS_elin4_2d(&AU, &AV, &U_in, &V_in, &M, &Du, &Dv, &wW, &wN, &wE, &wS, Mparam);

}
