/*
function [dU_out0 dU_out1] = dispsor_nlin_d(
      U_in0, 
      dU_in0, 
      Cu_in0, 
      Du0, 
      wW0, 
      wN0, 
      wE0, 
      wS0,
      U_in1, 
      dU_in1, 
      Cu_in1, 
      Du1, 
      wW1, 
      wN1, 
      wE1, 
      wS1,  
      iter,
      omega)

Calculates disparity using nonlinear diffusion as smoothness term and early linearization of the constancy terms. Imposes symmetry.

A MEX c-program for Matlab.

INPUT
U_in0		=		Disparity.
dU_in0		=		Change in disparity.
Cu_in0		=		Right hand side for U (A*dU = Cu).
Du0		=		Diagonal component of A from the data term.
wW0		=		Diffusion weight (west).
wN0		=		Diffusion weight (north).
wE0		=		Diffusion weight (east).
wS0		=		Diffusion weight (south).
U_in1		=		Disparity
dU_in1		=		Change in disparity.
Cu_in1		=		Right hand side for U (A*dU = Cu).
Du1		=		Diagonal component of A from the data term.
wW1		=		Diffusion weight (west).
wN1		=		Diffusion weight (north).
wE1		=		Diffusion weight (east).
wS1		=		Diffusion weight (south).
iter		=		Number of iterations.
omaga		=

OUTPUT
dU_out0	=		New estimation of dU after 'iter' iterations.
dU_out1	=		New estimation of dU after 'iter' iterations.

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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*Structs*/
	struct matrixM U_in0 = {0,NULL,NULL,0,0,0,0};	
	struct matrixM dU_in0 = {0,NULL,NULL,0,0,0,0};	
	struct matrixM dU_out0 = {0,NULL,NULL,0,0,0,0};	

	struct matrixM U_in1 = {0,NULL,NULL,0,0,0,0};	
	struct matrixM dU_in1 = {0,NULL,NULL,0,0,0,0};	
	struct matrixM dU_out1 = {0,NULL,NULL,0,0,0,0};	
	/*Constant terms*/
	struct matrixM Cu0 = {0,NULL,NULL,0,0,0,0};
	struct matrixM Du0 = {0,NULL,NULL,0,0,0,0};

	struct matrixM Cu1 = {0,NULL,NULL,0,0,0,0};
	struct matrixM Du1 = {0,NULL,NULL,0,0,0,0};
	/*Regularization diffusion terms*/
	struct matrixM wW0 = {0,NULL,NULL,0,0,0,0};
	struct matrixM wN0 = {0,NULL,NULL,0,0,0,0};
	struct matrixM wE0 = {0,NULL,NULL,0,0,0,0};
	struct matrixM wS0 = {0,NULL,NULL,0,0,0,0};

	struct matrixM wW1 = {0,NULL,NULL,0,0,0,0};
	struct matrixM wN1 = {0,NULL,NULL,0,0,0,0};
	struct matrixM wE1 = {0,NULL,NULL,0,0,0,0};
	struct matrixM wS1 = {0,NULL,NULL,0,0,0,0};

	struct mparam Mparam = {0,0,0,0,0,0,1,1,1};

	int solver = 0;
	void (*solveME)(struct matrixM *U_in0, 
			struct matrixM *dU_out0, 
			struct matrixM *Cu0, 
			struct matrixM *Du0, 
			struct matrixM *wW0, 
			struct matrixM *wN0, 
			struct matrixM *wE0, 
			struct matrixM *wS0, 
			struct matrixM *U_in1, 
			struct matrixM *dU_out1, 
			struct matrixM *Cu1, 
			struct matrixM *Du1, 
			struct matrixM *wW1, 
			struct matrixM *wN1, 
			struct matrixM *wE1, 
			struct matrixM *wS1,
			struct mparam Mparam);
  
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 19) {
		mexErrMsgTxt("Disp_sor_llin_sym4_2d parameter error: wrong number of input parameters!");
	} 
	/*
	//------------
	//--- U in ---
	//------------
	/*
	/* U_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'U_in0' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	U_in0.ndims = mxGetNumberOfDimensions(prhs[0]);
	U_in0.dimElems = mxGetDimensions(prhs[0]);
	U_in0.data = (float*)mxGetPr( prhs[0] );
	/*
	//------------
	//--- dU in ---
	//------------
	/*
	/* dU_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'dU_in0' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	dU_in0.ndims = mxGetNumberOfDimensions(prhs[1]);
	dU_in0.dimElems = mxGetDimensions(prhs[1]);
	dU_in0.data = (float*)mxGetPr( prhs[1] );
	/*
	//------------
	//--- Cu in ---
	//------------
	/*
	/* Cu must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[2])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'Cu0' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Cu0.ndims = mxGetNumberOfDimensions(prhs[2]);
	Cu0.dimElems = mxGetDimensions(prhs[2]);
	Cu0.data = (float*)mxGetPr( prhs[2] );
	/*
	//------------
	//--- Du in ---
	//------------
	/*
	/* Du must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[3])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'Du0' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Du0.ndims = mxGetNumberOfDimensions(prhs[3]);
	Du0.dimElems = mxGetDimensions(prhs[3]);
	Du0.data = (float*)mxGetPr( prhs[3] );
	/*
	//----------
	//--- wW ---
	//----------
	/*
	/* wW must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[4])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'wW0' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wW0.ndims = mxGetNumberOfDimensions(prhs[4]);
	wW0.dimElems = mxGetDimensions(prhs[4]);
	wW0.data = (float*)mxGetPr( prhs[4] );
	/*
	//----------
	//--- wN ---
	//----------
	/*
	/* wN must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[5])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'wN0' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wN0.ndims = mxGetNumberOfDimensions(prhs[5]);
	wN0.dimElems = mxGetDimensions(prhs[5]);
	wN0.data = (float*)mxGetPr( prhs[5] );
	/*
	//----------
	//--- wE ---
	//----------
	/*
	/* wE must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[6])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'wE0' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wE0.ndims = mxGetNumberOfDimensions(prhs[6]);
	wE0.dimElems = mxGetDimensions(prhs[6]);
	wE0.data = (float*)mxGetPr( prhs[6] );
	/*
	//----------
	//--- wS ---
	//----------
	/*
	/* wS must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[7])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'wS0' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wS0.ndims = mxGetNumberOfDimensions(prhs[7]);
	wS0.dimElems = mxGetDimensions(prhs[7]);
	wS0.data = (float*)mxGetPr( prhs[7] );
	/*
	//------------
	//--- U in ---
	//------------
	/*
	/* U_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[8])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'U_in1' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	U_in1.ndims = mxGetNumberOfDimensions(prhs[8]);
	U_in1.dimElems = mxGetDimensions(prhs[8]);
	U_in1.data = (float*)mxGetPr( prhs[8] );
	/*
	//------------
	//--- dU in ---
	//------------
	/*
	/* dU_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[9])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'dU_in1' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	dU_in1.ndims = mxGetNumberOfDimensions(prhs[9]);
	dU_in1.dimElems = mxGetDimensions(prhs[9]);
	dU_in1.data = (float*)mxGetPr( prhs[9] );
	/*
	//------------
	//--- Cu in ---
	//------------
	/*
	/* Cu must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[10])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'Cu1' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Cu1.ndims = mxGetNumberOfDimensions(prhs[10]);
	Cu1.dimElems = mxGetDimensions(prhs[10]);
	Cu1.data = (float*)mxGetPr( prhs[10] );
	/*
	//------------
	//--- Du in ---
	//------------
	/*
	/* Du must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[11])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'Du1' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Du1.ndims = mxGetNumberOfDimensions(prhs[11]);
	Du1.dimElems = mxGetDimensions(prhs[11]);
	Du1.data = (float*)mxGetPr( prhs[11] );
	/*
	//----------
	//--- wW ---
	//----------
	/*
	/* wW must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[12])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'wW1' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wW1.ndims = mxGetNumberOfDimensions(prhs[12]);
	wW1.dimElems = mxGetDimensions(prhs[12]);
	wW1.data = (float*)mxGetPr( prhs[12] );
	/*
	//----------
	//--- wN ---
	//----------
	/*
	/* wN must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[13])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'wN1' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wN1.ndims = mxGetNumberOfDimensions(prhs[13]);
	wN1.dimElems = mxGetDimensions(prhs[13]);
	wN1.data = (float*)mxGetPr( prhs[13] );
	/*
	//----------
	//--- wE ---
	//----------
	/*
	/* wE must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[14])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'wE1' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wE1.ndims = mxGetNumberOfDimensions(prhs[14]);
	wE1.dimElems = mxGetDimensions(prhs[14]);
	wE1.data = (float*)mxGetPr( prhs[14] );
	/*
	//----------
	//--- wS ---
	//----------
	/*
	/* wS must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[15])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'wS1' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wS1.ndims = mxGetNumberOfDimensions(prhs[15]);
	wS1.dimElems = mxGetDimensions(prhs[15]);
	wS1.data = (float*)mxGetPr( prhs[15] );
	/*
	//-----------------------------------
	//--- Number of iterations 'iter' ---
	//-----------------------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[16])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'iter' must be a noncomplex, single-type scalar");
	}
	Mparam.iter = *((float*)mxGetPr( prhs[16] ) );
	/*
	//------------------------------------
	//--- Relaxation parameter 'omega' ---
	//------------------------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[17])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'omega' must be a noncomplex, single-type scalar");
	}
	Mparam.omega = *((float*)mxGetPr( prhs[17] ) );
	/*
	//-------------------
	//--- Solver type ---
	//-------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[18])){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d error: 'solver' must be a noncomplex, single-type scalar");
	}
	solver = (int)*((float*)mxGetPr( prhs[18] ) );

	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	dU_out0.ndims = dU_in0.ndims;
	dU_out0.dimElems = dU_in0.dimElems;
	
	dU_out1.ndims = dU_in1.ndims;
	dU_out1.dimElems = dU_in1.dimElems;
		
	if( nlhs<2 ){
		mexErrMsgTxt("Disp_sor_llin_sym4_2d insufficient number of outputs. Outputs from this function are 'dU0' and 'dU1'");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(dU_out0.ndims, dU_out0.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("Disp_sor_llin_sym4_2d: error reserving space for output variable 'dU_out0'");
		if( (dU_out0.data = (float*)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("Disp_sor_llin_sym4_2d: error obtaining pointer to output variable 'dU_out0'");
	
		if( (plhs[1] = mxCreateNumericArray(dU_out1.ndims, dU_out1.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("Disp_sor_llin_sym4_2d: error reserving space for output variable 'dU_out1'");
		if( (dU_out1.data = (float*)mxGetPr( plhs[1] ))==NULL )
			mexErrMsgTxt("Disp_sor_llin_sym4_2d: error obtaining pointer to output variable 'dU_out1'");
	
	}

	/*Choose solver*/
	switch (solver)
	{
		case 1:
			solveME = &GS_SOR_llinsym4_2d;
			break;
		case 2:	
			solveME = &GS_ALR_SOR_llinsym4_2d;
			break;
		default:
			mexErrMsgTxt("dispsor_llinsym4_2d: no such solver");
	}

	memcpy( dU_out0.data, dU_in0.data, dU_in0.dimElems[0]*dU_in0.dimElems[1]*sizeof(float) );
	memcpy( dU_out1.data, dU_in1.data, dU_in1.dimElems[0]*dU_in1.dimElems[1]*sizeof(float) );

	solveME(	&U_in0, 
			&dU_out0, 
			&Cu0, 
			&Du0, 
			&wW0, 
			&wN0, 
			&wE0, 
			&wS0, 
			&U_in1, 
			&dU_out1, 
			&Cu1, 
			&Du1, 
			&wW1, 
			&wN1, 
			&wE1, 
			&wS1,
			Mparam
		);

}
