/*
function [dU_out dV_out varargout] = opsor_llin_d(
      U_in, V_in, dU_in, dV_in, M, Cu_in, Cv_in, Du, Dv, wW, wNW, wN, wNE, wE, wSE, wS, wSW, iter, solver)

Calculates optical flow using nonlinear diffusion as smoothness term and late linearization of the constancy terms. 

A MEX c-program for Matlab.

INPUT
U_in		=		Optical flow, horizontal component.
V_in		=		Optical flow, vertical component.
dU_in		=		Change in optical flow, horizontal component.
dV_in		=		Change in optical flow, vertical component.
M		=		Common multiplier
Cu_in		=		Right hand side for U (Au*dU = Cu - M*v).
Cv_in		=		Right hand size for V (Av*dV = Cv - M*u).
Du		=		Diagonal component of Au from the data term.
Dv		=		Diagonal component of Av from the data term.
wW		=		Diffusion weight (west).
wNW		=		Diffusion weight (north-west).
wN		=		Diffusion weight (north).
wNE		=		Diffusion weight (north-east).
wE		=		Diffusion weight (east).
wSE		=		Diffusion weight (south-east).
wS		=		Diffusion weight (south).
wSW		=		Diffusion weight (south-west).
iter		=		Number of iterations.
solver		=		Type of solver used:
			      		1 = Normal point-wise Gauss-Seidel relaxation.
			      		2 = Alternating line relaxation (ALR) - Block Gauss-Seidel.

OUTPUT
dU_out		=		New estimation of dU after 'iter' iterations.
dV_out		=		New estimation of dV after 'iter' iterations.
varargout{1}	=		Residual for U. Au*dU = Cu - M*v => r = Cu - M*v - Au*dU
varargout{2}	=		Resitual for V. Av*dV = Cv - M*u => r = Cv - M*u - Av*dV


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
	struct matrixM dU_in = {0,NULL,NULL,0,0,0,0};	
	struct matrixM dV_in = {0,NULL,NULL,0,0,0,0};	
	struct matrixM dU_out = {0,NULL,NULL,0,0,0,0};	
	struct matrixM dV_out = {0,NULL,NULL,0,0,0,0};	
	/*Constant terms*/
	struct matrixM M = {0,NULL,NULL,0,0,0,0};
	struct matrixM Cu = {0,NULL,NULL,0,0,0,0};
	struct matrixM Cv = {0,NULL,NULL,0,0,0,0};
	struct matrixM Du = {0,NULL,NULL,0,0,0,0};
	struct matrixM Dv = {0,NULL,NULL,0,0,0,0};
	/*Regularization diffusion terms*/
	struct matrixM wW = {0,NULL,NULL,0,0,0,0};
	struct matrixM wNW = {0,NULL,NULL,0,0,0,0};
	struct matrixM wN = {0,NULL,NULL,0,0,0,0};
	struct matrixM wNE = {0,NULL,NULL,0,0,0,0};
	struct matrixM wE = {0,NULL,NULL,0,0,0,0};
	struct matrixM wSE = {0,NULL,NULL,0,0,0,0};
	struct matrixM wS = {0,NULL,NULL,0,0,0,0};
	struct matrixM wSW = {0,NULL,NULL,0,0,0,0};
	/*Residuals*/
	struct matrixM RU = {0,NULL,NULL,0,0,0,0};
	struct matrixM RV = {0,NULL,NULL,0,0,0,0};

	struct mparam Mparam = {0,0,0,0,0,0,1,1,1};

	int solver = 0;
	void (*solveME)(struct matrixM *U, 
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
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 20) {
		mexErrMsgTxt("Oflow_sor_llin8_2d parameter error: wrong number of input parameters!");
	} 
	/*
	//------------
	//--- U in ---
	//------------
	*/
	/* U_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'U_in' must be a noncomplex single-valued matrix.");
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
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'V_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	V_in.ndims = mxGetNumberOfDimensions(prhs[1]);
	V_in.dimElems = mxGetDimensions(prhs[1]);
	V_in.data = (float*)mxGetPr( prhs[1] );
	/*
	//------------
	//--- dU in ---
	//------------
	*/
	/* dU_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[2])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'dU_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	dU_in.ndims = mxGetNumberOfDimensions(prhs[2]);
	dU_in.dimElems = mxGetDimensions(prhs[2]);
	dU_in.data = (float*)mxGetPr( prhs[2] );
	/*
	//------------
	//--- dV in ---
	//------------
	*/
	/* v_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[3])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'dV_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	dV_in.ndims = mxGetNumberOfDimensions(prhs[3]);
	dV_in.dimElems = mxGetDimensions(prhs[3]);
	dV_in.data = (float*)mxGetPr( prhs[3] );
	/*
	//------------
	//--- M in ---
	//------------
	*/
	/* M must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[4])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'M' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	M.ndims = mxGetNumberOfDimensions(prhs[4]);
	M.dimElems = mxGetDimensions(prhs[4]);
	M.data = (float*)mxGetPr( prhs[4] );
	/*
	//------------
	//--- Cu in ---
	//------------
	*/
	/* Cu must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[5])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'Cu' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Cu.ndims = mxGetNumberOfDimensions(prhs[5]);
	Cu.dimElems = mxGetDimensions(prhs[5]);
	Cu.data = (float*)mxGetPr( prhs[5] );
	/*
	//------------
	//--- Cv in ---
	//------------
	*/
	/* Cv must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[6])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'Cv' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Cv.ndims = mxGetNumberOfDimensions(prhs[6]);
	Cv.dimElems = mxGetDimensions(prhs[6]);
	Cv.data = (float*)mxGetPr( prhs[6] );
	/*
	//------------
	//--- Du in ---
	//------------
	*/
	/* Du must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[7])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'Du' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Du.ndims = mxGetNumberOfDimensions(prhs[7]);
	Du.dimElems = mxGetDimensions(prhs[7]);
	Du.data = (float*)mxGetPr( prhs[7] );
	/*
	//------------
	//--- Dv in ---
	//------------
	*/
	/* Dv must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[8])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'Dv' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Dv.ndims = mxGetNumberOfDimensions(prhs[8]);
	Dv.dimElems = mxGetDimensions(prhs[8]);
	Dv.data = (float*)mxGetPr( prhs[8] );
	/*
	//----------
	//--- wW ---
	//----------
	*/
	/* wW must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[9])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'wW' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wW.ndims = mxGetNumberOfDimensions(prhs[9]);
	wW.dimElems = mxGetDimensions(prhs[9]);
	wW.data = (float*)mxGetPr( prhs[9] );
	/*
	//----------
	//--- wNW ---
	//----------
	*/
	/* wN must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[10])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'wNW' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wNW.ndims = mxGetNumberOfDimensions(prhs[10]);
	wNW.dimElems = mxGetDimensions(prhs[10]);
	wNW.data = (float*)mxGetPr( prhs[10] );
	/*
	//----------
	//--- wN ---
	//----------
	*/
	/* wN must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[11])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'wN' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wN.ndims = mxGetNumberOfDimensions(prhs[11]);
	wN.dimElems = mxGetDimensions(prhs[11]);
	wN.data = (float*)mxGetPr( prhs[11] );
	/*
	//----------
	//--- wNE ---
	//----------
	*/
	/* wNE must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[12])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'wNE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wNE.ndims = mxGetNumberOfDimensions(prhs[12]);
	wNE.dimElems = mxGetDimensions(prhs[12]);
	wNE.data = (float*)mxGetPr( prhs[12] );
	/*
	//----------
	//--- wE ---
	//----------
	*/
	/* wE must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[13])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'wE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wE.ndims = mxGetNumberOfDimensions(prhs[13]);
	wE.dimElems = mxGetDimensions(prhs[13]);
	wE.data = (float*)mxGetPr( prhs[13] );
	/*
	//----------
	//--- wSE ---
	//----------
	*/
	/* wSE must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[14])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'wSE' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wSE.ndims = mxGetNumberOfDimensions(prhs[14]);
	wSE.dimElems = mxGetDimensions(prhs[14]);
	wSE.data = (float*)mxGetPr( prhs[14] );
	/*
	//----------
	//--- wS ---
	//----------
	*/
	/* wS must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[15])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'wS' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wS.ndims = mxGetNumberOfDimensions(prhs[15]);
	wS.dimElems = mxGetDimensions(prhs[15]);
	wS.data = (float*)mxGetPr( prhs[15] );
	/*
	//----------
	//--- wSW ---
	//----------
	*/
	/* wSW must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[16])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'wSW' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	wSW.ndims = mxGetNumberOfDimensions(prhs[16]);
	wSW.dimElems = mxGetDimensions(prhs[16]);
	wSW.data = (float*)mxGetPr( prhs[16] );
	/*
	//-----------------------------------
	//--- Number of iterations 'iter' ---
	//-----------------------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[17])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'iter' must be a noncomplex, single-type scalar");
	}
	Mparam.iter = *((float*)mxGetPr( prhs[17] ) );
	/*
	//------------------------------------
	//--- Relaxation parameter 'omega' ---
	//------------------------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[18])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'omega' must be a noncomplex, single-type scalar");
	}
	Mparam.omega = *((float*)mxGetPr( prhs[18] ) );
	/*
	//-------------------
	//--- Solver type ---
	//-------------------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[19])){
		mexErrMsgTxt("Oflow_sor_llin8_2d error: 'solver' must be a noncomplex, single-type scalar");
	}
	solver = (int)*((float*)mxGetPr( prhs[19] ) );
	
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	dU_out.ndims = dU_in.ndims;
	dU_out.dimElems = dU_in.dimElems;
	dV_out.ndims = dV_in.ndims;
	dV_out.dimElems = dV_in.dimElems;
	
	if( nlhs<2 ){
		mexErrMsgTxt("Oflow_sor_llin8_2d, insufficient number of outputs. Outputs from this function are 'dU' and 'dV'");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(dU_out.ndims, dU_out.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("Oflow_sor_llin8_2d: error reserving space for output variable 'dU_out'");
		if( (dU_out.data = (float*)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("Oflow_sor_llin8_2d: error obtaining pointer to output variable 'dU_out'");
	
		if( (plhs[1] = mxCreateNumericArray(dV_out.ndims, dV_out.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("Oflow_sor_llin8_2d: error reserving space for output variable 'dV_out'");
		if( (dV_out.data = (float*)mxGetPr( plhs[1] ))==NULL )
			mexErrMsgTxt("Oflow_sor_llin8_2d: error obtaining pointer to output variable 'dV_out'");
	}

	if( nlhs>=4 )
	{
		RU.ndims = M.ndims;
		RU.dimElems = M.dimElems;
		RV.ndims = M.ndims;
		RV.dimElems = M.dimElems;

		if( (plhs[2] = mxCreateNumericArray(RU.ndims, RU.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("Oflow_sor_llin8_2d: error reserving space for output variable 'RU'");
		if( (RU.data = (float*)mxGetPr( plhs[2] ))==NULL )
			mexErrMsgTxt("Oflow_sor_llin8_2d: error obtaining pointer to output variable 'RU'");
	
		if( (plhs[3] = mxCreateNumericArray(RV.ndims, RV.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("Oflow_sor_llin8_2d: error reserving space for output variable 'RV'");
		if( (RV.data = (float*)mxGetPr( plhs[3] ))==NULL )
			mexErrMsgTxt("Oflow_sor_llin8_2d: error obtaining pointer to output variable 'RV'");

	}

	/*Choose solver*/
	switch (solver)
	{
		case 1:
			solveME = &GS_SOR_llin8_2d;
			break;
		case 2:	
			solveME = &GS_ALR_SOR_llin8_2d;
			break;
		default:
			mexErrMsgTxt("Oflow_sor_llin8_2d: no such solver");
	}
	
	/*Call solver*/
	if( Mparam.iter>0 )
	{
		memcpy( dU_out.data, dU_in.data, dU_in.dimElems[0]*dU_in.dimElems[1]*sizeof(float) );
		memcpy( dV_out.data, dV_in.data, dV_in.dimElems[0]*dV_in.dimElems[1]*sizeof(float) );
		solveME(	&U_in,
				&V_in, 
				&dU_out,
				&dV_out,
				&M,
				&Cu,
				&Cv,
				&Du,
				&Dv,
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

	/*Calculate residuals*/
	/*
	if( nlhs>=4 )
		Residuals_llin8_2d(	&RU,
					&RV,
					&U_in,
					&V_in,
					&dU_in,
					&dV_in,
					&M,
					&Cu,
					&Cv,
					&Du,
					&Dv,
					&wW,
					&wNW,
					&wN,
					&wNE,
					&wE,
					&wSE,
					&wS,
					&wSW,
					Mparam);
	*/
}
