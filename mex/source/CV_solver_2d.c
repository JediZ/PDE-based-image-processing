/*
function [PHI_out varargout] = CV_solver_2d(PHIin, D, DH, GradNorm, tau, nu)

AOS solver for Chan&Vese active region segmentation model.

A MEX c-program for Matlab.

INPUT
PHIin		=		Level set function PHI.
D		=		Data term.
DH		=		Blurred version of Heaviside derivative.
GradNorm	=		Euclidean gradient norm of Oin.
tau		=		Time step.
nu		=		Diffusion term influence coefficient.

OUTPUT

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

#include "levelsetSolvers.h"
#include <omp.h>

/*
//------------------
//--- MAIN (MEX) ---
//------------------
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
	/*Structs*/
	struct matrixM PHI_in = {0,NULL,NULL,0,0,0,0};
	struct matrixM PHI_out = {0,NULL,NULL,0,0,0,0};
	struct matrixM D_in = {0,NULL,NULL,0,0,0,0};	
	struct matrixM DH_in = {0,NULL,NULL,0,0,0,0};
	struct matrixM GradNorm_in = {0,NULL,NULL,0,0,0,0};

	float tau = 0.0f;
	float nu = 0.0f;
	
	void (*solveME)(	struct matrixM *PHI_out,
				struct matrixM *PHI_in,
				struct matrixM *D_in,
				struct matrixM *GradNorm_in,
				struct matrixM *Diff_in,
				float tau,
				float nu );
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 6) {
		mexErrMsgTxt("cv_solver_2D parameter error: wrong number of input parameters!");
	} 
	/*
	//------------
	//--- O in ---
	//------------
	*/
	/* PHI_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("cv_solver_2D error: 'PHI_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	PHI_in.ndims = mxGetNumberOfDimensions(prhs[0]);
	PHI_in.dimElems = mxGetDimensions(prhs[0]);
	PHI_in.data = (float*)mxGetPr( prhs[0] );
	/*
	//------------
	//--- D in ---
	//------------
	*/
	/* D_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("cv_solver_2D error: 'D_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	D_in.ndims = mxGetNumberOfDimensions(prhs[1]);
	D_in.dimElems = mxGetDimensions(prhs[1]);
	D_in.data = (float*)mxGetPr( prhs[1] );
	/*
	//------------
	//--- DH in ---
	//------------
	*/
	/* DH_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[2])){
		mexErrMsgTxt("cv_solver_2D error: 'DH_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	DH_in.ndims = mxGetNumberOfDimensions(prhs[2]);
	DH_in.dimElems = mxGetDimensions(prhs[2]);
	DH_in.data = (float*)mxGetPr(prhs[2]);
	/*
	//------------
	//--- GradNorm in ---
	//------------
	*/
	/* GradNorm_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[3])){
		mexErrMsgTxt("cv_solver_2D error: 'GradNorm_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	GradNorm_in.ndims = mxGetNumberOfDimensions(prhs[3]);
	GradNorm_in.dimElems = mxGetDimensions(prhs[3]);
	GradNorm_in.data = (float*)mxGetPr( prhs[3] );
	/*
	//-----------------------------------
	//--- Time step 'tau ---
	//-----------------------------------
	*/
	/* tau must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[4])){
		mexErrMsgTxt("cv_solver_2D error: 'tau' must be a noncomplex, single-type scalar");
	}
	tau = *((float*)mxGetPr( prhs[4] ) );
	/*
	//-----------------------------------
	//--- Diffusion strength 'nu' ---
	//-----------------------------------
	*/
	/* nu must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[5])){
		mexErrMsgTxt("cv_solver_2D error: 'nu' must be a noncomplex, single-type scalar");
	}
	nu = *((float*)mxGetPr( prhs[5] ) );
	
	/* Check if maximum buffer size is big enough */
	switch( PHI_in.ndims )
	{
	  case 1:
		if (PHI_in.dimElems[0]>MAX_BUF_SIZE)
		    mexErrMsgTxt("cv_solver_2D error: MAX_BUF_SIZE too small...increase the size in the header file!");
		break;
	  default :
		if ( (PHI_in.dimElems[0]>MAX_BUF_SIZE) || (PHI_in.dimElems[1]>MAX_BUF_SIZE) )
		  mexErrMsgTxt("cv_solver_2D error: MAX_BUF_SIZE too small...increase the size in the header file!");
		break;
	}
	
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	PHI_out.ndims = PHI_in.ndims;
	PHI_out.dimElems = PHI_in.dimElems;

	if( nlhs<1 ){
		mexErrMsgTxt("cv_solver_2D error insufficient number of outputs. Outputs from this function is 'PHI_out'");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(PHI_out.ndims,PHI_out.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("cv_solver_2D error: error reserving space for output variable 'PHI_out'");

		if( (PHI_out.data = (float *)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("cv_solver_2D error: error obtaining pointer to output variable 'PHI_out'");
	}
	
	solveME = &CV_AOSOMP_4_2d;	/*OpenMP version...needs to be compiled "against" OPENMP, -fopenmp in C- and LDFLAGS*/
	/*solveME = &CV_AOS_4_2d;	/*Basic version...-fopenmp not needed*/

	solveME(	&PHI_out,
			&PHI_in,
			&D_in,
			&DH_in,
			&GradNorm_in,
			tau, 
			nu);
}
