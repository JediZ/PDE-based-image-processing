/*
function [wW wN wE wS] = DdiffWeights(D, eps)

Calculate diffusion coefficients for disparity 

A MEX c-program for Matlab.

INPUT
D	=	Data for which diffusion coefficients are calculated.
eps	=	Small coefficient to avoid division by zero.

OUTPUT
wW	=	West diffusion weights.
wN	=	North diffusion weights.
wE	=	East diffusion weights.
wS	=	South diffusion weights.

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

#include "imageDiffusionWeights.h"

/*
//------------------
//--- MAIN (MEX) ---
//------------------
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
	/*Structs*/
	struct matrixM D = {0,NULL,NULL,0,0,0,0};	
	
	/*Regularization diffusion terms*/
	struct matrixM wW = {0,NULL,NULL,0,0,0,0};
	struct matrixM wN = {0,NULL,NULL,0,0,0,0};
	struct matrixM wE = {0,NULL,NULL,0,0,0,0};
	struct matrixM wS = {0,NULL,NULL,0,0,0,0};

	float eps;

	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 2) {
		mexErrMsgTxt("diffusion6_2d parameter error: wrong number of input parameters!");
	} 
	/*
	//------------
	//--- D in ---
	//------------
	*/
	/* D_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("diffusion6_2d error: 'U_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	D.ndims = mxGetNumberOfDimensions(prhs[0]);
	D.dimElems = mxGetDimensions(prhs[0]);
	D.data = (float*)mxGetPr( prhs[0] );
	/*
	//-----------
	//--- eps ---
	//-----------
	*/
	/* iter must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("diffusion6_2d error: 'eps' must be a noncomplex, single-type scalar");
	}
	eps = *((float*)mxGetPr( prhs[1] ) );
		
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	wW.ndims = D.ndims;
	wW.dimElems = D.dimElems;
	
	wN.ndims = D.ndims;
	wN.dimElems = D.dimElems;

	wE.ndims = D.ndims;
	wE.dimElems = D.dimElems;
	
	wS.ndims = D.ndims;
	wS.dimElems = D.dimElems;

	if( nlhs<4 ){
		mexErrMsgTxt("diffusion6_2d error insufficient number of outputs. Outputs from this function are 'wW', 'wN', 'wE' and 'wS'");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(wW.ndims, wW.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("diffusion6_2d error: error reserving space for output variable 'wW'");
		if( (wW.data = (float *)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("diffusion6_2d error: error obtaining pointer to output variable 'wW'");

		if( (plhs[1] = mxCreateNumericArray(wN.ndims, wN.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("diffusion6_2d error: error reserving space for output variable 'wN'");
		if( (wN.data = (float *)mxGetPr( plhs[1] ))==NULL )
			mexErrMsgTxt("diffusion6_2d error: error obtaining pointer to output variable 'wN'");
	
		if( (plhs[2] = mxCreateNumericArray(wE.ndims, wE.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("diffusion6_2d error: error reserving space for output variable 'wE'");
		if( (wE.data = (float *)mxGetPr( plhs[2] ))==NULL )
			mexErrMsgTxt("diffusion6_2d error: error obtaining pointer to output variable 'wE'");

		if( (plhs[3] = mxCreateNumericArray(wS.ndims, wS.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("diffusion6_2d error: error reserving space for output variable 'wS'");
		if( (wS.data = (float *)mxGetPr( plhs[3] ))==NULL )
			mexErrMsgTxt("diffusion6_2d error: error obtaining pointer to output variable 'wS'");
	}

	diffWeights6_2D_c( &wW, &wN, &wE, &wS, &D, eps );

}
