/*
function [Idxt Idyt Idxx Idyy Idxy] = SndDerivatives( It0, It1 )

Approximates Second order image derivatives using adaptation of Simoncelli-filters.

A MEX c-program for Matlab.

INPUT
It0		=		Image at time t=0
It1		=		Image at time t=1

OUTPUT
Idxt		=		Temporal derivative
Idyt		=		Temporal derivative
Idxx		=		Horizontal spatial derivative
Idyy		=		Vertical spatial derivative
Idxy		=		Spatial cross-derivative

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

#include"imageDerivatives.h"

/*
//------------------
//--- MAIN (MEX) ---
//------------------
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
	/*Structs*/
	struct matrixM It0 = {0,NULL,NULL,0,0,0,0};
	struct matrixM It1 = {0,NULL,NULL,0,0,0,0};
	struct matrixM Idxt = {0,NULL,NULL,0,0,0,0};
	struct matrixM Idyt = {0,NULL,NULL,0,0,0,0};
	struct matrixM Idxx = {0,NULL,NULL,0,0,0,0};
	struct matrixM Idyy = {0,NULL,NULL,0,0,0,0};
	struct matrixM Idxy = {0,NULL,NULL,0,0,0,0};
	struct matrixM temp1 = {0,NULL,NULL,0,0,0,0};
	struct matrixM temp2 = {0,NULL,NULL,0,0,0,0};

	/*Simoncelli filters*/
	float Smoother[5] = 		{0.037659f,	0.249724f,	0.439911f,	0.249724f,	0.037659f};
	float fstSpatialDerivator[5] = 	{-0.104550f,	-0.292315f,	0.0f,		0.292315f,	0.104550f};
	float sndSpatialDerivator[5] = 	{0.232905f,	0.002668f,	-0.471147f,	0.002668f,	0.232905f};
	float TemporalDerivator[2] =	{0.50f, -0.50f};
	
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 2) {
		mexErrMsgTxt("sndDerivatives: wrong number of input parameters!");
	} 
	/*
	//--------------
	//--- It0 in ---
	//--------------
	*/
	/* It0 must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("sndDerivatives: 'It0' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	It0.ndims = mxGetNumberOfDimensions(prhs[0]);
	It0.dimElems = mxGetDimensions(prhs[0]);
	It0.data = (float*)mxGetPr( prhs[0] );
	/*
	//--------------
	//--- It1 in ---
	//--------------
	*/
	/* It1 must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("sndDerivatives: 'It1' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	It1.ndims = mxGetNumberOfDimensions(prhs[1]);
	It1.dimElems = mxGetDimensions(prhs[1]);
	It1.data = (float*)mxGetPr( prhs[1] );
	
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	Idxt.ndims = It0.ndims;
	Idxt.dimElems = It0.dimElems;
	Idyt.ndims = It0.ndims;
	Idyt.dimElems = It0.dimElems;
	Idxx.ndims = It0.ndims;
	Idxx.dimElems = It0.dimElems;
	Idyy.ndims = It0.ndims;
	Idyy.dimElems = It0.dimElems;
	Idxy.ndims = It0.ndims;
	Idxy.dimElems = It0.dimElems;
	temp1.ndims = It0.ndims;
	temp1.dimElems = It0.dimElems;
	temp2.ndims = It0.ndims;
	temp2.dimElems = It0.dimElems;
	
		
	if( nlhs<5 ){
		mexErrMsgTxt("sndDerivatives: insufficient number of outputs. Outputs from this function are 'Idxt', 'Idyt' 'Idxx', 'Idyy' and 'Idxy'.");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(Idxt.ndims, Idxt.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("sndDerivatives: error reserving space for output variable 'Idxt'");
		if( (Idxt.data = (float *)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("sndDerivatives: error obtaining pointer to output variable 'Idxt'");

		if( (plhs[1] = mxCreateNumericArray(Idyt.ndims, Idyt.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("sndDerivatives: error reserving space for output variable 'Idyt'");
		if( (Idyt.data = (float *)mxGetPr( plhs[1] ))==NULL )
			mexErrMsgTxt("sndDerivatives: error obtaining pointer to output variable 'Idyt'");

		if( (plhs[2] = mxCreateNumericArray(Idxx.ndims, Idxx.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("sndDerivatives: error reserving space for output variable 'Idxx'");
		if( (Idxx.data = (float *)mxGetPr( plhs[2] ))==NULL )
			mexErrMsgTxt("sndDerivatives: error obtaining pointer to output variable 'Idxx'");

		if( (plhs[3] = mxCreateNumericArray(Idyy.ndims, Idyy.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("sndDerivatives: error reserving space for output variable 'Idyy'");
		if( (Idyy.data = (float *)mxGetPr( plhs[3] ))==NULL )
			mexErrMsgTxt("sndDerivatives: error obtaining pointer to output variable 'Idyy'");

		if( (plhs[4] = mxCreateNumericArray(Idxy.ndims, Idxy.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("sndDerivatives: error reserving space for output variable 'Idxy'");
		if( (Idxy.data = (float *)mxGetPr( plhs[4] ))==NULL )
			mexErrMsgTxt("sndDerivatives: error obtaining pointer to output variable 'Idxy'");
	}

	if(temp1.ndims==2)
		temp1.data = (float*)mxCalloc( temp1.dimElems[0]*temp1.dimElems[1], sizeof(float) );
	else
		temp1.data = (float*)mxCalloc( temp1.dimElems[0]*temp1.dimElems[1]*temp1.dimElems[2], sizeof(float) );

	if(temp2.ndims==2)
		temp2.data = (float*)mxCalloc( temp2.dimElems[0]*temp2.dimElems[1], sizeof(float) );
	else
		temp2.data = (float*)mxCalloc( temp2.dimElems[0]*temp2.dimElems[1]*temp2.dimElems[2], sizeof(float) );

	
	sndSimoncelli_c(	&Idxt, &Idyt, &Idxx, &Idyy, &Idxy, &It0, &It1, &temp1, &temp2, 5, Smoother, fstSpatialDerivator,
				sndSpatialDerivator, TemporalDerivator );
	
	mxFree( temp1.data );
	mxFree( temp2.data );

}
