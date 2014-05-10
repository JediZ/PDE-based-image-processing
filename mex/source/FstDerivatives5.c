/*
function [Idt Idx Idy] = FstDerivatives5( It0, It1 )

Approximates first order image derivatives using adaptation of Simoncelli-filters.

A MEX c-program for Matlab.

INPUT
It0		=		Image at time t=0
It1		=		Image at time t=1

OUTPUT
Idt		=		Temporal derivative
Idx		=		Horizontal spatial derivative
Idy		=		Vertical spatial derivative

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
	struct matrixM Idt = {0,NULL,NULL,0,0,0,0};
	struct matrixM Idx = {0,NULL,NULL,0,0,0,0};
	struct matrixM Idy = {0,NULL,NULL,0,0,0,0};
	struct matrixM temp = {0,NULL,NULL,0,0,0,0};

	/*Simoncelli filters*/
	float Smoother[5] = 		{0.037659f,	0.249724f,	0.439911f,	0.249724f,	0.037659f};
	float SpatialDerivator[5] = 	{-0.104550f,	-0.292315f,	0.0f,		0.292315f,	0.104550f};
	float TemporalDerivator[5] = 	{0.50f,	-0.50f};

	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 2) {
		mexErrMsgTxt("fstDerivatives: wrong number of input parameters!");
	} 
	/*
	//--------------
	//--- It0 in ---
	//--------------
	*/
	/* It0 must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("fstDerivatives: 'It0' must be a noncomplex single-valued matrix.");
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
		mexErrMsgTxt("fstDerivatives: 'It1' must be a noncomplex single-valued matrix.");
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
	Idt.ndims = It0.ndims;
	Idt.dimElems = It0.dimElems;
	Idx.ndims = It0.ndims;
	Idx.dimElems = It0.dimElems;
	Idy.ndims = It0.ndims;
	Idy.dimElems = It0.dimElems;
	temp.ndims = It0.ndims;
	temp.dimElems = It0.dimElems;
	
		
	if( nlhs<3 ){
		mexErrMsgTxt("fstDerivatives: insufficient number of outputs...outputs from this function are 'Idt', 'Idx' and 'Idy'.");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(Idt.ndims, Idt.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("fstDerivatives: error reserving space for output variable 'Idt'");
		if( (Idt.data = (float *)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("fstDerivatives: error obtaining pointer to output variable 'Idt'");

		if( (plhs[1] = mxCreateNumericArray(Idx.ndims, Idx.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("fstDerivatives: error reserving space for output variable 'Idx'");
		if( (Idx.data = (float *)mxGetPr( plhs[1] ))==NULL )
			mexErrMsgTxt("fstDerivatives: error obtaining pointer to output variable 'Idx'");

		if( (plhs[2] = mxCreateNumericArray(Idy.ndims, Idy.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("fstDerivatives: error reserving space for output variable 'Idy'");
		if( (Idy.data = (float *)mxGetPr( plhs[2] ))==NULL )
			mexErrMsgTxt("fstDerivatives: error obtaining pointer to output variable 'Idy'");
	}

	if(temp.ndims==2)
		temp.data = (float*)mxCalloc( temp.dimElems[0]*temp.dimElems[1], sizeof(float) );
	else
		temp.data = (float*)mxCalloc( temp.dimElems[0]*temp.dimElems[1]*temp.dimElems[2], sizeof(float) );

	
	fstSimoncelli_c( &Idt, &Idx, &Idy, &It0, &It1, &temp, 5, Smoother, SpatialDerivator, TemporalDerivator );
	
	mxFree( temp.data );

}
