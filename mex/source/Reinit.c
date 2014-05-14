/*
function [PHI_out varargout] = reInitC(PHIin, T)

Re-initializes a level-set function (or functions).

A MEX c-program for Matlab.

INPUT
PHIin		=		Level set function PHI.
T		=		Time, 0:0.25:T (i.e. how many steps to perform).

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

	float T = 0.0f;
	unsigned int elems;
	unsigned int i;
	
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
	if (nrhs != 2) {
		mexErrMsgTxt("reInitC parameter error: wrong number of input parameters!");
	} 
	/*
	//------------
	//--- PHI in ---
	//------------
	*/
	/* PHI_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("reInitC: 'PHI_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	PHI_in.ndims = mxGetNumberOfDimensions(prhs[0]);
	PHI_in.dimElems = mxGetDimensions(prhs[0]);
	PHI_in.data = (float*)mxGetPr( prhs[0] );

	/*
	//-----------------------------------
	//--- Time T ---
	//-----------------------------------
	*/
	/* T must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("reInitC error: 'T' must be a noncomplex, single-type scalar");
	}
	T = *((float*)mxGetPr( prhs[1] ) );

	
	/* Check if maximum buffer size is big enough */
	switch( PHI_in.ndims )
	{
	  case 1:
		if (PHI_in.dimElems[0]>MAX_BUF_SIZE)
		    mexErrMsgTxt("reInitC error: MAX_BUF_SIZE too small...increase the size in the header file!");
		break;
	  default :
		if ( (PHI_in.dimElems[0]>MAX_BUF_SIZE) || (PHI_in.dimElems[1]>MAX_BUF_SIZE) )
		  mexErrMsgTxt("reInitC error: MAX_BUF_SIZE too small...increase the size in the header file!");
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
		mexErrMsgTxt("reInitC error insufficient number of outputs. Outputs from this function is 'PHI_out'");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(PHI_out.ndims,PHI_out.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("reInitC error: error reserving space for output variable 'PHI_out'");

		if( (PHI_out.data = (float *)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("reInitC error: error obtaining pointer to output variable 'PHI_out'");
	}
	
	elems = PHI_in.dimElems[0];
	
	for(i=1;i<PHI_out.ndims;i++)
		elems *= PHI_in.dimElems[i];
	
	reinit(	&PHI_in, T);
	memcpy( PHI_out.data, PHI_in.data, elems*sizeof(float) );

}
