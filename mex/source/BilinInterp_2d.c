/*
function Iout = BilinInterp_2d(Iin, X, Y)

Arbitrary image rewarping using bilinear interpolation. A MEX c-program for Matlab.
This is a single-threaded 'simple' version.

INPUT
Iin		=		Original image, MxNxS matrix.
X		=		Remap x-coordinates, MxN matrix.
Y		=		Remap y-coordinates, MxN matrix.
OUTPUT
Iout		=		Interpolated ("warped") image

Author: Jarno Ralli
E-mail: jarno@ralli.fi

If you use this code, please reference (some) of my papers available at http://www.jarnoralli.fi

Copyright 2011-2014, Jarno Ralli

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
 #include "imageInterpolation.h"
 
/*
-------------------
--- mexFunction ---
-------------------
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	struct matrixM Iin;
	struct matrixM Iout;
	struct matrixM X;
	struct matrixM Y;
	
	/* 
	//-----------------------------------------------------------
	//--- CHECK FOR CORRECT PARAMETERS INTO/FROM THE FUNCTION ---
	//-----------------------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 3) {
		mexErrMsgTxt("proper function call is 'bilinInterp2( Iin, X, Y)'");
	} else if (nlhs > 3) {
		mexErrMsgTxt("proper function call is 'bilinInterp2( Iin, X, Y)'");
	}
	
	/*
	//------------
	//--- Iin ---
	//------------
	*/
	/* Iin must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("error: 'Iin' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Iin.ndims = mxGetNumberOfDimensions(prhs[0]);
	Iin.dimElems = mxGetDimensions(prhs[0]);
	Iin.data = (float*)mxGetPr( prhs[0] );
	
	/*
	//------------
	//--- X ---
	//------------
	*/
	/* X must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("error: 'X' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	X.ndims = mxGetNumberOfDimensions(prhs[1]);
	X.dimElems = mxGetDimensions(prhs[1]);
	X.data = (float*)mxGetPr( prhs[1] );
	
	/*
	//------------
	//--- Y ---
	//------------
	*/
	/* Y must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[2])){
		mexErrMsgTxt("error: 'Y' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	Y.ndims = mxGetNumberOfDimensions(prhs[2]);
	Y.dimElems = mxGetDimensions(prhs[2]);
	Y.data = (float*)mxGetPr( prhs[2] );
	
	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT OUTPUT PARAMETERS ---
	//------------------------------------------
	*/
	Iout.ndims = Iin.ndims;
	Iout.dimElems = Iin.dimElems;

	if( nlhs<1 ){
		mexErrMsgTxt("insufficient number of outputs. Outputs from this function is 'Iout'");
	}else{
	
		if( (plhs[0] = mxCreateNumericArray(Iout.ndims,Iout.dimElems, mxSINGLE_CLASS, mxREAL))==NULL )
			mexErrMsgTxt("error reserving space for output variable 'Iout'");
		if( (Iout.data = (float *)mxGetPr( plhs[0] ))==NULL )
			mexErrMsgTxt("error obtaining pointer to output variable 'Iout'");
	}
	
	bilinInterp2(	&Iout,
			&Iin,
			&X,
			&Y);
}
