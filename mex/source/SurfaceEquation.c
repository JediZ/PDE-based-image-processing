/*
function [M_out varargout] = SurfaceEquation(A_in, B_in, M_in, err_thr, min_set_size, iter )

Searches for parameters of a 2D surface using RANSAC. Surface is modeled as a first or second order
multivariate polynomial: [X Y 1] or [X^2 Y^2 XY X Y 1]. Order is automatically detected from size of A_in.

A_in * M_out = B_in

A MEX c-program for Matlab.

INPUT
A_in		=		Homogeneous coordinate matrix of the surface [X Y 1] or [X^2 Y^2 XY X Y 1]
B_in		=		'Height' of the surface [Z].
M_in		=		Model in. Previously found 'good' model to compete with the newly generated models.
err_thr		=		Error threshold to accept a point into inliers.
min_set_size	=		Minimum set size (0..1) of inliers in percentage
iter		=		Number of iterations.

OUTPUT
M_out		=		Model out.

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

#include "ransac.h"
#include "blas.h"
#include "lapack.h"

#define M 20		/*Maximum number of individuals used in RANSAC*/
#define N 6		/*Maximum model vector length*/
#define LWORK 198

/*Multivariate polynomial*/
void mvarPolynomial( struct matrixM *A_in,
		     struct matrixM *B_in,
		     struct matrixM *model_out,
		     struct matrixM *model_in,
		     struct matrixM *err_out, 
		     unsigned int *rand_set, 
		     unsigned int n );

/*
//------------------
//--- MAIN (MEX) ---
//------------------
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
	/*Structs*/
	struct matrixM A_in = {0,NULL,NULL,0,0,0,0};
	struct matrixM B_in = {0,NULL,NULL,0,0,0,0};
	struct matrixM error_out = {0,NULL,NULL,0,0,0,0};
	struct matrixM M_in = {0,NULL,NULL,0,0,0,0};
	struct matrixM M_out = {0,NULL,NULL,0,0,0,0};

	float err_thr;
	float min_set_size;
	float iter;
	
	mwSize temp[2] = {0, 1};

	/* 
	//------------------------------------------
	//--- CHECK FOR CORRECT INPUT PARAMETERS ---
	//------------------------------------------
	*/
	/* Check for proper number of arguments. */
	if (nrhs != 6) {
		mexErrMsgTxt("SurfaceEquation error: wrong number of input parameters!");
	} 
	/*
	//---------------
	//--- A_in ---
	//---------------
	*/
	/* A_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[0])){
		mexErrMsgTxt("SurfaceEquation error: 'A_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	A_in.ndims = mxGetNumberOfDimensions(prhs[0]);
	A_in.dimElems = mxGetDimensions(prhs[0]);
	A_in.data = (float*)mxGetPr( prhs[0] );
	
	/*
	//---------------
	//--- B_in ---
	//---------------
	*/
	/* B_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[1])){
		mexErrMsgTxt("SurfaceEquation error: 'B_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	B_in.ndims = mxGetNumberOfDimensions(prhs[1]);
	B_in.dimElems = mxGetDimensions(prhs[1]);
	B_in.data = (float*)mxGetPr( prhs[1] );
	/*
	//------------
	//--- M in ---
	//------------
	*/
	/* M_in must be a noncomplex matrix. */
	if (!mxIsSingle(prhs[2])){
		mexErrMsgTxt("SurfaceEquation error: 'M_in' must be a noncomplex single-valued matrix.");
	}
	/*Get number of dimensions and elements per dimension */
	M_in.ndims = mxGetNumberOfDimensions(prhs[2]);
	M_in.dimElems = mxGetDimensions(prhs[2]);
	M_in.data = (float*)mxGetPr( prhs[2] );
	/*
	//-----------------------
	//--- Error threshold ---
	//-----------------------
	*/
	/* err_thr must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[3])){
		mexErrMsgTxt("SurfaceEquation error: 'err_thr' must be a noncomplex, single-type scalar");
	}
	err_thr = *((float*)mxGetPr( prhs[3] ) );
	/*
	//------------------------
	//--- Minimum set size ---
	//------------------------
	*/
	/* min_set_size must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[4])){
		mexErrMsgTxt("SurfaceEquation error: 'min_set_size' must be a noncomplex, single-type scalar");
	}
	min_set_size = *((float*)mxGetPr( prhs[4] ) );
	/*
	//-----------------------------------
	//--- Number of iterations ---
	//-----------------------------------
	*/
	/* min_set_size must be a noncomplex scalar. */
	if (!mxIsSingle(prhs[5])){
		mexErrMsgTxt("SurfaceEquation error: 'iter' must be a noncomplex, single-type scalar");
	}
	iter = *((float*)mxGetPr( prhs[5] ) );

	/*INPUT VALIDITY CHECKING*/
	if( M_in.data != NULL )
	{
	  if( A_in.dimElems[1] != M_in.dimElems[0] )
	  {	mexErrMsgTxt("SurfaceEquation error: M_in is a column vector with as many row elements as A_in has columns!");
	  }
	}
	if( A_in.dimElems[0] != B_in.dimElems[0] )
	  mexErrMsgTxt("SurfaceEquation error: A_in and B_in have to have same amount of rows!");

	switch( A_in.dimElems[1] )
	{
	  case 3:
		temp[0] = 3;
		break;
		
	  case 6:
		temp[0] = 6;
		break;
	  
	  default:
	      mexPrintf("A_in.dimElems[0]:%d, A_in.dimElems[1]:%d\n",A_in.dimElems[0], A_in.dimElems[1]);
	      mexErrMsgTxt("SurfaceEquation error: only 1st and 2nd order polynomials are implemented!");
	}
	
	if( nlhs<2 )
	{
		mexErrMsgTxt("SurfaceEquation error: insufficient number of outputs. Outputs from this function is 'M_out' and 'error_out'");
	}else{
		if( (plhs[0] = mxCreateNumericArray(2, temp, mxSINGLE_CLASS, mxREAL))==NULL )
		  mexErrMsgTxt("SurfaceEquation error: error reserving space for output variable 'M_out'");
		if( (M_out.data = (float *)mxGetPr( plhs[0] ))==NULL )
		  mexErrMsgTxt("SurfaceEquation error: error obtaining pointer to output variable 'M_out'");
		M_out.ndims = mxGetNumberOfDimensions(plhs[0]);
		M_out.dimElems = mxGetDimensions(plhs[0]);
		
		temp[0] = A_in.dimElems[0];
		if( (plhs[1] = mxCreateNumericArray(2, temp, mxSINGLE_CLASS, mxREAL))==NULL )
		  mexErrMsgTxt("SurfaceEquation error: error reserving space for output variable 'M_out'");
		if( (error_out.data = (float *)mxGetPr( plhs[1] ))==NULL )
		  mexErrMsgTxt("SurfaceEquation error: error obtaining pointer to output variable 'error_out'");
		error_out.ndims = mxGetNumberOfDimensions(plhs[1]);
		error_out.dimElems = mxGetDimensions(plhs[1]);
	}
		
	RANSAC(	&A_in,					/*Homogeneous coordinate matrix (e.g [X Y 1] or [X^2 Y^2 XY X Y 1])*/
		&B_in,					/*Surface 'height' [Z]*/
		&M_out,					/*The best model found by RANSAC*/
		&M_in,					/*The best model found by RANSAC*/
		&error_out,				/*Residual (or fitting error) of each data for the best model*/
		&mvarPolynomial,
		err_thr, 				/*Error threshold for determining if the data fits the model*/
		min_set_size, 				/*Minimum number of close data (percentage; 0.0-1.0) values required to assert that a model fits well to data*/
		(unsigned int)iter, 			/*Number of iterations*/
		(unsigned int)M_out.dimElems[0]+1 );	/*Number of data used for generating the model*/

}

/*Multivariate polynomial*/
void mvarPolynomial( struct matrixM *A_in,
		     struct matrixM *B_in,
		     struct matrixM *model_out,
		     struct matrixM *model_in,
		     struct matrixM *err_out, 
		     unsigned int *rand_set, 
		     unsigned int nd )
{
	float A[M*N], B[M];
	float work[LWORK];
	mwSize lwork = LWORK;
	mwSize	m = nd,
		Arows = A_in->dimElems[0],
		Acols = A_in->dimElems[1],
		Brows = B_in->dimElems[0],
		Bcols = B_in->dimElems[1],
		Mrows = model_out->dimElems[0],
		Mcols = model_out->dimElems[1],
		ipiv[N],
		info;
	
	unsigned int i, j, colOffset, colOffset2, colOffset3, colOffset4, colOffset5;
	char *chn = "N";
	/*for dgemm */
	float alpha = 1.0f, beta = -1.0f;
	
	colOffset = A_in->dimElems[0];
	colOffset2 = 2*colOffset;
	colOffset3 = 3*colOffset;
	colOffset4 = 4*colOffset;
	colOffset5 = 5*colOffset;
	
	/*rand_set[0] = 0;
	rand_set[1] = 1;
	rand_set[2] = 2;
	rand_set[3] = 3;
	rand_set[4] = 4;
	rand_set[5] = 5;
	rand_set[6] = 6;*/
	
	if( model_in == NULL )
	{
	  /*First or second order multivariate polynomial*/
	  switch( model_out->dimElems[0] )
	  {
	    /*1st order*/
	    case 3:

	      for(i=0;i<nd;i++)
	      {
		  A[ i ] = A_in->data[ rand_set[i]  ];
		  A[ i + nd ] = A_in->data[ rand_set[i] + colOffset ];
		  A[ i + 2*nd ] = A_in->data[ rand_set[i] + colOffset2 ];
		  
		  B[i] = B_in->data[ rand_set[i] ];
	      }
	      
/*	      A[0] = A_in->data[rand_set[0]];
	      A[1] = A_in->data[rand_set[1]];
	      A[2] = A_in->data[rand_set[2]];
	      
	      A[3] = A_in->data[rand_set[0]+colOffset];
	      A[4] = A_in->data[rand_set[1]+colOffset];
	      A[5] = A_in->data[rand_set[2]+colOffset];
	      
	      A[6] = 1.0f;
	      A[7] = 1.0f;
	      A[8] = 1.0f;
    
	      model_out->data[0] = B_in->data[rand_set[0]];
	      model_out->data[1] = B_in->data[rand_set[1]];
	      model_out->data[2] = B_in->data[rand_set[2]];
*/
	    break;
	    /*2nd order*/
	    case 6:
	      for(i=0;i<nd;i++)
	      {
		  A[ i ] = A_in->data[ rand_set[i]  ];
		  A[ i + nd ] = A_in->data[ rand_set[i] + colOffset ];
		  A[ i + 2*nd ] = A_in->data[ rand_set[i] + colOffset2 ];
		  A[ i + 3*nd ] = A_in->data[ rand_set[i] + colOffset3 ];
		  A[ i + 4*nd ] = A_in->data[ rand_set[i] + colOffset4 ];
		  A[ i + 5*nd ] = A_in->data[ rand_set[i] + colOffset5 ];
		  
		  B[i] = B_in->data[ rand_set[i] ];
	      }

/*	      A[0] = A_in->data[ rand_set[0] ];
	      A[1] = A_in->data[ rand_set[1] ];
	      A[2] = A_in->data[ rand_set[2] ];
	      A[3] = A_in->data[ rand_set[3] ];
	      A[4] = A_in->data[ rand_set[4] ];
	      A[5] = A_in->data[ rand_set[5] ];
	      
	      A[6] = A_in->data[ rand_set[0] + colOffset ];
	      A[7] = A_in->data[ rand_set[1] + colOffset ];
	      A[8] = A_in->data[ rand_set[2] + colOffset ];
	      A[9] = A_in->data[ rand_set[3] + colOffset ];
	      A[10] = A_in->data[ rand_set[4] + colOffset ];
	      A[11] = A_in->data[ rand_set[5] + colOffset ];
	      
	      A[12] = A_in->data[ rand_set[0] + colOffset2 ];
	      A[13] = A_in->data[ rand_set[1] + colOffset2 ];
	      A[14] = A_in->data[ rand_set[2] + colOffset2 ];
	      A[15] = A_in->data[ rand_set[3] + colOffset2 ];
	      A[16] = A_in->data[ rand_set[4] + colOffset2 ];
	      A[17] = A_in->data[ rand_set[5] + colOffset2 ];
	      
	      A[18] = A_in->data[ rand_set[0] + colOffset3 ];
	      A[19] = A_in->data[ rand_set[1] + colOffset3 ];
	      A[20] = A_in->data[ rand_set[2] + colOffset3 ];
	      A[21] = A_in->data[ rand_set[3] + colOffset3 ];
	      A[22] = A_in->data[ rand_set[4] + colOffset3 ];
	      A[23] = A_in->data[ rand_set[5] + colOffset3 ];
	      
	      A[24] = A_in->data[ rand_set[0] + colOffset4 ];
	      A[25] = A_in->data[ rand_set[1] + colOffset4 ];
	      A[26] = A_in->data[ rand_set[2] + colOffset4 ];
	      A[27] = A_in->data[ rand_set[3] + colOffset4 ];
	      A[28] = A_in->data[ rand_set[4] + colOffset4 ];
	      A[29] = A_in->data[ rand_set[5] + colOffset4 ];
	      
	      A[30] = A_in->data[ rand_set[0] + colOffset5 ];
	      A[31] = A_in->data[ rand_set[1] + colOffset5 ];
	      A[32] = A_in->data[ rand_set[2] + colOffset5 ];
	      A[33] = A_in->data[ rand_set[3] + colOffset5 ];
	      A[34] = A_in->data[ rand_set[4] + colOffset5 ];
	      A[35] = A_in->data[ rand_set[5] + colOffset5 ];
	      
	      model_out->data[0] = B_in->data[rand_set[0]];
	      model_out->data[1] = B_in->data[rand_set[1]];
	      model_out->data[2] = B_in->data[rand_set[2]];
	      model_out->data[3] = B_in->data[rand_set[3]];
	      model_out->data[4] = B_in->data[rand_set[4]];
	      model_out->data[5] = B_in->data[rand_set[5]];
*/
	      break;
	    
	    default:
	      mexErrMsgTxt("mvarPolynomial: only 1st and 2nd order multivariate polynomials are implemented!!");
	  }
	      
/*	    Solve the linear equation A*model_out = B ( B is stored in model_out )	*/
/*	    sgesv(	&n,		/*the number of linear equations		*/
/*			&Bcols,		/*number of columns in B (nrhs)			*/
/*			A,
/*			&n,		/*leading dimension of A, lda = max(1,n)	*/
/*			ipiv,		/*pivot indices, size n				*/
/*			model_out->data,
/*			&n,		/*leading dimension of B, ldb = max(1,n)	*/
/*			&info );							*/

	    sgels(	chn,
			&m,
			&Acols,
			&Bcols,
			A,
			&m,
			B,
			&m,
			work,
			&lwork,
			&info );
	    memcpy( model_out->data, B, model_out->dimElems[0]*sizeof(float) );
		
	}else
	{
	  memcpy( model_out->data, model_in->data, model_in->dimElems[0]*sizeof(float) );
	  info = 0;
	}

	if( info==0 )
	{
	  /* C = alpha*A*B + beta*C */
	  /* Using out variables: err_out = alpha*A*model_out + beta*err_out */
	  /*Calculate the error between the model and the data*/
	  memcpy( err_out->data, B_in->data, B_in->dimElems[0]*B_in->dimElems[1]*sizeof(float) );
	  sgemm(	chn,			/*transA*/
			chn,			/*transB*/
			&Arows,			/*m, number of rows of A and C*/
			&Mcols,			/*n, number of columns of B*/
			&Acols,			/*k, number of columns of A and number of rows of B*/
			&alpha,			/*alpha*/
			A_in->data, 		/*A*/
			&Arows,			/*lda = max(1,m)*/
			model_out->data, 	/*B*/
			&Acols,			/*ldb = max(1,k)*/
			&beta, 			/*beta*/
			err_out->data,		/*C*/
			&Arows);		/*ldc = max(1,m)*/
	  /*Quadratic error*/
	  for(i=0;i<err_out->dimElems[0];i++)
	    err_out->data[i] *= err_out->data[i];
	}else
	{	
	    for(i=0;i<err_out->dimElems[0];i++)
	      err_out->data[i] = FLT_MAX;
	}
	
}


