 /*
Functions related to PDE solvers.
 
Solvers search for a solution to Ax = B, such that
-TRACE = TRACE(A)
-x is the solution

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

#include "pdeSolvers.h"

#if !defined( isnan )
  #define isnan(x)((x)!=(x))
#endif

/*
//---------------------------------------------------------------------------------------
//---Point-wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_4_2d(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam)
{
	float wNeigh, *INV_TRACE, *B_temp;
	unsigned int i,j,k,pos,iter,wpos,npos,epos,spos;
	unsigned int nrows = X->dimElems[0];
	unsigned int ncols = X->dimElems[1];
	unsigned int nframes = 1;
	unsigned int iterations = (int)Mparam.iter;

	if(X->ndims>2)
		nframes = X->dimElems[2];

	if( (INV_TRACE = (float*)malloc( nrows*ncols*nframes*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_4_2d: error reserving memory for 'INV_TRACE'\n");
		return;
	}
	if( (B_temp = (float*)malloc( nrows*ncols*nframes*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_4_2d: error reserving memory for 'B_temp'\n");
		return;
	}

	for(iter=0;iter<iterations;iter++)
	{
	for(k=0;k<nframes;k++)
	{
		for(j=1;j<ncols-1;j++)
		{
		for(i=1;i<nrows-1;i++)		
		{
			/*Spatial position indices*/
			pos = j*nrows+i+k*nrows*ncols;
			wpos = pos-nrows;
			epos = pos+nrows;
			npos = pos-1;
			spos = pos+1;

			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			wNeigh =	X->data[epos]*wE->data[ pos ] +
					X->data[wpos]*wW->data[ pos ];
			wNeigh +=	X->data[spos]*wS->data[ pos ] +
					X->data[npos]*wN->data[ pos ];

			if(iter==0)
			{
				if(!isnan(TRACE->data[pos]))
				{
					INV_TRACE[pos] = 1.0f/TRACE->data[pos];
					B_temp[pos] = B->data[pos];
				}
				else
				{	
					INV_TRACE[pos] = 	wE->data[ pos ] +
								wW->data[ pos ];
					INV_TRACE[pos] += 	wS->data[ pos ] +
								wN->data[ pos ];
					INV_TRACE[pos] = 1.0f/INV_TRACE[pos];
					B_temp[pos] = 0.0f;
				}
			}
			/*SOR (Successive Over Relaxation)*/
			X->data[pos] = (1.0f-Mparam.omega)*X->data[pos];			/*old solution*/
			X->data[pos] += Mparam.omega*(B_temp[pos]+wNeigh)*INV_TRACE[pos];	/*new solution*/
		}
		}

		/*
		-------------------
		--Fill the borders
		-------------------
		*/
		for(j=0;j<ncols;j++)
		{
			pos = j*nrows+k*nrows*ncols;;

			X->data[ pos ] = X->data[ pos+1 ];
			X->data[ pos+nrows-1 ] = X->data[ pos+nrows-2 ];
		}
		for(i=0;i<nrows;i++)
		{
			pos = i+k*nrows*ncols;

			X->data[ pos ] = X->data[ pos+nrows ];
			X->data[ pos+(ncols-1)*nrows ] = X->data[ pos+(ncols-2)*nrows ];
		}
	}
	}
	
	free(INV_TRACE);
	free(B_temp);
}

/*
//---------------------------------------------------------------------------------------
//---Point-wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_8_2d(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam)
{
	float wNeigh, *INV_TRACE, *B_temp;
	unsigned int i,j,k,pos,iter,wpos,npos,epos,spos;
	unsigned int nrows = X->dimElems[0];
	unsigned int ncols = X->dimElems[1];
	unsigned int nframes = 1;
	unsigned int iterations = (int)Mparam.iter;

	if(X->ndims>2)
		nframes = X->dimElems[2];

	if( (INV_TRACE = (float*)malloc( nrows*ncols*nframes*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_8_2d: error reserving memory for 'INV_TRACE'\n");
		return;
	}
	if( (B_temp = (float*)malloc( nrows*ncols*nframes*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_8_2d: error reserving memory for 'B_temp'\n");
		return;
	}

	for(iter=0;iter<iterations;iter++)
	{
	for(k=0;k<nframes;k++)
	{
		for(j=1;j<ncols-1;j++)
		{
		for(i=1;i<nrows-1;i++)		
		{
		  
			/*Spatial position indices*/
			pos = j*nrows+i+k*nrows*ncols;
			wpos = pos-nrows;
			epos = pos+nrows;
			npos = pos-1;
			spos = pos+1;
			
			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			wNeigh =	X->data[epos]*wE->data[ pos ] +
					X->data[wpos]*wW->data[ pos ];
			wNeigh +=	X->data[spos]*wS->data[ pos ] +
					X->data[npos]*wN->data[ pos ];
			wNeigh +=	X->data[wpos+1]*wSW->data[ pos ] +
					X->data[wpos-1]*wNW->data[ pos ];
			wNeigh +=	X->data[epos+1]*wSE->data[ pos ] +
					X->data[epos-1]*wNE->data[ pos ];

			if(iter==0)
			{
				if(!isnan(TRACE->data[pos]))
				{
					INV_TRACE[pos] = 1.0f/TRACE->data[pos];
					B_temp[pos] = B->data[pos];
				}
				else
				{	
					INV_TRACE[pos] = 	wE->data[ pos ] +
								wW->data[ pos ];
					INV_TRACE[pos] += 	wS->data[ pos ] +
								wN->data[ pos ];
					INV_TRACE[pos] += 	wSW->data[ pos ] +
								wNW->data[ pos ];
					INV_TRACE[pos] += 	wSE->data[ pos ] +
								wNE->data[ pos ];
					INV_TRACE[pos] = 1.0f/INV_TRACE[pos];
					B_temp[pos] = 0.0f;
				}
			}
			/*SOR (Successive Over Relaxation)*/
			X->data[pos] = (1.0f-Mparam.omega)*X->data[pos];
			X->data[pos] += Mparam.omega*(B_temp[pos]+wNeigh)*INV_TRACE[pos];
		}
		}

		/*
		-------------------
		--Fill the borders
		-------------------
		*/
		for(j=0;j<ncols;j++)
		{
			pos = j*nrows+k*nrows*ncols;;

			X->data[ pos ] = X->data[ pos+1 ];
			X->data[ pos+nrows-1 ] = X->data[ pos+nrows-2 ];
		}
		for(i=0;i<nrows;i++)
		{
			pos = i+k*nrows*ncols;

			X->data[ pos ] = X->data[ pos+nrows ];
			X->data[ pos+(ncols-1)*nrows ] = X->data[ pos+(ncols-2)*nrows ];
		}
	}
	}
	
	free(INV_TRACE);
	free(B_temp);
}

/*
//---------------------------------------------------------------------------------------
//---Alternating Line Relaxation (ALR) scheme
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 4-neighbours
//---------------------------------------------------------------------------------------
*/
void GS_ALR_SOR_4_2d(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam)
{
	float *Cp, *Dp;
	int nrows = X->dimElems[0];
	int ncols = X->dimElems[1];
	int nframes = 1;
	int maxLength;
	int iterations = (int)Mparam.iter;
	int iter;
	int i,j,k,pos;

	if(X->ndims>2)
		nframes = X->dimElems[2];
	
	maxLength = (nrows>ncols)?nrows:ncols;

	/*Reserve space for Cprime and Dprime temp variables used in TDMA*/
	if( (Cp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("Error reserving space for 'Cp'\n");
		return;
	}
	if( (Dp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("Error reserving space for 'Dp'\n");
		return;
	}

	for(iter=0;iter<iterations;iter++)
	{
		/*
		------------------------
		-Vertical line relaxing
		------------------------
		*/
		TDMA_wcolumn_ALR_4(	X, TRACE, B, wW, wN, wE, wS, Cp, Dp, ncols, nrows, nframes, Mparam.omega);
		TDMA_mcolumn_ALR_4(	X, TRACE, B, wW, wN, wE, wS, Cp, Dp, ncols, nrows, nframes, Mparam.omega);
		TDMA_ecolumn_ALR_4(	X, TRACE, B, wW, wN, wE, wS, Cp, Dp, ncols, nrows, nframes, Mparam.omega);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		TDMA_nrow_ALR_4(		X, TRACE, B, wW, wN, wE, wS, Cp, Dp, ncols, nrows, nframes, Mparam.omega);
		TDMA_mrow_ALR_4(		X, TRACE, B, wW, wN, wE, wS, Cp, Dp, ncols, nrows, nframes, Mparam.omega);
		TDMA_srow_ALR_4(		X, TRACE, B, wW, wN, wE, wS, Cp, Dp, ncols, nrows, nframes, Mparam.omega);
	}
	
	free( Cp );
	free( Dp );
}

/*
//---------------------------------------------------------------------------------------
//---Alternating Line Relaxation (ALR) scheme
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 8-neighbours
//---------------------------------------------------------------------------------------
*/
void	GS_ALR_SOR_8_2d(struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam)
{
	float *Cp, *Dp;
	int nrows = X->dimElems[0];
	int ncols = X->dimElems[1];
	int nframes = 1;
	int maxLength;
	int iterations = 1;
	int iter;

	if(X->ndims>2)
		nframes = X->dimElems[2];
	
	maxLength = (nrows>ncols)?nrows:ncols;

	/*Reserve space for Cprime and Dprime temp variables used in TDMA*/
	if( (Cp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("Error reserving space for 'Cp'\n");
		return;
	}
	if( (Dp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("Error reserving space for 'Dp'\n");
		return;
	}

	for(iter=0;iter<iterations;iter++)
	{
		/*
		------------------------
		-Vertical line relaxing
		------------------------
		*/
		TDMAcolumn_ALR_8(	X, TRACE, B, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Cp, Dp, ncols, nrows, nframes, Mparam.omega);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		TDMArow_ALR_8(		X, TRACE, B, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Cp, Dp, ncols, nrows, nframes, Mparam.omega);

	}

	free( Cp );
	free( Dp );
}

/*
//------------------------------------------
// TDMA vertical line solving - west column
//------------------------------------------
*/
void TDMA_wcolumn_ALR_4(	struct matrixM *X,
				struct matrixM *TRACE,
				struct matrixM *B,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				float *cp,
				float *dp,
				int ncols,
				int nrows,
				int nframes,
				float omega)
{
	int i,j,k,pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	for(k=0;k<nframes;k++)
	{
		j = 0;
		/*---First row---*/
		i = 0;
		pos = k*nrows*ncols+j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;
	
		/*a = 	0.0f;*/
		/*b = 	TRACE->data[pos];*/
		c = 	-wS->data[pos];
		d =	wE->data[pos]*X->data[epos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wS->data[pos] + wE->data[pos];
		}
		
		cp[i] = c/b;
		dp[i] = d/b;
	
		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos++;
			wpos++;
			epos++;
			npos++;
			spos++;
	
			a = 	-wN->data[pos];
			/*b = 	TRACE->data[pos];*/
			c = 	-wS->data[pos];
			d =	wE->data[pos]*X->data[epos];
			
			if(!isnan(TRACE->data[pos]))
			{
			    b = TRACE->data[pos];
			    d += B->data[pos];
			}
			else
			{
			    b = wN->data[pos] + wS->data[pos] + wE->data[pos];
			}
				
			div = 1 / (b-cp[i-1]*a);
			cp[i] = c*div;
			dp[i] = (d-dp[i-1]*a)*div;
		}
	
		/*---Last row---*/
		pos++;
		wpos++;
		epos++;
		npos++;
		spos++;
	
		a = 	-wN->data[pos];
		/*b = 	TRACE->data[pos];*/
		/*c = 	0.0f;*/
		d =	wE->data[pos]*X->data[epos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wN->data[pos] + wE->data[pos];
		}

		dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
		temp1 = X->data[pos];
		X->data[pos] = dp[i];
	
		for(i=nrows-2;i>=0;i--)
		{
			pos = k*nrows*ncols+j*nrows+i;
			temp2 = X->data[pos];
			X->data[pos] = dp[i] - cp[i]*X->data[pos+1];

			X->data[pos+1] = omega*X->data[pos+1] + (1.0f-omega)*temp1;
			temp1 = temp2;
		}
		X->data[pos] = omega*X->data[pos] + (1.0f-omega)*temp1;

	}

}

/*
//--------------------------------------------
// TDMA vertical line solving - middle columns
//--------------------------------------------
*/
void TDMA_mcolumn_ALR_4(	struct matrixM *X,
				struct matrixM *TRACE,
				struct matrixM *B,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				float *cp,
				float *dp,
				int ncols,
				int nrows,
				int nframes,
				float omega)
{
	int i,j,k,pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	for(k=0;k<nframes;k++)
	{
	for(j=1;j<=ncols-2;j++)
	{
		/*---First row---*/
		i = 0;
		pos = k*nrows*ncols+j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;
	
		/*a = 	0.0f;*/
		/*b = 	TRACE->data[pos];*/
		c = 	-wS->data[pos];
		d =	wW->data[pos]*X->data[wpos] + wE->data[pos]*X->data[epos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wS->data[pos] + wW->data[pos]+ wE->data[pos];
		}
		
		cp[i] = c/b;
		dp[i] = d/b;
	
		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos++;
			wpos++;
			epos++;
			npos++;
			spos++;
	
			a = 	-wN->data[pos];
			/*b = 	TRACE->data[pos];*/
			c = 	-wS->data[pos];
			d =	wW->data[pos]*X->data[wpos] + wE->data[pos]*X->data[epos];
			
			if(!isnan(TRACE->data[pos]))
			{
			    b = TRACE->data[pos];
			    d += B->data[pos];
			}
			else
			{
			    b = wN->data[pos] + wS->data[pos]  + wW->data[pos]+ wE->data[pos];
			}
				
			div = 1 / (b-cp[i-1]*a);
			cp[i] = c*div;
			dp[i] = (d-dp[i-1]*a)*div;
		}
	
		/*---Last row---*/
		pos++;
		wpos++;
		epos++;
		npos++;
		spos++;
	
		a = 	-wN->data[pos];
		/*b = 	TRACE->data[pos];*/
		/*c = 	0.0f;*/
		d =	wW->data[pos]*X->data[wpos] + wE->data[pos]*X->data[epos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wN->data[pos] + wW->data[pos]+ wE->data[pos];
		}

		dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
		temp1 = X->data[pos];
		X->data[pos] = dp[i];
	
		for(i=nrows-2;i>=0;i--)
		{
			pos = k*nrows*ncols+j*nrows+i;
			temp2 = X->data[pos];
			X->data[pos] = dp[i] - cp[i]*X->data[pos+1];

			X->data[pos+1] = omega*X->data[pos+1] + (1.0f-omega)*temp1;
			temp1 = temp2;
		}
		X->data[pos] = omega*X->data[pos] + (1.0f-omega)*temp1;

	}
	}
}

/*
//--------------------------------------------
// TDMA vertical line solving - middle columns
//--------------------------------------------
*/
void TDMA_ecolumn_ALR_4(	struct matrixM *X,
				struct matrixM *TRACE,
				struct matrixM *B,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				float *cp,
				float *dp,
				int ncols,
				int nrows,
				int nframes,
				float omega)
{
	int i,j,k,pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	for(k=0;k<nframes;k++)
	{
		j = ncols-1;
		/*---First row---*/
		i = 0;
		pos = k*nrows*ncols+j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;
	
		/*a = 	0.0f;*/
		/*b = 	TRACE->data[pos];*/
		c = 	-wS->data[pos];
		d =	wW->data[pos]*X->data[wpos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wS->data[pos] + wW->data[pos];
		}
		
		cp[i] = c/b;
		dp[i] = d/b;
	
		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos++;
			wpos++;
			epos++;
			npos++;
			spos++;
	
			a = 	-wN->data[pos];
			/*b = 	TRACE->data[pos];*/
			c = 	-wS->data[pos];
			d =	wW->data[pos]*X->data[wpos];
			
			if(!isnan(TRACE->data[pos]))
			{
			    b = TRACE->data[pos];
			    d += B->data[pos];
			}
			else
			{
			    b = wN->data[pos] + wS->data[pos]  + wW->data[pos];
			}
				
			div = 1 / (b-cp[i-1]*a);
			cp[i] = c*div;
			dp[i] = (d-dp[i-1]*a)*div;
		}
	
		/*---Last row---*/
		pos++;
		wpos++;
		epos++;
		npos++;
		spos++;
	
		a = 	-wN->data[pos];
		/*b = 	TRACE->data[pos];*/
		/*c = 	0.0f;*/
		d =	wW->data[pos]*X->data[wpos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wN->data[pos] + wW->data[pos];
		}

		dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
		temp1 = X->data[pos];
		X->data[pos] = dp[i];
	
		for(i=nrows-2;i>=0;i--)
		{
			pos = k*nrows*ncols+j*nrows+i;
			temp2 = X->data[pos];
			X->data[pos] = dp[i] - cp[i]*X->data[pos+1];

			X->data[pos+1] = omega*X->data[pos+1] + (1.0f-omega)*temp1;
			temp1 = temp2;
		}
		X->data[pos] = omega*X->data[pos] + (1.0f-omega)*temp1;

	}
}


/*
//------------------------------------------
// TDMA horizontal line solving - north row
//------------------------------------------
*/
void TDMA_nrow_ALR_4(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega)
{
	int i,j,k,pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	for(k=0;k<nframes;k++)
	{
		i = 0;
		/*---First column---*/
		j = 0;
		pos = k*nrows*ncols+j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		/*a =	0.0f;*/
		/*b =	TRACE->data[pos];*/
		c =	-wE->data[pos];
		d =	wS->data[pos]*X->data[spos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wS->data[pos] + wE->data[pos];
		}

		cp[j] = c/b;
		dp[j] = d/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos+=nrows;
			wpos+=nrows;
			epos+=nrows;
			npos+=nrows;
			spos+=nrows;

			a = 	-wW->data[pos];	
			/*b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos];*/
			c = 	-wE->data[pos];
			d =	wS->data[pos]*X->data[spos];
			
			if(!isnan(TRACE->data[pos]))
			{
			    b = TRACE->data[pos];
			    d += B->data[pos];
			}
			else
			{
			    b = wS->data[pos]  + wW->data[pos]+ wE->data[pos];
			}
			
			div = 1 / (b-cp[j-1]*a);
			cp[j] = c*div;
			dp[j] = (d-dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos+=nrows;
		wpos+=nrows;
		epos+=nrows;
		npos+=nrows;
		spos+=nrows;

		a = 	-wW->data[pos];	
		/*b = 	wN->data[pos] + wS->data[pos] + wW->data[pos];*/
		/*c = 	0.0f;*/
		d =	wS->data[pos]*X->data[spos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wS->data[pos]  + wW->data[pos];
		}
	
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
		temp1 = X->data[pos];
		X->data[pos] = dp[j];

		for(j=ncols-2;j>=0;j--)
		{
			pos = k*nrows*ncols+j*nrows+i;
			temp2 = X->data[pos];
			X->data[pos] = dp[j] - cp[j]*X->data[pos+nrows];

			X->data[pos+nrows] = omega*X->data[pos+nrows] + (1.0f-omega)*temp1;
			temp1 = temp2;
		}
		X->data[pos] = omega*X->data[pos] + (1.0f-omega)*temp1;
	}
}

/*
//-------------------------------------------
// TDMA horizontal line solving - middle rows
//-------------------------------------------
*/
void TDMA_mrow_ALR_4(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega)
{
	int i,j,k,pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	for(k=0;k<nframes;k++)
	{
	for(i=1;i<=nrows-2;i++)
	{
		/*---First column---*/
		j = 0;
		pos = k*nrows*ncols+j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		/*a =	0.0f;*/
		/*b =	TRACE->data[pos];*/
		c =	-wE->data[pos];
		d =	wS->data[pos]*X->data[spos] + wN->data[pos]*X->data[npos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wN->data[pos] + wS->data[pos] + wE->data[pos];
		}

		cp[j] = c/b;
		dp[j] = d/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos+=nrows;
			wpos+=nrows;
			epos+=nrows;
			npos+=nrows;
			spos+=nrows;

			a = 	-wW->data[pos];	
			/*b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos];*/
			c = 	-wE->data[pos];
			d =	wS->data[pos]*X->data[spos] + wN->data[pos]*X->data[npos];
			
			if(!isnan(TRACE->data[pos]))
			{
			    b = TRACE->data[pos];
			    d += B->data[pos];
			}
			else
			{
			    b = wN->data[pos] + wS->data[pos]  + wW->data[pos]+ wE->data[pos];
			}
			
			div = 1 / (b-cp[j-1]*a);
			cp[j] = c*div;
			dp[j] = (d-dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos+=nrows;
		wpos+=nrows;
		epos+=nrows;
		npos+=nrows;
		spos+=nrows;

		a = 	-wW->data[pos];	
		/*b = 	wN->data[pos] + wS->data[pos] + wW->data[pos];*/
		/*c = 	0.0f;*/
		d =	wS->data[pos]*X->data[spos] + wN->data[pos]*X->data[npos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wN->data[pos] + wS->data[pos]  + wW->data[pos];
		}
	
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
		temp1 = X->data[pos];
		X->data[pos] = dp[j];

		for(j=ncols-2;j>=0;j--)
		{
			pos = k*nrows*ncols+j*nrows+i;
			temp2 = X->data[pos];
			X->data[pos] = dp[j] - cp[j]*X->data[pos+nrows];

			X->data[pos+nrows] = omega*X->data[pos+nrows] + (1.0f-omega)*temp1;
			temp1 = temp2;
		}
		X->data[pos] = omega*X->data[pos] + (1.0f-omega)*temp1;
	}
	}
}

/*
//-------------------------------------------
// TDMA horizontal line solving - south row
//-------------------------------------------
*/
void TDMA_srow_ALR_4(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega)
{
	int i,j,k,pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	for(k=0;k<nframes;k++)
	{
		i = nrows-1;
		/*---First column---*/
		j = 0;
		pos = k*nrows*ncols+j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		/*a =	0.0f;*/
		/*b =	TRACE->data[pos];*/
		c =	-wE->data[pos];
		d =	wN->data[pos]*X->data[npos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wN->data[pos] + wE->data[pos];
		}

		cp[j] = c/b;
		dp[j] = d/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos+=nrows;
			wpos+=nrows;
			epos+=nrows;
			npos+=nrows;
			spos+=nrows;

			a = 	-wW->data[pos];	
			/*b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos];*/
			c = 	-wE->data[pos];
			d =	wN->data[pos]*X->data[npos];
			
			if(!isnan(TRACE->data[pos]))
			{
			    b = TRACE->data[pos];
			    d += B->data[pos];
			}
			else
			{
			    b = wN->data[pos] + wW->data[pos]+ wE->data[pos];
			}
			
			div = 1 / (b-cp[j-1]*a);
			cp[j] = c*div;
			dp[j] = (d-dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos+=nrows;
		wpos+=nrows;
		epos+=nrows;
		npos+=nrows;
		spos+=nrows;

		a = 	-wW->data[pos];	
		/*b = 	wN->data[pos] + wS->data[pos] + wW->data[pos];*/
		/*c = 	0.0f;*/
		d =	wN->data[pos]*X->data[npos];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
		    b = wN->data[pos] + wW->data[pos];
		}
	
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
		temp1 = X->data[pos];
		X->data[pos] = dp[j];

		for(j=ncols-2;j>=0;j--)
		{
			pos = k*nrows*ncols+j*nrows+i;
			temp2 = X->data[pos];
			X->data[pos] = dp[j] - cp[j]*X->data[pos+nrows];

			X->data[pos+nrows] = omega*X->data[pos+nrows] + (1.0f-omega)*temp1;
			temp1 = temp2;
		}
		X->data[pos] = omega*X->data[pos] + (1.0f-omega)*temp1;
	}
}

/*
//----------------------------
// TDMA vertical line solving
//----------------------------
*/
void TDMAcolumn_ALR_8(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega)
{
	int i,j,k,pos,wpos,npos,epos,spos;
	float a,b,c,d,div,temp1,temp2;
	
	for(k=0;k<nframes;k++)
	{
	for(j=1;j<=ncols-2;j++)
	{
		/*---First row---*/
		i = 0;
		pos = k*nrows*ncols+j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;
	
		/*a = 	0.0f;*/
		/*b = 	TRACE->data[pos];*/
		c = 	-wS->data[pos];
		d =	wW->data[pos]*X->data[wpos] + wE->data[pos]*X->data[epos];
		d +=	wSW->data[pos]*X->data[wpos+1] + wSE->data[pos]*X->data[epos+1];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
			b =	wN->data[pos] + wS->data[pos]  + wW->data[pos]+ wE->data[pos];
			b +=	wNW->data[pos] + wNW->data[pos]  + wSW->data[pos]+ wSE->data[pos];
		}
		
		cp[i] = c/b;
		dp[i] = d/b;
	
		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos++;
			wpos++;
			epos++;
			npos++;
			spos++;
	
			a = 	-wN->data[pos];
			/*b = 	TRACE->data[pos];*/
			c = 	-wS->data[pos];
			d =	wW->data[pos]*X->data[wpos] + wE->data[pos]*X->data[epos];
			d +=	wSW->data[pos]*X->data[wpos+1] + wSE->data[pos]*X->data[epos+1];
			d +=	wNW->data[pos]*X->data[wpos-1] + wNE->data[pos]*X->data[epos-1];
			
			if(!isnan(TRACE->data[pos]))
			{
			    b = TRACE->data[pos];
			    d += B->data[pos];
			}
			else
			{
				b = 	wN->data[pos] + wS->data[pos]  + wW->data[pos]+ wE->data[pos];
				b +=	wNW->data[pos] + wNW->data[pos]  + wSW->data[pos]+ wSE->data[pos];
			}
				
			div = 1 / (b-cp[i-1]*a);
			cp[i] = c*div;
			dp[i] = (d-dp[i-1]*a)*div;
		}
	
		/*---Last row---*/
		pos++;
		wpos++;
		epos++;
		npos++;
		spos++;
	
		a = 	-wN->data[pos];
		/*b = 	TRACE->data[pos];*/
		/*c = 	0.0f;*/
		d =	wW->data[pos]*X->data[wpos] + wE->data[pos]*X->data[epos];
		d +=	wNW->data[pos]*X->data[wpos-1] + wNE->data[pos]*X->data[epos-1];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
			b = 	wN->data[pos] + wS->data[pos]  + wW->data[pos]+ wE->data[pos];
			b +=	wNW->data[pos] + wNW->data[pos]  + wSW->data[pos]+ wSE->data[pos];
		}

		dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
		temp1 = X->data[pos];
		X->data[pos] = dp[i];
	
		for(i=nrows-2;i>=0;i--)
		{
			pos = k*nrows*ncols+j*nrows+i;
			temp2 = X->data[pos];
			X->data[pos] = dp[i] - cp[i]*X->data[pos+1];

			X->data[pos+1] = omega*X->data[pos+1] + (1.0f-omega)*temp1;
			temp1 = temp2;
		}
		X->data[pos] = omega*X->data[pos] + (1.0f-omega)*temp1;

	}
	}

	
}
			
/*
//------------------------------
// TDMA horizontal line solving 
//------------------------------
*/
void TDMArow_ALR_8(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes,
			float omega)
{
	int i,j,k,pos,wpos,npos,epos,spos,swpos;
	float a,b,c,d,div,temp1,temp2;
	
	for(k=0;k<nframes;k++)
	{
	for(i=1;i<=nrows-2;i++)
	{
		/*---First column---*/
		j = 0;
		pos = k*nrows*ncols+j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		/*a =	0.0f;*/
		/*b =	TRACE->data[pos];*/
		c =	-wE->data[pos];
		d =	wS->data[pos]*X->data[spos] + wN->data[pos]*X->data[npos];
		d +=	wSE->data[pos]*X->data[spos+nrows] + wNE->data[pos]*X->data[npos+nrows];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
			b = 	wN->data[pos] + wS->data[pos]  + wW->data[pos] + wE->data[pos];
			b +=	wNW->data[pos] + wNW->data[pos]  + wSW->data[pos]+ wSE->data[pos];
		}

		cp[j] = c/b;
		dp[j] = d/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos+=nrows;
			wpos+=nrows;
			epos+=nrows;
			npos+=nrows;
			spos+=nrows;

			a = 	-wW->data[pos];	
			/*b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos];*/
			c = 	-wE->data[pos];
			d =	wS->data[pos]*X->data[spos] + wN->data[pos]*X->data[npos];
			d +=	wSW->data[pos]*X->data[spos-nrows] + wNW->data[pos]*X->data[npos-nrows];
			d +=	wSE->data[pos]*X->data[spos+nrows] + wNE->data[pos]*X->data[npos+nrows];
			
			if(!isnan(TRACE->data[pos]))
			{
			    b = TRACE->data[pos];
			    d += B->data[pos];
			}
			else
			{
				b = 	wN->data[pos] + wS->data[pos]  + wW->data[pos] + wE->data[pos];
				b +=	wNW->data[pos] + wNW->data[pos]  + wSW->data[pos]+ wSE->data[pos];
			}
			
			div = 1 / (b-cp[j-1]*a);
			cp[j] = c*div;
			dp[j] = (d-dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos+=nrows;
		wpos+=nrows;
		epos+=nrows;
		npos+=nrows;
		spos+=nrows;

		a = 	-wW->data[pos];	
		/*b = 	wN->data[pos] + wS->data[pos] + wW->data[pos];*/
		/*c = 	0.0f;*/
		d =	wS->data[pos]*X->data[spos] + wN->data[pos]*X->data[npos];
		d +=	wSW->data[pos]*X->data[spos-nrows] + wNW->data[pos]*X->data[npos-nrows];
		
		if(!isnan(TRACE->data[pos]))
		{
		    b = TRACE->data[pos];
		    d += B->data[pos];
		}
		else
		{
			b = 	wN->data[pos] + wS->data[pos]  + wW->data[pos] + wE->data[pos];
			b +=	wNW->data[pos] + wNW->data[pos]  + wSW->data[pos]+ wSE->data[pos];
		}
	
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
		temp1 = X->data[pos];
		X->data[pos] = dp[j];

		for(j=ncols-2;j>=0;j--)
		{
			pos = k*nrows*ncols+j*nrows+i;
			temp2 = X->data[pos];
			X->data[pos] = dp[j] - cp[j]*X->data[pos+nrows];

			X->data[pos+nrows] = omega*X->data[pos+nrows] + (1.0f-omega)*temp1;
			temp1 = temp2;
		}
		X->data[pos] = omega*X->data[pos] + (1.0f-omega)*temp1;
	}
	}

}

/*
//------------------------------------------------------
// TDMA vertical line solving - 4 neighbours
//------------------------------------------------------
*/

/*void TDMA_column4(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes)
{
	int i,j,k, ipos,wpos;
	float a,b,c,d,div,temp1,temp2;

	for(k=0;k<nframes;k++)
	{
		for(j=0;j<=ncols-1;j++)
		{
			/*---First row---*/
/*			i = 0;
			wpos = j*nrows+i;
			ipos = wpos+k*nrows*ncols;
		
			/*a = 	0.0f;*/
/*			b = 	2.0f + wS->data[wpos];
			c = 	-wS->data[wpos];
			d = 	Iold->data[ipos];
			
			cp[i] = c/b;
			dp[i] = d/b;
		
			/*---Middle rows---*/
/*			for(i=1;i<=nrows-2;i++)
			{
				wpos = j*nrows+i;
				ipos = wpos+k*nrows*ncols;
		
				a = 	-wN->data[wpos];
				b = 	2.0f + wN->data[wpos] + wS->data[wpos];
				c = 	-wS->data[wpos];
				d =	Iold->data[ipos];
	
				div = 1 / (b-cp[i-1]*a);
				cp[i] = c*div;
				dp[i] = (d-dp[i-1]*a)*div;
			}
		
			/*---Last row---*/
/*			wpos = j*nrows+i;
			ipos = wpos+k*nrows*ncols;
		
			a = 	-wN->data[wpos];
			b = 	2.0f + wN->data[wpos];
			/*c = 	0.0f;*/
/*			d =	Iold->data[ipos];
			
			dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
			temp1 = Inew->data[ipos];
			Inew->data[ipos] = dp[i];

			for(i=nrows-2;i>=0;i--)
			{
				wpos = j*nrows+i;
				ipos = wpos+k*nrows*ncols;
				
				temp2 = Inew->data[ipos];
				Inew->data[ipos] = dp[i] - cp[i]*Inew->data[ipos+1];
				Inew->data[ipos+1] = Inew->data[ipos+1] + temp1;
				temp1 = temp2;
				
			}
			Inew->data[ipos] = Inew->data[ipos] + temp1;
		}
	}
}


/*
//------------------------------------------------------
// TDMA horizontal line solving - 4 neighbours
//------------------------------------------------------
*/
/*void TDMA_row4(		struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes)
{
	int i,j,k,ipos,wpos;
	float a,b,c,d,div,temp1,temp2;
	
	for(k=0;k<nframes;k++)
	{
		for(i=0;i<=nrows-1;i++)
		{
			/*---First column---*/
/*			j = 0;
			wpos = j*nrows+i;
			ipos = wpos+k*nrows*ncols;
				
			/*a =	0.0f;*/
/*			b =	2.0f + wE->data[wpos];
			c =	-wE->data[wpos];
			d =	Iold->data[ipos];
	
			cp[j] = c/b;
			dp[j] = d/b;
	
			/*---Middle columns---*/
/*			for(j=1;j<=ncols-2;j++)
			{
				wpos = j*nrows+i;
				ipos = wpos+k*nrows*ncols;
	
				a = 	-wW->data[wpos];	
				b = 	2.0f + wE->data[wpos] + wW->data[wpos];
				c = 	-wE->data[wpos];
				d =	Iold->data[ipos];
				
				div = 1 / (b-cp[j-1]*a);
				cp[j] = c*div;
				dp[j] = (d-dp[j-1]*a)*div;
			}
	
			/*---Last column---*/
/*			wpos = j*nrows+i;
			ipos = wpos+k*nrows*ncols;
	
			a = 	-wW->data[wpos];	
			b = 	2.0f + wW->data[wpos];
			/*c = 	0.0f;*/
/*			d =	Iold->data[ipos];
		
			dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
			temp1 = Inew->data[ipos];
			Inew->data[ipos] = dp[j];
	
			for(j=ncols-2;j>=0;j--)
			{
				wpos = j*nrows+i;
				ipos = wpos+k*nrows*ncols;

				temp2 = Inew->data[ipos];
				Inew->data[ipos] = dp[j] - cp[j]*Inew->data[ipos+nrows];
				Inew->data[ipos+nrows] = Inew->data[ipos+nrows] + temp1;
				temp1 = temp2;
			}
			Inew->data[ipos] = Inew->data[ipos] + temp1;
		}
	}
}


/*
//------------------------------------------------------
// TDMA vertical line solving - 8 neighbours
// !!!NOT IMPLEMENTED!!!
//------------------------------------------------------
*/
void TDMA_column8(	struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes)
{
	int i,j,k, ipos,wpos;
	float a,b,c,d,div,temp1,temp2;

	
}


/*
//------------------------------------------------------
// TDMA horizontal line solving - 8 neighbours
// !!!NOT IMPLEMENTED!!!
//------------------------------------------------------
*/
void TDMA_row8(		struct matrixM *X,
			struct matrixM *TRACE,
			struct matrixM *B,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			float *cp,
			float *dp,
			int ncols,
			int nrows,
			int nframes)
{
	int i,j,k,ipos,wpos;
	float a,b,c,d,div,temp1,temp2;
	
	
}

