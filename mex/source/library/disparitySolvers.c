/*
Functions related to variational disparity calcuation

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

#if !defined( isnan )
  #define isnan(x)((x)!=(x))
#endif

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin4_2d(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam)
{
	float wNeigh, A, B;
	float *data_div, *dividend;
	unsigned int i,j,pos,iter,wpos,npos,epos,spos;
	unsigned int nrows = Cu->dimElems[0];
	unsigned int ncols = Cu->dimElems[1];
	unsigned int iterations = (int)Mparam.iter;

	if( (data_div = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin4_2d: error reserving memory for 'data_div'\n");
		return;
	}

	if( (dividend = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin4_2d: error reserving memory for 'data_div'\n");
		free(data_div);
		return;
	}

	for(iter=0;iter<iterations;iter++)
	{
		for(j=1;j<ncols-1;j++)
		{
		for(i=1;i<nrows-1;i++)		
		{
			/*Spatial and temporal position indices*/
			pos = j*nrows+i;
			wpos = pos-nrows;
			epos = pos+nrows;
			npos = pos-1;
			spos = pos+1;

			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			wNeigh =	(U->data[epos] + dU->data[epos] - U->data[pos])*wE->data[ pos ] +
					(U->data[wpos] + dU->data[wpos] - U->data[pos])*wW->data[ pos ] + 
					(U->data[spos] + dU->data[spos] - U->data[pos])*wS->data[ pos ] +
					(U->data[npos] + dU->data[npos] - U->data[pos])*wN->data[ pos ];

			if(iter==0)
			{
				if(!isnan(Cu->data[pos]))
				{
					dividend[pos] = Cu->data[pos];
					data_div[pos] = 1.0f/(	Du->data[pos]+ 
								wE->data[ pos ]+
								wW->data[ pos ]+
								wS->data[ pos ]+
								wN->data[ pos ] );
				}
				else
				{	
					dividend[pos] = 0.0f;
					data_div[pos] = 1.0f/(	wE->data[ pos ]+
								wW->data[ pos ]+
								wS->data[ pos ]+
								wN->data[ pos ] );
				}
			}

			/*dU_newapprox = (wNeigh+dividend[pos])*data_div[pos];*/
			A = (1.0f-Mparam.omega)*dU->data[ pos ];
			B = Mparam.omega*(wNeigh+dividend[pos])*data_div[pos];
			dU->data[ pos ] = A + B;
		}
		}

		/*
		-------------------
		--Fill the borders
		-------------------
		*/
		for(j=0;j<ncols;j++)
		{
			pos = j*nrows;

			dU->data[ pos ] = dU->data[ pos+1 ];
			dU->data[ pos+nrows-1 ] = dU->data[ pos+nrows-2 ];
		}

		for(i=0;i<nrows;i++)
		{
			dU->data[ i ] = dU->data[ i+nrows ];
			dU->data[ i+(ncols-1)*nrows ] = dU->data[ i+(ncols-2)*nrows ];
		}
	}
	
	free(data_div);
	free(dividend);
}

/*
//---------------------------------------------------------------------------------------
//---Alternating line relaxation Gauss-Seidel (Gauss-Seidel block variant)
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization.
//---Late linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void	GS_ALR_SOR_llin4_2d(	struct matrixM *U,
				struct matrixM *dU,
				struct matrixM *Cu,
				struct matrixM *Du,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam)
{
	float *Cp, *Dp;
	int nrows = Cu->dimElems[0];
	int ncols = Cu->dimElems[1];
	int maxLength;
	int iterations = (int)Mparam.iter;
	int iter;

	maxLength = (nrows>ncols)?nrows:ncols;

	/*Reserve space for Cprime and Dprime temp variables used in TDMA*/
	if( (Cp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("GS_ALR_SOR_llin4_2d: error reserving space for 'Cp'\n");
		return;
	}
	if( (Dp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("GS_ALR_SOR_llin4_2d: error reserving space for 'Dp'\n");
		free(Cp);
		return;
	}

	for(iter=0;iter<iterations;iter++)
	{
		/*
		------------------------
		-Vertical line relaxing
		------------------------
		*/
		westColumn4(	U, dU, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleColumn4(	U, dU, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		eastColumn4(	U, dU, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		northRow4(	U, dU, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleRow4(	U, dU, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		southRow4(	U, dU, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);

	}

	free( Cp );
	free( Dp );

}

/*
//---------------------------------------------------------------------------------------------------
// Calculate residuals r (Ax=b => r=b-Ax) for non-linear (late linearization) case with 4-neighbours
//---------------------------------------------------------------------------------------------------
*/
void Residuals_llin4_2d(	struct matrixM *RU,
				struct matrixM *U,
				struct matrixM *dU,
				struct matrixM *Cu,
				struct matrixM *Du,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam)
{
	float wNeigh;
	unsigned int i,j,pos,iter,wpos,npos,epos,spos;
	unsigned int nrows = Cu->dimElems[0];
	unsigned int ncols = Cu->dimElems[1];

	for(j=1;j<ncols-1;j++)
	{
	for(i=1;i<nrows-1;i++)		
	{
		/*Spatial and temporal position indices*/
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		/*
		---------------------
		--Spatial diffusion--
		---------------------	
		*/
		wNeigh =	(U->data[epos] + dU->data[epos] - U->data[pos])*wE->data[ pos ] +
				(U->data[wpos] + dU->data[wpos] - U->data[pos])*wW->data[ pos ] + 
				(U->data[spos] + dU->data[spos] - U->data[pos])*wS->data[ pos ] +
				(U->data[npos] + dU->data[npos] - U->data[pos])*wN->data[ pos ];

		/*We don't like NaN:s, do we*/
		if( !isnan(Cu->data[pos]) )
		{
			RU->data[pos] = Cu->data[pos] + wNeigh - dU->data[pos]*(	Du->data[pos] + 
											wW->data[pos] + 
											wN->data[pos] + 
											wE->data[pos] +
											wS->data[pos] );
		}else{
		
			RU->data[pos] = wNeigh - dU->data[pos]*(	wW->data[pos] + 
									wN->data[pos] + 
									wE->data[pos] +
									wS->data[pos] );
		}

	}
	}

	/*
	-------------------
	--Fill the borders
	-------------------
	*/
	for(j=0;j<ncols;j++)
	{
		pos = j*nrows;

		RU->data[ pos ] = RU->data[ pos+1 ];
		RU->data[ pos+nrows-1 ] = RU->data[ pos+nrows-2 ];
	}

	for(i=0;i<nrows;i++)
	{
		RU->data[ i ] = RU->data[ i+nrows ];
		RU->data[ i+(ncols-1)*nrows ] = RU->data[ i+(ncols-2)*nrows ];
	}

}

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions with symmetry.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llinsym4_2d(	struct matrixM *U0, 
				struct matrixM *dU0,
				struct matrixM *Cu0,
				struct matrixM *Du0,
				struct matrixM *wW0,
				struct matrixM *wN0,
				struct matrixM *wE0,
				struct matrixM *wS0,
				struct matrixM *U1, 
				struct matrixM *dU1,
				struct matrixM *Cu1,
				struct matrixM *Du1,
				struct matrixM *wW1,
				struct matrixM *wN1,
				struct matrixM *wE1,
				struct matrixM *wS1,
				struct mparam Mparam)
{
	float wNeigh0, wNeigh1;
	float dU_newapprox0, dU_newapprox1;
	float *u_data_div0, *u_data_div1, *dividend0, *dividend1;
	unsigned int i,j,pos,wpos,npos,epos,spos,iter;
	unsigned int nrows = Cu0->dimElems[0];
	unsigned int ncols = Cu0->dimElems[1];
	unsigned int channel_offset = nrows*ncols;

	/*Reserve space for the divisor terms*/
	if( (u_data_div0 = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym4_2d: error reserving space for 'u_data_div0'\n");
		return;
	}
	if( (u_data_div1 = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym4_2d: error reserving space for 'u_data_div1'\n");
		free(u_data_div0);
		return;
	}
	if( (dividend0 = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym4_2d: error reserving space for 'dividend0'\n");
		free(u_data_div0);
		free(u_data_div1);
		return;
	}
	if( (dividend1 = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym4_2d: error reserving space for 'dividend1'\n");
		free(u_data_div0);
		free(u_data_div1);
		free(dividend0);
		return;
	}

	for(iter=0;iter<(int)Mparam.iter;iter++)
	{
		for(j=1;j<ncols-1;j++)
		{
		for(i=1;i<nrows-1;i++)		
		{
			/*Spatial position indices*/
			pos = j*nrows+i;
			wpos = pos - nrows;
			epos = pos + nrows;
			npos = pos - 1;
			spos = pos + 1;

			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			wNeigh0 =	(U0->data[epos]+dU0->data[epos]-U0->data[pos])*wE0->data[ pos ] +
					(U0->data[wpos]+dU0->data[wpos]-U0->data[pos])*wW0->data[ pos ] + 
					(U0->data[spos]+dU0->data[spos]-U0->data[pos])*wS0->data[ pos ] +
					(U0->data[npos]+dU0->data[npos]-U0->data[pos])*wN0->data[ pos ];

			wNeigh1 =	(U1->data[epos]+dU1->data[epos]-U1->data[pos])*wE1->data[ pos ] +
					(U1->data[wpos]+dU1->data[wpos]-U1->data[pos])*wW1->data[ pos ] + 
					(U1->data[spos]+dU1->data[spos]-U1->data[pos])*wS1->data[ pos ] +
					(U1->data[npos]+dU1->data[npos]-U1->data[pos])*wN1->data[ pos ];

			/*Calculate those elements only once that don't change*/
			if(iter==0)
			{
				if(!isnan(Cu0->data[pos]))
				{
					dividend0[pos] = Cu0->data[pos];
					u_data_div0[pos] = 1.0f/(	Du0->data[pos]+ 
									wE0->data[ pos ]+
									wW0->data[ pos ]+
									wS0->data[ pos ]+
									wN0->data[ pos ]);
				}
				else
				{	u_data_div0[pos] = 1.0f/(	wE0->data[ pos ]+
									wW0->data[ pos ]+
									wS0->data[ pos ]+
									wN0->data[ pos ]);
				}
				if(!isnan(Cu1->data[pos]))
				{
					dividend1[pos] = Cu1->data[pos];
					u_data_div1[pos] = 1.0f/(	Du1->data[pos]+ 
									wE1->data[ pos ]+
									wW1->data[ pos ]+
									wS1->data[ pos ]+
									wN1->data[ pos ]);
				}
				else
				{	u_data_div1[pos] = 1.0f/(	wE1->data[ pos ]+
									wW1->data[ pos ]+
									wS1->data[ pos ]+
									wN1->data[ pos ]);
				}
			}

			dU_newapprox0 = (wNeigh0+dividend0[pos])*u_data_div0[pos];
			dU0->data[ pos ] = (1.0f-Mparam.omega)*dU0->data[ pos ] + Mparam.omega*dU_newapprox0;

			dU_newapprox1 = (wNeigh1+dividend1[pos])*u_data_div1[pos];
			dU1->data[ pos ] = (1.0f-Mparam.omega)*dU1->data[ pos ] + Mparam.omega*dU_newapprox1;
		}
		}

		/*
		-------------------
		--Fill the borders
		-------------------
		*/
		for(j=0;j<ncols;j++)
		{
			pos = j*nrows;

			dU0->data[ pos ] = dU0->data[ pos+1 ];
			dU0->data[ pos+nrows-1 ] = dU0->data[ pos+nrows-2 ];
			dU1->data[ pos ] = dU1->data[ pos+1 ];
			dU1->data[ pos+nrows-1 ] = dU1->data[ pos+nrows-2 ];
		}

		for(i=0;i<nrows;i++)
		{
			dU0->data[ i ] = dU0->data[ i+nrows ];
			dU0->data[ i+(ncols-1)*nrows ] = dU0->data[ i+(ncols-2)*nrows ];
			dU1->data[ i ] = dU1->data[ i+nrows ];
			dU1->data[ i+(ncols-1)*nrows ] = dU1->data[ i+(ncols-2)*nrows ];
		}
	}

	free(u_data_div0);
	free(u_data_div1);
	free(dividend0);
	free(dividend1);
}

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions with symmetry.
//---------------------------------------------------------------------------------------
*/
void GS_ALR_SOR_llinsym4_2d(	struct matrixM *U0, 
				struct matrixM *dU0,
				struct matrixM *Cu0,
				struct matrixM *Du0,
				struct matrixM *wW0,
				struct matrixM *wN0,
				struct matrixM *wE0,
				struct matrixM *wS0,
				struct matrixM *U1, 
				struct matrixM *dU1,
				struct matrixM *Cu1,
				struct matrixM *Du1,
				struct matrixM *wW1,
				struct matrixM *wN1,
				struct matrixM *wE1,
				struct matrixM *wS1,
				struct mparam Mparam)
{
	float *Cp, *Dp;
	unsigned int nrows = Cu0->dimElems[0];
	unsigned int ncols = Cu0->dimElems[1];
	unsigned int maxLength;
	int iterations = (int)Mparam.iter;
	int iter;

	maxLength = (nrows>ncols)?nrows:ncols;

	/*Reserve space for Cprime and Dprime temp variables used in TDMA*/
	if( (Cp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("GS_ALR_SOR_llinsym4_2d: error reserving space for 'Cp'\n");
		return;
	}
	if( (Dp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("GS_ALR_SOR_llinsym4_2d: error reserving space for 'Dp'\n");
		free(Dp);
		return;
	}


	for(iter=0;iter<iterations;iter++)
	{
		/*
		------------------------
		-Vertical line relaxing
		------------------------
		*/
		westColumn4(	U0, dU0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		middleColumn4(	U0, dU0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		eastColumn4(	U0, dU0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		northRow4(	U0, dU0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		middleRow4(	U0, dU0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		southRow4(	U0, dU0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);

		/*
		------------------------
		-Vertical line relaxing
		------------------------
		*/
		westColumn4(	U1, dU1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		middleColumn4(	U1, dU1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		eastColumn4(	U1, dU1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		northRow4(	U1, dU1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		middleRow4(	U1, dU1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		southRow4(	U1, dU1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);

	}

	free( Cp );
	free( Dp );
}

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin8_2d(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
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
	float wNeighHor, wNeighVer, sumHor, sumVer, dU_newapprox;
	float *data_div, *dividend, neighbours;
	unsigned int i,j,k,pos,wpos,nwpos,npos,nepos,epos,sepos,spos,swpos,iter;
	unsigned int nrows = Cu->dimElems[0];
	unsigned int ncols = Cu->dimElems[1];
	unsigned int channel_offset = nrows*ncols;
	unsigned int nchannels = 1;

	if(Cu->ndims>2)
	nchannels = Cu->dimElems[2];

	if( (data_div = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin8_2d: error reserving memory for 'data_div'\n");
		return;
	}

	if( (dividend = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin8_2d: error reserving memory for 'data_div'\n");
		free(data_div);
		return;
	}

	for(iter=0;iter<(int)Mparam.iter;iter++)
	{
		for(j=1;j<ncols-1;j++)
		{
		for(i=1;i<nrows-1;i++)		
		{
			/*Spatial and temporal position indices*/
			pos = j*nrows+i;

			npos = pos-1;
			spos = pos+1;
			wpos = pos-nrows;
			epos = pos+nrows;

			nwpos = wpos-1;
			nepos = epos-1;
			swpos = wpos+1;
			sepos = epos+1;
			
			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			neighbours = (U->data[wpos]+dU->data[wpos]-U->data[pos])*wW->data[ pos ];
			neighbours += (U->data[nwpos]+dU->data[nwpos]-U->data[pos])*wNW->data[ pos ];
			neighbours += (U->data[npos]+dU->data[npos]-U->data[pos])*wN->data[ pos ];
			neighbours += (U->data[nepos]+dU->data[nepos]-U->data[pos])*wNE->data[ pos ];
			neighbours += (U->data[epos]+dU->data[epos]-U->data[pos])*wE->data[ pos ];
			neighbours += (U->data[sepos]+dU->data[sepos]-U->data[pos])*wSE->data[ pos ];
			neighbours += (U->data[spos]+dU->data[spos]-U->data[pos])*wS->data[ pos ];
			neighbours += (U->data[swpos]+dU->data[swpos]-U->data[pos])*wSW->data[ pos ];

			/*Calculate those elements only once that don't change*/
			if(iter==0)
			{
				data_div[pos] = wW->data[ pos ]+wNW->data[ pos ]+wN->data[ pos ]+wNE->data[ pos ]
						+wE->data[ pos ]+wSE->data[ pos ]+wS->data[ pos ]+wSW->data[ pos ];

				if(!isnan(Cu->data[pos]))
				{
					dividend[pos] = dividend[pos] + Cu->data[pos];
					data_div[pos] = data_div[pos] + Du->data[pos];
				}
			}

			dU_newapprox = (neighbours+dividend[pos])/(data_div[pos]);
			if(!isnan(dU_newapprox))
				dU->data[ pos ] = (1.0f-Mparam.omega)*dU->data[ pos ] + Mparam.omega*dU_newapprox;
		}
		}

		/*
		-------------------
		--Fill the borders
		-------------------
		*/
		for(j=0;j<ncols;j++)
		{
			pos = j*nrows;

			dU->data[ pos ] = dU->data[ pos+1 ];
			dU->data[ pos+nrows-1 ] = dU->data[ pos+nrows-2 ];

		}

		for(i=0;i<nrows;i++)
		{
			dU->data[ i ] = dU->data[ i+nrows ];
			dU->data[ i+(ncols-1)*nrows ] = dU->data[ i+(ncols-2)*nrows ];

			
		}
		
			
	}

	free(data_div);
	free(dividend);
}

/*
//---------------------------------------------------------------------------------------
//---Alternating line relaxation Gauss-Seidel (Gauss-Seidel block variant)
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization.
//---Late linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void	GS_ALR_SOR_llin8_2d(	struct matrixM *U,
				struct matrixM *dU,
				struct matrixM *Cu,
				struct matrixM *Du,
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
	int nrows = Cu->dimElems[0];
	int ncols = Cu->dimElems[1];
	int maxLength;
	int iterations = (int)Mparam.iter;
	int iter;

	maxLength = (nrows>ncols)?nrows:ncols;

	/*Reserve space for Cprime and Dprime temp variables used in TDMA*/
	if( (Cp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("GS_ALR_SOR_llin8_2d: error reserving space for 'Cp'\n");
		return;
	}
	if( (Dp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("GS_ALR_SOR_llin8_2d: error reserving space for 'Dp'\n");
		free(Dp);
		return;
	}

	for(iter=0;iter<iterations;iter++)
	{
		/*
		------------------------
		-Vertical line relaxing
		------------------------
		*/
		westColumn8(	U, dU, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		middleColumn8(	U, dU, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		eastColumn8(	U, dU, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		northRow8(	U, dU, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		middleRow8(	U, dU, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		southRow8(	U, dU, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);

	}

	free( Cp );
	free( Dp );
}

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Late lizearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llinsym8_2d(	struct matrixM *U0, 
				struct matrixM *dU0,
				struct matrixM *Cu0,
				struct matrixM *Du0,
				struct matrixM *wW0,
				struct matrixM *wNW0,
				struct matrixM *wN0,
				struct matrixM *wNE0,
				struct matrixM *wE0,
				struct matrixM *wSE0,
				struct matrixM *wS0,
				struct matrixM *wSW0,
				struct matrixM *U1, 
				struct matrixM *dU1,
				struct matrixM *Cu1,
				struct matrixM *Du1,
				struct matrixM *wW1,
				struct matrixM *wNW1,
				struct matrixM *wN1,
				struct matrixM *wNE1,
				struct matrixM *wE1,
				struct matrixM *wSE1,
				struct matrixM *wS1,
				struct matrixM *wSW1,
				struct mparam Mparam)
{
	float dU_newapprox0, dU_newapprox1, neighbours0, neighbours1;
	float *u_data_div0, *u_data_div1, *dividend0, *dividend1;
	unsigned int i,j,k,pos,wpos,nwpos,npos,nepos,epos,sepos,spos,swpos,iter;
	unsigned int nrows = Cu0->dimElems[0];
	unsigned int ncols = Cu0->dimElems[1];
	unsigned int nchannels = 1;
	unsigned int channel_offset = nrows*ncols;

	if( Cu0->ndims>2)
		nchannels = Cu0->dimElems[2];
	

	/*Reserve space for the divisor terms*/
	if( (u_data_div0 = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym8_2d: error reserving space for 'u_data_div0'\n");
		return;
	}
	if( (u_data_div1 = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym8_2d: error reserving space for 'u_data_div1'\n");
		free(u_data_div0);
		return;
	}
	if( (dividend0 = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym8_2d: error reserving space for 'dividend0'\n");
		free(u_data_div0);
		free(u_data_div1);
		return;
	}
	if( (dividend1 = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym8_2d: error reserving space for 'dividend1'\n");
		free(u_data_div0);
		free(u_data_div1);
		free(dividend1);
		return;
	}

	for(iter=0;iter<(int)Mparam.iter;iter++)
	{
		for(j=1;j<ncols-1;j++)
		{
		for(i=1;i<nrows-1;i++)		
		{
			/*Spatial and temporal position indices*/
			pos = j*nrows+i;

			npos = pos-1;
			spos = pos+1;
			wpos = pos-nrows;
			epos = pos+nrows;

			nwpos = wpos-1;
			nepos = epos-1;
			swpos = wpos+1;
			sepos = epos+1;

			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			neighbours0 = (U0->data[wpos]+dU0->data[wpos]-U0->data[pos])*wW0->data[ pos ];
			neighbours0 += (U0->data[nwpos]+dU0->data[nwpos]-U0->data[pos])*wNW0->data[ pos ];
			neighbours0 += (U0->data[npos]+dU0->data[npos]-U0->data[pos])*wN0->data[ pos ];
			neighbours0 += (U0->data[nepos]+dU0->data[nepos]-U0->data[pos])*wNE0->data[ pos ];
			neighbours0 += (U0->data[epos]+dU0->data[epos]-U0->data[pos])*wE0->data[ pos ];
			neighbours0 += (U0->data[sepos]+dU0->data[sepos]-U0->data[pos])*wSE0->data[ pos ];
			neighbours0 += (U0->data[spos]+dU0->data[spos]-U0->data[pos])*wS0->data[ pos ];
			neighbours0 += (U0->data[swpos]+dU0->data[swpos]-U0->data[pos])*wSW0->data[ pos ];

			neighbours1 = (U1->data[wpos]+dU1->data[wpos]-U1->data[pos])*wW1->data[ pos ];
			neighbours1 += (U1->data[nwpos]+dU1->data[nwpos]-U1->data[pos])*wNW1->data[ pos ];
			neighbours1 += (U1->data[npos]+dU1->data[npos]-U1->data[pos])*wN1->data[ pos ];
			neighbours1 += (U1->data[nepos]+dU1->data[nepos]-U1->data[pos])*wNE1->data[ pos ];
			neighbours1 += (U1->data[epos]+dU1->data[epos]-U1->data[pos])*wE1->data[ pos ];
			neighbours1 += (U1->data[sepos]+dU1->data[sepos]-U1->data[pos])*wSE1->data[ pos ];
			neighbours1 += (U1->data[spos]+dU1->data[spos]-U1->data[pos])*wS1->data[ pos ];
			neighbours1 += (U1->data[swpos]+dU1->data[swpos]-U1->data[pos])*wSW1->data[ pos ];

			/*Calculate those elements only once that don't change*/
			if(iter==0)
			{
				u_data_div0[pos] = wW0->data[ pos ]+wNW0->data[ pos ]+wN0->data[ pos ]+wNE0->data[ pos ]
						+wE0->data[ pos ]+wSE0->data[ pos ]+wS0->data[ pos ]+wSW0->data[ pos ];

				u_data_div1[pos] = wW1->data[ pos ]+wNW1->data[ pos ]+wN1->data[ pos ]+wNE1->data[ pos ]
						+wE1->data[ pos ]+wSE1->data[ pos ]+wS1->data[ pos ]+wSW1->data[ pos ];

				if(!isnan(Cu0->data[pos]))
				{
					dividend0[pos] = dividend0[pos] + Cu0->data[pos];
					u_data_div0[pos] = u_data_div0[pos] + Du0->data[pos];
				}

				if(!isnan(Cu1->data[pos]))
				{
					dividend1[pos] = dividend1[pos] + Cu1->data[pos];
					u_data_div1[pos] = u_data_div1[pos] + Du1->data[pos];
				}
			}

			dU_newapprox0 = (neighbours0+dividend0[pos])/(u_data_div0[pos]);
			dU_newapprox1 = (neighbours1+dividend1[pos])/(u_data_div1[pos]);

			dU0->data[ pos ] = (1.0f-Mparam.omega)*dU0->data[ pos ] + Mparam.omega*dU_newapprox0;
			dU1->data[ pos ] = (1.0f-Mparam.omega)*dU1->data[ pos ] + Mparam.omega*dU_newapprox1;
			
		}
		}

		/*
		-------------------
		--Fill the borders
		-------------------
		*/
		for(j=0;j<ncols;j++)
		{
			pos = j*nrows;

			dU0->data[ pos ] = dU0->data[ pos+1 ];
			dU0->data[ pos+nrows-1 ] = dU0->data[ pos+nrows-2 ];
			dU1->data[ pos ] = dU1->data[ pos+1 ];
			dU1->data[ pos+nrows-1 ] = dU1->data[ pos+nrows-2 ];
		}

		for(i=0;i<nrows;i++)
		{
			dU0->data[ i ] = dU0->data[ i+nrows ];
			dU0->data[ i+(ncols-1)*nrows ] = dU0->data[ i+(ncols-2)*nrows ];
			dU1->data[ i ] = dU1->data[ i+nrows ];
			dU1->data[ i+(ncols-1)*nrows ] = dU1->data[ i+(ncols-2)*nrows ];
		}
	}

	free(u_data_div0);
	free(u_data_div1);
	free(dividend0);
	free(dividend1);
}

/*
//---------------------------------------------------------------------------------------
//---Succesive over relaxation using both brightness and gradient constraints
//---and 2D first and second order regularization. Non-lizearized constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin4_2d12(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW1,
			struct matrixM *wN1,
			struct matrixM *wE1,
			struct matrixM *wS1,
			struct matrixM *wW2,
			struct matrixM *wN2,
			struct matrixM *wE2,
			struct matrixM *wS2,
			struct mparam Mparam)
{
	float wNeigh1, wNeigh2, dU_newapprox;
	float *data_div, *dividend;
	unsigned int i,j,iter,pos,wpos,npos,epos,spos, wpos2,npos2,epos2,spos2;
	unsigned int nrows = Cu->dimElems[0];
	unsigned int ncols = Cu->dimElems[1];
	unsigned int channel_offset = nrows*ncols;

	if( (data_div = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin4_2d12: error reserving memory for 'data_div'\n");
		return;
	}
	if( (dividend = (float*)calloc( nrows*ncols, sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin4_2d12: error reserving memory for 'dividend'\n");
		free(data_div);
		return;
	}

	for(iter=0;iter<(int)Mparam.iter;iter++)
	{
		
		for(j=2;j<ncols-2;j++)
		{
		for(i=2;i<nrows-2;i++)		
		{
			/*Spatial and temporal position indices*/
			pos = j*nrows+i;
			wpos = pos-nrows;
			epos = pos+nrows;
			npos = pos-1;
			spos = pos+1;

			wpos2 = pos-2*nrows;
			epos2 = pos+2*nrows;
			npos2 = pos-2;
			spos2 = pos+2;

			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			/*First order regularization weights*/
			wNeigh1 =	(U->data[epos]+dU->data[epos]-U->data[pos]) * wE1->data[ pos ] +
					(U->data[wpos]+dU->data[wpos]-U->data[pos]) * wW1->data[ pos ] + 
					(U->data[spos]+dU->data[spos]-U->data[pos]) * wS1->data[ pos ] +
					(U->data[npos]+dU->data[npos]-U->data[pos]) * wN1->data[ pos ];
			/*Second order regularization weights*/
			wNeigh2 =	(-U->data[wpos]+3*U->data[pos]-3*U->data[epos]+U->data[epos2] +
					-dU->data[wpos]-3*dU->data[epos]+dU->data[epos2]) * wE2->data[pos] +
					
					(U->data[wpos2]-3*U->data[wpos]+3*U->data[pos]-U->data[epos] +
					+dU->data[wpos2]-3*dU->data[wpos]-dU->data[epos]) * wW2->data[pos] +
					
					(-U->data[npos]+3*U->data[pos]-3*U->data[spos]+U->data[spos2] +
					-dU->data[npos]-3*dU->data[spos]+dU->data[spos2]) * wS2->data[pos] +
					
					(U->data[npos2]-3*U->data[npos]+3*U->data[pos]-U->data[spos] +
					+dU->data[npos2]-3*dU->data[npos]-dU->data[spos]) * wN2->data[pos];

			if(iter==0)
			{
				if(!isnan(Cu->data[pos]))
				{
					dividend[pos] = Cu->data[pos];
					data_div[pos] = 1.0f/(	Du->data[pos]+ 
								wE1->data[ pos ]+
								wW1->data[ pos ]+
								wS1->data[ pos ]+
								wN1->data[ pos ]+
								3*wE2->data[ pos ]+
								3*wW2->data[ pos ]+
								3*wS2->data[ pos ]+
								3*wN2->data[ pos ]);
				}
				else
				{	data_div[pos] = 1.0f/(	wE1->data[ pos ]+
								wW1->data[ pos ]+
								wS1->data[ pos ]+
								wN1->data[ pos ]+
								3*wE2->data[ pos ]+
								3*wW2->data[ pos ]+
								3*wS2->data[ pos ]+
								3*wN2->data[ pos ]);
				}
			}

			dU_newapprox = (wNeigh1+wNeigh2+dividend[pos])*data_div[pos];
			dU->data[ pos ] = (1.0f-Mparam.omega)*dU->data[ pos ] + Mparam.omega*dU_newapprox;
		}
		}

		/*
		-------------------
		--Fill the borders
		-------------------
		*/
		for(j=0;j<ncols;j++)
		{
			pos = j*nrows;

			dU->data[ pos ] = dU->data[ pos+2 ];
			dU->data[ pos+1 ] = dU->data[ pos+2 ];
			dU->data[ pos+nrows-1 ] = dU->data[ pos+nrows-3 ];
			dU->data[ pos+nrows-2 ] = dU->data[ pos+nrows-3 ];
		}

		for(i=0;i<nrows;i++)
		{
			dU->data[ i ] = dU->data[ i+2*nrows ];
			dU->data[ i+nrows ] = dU->data[ i+2*nrows ];
			dU->data[ i+(ncols-1)*nrows ] = dU->data[ i+(ncols-3)*nrows ];
			dU->data[ i+(ncols-2)*nrows ] = dU->data[ i+(ncols-3)*nrows ];
		}
			
	}

	free(data_div);
	free(dividend);
}

/*
//------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation.
//---2D (spatial) regularization (4 neighbours).
//---Non-lizearized constancy assumptions.
//---Uses inline assembler (AT&T syntax for Linux's GAS).
//------------------------------------------------------
*/
/*TODO: this should probably be removed completely*/
void GS_SOR_llin4ASM_2d(	float *U,
				float *dU,
				float *Cu,
				float *Du,
				float *wW,
				float *wN,
				float *wE,
				float *wS,
				unsigned int nrows,
				unsigned int ncols,
				unsigned int iterations,
				float omega)
{

		float wNeigh, A,B;
		float *data_div, *dividend;
		float *data_div_ptr, *dividend_ptr;
		int i,j,pos,iter;
		unsigned short int fstcw, fstsw;
		
		/*Is this the only way of passing the pointers to inline GAS???*/
		data_div_ptr = &data_div[0];
		dividend_ptr = &dividend[0];
		
#if defined(GASASM)

	if( (data_div = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin4ASM_2d: error reserving memory for 'data_div'\n");
		return;
	}
	if( (dividend = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin4ASM_2d: error reserving memory for 'data_div'\n");
		free(data_div);
		return;
	}

	/*Store the FPU control word and initialize FPU for calculations*/
	__asm__ ( 	
			"fstcw %0\n\t"
			"fstsw %1\n\t"
			"finit\n\t"
			:"=m"(fstcw),"=m"(fstsw)
		);

	for(iter=0;iter<iterations;iter++)
	{

		for(j=1;j<ncols-1;j++)
		{
		for(i=1;i<nrows-1;i++)		
		{
			/*Spatial and temporal position indices*/
			pos = j*nrows+i;
			/*wpos = pos-nrows;
			epos = pos+nrows;
			npos = pos-1;
			spos = pos+1;*/

			__asm (		
					/*
					---------------------
					--Spatial diffusion--
					---------------------	
			
					wNeigh =	(U->data[pos+nrows] + dU->data[pos+nrows] - U->data[pos])*wE->data[ pos ] +
							(U->data[pos-nrows] + dU->data[pos-nrows] - U->data[pos])*wW->data[ pos ] + 
							(U->data[pos+1] + dU->data[pos+1] - U->data[pos])*wS->data[ pos ] +
							(U->data[pos-1] + dU->data[pos-1] - U->data[pos])*wN->data[ pos ];
					*/
					/*Calculate pointers to U[pos+x]*/
					"movl %[U], 			%%eax\n\t"
					"movl %[nrows],			%%ecx\n\t"
					"leal (%%eax,%[pos],4),		%%eax\n\t"
					"movl %%eax,			%%edx\n\t"
					"neg %%ecx\n\t"
					"leal (%%eax,%%ecx,4),		%%ecx\n\t"
					"leal (%%eax,%[nrows],4),	%%eax\n\t"
					/*Load values into FPU*/
					"flds (%%eax)\n\t"
					"flds (%%ecx)\n\t"
					"flds 4(%%edx)\n\t"
					"flds -4(%%edx)\n\t"

					/*Calculate pointers to dU[pos+x]*/
					"movl %[dU],			%%eax\n\t"
					"movl %[nrows],			%%ecx\n\t"
					"leal (%%eax,%[pos],4),		%%eax\n\t"
					"movl %%eax,			%%edx\n\t"
					"neg %%ecx\n\t"
					"leal (%%eax,%%ecx,4),		%%ecx\n\t"
					"leal (%%eax,%[nrows],4),	%%eax\n\t"
					/*Load values into FPU*/
					"flds (%%eax)\n\t"
					"flds (%%ecx)\n\t"
					"flds 4(%%edx)\n\t"
					"flds -4(%%edx)\n\t"
					/*U[pos+x]+dU[pos+x]*/
					"faddp %%st(4)\n\t"
					"faddp %%st(4)\n\t"
					"faddp %%st(4)\n\t"
					"faddp %%st(4)\n\t"
							
					/*Load U[pos] into FPU*/
					"movl %[U],			%%eax\n\t"
					"leal (%%eax,%[pos],4),		%%eax\n\t"
					"flds (%%eax)\n\t"
					
					/*(U[pos+x]+dU[pos+x])-U[pos]*/
					"fsubr %%st,%%st(1)\n\t"
					"fsubr %%st,%%st(2)\n\t"
					"fsubr %%st,%%st(3)\n\t"
					"fsubrp %%st,%%st(4)\n\t"

					/*Calculate pointers to diffusion weights*/
					"movl %[wE],			%%eax\n\t"
					"movl %[wW],			%%ecx\n\t"
					"movl %[wS],			%%edx\n\t"
					"movl %[wN],			%%edi\n\t"
					/*Load values into FPU*/
					"flds (%%eax,%[pos],4)\n\t"
					"flds (%%ecx,%[pos],4)\n\t"
					"flds (%%edx,%[pos],4)\n\t"
					"flds (%%edi,%[pos],4)\n\t"

					/*Multiply by diffusion weights*/
					"fmulp %%st(4)\n\t"
					"fmulp %%st(4)\n\t"
					"fmulp %%st(4)\n\t"
					"fmulp %%st(4)\n\t"
					/*Sum the result and store*/
					"faddp %%st(1)\n\t"
					"faddp %%st(1)\n\t"
					"faddp %%st(1)\n\t"			/*Leave wNeigh in the stack*/
					/*"fstps %[wNeigh]\n\t"*/

					/*if(iter==0)
					{
						if(!isnan(Cu[pos]))
						{
							dividend[pos] = Cu[pos];
							data_div[pos] = 1.0f/(	Du[pos]+ 
										wE[ pos ]+
										wW[ pos ]+
										wS[ pos ]+
										wN[ pos ] );
						else
						{	/*dividend[pos] = 0.0f;
							data_div[pos] = 1.0f/(	wE[ pos ]+
										wW[ pos ]+
										wS[ pos ]+
										wN[ pos ] );
						}
					}*/
					"movl %[iter], %%eax\n\t"
					"mov $0, %%ecx\n\t"
					"cmp %%eax, %%ecx\n\t"
					"jne cont\n\t"

					"movl %[Cu], %%eax\n\t"
					"flds (%%eax,%[pos],4)\n\t"
					"fucom %%st(0)\n\t"			/*Compare and leave Cu in the stack*/
					"fnstsw %%ax\n\t"
					"sahf\n\t"
					"jp nan\n\t"
					"jmp notnan\n\t"

					"nan:\n\t"
						"fucomp %%st(0)\n\t"		/*Empty Cu from the stack*/
						"movl %[dividend],	%%eax\n\t"
						"fldz \n\t"
						"fstps (%%eax,%[pos],4)\n\t"

						"movl %[wE],		%%eax\n\t"
						"movl %[wW],		%%ecx\n\t"
						"movl %[wS],		%%edx\n\t"
						"movl %[wN],		%%edi\n\t"

						"flds (%%eax,%[pos],4)\n\t"
						"flds (%%ecx,%[pos],4)\n\t"
						"flds (%%edx,%[pos],4)\n\t"
						"flds (%%edi,%[pos],4)\n\t"

						"movl %[data_div],	%%ecx\n\t"

						"faddp %%st(1)\n\t"
						"faddp %%st(1)\n\t"
						"faddp %%st(1)\n\t"

						"fld1\n\t"
						"fdivp %%st(1)\n\t"
						"fstps (%%ecx,%[pos],4)\n\t"
						"jmp cont\n\t"

					"notnan:\n\t"

						"movl %[dividend],	%%eax\n\t"
						/*"movl %[Cu],		%%ecx\n\t"
						"flds (%%ecx,%[pos],4)\n\t"*/
						"fstps (%%eax,%[pos],4)\n\t"

						"movl %[wE],		%%eax\n\t"
						"movl %[wW],		%%ecx\n\t"
						"movl %[wS],		%%edx\n\t"
						"movl %[wN],		%%edi\n\t"

						"flds (%%eax,%[pos],4)\n\t"
						"flds (%%ecx,%[pos],4)\n\t"
						"flds (%%edx,%[pos],4)\n\t"
						"flds (%%edi,%[pos],4)\n\t"

						"movl %[Du], 		%%eax\n\t"
						"movl %[data_div],	%%ecx\n\t"
						"flds (%%eax,%[pos],4)\n\t"

						"faddp %%st(1)\n\t"
						"faddp %%st(1)\n\t"
						"faddp %%st(1)\n\t"
						"faddp %%st(1)\n\t"

						"fld1\n\t"
						"fdivp %%st(1)\n\t"
						"fstps (%%ecx,%[pos],4)\n\t"

					"cont:\n\t"
					
					"movl %[dividend], %%eax\n\t"
					"movl %[data_div], %%ecx\n\t"
					"movl %[dU], %%edx\n\t"

					"flds (%%eax,%[pos],4)\n\t"
					"flds (%%ecx,%[pos],4)\n\t"
					"flds (%%edx,%[pos],4)\n\t"
					"fld1\n\t"
					"flds %[omega]\n\t"

					"fsubr %%st(0), %%st(1)\n\t"
					"fmulp %%st(0), %%st(3)\n\t"
					"fmulp %%st(1)\n\t"
					"fxch %%st(3)\n\t"
					"faddp %%st(0), %%st(2)\n\t"
					"fmulp %%st(1)\n\t"
					"faddp %%st(1)\n\t"
					
					"fstps (%%edx,%[pos],4)\n\t"
				
				:
				:	[dU]		"m"(dU),
					[data_div]	"m"(data_div_ptr),
					[dividend]	"m"(dividend_ptr),
					[U]		"m"(U),		
					[Cu]		"m"(Cu),
					[Du]		"m"(Du),	
					[wW]		"m"(wW),	
					[wN]		"m"(wN),	
					[wE]		"m"(wE),	
					[wS]		"m"(wS),	
					[iter]		"m"(iter),
					[omega]		"m"(omega),
					[pos]		"r"(pos),	
					[nrows]		"r"(nrows)	
				:"%eax","%ecx","%edx","%edi"
			);


		}
		}

		/*
		-------------------
		--Fill the borders
		-------------------
		*/
		for(j=0;j<ncols;j++)
		{
			pos = j*nrows;

			dU[ pos ] = dU[ pos+1 ];
			dU[ pos+nrows-1 ] = dU[ pos+nrows-2 ];
		}

		for(i=0;i<nrows;i++)
		{
			dU[ i ] = dU[ i+nrows ];
			dU[ i+(ncols-1)*nrows ] = dU[ i+(ncols-2)*nrows ];
		}
	}
	
	free(data_div);
	free(dividend);

	/*Return FPU's initial control word*/
	__asm__ ( 	
			"fldcw %0\n\t"
			:"=m"(fstcw)
		);
#else
/*  #warning "GS_SOR_llin4ASM_2d not compiled!" */
#endif
		
}

/*
//------------------------------------------------------
// TDMA vertical line solving - west most column solver
//------------------------------------------------------
*/
void westColumn4(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j, pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	/*---First row---*/
	i = 0;
	j = 0;
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	/*a = 	0.0f;*/
	b = 	wS->data[pos] + wE->data[pos];
	c = 	-wS->data[pos];
	d =	  wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
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
		b = 	wN->data[pos] + wS->data[pos] + wE->data[pos];
		c = 	-wS->data[pos];
		d =	  wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}
		div = 1.0f / (b-cp[i-1]*a);
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
	b = 	wN->data[pos] + wE->data[pos];
	/*c = 	0.0f;*/
	d =	  wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]);
	
	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}

	dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
	temp1 = dU->data[pos];
	dU->data[pos] = dp[i];

	for(i=nrows-2;i>=0;i--)
	{
		pos = j*nrows+i;
		temp2 = dU->data[pos];
		dU->data[pos] = dp[i] - cp[i]*dU->data[pos+1];

		dU->data[pos+1] = Mparam.omega*dU->data[pos+1] + (1.0f-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	dU->data[pos] = Mparam.omega*dU->data[pos] + (1.0f-Mparam.omega)*temp1;
}

/*
//------------------------------------------------------
// TDMA vertical line solving - middle columns solver
//------------------------------------------------------
*/
void middleColumn4(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j, pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	for(j=1;j<=ncols-2;j++)
	{
		/*---First row---*/
		i = 0;
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;
	
		/*a = 	0.0f;*/
		b = 	wS->data[pos] + wE->data[pos] + wW->data[pos];
		c = 	-wS->data[pos];
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
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
			b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos];
			c = 	-wS->data[pos];
			d =	  wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
				+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
				+ wS->data[pos]*(U->data[spos]-U->data[pos])
				+ wN->data[pos]*(U->data[npos]-U->data[pos]);

			if(!isnan(Cu->data[pos]))
			{
				b += Du->data[pos];
				d += Cu->data[pos];
			}
			div = 1.0f / (b-cp[i-1]*a);
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
		b = 	wN->data[pos] + wE->data[pos] + wW->data[pos];
		/*c = 	0.0f;*/
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}

		dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
		temp1 = dU->data[pos];
		dU->data[pos] = dp[i];
	
		for(i=nrows-2;i>=0;i--)
		{
			pos = j*nrows+i;
			temp2 = dU->data[pos];
			dU->data[pos] = dp[i] - cp[i]*dU->data[pos+1];

			dU->data[pos+1] = Mparam.omega*dU->data[pos+1] + (1.0f-Mparam.omega)*temp1;
			temp1 = temp2;
		}
		dU->data[pos] = Mparam.omega*dU->data[pos] + (1.0f-Mparam.omega)*temp1;
	}
}

/*
//------------------------------------------------------
// TDMA vertical line solving - east most column solver
//------------------------------------------------------
*/
void eastColumn4(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j, pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	/*---First row---*/
	i = 0;
	j = ncols - 1;
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	/*a = 	0.0f;*/
	b = 	wS->data[pos] + wW->data[pos];
	c = 	-wS->data[pos];
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
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
		b = 	wN->data[pos] + wS->data[pos] + wW->data[pos];
		c = 	-wS->data[pos];
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}

		div   =	1.0f / (b-cp[i-1]*a);
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
	b = 	wN->data[pos] + wW->data[pos];
	/*c = 	0.0f;*/
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]);
	
	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}

	dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
	temp1 = dU->data[pos];
	dU->data[pos] = dp[i];

	for(i=nrows-2;i>=0;i--)
	{
		pos = j*nrows+i;
		temp2 = dU->data[pos];
		dU->data[pos] = dp[i] - cp[i]*dU->data[pos+1];

		
		dU->data[pos+1] = Mparam.omega*dU->data[pos+1] + (1.0f-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	dU->data[pos] = Mparam.omega*dU->data[pos] + (1.0f-Mparam.omega)*temp1;

}

/*
//------------------------------------------------------
// TDMA horizontal line solving - top most row solver
//------------------------------------------------------
*/
void northRow4(		struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j, pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;


	/*---First column---*/
	j = 0, i = 0;
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	/*a =	0.0f;*/
	b =	wS->data[pos] + wE->data[pos];
	c =	-wE->data[pos];
	d =	  wE->data[pos]*(U->data[epos]-U->data[pos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
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
		b = 	wS->data[pos] + wE->data[pos] + wW->data[pos];
		c = 	-wE->data[pos];
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}
		div = 1.0f / (b-cp[j-1]*a);
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
	b = 	wS->data[pos] + wW->data[pos];
	/*c = 	0.0f;*/
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}

	dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
	temp1 = dU->data[pos];
	dU->data[pos] = dp[j];

	for(j=ncols-2;j>=0;j--)
	{
		pos = j*nrows+i;
		temp2 = dU->data[pos];
		dU->data[pos] = dp[j] - cp[j]*dU->data[pos+nrows];

		dU->data[pos+nrows] = Mparam.omega*dU->data[pos+nrows] + (1.0f-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	dU->data[pos] = Mparam.omega*dU->data[pos] + (1.0f-Mparam.omega)*temp1;
}

/*
//------------------------------------------------------
// TDMA horizontal line solving - middle rows solver
//------------------------------------------------------
*/
void middleRow4(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j, pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	for(i=1;i<=nrows-2;i++)
	{
		/*---First column---*/
		j = 0;
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		/*a =	0.0f;*/
		b =	wN->data[pos] + wS->data[pos] + wE->data[pos];
		c =	-wE->data[pos];
		d =	  wE->data[pos]*(U->data[epos]-U->data[pos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
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
			b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos];
			c = 	-wE->data[pos];
			d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
				+ wE->data[pos]*(U->data[epos]-U->data[pos])
				+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos])
				+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);
			
			if(!isnan(Cu->data[pos]))
			{
				b += Du->data[pos];
				d += Cu->data[pos];
			}
			div = 1.0f / (b-cp[j-1]*a);
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
		b = 	wN->data[pos] + wS->data[pos] + wW->data[pos];
		/*c = 	0.0f;*/
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);
	
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}

		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
		temp1 = dU->data[pos];
		dU->data[pos] = dp[j];

		for(j=ncols-2;j>=0;j--)
		{
			pos = j*nrows+i;
			temp2 = dU->data[pos];
			dU->data[pos] = dp[j] - cp[j]*dU->data[pos+nrows];

			dU->data[pos+nrows] = Mparam.omega*dU->data[pos+nrows] + (1.0f-Mparam.omega)*temp1;
			temp1 = temp2;
		}
		dU->data[pos] = Mparam.omega*dU->data[pos] + (1.0f-Mparam.omega)*temp1;
	}
}

/*
//------------------------------------------------------
// TDMA horizontal line solving - bottom most row solver
//------------------------------------------------------
*/
void southRow4(		struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j, pos,wpos,epos,npos,spos;
	float a,b,c,d,div,temp1,temp2;

	/*---First column---*/
	j = 0, i=nrows-1;
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	/*a = 	-wW->data[pos];	*/
	b = 	wN->data[pos]  + wE->data[pos];
	c = 	-wE->data[pos];
	d =	  wE->data[pos]*(U->data[epos]-U->data[pos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
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
		b = 	wN->data[pos] + wE->data[pos] + wW->data[pos];
		c = 	-wE->data[pos];
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}

		cp[j] = c / (b-cp[j-1]*a);
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
	}

	/*---Last column---*/
	pos+=nrows;
	wpos+=nrows;
	epos+=nrows;
	npos+=nrows;
	spos+=nrows;

	a = 	-wW->data[pos];	
	b = 	wN->data[pos] +  wW->data[pos];
	/*c = 	-wE->data[pos];*/
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}

	dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
	temp1 = dU->data[pos];
	dU->data[pos] = dp[j];

	for(j=ncols-2;j>=0;j--)
	{
		pos = j*nrows+i;
		temp2 = dU->data[pos];
		dU->data[pos] = dp[j] - cp[j]*dU->data[pos+nrows];

		dU->data[pos+nrows] = Mparam.omega*dU->data[pos+nrows] + (1.0f-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	dU->data[pos] = Mparam.omega*dU->data[pos] + (1.0f-Mparam.omega)*temp1;
}

/*
//------------------------------------------------------
// TDMA vertical line solving - west most column solver
//------------------------------------------------------
*/
void westColumn8(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j,pos,wpos,nwpos,npos,nepos,epos,sepos,spos,swpos;
	float a,b,c,d,div,temp1,temp2;

	/*---First row---*/
	i = 0;
	j = 0;
	pos = j*nrows+i;
	wpos = pos-nrows;
	nwpos = wpos-1;
	swpos = wpos+1;
	epos = pos+nrows;
	nepos = epos-1;
	sepos = epos+1;
	npos = pos-1;
	spos = pos+1;

	/*a = -wN->data[pos];*/
	b = wS->data[pos] + wE->data[pos] + wSE->data[pos];
	c = -wS->data[pos];
	d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
		+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
		+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}
	
	cp[i] = c/b;
	dp[i] = d/b;

	/*---Middle rows---*/
	for(i=1;i<=nrows-2;i++)
	{
		pos++;
		wpos++;
		nwpos++;
		swpos++;
		epos++;
		nepos++;
		sepos++;
		npos++;
		spos++;

		a = -wN->data[pos];
		b = wN->data[pos] + wS->data[pos] + wE->data[pos] + wNE->data[pos] + wSE->data[pos];
		c = -wS->data[pos];
		d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos])
			+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
			+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}
		div = 1 / (b-cp[i-1]*a);
		cp[i] = c*div;
		dp[i] = (d-dp[i-1]*a)*div;
	}

	/*---Last row---*/
	pos++;
	wpos++;
	nwpos++;
	swpos++;
	epos++;
	nepos++;
	sepos++;
	npos++;
	spos++;

	a = -wN->data[pos];
	b = wN->data[pos] + wE->data[pos] + wNE->data[pos];
	/*c = -wS->data[pos];*/
	d =	  wN->data[pos]*(U->data[npos]-U->data[pos])
		+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
		+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos]);
	
	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}

	dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
	temp1 = dU->data[pos];
	dU->data[pos] = dp[i];

	for(i=nrows-2;i>=0;i--)
	{
		pos = j*nrows+i;
		temp2 = dU->data[pos];
		dU->data[pos] = dp[i] - cp[i]*dU->data[pos+1];

		dU->data[pos+1] = Mparam.omega*dU->data[pos+1] + (1-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	dU->data[pos] = Mparam.omega*dU->data[pos] + (1-Mparam.omega)*temp1;
}

/*
//------------------------------------------------------
// TDMA vertical line solving - middle columns solver
//------------------------------------------------------
*/
void middleColumn8(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j,pos,wpos,nwpos,npos,nepos,epos,sepos,spos,swpos;
	float a,b,c,d,div,temp1,temp2;

	for(j=1;j<=ncols-2;j++)
	{
		/*---First row---*/
		i = 0;
		pos = j*nrows+i;
		wpos = pos-nrows;
		nwpos = wpos-1;
		swpos = wpos+1;
		epos = pos+nrows;
		nepos = epos-1;
		sepos = epos+1;
		npos = pos-1;
		spos = pos+1;
	
		/*a = 0.0f;*/
		b = wS->data[pos] + wE->data[pos] + wW->data[pos] + wSE->data[pos] + wSW->data[pos];
		c = -wS->data[pos];
		d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
			+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
			+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos])
			+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}
		
		cp[i] = c/b;
		dp[i] = d/b;
	
		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos++;
			wpos++;
			nwpos++;
			swpos++;
			epos++;
			nepos++;
			sepos++;
			npos++;
			spos++;
	
			a = -wN->data[pos];
			b = wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos] 
			 + wNW->data[pos] + wNE->data[pos] + wSW->data[pos] + wSE->data[pos];
			c = -wS->data[pos];
			d =	  wN->data[pos]*(U->data[npos]-U->data[pos])
				+ wS->data[pos]*(U->data[spos]-U->data[pos])
				+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
				+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
				+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
				+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
				+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos])
				+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos]);

			if(!isnan(Cu->data[pos]))
			{
				b += Du->data[pos];
				d += Cu->data[pos];
			}
			div = 1 / (b-cp[i-1]*a);
			cp[i] = c*div;
			dp[i] = (d-dp[i-1]*a)*div;
		}
	
		/*---Last row---*/
		pos++;
		wpos++;
		nwpos++;
		swpos++;
		epos++;
		nepos++;
		sepos++;
		npos++;
		spos++;
	
		a = -wN->data[pos];
		b = wN->data[pos] + wE->data[pos] + wW->data[pos] + wNW->data[pos] + wNE->data[pos];
		/*c = -wS->data[pos];*/
		d =	  wN->data[pos]*(U->data[npos]-U->data[pos])
			+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
			+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
			+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}

		dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
		temp1 = dU->data[pos];
		dU->data[pos] = dp[i];
	
		for(i=nrows-2;i>=0;i--)
		{
			pos = j*nrows+i;
			temp2 = dU->data[pos];
			dU->data[pos] = dp[i] - cp[i]*dU->data[pos+1];

			dU->data[pos+1] = Mparam.omega*dU->data[pos+1] + (1-Mparam.omega)*temp1;
			temp1 = temp2;
		}
		dU->data[pos] = Mparam.omega*dU->data[pos] + (1-Mparam.omega)*temp1;
	}
}

/*
//------------------------------------------------------
// TDMA vertical line solving - east most column solver
//------------------------------------------------------
*/
void eastColumn8(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j,pos,wpos,nwpos,npos,nepos,epos,sepos,spos,swpos;
	float a,b,c,d,div,temp1,temp2;

	/*---First row---*/
	i = 0;
	j = ncols - 1;
	pos = j*nrows+i;
	wpos = pos-nrows;
	nwpos = wpos-1;
	swpos = wpos+1;
	epos = pos+nrows;
	nepos = epos-1;
	sepos = epos+1;
	npos = pos-1;
	spos = pos+1;

	/*a = 0.0f;*/
	b = wS->data[pos] + wW->data[pos] + wSW->data[pos];
	c = -wS->data[pos];
	d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
  		+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
		+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos]);
		
	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}

	cp[i] = c/b;
	dp[i] = d/b;

	/*---Middle rows---*/
	for(i=1;i<=nrows-2;i++)
	{
		pos++;
		wpos++;
		nwpos++;
		swpos++;
		epos++;
		nepos++;
		sepos++;
		npos++;
		spos++;

		a = -wN->data[pos];
		b = wN->data[pos] + wS->data[pos] + wW->data[pos] + wNW->data[pos] + wSW->data[pos];
		c = -wS->data[pos];
		d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos])
			+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
			+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
			+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}

		div   =	1 / (b-cp[i-1]*a);
		cp[i] = c*div;
		dp[i] = (d-dp[i-1]*a)*div;
	}

	/*---Last row---*/
	pos++;
	wpos++;
	nwpos++;
	swpos++;
	epos++;
	nepos++;
	sepos++;
	npos++;
	spos++;

	a = -wN->data[pos];
	b = wN->data[pos] + wW->data[pos] + wNW->data[pos];
	/*c = 	0.0f;*/
	d =	  wN->data[pos]*(U->data[npos]-U->data[pos])
		+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
		+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos]);
	
	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}

	dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
	temp1 = dU->data[pos];
	dU->data[pos] = dp[i];

	for(i=nrows-2;i>=0;i--)
	{
		pos = j*nrows+i;
		temp2 = dU->data[pos];
		dU->data[pos] = dp[i] - cp[i]*dU->data[pos+1];

		dU->data[pos+1] = Mparam.omega*dU->data[pos+1] + (1-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	dU->data[pos] = Mparam.omega*dU->data[pos] + (1-Mparam.omega)*temp1;

}

/*
//------------------------------------------------------
// TDMA horizontal line solving - top most row solver
//------------------------------------------------------
*/
void northRow8(		struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j,pos,wpos,nwpos,npos,nepos,epos,sepos,spos,swpos;
	float a,b,c,d,div,temp1,temp2;

	/*---First column---*/
	j = 0, i = 0;
	pos = j*nrows+i;
	wpos = pos-nrows;
	nwpos = wpos-1;
	swpos = wpos+1;
	epos = pos+nrows;
	nepos = epos-1;
	sepos = epos+1;
	npos = pos-1;
	spos = pos+1;

	/*a = -wW->data[pos];	*/
	b = wS->data[pos] + wE->data[pos] + wSE->data[pos];
	c = -wE->data[pos];
	d =   wE->data[pos]*(U->data[epos]-U->data[pos])
		+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}
	
	cp[j] = c/b;
	dp[j] = d/b;

	/*---Middle columns---*/
	for(j=1;j<=ncols-2;j++)
	{
		pos+=nrows;
		wpos+=nrows;
		nwpos+=nrows;
		swpos+=nrows;
		epos+=nrows;
		nepos+=nrows;
		sepos+=nrows;
		npos+=nrows;
		spos+=nrows;

		a = -wW->data[pos];	
		b = wS->data[pos] + wE->data[pos] + wW->data[pos] + wSW->data[pos] + wSE->data[pos];
		c = -wE->data[pos];
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos])
			+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos])
			+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}
		div = 1 / (b-cp[j-1]*a);
		cp[j] = c*div;
		dp[j] = (d-dp[j-1]*a)*div;
	}

	/*---Last column---*/
	pos+=nrows;
	wpos+=nrows;
	nwpos+=nrows;
	swpos+=nrows;
	epos+=nrows;
	nepos+=nrows;
	sepos+=nrows;
	npos+=nrows;
	spos+=nrows;

	a = -wW->data[pos];	
	b = wS->data[pos] + wW->data[pos] + wSW->data[pos];
	/*c = 	-wE->data[pos];*/
	d =   wW->data[pos]*(U->data[wpos]-U->data[pos])
		+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}

	dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
	temp1 = dU->data[pos];
	dU->data[pos] = dp[j];

	for(j=ncols-2;j>=0;j--)
	{
		pos = j*nrows+i;
		temp2 = dU->data[pos];
		dU->data[pos] = dp[j] - cp[j]*dU->data[pos+nrows];

		dU->data[pos+nrows] = Mparam.omega*dU->data[pos+nrows] + (1-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	dU->data[pos] = Mparam.omega*dU->data[pos] + (1-Mparam.omega)*temp1;
}

/*
//------------------------------------------------------
// TDMA horizontal line solving - middle rows solver
//------------------------------------------------------
*/
void middleRow8(	struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j,pos,wpos,nwpos,npos,nepos,epos,sepos,spos,swpos;
	float a,b,c,d,div,temp1,temp2;

	for(i=1;i<=nrows-2;i++)
	{
		/*---First column---*/
		j = 0;
		pos = j*nrows+i;
		wpos = pos-nrows;
		nwpos = wpos-1;
		swpos = wpos+1;
		epos = pos+nrows;
		nepos = epos-1;
		sepos = epos+1;
		npos = pos-1;
		spos = pos+1;

		/*a = -wW->data[pos];	*/
		b = wN->data[pos] + wS->data[pos] + wE->data[pos] + wNE->data[pos] + wSE->data[pos];
		c = -wE->data[pos];
		d =	  wE->data[pos]*(U->data[epos]-U->data[pos])
			+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
			+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}
		
		cp[j] = c/b;
		dp[j] = d/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos+=nrows;
			wpos+=nrows;
			nwpos+=nrows;
			swpos+=nrows;
			epos+=nrows;
			nepos+=nrows;
			sepos+=nrows;
			npos+=nrows;
			spos+=nrows;

			a = -wW->data[pos];	
			b = wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos]
				+ wNW->data[pos] + wNE->data[pos] + wSW->data[pos] + wSE->data[pos];
			c = -wE->data[pos];
			d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
				+ wE->data[pos]*(U->data[epos]-U->data[pos])
				+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
				+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
				+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos])
				+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos])
				+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos])
				+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);
			
			if(!isnan(Cu->data[pos]))
			{
				b += Du->data[pos];
				d += Cu->data[pos];
			}
			div = 1 / (b-cp[j-1]*a);
			cp[j] = c*div;
			dp[j] = (d-dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos+=nrows;
		wpos+=nrows;
		nwpos+=nrows;
		swpos+=nrows;
		epos+=nrows;
		nepos+=nrows;
		sepos+=nrows;
		npos+=nrows;
		spos+=nrows;

		a = 	-wW->data[pos];	
		b = 	wN->data[pos] + wS->data[pos] + wW->data[pos] + wNW->data[pos] + wSW->data[pos];
		/*c = 	-wE->data[pos];*/
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
			+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
			+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);
	
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}

		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
		temp1 = dU->data[pos];
		dU->data[pos] = dp[j];

		for(j=ncols-2;j>=0;j--)
		{
			pos = j*nrows+i;
			temp2 = dU->data[pos];
			dU->data[pos] = dp[j] - cp[j]*dU->data[pos+nrows];

			dU->data[pos+nrows] = Mparam.omega*dU->data[pos+nrows] + (1-Mparam.omega)*temp1;
			temp1 = temp2;
		}
		dU->data[pos] = Mparam.omega*dU->data[pos] + (1-Mparam.omega)*temp1;
	}
}

/*
//------------------------------------------------------
// TDMA horizontal line solving - bottom most row solver
//------------------------------------------------------
*/
void southRow8(		struct matrixM *U,
			struct matrixM *dU,
			struct matrixM *Cu,
			struct matrixM *Du,
			struct matrixM *wW,
			struct matrixM *wNW,
			struct matrixM *wN,
			struct matrixM *wNE,
			struct matrixM *wE,
			struct matrixM *wSE,
			struct matrixM *wS,
			struct matrixM *wSW,
			struct mparam Mparam,
			float *cp,
			float *dp,
			int ncols,
			int nrows)
{
	int i,j,pos,wpos,nwpos,npos,nepos,epos,sepos,spos,swpos;
	float a,b,c,d,div,temp1,temp2;

	/*---First column---*/
	j = 0, i=nrows-1;
	pos = j*nrows+i;
	wpos = pos-nrows;
	nwpos = wpos-1;
	swpos = wpos+1;
	epos = pos+nrows;
	nepos = epos-1;
	sepos = epos+1;
	npos = pos-1;
	spos = pos+1;

	/*a = -wW->data[pos];	*/
	b = wN->data[pos] + wE->data[pos] + wNE->data[pos];
	c = -wE->data[pos];
	d =	  wE->data[pos]*(U->data[epos]-U->data[pos])
		+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}
	
	cp[j] = c/b;
	dp[j] = d/b;

	/*---Middle columns---*/
	for(j=1;j<=ncols-2;j++)
	{
		pos+=nrows;
		wpos+=nrows;
		nwpos+=nrows;
		swpos+=nrows;
		epos+=nrows;
		nepos+=nrows;
		sepos+=nrows;
		npos+=nrows;
		spos+=nrows;

		a = -wW->data[pos];	
		b = wN->data[pos] + wE->data[pos] + wW->data[pos] + wNW->data[pos] + wNE->data[pos];
		c = -wE->data[pos];
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos])
			+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
			+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
		}

		cp[j] = c / (b-cp[j-1]*a);
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
	}

	/*---Last column---*/
	pos+=nrows;
	wpos+=nrows;
	nwpos+=nrows;
	swpos+=nrows;
	epos+=nrows;
	nepos+=nrows;
	sepos+=nrows;
	npos+=nrows;
	spos+=nrows;

	a = -wW->data[pos];	
	b = wN->data[pos] + wW->data[pos] + wNW->data[pos];
	/*c = -wE->data[pos];*/
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
		+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
	}

	dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
	temp1 = dU->data[pos];
	dU->data[pos] = dp[j];

	for(j=ncols-2;j>=0;j--)
	{
		pos = j*nrows+i;
		temp2 = dU->data[pos];
		dU->data[pos] = dp[j] - cp[j]*dU->data[pos+nrows];

		dU->data[pos+nrows] = Mparam.omega*dU->data[pos+nrows] + (1-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	dU->data[pos] = Mparam.omega*dU->data[pos] + (1-Mparam.omega)*temp1;
}
