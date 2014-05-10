/*
Functions related to variational optical-flow calculation.

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

#if !defined( isnan )
  #define isnan(x)((x)!=(x))
#endif

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel succesive over relaxation with 2D (spatial) regularization.
//---Early lizearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_elin4_2d(	struct matrixM *U, 
			struct matrixM *V, 
			struct matrixM *M,
			struct matrixM *Cu,
			struct matrixM *Cv,
			struct matrixM *Du,
			struct matrixM *Dv,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam)
{
	float wNeighU, wNeighV, float_temp1, float_temp2, float_temp3;
	float U_newapprox, V_newapprox;
	float *u_data_div, *v_data_div;
	unsigned int i,j,k,pos,iter;
	unsigned int nrows = U->dimElems[0];
	unsigned int ncols = U->dimElems[1];
	
	/*Reserve space for the divisor terms*/
	if( (u_data_div = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_elin4_2d: error reserving space for 'u_data_div'\n");
		return;
	}
	if( (v_data_div = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_elin4_2d: error reserving space for 'v_data_div'\n");
		free(u_data_div);
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

			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			/*Weights for U*/
			wNeighU 	= U->data[ pos - nrows ]*wW->data[ pos ];
			float_temp1 	= U->data[ pos + nrows ]*wE->data[ pos ];
			wNeighU 	+= float_temp1;
			
			float_temp2 	= U->data[ pos - 1 ]*wN->data[ pos ];
			float_temp3 	= U->data[ pos + 1 ]*wS->data[ pos ];
			float_temp2 	+= float_temp3;
			
			wNeighU 	+= float_temp2;			

			/*Weights for V*/
			wNeighV 	= V->data[ pos - nrows ]*wW->data[ pos ];
			float_temp1 	= V->data[ pos + nrows ]*wE->data[ pos ];
			wNeighV 	+= float_temp1;

			float_temp2	= V->data[ pos - 1 ]*wN->data[ pos ];
			float_temp3	= V->data[ pos + 1 ]*wS->data[ pos ];
			float_temp2 	+= float_temp3;
			
			wNeighV 	+= float_temp2;


			if(iter==0)
			{
				/*Sum of weights*/
				float_temp1 = wW->data[ pos ] + wE->data[ pos ];
				float_temp2 = wN->data[ pos ] + wS->data[ pos ];
				float_temp1 += float_temp2;

				if(isnan(Du->data[pos]))
					u_data_div[pos] = 1.0f/( float_temp1 );
				else
					u_data_div[pos] = 1.0f/( float_temp1 + Du->data[pos] );

				if(isnan(Dv->data[pos]))
					v_data_div[pos] = 1.0f/( float_temp1 );
				else
					v_data_div[pos] = 1.0f/( float_temp1 + Dv->data[pos] );
			}

			if(isnan(Cu->data[pos]))
				U_newapprox = (wNeighU)*u_data_div[pos];
			else{
				float_temp1 = wNeighU + Cu->data[pos];
				float_temp2 = M->data[pos]*V->data[pos];
				float_temp1 = float_temp1 - float_temp2;
				U_newapprox = float_temp1 * u_data_div[pos];

				/*U_newapprox = (wNeighU - M->data[pos]*V->data[pos] + Cu->data[pos])*u_data_div[pos];*/
			}

			if(isnan(Cv->data[pos]))
				V_newapprox = (wNeighV)*v_data_div[pos];
			else{
				float_temp1 = wNeighV + Cv->data[pos];
				float_temp2 = M->data[pos]*U->data[pos];
				float_temp1 = float_temp1 - float_temp2;
				V_newapprox = float_temp1 * v_data_div[pos];

				/*V_newapprox = (wNeighV - M->data[pos]*U->data[pos] + Cv->data[pos])*v_data_div[pos];*/
			}

			U->data[ pos ] = (1.0f-Mparam.omega)*U->data[ pos ] + Mparam.omega*U_newapprox;
			V->data[ pos ] = (1.0f-Mparam.omega)*V->data[ pos ] + Mparam.omega*V_newapprox;
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
	
			U->data[ pos ] = U->data[ pos+1 ];
			U->data[ pos+nrows-1 ] = U->data[ pos+nrows-2 ];
	
			V->data[ pos ] = V->data[ pos+1 ];
			V->data[ pos+nrows-1 ] = V->data[ pos+nrows-2 ];
		}
	
		for(i=0;i<nrows;i++)
		{
			U->data[ i ] = U->data[ i+nrows ];
			U->data[ i+(ncols-1)*nrows ] = U->data[ i+(ncols-2)*nrows ];
	
			V->data[ i ] = V->data[ i+nrows ];
			V->data[ i+(ncols-1)*nrows ] = V->data[ i+(ncols-2)*nrows ];
		}
	
	}
	

	free(u_data_div);
	free(v_data_div);
}

/*
//---------------------------------------------------------------------------------------
//---Alternating line relaxation Gauss-Seidel (Gauss-Seidel block variant)
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization.
//---Early linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void	GS_ALR_SOR_elin4_2d(	struct matrixM *U, 
				struct matrixM *V, 
				struct matrixM *M,
				struct matrixM *Cu,
				struct matrixM *Cv,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam)
{
	float *Cp, *Dp;
	int iter;
	int nrows = Cu->dimElems[0];
	int ncols = Cu->dimElems[1];
	int maxLength;
	int iterations = (int)Mparam.iter;

	maxLength = (nrows>ncols)?nrows:ncols;

	/*Reserve space for Cprime and Dprime temp variables used in TDMA*/
	if( (Cp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("GS_ALR_SOR_elin4_2d: error reserving space for 'Cp'\n");
		return;
	}
	if( (Dp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("GS_ALR_SOR_elin4_2d: error reserving space for 'Dp'\n");
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
		westColumn_elin4(	U, V, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_elin4(	U, V, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_elin4(	U, V, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);

		westColumn_elin4(	V, U, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_elin4(	V, U, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_elin4(	V, U, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		northRow_elin4(	V, U, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleRow_elin4(V, U, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		southRow_elin4(	V, U, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		
		northRow_elin4(	U, V, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleRow_elin4(U, V, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		southRow_elin4(	U, V, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
	}

	free( Cp );
	free( Dp );
}

/*
//---------------------------------------------------------------------
// Calculate residuals r (Ax=b => r=b-Ax) for early linearization case
//---------------------------------------------------------------------
*/
void Residuals_elin4_2d(	struct matrixM *RU, 
				struct matrixM *RV, 
				struct matrixM *U,
				struct matrixM *V,
				struct matrixM *M,
				struct matrixM *Cu,
				struct matrixM *Cv,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam)
{

	float wNeighU, wNeighV, float_temp1, float_temp2, float_temp3;
	unsigned int i,j,k,pos,posOS,wpos,epos,npos,spos;
	unsigned int nrows = M->dimElems[0];
	unsigned int ncols = M->dimElems[1];
	unsigned int nframes = 1, frameoffset;

	if(M->ndims>2)
		nframes = M->dimElems[2];
	
	/*
	--------------------------------------
	--Calculate residuals Ax=b => r=b-Ax
	--------------------------------------
	*/
	for(k=0;k<nframes;k++)
	{
		frameoffset = k*nrows*ncols;
		for(j=1;j<ncols-1;j++)
		{
			for(i=1;i<nrows-1;i++)		
			{
				/*Spatial position indices*/
				pos = j*nrows+i;
				posOS = pos + frameoffset;	/*Offsetted position*/
				wpos = pos - nrows;
				epos = pos + nrows;
				npos = pos - 1;
				spos = pos + 1;

				/*
				---------------------
				--Spatial diffusion--
				---------------------	
				*/
		
				wNeighU =	U->data[ wpos ]*wW->data[ pos ] + 
						U->data[ epos ]*wE->data[ pos ] +
						U->data[ npos ]*wN->data[ pos ] + 
						U->data[ spos ]*wS->data[ pos ];
		
				wNeighV =	V->data[ wpos ]*wW->data[ pos ] + 
						V->data[ epos ]*wE->data[ pos ] +
						V->data[ npos ]*wN->data[ pos ] + 
						V->data[ spos ]*wS->data[ pos ];
		
				/*Sum of weights*/
				float_temp1 = wW->data[ pos ] + wE->data[ pos ];
				float_temp2 = wN->data[ pos ] + wS->data[ pos ];
				float_temp1 += float_temp2;

				/*We don't like NaN:s, do we*/
				if( !isnan(Cu->data[posOS]) ){
					RU->data[posOS] = Cu->data[posOS] - M->data[posOS]*V->data[pos] + wNeighU -(Du->data[posOS]+float_temp1)*U->data[pos];
				}else{
					RU->data[posOS] = wNeighU -(float_temp1)*U->data[pos];
				}
				if( !isnan(Cv->data[posOS]) ){
					RV->data[posOS] = Cv->data[posOS] - M->data[posOS]*U->data[pos] + wNeighV -(Dv->data[posOS]+float_temp1)*V->data[pos];
				}else{
					RV->data[posOS] = wNeighV -(float_temp1)*V->data[pos];
				}
			}
		}
	}

	/*
	-------------------
	--Fill the borders
	-------------------
	*/
	for(k=0;k<nframes;k++)
	{
		frameoffset = k*nrows*ncols;
		for(j=0;j<ncols;j++)
		{
			posOS = j*nrows + frameoffset;
	
			RU->data[ posOS ] = RU->data[ posOS+1 ];
			RU->data[ posOS+nrows-1 ] = RU->data[ posOS+nrows-2 ];
	
			RV->data[ posOS ] = RV->data[ posOS+1 ];
			RV->data[ posOS+nrows-1 ] = RV->data[ posOS+nrows-2 ];
		}
	
		for(i=0;i<nrows;i++)
		{
			posOS = i +  frameoffset;
	
			RU->data[ posOS ] = RU->data[ posOS+nrows ];
			RU->data[ posOS+(ncols-1)*nrows ] = RU->data[ posOS+(ncols-2)*nrows ];
	
			RV->data[ posOS ] = RV->data[ posOS+nrows ];
			RV->data[ posOS+(ncols-1)*nrows ] = RV->data[ posOS+(ncols-2)*nrows ];
		}
	}
}

/*
//----------------------------------------------------------------------
// Calculates LHS (Left hand size of Ax=b) for early linearization case
//----------------------------------------------------------------------
*/
void LHS_elin4_2d(	struct matrixM *AU, 
			struct matrixM *AV, 
			struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
			struct matrixM *Du,
			struct matrixM *Dv,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam)
{

	float wNeighU, wNeighV, float_temp1, float_temp2, float_temp3;
	unsigned int i,j,k,pos,posOS,wpos,epos,npos,spos;
	unsigned int nrows = M->dimElems[0];
	unsigned int ncols = M->dimElems[1];
	unsigned int nframes = 1, frameoffset;

	if(M->ndims>2)
		nframes = M->dimElems[2];
	
	/*
	----------------
	--Calculates Ax
	----------------
	*/
	for(k=0;k<nframes;k++)
	{
		frameoffset = k*nrows*ncols;
		for(j=1;j<ncols-1;j++)
		{
			for(i=1;i<nrows-1;i++)		
			{
				/*Spatial position indices*/
				pos = j*nrows+i;
				posOS = pos + frameoffset;	/*Offsetted position*/
				wpos = pos - nrows;
				epos = pos + nrows;
				npos = pos - 1;
				spos = pos + 1;

				/*
				---------------------
				--Spatial diffusion--
				---------------------	
				*/
		
				wNeighU =	U->data[ wpos ]*wW->data[ pos ] + 
						U->data[ epos ]*wE->data[ pos ] +
						U->data[ npos ]*wN->data[ pos ] + 
						U->data[ spos ]*wS->data[ pos ];
		
				wNeighV =	V->data[ wpos ]*wW->data[ pos ] + 
						V->data[ epos ]*wE->data[ pos ] +
						V->data[ npos ]*wN->data[ pos ] + 
						V->data[ spos ]*wS->data[ pos ];
		
				/*Sum of weights*/
				float_temp1 = wW->data[ pos ] + wE->data[ pos ];
				float_temp2 = wN->data[ pos ] + wS->data[ pos ];
				float_temp1 += float_temp2;

				/*We don't like NaN:s, do we*/
				if( !isnan(Du->data[posOS]) ){
					AU->data[posOS] = M->data[posOS]*V->data[pos] - wNeighU + (Du->data[posOS]+float_temp1)*U->data[pos];
				}else{
					AU->data[posOS] = - wNeighU + float_temp1*U->data[pos];
				}
				if( !isnan(Dv->data[posOS]) ){
					AV->data[posOS] = M->data[posOS]*U->data[pos] - wNeighV + (Dv->data[posOS]+float_temp1)*V->data[pos];
				}else{
					AV->data[posOS] = - wNeighV + float_temp1*V->data[pos];
				}
			}
		}
	}

	/*
	-------------------
	--Fill the borders
	-------------------
	*/
	for(k=0;k<nframes;k++)
	{
		frameoffset = k*nrows*ncols;
		for(j=0;j<ncols;j++)
		{
			posOS = j*nrows + frameoffset;
	
			AU->data[ posOS ] = AU->data[ posOS+1 ];
			AU->data[ posOS+nrows-1 ] = AU->data[ posOS+nrows-2 ];
	
			AV->data[ posOS ] = AV->data[ posOS+1 ];
			AV->data[ posOS+nrows-1 ] = AV->data[ posOS+nrows-2 ];
		}
		for(i=0;i<nrows;i++)
		{
			posOS = i + frameoffset;
	
			AU->data[ posOS ] = AU->data[ posOS+nrows ];
			AU->data[ posOS+(ncols-1)*nrows ] = AU->data[ posOS+(ncols-2)*nrows ];
	
			AV->data[ posOS ] = AV->data[ posOS+nrows ];
			AV->data[ posOS+(ncols-1)*nrows ] = AV->data[ posOS+(ncols-2)*nrows ];
		}
	}

}

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel successive over relaxation with 2D (spatial) regularization
//---Late linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin4_2d(	struct matrixM *U, 
			struct matrixM *V, 
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
			struct matrixM *Cu,
			struct matrixM *Cv,
			struct matrixM *Du,
			struct matrixM *Dv,
			struct matrixM *wW,
			struct matrixM *wN,
			struct matrixM *wE,
			struct matrixM *wS,
			struct mparam Mparam)
{
	float wNeighU, wNeighV, float_temp1, float_temp2, float_temp3;
	float dU_newapprox, dV_newapprox;
	float *u_data_div, *v_data_div;
	unsigned int i,j,k,pos,iter;
	unsigned int nrows = M->dimElems[0];
	unsigned int ncols = M->dimElems[1];
		
	/*Reserve space for the divisor terms*/
	if( (u_data_div = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin4_2d: error reserving space for 'u_data_div'\n");
		return;
	}
	if( (v_data_div = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin4_2d: error reserving space for 'v_data_div'\n");
		free(u_data_div);
		return;
	}

	for(iter=0;iter<(unsigned int)Mparam.iter;iter++)
	{
		for(j=1;j<ncols-1;j++)
		{
			for(i=1;i<nrows-1;i++)		
			{
				/*Spatial position indices*/
				pos = j*nrows+i;
				/*wpos = pos - nrows;
				epos = pos + nrows;
				npos = pos - 1;
				spos = pos + 1;*/

				/*
				---------------------
				--Spatial diffusion--
				---------------------	
				*/
				/*
				wNeighU =	(dU->data[ wpos ] + U->data[ wpos ] - U->data[ pos ])*wW->data[ pos ] + 
						(dU->data[ epos ] + U->data[ epos ] - U->data[ pos ])*wE->data[ pos ] +
						(dU->data[ npos ] + U->data[ npos ] - U->data[ pos ])*wN->data[ pos ] + 
						(dU->data[ spos ] + U->data[ spos ] - U->data[ pos ])*wS->data[ pos ];
				*/
				wNeighU		= dU->data[ pos - nrows ] + U->data[ pos - nrows ];
				float_temp1	= dU->data[ pos + nrows ] + U->data[ pos + nrows ];
				float_temp2	= dU->data[ pos - 1 ] + U->data[ pos - 1 ];
				float_temp3	= dU->data[ pos + 1 ] + U->data[ pos + 1 ];
				
				wNeighU		-= U->data[ pos ];
				float_temp1	-= U->data[ pos ];
				float_temp2	-= U->data[ pos ];
				float_temp3	-= U->data[ pos ];

				wNeighU		*= wW->data[ pos ];
				float_temp1	*= wE->data[ pos ];
				float_temp2	*= wN->data[ pos ];
				float_temp3	*= wS->data[ pos ];

				wNeighU		+= float_temp1;
				float_temp2	+= float_temp3;
				wNeighU		+= float_temp2;

				/*wNeighV =	(dV->data[ wpos ] + V->data[ wpos ] - V->data[ pos ])*wW->data[ pos ] + 
						(dV->data[ epos ] + V->data[ epos ] - V->data[ pos ])*wE->data[ pos ] +
						(dV->data[ npos ] + V->data[ npos ] - V->data[ pos ])*wN->data[ pos ] + 
						(dV->data[ spos ] + V->data[ spos ] - V->data[ pos ])*wS->data[ pos ];
				*/
				wNeighV		= dV->data[ pos - nrows ] + V->data[ pos - nrows ];
				float_temp1	= dV->data[ pos + nrows ] + V->data[ pos + nrows ];
				float_temp2	= dV->data[ pos - 1 ] + V->data[ pos - 1 ];
				float_temp3	= dV->data[ pos + 1 ] + V->data[ pos + 1 ];
				
				wNeighV		-= V->data[ pos ];
				float_temp1	-= V->data[ pos ];
				float_temp2	-= V->data[ pos ];
				float_temp3	-= V->data[ pos ];

				wNeighV		*= wW->data[ pos ];
				float_temp1	*= wE->data[ pos ];
				float_temp2	*= wN->data[ pos ];
				float_temp3	*= wS->data[ pos ];

				wNeighV		+= float_temp1;
				float_temp2	+= float_temp3;
				wNeighV		+= float_temp2;

				if(iter==0)
				{
					/*Sum of weights*/
					float_temp1 = wW->data[ pos ] + wE->data[ pos ];
					float_temp2 = wN->data[ pos ] + wS->data[ pos ];
					float_temp1 += float_temp2;

					if(isnan(Du->data[pos]))
						u_data_div[pos] = 1.0f/( float_temp1 );
					else
						u_data_div[pos] = 1.0f/( float_temp1 + Du->data[pos] );

					if(isnan(Dv->data[pos]))
						v_data_div[pos] = 1.0f/( float_temp1 );
					else
						v_data_div[pos] = 1.0f/( float_temp1 + Dv->data[pos] );
				}

				if(isnan(Cu->data[pos]))
					dU_newapprox = wNeighU*u_data_div[pos];
				else{
					float_temp1 = wNeighU + Cu->data[pos];
					float_temp2 = M->data[pos]*dV->data[pos];
					float_temp1 -= float_temp2;
					dU_newapprox = float_temp1 * u_data_div[pos];
								
					/*dU_newapprox = ( wNeighU - M->data[pos]*dV->data[pos] + Cu->data[pos] )*u_data_div[pos];*/
				}
				
				if(isnan(Cv->data[pos]))
					dV_newapprox = wNeighV*v_data_div[pos];
				else{
					float_temp1 = wNeighV + Cv->data[pos];
					float_temp2 = M->data[pos]*dU->data[pos];
					float_temp1 -= float_temp2;
					dV_newapprox = float_temp1 * v_data_div[pos];
				
					/*dV_newapprox = ( wNeighV - M->data[pos]*dU->data[pos] + Cv->data[pos] )*v_data_div[pos];*/
				}

				dU->data[ pos ] = (1.0f-Mparam.omega)*dU->data[ pos ] + Mparam.omega*dU_newapprox;
				dV->data[ pos ] = (1.0f-Mparam.omega)*dV->data[ pos ] + Mparam.omega*dV_newapprox;
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
	
			dV->data[ pos ] = dV->data[ pos+1 ];
			dV->data[ pos+nrows-1 ] = dV->data[ pos+nrows-2 ];
		}
	
		for(i=0;i<nrows;i++)
		{
			dU->data[ i ] = dU->data[ i+nrows ];
			dU->data[ i+(ncols-1)*nrows ] = dU->data[ i+(ncols-2)*nrows ];
	
			dV->data[ i ] = dV->data[ i+nrows ];
			dV->data[ i+(ncols-1)*nrows ] = dV->data[ i+(ncols-2)*nrows ];
		}
	}
	
	
	free(u_data_div);
	free(v_data_div);
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
				struct matrixM *V, 
				struct matrixM *dU,
				struct matrixM *dV,
				struct matrixM *M,
				struct matrixM *Cu,
				struct matrixM *Cv,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam)
{

	float *Cp, *Dp;
	int iter;
	int nrows = Cu->dimElems[0];
	int ncols = Cu->dimElems[1];
	int maxLength;
	int iterations = (int)Mparam.iter;

	maxLength = (nrows>ncols)?nrows:ncols;

	/*Reserve space for Cprime and Dprime temp variables used in TDMA*/
	if( (Cp = (float*)malloc( maxLength*sizeof(float) ))==NULL )
	{
		printf("GS_ALR_SOR_llin4_2d: Error reserving space for 'Cp'\n");
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
		westColumn_llin4(	U, V, dU, dV, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_llin4(	U, V, dU, dV, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_llin4(	U, V, dU, dV, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);

		westColumn_llin4(	V, U, dV, dU, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_llin4(	V, U, dV, dU, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_llin4(	V, U, dV, dU, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		northRow_llin4(	V, U, dV, dU, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleRow_llin4(V, U, dV, dU, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		southRow_llin4(	V, U, dV, dU, M, Cv, Dv, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		
		northRow_llin4(	U, V, dU, dV, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		middleRow_llin4(U, V, dU, dV, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
		southRow_llin4(	U, V, dU, dV, M, Cu, Du, wW, wN, wE, wS, Mparam, Cp, Dp, ncols, nrows);
	}

	free( Cp );
	free( Dp );
}

/*
//---------------------------------------------------------------------------------------------------
// Calculate residuals r (Ax=b => r=b-Ax) for late linearization case with 4-neighbours
//---------------------------------------------------------------------------------------------------
*/
void Residuals_llin4_2d(	struct matrixM *RU,
				struct matrixM *RV,
				struct matrixM *U, 
				struct matrixM *V, 
				struct matrixM *dU,
				struct matrixM *dV,
				struct matrixM *M,
				struct matrixM *Cu,
				struct matrixM *Cv,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam)
{

	float wNeighU, wNeighV, float_temp1, float_temp2, float_temp3;
	unsigned int i,j,k,pos,posOS,iter;
	unsigned int nrows = M->dimElems[0];
	unsigned int ncols = M->dimElems[1];
	unsigned int nframes = 1, frameoffset;

	if(M->ndims>2)
		nframes = M->dimElems[2];

	/*
	--------------------------------------
	--Calculate residuals Ax=b => r=b-Ax
	--------------------------------------
	*/
	for(k=0;k<nframes;k++)
	{
		frameoffset = k*nrows*ncols;
		for(j=1;j<ncols-1;j++)
		{
			for(i=1;i<nrows-1;i++)		
			{
				/*Spatial position indices*/
				pos = j*nrows+i;
				posOS = pos + frameoffset;
				/*wpos = pos - nrows;
				epos = pos + nrows;
				npos = pos - 1;
				spos = pos + 1;*/

				/*
				---------------------
				--Spatial diffusion--
				---------------------	
				*/
				/*
				wNeighU =	(dU->data[ wpos ] + U->data[ wpos ] - U->data[ pos ])*wW->data[ pos ] + 
						(dU->data[ epos ] + U->data[ epos ] - U->data[ pos ])*wE->data[ pos ] +
						(dU->data[ npos ] + U->data[ npos ] - U->data[ pos ])*wN->data[ pos ] + 
						(dU->data[ spos ] + U->data[ spos ] - U->data[ pos ])*wS->data[ pos ];
				*/
				wNeighU		= dU->data[ pos - nrows ] + U->data[ pos - nrows ];
				float_temp1	= dU->data[ pos + nrows ] + U->data[ pos + nrows ];
				float_temp2	= dU->data[ pos - 1 ] + U->data[ pos - 1 ];
				float_temp3	= dU->data[ pos + 1 ] + U->data[ pos + 1 ];
			
				wNeighU		-= U->data[ pos ];
				float_temp1	-= U->data[ pos ];
				float_temp2	-= U->data[ pos ];
				float_temp3	-= U->data[ pos ];

				wNeighU		*= wW->data[ pos ];
				float_temp1	*= wE->data[ pos ];
				float_temp2	*= wN->data[ pos ];
				float_temp3	*= wS->data[ pos ];

				wNeighU		+= float_temp1;
				float_temp2	+= float_temp3;
				wNeighU		+= float_temp2;

				/*wNeighV =	(dV->data[ wpos ] + V->data[ wpos ] - V->data[ pos ])*wW->data[ pos ] + 
						(dV->data[ epos ] + V->data[ epos ] - V->data[ pos ])*wE->data[ pos ] +
						(dV->data[ npos ] + V->data[ npos ] - V->data[ pos ])*wN->data[ pos ] + 
						(dV->data[ spos ] + V->data[ spos ] - V->data[ pos ])*wS->data[ pos ];
				*/
				wNeighV		= dV->data[ pos - nrows ] + V->data[ pos - nrows ];
				float_temp1	= dV->data[ pos + nrows ] + V->data[ pos + nrows ];
				float_temp2	= dV->data[ pos - 1 ] + V->data[ pos - 1 ];
				float_temp3	= dV->data[ pos + 1 ] + V->data[ pos + 1 ];
			
				wNeighV		-= V->data[ pos ];
				float_temp1	-= V->data[ pos ];
				float_temp2	-= V->data[ pos ];
				float_temp3	-= V->data[ pos ];

				wNeighV		*= wW->data[ pos ];
				float_temp1	*= wE->data[ pos ];
				float_temp2	*= wN->data[ pos ];
				float_temp3	*= wS->data[ pos ];

				wNeighV		+= float_temp1;
				float_temp2	+= float_temp3;
				wNeighV		+= float_temp2;

				/*Sum of weights*/
				float_temp1 = wW->data[ pos ] + wE->data[ pos ];
				float_temp2 = wN->data[ pos ] + wS->data[ pos ];
				float_temp1 += float_temp2;

				/*We don't like NaN:s, do we*/
				if( !isnan(Cu->data[posOS]) ){
					RU->data[posOS] = Cu->data[posOS] - M->data[posOS]*dV->data[pos] + wNeighU -(Du->data[posOS]+float_temp1)*dU->data[pos];
				}else{
					RU->data[posOS] = wNeighU -(float_temp1)*dU->data[pos];
				}
				if( !isnan(Cv->data[posOS]) ){
					RV->data[posOS] = Cv->data[posOS] - M->data[posOS]*dU->data[pos] + wNeighV -(Dv->data[posOS]+float_temp1)*dV->data[pos];
				}else{
					RV->data[posOS] = wNeighV -(float_temp1)*dV->data[pos];
				}
			}
		}
	}
	/*
	-------------------
	--Fill the borders
	-------------------
	*/
	for(k=0;k<nframes;k++)
	{
		frameoffset = k*nrows*ncols;
		for(j=0;j<ncols;j++)
		{
			posOS = j*nrows + frameoffset;
	
			RU->data[ posOS ] = RU->data[ posOS+1 ];
			RU->data[ posOS+nrows-1 ] = RU->data[ posOS+nrows-2 ];
	
			RV->data[ posOS ] = RV->data[ posOS+1 ];
			RV->data[ posOS+nrows-1 ] = RV->data[ posOS+nrows-2 ];
		}
	
		for(i=0;i<nrows;i++)
		{
			posOS = i + frameoffset;
	
			RU->data[ posOS ] = RU->data[ posOS+nrows ];
			RU->data[ posOS+(ncols-1)*nrows ] = RU->data[ posOS+(ncols-2)*nrows ];
	
			RV->data[ posOS ] = RV->data[ i+nrows ];
			RV->data[ posOS+(ncols-1)*nrows ] = RV->data[ posOS+(ncols-2)*nrows ];
		}
	}
}

/*
//---------------------------------------------------------------------
// Calculates LHS (Left hand size of Ax=b) for late linearization case
//---------------------------------------------------------------------
*/
void LHS_llin4_2d(		struct matrixM *AU,
				struct matrixM *AV,
				struct matrixM *U, 
				struct matrixM *V, 
				struct matrixM *dU,
				struct matrixM *dV,
				struct matrixM *M,
				struct matrixM *Du,
				struct matrixM *Dv,
				struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct mparam Mparam)
{

	float wNeighU, wNeighV, float_temp1, float_temp2, float_temp3;
	unsigned int i,j,k,pos,posOS,iter;
	unsigned int nrows = M->dimElems[0];
	unsigned int ncols = M->dimElems[1];
	unsigned int nframes = 1, frameoffset;

	if(M->ndims>2)
		nframes = M->dimElems[2];

	/*
	----------------
	--Calculates Ax
	----------------
	*/
	for(k=0;k<nframes;k++)
	{
		frameoffset = k*nrows*ncols;
		for(j=1;j<ncols-1;j++)
		{
			for(i=1;i<nrows-1;i++)		
			{
				/*Spatial position indices*/
				pos = j*nrows+i;
				posOS = pos + frameoffset;
				/*wpos = pos - nrows;
				epos = pos + nrows;
				npos = pos - 1;
				spos = pos + 1;*/

				/*
				---------------------
				--Spatial diffusion--
				---------------------	
				*/
				/*
				wNeighU =	(dU->data[ wpos ] + U->data[ wpos ] - U->data[ pos ])*wW->data[ pos ] + 
						(dU->data[ epos ] + U->data[ epos ] - U->data[ pos ])*wE->data[ pos ] +
						(dU->data[ npos ] + U->data[ npos ] - U->data[ pos ])*wN->data[ pos ] + 
						(dU->data[ spos ] + U->data[ spos ] - U->data[ pos ])*wS->data[ pos ];
				*/
				wNeighU		= dU->data[ pos - nrows ] + U->data[ pos - nrows ];
				float_temp1	= dU->data[ pos + nrows ] + U->data[ pos + nrows ];
				float_temp2	= dU->data[ pos - 1 ] + U->data[ pos - 1 ];
				float_temp3	= dU->data[ pos + 1 ] + U->data[ pos + 1 ];
			
				wNeighU		-= U->data[ pos ];
				float_temp1	-= U->data[ pos ];
				float_temp2	-= U->data[ pos ];
				float_temp3	-= U->data[ pos ];

				wNeighU		*= wW->data[ pos ];
				float_temp1	*= wE->data[ pos ];
				float_temp2	*= wN->data[ pos ];
				float_temp3	*= wS->data[ pos ];

				wNeighU		+= float_temp1;
				float_temp2	+= float_temp3;
				wNeighU		+= float_temp2;

				/*wNeighV =	(dV->data[ wpos ] + V->data[ wpos ] - V->data[ pos ])*wW->data[ pos ] + 
						(dV->data[ epos ] + V->data[ epos ] - V->data[ pos ])*wE->data[ pos ] +
						(dV->data[ npos ] + V->data[ npos ] - V->data[ pos ])*wN->data[ pos ] + 
						(dV->data[ spos ] + V->data[ spos ] - V->data[ pos ])*wS->data[ pos ];
				*/
				wNeighV		= dV->data[ pos - nrows ] + V->data[ pos - nrows ];
				float_temp1	= dV->data[ pos + nrows ] + V->data[ pos + nrows ];
				float_temp2	= dV->data[ pos - 1 ] + V->data[ pos - 1 ];
				float_temp3	= dV->data[ pos + 1 ] + V->data[ pos + 1 ];
			
				wNeighV		-= V->data[ pos ];
				float_temp1	-= V->data[ pos ];
				float_temp2	-= V->data[ pos ];
				float_temp3	-= V->data[ pos ];

				wNeighV		*= wW->data[ pos ];
				float_temp1	*= wE->data[ pos ];
				float_temp2	*= wN->data[ pos ];
				float_temp3	*= wS->data[ pos ];

				wNeighV		+= float_temp1;
				float_temp2	+= float_temp3;
				wNeighV		+= float_temp2;

				/*Sum of weights*/
				float_temp1 = wW->data[ pos ] + wE->data[ pos ];
				float_temp2 = wN->data[ pos ] + wS->data[ pos ];
				float_temp1 += float_temp2;

				/*We don't like NaN:s, do we*/
				if( !isnan(Du->data[posOS]) ){
					AU->data[posOS] = M->data[posOS]*dV->data[pos] - wNeighU + (Du->data[posOS]+float_temp1)*dU->data[pos];
				}else{
					AU->data[posOS] = - wNeighU + float_temp1*dU->data[pos];
				}
				if( !isnan(Dv->data[posOS]) ){
					AV->data[posOS] = M->data[posOS]*dU->data[pos] - wNeighV + (Dv->data[posOS]+float_temp1)*dV->data[pos];
				}else{
					AV->data[posOS] = - wNeighV + float_temp1*dV->data[pos];
				}
			}
		}
	}
	/*
	-------------------
	--Fill the borders
	-------------------
	*/
	for(k=0;k<nframes;k++)
	{
		frameoffset = k*nrows*ncols;
		for(j=0;j<ncols;j++)
		{
			posOS = j*nrows + frameoffset;
	
			AU->data[ posOS ] = AU->data[ posOS+1 ];
			AU->data[ posOS+nrows-1 ] = AU->data[ posOS+nrows-2 ];
	
			AV->data[ posOS ] = AV->data[ pos+1 ];
			AV->data[ posOS+nrows-1 ] = AV->data[ posOS+nrows-2 ];
		}
	
		for(i=0;i<nrows;i++)
		{
			posOS = i + frameoffset;
			AU->data[ posOS ] = AU->data[ posOS+nrows ];
			AU->data[ posOS+(ncols-1)*nrows ] = AU->data[ posOS+(ncols-2)*nrows ];
	
			AV->data[ posOS ] = AV->data[ posOS+nrows ];
			AV->data[ posOS+(ncols-1)*nrows ] = AV->data[ posOS+(ncols-2)*nrows ];
		}
	}
}

/*
//----------------------------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel successive over relaxation with 2D (spatial) regularization for symmetric case
//---Late linearization of constancy assumptions with symmetry.
//----------------------------------------------------------------------------------------------------------
*/
void GS_SOR_llinsym4_2d(	struct matrixM *U0, 
				struct matrixM *V0, 
				struct matrixM *dU0,
				struct matrixM *dV0,
				struct matrixM *M0,
				struct matrixM *Cu0,
				struct matrixM *Cv0,
				struct matrixM *Du0,
				struct matrixM *Dv0,
				struct matrixM *wW0,
				struct matrixM *wN0,
				struct matrixM *wE0,
				struct matrixM *wS0,
				struct matrixM *U1, 
				struct matrixM *V1, 
				struct matrixM *dU1,
				struct matrixM *dV1,
				struct matrixM *M1,
				struct matrixM *Cu1,
				struct matrixM *Cv1,
				struct matrixM *Du1,
				struct matrixM *Dv1,
				struct matrixM *wW1,
				struct matrixM *wN1,
				struct matrixM *wE1,
				struct matrixM *wS1,
				struct mparam Mparam)
{
	float wNeighU0, wNeighV0, wNeighU1, wNeighV1, float_temp1, float_temp2, float_temp3;
	float dU_newapprox0, dV_newapprox0;
	float dU_newapprox1, dV_newapprox1;
	float *u_data_div0, *v_data_div0;
	float *u_data_div1, *v_data_div1;
	unsigned int i,j,k,pos,iter;
	unsigned int nrows = M0->dimElems[0];
	unsigned int ncols = M0->dimElems[1];
	
	
	/*Reserve space for the divisor terms*/
	if( (u_data_div0 = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym4_2d: error reserving space for 'u0_data_div'\n");
		return;
	}
	if( (v_data_div0 = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym4_2d: error reserving space for 'v0_data_div'\n");
		free(u_data_div0);
		return;
	}
	if( (u_data_div1 = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym4_2d: error reserving space for 'u1_data_div'\n");
		free(u_data_div0);
		free(v_data_div0);
		return;
	}
	if( (v_data_div1 = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llinsym4_2d: error reserving space for 'v1_data_div'\n");
		free(u_data_div0);
		free(v_data_div0);
		free(u_data_div1);
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
			/*wpos = pos - nrows;
			epos = pos + nrows;
			npos = pos - 1;
			spos = pos + 1;*/

			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			/*wNeighU0 =	(dU0->data[ wpos ] + U0->data[ wpos ] - U0->data[ pos ])*wW0->data[ pos ] + 
					(dU0->data[ epos ] + U0->data[ epos ] - U0->data[ pos ])*wE0->data[ pos ] +
					(dU0->data[ npos ] + U0->data[ npos ] - U0->data[ pos ])*wN0->data[ pos ] + 
					(dU0->data[ spos ] + U0->data[ spos ] - U0->data[ pos ])*wS0->data[ pos ];*/
			wNeighU0	= dU0->data[ pos - nrows ] + U0->data[ pos - nrows ];
			float_temp1	= dU0->data[ pos + nrows ] + U0->data[ pos + nrows ];
			float_temp2	= dU0->data[ pos - 1 ] + U0->data[ pos - 1 ];
			float_temp3	= dU0->data[ pos + 1 ] + U0->data[ pos + 1 ];
				
			wNeighU0	-= U0->data[ pos ];
			float_temp1	-= U0->data[ pos ];
			float_temp2	-= U0->data[ pos ];
			float_temp3	-= U0->data[ pos ];

			wNeighU0	*= wW0->data[ pos ];
			float_temp1	*= wE0->data[ pos ];
			float_temp2	*= wN0->data[ pos ];
			float_temp3	*= wS0->data[ pos ];

			wNeighU0	+= float_temp1;
			float_temp2	+= float_temp3;
			wNeighU0	+= float_temp2;
				
			/*wNeighV0 =	(dV0->data[ wpos ] + V0->data[ wpos ] - V0->data[ pos ])*wW0->data[ pos ] + 
					(dV0->data[ epos ] + V0->data[ epos ] - V0->data[ pos ])*wE0->data[ pos ] +
					(dV0->data[ npos ] + V0->data[ npos ] - V0->data[ pos ])*wN0->data[ pos ] + 
					(dV0->data[ spos ] + V0->data[ spos ] - V0->data[ pos ])*wS0->data[ pos ];*/
			wNeighV0	= dV0->data[ pos - nrows ] + V0->data[ pos - nrows ];
			float_temp1	= dV0->data[ pos + nrows ] + V0->data[ pos + nrows ];
			float_temp2	= dV0->data[ pos - 1 ] + V0->data[ pos - 1 ];
			float_temp3	= dV0->data[ pos + 1 ] + V0->data[ pos + 1 ];
				
			wNeighV0	-= V0->data[ pos ];
			float_temp1	-= V0->data[ pos ];
			float_temp2	-= V0->data[ pos ];
			float_temp3	-= V0->data[ pos ];

			wNeighV0	*= wW0->data[ pos ];
			float_temp1	*= wE0->data[ pos ];
			float_temp2	*= wN0->data[ pos ];
			float_temp3	*= wS0->data[ pos ];

			wNeighV0	+= float_temp1;
			float_temp2	+= float_temp3;
			wNeighV0	+= float_temp2;
				
			/*wNeighHorU1 =	(dU1->data[ wpos ] + U1->data[ wpos ] - U1->data[ pos ])*wW1->data[ pos ] + 
					(dU1->data[ epos ] + U1->data[ epos ] - U1->data[ pos ])*wE1->data[ pos ] +
					(dU1->data[ npos ] + U1->data[ npos ] - U1->data[ pos ])*wN1->data[ pos ] + 
					(dU1->data[ spos ] + U1->data[ spos ] - U1->data[ pos ])*wS1->data[ pos ];*/
			wNeighU1	= dU1->data[ pos - nrows ] + U1->data[ pos - nrows ];
			float_temp1	= dU1->data[ pos + nrows ] + U1->data[ pos + nrows ];
			float_temp2	= dU1->data[ pos - 1 ] + U1->data[ pos - 1 ];
			float_temp3	= dU1->data[ pos + 1 ] + U1->data[ pos + 1 ];
				
			wNeighU1	-= U1->data[ pos ];
			float_temp1	-= U1->data[ pos ];
			float_temp2	-= U1->data[ pos ];
			float_temp3	-= U1->data[ pos ];

			wNeighU1	*= wW1->data[ pos ];
			float_temp1	*= wE1->data[ pos ];
			float_temp2	*= wN1->data[ pos ];
			float_temp3	*= wS1->data[ pos ];

			wNeighU1	+= float_temp1;
			float_temp2	+= float_temp3;
			wNeighU1	+= float_temp2;

			/*wNeighHorV1 =	(dV1->data[ wpos ] + V1->data[ wpos ] - V1->data[ pos ])*wW1->data[ pos ] + 
					(dV1->data[ epos ] + V1->data[ epos ] - V1->data[ pos ])*wE1->data[ pos ] +
					(dV1->data[ npos ] + V1->data[ npos ] - V1->data[ pos ])*wN1->data[ pos ] + 
					(dV1->data[ spos ] + V1->data[ spos ] - V1->data[ pos ])*wS1->data[ pos ];*/
			wNeighV1	= dV1->data[ pos - nrows ] + V1->data[ pos - nrows ];
			float_temp1	= dV1->data[ pos + nrows ] + V1->data[ pos + nrows ];
			float_temp2	= dV1->data[ pos - 1 ] + V1->data[ pos - 1 ];
			float_temp3	= dV1->data[ pos + 1 ] + V1->data[ pos + 1 ];
				
			wNeighV1	-= V1->data[ pos ];
			float_temp1	-= V1->data[ pos ];
			float_temp2	-= V1->data[ pos ];
			float_temp3	-= V1->data[ pos ];

			wNeighV1	*= wW1->data[ pos ];
			float_temp1	*= wE1->data[ pos ];
			float_temp2	*= wN1->data[ pos ];
			float_temp3	*= wS1->data[ pos ];

			wNeighV1	+= float_temp1;
			float_temp2	+= float_temp3;
			wNeighV1	+= float_temp2;
	
			if(iter==0)
			{
				/*Sum of weights*/
				float_temp1 = wW0->data[ pos ] + wE0->data[ pos ];
				float_temp2 = wN0->data[ pos ] + wS0->data[ pos ];
				float_temp1 += float_temp2;

				if(isnan(Du0->data[pos]))
					u_data_div0[pos] = 1.0f/( float_temp1 );
				else
					u_data_div0[pos] = 1.0f/( float_temp1 + Du0->data[pos] );

				if(isnan(Dv0->data[pos]))
					v_data_div0[pos] = 1.0f/( float_temp1 );
				else
					v_data_div0[pos] = 1.0f/( float_temp1 + Dv0->data[pos] );

				/*Sum of weights*/
				float_temp1 = wW1->data[ pos ] + wE1->data[ pos ];
				float_temp2 = wN1->data[ pos ] + wS1->data[ pos ];
				float_temp1 = float_temp1 + float_temp2;

				if(isnan(Du1->data[pos]))
					u_data_div1[pos] = 1.0f/( float_temp1 );
				else
					u_data_div1[pos] = 1.0f/( float_temp1 + Du1->data[pos] );

				if(isnan(Dv1->data[pos]))
					v_data_div1[pos] = 1.0f/( float_temp1 );
				else
					v_data_div1[pos] = 1.0f/( float_temp1 + Dv1->data[pos] );
			}
	
			if(isnan(Cu0->data[pos]))
				dU_newapprox0 = wNeighU0*u_data_div0[pos];
			else
				dU_newapprox0 = ( wNeighU0 + M0->data[pos]*dV0->data[pos] + Cu0->data[pos] )*u_data_div0[pos];

			if(isnan(Cv0->data[pos]))
				dV_newapprox0 = wNeighV0*v_data_div0[pos];
			else
				dV_newapprox0 = ( wNeighV0 - M0->data[pos]*dU0->data[pos] + Cv0->data[pos] )*v_data_div0[pos];

			if(isnan(Cu1->data[pos]))
				dU_newapprox1 = wNeighU1*u_data_div1[pos];
			else
				dU_newapprox1 = ( wNeighU1 - M1->data[pos]*dV1->data[pos] + Cu1->data[pos] )*u_data_div1[pos];

			if(isnan(Cv1->data[pos]))
				dV_newapprox1 = wNeighV1*v_data_div1[pos];
			else
				dV_newapprox1 = ( wNeighV1 + M1->data[pos]*dU1->data[pos]+ Cv1->data[pos] )*v_data_div1[pos];

			if(!isnan( dU_newapprox0 ))
				dU0->data[ pos ] = (1.0f-Mparam.omega)*dU0->data[ pos ] + Mparam.omega*dU_newapprox0;
			if(!isnan( dV_newapprox0 ))
				dV0->data[ pos ] = (1.0f-Mparam.omega)*dV0->data[ pos ] + Mparam.omega*dV_newapprox0;
			if(!isnan( dU_newapprox1 ))
				dU1->data[ pos ] = (1.0f-Mparam.omega)*dU1->data[ pos ] + Mparam.omega*dU_newapprox1;
			if(!isnan( dV_newapprox1 ))
				dV1->data[ pos ] = (1.0f-Mparam.omega)*dV1->data[ pos ] + Mparam.omega*dV_newapprox1;
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
	
			dV0->data[ pos ] = dV0->data[ pos+1 ];
			dV0->data[ pos+nrows-1 ] = dV0->data[ pos+nrows-2 ];
	
			dU1->data[ pos ] = dU1->data[ pos+1 ];
			dU1->data[ pos+nrows-1 ] = dU1->data[ pos+nrows-2 ];
	
			dV1->data[ pos ] = dV1->data[ pos+1 ];
			dV1->data[ pos+nrows-1 ] = dV1->data[ pos+nrows-2 ];
		}
	
		for(i=0;i<nrows;i++)
		{
			dU0->data[ i ] = dU0->data[ i+nrows ];
			dU0->data[ i+(ncols-1)*nrows ] = dU0->data[ i+(ncols-2)*nrows ];
	
			dV0->data[ i ] = dV0->data[ i+nrows ];
			dV0->data[ i+(ncols-1)*nrows ] = dV0->data[ i+(ncols-2)*nrows ];
	
			dU1->data[ i ] = dU1->data[ i+nrows ];
			dU1->data[ i+(ncols-1)*nrows ] = dU1->data[ i+(ncols-2)*nrows ];
	
			dV1->data[ i ] = dV1->data[ i+nrows ];
			dV1->data[ i+(ncols-1)*nrows ] = dV1->data[ i+(ncols-2)*nrows ];
		}

	}
	
	free(u_data_div0);
	free(v_data_div0);
	free(u_data_div1);
	free(v_data_div1);
}

/*
//---------------------------------------------------------------------------------------
//---Alternating line relaxation Gauss-Seidel (Gauss-Seidel block variant)
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization.
//---Late linearization of constancy assumptions imposing symmetry
//---------------------------------------------------------------------------------------
*/
void GS_ALR_SOR_llinsym4_2d(	struct matrixM *U0, 
				struct matrixM *V0, 
				struct matrixM *dU0,
				struct matrixM *dV0,
				struct matrixM *M0,
				struct matrixM *Cu0,
				struct matrixM *Cv0,
				struct matrixM *Du0,
				struct matrixM *Dv0,
				struct matrixM *wW0,
				struct matrixM *wN0,
				struct matrixM *wE0,
				struct matrixM *wS0,
				struct matrixM *U1, 
				struct matrixM *V1, 
				struct matrixM *dU1,
				struct matrixM *dV1,
				struct matrixM *M1,
				struct matrixM *Cu1,
				struct matrixM *Cv1,
				struct matrixM *Du1,
				struct matrixM *Dv1,
				struct matrixM *wW1,
				struct matrixM *wN1,
				struct matrixM *wE1,
				struct matrixM *wS1,
				struct mparam Mparam)
{
	float wNeigh, *Cp, *Dp;
	int iter;
	int nrows = Cu0->dimElems[0];
	int ncols = Cu0->dimElems[1];
	int maxLength;
	int iterations = (int)Mparam.iter;

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
		westColumn_llin4(	U0, V0, dU0, dV0, M0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_llin4(	U0, V0, dU0, dV0, M0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_llin4(	U0, V0, dU0, dV0, M0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);

		westColumn_llin4(	V0, U0, dV0, dU0, M0, Cv0, Dv0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_llin4(	V0, U0, dV0, dU0, M0, Cv0, Dv0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_llin4(	V0, U0, dV0, dU0, M0, Cv0, Dv0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		northRow_llin4(	U0, V0, dU0, dV0, M0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		middleRow_llin4(U0, V0, dU0, dV0, M0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		southRow_llin4(	U0, V0, dU0, dV0, M0, Cu0, Du0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);

		northRow_llin4(	V0, U0, dV0, dU0, M0, Cv0, Dv0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		middleRow_llin4(V0, U0, dV0, dU0, M0, Cv0, Dv0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);
		southRow_llin4(	V0, U0, dV0, dU0, M0, Cv0, Dv0, wW0, wN0, wE0, wS0, Mparam, Cp, Dp, ncols, nrows);

		/*
		------------------------
		-Vertical line relaxing
		------------------------
		*/
		westColumn_llin4(	U1, V1, dU1, dV1, M1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_llin4(	U1, V1, dU1, dV1, M1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_llin4(	U1, V1, dU1, dV1, M1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);

		westColumn_llin4(	V1, U1, dV1, dU1, M1, Cv1, Dv1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_llin4(	V1, U1, dV1, dU1, M1, Cv1, Dv1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_llin4(	V1, U1, dV1, dU1, M1, Cv1, Dv1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		northRow_llin4(	U1, V1, dU1, dV1, M1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		middleRow_llin4(U1, V1, dU1, dV1, M1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		southRow_llin4(	U1, V1, dU1, dV1, M1, Cu1, Du1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);

		northRow_llin4(	V1, U1, dV1, dU1, M1, Cv1, Dv1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		middleRow_llin4(V1, U1, dV1, dU1, M1, Cv1, Dv1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
		southRow_llin4(	V1, U1, dV1, dU1, M1, Cv1, Dv1, wW1, wN1, wE1, wS1, Mparam, Cp, Dp, ncols, nrows);
	
	}

	free( Cp );
	free( Dp );
}

/*
//---------------------------------------------------------------------------------------
//---Point wise Gauss-Seidel successive over relaxation with 2D (spatial) regularization
//---Late linearization of constancy assumptions.
//---------------------------------------------------------------------------------------
*/
void GS_SOR_llin8_2d(	struct matrixM *U, 
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
			struct mparam Mparam)
{
	float wNeighU, wNeighV, float_temp1, float_temp2, float_temp3;
	float dU_newapprox, dV_newapprox;
	float *u_data_div, *v_data_div;
	unsigned int i,j,k,pos,iter;
	unsigned int nrows = M->dimElems[0];
	unsigned int ncols = M->dimElems[1];
		
	/*Reserve space for the divisor terms*/
	if( (u_data_div = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin8_2d: error reserving space for 'u_data_div'\n");
		return;
	}
	if( (v_data_div = (float*)malloc( nrows*ncols*sizeof(float) ))==NULL )
	{
		printf("GS_SOR_llin8_2d: error reserving space for 'v_data_div'\n");
		free(v_data_div);
		return;
	}

	for(iter=0;iter<(unsigned int)Mparam.iter;iter++)
	{
		for(j=1;j<ncols-1;j++)
		{
		for(i=1;i<nrows-1;i++)		
		{
			/*Spatial position indices*/
			pos = j*nrows+i;
			/*wpos = pos - nrows;
			epos = pos + nrows;
			npos = pos - 1;
			spos = pos + 1;*/

			/*
			---------------------
			--Spatial diffusion--
			---------------------	
			*/
			/*
			wNeighU =	(dU->data[ wpos ] + U->data[ wpos ] - U->data[ pos ])*wW->data[ pos ] + 
					(dU->data[ epos ] + U->data[ epos ] - U->data[ pos ])*wE->data[ pos ] +
					(dU->data[ npos ] + U->data[ npos ] - U->data[ pos ])*wN->data[ pos ] + 
					(dU->data[ spos ] + U->data[ spos ] - U->data[ pos ])*wS->data[ pos ];
			*/
			wNeighU		= dU->data[ pos - nrows ] + U->data[ pos - nrows ];
			float_temp1	= dU->data[ pos + nrows ] + U->data[ pos + nrows ];
			float_temp2	= dU->data[ pos - 1 ] + U->data[ pos - 1 ];
			float_temp3	= dU->data[ pos + 1 ] + U->data[ pos + 1 ];
				
			wNeighU		-= U->data[ pos ];
			float_temp1	-= U->data[ pos ];
			float_temp2	-= U->data[ pos ];
			float_temp3	-= U->data[ pos ];

			wNeighU		*= wW->data[ pos ];
			float_temp1	*= wE->data[ pos ];
			float_temp2	*= wN->data[ pos ];
			float_temp3	*= wS->data[ pos ];

			wNeighU		+= float_temp1;
			float_temp2	+= float_temp3;
			wNeighU		+= float_temp2;

			/*wNeighV =	(dV->data[ wpos ] + V->data[ wpos ] - V->data[ pos ])*wW->data[ pos ] + 
					(dV->data[ epos ] + V->data[ epos ] - V->data[ pos ])*wE->data[ pos ] +
					(dV->data[ npos ] + V->data[ npos ] - V->data[ pos ])*wN->data[ pos ] + 
					(dV->data[ spos ] + V->data[ spos ] - V->data[ pos ])*wS->data[ pos ];
			*/
			wNeighV		= dV->data[ pos - nrows ] + V->data[ pos - nrows ];
			float_temp1	= dV->data[ pos + nrows ] + V->data[ pos + nrows ];
			float_temp2	= dV->data[ pos - 1 ] + V->data[ pos - 1 ];
			float_temp3	= dV->data[ pos + 1 ] + V->data[ pos + 1 ];
				
			wNeighV		-= V->data[ pos ];
			float_temp1	-= V->data[ pos ];
			float_temp2	-= V->data[ pos ];
			float_temp3	-= V->data[ pos ];

			wNeighV		*= wW->data[ pos ];
			float_temp1	*= wE->data[ pos ];
			float_temp2	*= wN->data[ pos ];
			float_temp3	*= wS->data[ pos ];

			wNeighV		+= float_temp1;
			float_temp2	+= float_temp3;
			wNeighV		+= float_temp2;

			if(iter==0)
			{
				/*Sum of weights*/
				float_temp1 = wW->data[ pos ] + wE->data[ pos ];
				float_temp2 = wN->data[ pos ] + wS->data[ pos ];
				float_temp1 += float_temp2;

				if(isnan(Du->data[pos]))
					u_data_div[pos] = 1.0f/( float_temp1 );
				else
					u_data_div[pos] = 1.0f/( float_temp1 + Du->data[pos] );

				if(isnan(Dv->data[pos]))
					v_data_div[pos] = 1.0f/( float_temp1 );
				else
					v_data_div[pos] = 1.0f/( float_temp1 + Dv->data[pos] );
			}

			if(isnan(Cu->data[pos]))
				dU_newapprox = wNeighU*u_data_div[pos];
			else{
				float_temp1 = wNeighU + Cu->data[pos];
				float_temp2 = M->data[pos]*dV->data[pos];
				float_temp1 -= float_temp2;
				dU_newapprox = float_temp1 * u_data_div[pos];
								
				/*dU_newapprox = ( wNeighU - M->data[pos]*dV->data[pos] + Cu->data[pos] )*u_data_div[pos];*/
			}
				
			if(isnan(Cv->data[pos]))
				dV_newapprox = wNeighV*v_data_div[pos];
			else{
				float_temp1 = wNeighV + Cv->data[pos];
				float_temp2 = M->data[pos]*dU->data[pos];
				float_temp1 -= float_temp2;
				dV_newapprox = float_temp1 * v_data_div[pos];
				
				/*dV_newapprox = ( wNeighV - M->data[pos]*dU->data[pos] + Cv->data[pos] )*v_data_div[pos];*/
			}

			dU->data[ pos ] = (1.0f-Mparam.omega)*dU->data[ pos ] + Mparam.omega*dU_newapprox;
			dV->data[ pos ] = (1.0f-Mparam.omega)*dV->data[ pos ] + Mparam.omega*dV_newapprox;
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
	
			dV->data[ pos ] = dV->data[ pos+1 ];
			dV->data[ pos+nrows-1 ] = dV->data[ pos+nrows-2 ];
		}
	
		for(i=0;i<nrows;i++)
		{
			dU->data[ i ] = dU->data[ i+nrows ];
			dU->data[ i+(ncols-1)*nrows ] = dU->data[ i+(ncols-2)*nrows ];
	
			dV->data[ i ] = dV->data[ i+nrows ];
			dV->data[ i+(ncols-1)*nrows ] = dV->data[ i+(ncols-2)*nrows ];
		}
	}
	
	
	free(u_data_div);
	free(v_data_div);
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
				struct mparam Mparam)
{
	float *Cp, *Dp;
	int iter;
	int nrows = Cu->dimElems[0];
	int ncols = Cu->dimElems[1];
	int maxLength;
	int iterations = (int)Mparam.iter;

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
		westColumn_llin8(	U, V, dU, dV, M, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_llin8(	U, V, dU, dV, M, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_llin8(	U, V, dU, dV, M, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);

		westColumn_llin8(	V, U, dV, dU, M, Cv, Dv, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		middleColumn_llin8(	V, U, dV, dU, M, Cv, Dv, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		eastColumn_llin8(	V, U, dV, dU, M, Cv, Dv, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);

		/*
		--------------------------
		-Horizontal line relaxing
		--------------------------
		*/
		northRow_llin8(	V, U, dV, dU, M, Cv, Dv, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		middleRow_llin8(V, U, dV, dU, M, Cv, Dv, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		southRow_llin8(	V, U, dV, dU, M, Cv, Dv, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		
		northRow_llin8(	U, V, dU, dV, M, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		middleRow_llin8(U, V, dU, dV, M, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);
		southRow_llin8(	U, V, dU, dV, M, Cu, Du, wW, wNW, wN, wNE, wE, wSE, wS, wSW, Mparam, Cp, Dp, ncols, nrows);

	}

	free( Cp );
	free( Dp );
}

/*
//*********************************************************************************
// LINE SOLVERS FOR EARLY LINEARIZATION
//*********************************************************************************
*/

/*
//------------------------------------------------------
// TDMA vertical line solving - west most columns solver
//------------------------------------------------------
*/
void westColumn_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
	i = 0; j=0;
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	/*a = 	0.0f;*/
	b = 	wS->data[pos] + wE->data[pos];
	c = 	-wS->data[pos];
	d =	wE->data[pos]*U->data[epos];
	
	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*V->data[pos];
	}
	
	cp[i] = c/b;
	dp[i] = d/b;

	/*---Middle rows---*/
	for(i=1;i<=nrows-2;i++)
	{
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		a = 	-wN->data[pos];
		b = 	wN->data[pos] + wS->data[pos] + wE->data[pos];
		c = 	-wS->data[pos];
		d =	wE->data[pos]*U->data[epos];
			
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*V->data[pos];
		}
		div = 1 / (b-cp[i-1]*a);
		cp[i] = c*div;
		dp[i] = (d-dp[i-1]*a)*div;
	}

	/*---Last row---*/
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	a = 	-wN->data[pos];
	b = 	wN->data[pos] + wE->data[pos];
	/*c = 	0.0f;*/
	d =	wE->data[pos]*U->data[epos];

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*V->data[pos];
	}

	dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
	temp1 = U->data[pos];
	U->data[pos] = dp[i];

	for(i=nrows-2;i>=0;i--)
	{
		pos = j*nrows+i;
		temp2 = U->data[pos];
		U->data[pos] = dp[i] - cp[i]*U->data[pos+1];
		
		U->data[pos+1] = Mparam.omega*U->data[pos+1] + (1.0f-Mparam.omega)*temp1;
		temp1 = temp2;
		
	}
	U->data[pos] = Mparam.omega*U->data[pos] + (1.0f-Mparam.omega)*temp1;
}

/*
//------------------------------------------------------
// TDMA vertical line solving - middle columns solver
//------------------------------------------------------
*/
void middleColumn_elin4(struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
		d =	wW->data[pos]*U->data[wpos] + wE->data[pos]*U->data[epos];
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*V->data[pos];
		}
		
		cp[i] = c/b;
		dp[i] = d/b;
	
		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos = j*nrows+i;
			wpos = pos-nrows;
			epos = pos+nrows;
			npos = pos-1;
			spos = pos+1;
	
			a = 	-wN->data[pos];
			b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos];
			c = 	-wS->data[pos];
			d =	wW->data[pos]*U->data[wpos] + wE->data[pos]*U->data[epos];
				
			if(!isnan(Cu->data[pos]))
			{
				b += Du->data[pos];
				d += Cu->data[pos];
				d -= M->data[pos]*V->data[pos];
			}
			div = 1 / (b-cp[i-1]*a);
			cp[i] = c*div;
			dp[i] = (d-dp[i-1]*a)*div;
		}
	
		/*---Last row---*/
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;
	
		a = 	-wN->data[pos];
		b = 	wN->data[pos] + wE->data[pos] + wW->data[pos];
		/*c = 	0.0f;*/
		d =	wW->data[pos]*U->data[wpos] + wE->data[pos]*U->data[epos];

		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*V->data[pos];
		}

		dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
		temp1 = U->data[pos];
		U->data[pos] = dp[i];
	
		for(i=nrows-2;i>=0;i--)
		{
			pos = j*nrows+i;
			temp2 = U->data[pos];
			U->data[pos] = dp[i] - cp[i]*U->data[pos+1];

			U->data[pos+1] = Mparam.omega*U->data[pos+1] + (1.0f-Mparam.omega)*temp1;
			temp1 = temp2;
		}
		U->data[pos] = Mparam.omega*U->data[pos] + (1.0f-Mparam.omega)*temp1;

	}
}

/*
//------------------------------------------------------
// TDMA vertical line solving - east most column solver
//------------------------------------------------------
*/
void eastColumn_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
	i = 0; j=ncols-1;
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	/*a = 	0.0f;*/
	b = 	wS->data[pos] + wW->data[pos];
	c = 	-wS->data[pos];
	d =	wW->data[pos]*U->data[wpos];
	
	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*V->data[pos];
	}
	
	cp[i] = c/b;
	dp[i] = d/b;

	/*---Middle rows---*/
	for(i=1;i<=nrows-2;i++)
	{
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		a = 	-wN->data[pos];
		b = 	wN->data[pos] + wS->data[pos] + wW->data[pos];
		c = 	-wS->data[pos];
		d =	wW->data[pos]*U->data[wpos];
			
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*V->data[pos];
		}
		div = 1 / (b-cp[i-1]*a);
		cp[i] = c*div;
		dp[i] = (d-dp[i-1]*a)*div;
	}

	/*---Last row---*/
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	a = 	-wN->data[pos];
	b = 	wN->data[pos] + wW->data[pos];
	/*c = 	0.0f;*/
	d =	wW->data[pos]*U->data[wpos];

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*V->data[pos];
	}

	dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
	temp1 = U->data[pos];
	U->data[pos] = dp[i];

	for(i=nrows-2;i>=0;i--)
	{
		pos = j*nrows+i;
		temp2 = U->data[pos];
		U->data[pos] = dp[i] - cp[i]*U->data[pos+1];

		U->data[pos+1] = Mparam.omega*U->data[pos+1] + (1.0f-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	U->data[pos] = Mparam.omega*U->data[pos] + (1.0f-Mparam.omega)*temp1;

}

/*
//------------------------------------------------------
// TDMA horizontal line solving - north most row solver
//------------------------------------------------------
*/
void northRow_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
	j = 0; i=0;
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	/*a =	0.0f;*/
	b =	wS->data[pos] + wE->data[pos];
	c =	-wE->data[pos];
	d =	wS->data[pos]*U->data[spos];

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*V->data[pos];
	}
	
	cp[j] = c/b;
	dp[j] = d/b;

	/*---Middle columns---*/
	for(j=1;j<=ncols-2;j++)
	{
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		a = 	-wW->data[pos];	
		b = 	wS->data[pos] + wE->data[pos] + wW->data[pos];
		c = 	-wE->data[pos];
		d =	wS->data[pos]*U->data[spos];
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*V->data[pos];
		}
		div = 1 / (b-cp[j-1]*a);
		cp[j] = c*div;
		dp[j] = (d-dp[j-1]*a)*div;
	}

	/*---Last column---*/
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	a = 	-wW->data[pos];	
	b = 	wS->data[pos] + wW->data[pos];
	/*c = 	0.0f;*/
	d =	wS->data[pos]*U->data[spos];

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*V->data[pos];
	}

	dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
	temp1 = U->data[pos];
	U->data[pos] = dp[j];

	for(j=ncols-2;j>=0;j--)
	{
		pos = j*nrows+i;
		temp2 = U->data[pos];
		U->data[pos] = dp[j] - cp[j]*U->data[pos+nrows];

		U->data[pos+nrows] = Mparam.omega*U->data[pos+nrows] + (1.0f-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	U->data[pos] = Mparam.omega*U->data[pos] + (1.0f-Mparam.omega)*temp1;

}

/*
//------------------------------------------------------
// TDMA horizontal line solving - middle rows solver
//------------------------------------------------------
*/
void middleRow_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
		d =	wS->data[pos]*U->data[spos] + wN->data[pos]*U->data[npos];

		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*V->data[pos];
		}
		
		cp[j] = c/b;
		dp[j] = d/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos = j*nrows+i;
			wpos = pos-nrows;
			epos = pos+nrows;
			npos = pos-1;
			spos = pos+1;

			a = 	-wW->data[pos];	
			b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos];
			c = 	-wE->data[pos];
			d =	wS->data[pos]*U->data[spos] + wN->data[pos]*U->data[npos];
			
			if(!isnan(Cu->data[pos]))
			{
				b += Du->data[pos];
				d += Cu->data[pos];
				d -= M->data[pos]*V->data[pos];
			}
			div = 1 / (b-cp[j-1]*a);
			cp[j] = c*div;
			dp[j] = (d-dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		a = 	-wW->data[pos];	
		b = 	wN->data[pos] + wS->data[pos] + wW->data[pos];
		/*c = 	0.0f;*/
		d =	wS->data[pos]*U->data[spos] + wN->data[pos]*U->data[npos];
	
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*V->data[pos];
		}

		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
		temp1 = U->data[pos];
		U->data[pos] = dp[j];

		for(j=ncols-2;j>=0;j--)
		{
			pos = j*nrows+i;
			temp2 = U->data[pos];
			U->data[pos] = dp[j] - cp[j]*U->data[pos+nrows];

			U->data[pos+nrows] = Mparam.omega*U->data[pos+nrows] + (1.0f-Mparam.omega)*temp1;
			temp1 = temp2;
		}
		U->data[pos] = Mparam.omega*U->data[pos] + (1.0f-Mparam.omega)*temp1;
	}
}

/*
//------------------------------------------------------
// TDMA horizontal line solving - south most rows solver
//------------------------------------------------------
*/
void southRow_elin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *M,
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
	j = 0; i = nrows-1;
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	/*a =	0.0f;*/
	b =	wN->data[pos] + wE->data[pos];
	c =	-wE->data[pos];
	d =	wN->data[pos]*U->data[npos];

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*V->data[pos];
	}
	
	cp[j] = c/b;
	dp[j] = d/b;

	/*---Middle columns---*/
	for(j=1;j<=ncols-2;j++)
	{
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

		a = 	-wW->data[pos];	
		b = 	wN->data[pos] + wE->data[pos] + wW->data[pos];
		c = 	-wE->data[pos];
		d =	wN->data[pos]*U->data[npos];
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*V->data[pos];
		}
		div = 1 / (b-cp[j-1]*a);
		cp[j] = c*div;
		dp[j] = (d-dp[j-1]*a)*div;
	}

	/*---Last column---*/
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	a = 	-wW->data[pos];	
	b = 	wN->data[pos] + wW->data[pos];
	/*c = 	0.0f;*/
	d =	wN->data[pos]*U->data[npos];

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*V->data[pos];
	}

	dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
	temp1 = U->data[pos];
	U->data[pos] = dp[j];

	for(j=ncols-2;j>=0;j--)
	{
		pos = j*nrows+i;
		temp2 = U->data[pos];
		U->data[pos] = dp[j] - cp[j]*U->data[pos+nrows];

		U->data[pos+nrows] = Mparam.omega*U->data[pos+nrows] + (1.0f-Mparam.omega)*temp1;
		temp1 = temp2;
	}
	U->data[pos] = Mparam.omega*U->data[pos] + (1.0f-Mparam.omega)*temp1;
}

/*
//*********************************************************************************
// LINE SOLVERS FOR LATE LINEARIZATION CASE
//*********************************************************************************
*/

/*
//------------------------------------------------------
// TDMA vertical line solving - west most columns solver
//------------------------------------------------------
*/
void westColumn_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
		d -= M->data[pos]*dV->data[pos];
	}
	
	cp[i] = c/b;
	dp[i] = d/b;

	/*---Middle rows---*/
	for(i=1;i<=nrows-2;i++)
	{
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

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
			d -= M->data[pos]*dV->data[pos];
		}
		div = 1 / (b-cp[i-1]*a);
		cp[i] = c*div;
		dp[i] = (d-dp[i-1]*a)*div;
	}

	/*---Last row---*/
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	a = 	-wN->data[pos];
	b = 	wN->data[pos] + wE->data[pos];
	/*c = 	0.0f;*/
	d =	  wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]);
	
	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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
void middleColumn_llin4(struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
			d -= M->data[pos]*dV->data[pos];
		}
		
		cp[i] = c/b;
		dp[i] = d/b;
	
		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos = j*nrows+i;
			wpos = pos-nrows;
			epos = pos+nrows;
			npos = pos-1;
			spos = pos+1;
	
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
				d -= M->data[pos]*dV->data[pos];
			}
			div = 1 / (b-cp[i-1]*a);
			cp[i] = c*div;
			dp[i] = (d-dp[i-1]*a)*div;
		}
	
		/*---Last row---*/
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;
	
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
			d -= M->data[pos]*dV->data[pos];
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
void eastColumn_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
		d -= M->data[pos]*dV->data[pos];
	}

	cp[i] = c/b;
	dp[i] = d/b;

	/*---Middle rows---*/
	for(i=1;i<=nrows-2;i++)
	{
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

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
			d -= M->data[pos]*dV->data[pos];
		}

		div   =	1 / (b-cp[i-1]*a);
		cp[i] = c*div;
		dp[i] = (d-dp[i-1]*a)*div;
	}

	/*---Last row---*/
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	a = 	-wN->data[pos];
	b = 	wN->data[pos] + wW->data[pos];
	/*c = 	0.0f;*/
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]);
	
	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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
// TDMA horizontal line solving - north most row solver
//------------------------------------------------------
*/
void northRow_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
	d =	wE->data[pos]*(U->data[epos]-U->data[pos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
	}
	
	cp[j] = c/b;
	dp[j] = d/b;

	/*---Middle columns---*/
	for(j=1;j<=ncols-2;j++)
	{
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

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
			d -= M->data[pos]*dV->data[pos];
		}
		div = 1 / (b-cp[j-1]*a);
		cp[j] = c*div;
		dp[j] = (d-dp[j-1]*a)*div;
	}

	/*---Last column---*/
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	a = 	-wW->data[pos];	
	b = 	wS->data[pos] + wW->data[pos];
	/*c = 	0.0f;*/
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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
void middleRow_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
			d -= M->data[pos]*dV->data[pos];
		}
		
		cp[j] = c/b;
		dp[j] = d/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos = j*nrows+i;
			wpos = pos-nrows;
			epos = pos+nrows;
			npos = pos-1;
			spos = pos+1;

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
				d -= M->data[pos]*dV->data[pos];
			}
			div = 1 / (b-cp[j-1]*a);
			cp[j] = c*div;
			dp[j] = (d-dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

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
			d -= M->data[pos]*dV->data[pos];
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
// TDMA horizontal line solving - south most rows solver
//------------------------------------------------------
*/
void southRow_llin4(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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

	/*a = 	-wW->data[pos];*/
	b = 	wN->data[pos] + wE->data[pos];
	c = 	-wE->data[pos];
	d =	  wE->data[pos]*(U->data[epos]-U->data[pos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
	}
	
	cp[j] = c/b;
	dp[j] = d/b;

	/*---Middle columns---*/
	for(j=1;j<=ncols-2;j++)
	{
		pos = j*nrows+i;
		wpos = pos-nrows;
		epos = pos+nrows;
		npos = pos-1;
		spos = pos+1;

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
			d -= M->data[pos]*dV->data[pos];
		}

		cp[j] = c / (b-cp[j-1]*a);
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);
	}

	/*---Last column---*/
	pos = j*nrows+i;
	wpos = pos-nrows;
	epos = pos+nrows;
	npos = pos-1;
	spos = pos+1;

	a = 	-wW->data[pos];	
	b = 	wN->data[pos] + wW->data[pos];
	/*c = 	-wE->data[pos];*/
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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
// TDMA vertical line solving - west most columns solver
//------------------------------------------------------
*/
void westColumn_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
	i = 0;j = 0;
	pos = j*nrows+i;
	wpos = pos-nrows;
	nwpos = wpos-1;
	swpos = wpos+1;
	epos = pos+nrows;
	nepos = epos-1;
	sepos = epos+1;
	npos = pos-1;
	spos = pos+1;


	/*a = 	-wN->data[pos];*/
	b = 	wS->data[pos] + wE->data[pos] + wSE->data[pos];
	c = 	-wS->data[pos];
	d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
		+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
		+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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

		a = 	-wN->data[pos];
		b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wNE->data[pos] + wSE->data[pos];
		c = 	-wS->data[pos];
		d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos])
			+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
			+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos]);

		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*dV->data[pos];
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

	a = 	-wN->data[pos];
	b = 	wN->data[pos] + wE->data[pos] + wNE->data[pos];
	/*c = 	-wS->data[pos];*/
	d =	  wN->data[pos]*(U->data[npos]-U->data[pos])
		+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
		+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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
void middleColumn_llin8(struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
	
		/*a = 	0.0f;*/
		b = 	wS->data[pos] + wE->data[pos] + wW->data[pos] + wSE->data[pos] + wSW->data[pos];
		c = 	-wS->data[pos];
		d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
			+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos])
			+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos])
			+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos]);

		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*dV->data[pos];
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
	
			a = 	-wN->data[pos];
			b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos]
				+ wNW->data[pos] + wNE->data[pos] + wSW->data[pos] + wSE->data[pos];
			c = 	-wS->data[pos];
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
				d -= M->data[pos]*dV->data[pos];
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
	
		a = 	-wN->data[pos];
		b = 	wN->data[pos] + wE->data[pos] + wW->data[pos] + wNW->data[pos] + wNE->data[pos];
		/*c = 	-wS->data[pos];*/
		d =	  wN->data[pos]*(U->data[npos]-U->data[pos])
			+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
			+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
			+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos]+dU->data[epos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*dV->data[pos];
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
void eastColumn_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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
	i = 0;j = ncols - 1;
	pos = j*nrows+i;
	wpos = pos-nrows;
	nwpos = wpos-1;
	swpos = wpos+1;
	epos = pos+nrows;
	nepos = epos-1;
	sepos = epos+1;
	npos = pos-1;
	spos = pos+1;

	/*a = 	0.0f;*/
	b = 	wS->data[pos] + wW->data[pos] + wSW->data[pos];
	c = 	-wS->data[pos];
	d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
  		+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
		+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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

		a = 	-wN->data[pos];
		b = 	wN->data[pos] + wS->data[pos] + wW->data[pos] + wNW->data[pos] + wSW->data[pos];
		c = 	-wS->data[pos];
		d =	  wS->data[pos]*(U->data[spos]-U->data[pos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos])
			+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
			+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
			+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*dV->data[pos];
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

	a = 	-wN->data[pos];
	b = 	wN->data[pos] + wW->data[pos] + wNW->data[pos];
	/*c = 	0.0f;*/
	d =	  wN->data[pos]*(U->data[npos]-U->data[pos])
		+ wW->data[pos]*(U->data[wpos]-U->data[pos]+dU->data[wpos])
		+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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
// TDMA horizontal line solving - north most row solver
//------------------------------------------------------
*/
void northRow_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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

	/*a = 	-wW->data[pos];	*/
	b = 	wS->data[pos] + wE->data[pos] + wSE->data[pos];
	c = 	-wE->data[pos];
	d =	  wE->data[pos]*(U->data[epos]-U->data[pos])
		+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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

		a = 	-wW->data[pos];	
		b = 	wS->data[pos] + wE->data[pos] + wW->data[pos] + wSW->data[pos] + wSE->data[pos];
		c = 	-wE->data[pos];
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos])
			+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos])
			+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*dV->data[pos];
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
	b = 	wS->data[pos] + wW->data[pos] + wSW->data[pos];
	/*c = 	-wE->data[pos];*/
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
		+ wSW->data[pos]*(U->data[swpos]-U->data[pos]+dU->data[swpos])
		+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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
void middleRow_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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

		/*a = 	-wW->data[pos];	*/
		b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wNE->data[pos] + wSE->data[pos];
		c = 	-wE->data[pos];
		d =	wE->data[pos]*(U->data[epos]-U->data[pos])
			+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
			+ wSE->data[pos]*(U->data[sepos]-U->data[pos]+dU->data[sepos])
			+ wS->data[pos]*(U->data[spos]-U->data[pos]+dU->data[spos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*dV->data[pos];
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

			a = 	-wW->data[pos];	
			b = 	wN->data[pos] + wS->data[pos] + wE->data[pos] + wW->data[pos]
				+ wNW->data[pos] + wNE->data[pos] + wSW->data[pos] + wSE->data[pos];
			c = 	-wE->data[pos];
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
				d -= M->data[pos]*dV->data[pos];
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
			d -= M->data[pos]*dV->data[pos];
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
// TDMA horizontal line solving - south most rows solver
//------------------------------------------------------
*/
void southRow_llin8(	struct matrixM *U,
			struct matrixM *V,
			struct matrixM *dU,
			struct matrixM *dV,
			struct matrixM *M,
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

	/*a = 	-wW->data[pos];	*/
	b = 	wN->data[pos] + wE->data[pos] + wNE->data[pos];
	c = 	-wE->data[pos];
	d =	wE->data[pos]*(U->data[epos]-U->data[pos])
		+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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

		a = 	-wW->data[pos];	
		b = 	wN->data[pos] + wE->data[pos] + wW->data[pos] + wNW->data[pos] + wNE->data[pos];
		c = 	-wE->data[pos];
		d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
			+ wE->data[pos]*(U->data[epos]-U->data[pos])
			+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
			+ wNE->data[pos]*(U->data[nepos]-U->data[pos]+dU->data[nepos])
			+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);
		
		if(!isnan(Cu->data[pos]))
		{
			b += Du->data[pos];
			d += Cu->data[pos];
			d -= M->data[pos]*dV->data[pos];
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

	a = 	-wW->data[pos];	
	b = 	wN->data[pos] + wW->data[pos] + wNW->data[pos];
	/*c = 	-wE->data[pos];*/
	d =	  wW->data[pos]*(U->data[wpos]-U->data[pos])
		+ wNW->data[pos]*(U->data[nwpos]-U->data[pos]+dU->data[nwpos])
		+ wN->data[pos]*(U->data[npos]-U->data[pos]+dU->data[npos]);

	if(!isnan(Cu->data[pos]))
	{
		b += Du->data[pos];
		d += Cu->data[pos];
		d -= M->data[pos]*dV->data[pos];
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
