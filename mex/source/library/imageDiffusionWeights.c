/*
Functions related to diffusion.

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
#define isnanC(x)((x)!=(x))

void Dver( struct matrixM *D, float *result )
{
	unsigned int frames = 1;
	unsigned int rows = D->dimElems[0];
	unsigned int cols = D->dimElems[1];
	unsigned int i,j,k,pos;

	float A,B;

	if(D->ndims>2)
		frames = D->dimElems[2];

	for(k=0;k<frames;k++)
	{
		for(i=0;i<cols;i++)
		{
			pos = i*rows+k*cols*rows;
		
			A = 0.25f*D->data[pos];
			B = -0.25f*D->data[pos+1];
			result[pos] = A+B;
	
			for(j=1;j<rows-1;j++)
			{
				pos++;
				A = 0.25*D->data[pos-1];
				B = -0.25*D->data[pos+1];
				result[pos] = A+B;
	
			}
	
			pos++;
			A = 0.25*D->data[pos-1];
			B = -0.25*D->data[pos];
			result[pos] = A+B;
	
		}
	}

}

void Dhor( struct matrixM *D, float *result )
{
	unsigned int frames = 1;
	unsigned int rows = D->dimElems[0];
	unsigned int cols = D->dimElems[1];
	unsigned int i,j,k,pos;

	float A,B;

	if(D->ndims>2)
		frames = D->dimElems[2];

	for(k=0;k<frames;k++)
	{
		for(j=0;j<rows;j++)
		{
			pos = j+k*rows*cols;
			A = 0.25*D->data[pos];
			B = -0.25*D->data[pos+rows];
			result[pos] = A + B;
		
			for(i=1;i<cols-1;i++)
			{
				pos+=rows;
	
				A = 0.25*D->data[pos-rows];
				B = -0.25*D->data[pos+rows];
				result[pos] = A + B;
			}
	
			pos+=rows;
			A = 0.25*D->data[pos-rows];
			B = -0.25*D->data[pos];
			result[pos] = A + B;
		}
	}
}

void Calc_wW( struct matrixM *Result, struct matrixM *D, float *ver, float *hor, float *temp, float eps )
{
	unsigned int frames = 1;
	unsigned int rows = D->dimElems[0];
	unsigned int cols = D->dimElems[1];
	unsigned int i,j,k,pos,ppos;

	float A,B;

	if( D->ndims>2 )
		frames = D->dimElems[2];

	for(k=0;k<frames;k++)
	{
		for(j=0;j<rows;j++)
		{
			pos = j+k*rows*cols;
			ppos = j;
			
			for(i=1;i<cols;i++)
			{
				pos+=rows;
				ppos+=rows;
				
				A = D->data[pos] - D->data[pos-rows];
				B = ver[pos] + ver[pos-rows];
				A = A*A;
				B = B*B;

				temp[pos] = A + B;

				/*In case of several "frames" find the maximum*/
				if( (k>=1)&&(temp[pos]>temp[ppos]) )
					temp[ppos] = temp[pos];
			}
		}
	}
	/*Finally calculate the diffusion weight*/
	for(j=0;j<rows;j++)
	{
		pos = j;

		for(i=1;i<cols;i++)
		{
			pos+=rows;
			Result->data[pos] = 1.0f/(float)sqrt(temp[pos]+eps);
		}
		/*Image edges are prone to errors thus more smoothing is needed*/
		/*pos+=rows;
		Result->data[pos] = 4.0f/(float)sqrt(temp[pos]+eps);*/
	}
	
}

void Calc_wN( struct matrixM *Result, struct matrixM *D, float *ver, float *hor, float *temp, float eps )
{
	unsigned int frames = 1;
	unsigned int rows = D->dimElems[0];
	unsigned int cols = D->dimElems[1];
	unsigned int i,j,k,pos,ppos;

	float A,B;

	if( D->ndims>2 )
		frames = D->dimElems[2];

	for(k=0;k<frames;k++)
	{
		for(i=0;i<cols;i++)
		{
			pos = i*rows+k*cols*rows;
			ppos = i*rows;
	
			for(j=1;j<rows;j++)
			{
				pos++;
				ppos++;
				
				A = D->data[pos] - D->data[pos-1];
				B = hor[pos] + hor[pos-1];
				A = A*A;
				B = B*B;

				temp[pos] = A + B;
				/*In case of several "frames" find the maximum*/
				if( (k>=1)&&(temp[pos]>temp[ppos]) )
					temp[ppos]=temp[pos];
			}
		}
	}
	/*Finally calculate the diffusion weight*/
	for(i=0;i<cols;i++)
	{
		pos = i*rows;
		for(j=1;j<rows;j++)
		{
			pos++;
			Result->data[pos] = 1.0f/(float)sqrt(temp[pos]+eps);
		}
		/*Image edges are prone to errors thus more smoothing is needed*/
		/*pos++;
		Result->data[pos] = 4.0f/(float)sqrt(temp[pos]+eps);*/
	}

	/*East most column especially is prone to error*/
	/*pos = (cols-1)*rows;
	for(j=1;j<rows;j++)
	{
		pos++;
		Result->data[pos] = 4.0f/(float)sqrt(temp[pos]+eps);
	}*/
}

void Calc_wE( struct matrixM *Result, struct matrixM *D, float *ver, float *hor, float *temp, float eps )
{
	unsigned int frames = 1;
	unsigned int rows = D->dimElems[0];
	unsigned int cols = D->dimElems[1];
	unsigned int i,j,k,pos,ppos;

	float A, B;

	if( D->ndims>2 )
		frames = D->dimElems[2];

	for(k=0;k<frames;k++)
	{
		for(j=0;j<rows;j++)
		{
			pos = j+k*rows*cols;
			ppos = j;
			
			for(i=0;i<cols-1;i++)
			{
				A = D->data[pos] - D->data[pos+rows];
				B = ver[pos] + ver[pos+rows];
				A = A*A;
				B = B*B;

				temp[pos] =  A + B;
				/*In case of several "frames" find the maximum*/
				if( (k>=1)&&(temp[pos]>temp[ppos]) )
					temp[ppos] = temp[pos];

				pos+=rows;
				ppos+=rows;
			}
		}
	}
	/*Finally calculate the diffusion weight*/
	for(j=0;j<rows;j++)
	{
		pos = j;
		/*Image edges are prone to errors thus more smoothing is needed*/
		/*Result->data[pos] = 4.0f/(float)sqrt(temp[pos]+eps);
		pos+=rows;*/
				
		for(i=0;i<cols-1;i++)
		{
			Result->data[pos] = 1.0f/(float)sqrt(temp[pos]+eps);
			pos+=rows;
		}
	}

}

void Calc_wS( struct matrixM *Result, struct matrixM *D, float *ver, float *hor, float *temp, float eps )
{
	unsigned int frames = 1;
	unsigned int rows = D->dimElems[0];
	unsigned int cols = D->dimElems[1];
	unsigned int i,j,k,pos,ppos;

	float A,B;

	if( D->ndims>2 )
		frames = D->dimElems[2];

	for(k=0;k<frames;k++)
	{
		for(i=0;i<cols;i++)
		{
			pos = i*rows+k*cols*rows;
			ppos = i*rows;
	
			for(j=0;j<rows-1;j++)
			{
				A = D->data[pos] - D->data[pos+1];
				B = hor[pos] + hor[pos+1];
				A = A*A;
				B = B*B;

				temp[pos] = A + B; 
				/*In case of several "frames" find the maximum*/
				if( (k>=1)&&(temp[pos]>temp[ppos]) )
					temp[ppos] = temp[pos];

				pos++;
				ppos++;
			}
		}
	}
	/*Finally calculate the diffusion weight*/
	for(i=0;i<cols;i++)
	{
		pos = i*rows;
		/*Image edges are prone to errors thus more smoothing is needed*/
		/*Result->data[pos] = 4.0f/(float)sqrt(temp[pos]+eps);
		pos++;*/
	
		for(j=0;j<rows-1;j++)
		{
			Result->data[pos] = 1.0f/(float)sqrt(temp[pos]+eps);
			pos++;
		}
	}
	/*East most column especially is prone to error*/
	/*pos = (cols-1)*rows;
	for(j=1;j<rows;j++)
	{
		Result->data[pos] = 4.0f/(float)sqrt(temp[pos]+eps);
		pos++;
	}*/
}

/*
//--------------------------------
//---Calculates diffusion weights 
//--------------------------------
*/
void diffWeights6_2D_c(		struct matrixM *wW,
				struct matrixM *wN,
				struct matrixM *wE,
				struct matrixM *wS,
				struct matrixM *D,
				float eps )
{
	unsigned int frames = 1;
	unsigned int rows = D->dimElems[0];
	unsigned int cols = D->dimElems[1];
	float *ver, *hor, *temp;

	if( D->ndims>2 )
		frames = D->dimElems[2];

	if( (ver = (float*)malloc( rows*cols*frames*sizeof(float) ))==NULL )
		mexErrMsgTxt("Error reserving space for Dhor:ver!");

	if( (hor = (float*)malloc( rows*cols*frames*sizeof(float) ))==NULL )
		mexErrMsgTxt("Error reserving space for Dhor:hor!");

	if( (temp = (float*)malloc( rows*cols*frames*sizeof(float) ))==NULL )
		mexErrMsgTxt("Error reserving space for Dhor:temp!");

	/*Intermediary results*/
	Dver( D, ver );
	Dhor( D, hor );

	/*Finally, the weights*/
	Calc_wW( wW, D, ver, hor, temp, eps );
	Calc_wN( wN, D, ver, hor, temp, eps );
	Calc_wE( wE, D, ver, hor, temp, eps );
	Calc_wS( wS, D, ver, hor, temp, eps );

	free( ver );
	free( hor );
	free( temp );
}

